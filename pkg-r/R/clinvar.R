# ClinVar: clinical significance, via NCBI E-utilities.
#
# Ported from variant-reviewer/R/api_clinvar.R.
#
# Two requests, not one: esearch resolves a term to a UID, esummary turns that
# UID into a record. There is no single-call endpoint, so the client makes both
# and the caller sees one envelope.
#
# THE NCBI KEY GOES IN THE QUERY STRING, NOT A HEADER.
#
# E-utilities raises a caller from 3 to 10 requests a second with an `api_key`
# parameter, and offers no header form. A query-string credential needs handling
# a header one does not: it would land in the cache key, it would print with the
# request, and it would ride along in the URL inside a transport error message.
#
# biohttp::secret_query handles all three. It attaches the key at dispatch, so
# nothing built from the request beforehand carries it, and redacts its value
# from any message built from a failure. Ported from genescout's
# fix/secret-redaction-and-ncbi-key.
#
# Docs: https://www.ncbi.nlm.nih.gov/books/NBK25500/

EUTILS_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# The optional key, read from the environment at call time rather than at load,
# so setting it in a running session takes effect without reloading the package.
clinvar_secret_query <- function() {
  key <- Sys.getenv("NCBI_API_KEY")
  if (nzchar(key)) list(api_key = key) else NULL
}

#' Collapse a ClinVar trait set into one condition string
#'
#' Pure.
#'
#' @param germline The `germline_classification` block of an esummary record.
#'
#' @return A single semicolon-separated string, or `NA_character_`.
#'
#' @examples
#' clinvar_conditions(list(trait_set = list(list(trait_name = "RASopathy"))))
#'
#' @export
clinvar_conditions <- function(germline) {
  traits <- biohttp::pluck_at(germline, "trait_set")
  if (is.null(traits) || length(traits) == 0) {
    return(NA_character_)
  }
  names <- vapply(traits, function(trait) chr_at(trait, "trait_name"), "")
  names <- unique(names[!is.na(names) & nzchar(names)])
  if (length(names) == 0) NA_character_ else paste(names, collapse = "; ")
}

#' Turn a ClinVar esummary record into a table
#'
#' Pure. Takes one record from the `result` block of an esummary response, not
#' the whole response, because the record is keyed by UID and the caller already
#' knows which UID it asked for.
#'
#' @param record A parsed ClinVar esummary record.
#' @param uid The ClinVar UID the record came from.
#'
#' @return A one-row tibble with `uid`, `accession`, `title`, `significance`,
#'   `review_status`, `last_evaluated`, and `conditions`. `NULL` when `record`
#'   is absent.
#'
#' @examples
#' record <- list(
#'   accession = "VCV000040389",
#'   title = "NM_004333.6(BRAF):c.1799T>G (p.Val600Gly)",
#'   germline_classification = list(
#'     description = "Pathogenic",
#'     review_status = "reviewed by expert panel"
#'   )
#' )
#' clinvar_parse_record(record, "40389")
#'
#' @export
clinvar_parse_record <- function(record, uid = NA_character_) {
  if (is.null(record)) {
    return(NULL)
  }
  germline <- biohttp::pluck_at(record, "germline_classification")
  tibble::tibble(
    uid = as.character(uid),
    accession = chr_at(record, "accession"),
    title = chr_at(record, "title"),
    significance = chr_at(germline, "description"),
    review_status = chr_at(germline, "review_status"),
    last_evaluated = chr_at(germline, "last_evaluated"),
    conditions = clinvar_conditions(germline)
  )
}

#' Bucket a ClinVar significance string into a coarse category
#'
#' Pure. For colouring and grouping, where the full free-text significance is
#' too granular to plot.
#'
#' The order of the checks matters. A conflicting record usually contains the
#' word "pathogenic" too, so conflicting has to be tested first or every
#' conflicting record reads as pathogenic.
#'
#' @param significance A ClinVar clinical significance string.
#'
#' @return One of `"Conflicting"`, `"Pathogenic / likely"`,
#'   `"Benign / likely"`, `"Uncertain"`, or `"Other"`.
#'
#' @examples
#' clinvar_category("Pathogenic")
#' clinvar_category("Conflicting interpretations of pathogenicity")
#'
#' @export
clinvar_category <- function(significance) {
  value <- tolower(as.character(significance %||% ""))
  if (grepl("conflict", value)) {
    "Conflicting"
  } else if (grepl("pathogenic", value)) {
    "Pathogenic / likely"
  } else if (grepl("benign", value)) {
    "Benign / likely"
  } else if (grepl("uncertain", value)) {
    "Uncertain"
  } else {
    "Other"
  }
}

#' Look up the ClinVar classification for a variant
#'
#' Resolves `term` to a ClinVar UID, then fetches that record. An rsID is the
#' term that works best.
#'
#' @section An rsID does not identify an allele:
#' The same caveat as [gnomad_frequency()]. `7:g.140753336A>T` and
#' `7:g.140753336A>C` share `rs113488022`, so a term that is only an rsID can
#' resolve to a record for the other allele. ClinVar returns the first matching
#' UID, and this function returns that record.
#'
#' @section The NCBI API key:
#' Set `NCBI_API_KEY` and both requests carry it, which raises the rate limit
#' from 3 to 10 requests a second. It is passed as a `secret_query`, so it stays
#' out of the cache key and out of every message. Without one the client works
#' at the lower limit.
#'
#' @param term A search term, usually an rsID or an accession.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [clinvar_parse_record()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(clinvar_classification("rs113488022"))
#' }
#'
#' @export
clinvar_classification <- function(term, ...) {
  if (biohttp::is_blank(term)) {
    return(biohttp::status_no_data(
      source = "ClinVar",
      detail = "no variant identifier was supplied"
    ))
  }
  search <- biohttp::get_json(
    EUTILS_BASE,
    path = "esearch.fcgi",
    query = list(db = "clinvar", term = as.character(term), retmode = "json"),
    source = "ClinVar",
    secret_query = clinvar_secret_query(),
    ...
  )
  if (!isTRUE(search$ok)) {
    return(search)
  }
  ids <- biohttp::pluck_at(search$data, "esearchresult", "idlist")
  if (is.null(ids) || length(ids) == 0) {
    return(biohttp::status_no_data(
      source = "ClinVar",
      http = search$http,
      detail = paste0("no ClinVar record for ", term)
    ))
  }
  uid <- as.character(ids[[1]])

  summary <- biohttp::get_json(
    EUTILS_BASE,
    path = "esummary.fcgi",
    query = list(db = "clinvar", id = uid, retmode = "json"),
    source = "ClinVar",
    secret_query = clinvar_secret_query(),
    ...
  )
  if (!isTRUE(summary$ok)) {
    return(summary)
  }
  parsed <- clinvar_parse_record(
    biohttp::pluck_at(summary$data, "result", uid),
    uid
  )
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "ClinVar",
      http = summary$http,
      detail = paste0("ClinVar returned no summary for UID ", uid)
    ))
  }
  biohttp::status_ok(data = parsed, source = "ClinVar", http = summary$http)
}
