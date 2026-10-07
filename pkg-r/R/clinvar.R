# ClinVar: clinical significance, via NCBI E-utilities.
#
# Ported from variant-reviewer/R/api_clinvar.R.
#
# Two requests, not one: esearch resolves a term to UIDs, esummary turns them
# into records. There is no single-call endpoint, so the client makes both and
# the caller sees one envelope. esummary takes every UID in one call, so a term
# that matches several records still costs two requests.
#
# A TERM CAN MATCH SEVERAL RECORDS, AND THE FIRST IS OFTEN THE WRONG ONE.
#
# esearch is a text search. rs113488022 matches BRAF V600G and V600E, and lists
# V600G first. NM_004333.6:c.1799T>A, the HGVS name of V600E, lists BRAF I208V
# first, a different variant at a different position. rs121913343 lists TP53
# R273G before R273C. Keeping the first UID returned the wrong record in all
# three. That was issue #40.
#
# An SPDI or a VCV accession matches one record. For any other term the caller
# can pass `allele`, a protein change, and the client keeps the record whose own
# `protein_change` field has it. That field lists the change once per transcript
# ("V600E, V512E, ..."), so it is split and matched as a list, not searched for
# in the title. Without `allele` the first record is still returned, and
# `n_matches` says how many records the term matched.
#
# THE NCBI KEY GOES IN THE QUERY STRING, NOT A HEADER.
#
# E-utilities raises a caller from 3 to 10 requests a second with an `api_key`
# parameter, and offers no header form. A query-string credential needs handling
# a header one does not: it would land in the cache key, it would print with the
# request, and it would ride along in the URL inside a transport error message.
#
# The `secret_query` argument biohttp's wrappers take handles all three. It
# attaches the key at dispatch, so nothing built from the request beforehand
# carries it, and redacts its value from any message built from a failure, in
# both the raw and the percent-encoded form a URL carries. Ported from
# genescout's fix/secret-redaction-and-ncbi-key.
#
# Docs: https://www.ncbi.nlm.nih.gov/books/NBK25500/

EUTILS_BASE <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# The optional key, read from the environment at call time rather than at load,
# so setting it in a running session takes effect without reloading the package.
clinvar_secret_query <- function() {
  key <- Sys.getenv("NCBI_API_KEY")
  if (nzchar(key)) list(api_key = key) else NULL
}

# The documented E-utilities rate: 3 requests a second without a key, 10 with
# one. Supplied as the default throttle so a caller gets it right without
# reading the docs, the way CLINGEN_THROTTLE does for the Allele Registry. A
# function rather than a constant because the rate follows the key, and the
# key is read at call time.
clinvar_throttle <- function() {
  per_second <- if (nzchar(Sys.getenv("NCBI_API_KEY"))) 10 else 3
  list(capacity = per_second, fill_time_s = 1)
}

# NCBI asks a caller to say who it is with `tool` and `email` on every
# request, so that a misbehaving client can be contacted before it is blocked.
# Both come from the variables biohttp already builds its User-Agent from, so
# a caller that has identified itself once is identified here too. A blank
# value is omitted rather than sent empty.
clinvar_identity_query <- function() {
  tool <- Sys.getenv("BIOHTTP_CALLER_IDENTITY", "")
  email <- Sys.getenv("BIOHTTP_CONTACT_EMAIL", "")
  out <- list()
  if (nzchar(tool)) {
    out$tool <- tool
  }
  if (nzchar(email)) {
    out$email <- email
  }
  out
}

# How many UIDs esearch is asked for. This is also its default, and it is sent
# so the limit is visible here. An rsID or an HGVS name matches a handful.
CLINVAR_MAX_RECORDS <- 20L

# Three-letter amino acid codes to the one-letter codes ClinVar uses in
# `protein_change`. Ter is a stop, which ClinVar writes as *.
CLINVAR_AMINO_ACIDS <- c(
  ala = "A",
  arg = "R",
  asn = "N",
  asp = "D",
  cys = "C",
  gln = "Q",
  glu = "E",
  gly = "G",
  his = "H",
  ile = "I",
  leu = "L",
  lys = "K",
  met = "M",
  phe = "F",
  pro = "P",
  ser = "S",
  thr = "T",
  trp = "W",
  tyr = "Y",
  val = "V",
  sec = "U",
  pyl = "O",
  ter = "*"
)

# Bring a protein change into the form ClinVar writes: "p.Val600Glu",
# "p.(Val600Glu)" and "v600e" all become "V600E", and any frameshift ends in a
# bare "fs", as in "Q1756fs". A form this does not recognise is returned
# trimmed, and then simply matches no record.
clinvar_protein_change <- function(allele) {
  x <- gsub("[()]", "", sub("^p\\.", "", trimws(as.character(allele))))
  three <- regmatches(
    x,
    regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|fs.*|\\*)$", x)
  )[[1]]
  if (length(three) == 4) {
    from <- CLINVAR_AMINO_ACIDS[tolower(three[[2]])]
    to <- if (startsWith(three[[4]], "fs")) {
      "fs"
    } else if (three[[4]] == "*") {
      "*"
    } else {
      CLINVAR_AMINO_ACIDS[tolower(three[[4]])]
    }
    if (!anyNA(c(from, to))) {
      return(unname(paste0(from, three[[3]], to)))
    }
  }
  one <- regmatches(x, regexec("^([A-Za-z])([0-9]+)([A-Za-z*]|fs.*)$", x))[[1]]
  if (length(one) == 4) {
    to <- if (startsWith(one[[4]], "fs")) "fs" else toupper(one[[4]])
    return(paste0(toupper(one[[2]]), one[[3]], to))
  }
  x
}

# TRUE when an esummary record lists `change` among its protein changes.
clinvar_has_protein_change <- function(record, change) {
  listed <- chr_at(record, "protein_change")
  if (is.na(listed)) {
    return(FALSE)
  }
  change %in% trimws(strsplit(listed, ",", fixed = TRUE)[[1]])
}

#' Collapse a ClinVar trait set into one condition string
#'
#' Pure.
#'
#' @param germline The `germline_classification` block of an esummary record.
#'
#' @return A single semicolon-separated string, or `NA_character_`.
#'
#' @inherit clinvar_classification references
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
#' @inherit clinvar_classification references
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
#' @inherit clinvar_classification references
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
#' Searches ClinVar for `term` and returns the classification of one record.
#'
#' @section A term can match several records:
#' ClinVar's search is a text search, and an rsID or an HGVS name often matches
#' more than one record. `rs113488022` matches BRAF V600G and V600E, and ClinVar
#' lists V600G first. `NM_004333.6:c.1799T>A`, the HGVS name of V600E, lists
#' BRAF I208V first, a different variant at a different position.
#'
#' There are two ways to get the record you mean:
#'
#' * Pass a precise `term`. An SPDI such as `NC_000007.14:140753335:A:T`, or a
#'   VCV accession such as `VCV000013961`, matches one record.
#' * Pass `allele`, the protein change you mean, such as `"V600E"` or
#'   `"p.Val600Glu"`. It is checked against each record's own list of protein
#'   changes, and the first record that has it is returned.
#'
#' Without either, the first record ClinVar lists is returned, and `n_matches`
#' says how many records the term matched. A value above 1 means the row may be
#' for a different variant from the one you meant.
#'
#' Only the first 20 records ClinVar lists are fetched, so `allele` is checked
#' against those. An rsID or an HGVS name matches far fewer.
#'
#' @section The NCBI API key:
#' Set `NCBI_API_KEY` and both requests carry it, which raises the rate limit
#' from 3 to 10 requests a second. It is passed as a `secret_query`, so it stays
#' out of the cache key and out of every message. Without one the client works
#' at the lower limit, and the default `throttle` follows: 3 a second without
#' a key, 10 with one.
#'
#' @section Identifying the caller:
#' NCBI asks every client to send `tool` and `email`. They are read from
#' `BIOHTTP_CALLER_IDENTITY` and `BIOHTTP_CONTACT_EMAIL`, the same variables
#' biohttp builds its User-Agent from, and a blank one is omitted. They are
#' ordinary query parameters, so they are part of the cache key.
#'
#' @param term A search term. An SPDI or a VCV accession matches one record; an
#'   rsID or an HGVS name can match several.
#' @param throttle A throttle spec, see [biohttp::req_defaults()]. Defaults to
#'   the documented E-utilities rate for the key in use.
#' @param allele Optional. The protein change you mean, such as `"V600E"` or
#'   `"p.Val600Glu"`, for a `term` that can match several records. `NULL`, `NA`
#'   and `""` mean no filter.
#' @param ... Passed to [biohttp::get_json()].
#'
#' @return A biohttp envelope whose `data` is a one-row tibble with the columns
#'   of [clinvar_parse_record()] and `n_matches`, the number of records that fit
#'   the request. `no_data` when nothing matches, including when no record has
#'   the protein change in `allele`.
#'
#' @references
#' Landrum et al. (2018). ClinVar: improving access to variant
#' interpretations and supporting evidence. Nucleic Acids Research 46(D1),
#' D1062-D1067. \doi{10.1093/nar/gkx1153}
#'
#' Service documentation: <https://www.ncbi.nlm.nih.gov/clinvar/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(clinvar_classification("NC_000007.14:140753335:A:T"))
#' biohttp::body_or_null(clinvar_classification("rs113488022", allele = "V600E"))
#' }
#'
#' @export
clinvar_classification <- function(
  term,
  throttle = clinvar_throttle(),
  allele = NULL,
  ...
) {
  if (biohttp::is_blank(term)) {
    return(biohttp::status_no_data(
      source = "ClinVar",
      detail = "no variant identifier was supplied"
    ))
  }
  if (length(allele) > 1) {
    stop("allele must be a single protein change", call. = FALSE)
  }
  wanted <- if (biohttp::is_blank(allele)) {
    NULL
  } else {
    clinvar_protein_change(allele)
  }
  identity <- clinvar_identity_query()
  search <- biohttp::get_json(
    EUTILS_BASE,
    path = "esearch.fcgi",
    query = c(
      list(
        db = "clinvar",
        term = as.character(term),
        retmode = "json",
        retmax = CLINVAR_MAX_RECORDS
      ),
      identity
    ),
    source = "ClinVar",
    throttle = throttle,
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
  ids <- as.character(unlist(ids, use.names = FALSE))
  matched <- suppressWarnings(as.integer(
    biohttp::pluck_at(search$data, "esearchresult", "count", default = NA)
  ))
  if (is.na(matched)) {
    matched <- length(ids)
  }

  summary <- biohttp::get_json(
    EUTILS_BASE,
    path = "esummary.fcgi",
    query = c(
      list(db = "clinvar", id = paste(ids, collapse = ","), retmode = "json"),
      identity
    ),
    source = "ClinVar",
    throttle = throttle,
    secret_query = clinvar_secret_query(),
    ...
  )
  if (!isTRUE(summary$ok)) {
    return(summary)
  }
  records <- lapply(ids, function(uid) {
    biohttp::pluck_at(summary$data, "result", uid)
  })
  names(records) <- ids
  records <- Filter(Negate(is.null), records)
  if (length(records) == 0) {
    return(biohttp::status_no_data(
      source = "ClinVar",
      http = summary$http,
      detail = paste0(
        "ClinVar returned no summary for UID ",
        paste(ids, collapse = ", ")
      )
    ))
  }
  if (!is.null(wanted)) {
    has_it <- vapply(
      records,
      clinvar_has_protein_change,
      logical(1),
      change = wanted
    )
    if (!any(has_it)) {
      return(biohttp::status_no_data(
        source = "ClinVar",
        http = summary$http,
        detail = paste0(
          "no ClinVar record for ",
          term,
          " has the protein change ",
          wanted
        )
      ))
    }
    records <- records[has_it]
    matched <- length(records)
  }

  parsed <- clinvar_parse_record(records[[1]], names(records)[[1]])
  parsed$n_matches <- as.integer(matched)
  biohttp::status_ok(data = parsed, source = "ClinVar", http = summary$http)
}
