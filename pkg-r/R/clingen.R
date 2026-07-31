# ClinGen Allele Registry: canonical allele identity.
#
# Ported from multi-variant-reviewer/R/api_clingen_registry.R.
#
# One batch POST of newline-delimited HGVS returns a JSON array in input order,
# each element carrying the Canonical Allele id (CAid) plus the cross-references
# that would otherwise need one lookup each: dbSNP, ClinVar allele and variation
# ids, the MyVariant hg38 id, and GRCh38 coordinates.
#
# NO KEY IS NEEDED. The batch POST is open, which is worth knowing because most
# identity services of this kind are not.
#
# Endpoint: https://reg.clinicalgenome.org/alleles

CLINGEN_REGISTRY_URL <- "https://reg.clinicalgenome.org/alleles"

# The documented rate. Supplied as the default so a caller gets it right without
# having to read the docs.
CLINGEN_THROTTLE <- list(capacity = 3, fill_time_s = 1)

# Pull the GRCh38 chromosomal locus out of the genomic alleles.
clingen_grch38_locus <- function(genomic_alleles) {
  empty <- list(
    chrom = NA_character_,
    pos = NA_integer_,
    ref = NA_character_,
    alt = NA_character_
  )
  for (allele in genomic_alleles) {
    is_grch38 <- identical(allele$referenceGenome, "GRCh38")
    # NC_ is the chromosomal accession. Without it this could be a scaffold.
    is_chromosomal <- any(grepl("^NC_", unlist(allele$hgvs %||% character())))
    if (!is_grch38 || !is_chromosomal || length(allele$coordinates) == 0) {
      next
    }
    coord <- allele$coordinates[[1]]
    return(list(
      chrom = sub(
        "^chr",
        "",
        as.character(allele$chromosome %||% NA_character_),
        ignore.case = TRUE
      ),
      pos = suppressWarnings(as.integer(coord$start %||% NA)),
      ref = as.character(coord$referenceAllele %||% NA_character_),
      alt = as.character(coord$allele %||% NA_character_)
    ))
  }
  empty
}

# The first id out of an externalRecords group, or NA.
clingen_first_record <- function(records, group, field) {
  entries <- records[[group]]
  if (is.null(entries) || length(entries) == 0) {
    return(NA)
  }
  entries[[1]][[field]] %||% NA
}

#' Turn one Allele Registry element into a table row
#'
#' Pure.
#'
#' @section An unresolved input is a row, not a gap:
#' ClinGen answers an input it could not resolve with an object carrying no
#' `@id`, and an `errorType`/`message` instead. That becomes a row with
#' `resolved = FALSE` and the reason, rather than being dropped. Dropping it
#' would shift every later row onto the wrong variant, and would lose the reason
#' the input was rejected.
#'
#' @param element One parsed Allele Registry element.
#'
#' @return A one-row tibble.
#'
#' @examples
#' clingen_parse_allele(list(errorType = "IncorrectHgvsPosition"))
#'
#' @export
clingen_parse_allele <- function(element) {
  caid_url <- element[["@id"]] %||% NA_character_
  if (is.na(caid_url)) {
    return(tibble::tibble(
      resolved = FALSE,
      caid = NA_character_,
      reason = as.character(
        element$message %||% element$errorType %||% "not resolved"
      ),
      rsid = NA_character_,
      clinvar_allele_id = NA_integer_,
      clinvar_variation_id = NA_integer_,
      myvariant_hg38 = NA_character_,
      chrom = NA_character_,
      pos = NA_integer_,
      ref = NA_character_,
      alt = NA_character_,
      title = NA_character_
    ))
  }
  records <- element$externalRecords %||% list()
  rs <- clingen_first_record(records, "dbSNP", "rs")
  locus <- clingen_grch38_locus(element$genomicAlleles %||% list())
  tibble::tibble(
    resolved = TRUE,
    caid = basename(as.character(caid_url)),
    reason = NA_character_,
    rsid = if (is.na(rs)) NA_character_ else paste0("rs", rs),
    clinvar_allele_id = as.integer(
      clingen_first_record(records, "ClinVarAlleles", "alleleId")
    ),
    clinvar_variation_id = as.integer(
      clingen_first_record(records, "ClinVarVariations", "variationId")
    ),
    myvariant_hg38 = as.character(
      clingen_first_record(records, "MyVariantInfo_hg38", "id")
    ),
    chrom = locus$chrom,
    pos = locus$pos,
    ref = locus$ref,
    alt = locus$alt,
    title = as.character(
      (element$communityStandardTitle %||% list(NA_character_))[[1]]
    )
  )
}

#' Turn an Allele Registry batch response into a table
#'
#' Pure. One row per element, in the order returned, which is the order asked.
#'
#' @param body A parsed Allele Registry response, an array.
#'
#' @return A tibble with one row per element.
#'
#' @export
clingen_parse_batch <- function(body) {
  if (is.null(body) || length(body) == 0) {
    return(NULL)
  }
  tibble::as_tibble(do.call(rbind, lapply(body, clingen_parse_allele)))
}

#' Resolve HGVS to canonical allele ids
#'
#' One batched POST. The response comes back in input order.
#'
#' @param hgvs HGVS strings.
#' @param throttle A throttle spec. Defaults to the documented 3 per second.
#' @param ... Passed to [biohttp::perform()] indirectly through the request.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [clingen_parse_allele()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(clingen_alleles("NM_000546.6:c.215C>G"))
#' }
#'
#' @export
clingen_alleles <- function(hgvs, throttle = CLINGEN_THROTTLE, ...) {
  hgvs <- as.character(hgvs %||% character())
  hgvs <- hgvs[!is.na(hgvs) & nzchar(hgvs)]
  if (length(hgvs) == 0) {
    return(biohttp::status_no_data(
      source = "ClinGen",
      detail = "no HGVS strings were supplied"
    ))
  }
  # A raw newline-delimited text body, not JSON, which is what file=hgvs means.
  #
  # This is the one client that builds its own request. biohttp has get_json(),
  # post_json() and get_text(), none of which send a text/plain body, so the
  # request is assembled with httr2 and then handed to biohttp::req_defaults()
  # and biohttp::perform(). Every transport rule still applies: retries, the
  # circuit breaker, the envelope, redacted headers. What is NOT happening is a
  # call to httr2::req_perform(), which would bypass all of it.
  #
  # A post_raw() in biohttp would remove even this.
  req <- httr2::request(CLINGEN_REGISTRY_URL)
  req <- httr2::req_url_query(req, file = "hgvs")
  req <- httr2::req_body_raw(
    req,
    paste(hgvs, collapse = "\n"),
    "text/plain"
  )
  req <- biohttp::req_defaults(req, timeout = 30, throttle = throttle, ...)
  res <- biohttp::perform(req, "ClinGen")
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- clingen_parse_batch(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "ClinGen",
      http = res$http,
      detail = "the Allele Registry returned no records"
    ))
  }
  biohttp::status_ok(data = parsed, source = "ClinGen", http = res$http)
}
