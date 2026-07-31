# VariantValidator: HGVS validation and normalization.
#
# Ported from multi-variant-reviewer/R/api_variantvalidator.R.
#
# Used as an independent second opinion on an HGVS string. Its validation
# warnings are the real rejection reasons, so they are returned rather than
# collapsed into a pass/fail.
#
# Single-variant endpoint, throttled to the documented 4 per second, so this is a
# targeted check on the rows that need one rather than a bulk pass.
#
# Endpoint: https://rest.variantvalidator.org

VARIANTVALIDATOR_URL <- "https://rest.variantvalidator.org"

# The documented rate, supplied as the default.
VARIANTVALIDATOR_THROTTLE <- list(capacity = 4, fill_time_s = 1)

#' Turn a VariantValidator response into a table
#'
#' Pure.
#'
#' @section The response is keyed by the answer, not by a fixed name:
#' VariantValidator returns an object whose key is the **normalized** variant,
#' which is not known before the call. `flag` and `metadata` are the only fixed
#' keys, so the record is found by removing those rather than by looking up a
#' name. A parser expecting a fixed key finds nothing.
#'
#' @param body A parsed VariantValidator response.
#' @param submitted The HGVS string that was sent.
#'
#' @return A one-row tibble of `resolved`, `submitted`, `normalized`, `gene`,
#'   `protein`, `chrom`, `pos`, `ref`, `alt`, and `warnings`, where `warnings` is
#'   a list column. `NULL` when there is no record at all.
#'
#' @examples
#' body <- list(
#'   flag = "gene_variant",
#'   metadata = list(),
#'   `NM_000546.6:c.215C>G` = list(gene_symbol = "TP53")
#' )
#' variantvalidator_parse(body, "NM_000546.6:c.215C>G")
#'
#' @export
variantvalidator_parse <- function(body, submitted = NA_character_) {
  keys <- setdiff(names(body), c("flag", "metadata"))
  if (length(keys) == 0) {
    return(NULL)
  }
  record <- body[[keys[1]]]
  vcf <- biohttp::pluck_at(
    record,
    "primary_assembly_loci",
    "grch38",
    "vcf",
    default = list()
  )
  tibble::tibble(
    resolved = TRUE,
    submitted = as.character(
      submitted %||% record$submitted_variant %||% keys[1]
    ),
    normalized = keys[1],
    gene = as.character(record$gene_symbol %||% NA_character_),
    protein = as.character(
      biohttp::pluck_at(
        record,
        "hgvs_predicted_protein_consequence",
        "tlr",
        default = NA_character_
      )
    ),
    chrom = sub(
      "^chr",
      "",
      as.character(vcf$chr %||% NA_character_),
      ignore.case = TRUE
    ),
    pos = suppressWarnings(as.integer(vcf$pos %||% NA)),
    ref = as.character(vcf$ref %||% NA_character_),
    alt = as.character(vcf$alt %||% NA_character_),
    # A list column, because the warnings are the useful part and joining them
    # into one string would make them unreadable.
    warnings = list(as.character(
      unlist(record$validation_warnings %||% list())
    ))
  )
}

#' Validate and normalize one HGVS variant
#'
#' @param hgvs An HGVS transcript variant.
#' @param build A genome build.
#' @param throttle A throttle spec. Defaults to the documented 4 per second.
#' @param ... Passed to [biohttp::get_json()].
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [variantvalidator_parse()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(variantvalidator_normalize("NM_000546.6:c.215C>G"))
#' }
#'
#' @export
variantvalidator_normalize <- function(
  hgvs,
  build = "GRCh38",
  throttle = VARIANTVALIDATOR_THROTTLE,
  ...
) {
  if (biohttp::is_blank(hgvs)) {
    return(biohttp::status_no_data(
      source = "VariantValidator",
      detail = "no HGVS string was supplied"
    ))
  }
  res <- biohttp::get_json(
    VARIANTVALIDATOR_URL,
    path = paste(
      "VariantValidator",
      "variantvalidator",
      build,
      as.character(hgvs),
      "all",
      sep = "/"
    ),
    query = list(`content-type` = "application/json"),
    source = "VariantValidator",
    timeout = 30,
    throttle = throttle,
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- variantvalidator_parse(res$data, hgvs)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "VariantValidator",
      http = res$http,
      detail = paste0("VariantValidator returned no record for ", hgvs)
    ))
  }
  biohttp::status_ok(
    data = parsed,
    source = "VariantValidator",
    http = res$http
  )
}
