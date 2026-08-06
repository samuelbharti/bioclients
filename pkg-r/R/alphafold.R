# AlphaFold DB: predicted protein structure.
#
# Ported from variant-reviewer/R/api_alphafold.R.
#
# The metadata lookup is kept separate from the coordinate download on purpose.
# A caller that wants to know whether a model exists, and how confident it is,
# should not pull a multi-megabyte structure file to find out. Only a viewer
# needs the coordinates.
#
# Endpoint: https://alphafold.ebi.ac.uk/api/prediction

ALPHAFOLD_URL <- "https://alphafold.ebi.ac.uk/api/prediction"
ALPHAFOLD_WEB <- "https://alphafold.ebi.ac.uk/entry"

#' Turn an AlphaFold prediction response into a table
#'
#' Pure.
#'
#' @section The one-element array:
#' The API returns a JSON **array** of model records, not an object, even for a
#' single accession. Reading it as an object yields `NULL` for every field while
#' looking like a successful parse. This handles either shape.
#'
#' @param body A parsed AlphaFold prediction response.
#' @param accession The UniProt accession that was queried.
#'
#' @return A one-row tibble of `accession`, `pdb_url`, `cif_url`, `version`,
#'   `mean_plddt`, and `source_url`. `NULL` when there is no model.
#'
#' @examples
#' body <- list(list(
#'   pdbUrl = "https://alphafold.ebi.ac.uk/files/AF-P15056-F1-model_v6.pdb",
#'   latestVersion = 6, globalMetricValue = 66.38
#' ))
#' alphafold_parse_model(body, "P15056")
#'
#' @export
alphafold_parse_model <- function(body, accession = NA_character_) {
  # An unnamed list is the array shape; a named one is already the record.
  entry <- if (is.list(body) && length(body) >= 1 && is.null(names(body))) {
    body[[1]]
  } else {
    body
  }
  pdb_url <- biohttp::pluck_at(entry, "pdbUrl")
  if (biohttp::is_blank(pdb_url)) {
    return(NULL)
  }
  tibble::tibble(
    accession = as.character(accession),
    pdb_url = as.character(pdb_url),
    cif_url = chr_at(entry, "cifUrl"),
    version = num_at(entry, "latestVersion"),
    # AlphaFold's global metric is mean pLDDT, 0-100. Higher is more confident.
    mean_plddt = num_at(entry, "globalMetricValue"),
    source_url = paste0(ALPHAFOLD_WEB, "/", accession)
  )
}

#' The predicted structure for a UniProt accession
#'
#' Metadata only. No coordinates are downloaded; use `pdb_url` from the result
#' for that, and only when something is actually going to render it.
#'
#' @param accession A UniProt accession, for example `"P15056"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [alphafold_parse_model()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(alphafold_model("P15056"))
#' }
#'
#' @export
alphafold_model <- function(accession, ...) {
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "AlphaFold",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  res <- biohttp::get_json(
    ALPHAFOLD_URL,
    path = cleaned,
    source = "AlphaFold",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- alphafold_parse_model(res$data, cleaned)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "AlphaFold",
      http = res$http,
      detail = paste0("no AlphaFold model for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "AlphaFold", http = res$http)
}
