# PDBe: experimentally determined structures covering a protein.
#
# Ported from genescout/R/tools/pdbe.R.
#
# One GET to the SIFTS best-structures mapping, keyed by UniProt accession. Every
# row is a real PDB accession, so the result is fully grounded.
#
# Endpoint: https://www.ebi.ac.uk/pdbe/api

PDBE_URL <- "https://www.ebi.ac.uk/pdbe/api"
PDBE_WEB <- "https://www.ebi.ac.uk/pdbe/entry/pdb"

#' Turn a PDBe best-structures response into a table
#'
#' Pure.
#'
#' @section One entry per structure, not per chain:
#' The response is keyed by the accession, and each value is a list of **per
#' chain** mappings. A single PDB entry therefore recurs once for every chain it
#' resolved, so 3 structures across 12 chains arrive as 12 records. Counting them
#' raw overstates structural coverage several-fold.
#'
#' Rows are collapsed to distinct `pdb_id`, keeping the first occurrence, which
#' is PDBe's own best (highest coverage) chain for that entry.
#'
#' @param body A parsed PDBe `best_structures` response.
#' @param accession The UniProt accession that was queried.
#'
#' @return A tibble of `pdb_id`, `method`, `resolution`, `coverage`, and
#'   `source_url`, one row per distinct structure. `NULL` when there are none.
#'
#' @examples
#' body <- list(P21359 = list(
#'   list(pdb_id = "7pgp", chain_id = "F", experimental_method = "Electron Microscopy",
#'        resolution = 3.1, coverage = 1),
#'   list(pdb_id = "7pgp", chain_id = "N", experimental_method = "Electron Microscopy",
#'        resolution = 3.1, coverage = 1)
#' ))
#' pdbe_parse_structures(body, "P21359")
#'
#' @export
pdbe_parse_structures <- function(body, accession = NA_character_) {
  records <- biohttp::pluck_at(body, accession)
  # PDBe keys the object by the accession it resolved, which is not always the
  # exact string that was asked for. Fall back to the only value present.
  if (is.null(records) && length(body) > 0) {
    records <- body[[1]]
  }
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  pdb_id <- col_chr(records, "pdb_id")
  keep <- !is.na(pdb_id) & nzchar(pdb_id) & !duplicated(pdb_id)
  if (!any(keep)) {
    return(NULL)
  }
  tibble::tibble(
    pdb_id = pdb_id[keep],
    method = col_chr(records, "experimental_method")[keep],
    resolution = col_num(records, "resolution")[keep],
    coverage = col_num(records, "coverage")[keep],
    source_url = paste0(PDBE_WEB, "/", pdb_id[keep])
  )
}

#' Experimental structures for a UniProt accession
#'
#' @param accession A UniProt accession, for example `"P21359"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [pdbe_parse_structures()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(pdbe_structures("P21359"))
#' }
#'
#' @export
pdbe_structures <- function(accession, ...) {
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "PDBe",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  res <- biohttp::get_json(
    PDBE_URL,
    path = paste0("mappings/best_structures/", cleaned),
    source = "PDBe",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- pdbe_parse_structures(res$data, cleaned)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "PDBe",
      http = res$http,
      detail = paste0("PDBe has no structures for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "PDBe", http = res$http)
}
