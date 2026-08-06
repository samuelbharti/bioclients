# QuickGO: Gene Ontology annotations for a gene product.
#
# Ported from genescout/R/tools/quickgo.R.
#
# One GET to annotation/search, keyed by UniProt accession as
# `geneProductId=UniProtKB:<acc>`. Resolve a symbol first with [mygene_gene()]
# and read its `uniprot` column.
#
# ONE ROW PER EVIDENCE LINE, NOT PER TERM.
#
# GO repeats a term once for every piece of evidence supporting it, so a gene
# with 20 distinct functions can arrive as 100 rows. Counting the raw rows
# overstates how many things a gene is annotated to. Rows are collapsed to
# distinct `go_id`, keeping the first evidence line for each.
#
# `includeFields=goName` is not optional in practice. Without it QuickGO returns
# the GO id and no label, which is enough to count terms and not enough to read
# them.
#
# Endpoint: https://www.ebi.ac.uk/QuickGO/services

QUICKGO_URL <- "https://www.ebi.ac.uk/QuickGO/services"
QUICKGO_WEB <- "https://www.ebi.ac.uk/QuickGO/term"

# The three GO aspects. biological_process is the default because it is the axis
# that says what a gene does rather than where it is or what it binds.
QUICKGO_ASPECTS <- c(
  "biological_process",
  "molecular_function",
  "cellular_component"
)

#' Turn a QuickGO annotation response into a table
#'
#' Pure.
#'
#' @param body A parsed QuickGO `annotation/search` response.
#'
#' @return A tibble of `go_id`, `go_name`, `evidence`, `reference`, `aspect`,
#'   and `source_url`, one row per distinct GO term. `NULL` when there are none.
#'
#' @inherit quickgo_annotations references
#'
#' @examples
#' body <- list(results = list(list(
#'   goId = "GO:0001937",
#'   goName = "negative regulation of endothelial cell proliferation",
#'   goEvidence = "IMP",
#'   reference = "PMID:16648142"
#' )))
#' quickgo_parse_annotations(body)
#'
#' @export
quickgo_parse_annotations <- function(body) {
  records <- biohttp::pluck_at(body, "results")
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  go_id <- col_chr(records, "goId")
  # A row whose goId is not a GO id is not an annotation, whatever else it is.
  keep <- !is.na(go_id) & grepl("^GO:", go_id) & !duplicated(go_id)
  if (!any(keep)) {
    return(NULL)
  }
  tibble::tibble(
    go_id = go_id[keep],
    go_name = trimws(col_chr(records, "goName")[keep]),
    evidence = col_chr(records, "goEvidence")[keep],
    reference = col_chr(records, "reference")[keep],
    aspect = col_chr(records, "goAspect")[keep],
    source_url = paste0(QUICKGO_WEB, "/", go_id[keep])
  )
}

#' GO annotations for a UniProt accession
#'
#' @param accession A UniProt accession, for example `"P21359"`.
#' @param aspect One of `"biological_process"`, `"molecular_function"`, or
#'   `"cellular_component"`.
#' @param limit Maximum annotation rows to ask for. Remember these are evidence
#'   lines rather than terms, so this is not a term count.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [quickgo_parse_annotations()].
#'
#' @references
#' Binns et al. (2009). QuickGO: a web-based tool for Gene Ontology
#' searching. Bioinformatics 25(22), 3045-3046.
#' \doi{10.1093/bioinformatics/btp536}
#'
#' Service documentation: <https://www.ebi.ac.uk/QuickGO/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(quickgo_annotations("P21359"))
#' }
#'
#' @export
quickgo_annotations <- function(
  accession,
  aspect = c(
    "biological_process",
    "molecular_function",
    "cellular_component"
  ),
  limit = 100,
  ...
) {
  aspect <- match.arg(aspect)
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "QuickGO",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  res <- biohttp::get_json(
    QUICKGO_URL,
    path = "annotation/search",
    query = list(
      geneProductId = paste0("UniProtKB:", cleaned),
      aspect = aspect,
      # Without this QuickGO sends ids and no labels.
      includeFields = "goName",
      limit = as.integer(limit)
    ),
    source = "QuickGO",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- quickgo_parse_annotations(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "QuickGO",
      http = res$http,
      detail = paste0("QuickGO has no ", aspect, " annotations for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "QuickGO", http = res$http)
}
