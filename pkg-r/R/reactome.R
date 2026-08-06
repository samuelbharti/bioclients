# Reactome: the curated pathways a gene participates in.
#
# Ported from genescout/R/tools/pathways.R.
#
# One GET to the HGNC mapping, which takes the gene symbol directly rather than
# an accession, so no identifier resolution is needed first.
#
# THE RESPONSE IS A BARE ARRAY.
#
# Reactome answers with a top-level JSON array, not an object wrapping one. A
# parser reaching for a `results` or `pathways` key finds nothing and reports
# every gene as having no pathways.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout narrows the list to the pathways that count toward disease
# relevance: those Reactome flags as disease-associated, plus any whose name
# shares a token with the pathways the review already cares about. Deciding that
# "Regulation of RAS by GAPs" answers a question about "RAS/MAPK" is the
# review's judgement, so `in_disease` is returned as the flag Reactome sets and
# the narrowing stays in the app.
#
# Endpoint: https://reactome.org/ContentService

REACTOME_URL <- "https://reactome.org/ContentService"
REACTOME_WEB <- "https://reactome.org/content/detail"

#' Turn a Reactome pathway array into a table
#'
#' Pure.
#'
#' @param body A parsed Reactome pathways response, which is a bare array.
#'
#' @return A tibble of `pathway_id`, `name`, `in_disease`, and `source_url`, one
#'   row per pathway. `NULL` when the gene has none.
#'
#' @inherit reactome_pathways references
#'
#' @examples
#' body <- list(list(
#'   stId = "R-HSA-5658442",
#'   displayName = "Regulation of RAS by GAPs",
#'   isInDisease = FALSE
#' ))
#' reactome_parse_pathways(body)
#'
#' @export
reactome_parse_pathways <- function(body) {
  if (is.null(body) || length(body) == 0) {
    return(NULL)
  }
  pathway_id <- col_chr(body, "stId")
  keep <- !is.na(pathway_id) & nzchar(pathway_id)
  if (!any(keep)) {
    return(NULL)
  }
  records <- body[keep]
  tibble::tibble(
    pathway_id = pathway_id[keep],
    name = trimws(col_chr(records, "displayName")),
    # Reactome's own flag for a pathway that represents disease biology. It is
    # a property of the pathway, not a judgement about this gene.
    in_disease = vapply(
      records,
      function(rec) {
        isTRUE(biohttp::pluck_at(rec, "isInDisease", default = FALSE))
      },
      logical(1)
    ),
    source_url = paste0(REACTOME_WEB, "/", pathway_id[keep])
  )
}

#' Reactome pathways for a gene symbol
#'
#' @param symbol A gene symbol, for example `"NF1"`.
#' @param species The NCBI taxon id. Defaults to human.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [reactome_parse_pathways()].
#'
#' @references
#' Milacic et al. (2024). The Reactome Pathway Knowledgebase 2024.
#' Nucleic Acids Research 52(D1), D672-D678. \doi{10.1093/nar/gkad1025}
#'
#' Service documentation: <https://reactome.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(reactome_pathways("NF1"))
#' }
#'
#' @export
reactome_pathways <- function(symbol, species = 9606, ...) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "Reactome",
      detail = "no usable gene symbol was supplied"
    ))
  }
  sym <- toupper(cleaned)
  res <- biohttp::get_json(
    REACTOME_URL,
    path = paste0("data/mapping/HGNC/", sym, "/pathways"),
    query = list(species = species),
    source = "Reactome",
    ...
  )
  # Reactome answers 404 for a gene it maps but has no pathways for, which is
  # the source saying it has nothing rather than the call failing.
  if (identical(res$http, 404L)) {
    return(biohttp::status_no_data(
      source = "Reactome",
      http = res$http,
      detail = paste0("Reactome has no pathways for ", sym)
    ))
  }
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- reactome_parse_pathways(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Reactome",
      http = res$http,
      detail = paste0("Reactome has no pathways for ", sym)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Reactome", http = res$http)
}
