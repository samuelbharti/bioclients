# Human Phenotype Ontology: the diseases and phenotypes HPO associates with a
# gene.
#
# Ported from genescout/R/tools/hpo.R.
#
# One GET to the JAX ontology API's gene annotation endpoint. The response
# carries two independent arrays, diseases and phenotypes, so both are returned.
#
# THE LOOKUP KEY IS AN NCBI GENE ID.
#
# The path is `network/annotation/NCBIGene:<id>`, and the id has to be the
# numeric Entrez gene id. A symbol, an Ensembl id, or an HGNC id in that slot is
# a request for an entity HPO does not have. Resolve a symbol first with
# [mygene_gene()], which returns `entrez`.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout also carries a list of generic disease words, a tokenizer that drops
# them, and a function that decides whether a gene's HPO diseases are relevant to
# the disease under review. Deciding that "Li-Fraumeni syndrome" answers a
# question about "adrenocortical carcinoma" is the review's judgement, not a
# property of HPO's response, so it stays in the app making the judgement.
#
# Endpoint: https://ontology.jax.org/api

HPO_URL <- "https://ontology.jax.org/api"
HPO_WEB <- "https://ontology.jax.org/app/browse/gene"

#' Turn an HPO gene annotation into a table of diseases
#'
#' Pure.
#'
#' Every row is grounded by an OMIM, ORPHA, or DECIPHER id in `id`, and by a
#' MONDO id in `mondo` where HPO has one.
#'
#' @param body A parsed HPO `network/annotation` response.
#'
#' @return A tibble of `id`, `name`, `mondo`, and `description`, one row per
#'   associated disease. `NULL` when there are none.
#'
#' @examples
#' body <- list(diseases = list(
#'   list(
#'     id = "OMIM:151623",
#'     name = "Li-Fraumeni syndrome",
#'     mondoId = "MONDO:0018875",
#'     description = NULL
#'   )
#' ))
#' hpo_parse_diseases(body)
#'
#' @export
hpo_parse_diseases <- function(body) {
  records <- biohttp::pluck_at(body, "diseases")
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  tibble::tibble(
    id = col_chr(records, "id"),
    name = col_chr(records, "name"),
    mondo = col_chr(records, "mondoId"),
    # HPO sends an explicit null here far more often than it sends prose.
    description = col_chr(records, "description")
  )
}

#' Turn an HPO gene annotation into a table of phenotypes
#'
#' Pure.
#'
#' These are the HPO terms observed in the diseases the gene is associated with,
#' not a separate assertion about the gene.
#'
#' @param body A parsed HPO `network/annotation` response.
#'
#' @return A tibble of `id` and `name`, one row per phenotype term. `NULL` when
#'   there are none.
#'
#' @examples
#' body <- list(phenotypes = list(list(id = "HP:0003002", name = "Breast carcinoma")))
#' hpo_parse_phenotypes(body)
#'
#' @export
hpo_parse_phenotypes <- function(body) {
  records <- biohttp::pluck_at(body, "phenotypes")
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  tibble::tibble(
    id = col_chr(records, "id"),
    name = col_chr(records, "name")
  )
}

#' HPO's annotation for a gene
#'
#' @param entrez An NCBI Gene id, for example `7157`. See the note on the
#'   lookup key in the file header: nothing else works here.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of two tibbles, `diseases`
#'   and `phenotypes`, either of which may be `NULL`. Two tables rather than one,
#'   because they are independent arrays in the response and joining them would
#'   invent a disease-to-phenotype mapping the response does not carry. It also
#'   carries `source_url`. See [hpo_parse_diseases()] and
#'   [hpo_parse_phenotypes()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(hpo_gene_annotation(7157))$diseases
#' }
#'
#' @export
hpo_gene_annotation <- function(entrez, ...) {
  id <- trimws(as.character(entrez %||% ""))
  if (!grepl("^[0-9]+$", id)) {
    return(biohttp::status_no_data(
      source = "HPO",
      detail = "an NCBI Gene id is required for an HPO lookup"
    ))
  }
  res <- biohttp::get_json(
    HPO_URL,
    path = paste0("network/annotation/NCBIGene:", id),
    source = "HPO",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  diseases <- hpo_parse_diseases(res$data)
  phenotypes <- hpo_parse_phenotypes(res$data)
  if (is.null(diseases) && is.null(phenotypes)) {
    return(biohttp::status_no_data(
      source = "HPO",
      http = res$http,
      detail = paste0("HPO has no annotation for NCBIGene:", id)
    ))
  }
  biohttp::status_ok(
    data = list(
      diseases = diseases,
      phenotypes = phenotypes,
      source_url = paste0(HPO_WEB, "/NCBIGene:", id)
    ),
    source = "HPO",
    http = res$http
  )
}
