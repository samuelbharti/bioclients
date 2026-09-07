# Human Phenotype Ontology: term search, term lookup, and the diseases and
# phenotypes HPO associates with a gene.
#
# Ported from genescout/R/tools/hpo.R (gene annotation) and
# a sibling app's api_hpo.R (search and term).
#
# Three endpoints, all keyless. JAX takes HP ids in the path directly and
# answers one term per request, so there is nothing here to batch.
#
# The gene annotation response carries two independent arrays, diseases and
# phenotypes, so both are returned.
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
# genescout carries a list of generic disease words, a tokenizer that drops them,
# and a function that decides whether a gene's HPO diseases are relevant to the
# disease under review. Deciding that "Li-Fraumeni syndrome" answers a question
# about "adrenocortical carcinoma" is the review's judgement, not a property of
# HPO's response, so it stays in the app making the judgement.
#
# A sibling app carries the ontology work built on top of these terms:
# ancestor and descendant walks, Resnik information content, and the propagated
# per-term gene counts. All of it needs the whole DAG rather than an API
# response, and none of it is a client concern.
#
# Endpoint: https://ontology.jax.org/api

HPO_URL <- "https://ontology.jax.org/api"
HPO_WEB <- "https://ontology.jax.org/app/browse/gene"
HPO_TERM_WEB <- "https://hpo.jax.org/browse/term"

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
#' @inherit hpo_search references
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
#' @inherit hpo_search references
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

#' Turn an HPO search response into a table
#'
#' Pure.
#'
#' @param body A parsed HPO `hp/search` response.
#'
#' @return A tibble of `id`, `name`, `definition`, and `descendant_count`, best
#'   match first, which is the order JAX returns. `NULL` when nothing matched.
#'
#' @inherit hpo_search references
#'
#' @examples
#' body <- list(terms = list(list(id = "HP:0001250", name = "Seizure")))
#' hpo_parse_search(body)
#'
#' @export
hpo_parse_search <- function(body) {
  terms <- biohttp::pluck_at(body, "terms")
  if (is.null(terms) || length(terms) == 0) {
    return(NULL)
  }
  tibble::tibble(
    id = col_chr(terms, "id"),
    name = col_chr(terms, "name"),
    definition = col_chr(terms, "definition"),
    # How many terms sit below this one, which is what tells a broad term from
    # a specific one without walking the ontology.
    descendant_count = col_num(terms, "descendantCount")
  )
}

#' Turn an HPO term response into a table
#'
#' Pure.
#'
#' @param body A parsed HPO `hp/terms/<id>` response.
#'
#' @return A one-row tibble of `id`, `name`, `definition`, `comment`,
#'   `descendant_count`, and the list columns `synonyms` and `xrefs`. `NULL`
#'   when the body carries no term.
#'
#' @inherit hpo_search references
#'
#' @examples
#' body <- list(id = "HP:0001250", name = "Seizure", synonyms = list("Seizures"))
#' hpo_parse_term(body)
#'
#' @export
hpo_parse_term <- function(body) {
  id <- chr_at(body, "id")
  if (is.na(id)) {
    return(NULL)
  }
  as_chr_list <- function(key) {
    list(as.character(unlist(
      biohttp::pluck_at(body, key, default = list()),
      use.names = FALSE
    )))
  }
  tibble::tibble(
    id = id,
    name = chr_at(body, "name"),
    definition = chr_at(body, "definition"),
    comment = chr_at(body, "comment"),
    descendant_count = num_at(body, "descendantCount"),
    synonyms = as_chr_list("synonyms"),
    xrefs = as_chr_list("xrefs"),
    source_url = paste0(HPO_TERM_WEB, "/", id)
  )
}

#' Search HPO terms by free text
#'
#' @param query Free text, for example `"seizure"`.
#' @param limit Maximum terms to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [hpo_parse_search()].
#'
#' @references
#' Gargano et al. (2024). The Human Phenotype Ontology in 2024: phenotypes
#' around the world. Nucleic Acids Research 52(D1), D1333-D1346.
#' \doi{10.1093/nar/gkad1005}
#'
#' Service documentation: <https://hpo.jax.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(hpo_search("seizure"))
#' }
#'
#' @export
hpo_search <- function(query, limit = 10, ...) {
  if (biohttp::is_blank(query)) {
    return(biohttp::status_no_data(
      source = "HPO",
      detail = "no search text was supplied"
    ))
  }
  res <- biohttp::get_json(
    HPO_URL,
    path = "hp/search",
    query = list(q = as.character(query), limit = as.integer(limit)),
    source = "HPO",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- hpo_parse_search(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "HPO",
      http = res$http,
      detail = paste0("HPO matched no term for ", query)
    ))
  }
  biohttp::status_ok(data = parsed, source = "HPO", http = res$http)
}

#' Resolve one HP id to its term
#'
#' @param id An HP id, for example `"HP:0001250"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the one-row tibble described in
#'   [hpo_parse_term()].
#'
#' @inherit hpo_search references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(hpo_term("HP:0001250"))
#' }
#'
#' @export
hpo_term <- function(id, ...) {
  term <- toupper(trimws(as.character(id %||% "")))
  # The id goes straight into the path, so anything that is not an HP id is
  # refused rather than turned into a request for some other resource.
  if (!grepl("^HP:[0-9]+$", term)) {
    return(biohttp::status_no_data(
      source = "HPO",
      detail = "an HP term id is required, for example HP:0001250"
    ))
  }
  res <- biohttp::get_json(
    HPO_URL,
    path = paste0("hp/terms/", term),
    source = "HPO",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- hpo_parse_term(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "HPO",
      http = res$http,
      detail = paste0("HPO has no term ", term)
    ))
  }
  biohttp::status_ok(data = parsed, source = "HPO", http = res$http)
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
#' @inherit hpo_search references
#'
#' @examples
#' \donttest{
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
