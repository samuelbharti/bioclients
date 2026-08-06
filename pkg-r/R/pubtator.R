# PubTator3: how many articles have a gene tagged as an entity.
#
# Ported from genescout/R/tools/pubtator.R.
#
# ENTITY TAGGING, NOT STRING MATCHING.
#
# This is the point of using PubTator over a Europe PMC symbol count. Europe PMC
# matches the raw text, so counting "MET" or "PIGS" or "SET" returns every paper
# using the English word. PubTator3 counts articles where its pipeline has
# tagged the actual gene, so an ambiguous symbol stops inflating its own count.
#
# PREFER THE ENTREZ ID OVER THE SYMBOL.
#
# `@GENE_7157` names TP53 exactly. `@GENE_TP53` relies on PubTator resolving the
# symbol, which is the thing that goes wrong for an alias. Pass `entrez` when
# there is one; [mygene_gene()] returns it.
#
# A COUNT OF ZERO IS AN ANSWER.
#
# Same rule as europepmc.R: a gene PubTator knows and has tagged nowhere is
# `ok` with a count of 0. Only a missing count is a failure. Note that this is
# the opposite call from impc.R, where 0 is refused, and the difference is
# whether the source can tell "nothing found" apart from "never looked".
# PubTator has searched the whole corpus either way; IMPC may never have
# phenotyped the gene.
#
# Endpoint: https://www.ncbi.nlm.nih.gov/research/pubtator3-api

PUBTATOR_URL <- "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
PUBTATOR_WEB <- "https://www.ncbi.nlm.nih.gov/research/pubtator3"

#' Build the PubTator3 entity token for a gene
#'
#' Prefers the Entrez id, which names the gene exactly, and falls back to the
#' symbol. Exported because the token is what appears in `source_id`, so a
#' caller reconstructing a citation needs the same rule.
#'
#' @param symbol A gene symbol.
#' @param entrez An NCBI Gene id. Used in preference to `symbol` when present.
#'
#' @return A single string such as `"@GENE_7157"`, or `NULL` when neither
#'   argument is usable.
#'
#' @inherit pubtator_gene_literature references
#'
#' @examples
#' pubtator_entity("TP53", 7157)
#' pubtator_entity("TP53")
#'
#' @export
pubtator_entity <- function(symbol = NULL, entrez = NULL) {
  id <- trimws(as.character(entrez %||% ""))
  if (length(id) == 1 && !is.na(id) && grepl("^[0-9]+$", id)) {
    return(paste0("@GENE_", id))
  }
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(NULL)
  }
  paste0("@GENE_", toupper(cleaned))
}

#' Read the article count off a PubTator3 search response
#'
#' Pure. `0` is a real count and comes back as `0L`; only an absent `count` is
#' `NA`.
#'
#' @param body A parsed PubTator3 `search/` response.
#'
#' @return A single integer, or `NA_integer_`.
#'
#' @inherit pubtator_gene_literature references
#'
#' @examples
#' pubtator_parse_count(list(count = 4180))
#' pubtator_parse_count(list(count = 0))
#'
#' @export
pubtator_parse_count <- function(body) {
  raw <- biohttp::pluck_at(body, "count", default = NA)
  if (biohttp::is_blank(raw)) {
    return(NA_integer_)
  }
  suppressWarnings(as.integer(raw))
}

#' Turn PubTator3 search results into a table
#'
#' Pure.
#'
#' @param body A parsed PubTator3 `search/` response.
#'
#' @return A tibble of `pmid`, `title`, `journal`, `authors` (a list column, as
#'   an article has any number), `year`, and `source_url`. `NULL` when there are
#'   no results.
#'
#' @inherit pubtator_gene_literature references
#'
#' @examples
#' body <- list(results = list(list(
#'   pmid = 36197410,
#'   title = "TP53 or Not TP53",
#'   journal = "Clin Cancer Res"
#' )))
#' pubtator_parse_results(body)
#'
#' @export
pubtator_parse_results <- function(body) {
  records <- biohttp::pluck_at(body, "results")
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  pmid <- col_chr(records, "pmid")
  tibble::tibble(
    pmid = pmid,
    title = col_chr(records, "title"),
    journal = col_chr(records, "journal"),
    authors = lapply(records, function(rec) {
      as.character(unlist(
        biohttp::pluck_at(rec, "authors", default = list()),
        use.names = FALSE
      ))
    }),
    year = col_chr(records, "date"),
    source_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid)
  )
}

#' Articles PubTator3 has tagged with a gene
#'
#' @param symbol A gene symbol.
#' @param entrez An NCBI Gene id, preferred over `symbol`. See the note in the
#'   file header.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `count`, `entity`,
#'   `results` (the tibble from [pubtator_parse_results()], which may be
#'   `NULL`), and `source_url`. A count of 0 is `ok`, not `no_data`.
#'
#' @references
#' Wei et al. (2024). PubTator 3.0: an AI-powered literature resource for
#' unlocking biomedical knowledge. Nucleic Acids Research 52(W1), W540-W546.
#' \doi{10.1093/nar/gkae235}
#'
#' Service documentation: <https://www.ncbi.nlm.nih.gov/research/pubtator3/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(pubtator_gene_literature("TP53", 7157))$count
#' }
#'
#' @export
pubtator_gene_literature <- function(symbol = NULL, entrez = NULL, ...) {
  entity <- pubtator_entity(symbol, entrez)
  if (is.null(entity)) {
    return(biohttp::status_no_data(
      source = "PubTator3",
      detail = "no usable gene symbol or Entrez id was supplied"
    ))
  }
  res <- biohttp::get_json(
    PUBTATOR_URL,
    # The trailing slash is load-bearing: without it NCBI redirects, which costs
    # a round trip on every lookup.
    path = "search/",
    query = list(text = entity),
    source = "PubTator3",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  count <- pubtator_parse_count(res$data)
  if (is.na(count)) {
    return(biohttp::status_no_data(
      source = "PubTator3",
      http = res$http,
      detail = "PubTator3 returned no count"
    ))
  }
  biohttp::status_ok(
    data = list(
      count = count,
      entity = entity,
      results = pubtator_parse_results(res$data),
      source_url = paste0(
        PUBTATOR_WEB,
        "?query=",
        utils::URLencode(entity, reserved = TRUE)
      )
    ),
    source = "PubTator3",
    http = res$http
  )
}
