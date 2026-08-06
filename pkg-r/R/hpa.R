# Human Protein Atlas: curated per-gene annotation.
#
# Ported from genescout/R/tools/hpa.R.
#
# HPA's per-gene record carries curated "Protein class" and "Disease involvement"
# tags. This client returns them as they are.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout also keeps a list of which of those tags it counts as
# disease-relevant, and a function that counts the matches into a signal. Which
# tags count is an editorial judgement about what makes a gene interesting, and
# that belongs to the app making the judgement. Porting the list down here would
# mean every consumer inherited genescout's opinion, and changing that opinion
# would need a release of this package.
#
# Endpoint: https://www.proteinatlas.org/<ensembl>.json

HPA_URL <- "https://www.proteinatlas.org"

# HPA fields are string-or-array depending on how many values there are, so a
# reader that assumes either shape breaks on the other.
hpa_as_character <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  as.character(unlist(x, use.names = FALSE))
}

#' Turn an HPA gene record into a table
#'
#' Pure.
#'
#' `protein_class` and `disease_involvement` are list columns, because a gene
#' carries any number of tags and flattening them to one string would make them
#' unusable without re-splitting.
#'
#' @param body A parsed HPA gene record.
#' @param ensembl The Ensembl gene id that was queried.
#'
#' @return A one-row tibble of `symbol`, `ensembl`, `uniprot`, `protein_class`,
#'   `disease_involvement`, and `source_url`. `NULL` when the record is empty.
#'
#' @inherit hpa_gene references
#'
#' @examples
#' body <- list(
#'   Gene = "TP53",
#'   Ensembl = "ENSG00000141510",
#'   Uniprot = list("P04637"),
#'   `Protein class` = list("Cancer-related genes", "Transcription factors"),
#'   `Disease involvement` = list("Tumor suppressor")
#' )
#' hpa_parse_gene(body, "ENSG00000141510")
#'
#' @export
hpa_parse_gene <- function(body, ensembl = NA_character_) {
  if (is.null(body) || length(body) == 0) {
    return(NULL)
  }
  symbol <- chr_at(body, "Gene")
  if (is.na(symbol)) {
    return(NULL)
  }
  tibble::tibble(
    symbol = symbol,
    ensembl = as.character(
      biohttp::pluck_at(body, "Ensembl", default = ensembl)
    ),
    uniprot = list(hpa_as_character(biohttp::pluck_at(body, "Uniprot"))),
    protein_class = list(
      hpa_as_character(biohttp::pluck_at(body, "Protein class"))
    ),
    disease_involvement = list(
      hpa_as_character(biohttp::pluck_at(body, "Disease involvement"))
    ),
    source_url = paste0(HPA_URL, "/", ensembl)
  )
}

#' The Human Protein Atlas record for a gene
#'
#' Keyed by Ensembl gene id, which is what HPA's per-gene JSON path takes.
#' Resolve a symbol first with [mygene_gene()].
#'
#' @param ensembl An Ensembl gene id, for example `"ENSG00000141510"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [hpa_parse_gene()].
#'
#' @references
#' Uhlen et al. (2015). Tissue-based map of the human proteome.
#' Science 347(6220), 1260419. \doi{10.1126/science.1260419}
#'
#' Service documentation: <https://www.proteinatlas.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(hpa_gene("ENSG00000141510"))
#' }
#'
#' @export
hpa_gene <- function(ensembl, ...) {
  id <- trimws(as.character(ensembl %||% ""))
  # HPA's path is the id with .json appended, so a malformed id would be a
  # request for an arbitrary page rather than a lookup.
  if (!grepl("^ENSG[0-9]+$", id)) {
    return(biohttp::status_no_data(
      source = "HPA",
      detail = "an Ensembl gene id is required for an HPA lookup"
    ))
  }
  res <- biohttp::get_json(
    HPA_URL,
    path = paste0(id, ".json"),
    source = "HPA",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- hpa_parse_gene(res$data, id)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "HPA",
      http = res$http,
      detail = paste0("HPA has no record for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "HPA", http = res$http)
}
