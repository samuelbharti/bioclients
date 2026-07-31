# DGIdb: curated drug-gene interactions.
#
# Ported from genescout/R/tools/dgidb.R and gene-list-builder/R/source_dgidb.R.
#
# The count of curated interactions is a druggability SIGNAL for prioritization,
# never a clinical call. This package reports the count; what an app does with it
# is the app's business.
#
# Endpoint: https://dgidb.org/api/graphql

DGIDB_URL <- "https://dgidb.org/api/graphql"
DGIDB_WEB <- "https://dgidb.org/genes"

# The query is already batch-shaped: `names` takes a list. genescout and
# gene-list-builder both pass a single symbol, which spends one round trip per
# gene for no reason. dgidb_genes() passes the whole list.
DGIDB_QUERY <- paste(
  "query($names: [String!]) {",
  "  genes(names: $names) {",
  "    nodes { name conceptId interactions { interactionScore } }",
  "  }",
  "}",
  sep = "\n"
)

#' Turn a DGIdb genes response into a table
#'
#' Pure. Rows come back in `symbols` order, one per input.
#'
#' @section A real zero is not a miss:
#' The distinction this parser exists to preserve. A gene DGIdb knows about with
#' no recorded interactions has `interaction_count = 0`: that is an answer, and
#' it means the gene is not currently druggable. A gene DGIdb has never heard of
#' has `interaction_count = NA`: that is an absence of evidence.
#'
#' Collapsing the two would tell a caller that an unknown gene is known not to be
#' druggable, which is a different and much stronger claim than the data
#' supports.
#'
#' @param body A parsed DGIdb GraphQL response body.
#' @param symbols The gene symbols that were queried, in the order asked.
#'
#' @return A tibble of `symbol`, `concept_id`, `interaction_count`, and
#'   `source_url`, one row per entry in `symbols`.
#'
#' @examples
#' body <- list(data = list(genes = list(nodes = list(
#'   list(name = "NF1", conceptId = "hgnc:7765", interactions = list(list(), list()))
#' ))))
#' dgidb_parse_genes(body, "NF1")
#'
#' @export
dgidb_parse_genes <- function(body, symbols) {
  nodes <- biohttp::pluck_at(body, "data", "genes", "nodes")
  wanted <- toupper(trimws(as.character(symbols)))
  found <- if (is.null(nodes)) character() else toupper(col_chr(nodes, "name"))

  rows <- lapply(seq_along(wanted), function(i) {
    at <- match(wanted[i], found)
    if (is.na(at)) {
      # DGIdb has never heard of this gene. NA, not 0.
      return(tibble::tibble(
        symbol = wanted[i],
        concept_id = NA_character_,
        interaction_count = NA_integer_,
        source_url = NA_character_
      ))
    }
    node <- nodes[[at]]
    concept <- biohttp::pluck_at(node, "conceptId", default = wanted[i])
    tibble::tibble(
      symbol = wanted[i],
      concept_id = as.character(concept),
      # A gene present with an empty interactions list is a real zero.
      interaction_count = length(
        biohttp::pluck_at(node, "interactions", default = list())
      ),
      source_url = paste0(DGIDB_WEB, "/", concept)
    )
  })
  do.call(rbind, rows)
}

#' Drug-gene interaction counts for many genes
#'
#' One request for the whole list, because DGIdb's `genes` query already takes an
#' array.
#'
#' @param symbols Gene symbols.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry in
#'   `symbols`, in the same order. See [dgidb_parse_genes()] for the columns and
#'   for why an unknown gene is `NA` rather than `0`.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(dgidb_genes(c("NF1", "BRAF")))
#' }
#'
#' @export
dgidb_genes <- function(symbols, ...) {
  cleaned <- vapply(
    symbols,
    function(symbol) clean_symbol(symbol) %||% NA_character_,
    character(1),
    USE.NAMES = FALSE
  )
  usable <- toupper(cleaned[!is.na(cleaned)])
  if (length(usable) == 0) {
    return(biohttp::status_no_data(
      source = "DGIdb",
      detail = "no usable gene symbols were supplied"
    ))
  }
  res <- biohttp::post_json(
    DGIDB_URL,
    body = list(query = DGIDB_QUERY, variables = list(names = as.list(usable))),
    source = "DGIdb",
    ...
  )
  bad <- biohttp::graphql_error(res, "DGIdb")
  if (!is.null(bad)) {
    return(bad)
  }
  biohttp::status_ok(
    data = dgidb_parse_genes(res$data, toupper(cleaned)),
    source = "DGIdb",
    http = res$http
  )
}

#' Drug-gene interaction count for one gene
#'
#' A thin wrapper over [dgidb_genes()].
#'
#' @param symbol A gene symbol.
#' @inheritParams dgidb_genes
#'
#' @return A biohttp envelope whose `data` is a one-row tibble.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(dgidb_gene("NF1"))
#' }
#'
#' @export
dgidb_gene <- function(symbol, ...) {
  dgidb_genes(symbol, ...)
}
