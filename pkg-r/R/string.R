# STRING: protein-protein interaction.
#
# Ported from variant-reviewer/R/api_string.R (interaction partners for one gene)
# and genescout/R/tools/string.R (the within-list network, plus the identifier
# map). Both are here because they answer different questions.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout's string_connectivity(), which counts how many other genes in a list
# each gene connects to, is not ported. It is cohort-relative prioritization
# maths over edges this client already returns, and its own comment describes it
# as a signal that "only nudges a connected gene up". That is ranking, and
# ranking stays in the app.
#
# Endpoint: https://string-db.org/api

STRING_URL <- "https://string-db.org/api"
STRING_WEB <- "https://string-db.org/cgi/network"
STRING_HUMAN <- 9606L

# required_score is on STRING's 0-1000 scale. 700 is its own "high confidence"
# threshold.
STRING_MIN_SCORE <- 700L

# Identifiers per network call, so the query URL stays bounded. A caller that
# exceeds it gets `truncated = TRUE` and a count rather than a silently short
# answer, because a gene that was never queried must not be mistaken for a gene
# measured to have no partners.
STRING_MAX_NODES <- 500L

#' Turn STRING interaction-partner rows into a table
#'
#' Pure. Sorted by combined score, strongest first.
#'
#' @param body A parsed STRING `interaction_partners` response, a JSON array.
#'
#' @return A tibble of `partner`, `score`, and the four evidence channels
#'   `experimental`, `database`, `coexpression`, and `textmining`. `NULL` when
#'   there are no partners.
#'
#' @inherit string_partners references
#'
#' @examples
#' body <- list(list(
#'   preferredName_B = "SFN", score = 0.999,
#'   escore = 0.981, dscore = 0.75, ascore = 0, tscore = 0.859
#' ))
#' string_parse_partners(body)
#'
#' @export
string_parse_partners <- function(body) {
  if (is.null(body) || !is.list(body) || length(body) == 0) {
    return(NULL)
  }
  out <- tibble::tibble(
    partner = col_chr(body, "preferredName_B"),
    score = col_num(body, "score"),
    experimental = col_num(body, "escore"),
    database = col_num(body, "dscore"),
    coexpression = col_num(body, "ascore"),
    textmining = col_num(body, "tscore")
  )
  out[order(-out$score), , drop = FALSE]
}

#' Turn a STRING network response into an edge table
#'
#' Pure. Edge endpoints come back in STRING's `preferredName` space; see
#' [string_reconcile_edges()] for why that matters.
#'
#' @param body A parsed STRING `network` response, a JSON array of edges.
#'
#' @return A tibble of `gene_a`, `gene_b`, and `score`, upper-cased. A zero-row
#'   tibble when there are no edges, because "no high-confidence edges" is a real
#'   answer about a set rather than an absence of one.
#'
#' @inherit string_partners references
#'
#' @examples
#' body <- list(list(
#'   preferredName_A = "TP53", preferredName_B = "NF1", score = 0.88
#' ))
#' string_parse_network(body)
#'
#' @export
string_parse_network <- function(body) {
  empty <- tibble::tibble(
    gene_a = character(),
    gene_b = character(),
    score = numeric()
  )
  if (is.null(body) || !is.list(body) || length(body) == 0) {
    return(empty)
  }
  gene_a <- toupper(col_chr(body, "preferredName_A"))
  gene_b <- toupper(col_chr(body, "preferredName_B"))
  score <- col_num(body, "score")
  keep <- !is.na(gene_a) & !is.na(gene_b) & nzchar(gene_a) & nzchar(gene_b)
  tibble::tibble(
    gene_a = gene_a[keep],
    gene_b = gene_b[keep],
    score = score[keep]
  )
}

#' Turn a STRING identifier-map response into a table
#'
#' Pure. One row per queried identifier.
#'
#' @param body A parsed STRING `get_string_ids` response, a JSON array.
#'
#' @return A tibble of `query`, `preferred`, and `string_id`, with `query` and
#'   `preferred` upper-cased.
#'
#' @inherit string_partners references
#'
#' @examples
#' body <- list(list(
#'   queryItem = "SEPTIN9", preferredName = "SEPT9",
#'   stringId = "9606.ENSP00000329125"
#' ))
#' string_parse_ids(body)
#'
#' @export
string_parse_ids <- function(body) {
  empty <- tibble::tibble(
    query = character(),
    preferred = character(),
    string_id = character()
  )
  if (is.null(body) || !is.list(body) || length(body) == 0) {
    return(empty)
  }
  query <- toupper(col_chr(body, "queryItem"))
  preferred <- toupper(col_chr(body, "preferredName"))
  keep <- !is.na(query) & nzchar(query) & !is.na(preferred) & nzchar(preferred)
  out <- tibble::tibble(
    query = query[keep],
    preferred = preferred[keep],
    string_id = col_chr(body, "stringId")[keep]
  )
  out[!duplicated(out$query), , drop = FALSE]
}

#' Rewrite edge endpoints back into the queried symbol space
#'
#' Pure, and it exists because of a real and quiet failure.
#'
#' STRING's canonical `preferredName` lags HGNC. Query `SEPTIN9` and the edges
#' come back naming `SEPT9`, and the `network` endpoint does **not** echo the
#' query term. A caller matching edges against the symbols it asked about
#' therefore finds nothing for that gene and records it as having no partners,
#' which is indistinguishable from a real isolate.
#'
#' Passing the map from [string_map_ids()] through here translates the endpoints
#' back. An endpoint with no mapping is left alone, and an empty map leaves the
#' edges untouched, so reconciliation can never turn a good answer into a worse
#' one.
#'
#' @param edges An edge tibble from [string_parse_network()].
#' @param id_map An identifier map from [string_parse_ids()], or `NULL`.
#'
#' @return The edge tibble, with endpoints translated where a mapping existed.
#'
#' @inherit string_partners references
#'
#' @examples
#' edges <- tibble::tibble(gene_a = "SEPT9", gene_b = "TP53", score = 0.9)
#' map <- tibble::tibble(
#'   query = "SEPTIN9", preferred = "SEPT9", string_id = "9606.ENSP00000329125"
#' )
#' string_reconcile_edges(edges, map)
#'
#' @export
string_reconcile_edges <- function(edges, id_map) {
  if (is.null(id_map) || nrow(id_map) == 0 || nrow(edges) == 0) {
    return(edges)
  }
  lookup <- stats::setNames(id_map$query, id_map$preferred)
  translate <- function(x) {
    hit <- unname(lookup[x])
    ifelse(is.na(hit), x, hit)
  }
  edges$gene_a <- translate(edges$gene_a)
  edges$gene_b <- translate(edges$gene_b)
  edges
}

# Clean a symbol set down to what STRING will accept, upper-cased and unique.
string_symbols <- function(symbols) {
  syms <- toupper(trimws(as.character(symbols %||% character())))
  unique(syms[
    !is.na(syms) & nzchar(syms) & grepl("^[A-Za-z0-9._-]+$", syms)
  ])
}

#' Interaction partners for one gene
#'
#' @section On STRING's content type:
#' STRING serves JSON as `text/json` rather than `application/json`. A client
#' that trusts the content-type header rejects a perfectly good body. `biohttp`
#' parses with `check_type = FALSE`, so this works, but it is the reason not to
#' "tidy" that up.
#'
#' @param symbol A gene symbol.
#' @param species An NCBI taxon id. Defaults to human.
#' @param limit How many partners to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [string_parse_partners()].
#'
#' @references
#' Szklarczyk et al. (2023). The STRING database in 2023: protein-protein
#' association networks and functional enrichment analyses for any sequenced
#' genome of interest. Nucleic Acids Research 51(D1), D638-D646.
#' \doi{10.1093/nar/gkac1000}
#'
#' Service documentation: <https://string-db.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(string_partners("TP53"))
#' }
#'
#' @export
string_partners <- function(
  symbol,
  species = STRING_HUMAN,
  limit = 25,
  ...
) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "STRING",
      detail = "no usable gene symbol was supplied"
    ))
  }
  res <- biohttp::get_json(
    STRING_URL,
    path = "json/interaction_partners",
    query = list(identifiers = cleaned, species = species, limit = limit),
    source = "STRING",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- string_parse_partners(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "STRING",
      http = res$http,
      detail = paste0("STRING found no partners for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "STRING", http = res$http)
}

#' The STRING identifier map for a set of symbols
#'
#' Maps each queried symbol to STRING's own `preferredName`. Needed to interpret
#' network edges; see [string_reconcile_edges()].
#'
#' @param symbols Gene symbols.
#' @inheritParams string_partners
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [string_parse_ids()].
#'
#' @inherit string_partners references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(string_map_ids(c("SEPTIN9", "TP53")))
#' }
#'
#' @export
string_map_ids <- function(symbols, species = STRING_HUMAN, ...) {
  syms <- string_symbols(symbols)
  if (length(syms) == 0) {
    return(biohttp::status_no_data(
      source = "STRING",
      detail = "no usable gene symbols were supplied"
    ))
  }
  res <- biohttp::get_json(
    STRING_URL,
    path = "json/get_string_ids",
    query = list(
      identifiers = paste(syms, collapse = "\r"),
      species = species,
      limit = 1
    ),
    source = "STRING",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  biohttp::status_ok(
    data = string_parse_ids(res$data),
    source = "STRING",
    http = res$http
  )
}

#' The interaction network within a set of genes
#'
#' High-confidence edges among the genes you asked about, with the endpoints
#' reconciled back to your symbols.
#'
#' @section Knowing what was asked:
#' The result carries `queried`, the exact symbol set STRING was asked about
#' after the 500-identifier cap, plus `truncated` and `n_dropped`. That is
#' what lets a caller tell a gene measured to have no partners from a gene that
#' was never sent. Reporting the second as the first would invent a negative
#' result.
#'
#' @param symbols Gene symbols. At least two, since one gene has no network.
#' @param required_score STRING's 0-1000 combined-score threshold.
#' @param reconcile Whether to fetch the identifier map and translate edge
#'   endpoints back into your symbols. Costs one extra request and is skipped
#'   automatically when there are no edges to translate.
#' @inheritParams string_partners
#'
#' @return A biohttp envelope whose `data` is a list of `edges` (see
#'   [string_parse_network()]), `queried`, `n_query`, `truncated`, and
#'   `n_dropped`.
#'
#' @inherit string_partners references
#'
#' @examples
#' \donttest{
#' res <- string_network(c("TP53", "NF1", "EGFR"))
#' biohttp::body_or_null(res)$edges
#' }
#'
#' @export
string_network <- function(
  symbols,
  species = STRING_HUMAN,
  required_score = STRING_MIN_SCORE,
  reconcile = TRUE,
  ...
) {
  syms <- string_symbols(symbols)
  if (length(syms) < 2) {
    return(biohttp::status_no_data(
      source = "STRING",
      detail = "a network needs at least two usable gene symbols"
    ))
  }
  n_all <- length(syms)
  truncated <- n_all > STRING_MAX_NODES
  if (truncated) {
    syms <- syms[seq_len(STRING_MAX_NODES)]
  }
  res <- biohttp::get_json(
    STRING_URL,
    path = "json/network",
    query = list(
      identifiers = paste(syms, collapse = "\r"),
      species = species,
      required_score = required_score
    ),
    source = "STRING",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  edges <- string_parse_network(res$data)
  # Only worth a second request when there is something to translate.
  if (isTRUE(reconcile) && nrow(edges) > 0) {
    mapped <- string_map_ids(syms, species, ...)
    if (isTRUE(mapped$ok)) {
      edges <- string_reconcile_edges(edges, mapped$data)
    }
  }
  biohttp::status_ok(
    data = list(
      edges = edges,
      queried = syms,
      n_query = length(syms),
      truncated = truncated,
      n_dropped = if (truncated) n_all - STRING_MAX_NODES else 0L
    ),
    source = "STRING",
    http = res$http
  )
}
