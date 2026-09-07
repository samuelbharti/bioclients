# MyGene.info: gene identity and annotation.
#
# Ported from genescout/R/tools/mygene.R, which is the better of the family's two
# copies. variant-reviewer/R/api_mygene.R is the other and had no batch path.
#
# Endpoint: https://mygene.info/v3

MYGENE_BASE <- "https://mygene.info/v3"

# The fields fetched by both the single and the batch path, kept in one place so
# the two cannot drift.
#
# HGNC is upper case, and that is not a typo. Every other field here is lower
# case, but MyGene names this one `HGNC` in both the request and the response,
# so asking for `hgnc` returns nothing and reads as "this gene has no HGNC id".
# The id matters because Monarch's gene endpoints are keyed on it and nothing
# else in this table substitutes.
MYGENE_FIELDS <- paste(
  "name",
  "symbol",
  "entrezgene",
  "ensembl.gene",
  "uniprot",
  "HGNC",
  "type_of_gene",
  "summary",
  sep = ","
)

# The fields the batch POST matches each query against, so mixed identifier
# types resolve in one request without a per-token prefix.
MYGENE_BATCH_SCOPES <- "symbol,alias,ensembl.gene,entrezgene,retired"

# The documented ceiling on identifiers per batch POST. MyGene answers a
# larger body with an error rather than truncating it, so a gene list longer
# than this has to be chunked and the chunks merged back in input order.
MYGENE_BATCH <- 1000L

# Build the query term, detecting Ensembl gene and Entrez ids so they resolve
# precisely rather than as free-text symbol matches.
mygene_query_term <- function(symbol) {
  if (grepl("^ENSG\\d+", symbol, ignore.case = TRUE)) {
    paste0("ensembl.gene:", symbol)
  } else if (grepl("^\\d+$", symbol)) {
    paste0("entrezgene:", symbol)
  } else {
    symbol
  }
}

# A MyGene field can be a scalar or a list of several mappings. Take the first.
mygene_first <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  if (length(x) == 0) NA_character_ else as.character(x[[1]])
}

#' Choose the best MyGene hit for a queried token
#'
#' MyGene's `_score` can rank a fuzzy alias or retired match to a **different**
#' gene above the exact symbol match. Querying `"TTN"` with the alias and retired
#' scopes returns TTR (entrez 7276, score 19.07) ahead of TTN (7273, score
#' 18.29). Taking the top-scored hit therefore returns the wrong gene, quietly.
#'
#' So a hit whose official symbol equals the token wins, case-insensitively.
#' Only when none does, which is the case for a deliberate alias such as
#' `"p53"`, does MyGene's own best-scored hit get used.
#'
#' This is the single most important thing in this file. It is exported so a
#' caller assembling its own hits can apply the same rule.
#'
#' @param hits A list of hit records from a MyGene response.
#' @param token The identifier that was queried.
#'
#' @return One hit record, or `NULL` when `hits` is empty.
#'
#' @inherit mygene_gene references
#'
#' @examples
#' hits <- list(
#'   list(symbol = "TTR", entrezgene = "7276"),
#'   list(symbol = "TTN", entrezgene = "7273")
#' )
#' mygene_pick_hit(hits, "TTN")$entrezgene
#'
#' @export
mygene_pick_hit <- function(hits, token) {
  if (is.null(hits) || length(hits) == 0) {
    return(NULL)
  }
  for (hit in hits) {
    symbol <- biohttp::pluck_at(hit, "symbol", default = "")
    if (nzchar(symbol) && identical(toupper(symbol), toupper(token))) {
      return(hit)
    }
  }
  hits[[1]]
}

#' Turn MyGene hits into a gene table
#'
#' Pure. Takes an already-parsed response body and never touches the network, so
#' it is tested directly against a stored response.
#'
#' @param body A parsed MyGene `/query` response, the whole body including
#'   `hits`.
#' @param symbol The identifier that was queried, used to break the scoring tie
#'   described in [mygene_pick_hit()] and as the fallback symbol.
#'
#' @return A one-row tibble with `symbol`, `name`, `summary`, `entrez`,
#'   `ensembl_gene`, `uniprot`, `hgnc`, and `type_of_gene`. `NULL` when the body
#'   carries no usable hit.
#'
#' @inherit mygene_gene references
#'
#' @examples
#' body <- list(hits = list(list(
#'   symbol = "TP53",
#'   name = "tumor protein p53",
#'   entrezgene = "7157"
#' )))
#' mygene_parse_hits(body, "TP53")
#'
#' @export
mygene_parse_hits <- function(body, symbol = NA_character_) {
  hits <- biohttp::pluck_at(body, "hits")
  hit <- mygene_pick_hit(hits, symbol)
  mygene_row(hit, symbol)
}

# One hit record to one tibble row. NULL for a missing or explicitly notfound
# hit, so a caller never has to tell "no gene" apart from "wrong gene".
mygene_row <- function(hit, fallback_symbol = NA_character_) {
  if (is.null(hit)) {
    return(NULL)
  }
  if (isTRUE(biohttp::pluck_at(hit, "notfound", default = FALSE))) {
    return(NULL)
  }
  tibble::tibble(
    symbol = as.character(
      biohttp::pluck_at(hit, "symbol", default = fallback_symbol)
    ),
    name = chr_at(hit, "name"),
    summary = chr_at(hit, "summary"),
    entrez = as.character(
      biohttp::pluck_at(hit, "entrezgene", default = NA_character_)
    ),
    ensembl_gene = mygene_first(biohttp::pluck_at(hit, "ensembl", "gene")),
    uniprot = mygene_first(
      biohttp::pluck_at(hit, "uniprot", "Swiss-Prot")
    ),
    # Bare digits, the way MyGene sends it. Monarch wants the CURIE form, which
    # monarch_hgnc_id() builds; keeping the raw value here means this column
    # does not bake in one consumer's formatting.
    hgnc = chr_at(hit, "HGNC"),
    type_of_gene = chr_at(hit, "type_of_gene")
  )
}

#' Turn a MyGene batch response into a gene table
#'
#' Pure. The batch `/query` POST returns a flat array where each element echoes
#' its input `query`. An ambiguous query yields several elements, best `_score`
#' first, and an unmatched one yields an element with `notfound = true`.
#'
#' Hits are grouped by the echoed query and resolved through [mygene_pick_hit()],
#' then mapped back onto `symbols` **in input order**, so a caller zips the
#' result onto its input by position. An unmatched or invalid token yields a row
#' of `NA` rather than being dropped, because a shorter table would silently
#' shift every row after it.
#'
#' @param body A parsed MyGene batch response, a flat list of hit records.
#' @param symbols The identifiers that were queried, in the order asked.
#'
#' @return A tibble with one row per entry in `symbols`, same order.
#'
#' @inherit mygene_gene references
#'
#' @examples
#' body <- list(
#'   list(query = "TP53", symbol = "TP53", entrezgene = "7157"),
#'   list(query = "NOPE", notfound = TRUE)
#' )
#' mygene_parse_batch(body, c("TP53", "NOPE"))
#'
#' @export
mygene_parse_batch <- function(body, symbols) {
  by_query <- list()
  for (hit in body) {
    query <- biohttp::pluck_at(hit, "query")
    if (biohttp::is_blank(query)) {
      next
    }
    # Keep EVERY hit for a query, in MyGene's best-score-first order, so
    # mygene_pick_hit() can override the score with an exact symbol match.
    by_query[[query]] <- c(by_query[[query]], list(hit))
  }
  rows <- lapply(symbols, function(symbol) {
    cleaned <- clean_symbol(symbol)
    if (is.null(cleaned)) {
      return(mygene_empty_row(as.character(symbol)))
    }
    row <- mygene_row(mygene_pick_hit(by_query[[cleaned]], cleaned), cleaned)
    row %||% mygene_empty_row(cleaned)
  })
  do.call(rbind, rows)
}

# A placeholder row for a token that resolved to nothing. Keeps the output the
# same length as the input.
mygene_empty_row <- function(symbol) {
  tibble::tibble(
    symbol = as.character(symbol),
    name = NA_character_,
    summary = NA_character_,
    entrez = NA_character_,
    ensembl_gene = NA_character_,
    uniprot = NA_character_,
    hgnc = NA_character_,
    type_of_gene = NA_character_
  )
}

#' Look up one gene
#'
#' @param symbol A gene symbol, Ensembl gene id, or Entrez id.
#' @param species Passed through to MyGene.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [mygene_parse_hits()] for the columns.
#'
#' @references
#' Xin et al. (2016). High-performance web services for querying gene and
#' variant annotation. Genome Biology 17, 91.
#' \doi{10.1186/s13059-016-0953-9}
#'
#' Service documentation: <https://mygene.info/>
#'
#' @examples
#' \donttest{
#' res <- mygene_gene("TP53")
#' biohttp::body_or_null(res)
#' }
#'
#' @export
mygene_gene <- function(symbol, species = "human", ...) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "MyGene",
      detail = "no usable gene identifier was supplied"
    ))
  }
  res <- biohttp::get_json(
    MYGENE_BASE,
    path = "query",
    query = list(
      q = mygene_query_term(cleaned),
      species = species,
      # Ask for several candidates rather than the top score alone, so
      # mygene_pick_hit() has something to choose between.
      size = 5,
      fields = MYGENE_FIELDS
    ),
    source = "MyGene",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- mygene_parse_hits(res$data, cleaned)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "MyGene",
      http = res$http,
      detail = paste0("no MyGene hit for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "MyGene", http = res$http)
}

#' Look up many genes in one request
#'
#' The batch POST, which is the reason this client is worth installing. An
#' N-symbol list is one round trip rather than N, and the round trip is where
#' essentially all the time goes.
#'
#' MyGene takes at most 1000 identifiers per POST, see `MYGENE_BATCH`. A
#' longer list is chunked, the chunks are dispatched through
#' [biohttp::post_json_many()], and the hits are merged back onto `symbols`
#' in input order. A chunk that failed yields a row of `NA` per identifier
#' rather than taking the whole call down, following [gnomad_constraints()];
#' only when every chunk failed is the failing envelope returned.
#'
#' @param symbols Gene symbols, Ensembl gene ids, or Entrez ids.
#' @inheritParams mygene_gene
#' @param chunk_size Identifiers per request, at most `MYGENE_BATCH`.
#' @param ... Passed to [biohttp::post_json_many()].
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry in
#'   `symbols`, in the same order.
#'
#' @inherit mygene_gene references
#'
#' @examples
#' \donttest{
#' res <- mygene_genes(c("TP53", "BRCA1", "EGFR"))
#' biohttp::body_or_null(res)
#' }
#'
#' @export
mygene_genes <- function(
  symbols,
  species = "human",
  chunk_size = MYGENE_BATCH,
  ...
) {
  cleaned <- vapply(
    symbols,
    function(symbol) clean_symbol(symbol) %||% NA_character_,
    character(1),
    USE.NAMES = FALSE
  )
  usable <- unique(cleaned[!is.na(cleaned)])
  if (length(usable) == 0) {
    return(biohttp::status_no_data(
      source = "MyGene",
      detail = "no usable gene identifiers were supplied"
    ))
  }
  if (chunk_size < 1 || chunk_size > MYGENE_BATCH) {
    stop(
      "chunk_size must be between 1 and ",
      MYGENE_BATCH,
      ", the MyGene limit per request",
      call. = FALSE
    )
  }
  chunks <- split(usable, ceiling(seq_along(usable) / chunk_size))
  bodies <- lapply(chunks, function(chunk) {
    list(
      q = as.list(chunk),
      scopes = MYGENE_BATCH_SCOPES,
      fields = MYGENE_FIELDS,
      species = species
    )
  })
  results <- biohttp::post_json_many(
    paste0(MYGENE_BASE, "/query"),
    bodies = bodies,
    source = "MyGene",
    ...
  )
  answered <- Filter(function(res) isTRUE(res$ok), results)
  if (length(answered) == 0) {
    return(results[[1]])
  }
  # Every answered chunk is a flat array of hits echoing its query, so the
  # arrays concatenate into one and the batch parser maps them back onto the
  # caller's input. A hit from a failed chunk is simply absent, and its
  # identifier gets the NA row the parser gives any unmatched token.
  hits <- unlist(
    lapply(answered, function(res) res$data),
    recursive = FALSE,
    use.names = FALSE
  )
  biohttp::status_ok(
    data = mygene_parse_batch(hits, symbols),
    source = "MyGene"
  )
}
