# JensenLab DISEASES: disease-to-gene associations with a confidence score.
#
# Ported from genescout/R/tools/diseases.R.
#
# DISEASES publishes its associations in separate channels. Two are queried here:
# `Knowledge`, which is curated, and `Textmining`, which is mined from the
# literature. Each answers with its own confidence score on a 0 to 5 scale.
#
# THE RESPONSE IS AN ARRAY WHOSE FIRST ELEMENT IS THE DATA.
#
# The body is a JSON array, and the associations are a map inside `[[1]]`, keyed
# by Ensembl protein id. Element two is an empty object. Reading the top level as
# the list of records therefore finds two things, neither of which is an
# association, and reports zero genes for every disease without failing.
#
# THE QUERY TYPE CODES ARE NOT SELF-DESCRIBING.
#
# `type1 = -26` means the id is a Disease Ontology term and `type2 = 9606` means
# only human genes are wanted. Both are positional magic numbers in the API, so
# they are named constants here.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# The DOID itself. genescout resolves it upstream from the disease context, and
# turning a disease name into an ontology id is a resolution step with its own
# judgement in it.
#
# Endpoint: https://api.jensenlab.org

DISEASES_URL <- "https://api.jensenlab.org"
DISEASES_WEB <- "https://diseases.jensenlab.org"

# Curated knowledge and mined literature. Named, because a caller reading the
# result needs to know which channel a score came from.
DISEASES_CHANNELS <- c("Knowledge", "Textmining")

# -26 selects a Disease Ontology term, 9606 selects human.
DISEASES_DISEASE_TYPE <- -26L
DISEASES_HUMAN_TAXON <- 9606L

#' Turn one DISEASES channel response into a table
#'
#' Pure.
#'
#' @param body A parsed response from one DISEASES channel.
#' @param channel The channel the body came from, recorded in the result.
#'
#' @return A tibble of `symbol`, `protein`, `score`, and `channel`, one row per
#'   associated gene. `NULL` when the channel has no associations. Rows with no
#'   symbol or a non-finite score are dropped, because neither is usable.
#'
#' @inherit diseases_channel references
#'
#' @examples
#' body <- list(
#'   list(ENSP00000351015 = list(name = "NF1", score = 5)),
#'   list()
#' )
#' diseases_parse_channel(body, "Knowledge")
#'
#' @export
diseases_parse_channel <- function(body, channel = NA_character_) {
  # See the file header: the associations are one level down, not at the top.
  entries <- if (is.list(body) && length(body) >= 1) body[[1]] else NULL
  if (is.null(entries) || length(entries) == 0) {
    return(NULL)
  }
  symbol <- col_chr(entries, "name")
  # A score that is not a number is dropped by the is.finite() guard below. The
  # coercion warning that gets there first names neither the field nor the
  # service, so it is noise on a row that is already being discarded.
  score <- suppressWarnings(col_num(entries, "score"))
  keep <- !is.na(symbol) & nzchar(symbol) & is.finite(score)
  if (!any(keep)) {
    return(NULL)
  }
  # The map is keyed by Ensembl protein id, which is the only grounding the
  # response carries for a row. vapply over a named list keeps those names on
  # every column, so they are stripped here rather than in each caller.
  protein <- names(entries) %||% rep(NA_character_, length(entries))
  tibble::tibble(
    symbol = unname(symbol[keep]),
    protein = unname(protein[keep]),
    score = unname(score[keep]),
    channel = channel
  )
}

#' Combine DISEASES channels, keeping the strongest score per gene
#'
#' Pure.
#'
#' @section This is a combining convention, not the source's own answer:
#' DISEASES scores each channel separately and does not publish a combined
#' figure. Taking the maximum treats curated and mined evidence as
#' interchangeable, which suits a screen looking for any support at all and
#' suits nothing that cares where the support came from. Use
#' [diseases_channel()] and keep the channels apart if that distinction
#' matters.
#'
#' @param tables A list of tibbles from [diseases_parse_channel()]. `NULL`
#'   entries are ignored.
#'
#' @return A tibble of `symbol`, `score`, and `channel`, ordered by score,
#'   highest first. `channel` names the channel the winning score came from.
#'   `NULL` when there is nothing to combine.
#'
#' @inherit diseases_channel references
#'
#' @examples
#' diseases_merge_channels(list(
#'   tibble::tibble(symbol = "NF1", protein = NA, score = 5, channel = "Knowledge"),
#'   tibble::tibble(symbol = "NF1", protein = NA, score = 2, channel = "Textmining")
#' ))
#'
#' @export
diseases_merge_channels <- function(tables) {
  tables <- Filter(function(x) !is.null(x) && nrow(x) > 0, tables)
  if (length(tables) == 0) {
    return(NULL)
  }
  symbol <- unlist(lapply(tables, function(x) x$symbol), use.names = FALSE)
  score <- unlist(lapply(tables, function(x) x$score), use.names = FALSE)
  channel <- unlist(lapply(tables, function(x) x$channel), use.names = FALSE)
  best <- vapply(
    split(seq_along(symbol), symbol),
    function(idx) {
      idx[[which.max(score[idx])]]
    },
    integer(1)
  )
  ord <- order(-score[best])
  tibble::tibble(
    symbol = unname(symbol[best][ord]),
    score = unname(score[best][ord]),
    channel = unname(channel[best][ord])
  )
}

#' Genes DISEASES associates with a disease, from one channel
#'
#' @param doid A Disease Ontology id, for example `"DOID:0060293"`.
#' @param channel `"Knowledge"` for curated associations or `"Textmining"` for
#'   mined ones.
#' @param limit Maximum genes to ask for.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [diseases_parse_channel()].
#'
#' @references
#' Pletscher-Frankild et al. (2015). DISEASES: text mining and data
#' integration of disease-gene associations. Methods 74, 83-89.
#' \doi{10.1016/j.ymeth.2014.11.020}
#'
#' Service documentation: <https://diseases.jensenlab.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(diseases_channel("DOID:0060293", "Knowledge"))
#' }
#'
#' @export
diseases_channel <- function(
  doid,
  channel = c("Knowledge", "Textmining"),
  limit = 300,
  ...
) {
  channel <- match.arg(channel)
  if (biohttp::is_blank(doid)) {
    return(biohttp::status_no_data(
      source = "DISEASES",
      detail = "a Disease Ontology id is required for a DISEASES lookup"
    ))
  }
  res <- biohttp::get_json(
    DISEASES_URL,
    path = channel,
    query = diseases_query(doid, limit),
    source = "DISEASES",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- diseases_parse_channel(res$data, channel)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "DISEASES",
      http = res$http,
      detail = paste0("DISEASES has no ", channel, " genes for ", doid)
    ))
  }
  biohttp::status_ok(data = parsed, source = "DISEASES", http = res$http)
}

diseases_query <- function(doid, limit) {
  list(
    type1 = DISEASES_DISEASE_TYPE,
    id1 = as.character(doid),
    type2 = DISEASES_HUMAN_TAXON,
    limit = limit,
    format = "json"
  )
}

#' Genes DISEASES associates with a disease, across both channels
#'
#' Queries `Knowledge` and `Textmining` together and keeps the strongest score
#' per gene. Read the note on [diseases_merge_channels()] before relying on the
#' combined score.
#'
#' @inheritParams diseases_channel
#'
#' @section A failed channel fails the call:
#' If either channel errors, its envelope is returned rather than a combined
#' score built from the other one. A result that silently dropped the curated
#' channel would look like a weaker association rather than a partial answer,
#' and the envelope has no status that says "half of this is missing". Call
#' [diseases_channel()] per channel to handle the halves separately.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [diseases_merge_channels()], with `source_url` added.
#'
#' @inherit diseases_channel references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(diseases_gene_associations("DOID:0060293"))
#' }
#'
#' @export
diseases_gene_associations <- function(doid, limit = 300, ...) {
  if (biohttp::is_blank(doid)) {
    return(biohttp::status_no_data(
      source = "DISEASES",
      detail = "a Disease Ontology id is required for a DISEASES lookup"
    ))
  }
  results <- biohttp::get_json_many(
    DISEASES_URL,
    path = as.list(DISEASES_CHANNELS),
    queries = rep(list(diseases_query(doid, limit)), length(DISEASES_CHANNELS)),
    source = "DISEASES",
    ...
  )
  failed <- Filter(function(r) !isTRUE(r$ok), results)
  if (length(failed) > 0) {
    return(failed[[1]])
  }
  tables <- lapply(seq_along(results), function(i) {
    diseases_parse_channel(results[[i]]$data, DISEASES_CHANNELS[[i]])
  })
  merged <- diseases_merge_channels(tables)
  http <- results[[1]]$http
  if (is.null(merged)) {
    return(biohttp::status_no_data(
      source = "DISEASES",
      http = http,
      detail = paste0("DISEASES has no genes for ", doid)
    ))
  }
  merged$source_url <- paste0(
    DISEASES_WEB,
    "/Entity?type1=",
    DISEASES_DISEASE_TYPE,
    "&type2=",
    DISEASES_HUMAN_TAXON,
    "&id1=",
    utils::URLencode(as.character(doid), reserved = TRUE)
  )
  biohttp::status_ok(data = merged, source = "DISEASES", http = http)
}
