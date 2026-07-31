# Europe PMC: literature search and hit counts.
#
# Ported from genescout/R/tools/literature.R (search and count) and
# variant-reviewer/R/api_europepmc.R (the quoted query builder and the citation
# fields), which had diverged into two clients against one endpoint.
#
# A COUNT OF ZERO IS AN ANSWER, NOT A MISS.
#
# [europepmc_count()] returns `ok` with a count of 0 when a query genuinely
# matches nothing. That is the source saying "no papers", which is a real
# measurement and different from the source being unreachable. Folding it into
# `no_data` would make a gene with no literature indistinguishable from an
# outage, and a caller that treats absence as evidence would then be reading an
# outage as evidence.
#
# [europepmc_search()] does report `no_data` for an empty result list, because
# there are no rows to hand back. The two are not inconsistent: one is being
# asked how many, the other is being asked for which.
#
# Compare the opposite case in impc.R, where a count of 0 is refused precisely
# because it cannot be told apart from "never tested".
#
# Endpoint: https://www.ebi.ac.uk/europepmc/webservices/rest

EUROPEPMC_URL <- "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUROPEPMC_WEB <- "https://europepmc.org"

#' Build a Europe PMC query from terms
#'
#' Each term is quoted so a multi-word term matches as a phrase rather than as
#' loose words, and the terms are ANDed. Quoting also stops a term containing a
#' space or a colon from being read as query syntax.
#'
#' @param ... Terms, for example a gene symbol and an rsID. Blank terms are
#'   dropped.
#'
#' @return A single query string, or `NULL` when nothing usable was given.
#'
#' @examples
#' europepmc_query("BRAF")
#' europepmc_query("BRAF", "V600E")
#'
#' @export
europepmc_query <- function(...) {
  terms <- trimws(as.character(c(...)))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0) {
    return(NULL)
  }
  # A quote inside a term would close the phrase early and change the query.
  terms <- gsub("\"", "", terms, fixed = TRUE)
  paste(paste0("\"", terms, "\""), collapse = " AND ")
}

#' Turn Europe PMC search results into a table
#'
#' Pure.
#'
#' `source_id` is the citable identifier: `PMID:<n>` where there is a PMID, and
#' `<source>:<id>` otherwise, because a preprint or a patent record has no PMID
#' but is still groundable.
#'
#' @param body A parsed Europe PMC `search` response.
#'
#' @return A tibble of `title`, `authors`, `year`, `journal`, `pmid`, `doi`,
#'   `cited_by`, `source`, `source_id`, and `source_url`, newest first, which is
#'   the order Europe PMC returns. `NULL` when nothing matched.
#'
#' @examples
#' body <- list(resultList = list(result = list(list(
#'   id = "36197410",
#'   source = "MED",
#'   pmid = "36197410",
#'   title = "TP53 or Not TP53",
#'   authorString = "Green SD.",
#'   journalTitle = "Clin Cancer Res",
#'   pubYear = "2022"
#' ))))
#' europepmc_parse_results(body)
#'
#' @export
europepmc_parse_results <- function(body) {
  records <- biohttp::pluck_at(body, "resultList", "result") %||% body
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  pmid <- col_chr(records, "pmid")
  source <- col_chr(records, "source")
  id <- col_chr(records, "id")
  tibble::tibble(
    title = col_chr(records, "title"),
    authors = col_chr(records, "authorString"),
    year = col_chr(records, "pubYear"),
    journal = col_chr(records, "journalTitle"),
    pmid = pmid,
    doi = col_chr(records, "doi"),
    cited_by = col_num(records, "citedByCount"),
    source = source,
    source_id = ifelse(
      !is.na(pmid) & nzchar(pmid),
      paste0("PMID:", pmid),
      paste0(source, ":", id)
    ),
    source_url = paste0(EUROPEPMC_WEB, "/article/", source, "/", id)
  )
}

#' Read the hit count off a Europe PMC response
#'
#' Pure. `0` is a real count and comes back as `0L`; only an absent `hitCount`
#' is `NA`.
#'
#' @param body A parsed Europe PMC `search` response.
#'
#' @return A single integer, or `NA_integer_`.
#'
#' @examples
#' europepmc_parse_count(list(hitCount = 2530))
#' europepmc_parse_count(list(hitCount = 0))
#'
#' @export
europepmc_parse_count <- function(body) {
  raw <- biohttp::pluck_at(body, "hitCount", default = NA)
  if (biohttp::is_blank(raw)) {
    return(NA_integer_)
  }
  suppressWarnings(as.integer(raw))
}

#' Search Europe PMC
#'
#' @param query A query string. Build one with [europepmc_query()] to get the
#'   quoting right.
#' @param limit Maximum results to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `results` (the tibble
#'   from [europepmc_parse_results()]), `count` (Europe PMC's total, not the
#'   page size), and `query`.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(europepmc_search(europepmc_query("BRAF", "V600E")))
#' }
#'
#' @export
europepmc_search <- function(query, limit = 15, ...) {
  if (biohttp::is_blank(query)) {
    return(biohttp::status_no_data(
      source = "Europe PMC",
      detail = "no literature query was supplied"
    ))
  }
  res <- biohttp::get_json(
    EUROPEPMC_URL,
    path = "search",
    query = list(
      query = as.character(query),
      format = "json",
      resultType = "lite",
      pageSize = as.integer(limit)
    ),
    source = "Europe PMC",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- europepmc_parse_results(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Europe PMC",
      http = res$http,
      detail = paste0("Europe PMC found no publications for ", query)
    ))
  }
  biohttp::status_ok(
    data = list(
      results = parsed,
      count = europepmc_parse_count(res$data),
      query = as.character(query)
    ),
    source = "Europe PMC",
    http = res$http
  )
}

#' How many publications Europe PMC has for a query
#'
#' Asks for one row with `resultType = "idlist"`, because only the count is
#' wanted. Read the note in the file header: a count of 0 comes back as `ok`,
#' not as `no_data`.
#'
#' @inheritParams europepmc_search
#'
#' @return A biohttp envelope whose `data` is a list of `count`, `query`, and
#'   `source_url`.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(europepmc_count(europepmc_query("NF1")))$count
#' }
#'
#' @export
europepmc_count <- function(query, ...) {
  if (biohttp::is_blank(query)) {
    return(biohttp::status_no_data(
      source = "Europe PMC",
      detail = "no literature query was supplied"
    ))
  }
  res <- biohttp::get_json(
    EUROPEPMC_URL,
    path = "search",
    query = list(
      query = as.character(query),
      format = "json",
      pageSize = 1,
      resultType = "idlist"
    ),
    source = "Europe PMC",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  count <- europepmc_parse_count(res$data)
  # An absent hitCount is a malformed answer, unlike a hitCount of 0.
  if (is.na(count)) {
    return(biohttp::status_no_data(
      source = "Europe PMC",
      http = res$http,
      detail = "Europe PMC returned no hit count"
    ))
  }
  biohttp::status_ok(
    data = list(
      count = count,
      query = as.character(query),
      source_url = paste0(
        EUROPEPMC_WEB,
        "/search?query=",
        utils::URLencode(as.character(query), reserved = TRUE)
      )
    ),
    source = "Europe PMC",
    http = res$http
  )
}
