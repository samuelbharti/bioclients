# Pharos / IDG: Target Development Level.
#
# Ported from genescout/R/tools/pharos.R and a sibling app's source_pharos.R.
#
# TDL is IDG's tractability classification, ordered Tclin > Tchem > Tbio > Tdark.
# It is a druggability SIGNAL for prioritization, never a clinical call.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# Both source clients also carry a TDL -> 0-1 numeric map (Tclin = 1.0,
# Tchem = 0.75, Tbio = 0.5, Tdark = 0.25). That map is not ported, because
# turning a category into a weight is scoring, and scoring is the consuming app's
# job: a ranking model is the intellectual core of the app that holds it, and
# section 3.2 of the build plan keeps it there. This client reports the TDL that
# Pharos assigned. An app that wants a number applies its own.
#
# Endpoint: https://pharos-api.ncats.io/graphql

PHAROS_URL <- "https://pharos-api.ncats.io/graphql"
PHAROS_WEB <- "https://pharos.nih.gov/targets"

# The known TDL values, best-characterized first. Used to validate what came
# back, not to score it.
PHAROS_TDL_LEVELS <- c("Tclin", "Tchem", "Tbio", "Tdark")

# Batch-shaped already, like DGIdb's. Both source clients pass one symbol.
PHAROS_QUERY <- paste(
  "query($syms: [String!]) {",
  "  targets(targets: $syms) { targets { sym tdl } }",
  "}",
  sep = "\n"
)

#' Turn a Pharos targets response into a table
#'
#' Pure. Rows come back in `symbols` order, one per input.
#'
#' A TDL Pharos does not recognise is reported as `NA` rather than passed
#' through, so a caller never has to guess whether an unexpected string is a new
#' level or a typo.
#'
#' @param body A parsed Pharos GraphQL response body.
#' @param symbols The gene symbols that were queried, in the order asked.
#'
#' @return A tibble of `symbol`, `tdl`, and `source_url`, one row per entry in
#'   `symbols`. `tdl` is one of `Tclin`, `Tchem`, `Tbio`, `Tdark`, or `NA`.
#'
#' @examples
#' body <- list(data = list(targets = list(targets = list(
#'   list(sym = "NF1", tdl = "Tbio")
#' ))))
#' pharos_parse_targets(body, "NF1")
#'
#' @export
pharos_parse_targets <- function(body, symbols) {
  targets <- biohttp::pluck_at(body, "data", "targets", "targets")
  wanted <- toupper(trimws(as.character(symbols)))
  found <- if (is.null(targets)) {
    character()
  } else {
    toupper(col_chr(targets, "sym"))
  }

  rows <- lapply(seq_along(wanted), function(i) {
    at <- match(wanted[i], found)
    tdl <- if (is.na(at)) NA_character_ else chr_at(targets[[at]], "tdl")
    # An unrecognized level is NA. Passing an unknown string through would let a
    # caller branch on something this package has never seen.
    if (!is.na(tdl) && !(tdl %in% PHAROS_TDL_LEVELS)) {
      tdl <- NA_character_
    }
    tibble::tibble(
      symbol = wanted[i],
      tdl = tdl,
      source_url = if (is.na(at)) {
        NA_character_
      } else {
        paste0(PHAROS_WEB, "/", wanted[i])
      }
    )
  })
  do.call(rbind, rows)
}

#' Target Development Level for many genes
#'
#' One request for the whole list, because Pharos's `targets` query already takes
#' an array.
#'
#' @param symbols Gene symbols.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry in
#'   `symbols`, in the same order. See [pharos_parse_targets()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(pharos_targets(c("NF1", "EGFR")))
#' }
#'
#' @export
pharos_targets <- function(symbols, ...) {
  cleaned <- vapply(
    symbols,
    function(symbol) clean_symbol(symbol) %||% NA_character_,
    character(1),
    USE.NAMES = FALSE
  )
  usable <- toupper(cleaned[!is.na(cleaned)])
  if (length(usable) == 0) {
    return(biohttp::status_no_data(
      source = "Pharos",
      detail = "no usable gene symbols were supplied"
    ))
  }
  res <- biohttp::post_json(
    PHAROS_URL,
    body = list(query = PHAROS_QUERY, variables = list(syms = as.list(usable))),
    source = "Pharos",
    ...
  )
  bad <- biohttp::graphql_error(res, "Pharos")
  if (!is.null(bad)) {
    return(bad)
  }
  biohttp::status_ok(
    data = pharos_parse_targets(res$data, toupper(cleaned)),
    source = "Pharos",
    http = res$http
  )
}

#' Target Development Level for one gene
#'
#' A thin wrapper over [pharos_targets()].
#'
#' @param symbol A gene symbol.
#' @inheritParams pharos_targets
#'
#' @return A biohttp envelope whose `data` is a one-row tibble.
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(pharos_target("NF1"))
#' }
#'
#' @export
pharos_target <- function(symbol, ...) {
  pharos_targets(symbol, ...)
}
