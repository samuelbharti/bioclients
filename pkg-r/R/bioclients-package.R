#' @keywords internal
#'
#' @details
#' Every service module ships two halves, and the split is the point.
#'
#' The **client** builds a request, calls into `biohttp`, and returns its
#' envelope. It knows URLs, parameters, and rate limits, and it touches the
#' network. Start at [mygene_gene()], [gnomad_constraint()],
#' [gnomad_frequency()], or [clinvar_classification()].
#'
#' The **parser** takes an already-parsed response body and returns a canonical
#' structure. It is pure, it never touches the network, and it is what the
#' offline tests exercise. See [mygene_parse_hits()],
#' [gnomad_parse_constraint()], and [clinvar_parse_record()].
#'
#' A caller that already has a response body can use a parser on its own.
#'
#' @section What comes back:
#' Every client returns a `biohttp` envelope rather than raising, so a caller
#' branches on `res$status` and writes no `tryCatch()` of its own. Reach for
#' `biohttp::body_or_null()` when the reason for a failure is genuinely not
#' actionable.
#'
#' Reaching a source that has nothing for a query is `no_data`, not an error.
#' It is an answer.
#'
#' @section Batching:
#' [mygene_genes()] and [gnomad_constraints()] ask about many genes in one or a
#' few requests rather than one request per gene. The round trip is where
#' essentially all the time in these calls goes, so this is the difference that
#' matters. Both return rows in input order, one row per input, so a caller can
#' zip results onto its inputs by position.
"_PACKAGE"
