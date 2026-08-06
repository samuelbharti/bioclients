# ProtVar (EBI): protein-level context at a residue.
#
# Ported from variant-reviewer/R/api_protvar.R.
#
# Two per-position endpoints, both keyed by UniProt accession plus protein
# position: /function/{acc}/{pos} and /population/{acc}/{pos}. They are separate
# entry points here rather than one combined call, because a caller often wants
# only one of them and the combined version made two requests either way.
#
# Endpoint: https://www.ebi.ac.uk/ProtVar/api

PROTVAR_URL <- "https://www.ebi.ac.uk/ProtVar/api"

#' Pull a residue position out of a protein-change string
#'
#' Pure. Accepts the several forms these arrive in: `"R175H"`, `"p.Arg175His"`,
#' `"Arg175His"`, or a bare `"175"`.
#'
#' @param variant A protein-change string.
#'
#' @return An integer position, or `NULL` when there is no number in it.
#'
#' @examples
#' protvar_position("p.Arg175His")
#' protvar_position("R175H")
#' protvar_position("no digits here")
#'
#' @export
protvar_position <- function(variant) {
  if (biohttp::is_blank(variant)) {
    return(NULL)
  }
  found <- regmatches(
    as.character(variant),
    regexpr("\\d+", as.character(variant))
  )
  if (length(found) == 0 || !nzchar(found)) {
    return(NULL)
  }
  as.integer(found)
}

#' Strip inline citations out of a UniProt function comment
#'
#' Pure.
#'
#' UniProt FUNCTION comments carry their evidence inline, as
#' `(PubMed:11025664, PubMed:12524540, ...)`. For a well-studied protein the
#' citations run longer than the prose they support. This drops the citation
#' groups and tidies the punctuation left behind, while keeping non-citation
#' parentheticals such as `(By similarity)`.
#'
#' @param text A function comment.
#'
#' @return The text with citation groups removed.
#'
#' @examples
#' protvar_strip_citations("Induces arrest (PubMed:11025664, PubMed:12524540).")
#' protvar_strip_citations("Binds DNA (By similarity).")
#'
#' @export
protvar_strip_citations <- function(text) {
  if (biohttp::is_blank(text)) {
    return(text)
  }
  ref <- "(?:PubMed:\\d+|Ref\\.\\s*\\d+|ECO:[0-9|.A-Za-z:-]+)"
  out <- as.character(text)
  # Parentheticals that are nothing but citations.
  out <- gsub(
    sprintf("\\s*\\(%s(?:\\s*[,;]\\s*%s)*\\)", ref, ref),
    "",
    out,
    perl = TRUE
  )
  # Citations mixed into a parenthetical that also says something else.
  out <- gsub(sprintf("\\s*[,;]?\\s*%s", ref), "", out, perl = TRUE)
  # Tidy what is left: empty brackets, space before punctuation, doubled spaces.
  out <- gsub("\\(\\s*[,;]*\\s*\\)", "", out, perl = TRUE)
  out <- gsub("\\s+([.,;:])", "\\1", out, perl = TRUE)
  out <- gsub("\\s{2,}", " ", out, perl = TRUE)
  trimws(out)
}

#' Turn a ProtVar function response into its text
#'
#' Pure. Returns the first `FUNCTION` comment, with citations stripped.
#'
#' @param body A parsed ProtVar `/function` response.
#'
#' @return A single string, or `NA_character_` when there is no function comment.
#'
#' @examples
#' body <- list(comments = list(list(
#'   type = "FUNCTION",
#'   text = list(list(value = "Induces arrest (PubMed:11025664)."))
#' )))
#' protvar_parse_function(body)
#'
#' @export
protvar_parse_function <- function(body) {
  comments <- biohttp::pluck_at(body, "comments")
  if (is.null(comments)) {
    return(NA_character_)
  }
  for (comment in comments) {
    if (!identical(biohttp::pluck_at(comment, "type"), "FUNCTION")) {
      next
    }
    text <- biohttp::pluck_at(comment, "text")
    if (!is.null(text) && length(text) > 0) {
      return(protvar_strip_citations(chr_at(text[[1]], "value")))
    }
  }
  NA_character_
}

#' Turn a ProtVar population response into a table
#'
#' Pure. Known variants at the residue.
#'
#' @section Sources are deduplicated:
#' A variant carries one cross-reference per supporting record, so the same
#' source name recurs many times: the TP53 175 fixture lists `NCI-TCGA` four
#' times for a single change. Reporting them raw makes a single database look
#' like corroboration from several.
#'
#' @param body A parsed ProtVar `/population` response.
#'
#' @return A tibble of `change` and `sources`, where `sources` is a
#'   comma-separated list of distinct source names. `NULL` when there are no
#'   variants.
#'
#' @examples
#' body <- list(variants = list(list(
#'   alternativeSequence = "Cys",
#'   xrefs = list(list(name = "NCI-TCGA"), list(name = "NCI-TCGA"))
#' )))
#' protvar_parse_population(body)
#'
#' @export
protvar_parse_population <- function(body) {
  variants <- biohttp::pluck_at(body, "variants")
  if (is.null(variants) || length(variants) == 0) {
    return(NULL)
  }
  change <- col_chr(variants, "alternativeSequence")
  sources <- vapply(
    variants,
    function(variant) {
      xrefs <- biohttp::pluck_at(variant, "xrefs")
      if (is.null(xrefs) || length(xrefs) == 0) {
        return(NA_character_)
      }
      names <- unique(col_chr(xrefs, "name"))
      names <- names[!is.na(names)]
      if (length(names) == 0) NA_character_ else paste(names, collapse = ", ")
    },
    character(1)
  )
  out <- tibble::tibble(change = change, sources = sources)
  out <- out[!is.na(out$change), , drop = FALSE]
  if (nrow(out) == 0) {
    return(NULL)
  }
  out
}

# Both entry points are the same GET with a different path prefix.
protvar_get <- function(endpoint, accession, position, ...) {
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "ProtVar",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  pos <- suppressWarnings(as.integer(position))
  if (length(pos) != 1L || is.na(pos)) {
    return(biohttp::status_no_data(
      source = "ProtVar",
      detail = "no usable protein position was supplied"
    ))
  }
  biohttp::get_json(
    PROTVAR_URL,
    path = paste0(endpoint, "/", cleaned, "/", pos),
    source = "ProtVar",
    ...
  )
}

#' Functional context for a residue
#'
#' @param accession A UniProt accession, for example `"P04637"`.
#' @param position A protein position. Use [protvar_position()] to get one out
#'   of a protein-change string.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a single string. See
#'   [protvar_parse_function()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(protvar_function("P04637", 175))
#' }
#'
#' @export
protvar_function <- function(accession, position, ...) {
  res <- protvar_get("function", accession, position, ...)
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- protvar_parse_function(res$data)
  if (is.na(parsed)) {
    return(biohttp::status_no_data(
      source = "ProtVar",
      http = res$http,
      detail = "ProtVar returned no function comment for this residue"
    ))
  }
  biohttp::status_ok(data = parsed, source = "ProtVar", http = res$http)
}

#' Known variants at a residue
#'
#' @inheritParams protvar_function
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [protvar_parse_population()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(protvar_population("P04637", 175))
#' }
#'
#' @export
protvar_population <- function(accession, position, ...) {
  res <- protvar_get("population", accession, position, ...)
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- protvar_parse_population(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "ProtVar",
      http = res$http,
      detail = "ProtVar lists no variants at this residue"
    ))
  }
  biohttp::status_ok(data = parsed, source = "ProtVar", http = res$http)
}
