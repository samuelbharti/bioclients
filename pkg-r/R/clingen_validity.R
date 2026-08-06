# ClinGen gene-disease validity: how strong the curated evidence is that a gene
# causes a disease.
#
# Ported from genescout/R/tools/clingen.R.
#
# A different service from the Allele Registry in clingen.R, on a different
# host, so the two are kept apart and biohttp gives each its own breaker.
#
# THERE IS NO PER-GENE API.
#
# ClinGen publishes this as one public CC0 bulk CSV and nothing else, so the
# client fetches the whole file and the caller filters it. biohttp caches on the
# URL, so the download happens once per process and every later lookup is served
# in memory. That is why this is one fetch plus a pure filter rather than a
# request per gene.
#
# THE HEADER IS NOT THE FIRST ROW, AND THERE ARE TWO DECORATION ROWS.
#
# The file opens with a three-line banner, then a row of "++++" separators, then
# the real header, then ANOTHER row of "++++", and only then the data. Reading
# it with any ordinary CSV reader gives a table whose column names are
# "CLINGEN GENE DISEASE VALIDITY CURATIONS" and whose first two data rows are
# rows of plus signs. So the header is found by looking for it, and the
# decoration rows are dropped.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout maps the classification to an ordinal strength, Definitive 4 through
# Limited 1, with Disputed and Refuted at 0. That is a weight, and turning a
# curation verdict into a number is scoring. The classification string is
# returned as ClinGen wrote it.
#
# ON NOT ADDING A DEPENDENCY.
#
# utils::read.csv reads this in well under the time the download takes, so
# neither data.table nor vroom earns a place in Suggests for it. If a real
# bottleneck ever shows up here, CONTRIBUTING says which to reach for.
#
# Endpoint: https://search.clinicalgenome.org/kb/gene-validity/download

CLINGEN_VALIDITY_URL <- "https://search.clinicalgenome.org"
CLINGEN_VALIDITY_PATH <- "kb/gene-validity/download"

# The marker that finds the real header row.
CLINGEN_HEADER_MARKER <- "\"GENE SYMBOL\""

#' Parse the ClinGen gene-validity CSV
#'
#' Pure. Takes the file as a single string, the way [biohttp::get_text()]
#' returns it.
#'
#' @param text The CSV file contents.
#'
#' @return A tibble of `gene`, `hgnc`, `disease`, `mondo`, `moi`,
#'   `classification`, `sop`, `panel`, `date`, and `source_url`, one row per
#'   curation. A gene appears once per disease it has been curated against.
#'   `NULL` when the file is empty or carries no header.
#'
#' @examples
#' text <- paste(
#'   '"CLINGEN GENE DISEASE VALIDITY CURATIONS","",""',
#'   '"+++","+++","+++"',
#'   '"GENE SYMBOL","DISEASE LABEL","CLASSIFICATION"',
#'   '"+++","+++","+++"',
#'   '"NF1","neurofibromatosis type 1","Definitive"',
#'   sep = "\n"
#' )
#' clingen_parse_validity(text)
#'
#' @export
clingen_parse_validity <- function(text) {
  if (is.null(text) || length(text) != 1 || !nzchar(text)) {
    return(NULL)
  }
  lines <- strsplit(text, "\r?\n")[[1]]
  header <- which(grepl(CLINGEN_HEADER_MARKER, lines, fixed = TRUE))
  if (length(header) == 0) {
    return(NULL)
  }
  body <- lines[header[1]:length(lines)]
  # Drop the separator rows and any blank line. A row of plus signs would
  # otherwise arrive as a curation.
  body <- body[!grepl("^\"\\++\"", body) & nzchar(trimws(body))]
  if (length(body) <= 1) {
    return(NULL)
  }
  frame <- tryCatch(
    utils::read.csv(
      text = paste(body, collapse = "\n"),
      check.names = FALSE,
      colClasses = "character"
    ),
    error = function(e) NULL
  )
  if (is.null(frame) || nrow(frame) == 0) {
    return(NULL)
  }
  column <- function(name) {
    if (name %in% names(frame)) {
      as.character(frame[[name]])
    } else {
      rep(NA_character_, nrow(frame))
    }
  }
  tibble::tibble(
    gene = column("GENE SYMBOL"),
    hgnc = column("GENE ID (HGNC)"),
    disease = column("DISEASE LABEL"),
    mondo = column("DISEASE ID (MONDO)"),
    moi = column("MOI"),
    classification = column("CLASSIFICATION"),
    sop = column("SOP"),
    panel = column("GCEP"),
    date = column("CLASSIFICATION DATE"),
    source_url = column("ONLINE REPORT")
  )
}

#' Filter a parsed validity table to one or more genes
#'
#' Pure. Matching is case-insensitive on the symbol.
#'
#' @param table A tibble from [clingen_parse_validity()].
#' @param symbols Gene symbols.
#'
#' @return The matching rows, or `NULL` when none match.
#'
#' @examples
#' table <- tibble::tibble(
#'   gene = c("NF1", "TP53"),
#'   classification = c("Definitive", "Definitive")
#' )
#' clingen_validity_for(table, "nf1")
#'
#' @export
clingen_validity_for <- function(table, symbols) {
  if (is.null(table) || nrow(table) == 0) {
    return(NULL)
  }
  wanted <- toupper(trimws(as.character(symbols)))
  wanted <- wanted[!is.na(wanted) & nzchar(wanted)]
  if (length(wanted) == 0) {
    return(NULL)
  }
  hit <- toupper(table$gene) %in% wanted
  if (!any(hit)) {
    return(NULL)
  }
  table[hit, , drop = FALSE]
}

#' The ClinGen gene-disease validity table
#'
#' Fetches the whole published CSV. See the note in the file header on why this
#' is one download rather than a request per gene, and use
#' [clingen_validity_for()] to narrow the result.
#'
#' @param ... Passed to [biohttp::get_text()], for example `timeout`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [clingen_parse_validity()].
#'
#' @examples
#' \donttest{
#' table <- biohttp::body_or_null(clingen_gene_validity())
#' clingen_validity_for(table, c("NF1", "TP53"))
#' }
#'
#' @export
clingen_gene_validity <- function(...) {
  res <- biohttp::get_text(
    CLINGEN_VALIDITY_URL,
    path = CLINGEN_VALIDITY_PATH,
    source = "ClinGen",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- clingen_parse_validity(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "ClinGen",
      http = res$http,
      detail = "the ClinGen validity file carried no curations"
    ))
  }
  biohttp::status_ok(data = parsed, source = "ClinGen", http = res$http)
}
