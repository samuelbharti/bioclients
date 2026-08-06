# Genomics England PanelApp: curated diagnostic gene panels.
#
# Ported from genescout/R/tools/panelapp.R.
#
# THE PANEL INDEX HAS NO WORKING SEARCH.
#
# `panels/` accepts a `search` parameter and answers HTTP 200 to it, but the
# result is the unfiltered index. Sending a disease name there and reading the
# first result gives whichever panel sorts first overall, for every query. That
# is why [panelapp_all_panels()] exists at all: the index has to be walked and
# matched by the caller.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# Two things, and both are the same kind of thing.
#
# genescout picks the panel matching a disease name, using a stop-word list and
# a rule requiring a strict majority of the disease's tokens to match. Which
# panel answers a clinical question is the review's judgement.
#
# genescout also maps the confidence level to a 0 to 1 weight, green to 1.0 and
# amber to 0.5, and drops red. This client returns the level PanelApp publishes,
# for every gene including red. Turning a curation rating into a weight, or
# deciding red is not worth seeing, is scoring, and it belongs to whatever is
# doing the scoring.
#
# Endpoint: https://panelapp.genomicsengland.co.uk/api/v1

PANELAPP_URL <- "https://panelapp.genomicsengland.co.uk/api/v1"
PANELAPP_WEB <- "https://panelapp.genomicsengland.co.uk/panels"

# PanelApp's traffic-light rating. The API sends the number; the colour is what
# the curators and the web interface actually use.
PANELAPP_LEVELS <- c("3" = "green", "2" = "amber", "1" = "red")

#' Turn a PanelApp panel index page into a table
#'
#' Pure.
#'
#' @param body A parsed PanelApp `panels/` response.
#'
#' @return A tibble of `id`, `name`, `version`, `disorders` (a list column, as a
#'   panel carries any number), and `source_url`, one row per panel. `NULL` when
#'   the page is empty.
#'
#' @inherit panelapp_panels references
#'
#' @examples
#' body <- list(results = list(
#'   list(id = 255, name = "Neurofibromatosis Type 1",
#'        relevant_disorders = list("NF1"))
#' ))
#' panelapp_parse_index(body)
#'
#' @export
panelapp_parse_index <- function(body) {
  # An index page nests its panels under `results`; a bare list of panels is
  # accepted too, so a caller holding an already-flattened list can reuse this.
  records <- biohttp::pluck_at(body, "results") %||% body
  if (is.null(records) || length(records) == 0) {
    return(NULL)
  }
  id <- col_chr(records, "id")
  keep <- !is.na(id) & nzchar(id)
  if (!any(keep)) {
    return(NULL)
  }
  records <- records[keep]
  tibble::tibble(
    id = id[keep],
    name = col_chr(records, "name"),
    version = col_chr(records, "version"),
    disorders = lapply(records, function(rec) {
      as.character(unlist(
        biohttp::pluck_at(rec, "relevant_disorders", default = list()),
        use.names = FALSE
      ))
    }),
    source_url = paste0(PANELAPP_WEB, "/", id[keep], "/")
  )
}

#' Turn a PanelApp panel detail into a table of genes
#'
#' Pure.
#'
#' Every gene on the panel is returned, red included. See the note in the file
#' header on why the confidence level is not converted to a weight here.
#'
#' @param body A parsed PanelApp panel detail response.
#'
#' @return A tibble of `symbol`, `confidence` (the raw level PanelApp sent),
#'   `level` (`"green"`, `"amber"`, `"red"`, or `NA` for anything else),
#'   `hgnc_id`, and `source_url`. `NULL` when the panel has no genes.
#'
#' @inherit panelapp_panels references
#'
#' @examples
#' body <- list(id = 255, genes = list(
#'   list(
#'     entity_name = "NF1",
#'     gene_data = list(gene_symbol = "NF1", hgnc_id = "HGNC:7765"),
#'     confidence_level = "3"
#'   )
#' ))
#' panelapp_parse_panel(body)
#'
#' @export
panelapp_parse_panel <- function(body) {
  genes <- biohttp::pluck_at(body, "genes")
  if (is.null(genes) || length(genes) == 0) {
    return(NULL)
  }
  # `gene_data` is absent on a panel entry for a region or an STR, where
  # `entity_name` is the only name there is.
  symbol <- col_chr(genes, "gene_data", "gene_symbol")
  fallback <- col_chr(genes, "entity_name")
  symbol[is.na(symbol)] <- fallback[is.na(symbol)]
  keep <- !is.na(symbol) & nzchar(symbol)
  if (!any(keep)) {
    return(NULL)
  }
  confidence <- col_chr(genes, "confidence_level")[keep]
  panel_id <- chr_at(body, "id")
  tibble::tibble(
    symbol = symbol[keep],
    confidence = confidence,
    level = unname(PANELAPP_LEVELS[confidence]),
    hgnc_id = col_chr(genes, "gene_data", "hgnc_id")[keep],
    source_url = paste0(PANELAPP_WEB, "/", panel_id, "/")
  )
}

#' One page of the PanelApp panel index
#'
#' @param page The 1-based page number.
#' @param page_size Panels per page.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `panels` (the tibble
#'   from [panelapp_parse_index()]) and `has_more`, which is `TRUE` when
#'   PanelApp sent a `next` link.
#'
#' @references
#' Martin et al. (2019). PanelApp crowdsources expert knowledge to establish
#' consensus diagnostic gene panels. Nature Genetics 51(11), 1560-1565.
#' \doi{10.1038/s41588-019-0528-2}
#'
#' Service documentation: <https://panelapp.genomicsengland.co.uk/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(panelapp_panels())$panels
#' }
#'
#' @export
panelapp_panels <- function(page = 1, page_size = 100, ...) {
  res <- biohttp::get_json(
    PANELAPP_URL,
    path = "panels/",
    query = list(page = page, page_size = page_size),
    source = "PanelApp",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- panelapp_parse_index(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "PanelApp",
      http = res$http,
      detail = paste0("PanelApp returned no panels on page ", page)
    ))
  }
  biohttp::status_ok(
    data = list(
      panels = parsed,
      has_more = !biohttp::is_blank(
        biohttp::pluck_at(res$data, "next", default = NULL)
      )
    ),
    source = "PanelApp",
    http = res$http
  )
}

#' The whole PanelApp panel index
#'
#' Walks the paginated index and stacks the pages. Stops at `max_pages`, which
#' is a guard against an unbounded walk rather than a tuned number.
#'
#' @param max_pages Maximum pages to fetch.
#' @param page_size Panels per page.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble from
#'   [panelapp_parse_index()], stacked across pages. A page that fails after the
#'   first is where the walk stops, keeping what was already collected; a
#'   failure on the first page is returned as-is.
#'
#' @inherit panelapp_panels references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(panelapp_all_panels(max_pages = 2))
#' }
#'
#' @export
panelapp_all_panels <- function(max_pages = 6, page_size = 100, ...) {
  pages <- list()
  http <- NA_integer_
  for (page in seq_len(max_pages)) {
    res <- panelapp_panels(page = page, page_size = page_size, ...)
    if (!isTRUE(res$ok)) {
      # Nothing collected yet means the index itself is unreachable, which the
      # caller needs to see. A later failure still leaves a usable partial
      # index, and truncating it is what max_pages does anyway.
      if (page == 1) {
        return(res)
      }
      break
    }
    http <- res$http
    pages[[length(pages) + 1]] <- res$data$panels
    if (!isTRUE(res$data$has_more)) {
      break
    }
  }
  if (length(pages) == 0) {
    return(biohttp::status_no_data(
      source = "PanelApp",
      http = http,
      detail = "PanelApp returned no panels"
    ))
  }
  biohttp::status_ok(
    data = do.call(rbind, pages),
    source = "PanelApp",
    http = http
  )
}

#' The genes on one PanelApp panel
#'
#' @param panel_id A PanelApp panel id, for example `255`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [panelapp_parse_panel()].
#'
#' @inherit panelapp_panels references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(panelapp_panel(255))
#' }
#'
#' @export
panelapp_panel <- function(panel_id, ...) {
  id <- trimws(as.character(panel_id %||% ""))
  # The id goes straight into the path, so anything that is not one is refused
  # rather than turned into a request for some other resource.
  if (!grepl("^[0-9]+$", id)) {
    return(biohttp::status_no_data(
      source = "PanelApp",
      detail = "a numeric PanelApp panel id is required"
    ))
  }
  res <- biohttp::get_json(
    PANELAPP_URL,
    path = paste0("panels/", id, "/"),
    source = "PanelApp",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- panelapp_parse_panel(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "PanelApp",
      http = res$http,
      detail = paste0("PanelApp panel ", id, " has no genes")
    ))
  }
  biohttp::status_ok(data = parsed, source = "PanelApp", http = res$http)
}
