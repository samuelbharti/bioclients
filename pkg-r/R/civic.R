# CIViC: expert-curated clinical evidence for cancer variants.
#
# Ported from genescout/R/tools/civic.R.
#
# The counts are a measure of how much variant-interpretation literature the
# community has curated for a gene. That is RESEARCH evidence, never a clinical
# call about a patient. CIViC is CC0.
#
# Endpoint: https://civicdb.org/api/graphql

CIVIC_URL <- "https://civicdb.org/api/graphql"
CIVIC_WEB <- "https://civicdb.org"

#' Turn a CIViC gene response into a table
#'
#' Pure.
#'
#' @param body A parsed CIViC GraphQL response body.
#'
#' @return A one-row tibble of `id`, `symbol`, `entrez_id`, `evidence_items`,
#'   `assertions`, `variants`, and `source_url`. `NULL` when CIViC does not track
#'   the gene, which it reports as `data.gene = null`.
#'
#' @inherit civic_gene references
#'
#' @examples
#' body <- list(data = list(gene = list(
#'   id = 3867, name = "NF1", entrezId = 4763, link = "/features/3867",
#'   stats = list(evidenceItemCount = 56, assertionCount = 0, variantCount = 37)
#' )))
#' civic_parse_gene(body)
#'
#' @export
civic_parse_gene <- function(body) {
  gene <- biohttp::pluck_at(body, "data", "gene")
  if (is.null(gene)) {
    return(NULL)
  }
  link <- as.character(biohttp::pluck_at(gene, "link", default = ""))
  tibble::tibble(
    id = as.character(biohttp::pluck_at(gene, "id", default = NA_character_)),
    symbol = chr_at(gene, "name"),
    entrez_id = as.character(
      biohttp::pluck_at(gene, "entrezId", default = NA_character_)
    ),
    # A gene CIViC tracks but has curated nothing for is a real zero, so these
    # default to 0 rather than NA. Absence of the gene entirely is the NULL above.
    evidence_items = as.numeric(
      biohttp::pluck_at(gene, "stats", "evidenceItemCount", default = 0)
    ),
    assertions = as.numeric(
      biohttp::pluck_at(gene, "stats", "assertionCount", default = 0)
    ),
    variants = as.numeric(
      biohttp::pluck_at(gene, "stats", "variantCount", default = 0)
    ),
    source_url = if (nzchar(link)) paste0(CIVIC_WEB, link) else CIVIC_WEB
  )
}

#' Curated clinical evidence counts for a gene
#'
#' @section On the interpolated query:
#' CIViC's `gene` query takes the symbol as an inline string rather than a
#' GraphQL variable, so the symbol is substituted into the query text. It is run
#' through this package's internal `clean_symbol()` first, which strips
#' everything outside `[A-Za-z0-9._-]`. That is what stops a crafted symbol
#' closing the string and appending its own query. Do not remove it, and do not
#' switch to `paste0()` without it.
#'
#' @param symbol A HUGO gene symbol.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [civic_parse_gene()].
#'
#' @references
#' Griffith et al. (2017). CIViC is a community knowledgebase for expert
#' crowdsourcing the clinical interpretation of variants in cancer.
#' Nature Genetics 49(2), 170-174. \doi{10.1038/ng.3774}
#'
#' Service documentation: <https://civicdb.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(civic_gene("NF1"))
#' }
#'
#' @export
civic_gene <- function(symbol, ...) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "CIViC",
      detail = "no usable gene symbol was supplied"
    ))
  }
  cleaned <- toupper(cleaned)
  query <- sprintf(
    paste0(
      "{ gene(entrezSymbol: \"%s\") { id name entrezId link ",
      "stats { evidenceItemCount assertionCount variantCount } } }"
    ),
    cleaned
  )
  res <- biohttp::post_json(
    CIVIC_URL,
    body = list(query = query),
    source = "CIViC",
    ...
  )
  bad <- biohttp::graphql_error(res, "CIViC")
  if (!is.null(bad)) {
    return(bad)
  }
  parsed <- civic_parse_gene(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "CIViC",
      http = res$http,
      detail = paste0("CIViC does not track ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "CIViC", http = res$http)
}
