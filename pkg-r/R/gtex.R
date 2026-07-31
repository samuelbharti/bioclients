# GTEx Portal: median gene expression across tissues.
#
# Ported from variant-reviewer/R/api_gtex.R and genescout/R/tools/gtex.R.
#
# Endpoint: https://gtexportal.org/api/v2

GTEX_URL <- "https://gtexportal.org/api/v2"
GTEX_WEB <- "https://gtexportal.org/home/gene"

# GTEx v8 is what the versioned GENCODE ids in the reference endpoint belong to.
# Changing this without changing where the ids come from returns nothing.
GTEX_DATASET <- "gtex_v8"

#' Turn a GTEx gene-reference response into a table
#'
#' Pure.
#'
#' @param body A parsed GTEx `reference/gene` response.
#'
#' @return A tibble of `symbol`, `gencode_id`, `entrez`, `chromosome`, and
#'   `gene_type`. `NULL` when the gene is not in the reference.
#'
#' @examples
#' body <- list(data = list(list(
#'   geneSymbol = "TP53", gencodeId = "ENSG00000141510.16", entrezGeneId = 7157
#' )))
#' gtex_parse_reference(body)
#'
#' @export
gtex_parse_reference <- function(body) {
  rows <- biohttp::pluck_at(body, "data")
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  tibble::tibble(
    symbol = col_chr(rows, "geneSymbol"),
    gencode_id = col_chr(rows, "gencodeId"),
    entrez = col_chr(rows, "entrezGeneId"),
    chromosome = col_chr(rows, "chromosome"),
    gene_type = col_chr(rows, "geneType")
  )
}

#' Turn a GTEx median-expression response into a table
#'
#' Pure. Tissue ids arrive underscore-separated (`Adipose_Subcutaneous`) and are
#' returned both raw and as a readable label.
#'
#' @param body A parsed GTEx `expression/medianGeneExpression` response.
#'
#' @return A tibble of `tissue_id`, `tissue`, `median_tpm`, and `gencode_id`.
#'   Rows with no median are dropped. `NULL` when there is nothing usable.
#'
#' @examples
#' body <- list(data = list(list(
#'   median = 22.459, tissueSiteDetailId = "Adipose_Subcutaneous",
#'   gencodeId = "ENSG00000141510.16", unit = "TPM"
#' )))
#' gtex_parse_expression(body)
#'
#' @export
gtex_parse_expression <- function(body) {
  rows <- biohttp::pluck_at(body, "data")
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  tissue_id <- col_chr(rows, "tissueSiteDetailId")
  out <- tibble::tibble(
    tissue_id = tissue_id,
    tissue = gsub("_", " ", tissue_id),
    median_tpm = col_num(rows, "median"),
    gencode_id = col_chr(rows, "gencodeId")
  )
  out <- out[!is.na(out$median_tpm), , drop = FALSE]
  if (nrow(out) == 0) {
    return(NULL)
  }
  out
}

#' Resolve a gene to GTEx's versioned GENCODE id
#'
#' @section Why this call exists:
#' GTEx expression queries need a **versioned** GENCODE id such as
#' `ENSG00000141510.16`. MyGene and most other sources return the unversioned
#' `ENSG00000141510`, and GTEx returns nothing for it. There is no way to derive
#' the version, so it has to be looked up here first.
#'
#' This is the reason a GTEx expression lookup costs two requests.
#'
#' @param gene A gene symbol or an Ensembl gene id.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [gtex_parse_reference()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(gtex_gene_reference("TP53"))$gencode_id
#' }
#'
#' @export
gtex_gene_reference <- function(gene, ...) {
  cleaned <- clean_symbol(gene)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "GTEx",
      detail = "no usable gene identifier was supplied"
    ))
  }
  res <- biohttp::get_json(
    GTEX_URL,
    path = "reference/gene",
    query = list(geneId = cleaned),
    source = "GTEx",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- gtex_parse_reference(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "GTEx",
      http = res$http,
      detail = paste0("GTEx has no reference entry for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "GTEx", http = res$http)
}

#' Median expression across tissues
#'
#' Two requests: one to resolve the versioned GENCODE id, one for the expression.
#' See [gtex_gene_reference()] for why the first is unavoidable.
#'
#' Pass `gencode_id` directly to skip the lookup when you already have a
#' versioned id.
#'
#' @param gene A gene symbol or an Ensembl gene id. Ignored when `gencode_id` is
#'   given.
#' @param gencode_id A versioned GENCODE id, for example `"ENSG00000141510.16"`.
#' @param dataset The GTEx dataset. Must match the reference the id came from.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [gtex_parse_expression()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(gtex_median_expression("TP53"))
#' }
#'
#' @export
gtex_median_expression <- function(
  gene = NULL,
  gencode_id = NULL,
  dataset = GTEX_DATASET,
  ...
) {
  if (biohttp::is_blank(gencode_id)) {
    reference <- gtex_gene_reference(gene, ...)
    if (!isTRUE(reference$ok)) {
      return(reference)
    }
    gencode_id <- reference$data$gencode_id[1]
  }
  if (biohttp::is_blank(gencode_id)) {
    return(biohttp::status_no_data(
      source = "GTEx",
      detail = "could not resolve a versioned GENCODE id"
    ))
  }
  res <- biohttp::get_json(
    GTEX_URL,
    path = "expression/medianGeneExpression",
    query = list(
      gencodeId = gencode_id,
      # Both of these are required. Omitting datasetId returns nothing, and so
      # does an unversioned gencodeId.
      datasetId = dataset,
      itemsPerPage = 100
    ),
    source = "GTEx",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- gtex_parse_expression(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "GTEx",
      http = res$http,
      detail = paste0("GTEx has no expression data for ", gencode_id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "GTEx", http = res$http)
}
