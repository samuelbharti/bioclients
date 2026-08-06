# IMPC: the phenotypes a whole-gene knockout produces in mouse.
#
# Ported from genescout/R/tools/impc.R.
#
# Two GETs against the public IMPC Solr, because they are two different cores
# and the second cannot be reached from a human symbol.
#
# THE PHENOTYPE CORE HAS NO HUMAN SYMBOL FIELD.
#
# The `genotype-phenotype` core is keyed by MGI marker accession and carries no
# `human_gene_symbol` at all. Solr does not reject a query against a field it
# does not have, so asking it for `human_gene_symbol:NF1` returns HTTP 200 with
# zero documents. That reads exactly like "IMPC found no phenotype for this
# gene", which is why [impc_gene_phenotypes()] always resolves the ortholog
# through the `gene` core first.
#
# WHAT A MISS MEANS HERE.
#
# A gene IMPC never phenotyped and a gene phenotyped with no significant
# abnormality both come back with zero associations, and the response does not
# distinguish them. Both are `no_data`, never a count of 0, so an untested gene
# is never mistaken for a tested one that came back clean.
#
# Endpoint: https://www.ebi.ac.uk/mi/impc/solr

IMPC_URL <- "https://www.ebi.ac.uk/mi/impc/solr"
IMPC_WEB <- "https://www.mousephenotype.org/data/genes"

# IMPC reports one document per phenotype, sex, zygosity, and parameter, so a
# gene with 4 significant terms can arrive as several hundred rows.
IMPC_MAX_ROWS <- 500L

#' Read the mouse ortholog out of an IMPC gene-core response
#'
#' Pure.
#'
#' @param body A parsed IMPC `gene/select` response.
#'
#' @return A one-row tibble of `mgi` and `marker_symbol`. `NULL` when the gene
#'   has no IMPC mouse ortholog, which includes a document carrying no usable
#'   accession.
#'
#' @inherit impc_mouse_ortholog references
#'
#' @examples
#' body <- list(response = list(docs = list(
#'   list(mgi_accession_id = "MGI:97306", marker_symbol = "Nf1")
#' )))
#' impc_parse_ortholog(body)
#'
#' @export
impc_parse_ortholog <- function(body) {
  docs <- biohttp::pluck_at(body, "response", "docs")
  if (is.null(docs) || length(docs) == 0) {
    return(NULL)
  }
  mgi <- chr_at(docs[[1]], "mgi_accession_id")
  # Everything downstream interpolates this into a Solr query, so a value that
  # is not an accession is refused here rather than sent.
  if (is.na(mgi) || !grepl("^MGI:[0-9]+$", mgi)) {
    return(NULL)
  }
  tibble::tibble(
    mgi = mgi,
    marker_symbol = chr_at(docs[[1]], "marker_symbol")
  )
}

#' Turn an IMPC phenotype response into a table
#'
#' Pure.
#'
#' @section One row per term, not per observation:
#' IMPC reports a document per phenotype, sex, zygosity, and parameter
#' combination, so the same term recurs many times over. Rows are collapsed to
#' distinct `mp_id`, keeping the first occurrence. Counting the raw documents
#' would report a gene's phenotype breadth several times over.
#'
#' `mp_id` is an MP term for most phenotypes and an MPATH term for pathology
#' findings. Both appear in the same field.
#'
#' @param body A parsed IMPC `genotype-phenotype/select` response.
#' @param mgi The MGI marker accession that was queried.
#'
#' @return A tibble of `mp_id`, `mp_name`, `allele`, `allele_symbol`,
#'   `zygosity`, and `source_url`, one row per distinct phenotype term. `NULL`
#'   when there are none.
#'
#' @inherit impc_mouse_ortholog references
#'
#' @examples
#' body <- list(response = list(docs = list(
#'   list(
#'     mp_term_id = "MP:0011100",
#'     mp_term_name = "preweaning lethality, complete penetrance",
#'     allele_accession_id = "MGI:4364806",
#'     zygosity = "homozygote"
#'   )
#' )))
#' impc_parse_phenotypes(body, "MGI:97306")
#'
#' @export
impc_parse_phenotypes <- function(body, mgi = NA_character_) {
  docs <- biohttp::pluck_at(body, "response", "docs")
  if (is.null(docs) || length(docs) == 0) {
    return(NULL)
  }
  mp_id <- col_chr(docs, "mp_term_id")
  keep <- !is.na(mp_id) & nzchar(mp_id) & !duplicated(mp_id)
  if (!any(keep)) {
    return(NULL)
  }
  tibble::tibble(
    mp_id = mp_id[keep],
    mp_name = col_chr(docs, "mp_term_name")[keep],
    allele = col_chr(docs, "allele_accession_id")[keep],
    allele_symbol = col_chr(docs, "allele_symbol")[keep],
    zygosity = col_chr(docs, "zygosity")[keep],
    source_url = paste0(IMPC_WEB, "/", mgi)
  )
}

#' The mouse ortholog IMPC holds for a human gene
#'
#' @param symbol A human gene symbol, for example `"NF1"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the one-row tibble described in
#'   [impc_parse_ortholog()].
#'
#' @references
#' Groza et al. (2023). The International Mouse Phenotyping Consortium:
#' comprehensive knockout phenotyping underpinning the study of human
#' disease. Nucleic Acids Research 51(D1), D1038-D1045.
#' \doi{10.1093/nar/gkac972}
#'
#' Service documentation: <https://www.mousephenotype.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(impc_mouse_ortholog("NF1"))
#' }
#'
#' @export
impc_mouse_ortholog <- function(symbol, ...) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "IMPC",
      detail = "no usable gene symbol was supplied"
    ))
  }
  res <- biohttp::get_json(
    IMPC_URL,
    path = "gene/select",
    # The symbol goes into a Solr query expression, so it is stripped to
    # identifier characters first. A colon or a quote in it would otherwise
    # change which field is being searched.
    query = list(
      q = paste0("human_gene_symbol:", toupper(cleaned)),
      fl = "mgi_accession_id,marker_symbol",
      rows = 1,
      wt = "json"
    ),
    source = "IMPC",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- impc_parse_ortholog(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "IMPC",
      http = res$http,
      detail = paste0("IMPC has no mouse ortholog for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "IMPC", http = res$http)
}

#' Significant knockout phenotypes IMPC records for a human gene
#'
#' Resolves the mouse ortholog first, then queries the phenotype core by its MGI
#' accession. See the note in the file header on why the second query cannot be
#' made from the human symbol.
#'
#' @param symbol A human gene symbol, for example `"NF1"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `mgi`, `marker_symbol`,
#'   `phenotypes` (the tibble from [impc_parse_phenotypes()]), and `source_url`.
#'
#' @inherit impc_mouse_ortholog references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(impc_gene_phenotypes("NF1"))$phenotypes
#' }
#'
#' @export
impc_gene_phenotypes <- function(symbol, ...) {
  ortholog <- impc_mouse_ortholog(symbol, ...)
  if (!isTRUE(ortholog$ok)) {
    return(ortholog)
  }
  mgi <- ortholog$data$mgi[[1]]
  res <- biohttp::get_json(
    IMPC_URL,
    path = "genotype-phenotype/select",
    query = list(
      q = paste0("marker_accession_id:\"", mgi, "\""),
      fl = paste(
        "mp_term_id",
        "mp_term_name",
        "allele_accession_id",
        "allele_symbol",
        "zygosity",
        sep = ","
      ),
      rows = IMPC_MAX_ROWS,
      wt = "json"
    ),
    source = "IMPC",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- impc_parse_phenotypes(res$data, mgi)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "IMPC",
      http = res$http,
      detail = paste0("IMPC records no significant phenotype for ", mgi)
    ))
  }
  biohttp::status_ok(
    data = list(
      mgi = mgi,
      marker_symbol = ortholog$data$marker_symbol[[1]],
      phenotypes = parsed,
      source_url = paste0(IMPC_WEB, "/", mgi)
    ),
    source = "IMPC",
    http = res$http
  )
}
