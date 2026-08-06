# Ensembl REST: VEP by variant id, and the gene model.
#
# Ported from variant-reviewer/R/api_ensembl.R.
#
# TWO ENDPOINTS THAT ARE NOT THE ONE IN vep.R.
#
# vep.R posts genomic regions to `vep/homo_sapiens/region` in batches of 200.
# This file covers `vep/human/id/<rsid>`, which takes a variant id and answers
# for one variant, and `lookup/id/<id>?expand=1`, which returns a gene's
# transcripts and exons. Use vep.R when the coordinates are in hand; use this
# when an rsID or a gene id is.
#
# EXONS ARE NUMBERED IN TRANSCRIPTION ORDER, NOT COORDINATE ORDER.
#
# On the minus strand, exon 1 is the one with the HIGHEST genomic coordinate.
# Numbering by ascending coordinate regardless of strand puts exon 1 at the
# wrong end of every minus-strand gene, which is roughly half of them, and the
# result looks perfectly plausible.
#
# Endpoint: https://rest.ensembl.org

ENSEMBL_URL <- "https://rest.ensembl.org"

#' Turn Ensembl transcript consequences into a table
#'
#' Pure.
#'
#' Every consequence is returned, with its `biotype`, rather than only the
#' protein-coding ones. Which biotypes are worth showing is the caller's call,
#' and a filter here would hide the non-coding rows from a caller that wanted
#' them.
#'
#' @param consequences The `transcript_consequences` array from a VEP record.
#'
#' @return A tibble of `gene`, `gene_id`, `transcript`, `biotype`,
#'   `consequence`, `impact`, `sift`, and `polyphen`. `NULL` when there are
#'   none.
#'
#' @examples
#' ensembl_parse_consequences(list(list(
#'   gene_symbol = "BRAF",
#'   transcript_id = "ENST00000646891",
#'   biotype = "protein_coding",
#'   consequence_terms = list("missense_variant"),
#'   impact = "MODERATE"
#' )))
#'
#' @export
ensembl_parse_consequences <- function(consequences) {
  if (is.null(consequences) || length(consequences) == 0) {
    return(NULL)
  }
  tibble::tibble(
    gene = col_chr(consequences, "gene_symbol"),
    gene_id = col_chr(consequences, "gene_id"),
    transcript = col_chr(consequences, "transcript_id"),
    biotype = col_chr(consequences, "biotype"),
    # A transcript can carry several terms; they are joined rather than reduced
    # to the first, which would drop the rest with no sign it had happened.
    consequence = vapply(
      consequences,
      function(rec) {
        terms <- biohttp::pluck_at(rec, "consequence_terms")
        if (is.null(terms) || length(terms) == 0) {
          return(NA_character_)
        }
        paste(as.character(unlist(terms, use.names = FALSE)), collapse = ", ")
      },
      character(1)
    ),
    impact = col_chr(consequences, "impact"),
    sift = col_chr(consequences, "sift_prediction"),
    polyphen = col_chr(consequences, "polyphen_prediction")
  )
}

# The VEP id endpoint answers with an array of records. Taking [[1]] blindly
# breaks on a body that is already a single record, because [[1]] then returns
# the first FIELD of it, and a string parses to a row of NA rather than to an
# error. Named fields are what tells the two apart.
ensembl_first_record <- function(body) {
  if (is.null(body) || length(body) == 0) {
    return(NULL)
  }
  if (!is.null(names(body)) && any(nzchar(names(body)))) {
    return(body)
  }
  body[[1]]
}

#' Turn an Ensembl VEP record into a result
#'
#' Pure.
#'
#' @param record One element of a VEP response.
#'
#' @return A list of `most_severe`, `assembly`, and `consequences` (the tibble
#'   from [ensembl_parse_consequences()]). `NULL` when the record is empty.
#'
#' @examples
#' ensembl_parse_vep(list(
#'   most_severe_consequence = "missense_variant",
#'   assembly_name = "GRCh38"
#' ))
#'
#' @export
ensembl_parse_vep <- function(record) {
  if (is.null(record) || length(record) == 0) {
    return(NULL)
  }
  list(
    most_severe = chr_at(record, "most_severe_consequence"),
    assembly = chr_at(record, "assembly_name"),
    consequences = ensembl_parse_consequences(
      biohttp::pluck_at(record, "transcript_consequences")
    )
  )
}

#' Turn an Ensembl gene lookup into a gene model
#'
#' Pure.
#'
#' @section Exon numbering follows the strand:
#' Exons are sorted by genomic coordinate and then numbered in transcription
#' order. On the minus strand that means the highest-coordinate exon is exon 1.
#' Numbering by coordinate alone would put exon 1 at the wrong end of every
#' minus-strand gene.
#'
#' @param record A parsed `lookup/id` response fetched with `expand=1`.
#'
#' @return A list of `transcript`, `strand`, `region`, `gene_start`,
#'   `gene_end`, and `exons`, a tibble of `start`, `end`, and `number` sorted by
#'   `start`. `NULL` when the record carries no transcript with exons.
#'
#' @examples
#' record <- list(
#'   seq_region_name = "17", start = 1, end = 900,
#'   Transcript = list(list(
#'     id = "ENST1", is_canonical = 1, strand = -1,
#'     Exon = list(list(start = 800, end = 900), list(start = 1, end = 100))
#'   ))
#' )
#' ensembl_parse_gene_model(record)$exons
#'
#' @export
ensembl_parse_gene_model <- function(record) {
  transcripts <- biohttp::pluck_at(record, "Transcript")
  if (is.null(transcripts) || length(transcripts) == 0) {
    return(NULL)
  }
  canonical <- Filter(
    function(tx) isTRUE(as.logical(biohttp::pluck_at(tx, "is_canonical"))),
    transcripts
  )
  tx <- if (length(canonical) > 0) canonical[[1]] else transcripts[[1]]
  exons <- biohttp::pluck_at(tx, "Exon")
  if (is.null(exons) || length(exons) == 0) {
    return(NULL)
  }
  start <- col_num(exons, "start")
  end <- col_num(exons, "end")
  keep <- !is.na(start) & !is.na(end)
  if (!any(keep)) {
    return(NULL)
  }
  start <- start[keep]
  end <- end[keep]
  order_by_start <- order(start)
  start <- start[order_by_start]
  end <- end[order_by_start]
  strand <- num_at(tx, "strand")
  # See the section above: transcription order, not coordinate order.
  number <- if (isTRUE(strand < 0)) {
    rev(seq_along(start))
  } else {
    seq_along(start)
  }
  list(
    transcript = chr_at(tx, "id"),
    strand = strand,
    region = chr_at(record, "seq_region_name"),
    gene_start = num_at(record, "start", default = min(start)),
    gene_end = num_at(record, "end", default = max(end)),
    exons = tibble::tibble(start = start, end = end, number = number)
  )
}

#' Run VEP for a variant id
#'
#' @section An rsID does not identify an allele:
#' The same caveat as [gnomad_frequency()] and [clinvar_classification()].
#' `7:g.140753336A>T` and `7:g.140753336A>C` share `rs113488022`, so Ensembl may
#' answer for more than one allele. The first record is used. Post coordinates
#' with [vep_variants()] when the allele matters.
#'
#' @param rsid A variant id, for example `"rs113488022"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the list described in
#'   [ensembl_parse_vep()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(ensembl_vep_id("rs113488022"))$consequences
#' }
#'
#' @export
ensembl_vep_id <- function(rsid, ...) {
  id <- trimws(as.character(rsid %||% ""))
  # The id goes into the path, so anything that is not a variant id is refused
  # rather than turned into a request for some other resource.
  if (!grepl("^[A-Za-z0-9_.:-]+$", id)) {
    return(biohttp::status_no_data(
      source = "Ensembl VEP",
      detail = "a variant id is required, for example rs113488022"
    ))
  }
  res <- biohttp::get_json(
    ENSEMBL_URL,
    # `content-type` is a query parameter here, not a header. Drop it and
    # Ensembl serves its HTML browser page with HTTP 200, so the call looks like
    # a success right up to the point the JSON parser hits the first tag.
    # Confirmed live, see test-live.R.
    path = paste0("vep/human/id/", id),
    query = list(`content-type` = "application/json"),
    source = "Ensembl VEP",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- ensembl_parse_vep(ensembl_first_record(res$data))
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Ensembl VEP",
      http = res$http,
      detail = paste0("Ensembl VEP has no record for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Ensembl VEP", http = res$http)
}

#' The exon model for a gene
#'
#' @param gene_id An Ensembl gene id, for example `"ENSG00000157764"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the list described in
#'   [ensembl_parse_gene_model()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(ensembl_gene_model("ENSG00000157764"))$exons
#' }
#'
#' @export
ensembl_gene_model <- function(gene_id, ...) {
  id <- trimws(as.character(gene_id %||% ""))
  if (!grepl("^ENS[A-Z]*[GT][0-9.]+$", id)) {
    return(biohttp::status_no_data(
      source = "Ensembl",
      detail = "an Ensembl gene id is required for the gene model"
    ))
  }
  res <- biohttp::get_json(
    ENSEMBL_URL,
    path = paste0("lookup/id/", id),
    # expand=1 is what makes Ensembl send the transcripts and their exons. The
    # response is otherwise a bare gene record and the model comes out empty.
    query = list(expand = 1, `content-type` = "application/json"),
    source = "Ensembl",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- ensembl_parse_gene_model(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Ensembl",
      http = res$http,
      detail = paste0("Ensembl returned no exon model for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Ensembl", http = res$http)
}
