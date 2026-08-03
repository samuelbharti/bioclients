# MyVariant.info: dbNSFP predictors and ClinVar significance.
#
# Ported from a sibling app's api_myvariant.R, which is the better of the two
# copies. variant-reviewer/R/api_myvariant.R is the other.
#
# THIS FILE HAS TWO TRAPS THAT FAIL SILENTLY. Both are in the tests.
#
# 1. assembly = "hg38" is load-bearing. Omit it and MyVariant answers HTTP 200
#    with {notfound: true} for every GRCh38 variant. Nothing errors, nothing
#    warns, the variants simply vanish.
#
# 2. dbNSFP exposes CADD's rankscore as `raw_rankscore`. There is no plain
#    `cadd.rankscore` key, so a reader that follows the pattern of every other
#    predictor gets NA for CADD on every variant, forever, with no error.
#
# AND ONE THING NOT TO DO.
#
# Do not read gnomAD frequencies from here. MyVariant's snapshot is gnomAD 2.1.1
# from 2019. Query gnomAD directly with gnomad_frequency().
#
# Endpoint: https://myvariant.info/v1

MYVARIANT_URL <- "https://myvariant.info/v1"

# The verified batch limit, recorded in a sibling app's source registry as a
# hard limit rather than tuning.
MYVARIANT_BATCH <- 1000L

MYVARIANT_FIELDS <- paste(
  "dbnsfp.genename",
  "dbnsfp.revel",
  "dbnsfp.cadd",
  "dbnsfp.clinpred",
  "dbnsfp.alphamissense",
  "clinvar.rcv.clinical_significance",
  sep = ","
)

#' Build a MyVariant identifier from variant components
#'
#' Pure. MyVariant indexes variants by an HGVS `g.` string, so a caller holding
#' chromosome, position, ref and alt needs this to ask about one.
#'
#' @section This is not variant identity:
#' This formats a variant the way one service wants it written. It is not a
#' canonical variant key and it does not normalize anything. The canonical key
#' is `(assembly, chromosome, position, ref, alt)`, and turning arbitrary input
#' into one is `vcfcanon`'s job, not this package's.
#'
#' @param chrom A chromosome, with or without a `chr` prefix.
#' @param pos A 1-based position.
#' @param ref,alt Reference and alternate alleles.
#'
#' @return A single string.
#'
#' @examples
#' myvariant_id("17", 7676154, "G", "C")
#' myvariant_id("chr1", 100, "AT", "A")
#'
#' @export
myvariant_id <- function(chrom, pos, ref, alt) {
  chrom <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  pos <- as.integer(pos)
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  if (nchar(ref) == 1L && nchar(alt) == 1L) {
    sprintf("chr%s:g.%d%s>%s", chrom, pos, ref, alt)
  } else if (nchar(alt) == 1L && nchar(ref) > 1L) {
    sprintf("chr%s:g.%d_%ddel", chrom, pos + 1L, pos + nchar(ref) - 1L)
  } else if (nchar(ref) == 1L && nchar(alt) > 1L) {
    sprintf(
      "chr%s:g.%d_%dins%s",
      chrom,
      pos,
      pos + 1L,
      substr(alt, 2L, nchar(alt))
    )
  } else {
    # A delins falls back to the substitution form. If MyVariant does not index
    # it, the answer is no_data, which is honest.
    sprintf("chr%s:g.%d%s>%s", chrom, pos, ref, alt)
  }
}

# ClinVar significance arrives as a string, or as a list of them across
# submissions. Collapse to one representative.
myvariant_clinvar_sig <- function(record) {
  rcv <- biohttp::pluck_at(record, "clinvar", "rcv")
  if (is.null(rcv)) {
    return(NA_character_)
  }
  sigs <- if (!is.null(rcv$clinical_significance)) {
    rcv$clinical_significance
  } else {
    unlist(lapply(rcv, function(entry) entry$clinical_significance))
  }
  sigs <- unlist(sigs)
  if (length(sigs) == 0) NA_character_ else as.character(sigs)[[1]]
}

# dbNSFP's genename is one entry per transcript, usually the same symbol
# repeated. Reduce to one representative.
myvariant_gene <- function(record) {
  gene <- biohttp::pluck_at(record, "dbnsfp", "genename")
  if (is.null(gene) || length(gene) == 0) {
    return(NA_character_)
  }
  as.character(unlist(gene, use.names = FALSE))[[1]]
}

#' Turn one MyVariant record into a table row
#'
#' Pure.
#'
#' @section CADD's key is not what you would guess:
#' Every other dbNSFP predictor exposes `<name>.rankscore`. CADD does not: it is
#' `cadd.raw_rankscore`, and there is no plain `cadd.rankscore`. Following the
#' pattern gives `NA` for CADD on every variant with nothing to indicate it.
#'
#' @param record A parsed MyVariant record.
#'
#' @return A one-row tibble of `gene`, `revel`, `cadd`, `clinpred`,
#'   `alphamissense`, and `clinvar_sig`.
#'
#' @examples
#' record <- list(dbnsfp = list(
#'   genename = list("TP53"),
#'   cadd = list(raw_rankscore = 0.17018),
#'   revel = list(rankscore = 0.4)
#' ))
#' myvariant_parse_record(record)
#'
#' @export
myvariant_parse_record <- function(record) {
  tibble::tibble(
    gene = myvariant_gene(record),
    revel = num_at(biohttp::pluck_at(record, "dbnsfp", "revel"), "rankscore"),
    # raw_rankscore, NOT rankscore. See the section above.
    cadd = num_at(
      biohttp::pluck_at(record, "dbnsfp", "cadd"),
      "raw_rankscore"
    ),
    clinpred = num_at(
      biohttp::pluck_at(record, "dbnsfp", "clinpred"),
      "rankscore"
    ),
    alphamissense = num_at(
      biohttp::pluck_at(record, "dbnsfp", "alphamissense"),
      "rankscore"
    ),
    clinvar_sig = myvariant_clinvar_sig(record)
  )
}

#' Turn a MyVariant batch response into a table
#'
#' Pure. One row per requested id, in the order asked.
#'
#' @section notfound is an answer, not an error:
#' MyVariant reports an unknown variant as HTTP 200 with `notfound: true`. That
#' is the source saying it has nothing, which is different from the call
#' failing. Such rows come back as `NA` rather than being dropped, so a caller
#' zipping by position stays aligned.
#'
#' @param body A parsed MyVariant batch response, a flat array.
#' @param ids The MyVariant ids that were requested, in order.
#'
#' @return A tibble with one row per entry in `ids`, plus an `id` column.
#'
#' @examples
#' body <- list(list(query = "chr17:g.7676154G>C", dbnsfp = list(
#'   genename = list("TP53"), cadd = list(raw_rankscore = 0.17)
#' )))
#' myvariant_parse_batch(body, "chr17:g.7676154G>C")
#'
#' @export
myvariant_parse_batch <- function(body, ids) {
  by_query <- list()
  for (record in body) {
    query <- biohttp::pluck_at(record, "query")
    if (!biohttp::is_blank(query)) {
      by_query[[query]] <- record
    }
  }
  rows <- lapply(ids, function(id) {
    record <- by_query[[id]]
    row <- if (
      is.null(record) ||
        isTRUE(biohttp::pluck_at(record, "notfound", default = FALSE))
    ) {
      myvariant_empty_row()
    } else {
      myvariant_parse_record(record)
    }
    cbind(tibble::tibble(id = as.character(id)), row)
  })
  tibble::as_tibble(do.call(rbind, rows))
}

myvariant_empty_row <- function() {
  tibble::tibble(
    gene = NA_character_,
    revel = NA_real_,
    cadd = NA_real_,
    clinpred = NA_real_,
    alphamissense = NA_real_,
    clinvar_sig = NA_character_
  )
}

#' Annotate many variants in one request
#'
#' @section assembly is not optional:
#' `assembly = "hg38"` is always sent. Without it MyVariant answers 200 with
#' `notfound` for every GRCh38 variant, so a whole cohort disappears with no
#' error anywhere. There is no way to omit it through this function, and that is
#' deliberate.
#'
#' @param ids MyVariant ids. Build them with [myvariant_id()].
#' @param assembly The genome assembly. Changing this from `"hg38"` only makes
#'   sense for genuinely GRCh37 coordinates.
#' @param fields The MyVariant fields to request.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry in
#'   `ids`, in the same order. See [myvariant_parse_batch()].
#'
#' @examples
#' \dontrun{
#' ids <- myvariant_id("17", 7676154, "G", "C")
#' biohttp::body_or_null(myvariant_variants(ids))
#' }
#'
#' @export
myvariant_variants <- function(
  ids,
  assembly = "hg38",
  fields = MYVARIANT_FIELDS,
  ...
) {
  ids <- as.character(ids %||% character())
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(biohttp::status_no_data(
      source = "MyVariant",
      detail = "no variant ids were supplied"
    ))
  }
  if (length(ids) > MYVARIANT_BATCH) {
    stop(
      "`ids` must be at most ",
      MYVARIANT_BATCH,
      " per request, not ",
      length(ids),
      call. = FALSE
    )
  }
  res <- biohttp::post_json(
    paste0(
      MYVARIANT_URL,
      "/variant?assembly=",
      assembly,
      "&fields=",
      utils::URLencode(fields, reserved = TRUE)
    ),
    body = list(ids = paste(ids, collapse = ",")),
    source = "MyVariant",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  biohttp::status_ok(
    data = myvariant_parse_batch(res$data, ids),
    source = "MyVariant",
    http = res$http
  )
}
