# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported VEP fixture parses as expected", {
  body <- read_fixture("vep_region.json")
  keys <- vapply(body, vep_element_key, character(1))
  out <- vep_parse_batch(body, keys)

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), length(keys))
  expect_true(any(!is.na(out$consequence)))
})

# --- Trap: AlphaMissense is nested -------------------------------------------

test_that("AlphaMissense is read from inside the transcript", {
  # It lives at transcript_consequences[].alphamissense, never at the top level.
  # Reading it from the top finds nothing, which reads as "VEP does not serve
  # AlphaMissense" rather than as a bug.
  element <- list(
    most_severe_consequence = "missense_variant",
    transcript_consequences = list(list(
      gene_symbol = "BRAF",
      consequence_terms = list("missense_variant"),
      alphamissense = list(am_pathogenicity = 0.99, am_class = "pathogenic")
    ))
  )
  out <- vep_parse_element(element)

  expect_equal(out$alphamissense, 0.99)
  expect_identical(out$alphamissense_class, "pathogenic")
})

test_that("a top-level AlphaMissense is not where it is looked for", {
  # Proving the nesting matters: the same numbers at the top level yield NA.
  element <- list(
    most_severe_consequence = "missense_variant",
    alphamissense = list(am_pathogenicity = 0.99),
    transcript_consequences = list(list(
      consequence_terms = list("missense_variant")
    ))
  )
  expect_true(is.na(vep_parse_element(element)$alphamissense))
})

# --- Trap: the body must be a JSON array -------------------------------------

test_that("variants serialises as an array, not an object", {
  # A named R list becomes a JSON object and VEP answers 500. unname() is
  # load-bearing.
  reset_transport()
  sent <- NULL
  httr2::local_mocked_responses(function(req) {
    sent <<- req$body$data
    mock_json("[]")
  })

  vep_variants(c("7", "17"), c(140753336, 7676154), c("A", "G"), c("T", "C"))

  expect_null(names(sent$variants))
  expect_type(sent$variants, "list")
  expect_length(sent$variants, 2)
})

# --- Trap: the verified batch cap --------------------------------------------

test_that("more than 200 variants is refused rather than truncated", {
  reset_transport()
  n <- VEP_BATCH + 1
  expect_error(
    vep_variants(rep("1", n), seq_len(n), rep("A", n), rep("T", n)),
    "at most 200"
  )
})

test_that("the cap matches the verified limit", {
  expect_identical(VEP_BATCH, 200L)
})

# --- Trap: indels are renumbered ---------------------------------------------

test_that("every input row of the recorded indel batch comes back populated", {
  # Recorded live from rest.ensembl.org with AlphaMissense=1&mane=1&numbers=1&
  # vcf_string=1 for one SNV, one insertion, one deletion and one delins on
  # public sites: TP53 p.Pro72Arg, a PCSK9 intronic insertion, CFTR p.Phe508del
  # and a TP53 exonic delins. VEP renumbers three of the four, so a key rebuilt
  # from the reported coordinates matched only the SNV.
  body <- read_fixture("vep_region_indels.json")
  keys <- c(
    vep_key("17", 7676154, "G", "C"),
    vep_key("1", 55516888, "T", "TA"),
    vep_key("7", 117559590, "ATCT", "A"),
    vep_key("17", 7675088, "CGC", "TA")
  )
  out <- vep_parse_batch(body, keys)

  expect_identical(nrow(out), 4L)
  expect_identical(out$key, keys)
  expect_false(anyNA(out$consequence))
  expect_identical(
    out$consequence,
    c(
      "missense_variant",
      "intron_variant",
      "inframe_deletion",
      "frameshift_variant"
    )
  )
  expect_identical(out$gene[c(1, 3, 4)], c("TP53", "CFTR", "TP53"))
})

test_that("the key VEP rebuilds for an indel is not the key that was asked", {
  # This is the bug. The insertion sent as 1 55516888 T>TA is reported with
  # start 55516889 and allele_string -/A.
  body <- read_fixture("vep_region_indels.json")
  insertion <- body[[2]]

  expect_identical(insertion$start, 55516889L)
  expect_identical(insertion$allele_string, "-/A")
  expect_false(
    identical(vep_element_key(insertion), vep_key("1", 55516888, "T", "TA"))
  )
})

test_that("the echoed input is matched before the vcf_string", {
  # VEP re-anchors a delins in vcf_string: CGC>TA at 7675088 is reported as
  # 17-7675087-GCGC-GTA. Only the echoed input carries the key that was asked.
  element <- list(
    input = "17 7675088 . CGC TA . . .",
    vcf_string = "17-7675087-GCGC-GTA",
    seq_region_name = "17",
    start = 7675088,
    allele_string = "CGC/TA",
    most_severe_consequence = "frameshift_variant",
    transcript_consequences = list(list(
      gene_symbol = "TP53",
      consequence_terms = list("frameshift_variant")
    ))
  )
  out <- vep_parse_batch(list(element), vep_key("17", 7675088, "CGC", "TA"))

  expect_identical(out$gene, "TP53")
  expect_identical(
    vep_element_keys(element)[1],
    vep_key("17", 7675088, "CGC", "TA")
  )
})

test_that("an element without an input line still matches on vcf_string", {
  element <- list(
    vcf_string = "1-55516888-T-TA",
    seq_region_name = "1",
    start = 55516889,
    allele_string = "-/A",
    most_severe_consequence = "intron_variant",
    transcript_consequences = list()
  )
  out <- vep_parse_batch(list(element), vep_key("1", 55516888, "T", "TA"))

  expect_identical(out$consequence, "intron_variant")
})

test_that("vcf_string is requested so the second identity is present", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[]")
  })

  vep_variants("7", 140753336, "A", "T")

  expect_match(url, "vcf_string=1", fixed = TRUE)
})

# --- Request options ---------------------------------------------------------

test_that("the default options are the four flags the parsers rely on", {
  expect_identical(
    vep_default_options(),
    list(AlphaMissense = 1, mane = 1, numbers = 1, vcf_string = 1)
  )
})

test_that("options become a query string, with off flags omitted", {
  options <- c(
    vep_default_options(),
    list(CADD = TRUE, REVEL = FALSE, af = 0, hgvs = NULL, SpliceAI = 1)
  )
  expect_identical(
    vep_query_string(options),
    "AlphaMissense=1&mane=1&numbers=1&vcf_string=1&CADD=1&SpliceAI=1"
  )
  expect_identical(vep_query_string(list()), "")
  expect_identical(vep_query_string(NULL), "")
})

test_that("a valued option is carried through and encoded", {
  expect_identical(
    vep_query_string(list(pick_order = "mane_select,canonical")),
    "pick_order=mane_select%2Ccanonical"
  )
})

test_that("dbNSFP is refused", {
  # Trap 4 in the file header. dbNSFP values come back in dbNSFP's own
  # transcript order, not aligned to the transcript VEP reports, so they
  # silently describe a different transcript than the rest of the row.
  expect_error(vep_query_string(list(dbNSFP = "REVEL_score")), "dbNSFP")
  reset_transport()
  expect_error(
    vep_variants("7", 140753336, "A", "T", options = list(dbNSFP = "x")),
    "dbNSFP"
  )
})

test_that("unnamed options are refused rather than sent as garbage", {
  expect_error(vep_query_string(list(1, 2)), "named")
})

test_that("vep_variants sends the options on the query string", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[]")
  })

  vep_variants(
    "7",
    140753336,
    "A",
    "T",
    options = c(vep_default_options(), list(af_gnomadg = 1, CADD = 1))
  )

  expect_match(url, "AlphaMissense=1", fixed = TRUE)
  expect_match(url, "af_gnomadg=1", fixed = TRUE)
  expect_match(url, "CADD=1", fixed = TRUE)
})

test_that("the default request asks for exactly the default flags", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[]")
  })

  vep_variants("7", 140753336, "A", "T")

  expect_identical(
    sub("^.*\\?", "", url),
    "AlphaMissense=1&mane=1&numbers=1&vcf_string=1"
  )
})

# --- The transcript columns the optional flags add ---------------------------

test_that("the recorded frequency batch fills the predictor columns", {
  # Recorded live with the defaults plus af_gnomadg=1&af_gnomade=1&CADD=1&
  # SpliceAI=1&REVEL=1&hgvs=1 for the same four sites as vep_region_indels.
  body <- read_fixture("vep_region_indels_freq.json")
  keys <- c(
    vep_key("17", 7676154, "G", "C"),
    vep_key("1", 55516888, "T", "TA"),
    vep_key("7", 117559590, "ATCT", "A"),
    vep_key("17", 7675088, "CGC", "TA")
  )
  out <- vep_parse_batch(body, keys)

  expect_identical(out$key, keys)
  expect_identical(
    out$transcript,
    c(
      "ENST00000269305",
      "ENST00000643167",
      "ENST00000003084",
      "ENST00000269305"
    )
  )
  expect_identical(out$gene_id[1], "ENSG00000141510")
  expect_identical(out$biotype[2], "lncRNA")
  expect_identical(out$hgvsc[1], "ENST00000269305.9:c.215C>G")
  expect_identical(out$hgvsp[1], "ENSP00000269305.4:p.Pro72Arg")
  expect_identical(out$hgvsp[3], "ENSP00000003084.6:p.Phe508del")
  expect_identical(out$hgvsp[4], "ENSP00000269305.4:p.Arg174SerfsTer73")
  expect_identical(out$codons[1], "cCc/cGc")
  expect_identical(out$amino_acids[1], "P/R")
  expect_equal(out$cadd_phred[1], 12.91)
  expect_equal(out$cadd_raw[1], 1.299883)
  expect_equal(out$revel[1], 0.368)
  expect_equal(out$cadd_phred[3], 17.55)
  # SpliceAI reports four delta scores per transcript. All zero here, which is
  # an answer rather than an absence.
  expect_equal(out$spliceai_ds_ag[1], 0)
  expect_equal(out$spliceai_max[1], 0)
  expect_true(is.na(out$spliceai_max[2]))
  # An intronic lncRNA insertion has no protein annotation and no scores.
  expect_true(is.na(out$hgvsp[2]))
  expect_true(is.na(out$cadd_phred[2]))
})

test_that("spliceai_max is the largest of the four delta scores", {
  element <- list(
    most_severe_consequence = "splice_region_variant",
    transcript_consequences = list(list(
      consequence_terms = list("splice_region_variant"),
      spliceai = list(DS_AG = 0.01, DS_AL = 0.62, DS_DG = 0, DS_DL = 0.2)
    ))
  )
  out <- vep_parse_element(element)

  expect_equal(out$spliceai_ds_al, 0.62)
  expect_equal(out$spliceai_max, 0.62)
})

test_that("canonical is TRUE on the marked transcript and NA elsewhere", {
  marked <- list(
    most_severe_consequence = "missense_variant",
    transcript_consequences = list(list(
      consequence_terms = list("missense_variant"),
      canonical = 1L
    ))
  )
  unmarked <- list(
    most_severe_consequence = "missense_variant",
    transcript_consequences = list(list(
      consequence_terms = list("missense_variant")
    ))
  )
  expect_true(vep_parse_element(marked)$canonical)
  expect_true(is.na(vep_parse_element(unmarked)$canonical))
})

test_that("the optional columns are NA when the flags were not asked", {
  body <- read_fixture("vep_region.json")
  out <- vep_parse_element(body[[1]])

  expect_identical(out$transcript, "ENST00000269305")
  expect_true(is.na(out$hgvsc))
  expect_true(is.na(out$cadd_phred))
  expect_true(is.na(out$revel))
  expect_true(is.na(out$spliceai_max))
  expect_true(is.na(out$lof))
})

test_that("the LOFTEE call is read into lof", {
  element <- list(
    most_severe_consequence = "stop_gained",
    transcript_consequences = list(list(
      consequence_terms = list("stop_gained"),
      lof = "HC"
    ))
  )
  expect_identical(vep_parse_element(element)$lof, "HC")
})

# --- Colocated variants ------------------------------------------------------

test_that("the recorded frequency batch fills the colocated columns", {
  body <- read_fixture("vep_region_indels_freq.json")
  keys <- c(
    vep_key("17", 7676154, "G", "C"),
    vep_key("1", 55516888, "T", "TA"),
    vep_key("7", 117559590, "ATCT", "A"),
    vep_key("17", 7675088, "CGC", "TA")
  )
  out <- vep_parse_batch(body, keys)

  expect_identical(out$rsid, c("rs1042522", NA, "rs113993960", NA))
  expect_equal(out$gnomadg_af[1], 0.6268)
  expect_equal(out$gnomade_af[1], 0.7163)
  # The largest per-population value, not the overall one.
  expect_equal(out$gnomadg_af_max[1], 0.745)
  expect_equal(out$gnomade_af_max[1], 0.7477)
  expect_gt(out$gnomadg_af_max[1], out$gnomadg_af[1])
  expect_equal(out$gnomadg_af[3], 0.007884)
  expect_identical(
    out$clin_sig[3],
    "risk_factor;pathogenic;drug_response;likely_pathogenic"
  )
  expect_true(all(is.na(out$gnomadg_af[c(2, 4)])))
})

test_that("clin_sig is the significance of the element's own allele", {
  # rs1042522 is multi-allelic. The C allele is what was asked, and VEP lists
  # the significance of every allele in clin_sig, so reading that would hand
  # the T allele's conflicting call to the C allele.
  body <- read_fixture("vep_region_indels_freq.json")
  out <- vep_parse_colocated(body[[1]])

  expect_identical(out$clin_sig, "pathogenic;benign")
  expect_false(grepl("conflicting", out$clin_sig, fixed = TRUE))
})

test_that("frequencies are read for the element's allele", {
  element <- list(
    allele_string = "G/T",
    colocated_variants = list(list(
      id = "rs1",
      frequencies = list(
        C = list(gnomadg = 0.6),
        T = list(gnomadg = 0.01, gnomadg_afr = 0.02)
      )
    ))
  )
  out <- vep_parse_colocated(element)

  expect_equal(out$gnomadg_af, 0.01)
  expect_equal(out$gnomadg_af_max, 0.02)
})

test_that("a deletion's frequencies are keyed by the dash allele", {
  body <- read_fixture("vep_region_indels_freq.json")
  deletion <- body[[3]]

  expect_identical(deletion$allele_string, "TCT/-")
  expect_identical(names(deletion$colocated_variants[[1]]$frequencies), "-")
  expect_equal(vep_parse_colocated(deletion)$gnomadg_af, 0.007884)
})

test_that("the dbSNP record is read, not the COSMIC or HGMD entry beside it", {
  body <- read_fixture("vep_region_indels_freq.json")
  ids <- vapply(body[[1]]$colocated_variants, function(co) co$id, character(1))

  expect_true(any(grepl("^COSV", ids)))
  expect_identical(vep_parse_colocated(body[[1]])$rsid, "rs1042522")
})

test_that("an element with only somatic entries has no rsid", {
  element <- list(
    allele_string = "G/C",
    colocated_variants = list(list(id = "COSV123", somatic = 1))
  )
  out <- vep_parse_colocated(element)

  expect_true(is.na(out$rsid))
  expect_true(is.na(out$gnomadg_af))
})

test_that("no colocated variants is a row of NA", {
  out <- vep_parse_colocated(list(allele_string = "G/C"))
  expect_identical(nrow(out), 1L)
  expect_true(all(is.na(out)))
})

test_that("clin_sig falls back to the union when there is no per-allele form", {
  element <- list(
    allele_string = "G/C",
    colocated_variants = list(list(
      id = "rs1",
      clin_sig = list("benign", "likely_benign")
    ))
  )
  expect_identical(
    vep_parse_colocated(element)$clin_sig,
    "benign;likely_benign"
  )
})

test_that("the batch row carries both halves and the empty row matches", {
  body <- read_fixture("vep_region_indels_freq.json")
  out <- vep_parse_batch(body, c(vep_key("17", 7676154, "G", "C"), "9-1-A-T"))

  expect_true(all(c("hgvsc", "rsid", "gnomadg_af") %in% names(out)))
  expect_identical(nrow(out), 2L)
  expect_true(is.na(out$rsid[2]))
})

# --- Trap: order is not promised ---------------------------------------------

test_that("results are matched by identity, not by array position", {
  # VEP does not promise response order. Zipping by index assigns consequences
  # to the wrong variants, silently.
  body <- list(
    list(
      seq_region_name = "17",
      start = 7676154,
      allele_string = "G/C",
      most_severe_consequence = "missense_variant",
      transcript_consequences = list(list(
        gene_symbol = "TP53",
        consequence_terms = list("missense_variant")
      ))
    ),
    list(
      seq_region_name = "7",
      start = 140753336,
      allele_string = "A/T",
      most_severe_consequence = "missense_variant",
      transcript_consequences = list(list(
        gene_symbol = "BRAF",
        consequence_terms = list("missense_variant")
      ))
    )
  )
  # Asked in the opposite order to the response.
  keys <- c(vep_key("7", 140753336, "A", "T"), vep_key("17", 7676154, "G", "C"))
  out <- vep_parse_batch(body, keys)

  expect_identical(out$gene, c("BRAF", "TP53"))
})

test_that("a variant VEP returned nothing for keeps its row", {
  keys <- c(vep_key("1", 1, "A", "T"), vep_key("2", 2, "C", "G"))
  out <- vep_parse_batch(list(), keys)

  expect_identical(nrow(out), 2L)
  expect_true(all(is.na(out$gene)))
})

# --- Transcript choice -------------------------------------------------------

test_that("MANE Select wins over everything else", {
  element <- list(
    most_severe_consequence = "missense_variant",
    transcript_consequences = list(
      list(gene_symbol = "A", consequence_terms = list("missense_variant")),
      list(
        gene_symbol = "B",
        consequence_terms = list("intron_variant"),
        mane_select = "NM_1"
      )
    )
  )
  expect_identical(vep_parse_element(element)$gene, "B")
})

test_that("without MANE, the most severe consequence wins", {
  # Taking the first transcript outright reports an arbitrary one, often
  # non-coding, which makes the consequence look milder than it is.
  element <- list(
    most_severe_consequence = "stop_gained",
    transcript_consequences = list(
      list(gene_symbol = "A", consequence_terms = list("intron_variant")),
      list(gene_symbol = "B", consequence_terms = list("stop_gained"))
    )
  )
  expect_identical(vep_parse_element(element)$gene, "B")
})

test_that("a variant with no transcripts still reports its consequence", {
  element <- list(
    most_severe_consequence = "intergenic_variant",
    transcript_consequences = list()
  )
  out <- vep_parse_element(element)

  expect_identical(out$consequence, "intergenic_variant")
  expect_true(is.na(out$gene))
})

# --- Region strings ----------------------------------------------------------

test_that("a region string is the VCF-like form VEP expects", {
  expect_identical(
    vep_region("7", 140753336, "A", "T"),
    "7 140753336 . A T . . ."
  )
  expect_identical(
    vep_region("chr7", 140753336, "a", "t"),
    "7 140753336 . A T . . ."
  )
})

test_that("no variants is no_data and never reaches the network", {
  reset_transport()
  expect_identical(
    vep_variants(character(), integer(), character(), character())$status,
    "no_data"
  )
})

# --- Chunking through post_json_many -----------------------------------------

# A mock that answers each chunk with one element per variant it was sent,
# echoing the input line the way VEP does, so the match-back has something to
# match. `fail_when` makes a chunk fail when it returns TRUE for its inputs.
vep_chunk_mock <- function(seen, fail_when = function(inputs) FALSE) {
  function(req) {
    inputs <- unlist(req$body$data$variants, use.names = FALSE)
    seen$calls <- c(seen$calls, list(inputs))
    if (isTRUE(fail_when(inputs))) {
      return(httr2::response(status_code = 500))
    }
    elements <- lapply(inputs, function(input) {
      parts <- strsplit(input, " ", fixed = TRUE)[[1]]
      list(
        input = input,
        most_severe_consequence = "missense_variant",
        transcript_consequences = list(list(
          gene_symbol = paste0("GENE", parts[2]),
          consequence_terms = list("missense_variant")
        ))
      )
    })
    mock_json(jsonlite::toJSON(elements, auto_unbox = TRUE))
  }
}

test_that("vep_variants_all chunks the input and keeps input order", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(vep_chunk_mock(seen))

  res <- vep_variants_all(
    rep("1", 5),
    1:5,
    rep("A", 5),
    rep("T", 5),
    chunk_size = 2
  )

  expect_true(res$ok)
  expect_length(seen$calls, 3)
  expect_identical(lengths(seen$calls), c(2L, 2L, 1L))
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), 5L)
  expect_identical(out$key, vep_key("1", 1:5, "A", "T"))
  expect_identical(out$gene, paste0("GENE", 1:5))
  expect_identical(out$status, rep("ok", 5))
})

test_that("a failed chunk degrades to NA rows with its status", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(vep_chunk_mock(
    seen,
    fail_when = function(inputs) any(grepl("^1 3 ", inputs))
  ))

  res <- vep_variants_all(
    rep("1", 5),
    1:5,
    rep("A", 5),
    rep("T", 5),
    chunk_size = 2
  )

  expect_true(res$ok)
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), 5L)
  # Rows 3 and 4 were the failed chunk. They keep their place.
  expect_identical(out$gene, c("GENE1", "GENE2", NA, NA, "GENE5"))
  expect_identical(out$status, c("ok", "ok", "error", "error", "ok"))
  reset_transport()
})

test_that("when every chunk fails the failing envelope comes back", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 503)
  })

  res <- vep_variants_all("1", 1, "A", "T")

  expect_false(res$ok)
  expect_identical(res$status, "error")
  expect_identical(res$source, "VEP")
  reset_transport()
})

test_that("the chunk size may not exceed the verified limit", {
  reset_transport()
  expect_error(
    vep_variants_all("1", 1, "A", "T", chunk_size = VEP_BATCH + 1),
    "chunk_size"
  )
  expect_error(vep_variants_all("1", 1, "A", "T", chunk_size = 0), "chunk_size")
})

test_that("vep_variants_all sends the options and refuses dbNSFP", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[]")
  })

  vep_variants_all(
    "1",
    1,
    "A",
    "T",
    options = c(vep_default_options(), list(REVEL = 1))
  )
  expect_match(url, "REVEL=1", fixed = TRUE)

  expect_error(
    vep_variants_all("1", 1, "A", "T", options = list(dbNSFP = "x")),
    "dbNSFP"
  )
})

test_that("no variants is no_data without a request", {
  reset_transport()
  expect_identical(
    vep_variants_all(character(), integer(), character(), character())$status,
    "no_data"
  )
})

# --- Sorted, deduplicated, and throttled before dispatch ---------------------

test_that("vep_variants_all sorts into genome order before chunking", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(vep_chunk_mock(seen))

  # Deliberately out of genome order: chrX before chr1, descending position
  # within chr1.
  res <- vep_variants_all(
    c("X", "1", "1"),
    c(500, 300, 100),
    c("A", "A", "A"),
    c("T", "T", "T"),
    chunk_size = 2
  )

  expect_true(res$ok)
  # Chunk 1 dispatches chr1:100 then chr1:300, in genome order; chunk 2
  # dispatches chrX:500.
  expect_identical(
    seen$calls[[1]],
    c("1 100 . A T . . .", "1 300 . A T . . .")
  )
  expect_identical(seen$calls[[2]], "X 500 . A T . . .")

  out <- biohttp::body_or_null(res)
  # The output stays in the ORIGINAL input order, regardless of dispatch
  # order.
  expect_identical(
    out$key,
    vep_key(
      c("X", "1", "1"),
      c(500, 300, 100),
      c("A", "A", "A"),
      c("T", "T", "T")
    )
  )
})

test_that("a repeated variant is sent once and every occurrence gets the same answer", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(vep_chunk_mock(seen))

  # Positions 1 and 3 are the same variant, as a multiallelic split can
  # produce.
  res <- vep_variants_all(
    c("1", "1", "1"),
    c(100, 200, 100),
    c("A", "A", "A"),
    c("T", "G", "T"),
    chunk_size = 10
  )

  expect_true(res$ok)
  # Two distinct variants dispatched, not three.
  expect_length(unlist(seen$calls), 2)
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), 3L)
  expect_identical(out$gene[[1]], out$gene[[3]])
  expect_identical(out$status, rep("ok", 3))
})

test_that("vep_variants_all applies a default throttle for VEP's host", {
  reset_transport()
  captured <- "unset"
  httr2::local_mocked_responses(function(req) {
    captured <<- req$policies$throttle_realm
    mock_json("[]")
  })

  vep_variants_all("1", 1, "A", "T")
  expect_identical(captured, "rest.ensembl.org")

  # A fresh cache, or the second call (same URL and body as the first) would
  # be served from the success-only cache and never re-dispatch.
  reset_transport()
  captured <- "unset"
  vep_variants_all("1", 1, "A", "T", throttle = NULL)
  expect_null(captured)
})
