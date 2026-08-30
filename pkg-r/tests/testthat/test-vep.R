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
