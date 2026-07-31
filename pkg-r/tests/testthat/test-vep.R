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
