# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported MyVariant batch fixture parses as expected", {
  body <- read_fixture("myvariant_batch.json")
  ids <- vapply(body, function(record) record$query, character(1))
  out <- myvariant_parse_batch(body, ids)

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), length(ids))
  expect_identical(out$id[1], "chr17:g.7676154G>C")
  expect_identical(out$gene[1], "TP53")
})

# --- Trap one: CADD's rankscore key ------------------------------------------

test_that("CADD is read from raw_rankscore, which is the only key there is", {
  # Every other dbNSFP predictor is <name>.rankscore. CADD is
  # cadd.raw_rankscore, and there is no plain cadd.rankscore. Following the
  # pattern gives NA for CADD on every variant, forever, with no error.
  body <- read_fixture("myvariant_batch.json")
  raw <- body[[1]]$dbnsfp$cadd

  expect_true("raw_rankscore" %in% names(raw))
  expect_false("rankscore" %in% names(raw))

  out <- myvariant_parse_batch(body, body[[1]]$query)
  expect_equal(out$cadd[1], raw$raw_rankscore)
  expect_false(is.na(out$cadd[1]))
})

test_that("the other predictors do use plain rankscore", {
  # The asymmetry is real, not a mistake in the parser.
  record <- list(
    dbnsfp = list(
      clinpred = list(rankscore = 0.5),
      revel = list(rankscore = 0.6),
      alphamissense = list(rankscore = 0.7),
      cadd = list(raw_rankscore = 0.8)
    )
  )
  out <- myvariant_parse_record(record)

  expect_equal(out$clinpred, 0.5)
  expect_equal(out$revel, 0.6)
  expect_equal(out$alphamissense, 0.7)
  expect_equal(out$cadd, 0.8)
})

# --- Trap two: assembly=hg38 -------------------------------------------------

test_that("assembly=hg38 is always sent", {
  # Omit it and MyVariant answers 200 with notfound for every GRCh38 variant.
  # Nothing errors and nothing warns; a whole cohort simply disappears.
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json("[]")
  })

  myvariant_variants("chr17:g.7676154G>C")

  expect_match(urls[1], "assembly=hg38", fixed = TRUE)
})

test_that("notfound is no_data for that row, not an error for the call", {
  # MyVariant reports an unknown variant as 200 + notfound. That is the source
  # saying it has nothing, which is different from the call failing.
  body <- list(
    list(query = "chr1:g.1A>T", notfound = TRUE),
    list(query = "chr2:g.2C>G", dbnsfp = list(genename = list("BRCA2")))
  )
  out <- myvariant_parse_batch(body, c("chr1:g.1A>T", "chr2:g.2C>G"))

  expect_identical(nrow(out), 2L)
  expect_true(is.na(out$gene[1]))
  expect_identical(out$gene[2], "BRCA2")
})

test_that("a variant absent from the response keeps its row", {
  # Dropping it would shift every row after it onto the wrong variant.
  body <- list(list(query = "chr2:g.2C>G", dbnsfp = list(genename = list("X"))))
  out <- myvariant_parse_batch(body, c("chr1:g.1A>T", "chr2:g.2C>G"))

  expect_identical(out$id, c("chr1:g.1A>T", "chr2:g.2C>G"))
  expect_true(is.na(out$gene[1]))
})

# --- The verified batch cap --------------------------------------------------

test_that("more than 1000 ids is refused rather than silently truncated", {
  reset_transport()
  expect_error(
    myvariant_variants(rep("chr1:g.1A>T", MYVARIANT_BATCH + 1)),
    "at most 1000"
  )
})

test_that("the cap matches the verified limit", {
  expect_identical(MYVARIANT_BATCH, 1000L)
})

# --- Id building -------------------------------------------------------------

test_that("an SNV becomes a substitution id", {
  expect_identical(myvariant_id("17", 7676154, "G", "C"), "chr17:g.7676154G>C")
  expect_identical(
    myvariant_id("chr17", 7676154, "G", "C"),
    "chr17:g.7676154G>C"
  )
})

test_that("a deletion and an insertion get their own forms", {
  expect_identical(myvariant_id("1", 100, "AT", "A"), "chr1:g.101_101del")
  expect_identical(myvariant_id("1", 100, "A", "AT"), "chr1:g.100_101insT")
})

# --- Field shapes ------------------------------------------------------------

test_that("a per-transcript genename array reduces to one symbol", {
  record <- list(dbnsfp = list(genename = list("TP53", "TP53", "TP53")))
  expect_identical(myvariant_parse_record(record)$gene, "TP53")
})

test_that("clinvar significance reads both the list and the scalar shape", {
  many <- list(
    clinvar = list(
      rcv = list(
        list(clinical_significance = "Benign"),
        list(clinical_significance = "Pathogenic")
      )
    )
  )
  one <- list(clinvar = list(rcv = list(clinical_significance = "Pathogenic")))

  expect_identical(myvariant_parse_record(many)$clinvar_sig, "Benign")
  expect_identical(myvariant_parse_record(one)$clinvar_sig, "Pathogenic")
})

# --- The client half ---------------------------------------------------------

test_that("no ids is no_data and never reaches the network", {
  reset_transport()
  expect_identical(myvariant_variants(character())$status, "no_data")
  expect_identical(myvariant_variants(NULL)$status, "no_data")
})
