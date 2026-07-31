# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported HPA fixture parses to the curated tags", {
  body <- read_fixture("hpa_tp53.json")
  out <- hpa_parse_gene(body, "ENSG00000141510")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$symbol, "TP53")
  expect_identical(out$ensembl, "ENSG00000141510")
  expect_identical(out$uniprot[[1]], "P04637")
  expect_true("Cancer-related genes" %in% out$protein_class[[1]])
  expect_true("Tumor suppressor" %in% out$disease_involvement[[1]])
  expect_identical(
    out$source_url,
    "https://www.proteinatlas.org/ENSG00000141510"
  )
})

# --- Editorial judgement is not this package's job ---------------------------

test_that("tags come back as they are, uncounted and unfiltered", {
  # genescout keeps a list of which tags it treats as disease-relevant and
  # counts the matches. Which tags count is an opinion about what makes a gene
  # interesting, and it belongs to the app holding the opinion.
  body <- read_fixture("hpa_tp53.json")
  out <- hpa_parse_gene(body, "ENSG00000141510")

  expect_length(out$protein_class[[1]], 5)
  expect_true("Transcription factors" %in% out$protein_class[[1]])
  expect_false(any(c("n", "score", "present") %in% names(out)))
})

# --- The string-or-array shape -----------------------------------------------

test_that("a single-value field reads the same as a multi-value one", {
  # HPA returns a bare string when there is one value and an array when there
  # are several, so a reader that assumes either shape breaks on the other.
  single <- hpa_parse_gene(
    list(Gene = "X", Ensembl = "ENSG1", `Protein class` = "Solo"),
    "ENSG1"
  )
  many <- hpa_parse_gene(
    list(Gene = "X", Ensembl = "ENSG1", `Protein class` = list("A", "B")),
    "ENSG1"
  )

  expect_identical(single$protein_class[[1]], "Solo")
  expect_identical(many$protein_class[[1]], c("A", "B"))
})

test_that("a missing field is an empty vector, not NULL or NA", {
  out <- hpa_parse_gene(list(Gene = "X", Ensembl = "ENSG1"), "ENSG1")
  expect_identical(out$disease_involvement[[1]], character())
})

test_that("an empty record parses to NULL", {
  expect_null(hpa_parse_gene(list(), "ENSG1"))
  expect_null(hpa_parse_gene(list(Ensembl = "ENSG1"), "ENSG1"))
})

# --- The client half ---------------------------------------------------------

test_that("hpa_gene returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(testthat::test_path("fixtures", "hpa_tp53.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- hpa_gene("ENSG00000141510")

  expect_true(res$ok)
  expect_identical(res$source, "HPA")
  expect_identical(biohttp::body_or_null(res)$symbol, "TP53")
})

test_that("only an Ensembl gene id is accepted", {
  # HPA's path is the id with .json appended, so anything else is a request for
  # an arbitrary page on proteinatlas.org rather than a gene lookup.
  reset_transport()
  expect_identical(hpa_gene("TP53")$status, "no_data")
  expect_identical(hpa_gene("")$status, "no_data")
  expect_identical(hpa_gene("ENSG123/../search")$status, "no_data")
})
