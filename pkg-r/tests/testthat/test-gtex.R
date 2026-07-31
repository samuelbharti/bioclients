# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported GTEx reference fixture parses to the versioned id", {
  body <- read_fixture("gtex_reference_tp53.json")
  out <- gtex_parse_reference(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$symbol[1], "TP53")
  expect_identical(out$gencode_id[1], "ENSG00000141510.16")
  expect_identical(out$entrez[1], "7157")
})

test_that("the ported GTEx expression fixture parses to tissue medians", {
  body <- read_fixture("gtex_tp53.json")
  out <- gtex_parse_expression(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$tissue_id[1], "Adipose_Subcutaneous")
  expect_identical(out$tissue[1], "Adipose Subcutaneous")
  expect_equal(out$median_tpm[1], 22.459)
  expect_identical(out$gencode_id[1], "ENSG00000141510.16")
})

test_that("the second recorded expression fixture parses too", {
  body <- read_fixture("gtex_median_tp53.json")
  out <- gtex_parse_expression(body)

  expect_identical(out$tissue[1], "Nerve Tibial")
  expect_equal(out$median_tpm[1], 30.5)
})

test_that("a row with no median is dropped rather than carried as NA", {
  body <- list(
    data = list(
      list(median = 5, tissueSiteDetailId = "Liver"),
      list(tissueSiteDetailId = "Brain_Cortex")
    )
  )
  out <- gtex_parse_expression(body)

  expect_identical(nrow(out), 1L)
  expect_identical(out$tissue, "Liver")
})

test_that("an empty response parses to NULL", {
  expect_null(gtex_parse_reference(list(data = list())))
  expect_null(gtex_parse_expression(list(data = list())))
})

# --- The versioned-id requirement --------------------------------------------
# GTEx needs a versioned GENCODE id AND a datasetId. An unversioned Ensembl id,
# which is what MyGene and most other sources hand you, returns nothing. There
# is no way to derive the version, so it has to be looked up.

test_that("an expression lookup resolves the versioned id first", {
  reset_transport()
  urls <- character()
  reference <- paste(
    readLines(
      testthat::test_path("fixtures", "gtex_reference_tp53.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  expression <- paste(
    readLines(testthat::test_path("fixtures", "gtex_tp53.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    if (grepl("reference/gene", req$url, fixed = TRUE)) {
      return(mock_json(reference))
    }
    mock_json(expression)
  })

  res <- gtex_median_expression("TP53")

  expect_true(res$ok)
  expect_length(urls, 2)
  expect_match(urls[1], "reference/gene")
  # The versioned id from the first call is what the second one asks with.
  expect_match(urls[2], "ENSG00000141510\\.16")
  expect_match(urls[2], "datasetId=gtex_v8")
})

test_that("a known versioned id skips the reference lookup", {
  reset_transport()
  urls <- character()
  expression <- paste(
    readLines(testthat::test_path("fixtures", "gtex_tp53.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json(expression)
  })

  res <- gtex_median_expression(gencode_id = "ENSG00000141510.16")

  expect_true(res$ok)
  expect_length(urls, 1)
  expect_false(grepl("reference/gene", urls[1], fixed = TRUE))
})

test_that("a gene GTEx does not know stops before the expression call", {
  reset_transport()
  calls <- 0L
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    mock_json('{"data":[]}')
  })

  res <- gtex_median_expression("NOSUCHGENE")

  expect_identical(res$status, "no_data")
  # One request. Asking for expression with no id could not have succeeded.
  expect_identical(calls, 1L)
})

test_that("a blank gene is no_data and never reaches the network", {
  reset_transport()
  expect_identical(gtex_gene_reference("")$status, "no_data")
  expect_identical(gtex_median_expression("")$status, "no_data")
})
