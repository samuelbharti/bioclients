# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported DGIdb fixture parses to the interaction count", {
  body <- read_fixture("dgidb_nf1.json")
  out <- dgidb_parse_genes(body, "NF1")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$symbol, "NF1")
  expect_identical(out$concept_id, "hgnc:7765")
  expect_identical(out$interaction_count, 2L)
  expect_identical(out$source_url, "https://dgidb.org/genes/hgnc:7765")
})

# --- A real zero is not a miss -----------------------------------------------
# The distinction the whole parser exists to preserve. Collapsing them would
# report an unknown gene as known-not-druggable, which is a much stronger claim
# than the data supports.

test_that("a gene DGIdb knows with no interactions is a real zero", {
  body <- read_fixture("dgidb_zero.json")
  out <- dgidb_parse_genes(body, "XYZ")

  expect_identical(out$interaction_count, 0L)
  expect_false(is.na(out$interaction_count))
  expect_identical(out$concept_id, "hgnc:1")
})

test_that("a gene DGIdb has never heard of is NA, not zero", {
  body <- read_fixture("dgidb_empty.json")
  out <- dgidb_parse_genes(body, "NOSUCHGENE")

  expect_true(is.na(out$interaction_count))
  expect_true(is.na(out$concept_id))
  expect_true(is.na(out$source_url))
})

test_that("the two cases are distinguishable in one batch", {
  body <- read_fixture("dgidb_zero.json")
  out <- dgidb_parse_genes(body, c("XYZ", "NOSUCHGENE"))

  expect_identical(nrow(out), 2L)
  expect_identical(out$interaction_count[1], 0L)
  expect_true(is.na(out$interaction_count[2]))
})

# --- Order and alignment -----------------------------------------------------

test_that("rows come back in input order regardless of response order", {
  # DGIdb returns nodes in its own order, so a caller zipping by position needs
  # this. Matching is by name, not by index.
  body <- list(
    data = list(
      genes = list(
        nodes = list(
          list(
            name = "BRAF",
            conceptId = "hgnc:1097",
            interactions = list(list())
          ),
          list(name = "NF1", conceptId = "hgnc:7765", interactions = list())
        )
      )
    )
  )
  out <- dgidb_parse_genes(body, c("NF1", "BRAF"))

  expect_identical(out$symbol, c("NF1", "BRAF"))
  expect_identical(out$interaction_count, c(0L, 1L))
})

test_that("symbols are matched case-insensitively", {
  body <- read_fixture("dgidb_nf1.json")
  expect_identical(dgidb_parse_genes(body, "nf1")$interaction_count, 2L)
})

# --- The client half ---------------------------------------------------------

test_that("dgidb_genes sends the whole list in one request", {
  # The query is batch-shaped and both source apps passed one symbol at a time.
  # Sending N requests for N genes is the thing this client exists to stop.
  reset_transport()
  calls <- 0L
  sent <- NULL
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    sent <<- req$body$data$variables$names
    mock_json('{"data":{"genes":{"nodes":[]}}}')
  })

  dgidb_genes(c("NF1", "BRAF", "EGFR"))

  expect_identical(calls, 1L)
  expect_length(sent, 3)
})

test_that("dgidb_gene returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(testthat::test_path("fixtures", "dgidb_nf1.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- dgidb_gene("NF1")

  expect_true(res$ok)
  expect_identical(res$source, "DGIdb")
  expect_identical(biohttp::body_or_null(res)$interaction_count, 2L)
})

test_that("a blank symbol is no_data and never reaches the network", {
  reset_transport()
  expect_identical(dgidb_gene("")$status, "no_data")
})
