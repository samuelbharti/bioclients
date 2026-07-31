# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported QuickGO fixture parses as expected", {
  body <- read_fixture("quickgo_nf1.json")
  out <- quickgo_parse_annotations(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$go_id[1], "GO:0001937")
  expect_identical(
    out$go_name[1],
    "negative regulation of endothelial cell proliferation"
  )
  expect_identical(out$evidence[1], "IMP")
  expect_identical(out$reference[1], "PMID:16648142")
})

# --- One row per term, not per evidence line ---------------------------------

test_that("a term repeated across evidence lines collapses to one row", {
  # GO repeats a term once per supporting piece of evidence, so a gene with 20
  # distinct functions can arrive as 100 rows. Counting raw rows overstates how
  # many things a gene is annotated to.
  body <- read_fixture("quickgo_nf1.json")
  out <- quickgo_parse_annotations(body)

  expect_lt(nrow(out), length(body$results))
  expect_false(any(duplicated(out$go_id)))
})

test_that("the first evidence line for a term is the one kept", {
  body <- list(
    results = list(
      list(goId = "GO:1", goName = "a", goEvidence = "IMP"),
      list(goId = "GO:1", goName = "a", goEvidence = "IEA")
    )
  )
  expect_identical(quickgo_parse_annotations(body)$evidence, "IMP")
})

# --- Only GO ids count -------------------------------------------------------

test_that("a row whose goId is not a GO id is dropped", {
  body <- list(
    results = list(
      list(goId = "UniProtKB:P21359", goName = "not a term"),
      list(goId = "GO:0001937", goName = "real")
    )
  )
  expect_identical(quickgo_parse_annotations(body)$go_id, "GO:0001937")
})

test_that("no results is NULL", {
  expect_null(quickgo_parse_annotations(list(results = list())))
  expect_null(quickgo_parse_annotations(list()))
})

# --- The request -------------------------------------------------------------

test_that("goName is requested, because ids alone are unreadable", {
  # Without includeFields=goName, QuickGO sends the GO id and no label, which is
  # enough to count terms and not enough to read them.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"results":[]}')
  })

  quickgo_annotations("P21359")

  expect_match(url, "includeFields=goName", fixed = TRUE)
  expect_match(url, "geneProductId=UniProtKB%3AP21359")
  expect_match(url, "aspect=biological_process", fixed = TRUE)
})

test_that("the aspect is checked rather than pasted into the query", {
  reset_transport()
  expect_error(quickgo_annotations("P21359", aspect = "anything"))
})

test_that("another aspect can be asked for", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"results":[]}')
  })

  quickgo_annotations("P21359", aspect = "molecular_function")

  expect_match(url, "aspect=molecular_function", fixed = TRUE)
})

test_that("the accession is stripped before it goes into the query", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"results":[]}')
  })

  quickgo_annotations("p21359!bad")

  expect_match(url, "UniProtKB%3AP21359BAD")
})

test_that("no usable accession never reaches the network", {
  reset_transport()
  expect_identical(quickgo_annotations("")$status, "no_data")
  expect_identical(quickgo_annotations(NULL)$status, "no_data")
})
