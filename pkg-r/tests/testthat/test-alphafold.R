# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported AlphaFold fixture parses to the model metadata", {
  body <- read_fixture("alphafold_p15056.json")
  out <- alphafold_parse_model(body, "P15056")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$accession, "P15056")
  expect_match(out$pdb_url, "AF-P15056-F1")
  expect_equal(out$version, 6)
  expect_equal(out$mean_plddt, 66.38)
  expect_identical(out$source_url, "https://alphafold.ebi.ac.uk/entry/P15056")
})

# --- The one-element array ---------------------------------------------------

test_that("the response is read as an array, which is what it is", {
  # AlphaFold returns a JSON array of model records even for a single
  # accession. Reading it as an object yields NULL for every field while still
  # looking like a successful parse.
  body <- read_fixture("alphafold_p15056.json")
  expect_true(is.null(names(body)))
  expect_false(is.na(alphafold_parse_model(body, "P15056")$pdb_url))
})

test_that("an already-unwrapped record still parses", {
  body <- list(pdbUrl = "https://example.org/x.pdb", latestVersion = 4)
  out <- alphafold_parse_model(body, "P00000")

  expect_identical(out$pdb_url, "https://example.org/x.pdb")
  expect_equal(out$version, 4)
})

test_that("no model parses to NULL", {
  expect_null(alphafold_parse_model(list(), "P00000"))
  expect_null(alphafold_parse_model(list(list(latestVersion = 1)), "P00000"))
})

# --- The client half ---------------------------------------------------------

test_that("alphafold_model returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(
      testthat::test_path("fixtures", "alphafold_p15056.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- alphafold_model("P15056")

  expect_true(res$ok)
  expect_identical(res$source, "AlphaFold")
  expect_equal(biohttp::body_or_null(res)$mean_plddt, 66.38)
})

test_that("the metadata lookup downloads no coordinates", {
  # A caller asking whether a model exists must not pull a multi-megabyte
  # structure file to find out. One request, to the prediction endpoint.
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('[{"pdbUrl":"https://example.org/x.pdb","latestVersion":1}]')
  })

  alphafold_model("P15056")

  expect_length(urls, 1)
  expect_match(urls[1], "api/prediction")
  expect_false(grepl("\\.pdb$", urls[1]))
})

test_that("an accession is sanitized before it reaches the url path", {
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('[{"pdbUrl":"https://example.org/x.pdb"}]')
  })

  alphafold_model("p15056/../../etc")

  expect_match(urls[1], "/P15056ETC$")
})

test_that("a blank accession is no_data and never reaches the network", {
  reset_transport()
  expect_identical(alphafold_model("")$status, "no_data")
  expect_identical(alphafold_model("///")$status, "no_data")
})
