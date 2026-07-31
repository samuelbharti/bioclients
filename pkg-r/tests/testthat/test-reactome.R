# --- The ported fixtures, which must pass unchanged ---------------------------

test_that("the ported Reactome fixture parses as expected", {
  body <- read_fixture("reactome_nf1.json")
  out <- reactome_parse_pathways(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$pathway_id[1], "R-HSA-5658442")
  expect_identical(out$name[1], "Regulation of RAS by GAPs")
  expect_match(out$source_url[1], "R-HSA-5658442", fixed = TRUE)
})

test_that("the second ported fixture parses too", {
  out <- reactome_parse_pathways(read_fixture("reactome_ttn.json"))

  expect_s3_class(out, "tbl_df")
  expect_true(all(grepl("^R-HSA-", out$pathway_id)))
})

# --- Trap: the response is a bare array --------------------------------------

test_that("the pathways are read from the top level", {
  # Reactome answers with a JSON array, not an object wrapping one. A parser
  # reaching for a results or pathways key finds nothing and reports every gene
  # as having no pathways.
  body <- read_fixture("reactome_nf1.json")

  expect_null(names(body))
  expect_true(length(body) > 1)
  expect_identical(nrow(reactome_parse_pathways(body)), length(body))
})

test_that("the empty fixture is NULL, not an error", {
  expect_null(reactome_parse_pathways(read_fixture("reactome_empty.json")))
  expect_null(reactome_parse_pathways(list()))
  expect_null(reactome_parse_pathways(NULL))
})

test_that("a record with no stable id is dropped", {
  body <- list(
    list(displayName = "no id here"),
    list(stId = "R-HSA-1", displayName = "real")
  )
  expect_identical(reactome_parse_pathways(body)$pathway_id, "R-HSA-1")
})

# --- The disease flag --------------------------------------------------------

test_that("isInDisease is carried through as a logical", {
  # It is Reactome's flag on the pathway. Which pathways matter to a review is
  # the review's call, so the flag is reported rather than used to filter.
  body <- list(
    list(stId = "R-HSA-1", displayName = "a", isInDisease = TRUE),
    list(stId = "R-HSA-2", displayName = "b", isInDisease = FALSE),
    list(stId = "R-HSA-3", displayName = "c")
  )
  out <- reactome_parse_pathways(body)

  expect_type(out$in_disease, "logical")
  expect_identical(out$in_disease, c(TRUE, FALSE, FALSE))
})

test_that("no pathway is filtered out of the result", {
  body <- read_fixture("reactome_nf1.json")
  expect_identical(nrow(reactome_parse_pathways(body)), length(body))
})

# --- The client half ---------------------------------------------------------

test_that("the symbol goes into the HGNC mapping path, upper case", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('[{"stId":"R-HSA-1","displayName":"x"}]')
  })

  reactome_pathways("nf1")

  expect_match(url, "data/mapping/HGNC/NF1/pathways", fixed = TRUE)
  expect_match(url, "species=9606", fixed = TRUE)
})

test_that("a 404 is the source having nothing, not the call failing", {
  # Reactome answers 404 for a gene it maps but has no pathways for.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 404L)
  })

  res <- reactome_pathways("NF1")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "no pathways")
})

test_that("no usable symbol never reaches the network", {
  reset_transport()
  expect_identical(reactome_pathways("")$status, "no_data")
  expect_identical(reactome_pathways(NULL)$status, "no_data")
})
