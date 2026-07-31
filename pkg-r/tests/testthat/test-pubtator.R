# --- The ported fixtures, which must pass unchanged ---------------------------

test_that("the ported PubTator fixture parses as expected", {
  body <- read_fixture("pubtator_search_tp53.json")
  out <- pubtator_parse_results(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$pmid[1], "36197410")
  expect_identical(out$title[1], "TP53 or Not TP53: That Is the Question.")
  expect_identical(out$journal[1], "Clin Cancer Res")
})

test_that("authors is a list column", {
  # An article has any number, and flattening them to one string would make
  # them unusable without re-splitting.
  body <- read_fixture("pubtator_search_tp53.json")
  out <- pubtator_parse_results(body)

  expect_type(out$authors, "list")
  expect_true(length(out$authors[[1]]) > 1)
})

test_that("the count comes off the ported fixture", {
  body <- read_fixture("pubtator_search_tp53.json")
  expect_false(is.na(pubtator_parse_count(body)))
})

# --- A count of zero is an answer --------------------------------------------

test_that("the zero fixture is a real count with no results", {
  body <- read_fixture("pubtator_search_zero.json")

  expect_identical(pubtator_parse_count(body), 0L)
  expect_null(pubtator_parse_results(body))
})

test_that("a count of zero comes back ok, not no_data", {
  # PubTator has searched the whole corpus either way, so 0 is a measurement.
  # Contrast impc.R, where 0 is refused because "no phenotype found" cannot be
  # told apart from "never phenotyped".
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"results":[],"count":0}')
  })

  res <- pubtator_gene_literature("XYZ999")

  expect_true(res$ok)
  expect_identical(res$data$count, 0L)
  expect_null(res$data$results)
})

test_that("an absent count is no_data, unlike a zero one", {
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json('{"results":[]}'))

  expect_identical(pubtator_gene_literature("TP53")$status, "no_data")
})

# --- The entity token --------------------------------------------------------

test_that("the Entrez id is preferred over the symbol", {
  # @GENE_7157 names TP53 exactly. @GENE_TP53 relies on PubTator resolving the
  # symbol, which is what goes wrong for an alias.
  expect_identical(pubtator_entity("TP53", 7157), "@GENE_7157")
  expect_identical(pubtator_entity("TP53", "7157"), "@GENE_7157")
})

test_that("the symbol is the fallback when there is no Entrez id", {
  expect_identical(pubtator_entity("TP53"), "@GENE_TP53")
  expect_identical(pubtator_entity("tp53"), "@GENE_TP53")
  expect_identical(pubtator_entity("TP53", NULL), "@GENE_TP53")
})

test_that("a non-numeric entrez falls back rather than being sent", {
  expect_identical(pubtator_entity("TP53", "not-an-id"), "@GENE_TP53")
})

test_that("neither argument usable is NULL", {
  expect_null(pubtator_entity())
  expect_null(pubtator_entity(""))
  expect_null(pubtator_entity(NULL, NULL))
})

# --- The client half ---------------------------------------------------------

test_that("the entity token is what gets searched", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"results":[],"count":0}')
  })

  pubtator_gene_literature("TP53", 7157)

  expect_match(url, "search/", fixed = TRUE)
  expect_match(url, "text=%40GENE_7157")
})

test_that("the reported entity matches what was queried", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"results":[],"count":12}')
  })

  out <- biohttp::body_or_null(pubtator_gene_literature("TP53", 7157))

  expect_identical(out$entity, "@GENE_7157")
  expect_identical(out$count, 12L)
})

test_that("no usable gene never reaches the network", {
  reset_transport()
  expect_identical(pubtator_gene_literature("")$status, "no_data")
  expect_identical(pubtator_gene_literature(NULL, NULL)$status, "no_data")
})
