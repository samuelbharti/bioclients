# --- The ported fixtures, which must pass unchanged ---------------------------

test_that("the ported genescout fixture parses as expected", {
  body <- read_fixture("europepmc_nf1.json")
  out <- europepmc_parse_results(body)

  expect_s3_class(out, "tbl_df")
  expect_true(nrow(out) > 0)
  expect_false(is.na(out$title[1]))
  expect_identical(europepmc_parse_count(body), 2530L)
})

test_that("the ported variant-reviewer fixture parses as expected", {
  # The two apps had diverged into two clients against one endpoint, so both
  # stored responses have to pass the one parser that replaces them.
  body <- read_fixture("europepmc_braf_v600e.json")
  out <- europepmc_parse_results(body)

  expect_s3_class(out, "tbl_df")
  expect_true(nrow(out) > 0)
})

test_that("the citation fields from both clients are present", {
  # genescout carried pmid and the grounded source_id; variant-reviewer carried
  # doi and cited_by. The replacement has to keep all of them or one app loses
  # a column it was using.
  body <- read_fixture("europepmc_braf_v600e.json")
  out <- europepmc_parse_results(body)

  expect_true(all(
    c(
      "title",
      "authors",
      "year",
      "journal",
      "pmid",
      "doi",
      "cited_by",
      "source",
      "source_id",
      "source_url"
    ) %in%
      names(out)
  ))
  expect_type(out$cited_by, "double")
})

# --- Grounding ---------------------------------------------------------------

test_that("a record with a PMID is grounded on it", {
  body <- list(
    resultList = list(
      result = list(
        list(id = "36197410", source = "MED", pmid = "36197410", title = "x")
      )
    )
  )
  out <- europepmc_parse_results(body)

  expect_identical(out$source_id, "PMID:36197410")
  expect_identical(out$source_url, "https://europepmc.org/article/MED/36197410")
})

test_that("a record with no PMID falls back to source and id", {
  # A preprint or a patent has no PMID and is still citable.
  body <- list(
    resultList = list(
      result = list(
        list(id = "PPR123456", source = "PPR", title = "a preprint")
      )
    )
  )
  expect_identical(europepmc_parse_results(body)$source_id, "PPR:PPR123456")
})

test_that("no results is NULL", {
  expect_null(europepmc_parse_results(list(resultList = list(result = list()))))
})

# --- The query builder -------------------------------------------------------

test_that("terms are quoted so a multi-word term is one phrase", {
  expect_identical(europepmc_query("BRAF"), "\"BRAF\"")
  expect_identical(
    europepmc_query("BRAF", "V600E"),
    "\"BRAF\" AND \"V600E\""
  )
})

test_that("a quote inside a term cannot close the phrase early", {
  # An unescaped quote would end the phrase and turn the rest into query syntax.
  expect_identical(europepmc_query("BRAF\" OR x:\""), "\"BRAF OR x:\"")
})

test_that("blank terms are dropped rather than quoted empty", {
  expect_identical(europepmc_query("BRAF", NA, ""), "\"BRAF\"")
  expect_null(europepmc_query(""))
  expect_null(europepmc_query())
})

# --- A count of zero is an answer --------------------------------------------

test_that("the zero fixture is a real count, not a missing one", {
  body <- read_fixture("europepmc_zero.json")

  expect_identical(europepmc_parse_count(body), 0L)
  expect_false(is.na(europepmc_parse_count(body)))
})

test_that("a count of zero comes back ok, not no_data", {
  # A gene with no literature must not be indistinguishable from an outage. A
  # caller treating absence as evidence would otherwise read an outage as
  # evidence.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"hitCount":0,"resultList":{"result":[]}}')
  })

  res <- europepmc_count(europepmc_query("XYZ999"))

  expect_true(res$ok)
  expect_identical(res$status, "ok")
  expect_identical(res$data$count, 0L)
})

test_that("an absent hitCount is no_data, unlike a zero one", {
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json("{}"))

  expect_identical(europepmc_count("NF1")$status, "no_data")
})

test_that("an empty search is no_data, because there are no rows", {
  # Being asked how many and being asked for which are different questions.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"hitCount":0,"resultList":{"result":[]}}')
  })

  expect_identical(europepmc_search("XYZ999")$status, "no_data")
})

# --- The client half ---------------------------------------------------------

test_that("the count request asks for one idlist row, not a page of results", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"hitCount":7}')
  })

  europepmc_count("NF1")

  expect_match(url, "resultType=idlist", fixed = TRUE)
  expect_match(url, "pageSize=1", fixed = TRUE)
})

test_that("the search reports the total alongside the page", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json(
      '{"hitCount":2530,"resultList":{"result":[{"id":"1","source":"MED"}]}}'
    )
  })

  out <- biohttp::body_or_null(europepmc_search("NF1"))

  expect_identical(out$count, 2530L)
  expect_identical(nrow(out$results), 1L)
  expect_identical(out$query, "NF1")
})

test_that("a blank query never reaches the network", {
  reset_transport()
  expect_identical(europepmc_search("")$status, "no_data")
  expect_identical(europepmc_count(NULL)$status, "no_data")
})
