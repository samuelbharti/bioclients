# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported HPO fixture parses as expected", {
  body <- read_fixture("hpo_tp53.json")
  out <- hpo_parse_diseases(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 4L)
  expect_identical(out$id[1], "OMIM:151623")
  expect_identical(out$name[1], "Li-Fraumeni syndrome")
  expect_identical(out$mondo[1], "MONDO:0018875")
})

test_that("the same body carries phenotypes as a separate table", {
  body <- read_fixture("hpo_tp53.json")
  out <- hpo_parse_phenotypes(body)

  expect_identical(nrow(out), 2L)
  expect_identical(out$id, c("HP:0003002", "HP:0006744"))
})

test_that("an explicit null description reads as NA, not an error", {
  # Every disease in the fixture sends description: null, so a parser that
  # assumed prose here would fail on the common case rather than the rare one.
  body <- read_fixture("hpo_tp53.json")
  expect_true(all(is.na(hpo_parse_diseases(body)$description)))
})

test_that("both ids are kept, not just OMIM", {
  # HPO grounds a disease in OMIM, ORPHA, or DECIPHER, and the fixture has both
  # of the first two.
  body <- read_fixture("hpo_tp53.json")
  expect_true(any(grepl("^ORPHA:", hpo_parse_diseases(body)$id)))
})

# --- Empty shapes ------------------------------------------------------------

test_that("an absent or empty array is NULL, not a zero-row tibble", {
  expect_null(hpo_parse_diseases(list()))
  expect_null(hpo_parse_diseases(list(diseases = list())))
  expect_null(hpo_parse_phenotypes(list(phenotypes = list())))
})

test_that("a gene with diseases but no phenotypes still parses", {
  body <- list(diseases = list(list(id = "OMIM:1", name = "X")))
  expect_identical(nrow(hpo_parse_diseases(body)), 1L)
  expect_null(hpo_parse_phenotypes(body))
})

# --- The lookup key ----------------------------------------------------------

test_that("the request is keyed by NCBI Gene id", {
  # The path is network/annotation/NCBIGene:<id>. Anything but the numeric
  # Entrez id there is a request for an entity HPO does not have.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"diseases":[]}')
  })

  hpo_gene_annotation(7157)

  expect_match(url, "network/annotation/NCBIGene:7157", fixed = TRUE)
})

test_that("a symbol is refused rather than sent as an id", {
  # TP53 in the id slot would be a live request that comes back empty, which
  # reads as "HPO has nothing for this gene".
  reset_transport()
  res <- hpo_gene_annotation("TP53")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "NCBI Gene id")
})

test_that("a blank id never reaches the network", {
  reset_transport()
  expect_identical(hpo_gene_annotation("")$status, "no_data")
  expect_identical(hpo_gene_annotation(NULL)$status, "no_data")
})

# --- The client half ---------------------------------------------------------

test_that("both tables come back on one call", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json(paste0(
      '{"diseases":[{"id":"OMIM:151623","name":"Li-Fraumeni syndrome",',
      '"mondoId":"MONDO:0018875"}],',
      '"phenotypes":[{"id":"HP:0003002","name":"Breast carcinoma"}]}'
    ))
  })

  out <- biohttp::body_or_null(hpo_gene_annotation(7157))

  expect_identical(nrow(out$diseases), 1L)
  expect_identical(nrow(out$phenotypes), 1L)
  expect_match(out$source_url, "NCBIGene:7157", fixed = TRUE)
})

test_that("a gene HPO has no annotation for is no_data", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"diseases":[],"phenotypes":[]}')
  })

  expect_identical(hpo_gene_annotation(7157)$status, "no_data")
})

# --- Term search -------------------------------------------------------------

test_that("the ported HPO search fixture parses as expected", {
  body <- read_fixture("hpo_search_seizure.json")
  out <- hpo_parse_search(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$id[1], "HP:0007207")
  expect_identical(out$name[1], "Photosensitive tonic-clonic seizure")
  expect_match(out$definition[1], "flashing or flickering light")
})

test_that("the descendant count comes through", {
  # It is what tells a broad term from a specific one without walking the DAG.
  body <- read_fixture("hpo_search_seizure.json")
  expect_equal(hpo_parse_search(body)$descendant_count[1], 0)
})

test_that("no match is NULL, not a zero-row tibble", {
  expect_null(hpo_parse_search(list(terms = list())))
  expect_null(hpo_parse_search(list()))
})

test_that("the search sends the text and the limit", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"terms":[]}')
  })

  hpo_search("seizure", limit = 5)

  expect_match(url, "hp/search", fixed = TRUE)
  expect_match(url, "q=seizure", fixed = TRUE)
  expect_match(url, "limit=5", fixed = TRUE)
})

# --- Term lookup -------------------------------------------------------------

test_that("the ported HPO term fixture parses as expected", {
  body <- read_fixture("hpo_term_HP0001250.json")
  out <- hpo_parse_term(body)

  expect_identical(nrow(out), 1L)
  expect_identical(out$id, "HP:0001250")
  expect_identical(out$name, "Seizure")
  expect_equal(out$descendant_count, 346)
})

test_that("synonyms and xrefs are list columns", {
  # A term carries any number of each, and flattening them to one string would
  # make them unusable without re-splitting.
  body <- read_fixture("hpo_term_HP0001250.json")
  out <- hpo_parse_term(body)

  expect_type(out$synonyms, "list")
  expect_true("Epilepsy" %in% out$synonyms[[1]])
  expect_true(any(grepl("^SNOMEDCT_US:", out$xrefs[[1]])))
})

test_that("a term with no synonyms gets an empty vector, not NULL", {
  out <- hpo_parse_term(list(id = "HP:1", name = "X"))

  expect_identical(out$synonyms[[1]], character())
  expect_identical(out$xrefs[[1]], character())
})

test_that("a body with no id is NULL", {
  expect_null(hpo_parse_term(list(name = "Seizure")))
  expect_null(hpo_parse_term(list()))
})

test_that("the term id is upper-cased into the path, colon and all", {
  # JAX takes an HP id in the path directly, with the colon verbatim.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"id":"HP:0001250","name":"Seizure"}')
  })

  hpo_term("hp:0001250")

  expect_match(url, "hp/terms/HP:0001250", fixed = TRUE)
})

test_that("anything that is not an HP id never reaches the network", {
  reset_transport()
  expect_identical(hpo_term("seizure")$status, "no_data")
  expect_identical(hpo_term("HP:0001250/../../x")$status, "no_data")
  expect_identical(hpo_search("")$status, "no_data")
})
