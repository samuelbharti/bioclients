# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported DISEASES fixture parses as expected", {
  body <- read_fixture("diseases_nf1.json")
  out <- diseases_parse_channel(body, "Knowledge")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 4L)
  expect_identical(out$symbol[1], "NF1")
  expect_equal(out$score[1], 5.0)
  expect_identical(out$channel[1], "Knowledge")
})

test_that("the Ensembl protein key is kept as the row's grounding", {
  # It is the only identifier the response carries for an association.
  body <- read_fixture("diseases_nf1.json")
  out <- diseases_parse_channel(body, "Knowledge")

  expect_identical(out$protein[1], "ENSP00000351015")
})

# --- Trap: the data is one level down ----------------------------------------

test_that("the associations are read from the first array element", {
  # The body is a JSON array. Element one is the map of associations, element
  # two is an empty object. Treating the top level as the record list finds two
  # things, neither an association, and reports zero genes for every disease
  # without failing.
  body <- read_fixture("diseases_nf1.json")

  expect_length(body, 2)
  expect_length(body[[2]], 0)
  expect_identical(nrow(diseases_parse_channel(body)), 4L)
})

test_that("an empty first element is NULL", {
  expect_null(diseases_parse_channel(list(list(), list())))
  expect_null(diseases_parse_channel(list()))
})

# --- Column types ------------------------------------------------------------

test_that("the columns carry no names from the protein-keyed map", {
  # The map is keyed by ENSP, so vapply over it names every column. A named
  # numeric breaks equality against a plain one, silently, in whatever joins
  # the result later.
  body <- read_fixture("diseases_nf1.json")
  out <- diseases_parse_channel(body, "Knowledge")

  expect_null(names(out$score))
  expect_null(names(out$symbol))
  expect_identical(out$score, c(5.0, 4.2, 3.5, 2.7))
})

test_that("a row with no symbol or no usable score is dropped, quietly", {
  # A dropped row should not also emit a bare coercion warning naming neither
  # the field nor the service.
  body <- list(list(
    ENSP1 = list(name = "NF1", score = 5),
    ENSP2 = list(name = "", score = 3),
    ENSP3 = list(name = "TP53", score = NULL),
    ENSP4 = list(name = "BRCA1", score = "not a number")
  ))

  expect_silent(out <- diseases_parse_channel(body))
  expect_identical(out$symbol, "NF1")
})

# --- Merging channels --------------------------------------------------------

test_that("the strongest score per gene wins and names its channel", {
  knowledge <- diseases_parse_channel(
    list(list(E1 = list(name = "NF1", score = 5))),
    "Knowledge"
  )
  textmining <- diseases_parse_channel(
    list(list(
      E1 = list(name = "NF1", score = 2),
      E2 = list(name = "TP53", score = 4)
    )),
    "Textmining"
  )
  out <- diseases_merge_channels(list(knowledge, textmining))

  expect_identical(out$symbol, c("NF1", "TP53"))
  expect_equal(out$score, c(5, 4))
  expect_identical(out$channel, c("Knowledge", "Textmining"))
})

test_that("merging ignores NULL and empty channels", {
  one <- diseases_parse_channel(list(list(E1 = list(name = "NF1", score = 5))))

  expect_identical(nrow(diseases_merge_channels(list(one, NULL))), 1L)
  expect_null(diseases_merge_channels(list(NULL, NULL)))
  expect_null(diseases_merge_channels(list()))
})

test_that("the result is ordered by score, highest first", {
  tbl <- diseases_parse_channel(list(list(
    E1 = list(name = "A", score = 1),
    E2 = list(name = "B", score = 9),
    E3 = list(name = "C", score = 5)
  )))
  expect_identical(diseases_merge_channels(list(tbl))$symbol, c("B", "C", "A"))
})

# --- The client half ---------------------------------------------------------

test_that("the query carries the Disease Ontology and human type codes", {
  # type1 = -26 selects a Disease Ontology term and type2 = 9606 selects human.
  # Both are positional magic numbers, so getting either wrong returns a
  # perfectly well-formed answer about something else.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[{},{}]")
  })

  diseases_channel("DOID:0060293", "Knowledge")

  expect_match(url, "type1=-26", fixed = TRUE)
  expect_match(url, "type2=9606", fixed = TRUE)
  expect_match(url, "id1=DOID", fixed = TRUE)
  expect_match(url, "/Knowledge", fixed = TRUE)
})

test_that("the channel is checked rather than pasted into the path", {
  reset_transport()
  expect_error(diseases_channel("DOID:1", "Experiments"))
})

test_that("both channels are queried for a combined lookup", {
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('[{"E1":{"name":"NF1","score":5}},{}]')
  })

  out <- biohttp::body_or_null(diseases_gene_associations("DOID:0060293"))

  expect_length(urls, 2)
  expect_true(any(grepl("/Knowledge", urls, fixed = TRUE)))
  expect_true(any(grepl("/Textmining", urls, fixed = TRUE)))
  expect_identical(out$symbol, "NF1")
})

test_that("a failed channel fails the combined call", {
  # A combined score quietly missing the curated channel looks like a weaker
  # association, not a partial answer.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    if (grepl("Textmining", req$url, fixed = TRUE)) {
      return(mock_json("boom", status = 500L))
    }
    mock_json('[{"E1":{"name":"NF1","score":5}},{}]')
  })

  expect_false(isTRUE(diseases_gene_associations("DOID:0060293")$ok))
})

test_that("a blank DOID never reaches the network", {
  reset_transport()
  expect_identical(diseases_channel("")$status, "no_data")
  expect_identical(diseases_gene_associations(NULL)$status, "no_data")
})
