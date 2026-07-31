# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported PanelApp fixture parses as expected", {
  body <- read_fixture("panelapp_panel.json")
  out <- panelapp_parse_panel(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$symbol[1], "NF1")
  expect_identical(out$hgnc_id[1], "HGNC:7765")
  expect_match(out$source_url[1], "/panels/255/", fixed = TRUE)
})

# --- The confidence level is returned, not converted -------------------------

test_that("every gene comes back, red included", {
  # genescout keeps green and amber and drops red. Deciding red is not worth
  # seeing is the caller's call, so all four genes are returned here.
  body <- read_fixture("panelapp_panel.json")
  out <- panelapp_parse_panel(body)

  expect_identical(nrow(out), 4L)
  expect_true("LZTR1" %in% out$symbol)
})

test_that("the traffic-light level is named, not weighted", {
  # genescout maps green to 1.0 and amber to 0.5. That is a weight, and a
  # weight is scoring.
  body <- read_fixture("panelapp_panel.json")
  out <- panelapp_parse_panel(body)

  expect_identical(out$confidence, c("3", "3", "2", "1"))
  expect_identical(out$level, c("green", "green", "amber", "red"))
  expect_type(out$level, "character")
})

test_that("an unknown confidence level is NA rather than a guess", {
  body <- list(
    id = 1,
    genes = list(
      list(entity_name = "A", confidence_level = "0"),
      list(entity_name = "B", confidence_level = "3")
    )
  )
  out <- panelapp_parse_panel(body)

  expect_identical(out$level, c(NA_character_, "green"))
})

# --- Entity naming -----------------------------------------------------------

test_that("an entry with no gene_data falls back to entity_name", {
  # Panels carry regions and STRs alongside genes, and those have no gene_data.
  body <- list(
    id = 1,
    genes = list(
      list(entity_name = "ISCA-37421", confidence_level = "3"),
      list(
        entity_name = "NF1",
        gene_data = list(gene_symbol = "NF1"),
        confidence_level = "3"
      )
    )
  )
  out <- panelapp_parse_panel(body)

  expect_identical(out$symbol, c("ISCA-37421", "NF1"))
})

test_that("a panel with no genes is NULL", {
  expect_null(panelapp_parse_panel(list(id = 1, genes = list())))
  expect_null(panelapp_parse_panel(list()))
})

# --- The index ---------------------------------------------------------------

test_that("an index page parses, and a bare list of panels does too", {
  page <- list(
    results = list(
      list(
        id = 255,
        name = "Neurofibromatosis Type 1",
        version = "1.12",
        relevant_disorders = list("Neurofibromatosis type 1", "NF1")
      )
    )
  )
  out <- panelapp_parse_index(page)

  expect_identical(out$id, "255")
  expect_identical(out$disorders[[1]], c("Neurofibromatosis type 1", "NF1"))
  expect_identical(panelapp_parse_index(page$results)$id, "255")
})

test_that("a panel with no relevant disorders gets an empty vector", {
  out <- panelapp_parse_index(list(results = list(list(id = 1, name = "X"))))

  expect_identical(out$disorders[[1]], character())
})

# --- Trap: the index search parameter does not filter ------------------------

test_that("no search parameter is sent to the index", {
  # panels/ accepts `search` and answers 200 to it, but returns the unfiltered
  # index. A client that sent a disease name there and read the first result
  # would get whichever panel sorts first, for every query.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"results":[{"id":1,"name":"X"}],"next":null}')
  })

  panelapp_panels()

  expect_false(grepl("search", url, fixed = TRUE))
  expect_match(url, "page=1", fixed = TRUE)
})

test_that("the walk follows next and stops when it is absent", {
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    if (grepl("page=1", req$url, fixed = TRUE)) {
      return(mock_json(
        '{"results":[{"id":1,"name":"A"}],"next":"https://x/?page=2"}'
      ))
    }
    mock_json('{"results":[{"id":2,"name":"B"}],"next":null}')
  })

  out <- biohttp::body_or_null(panelapp_all_panels())

  expect_length(urls, 2)
  expect_identical(out$id, c("1", "2"))
})

test_that("the walk is capped even when next never runs out", {
  reset_transport()
  n <- 0L
  httr2::local_mocked_responses(function(req) {
    n <<- n + 1L
    mock_json('{"results":[{"id":1,"name":"A"}],"next":"https://x/?page=99"}')
  })

  panelapp_all_panels(max_pages = 3)

  expect_identical(n, 3L)
})

test_that("a failure on the first page is returned, not swallowed", {
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json("boom", status = 500L))

  expect_false(isTRUE(panelapp_all_panels()$ok))
})

# --- The detail client -------------------------------------------------------

test_that("a non-numeric panel id is refused rather than pasted into the path", {
  reset_transport()
  res <- panelapp_panel("255/../../admin")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "numeric")
})

test_that("the panel path keeps its trailing slash", {
  # PanelApp redirects a path without one, which costs a round trip on every
  # lookup.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"id":255,"genes":[]}')
  })

  panelapp_panel(255)

  expect_match(url, "/panels/255/", fixed = TRUE)
})

test_that("a blank id never reaches the network", {
  reset_transport()
  expect_identical(panelapp_panel("")$status, "no_data")
  expect_identical(panelapp_panel(NULL)$status, "no_data")
})
