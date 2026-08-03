# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported Pharos fixture parses to the TDL", {
  body <- read_fixture("pharos_nf1.json")
  out <- pharos_parse_targets(body, "NF1")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$symbol, "NF1")
  expect_identical(out$tdl, "Tbio")
  expect_identical(out$source_url, "https://pharos.nih.gov/targets/NF1")
})

# --- Scoring is not this package's job ---------------------------------------

test_that("no numeric score column is returned", {
  # Both source clients carried a TDL -> 0-1 map. Turning a category into a
  # weight is scoring, and scoring belongs to the consuming app. Adding a score
  # column here would put that ranking model in two places.
  body <- read_fixture("pharos_nf1.json")
  out <- pharos_parse_targets(body, "NF1")

  expect_false("score" %in% names(out))
  expect_identical(names(out), c("symbol", "tdl", "source_url"))
})

test_that("the recognized levels are the four IDG ones", {
  expect_identical(PHAROS_TDL_LEVELS, c("Tclin", "Tchem", "Tbio", "Tdark"))
})

# --- Unknown values ----------------------------------------------------------

test_that("a TDL Pharos does not recognize becomes NA", {
  # Passing an unexpected string through would let a caller branch on something
  # this package has never seen and cannot describe.
  body <- list(
    data = list(
      targets = list(
        targets = list(
          list(sym = "NF1", tdl = "Tsomethingnew")
        )
      )
    )
  )
  expect_true(is.na(pharos_parse_targets(body, "NF1")$tdl))
})

test_that("a gene Pharos does not return is NA with no url", {
  body <- list(data = list(targets = list(targets = list())))
  out <- pharos_parse_targets(body, "NOSUCHGENE")

  expect_identical(nrow(out), 1L)
  expect_true(is.na(out$tdl))
  expect_true(is.na(out$source_url))
})

# --- Order and alignment -----------------------------------------------------

test_that("rows come back in input order regardless of response order", {
  body <- list(
    data = list(
      targets = list(
        targets = list(
          list(sym = "EGFR", tdl = "Tclin"),
          list(sym = "NF1", tdl = "Tbio")
        )
      )
    )
  )
  out <- pharos_parse_targets(body, c("NF1", "EGFR"))

  expect_identical(out$symbol, c("NF1", "EGFR"))
  expect_identical(out$tdl, c("Tbio", "Tclin"))
})

# --- The client half ---------------------------------------------------------

test_that("pharos_targets sends the whole list in one request", {
  reset_transport()
  calls <- 0L
  sent <- NULL
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    sent <<- req$body$data$variables$syms
    mock_json('{"data":{"targets":{"targets":[]}}}')
  })

  pharos_targets(c("NF1", "EGFR", "BRAF"))

  expect_identical(calls, 1L)
  expect_length(sent, 3)
})

test_that("pharos_target returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(testthat::test_path("fixtures", "pharos_nf1.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- pharos_target("NF1")

  expect_true(res$ok)
  expect_identical(res$source, "Pharos")
  expect_identical(biohttp::body_or_null(res)$tdl, "Tbio")
})

test_that("a blank symbol is no_data and never reaches the network", {
  reset_transport()
  expect_identical(pharos_target("")$status, "no_data")
})
