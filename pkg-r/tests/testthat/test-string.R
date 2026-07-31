# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported interaction-partners fixture parses as expected", {
  body <- read_fixture("string_tp53.json")
  out <- string_parse_partners(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$partner[1], "SFN")
  expect_equal(out$score[1], 0.999)
  expect_equal(out$experimental[1], 0.981)
  expect_equal(out$database[1], 0.75)
  expect_equal(out$textmining[1], 0.859)
  expect_equal(out$coexpression[1], 0)
})

test_that("partners come back strongest first", {
  body <- read_fixture("string_tp53.json")
  out <- string_parse_partners(body)
  expect_true(all(diff(out$score) <= 0))
})

test_that("the ported network fixture parses into edges", {
  body <- read_fixture("string_network.json")
  out <- string_parse_network(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 7L)
  # Edges run between any two members of the queried set, not all from one hub.
  expect_setequal(
    unique(c(out$gene_a, out$gene_b)),
    c("TP53", "SUZ12", "NF1", "EGFR", "CDKN2A")
  )
  expect_equal(
    out$score[out$gene_a == "TP53" & out$gene_b == "CDKN2A"],
    0.999
  )
  expect_equal(
    out$score[out$gene_a == "NF1" & out$gene_b == "CDKN2A"],
    0.772
  )
})

test_that("the ported identifier-map fixture parses", {
  body <- read_fixture("string_ids.json")
  out <- string_parse_ids(body)

  expect_identical(nrow(out), 2L)
  expect_true("SEPTIN9" %in% out$query)
  expect_identical(out$preferred[out$query == "SEPTIN9"], "SEPT9")
})

# --- The renamed-gene trap ---------------------------------------------------

test_that("edges are translated back into the queried symbol space", {
  # STRING's preferredName lags HGNC. Query SEPTIN9 and the edges name SEPT9,
  # and the network endpoint does not echo the query term. A caller matching
  # edges against what it asked for finds nothing and records a real isolate.
  # The fixture pair exists precisely because SEPTIN9 is that case.
  edges <- string_parse_network(list(list(
    preferredName_A = "SEPT9",
    preferredName_B = "TP53",
    score = 0.9
  )))
  id_map <- string_parse_ids(read_fixture("string_ids.json"))

  reconciled <- string_reconcile_edges(edges, id_map)

  expect_identical(reconciled$gene_a, "SEPTIN9")
  expect_identical(reconciled$gene_b, "TP53")
})

test_that("reconciliation never makes a result worse", {
  edges <- tibble::tibble(gene_a = "A", gene_b = "B", score = 0.9)
  empty_map <- string_parse_ids(list())

  # No map, empty map, and an unmapped endpoint all leave the edges alone.
  expect_identical(string_reconcile_edges(edges, NULL), edges)
  expect_identical(string_reconcile_edges(edges, empty_map), edges)
  expect_identical(
    string_reconcile_edges(
      edges,
      tibble::tibble(
        query = "X",
        preferred = "Y",
        string_id = "1"
      )
    )$gene_a,
    "A"
  )
})

# --- Empty results -----------------------------------------------------------

test_that("no edges is a zero-row tibble, not NULL", {
  # "This set has no high-confidence edges" is a real answer about the set. NULL
  # would make it indistinguishable from a failed call.
  out <- string_parse_network(list())
  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 0L)
})

test_that("no partners parses to NULL", {
  expect_null(string_parse_partners(list()))
})

# --- The client half ---------------------------------------------------------

test_that("STRING's text/json content type is accepted", {
  # STRING serves JSON as text/json, not application/json. A client that trusts
  # the header rejects a perfectly good body. biohttp parses with
  # check_type = FALSE, and this is what pins that.
  reset_transport()
  fixture <- paste(
    readLines(
      testthat::test_path("fixtures", "string_tp53.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) {
    httr2::response(
      status_code = 200,
      headers = list(`content-type` = "text/json"),
      body = charToRaw(fixture)
    )
  })

  res <- string_partners("TP53")

  expect_true(res$ok)
  expect_identical(biohttp::body_or_null(res)$partner[1], "SFN")
})

test_that("a network reports what was actually queried", {
  # So a caller can tell a gene measured to have no partners from one that was
  # never sent. Reporting the second as the first invents a negative result.
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json("[]"))

  res <- string_network(c("TP53", "NF1"), reconcile = FALSE)
  out <- biohttp::body_or_null(res)

  expect_identical(out$queried, c("TP53", "NF1"))
  expect_identical(out$n_query, 2L)
  expect_false(out$truncated)
  expect_identical(out$n_dropped, 0L)
  expect_identical(nrow(out$edges), 0L)
})

test_that("an oversized set is capped and says so", {
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json("[]"))
  symbols <- paste0("GENE", seq_len(STRING_MAX_NODES + 25))

  out <- biohttp::body_or_null(string_network(symbols, reconcile = FALSE))

  expect_identical(out$n_query, STRING_MAX_NODES)
  expect_true(out$truncated)
  expect_identical(out$n_dropped, 25L)
})

test_that("reconciliation is skipped when there are no edges", {
  # It costs a second request, and there is nothing to translate.
  reset_transport()
  calls <- 0L
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    mock_json("[]")
  })

  string_network(c("TP53", "NF1"), reconcile = TRUE)

  expect_identical(calls, 1L)
})

test_that("fewer than two genes is no_data and never reaches the network", {
  reset_transport()
  expect_identical(string_network("TP53")$status, "no_data")
  expect_identical(string_network(character())$status, "no_data")
})

test_that("a blank symbol is no_data for partners too", {
  reset_transport()
  expect_identical(string_partners("")$status, "no_data")
})
