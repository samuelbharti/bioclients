# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported PDBe fixture parses to distinct structures", {
  body <- read_fixture("pdbe_nf1.json")
  out <- pdbe_parse_structures(body, "P21359")

  expect_s3_class(out, "tbl_df")
  expect_identical(out$pdb_id[1], "7pgp")
  expect_identical(out$method[1], "Electron Microscopy")
  expect_equal(out$resolution[1], 3.1)
  expect_identical(
    out$source_url[1],
    "https://www.ebi.ac.uk/pdbe/entry/pdb/7pgp"
  )
})

# --- One row per structure, not per chain ------------------------------------

test_that("a structure is counted once however many chains it resolved", {
  # The response is per-chain, so one PDB entry recurs once per chain. Counting
  # raw records overstates structural coverage several-fold: the fixture has
  # 7pgp on both chain F and chain N.
  body <- read_fixture("pdbe_nf1.json")
  out <- pdbe_parse_structures(body, "P21359")

  expect_false(any(duplicated(out$pdb_id)))
  # More raw records than distinct structures, which is the point.
  expect_true(length(body[["P21359"]]) > nrow(out))
})

test_that("the first occurrence is kept, which is PDBe's best chain", {
  body <- list(
    P00000 = list(
      list(pdb_id = "1abc", chain_id = "A", coverage = 1.0, resolution = 2.0),
      list(pdb_id = "1abc", chain_id = "B", coverage = 0.4, resolution = 2.0)
    )
  )
  out <- pdbe_parse_structures(body, "P00000")

  expect_identical(nrow(out), 1L)
  expect_equal(out$coverage, 1.0)
})

# --- Key handling ------------------------------------------------------------

test_that("the body is read even when keyed by a different accession", {
  # PDBe keys the object by the accession it resolved, which is not always the
  # exact string asked for.
  body <- list(P21359 = list(list(pdb_id = "7pgp", coverage = 1)))
  expect_identical(pdbe_parse_structures(body, "P99999")$pdb_id, "7pgp")
})

test_that("no structures parses to NULL", {
  expect_null(pdbe_parse_structures(list(), "P00000"))
  expect_null(pdbe_parse_structures(list(P00000 = list()), "P00000"))
})

# --- The client half ---------------------------------------------------------

test_that("pdbe_structures returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(testthat::test_path("fixtures", "pdbe_nf1.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- pdbe_structures("P21359")

  expect_true(res$ok)
  expect_identical(res$source, "PDBe")
  expect_identical(biohttp::body_or_null(res)$pdb_id[1], "7pgp")
})

test_that("a blank accession is no_data and never reaches the network", {
  reset_transport()
  expect_identical(pdbe_structures("")$status, "no_data")
})
