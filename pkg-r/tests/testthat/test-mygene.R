# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported MyGene fixture parses to the expected gene", {
  body <- read_fixture("mygene_tp53.json")
  out <- mygene_parse_hits(body, "TP53")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$symbol, "TP53")
  expect_identical(out$entrez, "7157")
  expect_identical(out$ensembl_gene, "ENSG00000141510")
  expect_identical(out$uniprot, "P04637")
  expect_identical(out$name, "tumor protein p53")
  expect_identical(out$type_of_gene, "protein-coding")
  expect_match(out$summary, "tumor suppressor")
})

test_that("uniprot takes Swiss-Prot, not the first TrEMBL accession", {
  # The fixture carries a Swiss-Prot scalar and a long TrEMBL array. Reading the
  # wrong one yields a plausible accession for the wrong protein record.
  body <- read_fixture("mygene_tp53.json")
  expect_identical(mygene_parse_hits(body, "TP53")$uniprot, "P04637")
})

# --- The scoring trap, which is the point of this client ----------------------

test_that("an exact symbol match beats a higher-scored match to another gene", {
  # Querying TTN with the alias and retired scopes returns TTR (entrez 7276,
  # score 19.07) ahead of TTN (7273, score 18.29). Taking the top hit returns
  # the wrong gene, and it looks entirely plausible while doing it.
  hits <- list(
    list(symbol = "TTR", entrezgene = "7276", `_score` = 19.07),
    list(symbol = "TTN", entrezgene = "7273", `_score` = 18.29)
  )
  expect_identical(mygene_pick_hit(hits, "TTN")$entrezgene, "7273")
})

test_that("a deliberate alias still falls back to the best-scored hit", {
  # "p53" is nobody's official symbol, so there is no exact match to prefer and
  # MyGene's own ranking is the right answer.
  hits <- list(
    list(symbol = "TP53", entrezgene = "7157"),
    list(symbol = "TP53BP1", entrezgene = "7158")
  )
  expect_identical(mygene_pick_hit(hits, "p53")$entrezgene, "7157")
})

test_that("picking from no hits is NULL rather than an error", {
  expect_null(mygene_pick_hit(list(), "TP53"))
  expect_null(mygene_pick_hit(NULL, "TP53"))
})

# --- The batch parser --------------------------------------------------------

test_that("batch results come back in input order, one row per input", {
  # MyGene returns the flat array in whatever order it likes, so a caller that
  # zips by position needs this guarantee. A dropped row would shift every row
  # after it onto the wrong gene.
  body <- list(
    list(query = "EGFR", symbol = "EGFR", entrezgene = "1956"),
    list(query = "TP53", symbol = "TP53", entrezgene = "7157")
  )
  out <- mygene_parse_batch(body, c("TP53", "EGFR"))

  expect_identical(nrow(out), 2L)
  expect_identical(out$symbol, c("TP53", "EGFR"))
  expect_identical(out$entrez, c("7157", "1956"))
})

test_that("an unmatched token becomes an NA row, not a missing row", {
  body <- list(
    list(query = "TP53", symbol = "TP53", entrezgene = "7157"),
    list(query = "NOPE", notfound = TRUE)
  )
  out <- mygene_parse_batch(body, c("TP53", "NOPE"))

  expect_identical(nrow(out), 2L)
  expect_identical(out$symbol, c("TP53", "NOPE"))
  expect_true(is.na(out$entrez[2]))
})

test_that("an unusable token still occupies its row", {
  out <- mygene_parse_batch(list(), c("TP53", "   "))
  expect_identical(nrow(out), 2L)
  expect_true(is.na(out$entrez[1]))
})

test_that("the batch parser applies the exact-symbol rule too", {
  body <- list(
    list(query = "TTN", symbol = "TTR", entrezgene = "7276"),
    list(query = "TTN", symbol = "TTN", entrezgene = "7273")
  )
  expect_identical(mygene_parse_batch(body, "TTN")$entrez, "7273")
})

# --- The client half ---------------------------------------------------------

test_that("mygene_gene returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- readLines(
    testthat::test_path("fixtures", "mygene_tp53.json"),
    warn = FALSE
  )
  httr2::local_mocked_responses(function(req) {
    mock_json(paste(fixture, collapse = ""))
  })

  res <- mygene_gene("TP53")

  expect_true(res$ok)
  expect_identical(res$status, "ok")
  expect_identical(res$source, "MyGene")
  expect_identical(biohttp::body_or_null(res)$symbol, "TP53")
})

test_that("a blank identifier is no_data and never reaches the network", {
  reset_transport()
  # No mock installed, so a dispatched request would attempt a real call.
  res <- mygene_gene("   ")
  expect_identical(res$status, "no_data")
})

test_that("an empty hit list is no_data rather than a wrong gene", {
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json('{"hits":[]}'))
  expect_identical(mygene_gene("NOSUCHGENE")$status, "no_data")
})

test_that("a transport failure passes the envelope straight through", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 503)
  })
  res <- mygene_gene("TP53")

  expect_false(res$ok)
  expect_identical(res$status, "error")
  expect_identical(res$source, "MyGene")
  reset_transport()
})

# --- The HGNC id -------------------------------------------------------------

test_that("HGNC is requested, and it is upper case", {
  # Every other field is lower case. MyGene names this one HGNC in both the
  # request and the response, so asking for `hgnc` returns nothing and reads as
  # "this gene has no HGNC id" rather than as a mistake.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"hits":[]}')
  })

  mygene_gene("TP53")

  expect_match(url, "HGNC", fixed = TRUE)
  expect_false(grepl("fields=[^&]*[^A-Z]hgnc", url))
})

test_that("the HGNC id is read off a hit", {
  # The stored response predates the field, so this is pinned against a record
  # rather than against the fixture.
  body <- list(hits = list(list(symbol = "TP53", HGNC = "11998")))
  out <- mygene_parse_hits(body, "TP53")

  expect_identical(out$hgnc, "11998")
})

test_that("the id stays in the bare form MyGene sends", {
  # Monarch wants HGNC:11998 and builds that itself. Baking one consumer's
  # formatting into the column would make every other consumer strip it again.
  body <- list(hits = list(list(symbol = "TP53", HGNC = "11998")))

  expect_false(grepl(
    "HGNC:",
    mygene_parse_hits(body, "TP53")$hgnc,
    fixed = TRUE
  ))
})

test_that("a gene with no HGNC id gets NA, not an error", {
  body <- list(hits = list(list(symbol = "TP53")))
  expect_true(is.na(mygene_parse_hits(body, "TP53")$hgnc))
})

test_that("the batch path carries the column too", {
  body <- list(
    list(query = "TP53", symbol = "TP53", HGNC = "11998"),
    list(query = "NOPE", notfound = TRUE)
  )
  out <- mygene_parse_batch(body, c("TP53", "NOPE"))

  expect_identical(out$hgnc, c("11998", NA_character_))
})
