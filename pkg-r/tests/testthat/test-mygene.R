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

# --- Genes with several Ensembl ids -----------------------------------------

test_that("PTEN keeps its reference id when MyGene also lists a patch id", {
  # Recorded live on 2026-10-07 with the client's own request: q=PTEN,
  # species=human, size=5 and MYGENE_FIELDS. MyGene sends `ensembl` as a list
  # of two objects, and the second id sits on HG2334_PATCH. Before this was
  # fixed, the column came back NA for every such gene.
  body <- read_fixture("mygene_pten.json")
  expect_null(names(body$hits[[1]]$ensembl))

  out <- mygene_parse_hits(body, "PTEN")
  expect_identical(out$ensembl_gene, "ENSG00000171862")
})

test_that("HLA-A keeps the chromosome 6 id, which MyGene lists last", {
  # Recorded live on 2026-10-07, same request as above with q=HLA-A. HLA-A has
  # eight ids, seven of them on alternate haplotypes. The first in the list is
  # one of those, so taking the first would give an id that looks fine and that
  # Open Targets and the Human Protein Atlas do not know.
  body <- read_fixture("mygene_hla_a.json")
  ensembl <- body$hits[[1]]$ensembl
  expect_length(ensembl, 8)
  expect_false(identical(ensembl[[1]]$gene, "ENSG00000206503"))

  out <- mygene_parse_hits(body, "HLA-A")
  expect_identical(out$ensembl_gene, "ENSG00000206503")
})

test_that("a gene with an id on X and one on Y keeps the X one", {
  # SHOX sits in the region X and Y share, with a separate id on each. Y is
  # listed first here on purpose, so the test shows the pick follows the
  # chromosome and not the list order.
  body <- list(
    hits = list(list(
      symbol = "SHOX",
      ensembl = list(
        list(gene = "ENSG00000292354"),
        list(gene = "ENSG00000185960")
      ),
      genomic_pos = list(
        list(chr = "Y", ensemblgene = "ENSG00000292354"),
        list(chr = "X", ensemblgene = "ENSG00000185960")
      )
    ))
  )

  expect_identical(
    mygene_parse_hits(body, "SHOX")$ensembl_gene,
    "ENSG00000185960"
  )
})

test_that("several ids and none on a reference chromosome give NA", {
  # HLA-DRB3 exists only on alternate haplotypes. Any of its three ids would
  # be a guess, and none of them is a reference id.
  body <- list(
    hits = list(list(
      symbol = "HLA-DRB3",
      ensembl = list(
        list(gene = "ENSG00000231679"),
        list(gene = "ENSG00000230463"),
        list(gene = "ENSG00000196101")
      ),
      genomic_pos = list(
        list(chr = "HSCHR6_MHC_QBL_CTG1", ensemblgene = "ENSG00000196101"),
        list(chr = "HSCHR6_MHC_COX_CTG1", ensemblgene = "ENSG00000231679"),
        list(chr = "HSCHR6_MHC_APD_CTG1", ensemblgene = "ENSG00000230463")
      )
    ))
  )

  expect_true(is.na(mygene_parse_hits(body, "HLA-DRB3")$ensembl_gene))
})

test_that("several ids with no genomic_pos give NA, not a guess", {
  body <- list(
    hits = list(list(
      symbol = "PTEN",
      ensembl = list(
        list(gene = "ENSG00000171862"),
        list(gene = "ENSG00000284792")
      )
    ))
  )

  expect_true(is.na(mygene_parse_hits(body, "PTEN")$ensembl_gene))
})

test_that("a single id is kept wherever it sits, as before", {
  # Only a choice between several ids needs a reference chromosome. One id is
  # the only answer MyGene has, so it is kept.
  body <- list(
    hits = list(list(
      symbol = "HLA-DRB3",
      ensembl = list(gene = "ENSG00000196101"),
      genomic_pos = list(
        chr = "HSCHR6_MHC_QBL_CTG1",
        ensemblgene = "ENSG00000196101"
      )
    ))
  )

  expect_identical(
    mygene_parse_hits(body, "HLA-DRB3")$ensembl_gene,
    "ENSG00000196101"
  )
})

test_that("the batch path picks the reference id too", {
  # Both paths build their row the same way, and the issue showed up in both.
  body <- list(
    list(
      query = "PTEN",
      symbol = "PTEN",
      ensembl = list(
        list(gene = "ENSG00000171862"),
        list(gene = "ENSG00000284792")
      ),
      genomic_pos = list(
        list(chr = "10", ensemblgene = "ENSG00000171862"),
        list(chr = "HG2334_PATCH", ensemblgene = "ENSG00000284792")
      )
    ),
    list(
      query = "NF1",
      symbol = "NF1",
      ensembl = list(gene = "ENSG00000196712")
    )
  )
  out <- mygene_parse_batch(body, c("PTEN", "NF1"))

  expect_identical(out$ensembl_gene, c("ENSG00000171862", "ENSG00000196712"))
})

test_that("genomic_pos is requested, both the chromosome and the id", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"hits":[]}')
  })

  mygene_gene("PTEN")

  expect_match(url, "genomic_pos.chr", fixed = TRUE)
  expect_match(url, "genomic_pos.ensemblgene", fixed = TRUE)
})

# --- The batch client and its chunking ---------------------------------------

# A mock that answers each batch POST with one hit per query it was sent,
# echoing the query the way MyGene does. `fail_when` makes a chunk fail.
mygene_chunk_mock <- function(seen, fail_when = function(queries) FALSE) {
  function(req) {
    queries <- unlist(req$body$data$q, use.names = FALSE)
    seen$calls <- c(seen$calls, list(queries))
    if (isTRUE(fail_when(queries))) {
      return(httr2::response(status_code = 503))
    }
    hits <- lapply(queries, function(query) {
      list(query = query, symbol = query, entrezgene = paste0("id_", query))
    })
    mock_json(jsonlite::toJSON(hits, auto_unbox = TRUE))
  }
}

test_that("mygene_genes returns one row per input in input order", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(mygene_chunk_mock(seen))

  res <- mygene_genes(c("TP53", "BRCA1", "EGFR"))

  expect_true(res$ok)
  expect_length(seen$calls, 1)
  out <- biohttp::body_or_null(res)
  expect_identical(out$symbol, c("TP53", "BRCA1", "EGFR"))
  expect_identical(out$entrez, c("id_TP53", "id_BRCA1", "id_EGFR"))
})

test_that("the batch is chunked at the MyGene limit", {
  # 1001 identifiers is two requests: 1000 and 1. MyGene answers a body over
  # the limit with an error rather than a truncated result.
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(mygene_chunk_mock(seen))
  symbols <- paste0("G", seq_len(MYGENE_BATCH + 1L))

  res <- mygene_genes(symbols)

  expect_true(res$ok)
  expect_identical(lengths(seen$calls), c(1000L, 1L))
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), MYGENE_BATCH + 1L)
  expect_identical(out$symbol, symbols)
  expect_identical(out$entrez, paste0("id_", symbols))
})

test_that("the limit matches the documented one", {
  expect_identical(MYGENE_BATCH, 1000L)
})

test_that("chunks merge in input order even when a later chunk answers first", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(mygene_chunk_mock(seen))

  res <- mygene_genes(c("A", "B", "C", "D", "E"), chunk_size = 2)

  expect_length(seen$calls, 3)
  expect_identical(
    biohttp::body_or_null(res)$symbol,
    c("A", "B", "C", "D", "E")
  )
})

test_that("a failed chunk becomes NA rows rather than failing the call", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(mygene_chunk_mock(
    seen,
    fail_when = function(queries) "C" %in% queries
  ))

  res <- mygene_genes(c("A", "B", "C", "D", "E"), chunk_size = 2)

  expect_true(res$ok)
  out <- biohttp::body_or_null(res)
  expect_identical(out$symbol, c("A", "B", "C", "D", "E"))
  expect_identical(out$entrez, c("id_A", "id_B", NA, NA, "id_E"))
  reset_transport()
})

test_that("when every chunk fails the envelope passes straight through", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 503)
  })

  res <- mygene_genes(c("TP53", "BRCA1"))

  expect_false(res$ok)
  expect_identical(res$status, "error")
  expect_identical(res$source, "MyGene")
  reset_transport()
})

test_that("a chunk size over the limit is refused", {
  reset_transport()
  expect_error(
    mygene_genes("TP53", chunk_size = MYGENE_BATCH + 1L),
    "chunk_size"
  )
})

test_that("duplicates and blanks are sent once and still fill every row", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  seen <- new.env()
  seen$calls <- list()
  httr2::local_mocked_responses(mygene_chunk_mock(seen))

  res <- mygene_genes(c("TP53", "  ", "TP53"))

  expect_identical(seen$calls[[1]], "TP53")
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), 3L)
  expect_identical(out$entrez, c("id_TP53", NA, "id_TP53"))
})
