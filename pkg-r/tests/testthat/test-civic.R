# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported CIViC fixture parses to the evidence counts", {
  body <- read_fixture("civic_nf1.json")
  out <- civic_parse_gene(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$id, "3867")
  expect_identical(out$symbol, "NF1")
  expect_identical(out$entrez_id, "4763")
  expect_equal(out$evidence_items, 56)
  expect_equal(out$assertions, 0)
  expect_equal(out$variants, 37)
  expect_identical(out$source_url, "https://civicdb.org/features/3867")
})

test_that("a gene CIViC does not track parses to NULL", {
  # CIViC reports it as data.gene = null. That is an absence of curation, and it
  # must not be reported as zero curated evidence, which is a claim about a gene
  # CIViC has actually looked at.
  expect_null(civic_parse_gene(list(data = list(gene = NULL))))
})

test_that("a tracked gene with nothing curated is a real zero", {
  body <- list(
    data = list(
      gene = list(
        id = 1,
        name = "X",
        entrezId = 2,
        link = "/features/1",
        stats = list(
          evidenceItemCount = 0,
          assertionCount = 0,
          variantCount = 0
        )
      )
    )
  )
  out <- civic_parse_gene(body)

  expect_equal(out$evidence_items, 0)
  expect_false(is.na(out$evidence_items))
})

test_that("a missing link falls back to the site root", {
  body <- list(data = list(gene = list(id = 1, name = "X")))
  expect_identical(civic_parse_gene(body)$source_url, "https://civicdb.org")
})

# --- The interpolated query --------------------------------------------------

test_that("a crafted symbol cannot escape the GraphQL string", {
  # CIViC takes the symbol inline rather than as a GraphQL variable, so the
  # symbol is substituted into the query text. clean_symbol() strips everything
  # outside [A-Za-z0-9._-], which is what stops a quote closing the string and
  # appending another query.
  reset_transport()
  sent <- NULL
  httr2::local_mocked_responses(function(req) {
    sent <<- req$body$data$query
    mock_json('{"data":{"gene":null}}')
  })

  civic_gene('NF1") { id } evil(x: "')

  expect_false(grepl("evil", sent, fixed = TRUE))
  expect_match(sent, 'entrezSymbol: "NF1IDEVILX"', fixed = TRUE)
  # Exactly two double quotes in the whole query: the pair around the symbol.
  # Any more means the payload opened a string of its own.
  expect_identical(lengths(regmatches(sent, gregexpr('"', sent))), 2L)
})

test_that("a symbol with only illegal characters never reaches the network", {
  reset_transport()
  expect_identical(civic_gene('"; }')$status, "no_data")
})

# --- The client half ---------------------------------------------------------

test_that("civic_gene returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(testthat::test_path("fixtures", "civic_nf1.json"), warn = FALSE),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- civic_gene("NF1")

  expect_true(res$ok)
  expect_identical(res$source, "CIViC")
  expect_equal(biohttp::body_or_null(res)$evidence_items, 56)
})

test_that("an untracked gene is no_data rather than an error", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"data":{"gene":null}}')
  })
  expect_identical(civic_gene("NOSUCHGENE")$status, "no_data")
})

test_that("a GraphQL errors array inside a 200 is a failure", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"errors":[{"message":"bad"}]}')
  })
  expect_identical(civic_gene("NF1")$status, "error")
})
