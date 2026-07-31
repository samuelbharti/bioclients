# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported ClinGen validity CSV parses as expected", {
  out <- clingen_parse_validity(read_fixture_text("clingen_gene_validity.csv"))

  expect_s3_class(out, "tbl_df")
  expect_identical(out$gene[1], "NF1")
  expect_identical(out$hgnc[1], "HGNC:7765")
  expect_identical(out$disease[1], "familial ovarian cancer")
  expect_identical(out$mondo[1], "MONDO:0016248")
  expect_identical(out$moi[1], "AD")
})

test_that("a gene curated against several diseases gets a row for each", {
  out <- clingen_parse_validity(read_fixture_text("clingen_gene_validity.csv"))
  nf1 <- out[out$gene == "NF1", ]

  expect_true(nrow(nf1) > 1)
  expect_true("neurofibromatosis type 1" %in% nf1$disease)
})

# --- Trap: the header is not the first row -----------------------------------

test_that("the banner and both separator rows are skipped", {
  # The file opens with a three-line banner, a row of "++++", the real header,
  # ANOTHER row of "++++", and only then the data. An ordinary read gives a
  # table named after the banner whose first two rows are plus signs.
  lines <- readLines(
    testthat::test_path("fixtures", "clingen_gene_validity.csv"),
    warn = FALSE
  )
  expect_match(lines[1], "CLINGEN GENE DISEASE VALIDITY CURATIONS")
  expect_match(lines[4], "^\"\\+")
  expect_match(lines[5], "GENE SYMBOL")
  expect_match(lines[6], "^\"\\+")

  out <- clingen_parse_validity(paste(lines, collapse = "\n"))

  expect_true("gene" %in% names(out))
  expect_false(any(grepl("^\\++$", out$gene)))
})

test_that("a file with no header row is NULL", {
  expect_null(clingen_parse_validity("\"BANNER\",\"\"\n\"+++\",\"+++\""))
})

test_that("a header with no data rows is NULL", {
  text <- paste(
    '"GENE SYMBOL","CLASSIFICATION"',
    '"+++","+++"',
    sep = "\n"
  )
  expect_null(clingen_parse_validity(text))
})

test_that("empty input is NULL, not an error", {
  expect_null(clingen_parse_validity(""))
  expect_null(clingen_parse_validity(NULL))
})

test_that("a missing column becomes NA rather than failing the parse", {
  # ClinGen has changed this file's columns before, and losing one should cost
  # that column rather than the whole table.
  text <- paste(
    '"GENE SYMBOL","CLASSIFICATION"',
    '"NF1","Definitive"',
    sep = "\n"
  )
  out <- clingen_parse_validity(text)

  expect_identical(out$gene, "NF1")
  expect_true(is.na(out$mondo))
})

# --- The classification is not converted to a score --------------------------

test_that("the classification is the string ClinGen wrote", {
  # genescout maps Definitive to 4 and Refuted to 0. That is a weight, and a
  # weight is scoring.
  out <- clingen_parse_validity(read_fixture_text("clingen_gene_validity.csv"))

  expect_type(out$classification, "character")
  expect_true("Definitive" %in% out$classification)
  expect_true("No Known Disease Relationship" %in% out$classification)
})

test_that("a negative curation is kept, not filtered out", {
  # "No Known Disease Relationship" is evidence, and dropping it would leave a
  # caller unable to tell it from a gene ClinGen never looked at.
  out <- clingen_parse_validity(read_fixture_text("clingen_gene_validity.csv"))

  expect_true(any(out$classification == "No Known Disease Relationship"))
})

# --- Filtering ---------------------------------------------------------------

test_that("filtering by symbol is case-insensitive", {
  table <- clingen_parse_validity(
    read_fixture_text("clingen_gene_validity.csv")
  )

  expect_identical(
    nrow(clingen_validity_for(table, "nf1")),
    nrow(clingen_validity_for(table, "NF1"))
  )
})

test_that("several genes can be asked for at once", {
  table <- tibble::tibble(gene = c("NF1", "TP53", "BRCA1"))
  out <- clingen_validity_for(table, c("NF1", "BRCA1"))

  expect_identical(out$gene, c("NF1", "BRCA1"))
})

test_that("no match is NULL", {
  table <- tibble::tibble(gene = c("NF1", "TP53"))

  expect_null(clingen_validity_for(table, "XYZ999"))
  expect_null(clingen_validity_for(table, ""))
  expect_null(clingen_validity_for(NULL, "NF1"))
})

# --- The client half ---------------------------------------------------------

test_that("the whole file is fetched once and served from cache after", {
  # There is no per-gene endpoint, so the download has to pay for itself across
  # every later lookup.
  reset_transport()
  calls <- 0L
  body <- read_fixture_text("clingen_gene_validity.csv")
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    httr2::response(status_code = 200L, body = charToRaw(body))
  })

  first <- clingen_gene_validity()
  second <- clingen_gene_validity()

  expect_true(first$ok)
  expect_identical(calls, 1L)
  expect_identical(nrow(first$data), nrow(second$data))
})

test_that("the download path is the published one", {
  reset_transport()
  url <- NULL
  body <- read_fixture_text("clingen_gene_validity.csv")
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    httr2::response(status_code = 200L, body = charToRaw(body))
  })

  clingen_gene_validity()

  expect_match(url, "kb/gene-validity/download", fixed = TRUE)
})

test_that("a file with no curations is no_data", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 200L, body = charToRaw("\"BANNER\",\"\""))
  })

  expect_identical(clingen_gene_validity()$status, "no_data")
})
