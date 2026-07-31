# --- The ported fixtures, which must pass unchanged ---------------------------

test_that("the ported IMPC gene fixture resolves the mouse ortholog", {
  body <- read_fixture("impc_gene_nf1.json")
  out <- impc_parse_ortholog(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$mgi, "MGI:97306")
  expect_identical(out$marker_symbol, "Nf1")
})

test_that("the ported IMPC phenotype fixture parses as expected", {
  body <- read_fixture("impc_pheno_nf1.json")
  out <- impc_parse_phenotypes(body, "MGI:97306")

  expect_s3_class(out, "tbl_df")
  expect_identical(out$allele[1], "MGI:4364806")
  expect_match(out$allele_symbol[1], "tm1a(KOMP)Wtsi", fixed = TRUE)
})

# --- One row per term, not per observation -----------------------------------

test_that("repeated observations of a term collapse to one row", {
  # IMPC reports a document per phenotype, sex, zygosity, and parameter, so the
  # fixture's 4 documents are 2 real phenotypes. Counting documents overstates
  # a gene's phenotype breadth several-fold.
  body <- read_fixture("impc_pheno_nf1.json")
  docs <- body$response$docs

  expect_length(docs, 4)
  expect_identical(nrow(impc_parse_phenotypes(body, "MGI:97306")), 2L)
})

test_that("MP and MPATH terms both count", {
  # mp_term_id carries MP terms for most phenotypes and MPATH terms for
  # pathology findings. Filtering to MP: would drop the second silently.
  body <- read_fixture("impc_pheno_nf1.json")
  out <- impc_parse_phenotypes(body, "MGI:97306")

  expect_setequal(out$mp_id, c("MP:0011100", "MPATH:212"))
})

test_that("the first observation of a term is the one kept", {
  body <- read_fixture("impc_pheno_nf1.json")
  out <- impc_parse_phenotypes(body, "MGI:97306")

  expect_identical(out$zygosity, c("homozygote", "heterozygote"))
})

# --- Empty and unusable shapes -----------------------------------------------

test_that("no documents is NULL for both parsers", {
  empty <- list(response = list(docs = list()))

  expect_null(impc_parse_ortholog(empty))
  expect_null(impc_parse_phenotypes(empty))
  expect_null(impc_parse_ortholog(list()))
})

test_that("a document with no usable accession is not an ortholog", {
  # The accession is interpolated into the next Solr query, so a blank or
  # malformed one would build a query for something else entirely.
  expect_null(impc_parse_ortholog(list(
    response = list(
      docs = list(
        list(mgi_accession_id = "", marker_symbol = "Nf1")
      )
    )
  )))
  expect_null(impc_parse_ortholog(list(
    response = list(
      docs = list(
        list(mgi_accession_id = "97306")
      )
    )
  )))
})

# --- Trap: the phenotype core has no human symbol ----------------------------

test_that("the phenotype core is queried by MGI accession, never by symbol", {
  # The genotype-phenotype core carries no human_gene_symbol field. Solr does
  # not reject a query against a field it lacks, so asking it for
  # human_gene_symbol:NF1 returns 200 with zero documents, which reads exactly
  # like "IMPC found no phenotype for this gene".
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    if (grepl("genotype-phenotype", req$url, fixed = TRUE)) {
      return(mock_json(
        '{"response":{"docs":[{"mp_term_id":"MP:1","mp_term_name":"x"}]}}'
      ))
    }
    mock_json(
      '{"response":{"docs":[{"mgi_accession_id":"MGI:97306","marker_symbol":"Nf1"}]}}'
    )
  })

  out <- biohttp::body_or_null(impc_gene_phenotypes("NF1"))

  expect_length(urls, 2)
  expect_match(urls[1], "human_gene_symbol%3ANF1")
  expect_match(urls[2], "marker_accession_id")
  expect_false(grepl("human_gene_symbol", urls[2], fixed = TRUE))
  expect_identical(out$mgi, "MGI:97306")
})

test_that("a gene with no mouse ortholog stops before the second query", {
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('{"response":{"docs":[]}}')
  })

  res <- impc_gene_phenotypes("NOTAGENE")

  expect_identical(res$status, "no_data")
  expect_length(urls, 1)
})

# --- A miss is never a zero --------------------------------------------------

test_that("no significant phenotype is no_data, not a count of zero", {
  # A gene IMPC never tested and a gene tested with no significant abnormality
  # both come back empty, and the response does not tell them apart. Reporting
  # 0 would let an untested gene be scored as if it had been tested.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    if (grepl("genotype-phenotype", req$url, fixed = TRUE)) {
      return(mock_json('{"response":{"docs":[]}}'))
    }
    mock_json(
      '{"response":{"docs":[{"mgi_accession_id":"MGI:97306","marker_symbol":"Nf1"}]}}'
    )
  })

  res <- impc_gene_phenotypes("NF1")

  expect_identical(res$status, "no_data")
  expect_null(res$data)
})

# --- The Solr query ----------------------------------------------------------

test_that("a symbol is stripped before it goes into a Solr expression", {
  # A colon or a quote in the symbol would change which field is searched.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"response":{"docs":[]}}')
  })

  impc_mouse_ortholog('NF1" OR marker_symbol:*')

  # The quote, the colon, the spaces, and the wildcard are all gone, so what is
  # left cannot escape the field it is searching. The underscore survives
  # because it is a legitimate gene-symbol character.
  expect_match(url, "human_gene_symbol%3ANF1ORMARKER_SYMBOL&")
})

test_that("no usable symbol never reaches the network", {
  reset_transport()
  expect_identical(impc_mouse_ortholog("")$status, "no_data")
  expect_identical(impc_gene_phenotypes(NULL)$status, "no_data")
})
