# --- The ported fixtures, which must pass unchanged ---------------------------

test_that("the ported Monarch disease search fixture parses as expected", {
  body <- read_fixture("monarch_search_marfan.json")
  out <- monarch_parse_search(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$id[1], "MONDO:0007947")
  expect_identical(out$name[1], "Marfan syndrome")
  expect_identical(out$category[1], "biolink:Disease")
})

test_that("the ported Monarch association fixture parses as expected", {
  body <- read_fixture("monarch_disease_phenotype.json")
  out <- monarch_parse_associations(body)

  expect_identical(nrow(out), 2L)
  expect_identical(out$subject[1], "MONDO:0017309")
  expect_identical(out$object[1], "HP:0001653")
  expect_identical(out$object_label[1], "Mitral regurgitation")
  expect_identical(out$primary_knowledge_source[1], "infores:orphanet")
})

# --- Trap: a gene symbol is not unique across species ------------------------

test_that("a gene search returns the orthologs under the same symbol", {
  # FBN1 returns the human gene and then the chicken and dog orthologs, all
  # three named FBN1. Taking the first result is right here and wrong in
  # general, so the taxon has to be on the row.
  body <- read_fixture("monarch_search_fbn1.json")
  out <- monarch_parse_search(body)

  expect_identical(out$name, c("FBN1", "FBN1", "FBN1"))
  expect_identical(
    out$taxon,
    c("Homo sapiens", "Gallus gallus", "Canis lupus familiaris")
  )
  expect_identical(out$id[1], "HGNC:3603")
})

test_that("a disease match has no taxon rather than a wrong one", {
  body <- read_fixture("monarch_search_marfan.json")
  expect_true(all(is.na(monarch_parse_search(body)$taxon)))
})

# --- The readable gloss ------------------------------------------------------

test_that("a gene falls back to full_name and a disease uses description", {
  # They are different fields, and a parser reading only one leaves half the
  # results with no description at all.
  gene <- read_fixture("monarch_search_fbn1.json")
  expect_identical(monarch_parse_search(gene)$description[1], "fibrillin 1")

  disease <- read_fixture("monarch_search_marfan.json")
  expect_false(is.na(monarch_parse_search(disease)$description[1]))
})

# --- Publications ------------------------------------------------------------

test_that("publications is a list column that survives null and an array", {
  # The fixture sends null; the gene-phenotype fixture sends an array. A parser
  # assuming either shape breaks on the other.
  none <- monarch_parse_associations(read_fixture(
    "monarch_disease_phenotype.json"
  ))
  expect_type(none$publications, "list")
  expect_identical(none$publications[[1]], character())

  some <- monarch_parse_associations(read_fixture(
    "monarch_phenotypes_tp53.json"
  ))
  expect_identical(some$publications[[1]], "PMID:30146126")
})

test_that("an association list can be passed without its wrapper", {
  items <- list(list(subject = "A", predicate = "p", object = "B"))
  expect_identical(nrow(monarch_parse_associations(items)), 1L)
})

test_that("no items is NULL for both parsers", {
  expect_null(monarch_parse_search(list(items = list())))
  expect_null(monarch_parse_associations(list(items = list())))
  expect_null(monarch_parse_search(list()))
})

# --- The HGNC CURIE ----------------------------------------------------------

test_that("an HGNC id is normalised to the CURIE Monarch expects", {
  # MyGene returns bare digits; Monarch's path wants HGNC:11998.
  expect_identical(monarch_hgnc_id("11998"), "HGNC:11998")
  expect_identical(monarch_hgnc_id("HGNC:11998"), "HGNC:11998")
  expect_identical(monarch_hgnc_id("hgnc:11998"), "HGNC:11998")
  expect_identical(monarch_hgnc_id(" 11998 "), "HGNC:11998")
})

test_that("anything that is not an HGNC id is refused", {
  # It goes into a URL path, so a symbol here would request some other entity.
  expect_null(monarch_hgnc_id("TP53"))
  expect_null(monarch_hgnc_id(""))
  expect_null(monarch_hgnc_id(NULL))
  expect_null(monarch_hgnc_id("11998/../../x"))
})

# --- The client half ---------------------------------------------------------

test_that("search reports Monarch's own total, not the page size", {
  # The fixture shows 3 of 26. Presenting 3 as everything would be wrong.
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"total":26,"items":[{"id":"MONDO:1","name":"X"}]}')
  })

  out <- biohttp::body_or_null(monarch_search("Marfan syndrome"))

  expect_identical(out$total, 26L)
  expect_identical(nrow(out$matches), 1L)
})

test_that("the association end is the parameter that changes", {
  # subject= and object= are the two directions of the same edge.
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('{"items":[{"subject":"A","object":"B"}]}')
  })

  monarch_associations("MONDO:0007947", end = "subject")
  monarch_associations("MONDO:0007947", end = "object")

  expect_match(urls[1], "subject=MONDO", fixed = TRUE)
  expect_match(urls[2], "object=MONDO", fixed = TRUE)
})

test_that("an unknown end is refused rather than sent", {
  reset_transport()
  expect_error(monarch_associations("MONDO:1", end = "predicate"))
})

test_that("the CURIE colon reaches Monarch unencoded", {
  # Both path segments are CURIEs, and Monarch wants the colon verbatim rather
  # than as %3A. This pins that the request builder leaves it alone.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"total":48,"items":[{"object":"HP:1"}]}')
  })

  monarch_gene_phenotypes("11998")

  expect_match(url, "entity/HGNC:11998/", fixed = TRUE)
  expect_match(url, "biolink:GeneToPhenotypicFeatureAssociation", fixed = TRUE)
  expect_false(grepl("%3A", url, fixed = TRUE))
})

test_that("a gene with no HGNC id never reaches the network", {
  reset_transport()
  res <- monarch_gene_phenotypes("TP53")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "HGNC")
})

test_that("blank input is no_data for every entry point", {
  reset_transport()
  expect_identical(monarch_search("")$status, "no_data")
  expect_identical(monarch_associations("")$status, "no_data")
  expect_identical(monarch_gene_phenotypes(NULL)$status, "no_data")
})
