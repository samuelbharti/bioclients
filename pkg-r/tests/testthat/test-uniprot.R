# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported UniProt disease fixture parses as expected", {
  body <- read_fixture("uniprot_disease_nf1.json")
  out <- uniprot_parse_diseases(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$id[1], "DI-02396")
  expect_identical(out$name[1], "Neurofibromatosis 1")
  expect_identical(out$acronym[1], "NF1")
  expect_identical(out$mim[1], "162200")
  expect_identical(
    out$source_url[1],
    "https://www.uniprot.org/diseases/DI-02396"
  )
})

test_that("the ported Proteins features fixture parses as expected", {
  body <- read_fixture("proteins_features_p15056.json")
  out <- uniprot_parse_features(body)

  expect_s3_class(out, "tbl_df")
  expect_true(nrow(out) > 0)
  expect_true(all(!is.na(out$begin)))
  # Sorted by start position.
  expect_true(all(diff(out$begin) >= 0))
  # Terse codes get a readable label.
  expect_true(all(out$label != "" & !is.na(out$label)))
})

# --- The causal distinction --------------------------------------------------

test_that("caused-by is told apart from may-be-involved", {
  # UniProt files both as DISEASE comments, but "caused by variants" is a
  # Mendelian claim and "may be involved in the pathogenesis" is far weaker.
  # Treating them alike promotes a speculative association to a causal one.
  causal <- list(
    comments = list(list(
      commentType = "DISEASE",
      disease = list(diseaseAccession = "DI-1", diseaseId = "A"),
      note = list(
        texts = list(list(
          value = "The disease is caused by variants affecting the gene"
        ))
      )
    ))
  )
  weak <- list(
    comments = list(list(
      commentType = "DISEASE",
      disease = list(diseaseAccession = "DI-2", diseaseId = "B"),
      note = list(
        texts = list(list(
          value = "The gene may be involved in the pathogenesis"
        ))
      )
    ))
  )

  expect_true(uniprot_parse_diseases(causal)$causal)
  expect_false(uniprot_parse_diseases(weak)$causal)
})

test_that("non-disease comments are ignored", {
  body <- list(
    comments = list(
      list(commentType = "FUNCTION", texts = list(list(value = "Kinase"))),
      list(
        commentType = "DISEASE",
        disease = list(diseaseAccession = "DI-1", diseaseId = "A")
      )
    )
  )
  expect_identical(nrow(uniprot_parse_diseases(body)), 1L)
})

test_that("a disease comment naming nothing is dropped", {
  # Nothing a caller could ground a claim on.
  body <- list(
    comments = list(list(
      commentType = "DISEASE",
      disease = list(diseaseAccession = "", diseaseId = "")
    ))
  )
  expect_null(uniprot_parse_diseases(body))
})

test_that("a MIM cross-reference is only read when it is a MIM one", {
  body <- list(
    comments = list(list(
      commentType = "DISEASE",
      disease = list(
        diseaseAccession = "DI-1",
        diseaseId = "A",
        diseaseCrossReference = list(database = "MedGen", id = "C0027831")
      )
    ))
  )
  expect_true(is.na(uniprot_parse_diseases(body)$mim))
})

test_that("an entry with no disease comments parses to NULL", {
  expect_null(uniprot_parse_diseases(list(comments = list())))
  expect_null(uniprot_parse_diseases(list()))
})

# --- Features ----------------------------------------------------------------

test_that("no features is a zero-row tibble, not NULL", {
  # "This protein has no annotated domains" is a real answer.
  out <- uniprot_parse_features(list(features = list()))
  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 0L)
})

test_that("an unmapped feature type falls back to its own code", {
  # A new UniProt type should still read sensibly rather than becoming NA.
  body <- list(
    features = list(
      list(type = "SOMETHING_NEW", description = "x", begin = "1", end = "2")
    )
  )
  expect_identical(uniprot_parse_features(body)$label, "SOMETHING_NEW")
})

test_that("begin and end arrive as strings and become integers", {
  body <- list(
    features = list(
      list(type = "DOMAIN", description = "K", begin = "457", end = "717")
    )
  )
  out <- uniprot_parse_features(body)

  expect_type(out$begin, "integer")
  expect_identical(out$begin, 457L)
  expect_identical(out$end, 717L)
})

test_that("features spanning a residue are selected", {
  features <- uniprot_parse_features(list(
    features = list(
      list(type = "DOMAIN", description = "Kinase", begin = "457", end = "717"),
      list(type = "SITE", description = "Active", begin = "576", end = "576")
    )
  ))

  expect_identical(nrow(uniprot_features_at(features, 600)), 1L)
  expect_identical(nrow(uniprot_features_at(features, 576)), 2L)
  expect_identical(nrow(uniprot_features_at(features, 100)), 0L)
  expect_identical(nrow(uniprot_features_at(features, NA)), 0L)
})

# --- The client half ---------------------------------------------------------

test_that("the two endpoints are labelled as the different services they are", {
  # rest.uniprot.org and the EBI Proteins API are separate services with
  # separate availability. One being down must not take the other out, which
  # means separate source labels and separate circuit breakers.
  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json('{"comments":[]}'))
  expect_identical(uniprot_diseases("P04637")$source, "UniProt")

  reset_transport()
  httr2::local_mocked_responses(function(req) mock_json('{"features":[]}'))
  expect_identical(uniprot_features("P15056")$source, "EBI Proteins")
})

test_that("uniprot_diseases returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- paste(
    readLines(
      testthat::test_path("fixtures", "uniprot_disease_nf1.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- uniprot_diseases("P21359")

  expect_true(res$ok)
  expect_identical(biohttp::body_or_null(res)$acronym[1], "NF1")
})

test_that("blank accessions are no_data and never reach the network", {
  reset_transport()
  expect_identical(uniprot_diseases("")$status, "no_data")
  expect_identical(uniprot_features("")$status, "no_data")
})
