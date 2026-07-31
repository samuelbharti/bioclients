# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported ClinVar fixture parses to the expected record", {
  record <- read_fixture("clinvar_40389.json")
  out <- clinvar_parse_record(record, "40389")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$uid, "40389")
  expect_identical(out$accession, "VCV000040389")
  expect_identical(out$title, "NM_004333.6(BRAF):c.1799T>G (p.Val600Gly)")
  expect_identical(out$significance, "Pathogenic")
  expect_identical(out$review_status, "reviewed by expert panel")
  expect_identical(out$last_evaluated, "2020/06/25 00:00")
  expect_identical(out$conditions, "RASopathy")
})

test_that("significance is read from the germline block, not the top level", {
  # ClinVar nests the classification under germline_classification. Reading a
  # top-level `description` would find nothing and report NA for every record.
  record <- read_fixture("clinvar_40389.json")
  expect_identical(
    clinvar_parse_record(record)$significance,
    record$germline_classification$description
  )
})

test_that("an absent record parses to NULL", {
  expect_null(clinvar_parse_record(NULL, "40389"))
})

# --- Conditions --------------------------------------------------------------

test_that("several traits collapse into one string", {
  germline <- list(
    trait_set = list(
      list(trait_name = "RASopathy"),
      list(trait_name = "Noonan syndrome")
    )
  )
  expect_identical(clinvar_conditions(germline), "RASopathy; Noonan syndrome")
})

test_that("duplicate trait names appear once", {
  germline <- list(
    trait_set = list(
      list(trait_name = "RASopathy"),
      list(trait_name = "RASopathy")
    )
  )
  expect_identical(clinvar_conditions(germline), "RASopathy")
})

test_that("no trait set is NA rather than an empty string", {
  expect_true(is.na(clinvar_conditions(list())))
  expect_true(is.na(clinvar_conditions(NULL)))
})

# --- Categories --------------------------------------------------------------

test_that("conflicting is tested before pathogenic", {
  # A conflicting record almost always contains the word "pathogenic" too, so
  # checking pathogenic first would file every conflicting record as pathogenic.
  # That turns a disputed call into a confident one, which is the wrong way to
  # be wrong.
  expect_identical(
    clinvar_category("Conflicting interpretations of pathogenicity"),
    "Conflicting"
  )
  expect_identical(clinvar_category("Pathogenic"), "Pathogenic / likely")
  expect_identical(clinvar_category("Likely pathogenic"), "Pathogenic / likely")
})

test_that("the remaining categories bucket as expected", {
  expect_identical(clinvar_category("Benign"), "Benign / likely")
  expect_identical(clinvar_category("Likely benign"), "Benign / likely")
  expect_identical(clinvar_category("Uncertain significance"), "Uncertain")
  expect_identical(clinvar_category("drug response"), "Other")
  expect_identical(clinvar_category(NULL), "Other")
})

# --- The client half ---------------------------------------------------------

test_that("clinvar_classification chains esearch into esummary", {
  reset_transport()
  record <- paste(
    readLines(
      testthat::test_path("fixtures", "clinvar_40389.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  seen <- character()
  httr2::local_mocked_responses(function(req) {
    seen <<- c(seen, req$url)
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json('{"esearchresult":{"idlist":["40389"]}}'))
    }
    mock_json(paste0('{"result":{"40389":', record, "}}"))
  })

  res <- clinvar_classification("rs113488022")

  expect_true(res$ok)
  expect_identical(res$source, "ClinVar")
  expect_identical(biohttp::body_or_null(res)$accession, "VCV000040389")
  expect_length(seen, 2)
  expect_match(seen[1], "esearch")
  expect_match(seen[2], "esummary")
})

test_that("an empty idlist stops before esummary", {
  reset_transport()
  calls <- 0L
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  res <- clinvar_classification("rs000000000")

  expect_identical(res$status, "no_data")
  # One request, not two. Asking esummary for a UID we do not have would be a
  # second round trip that cannot succeed.
  expect_identical(calls, 1L)
})

test_that("a blank term is no_data and never reaches the network", {
  reset_transport()
  expect_identical(clinvar_classification("")$status, "no_data")
})

test_that("a failed esearch passes its envelope straight through", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    httr2::response(status_code = 500)
  })
  res <- clinvar_classification("rs113488022")

  expect_false(res$ok)
  expect_identical(res$status, "error")
  reset_transport()
})
