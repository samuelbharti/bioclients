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

# --- A term that matches several records -------------------------------------

# Answers from the two recorded rs113488022 bodies and keeps every URL asked
# for. Both were recorded live on 2026-10-07 with the client's own requests:
# esearch with term=rs113488022 and retmax=20, then esummary with
# id=40389,13961. ClinVar lists V600G (40389) before V600E (13961).
clinvar_two_record_mock <- function(seen) {
  esearch <- read_fixture_text("clinvar_esearch_rs113488022.json")
  esummary <- read_fixture_text("clinvar_esummary_rs113488022.json")
  function(req) {
    seen$urls <- c(seen$urls, req$url)
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json(esearch))
    }
    mock_json(esummary)
  }
}

test_that("an rsID that matches two records returns the first and says so", {
  # This is the behaviour issue #40 found: the row is V600G, not V600E. It stays
  # the default so existing code keeps working, and n_matches now shows it.
  reset_transport()
  seen <- new.env()
  httr2::local_mocked_responses(clinvar_two_record_mock(seen))

  res <- clinvar_classification("rs113488022")
  out <- biohttp::body_or_null(res)

  expect_identical(out$accession, "VCV000040389")
  expect_identical(out$n_matches, 2L)
  expect_length(seen$urls, 2)
  expect_match(seen$urls[[2]], "40389(%2C|,)13961")
})

test_that("allele picks the record that has that protein change", {
  reset_transport()
  seen <- new.env()
  httr2::local_mocked_responses(clinvar_two_record_mock(seen))

  out <- biohttp::body_or_null(
    clinvar_classification("rs113488022", allele = "V600E")
  )

  expect_identical(out$accession, "VCV000013961")
  expect_identical(
    out$significance,
    "Conflicting classifications of pathogenicity"
  )
  expect_identical(out$n_matches, 1L)
  expect_length(seen$urls, 2)
})

test_that("allele takes the three-letter form too", {
  reset_transport()
  seen <- new.env()
  httr2::local_mocked_responses(clinvar_two_record_mock(seen))

  out <- biohttp::body_or_null(
    clinvar_classification("rs113488022", allele = "p.Val600Glu")
  )

  expect_identical(out$accession, "VCV000013961")
})

test_that("an allele no record has is no_data, not the first record", {
  # Falling back to the first record here would bring back the bug: the caller
  # asked for one variant and would get another.
  reset_transport()
  seen <- new.env()
  httr2::local_mocked_responses(clinvar_two_record_mock(seen))

  res <- clinvar_classification("rs113488022", allele = "V600K")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "V600K", fixed = TRUE)
  expect_length(seen$urls, 2)
})

test_that("a blank or NA allele means no filter", {
  reset_transport()
  seen <- new.env()
  httr2::local_mocked_responses(clinvar_two_record_mock(seen))

  out <- biohttp::body_or_null(
    clinvar_classification("rs113488022", allele = NA_character_)
  )

  expect_identical(out$accession, "VCV000040389")
  expect_identical(out$n_matches, 2L)
})

test_that("more than one allele is refused before any request", {
  reset_transport()
  calls <- 0L
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  expect_error(
    clinvar_classification("rs113488022", allele = c("V600E", "V600K")),
    "single"
  )
  expect_identical(calls, 0L)
})

test_that("a record esummary cannot summarise is skipped", {
  # ClinVar can list a UID it no longer has a summary for. esummary then sends
  # a record holding only `error`, and parsing that would give a row of NA.
  reset_transport()
  record <- paste(
    readLines(
      testthat::test_path("fixtures", "clinvar_40389.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) {
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json(
        '{"esearchresult":{"count":"2","idlist":["999","40389"]}}'
      ))
    }
    mock_json(paste0(
      '{"result":{"uids":["999","40389"],',
      '"999":{"error":"cannot get document summary","uid":"999"},',
      '"40389":',
      record,
      "}}"
    ))
  })

  out <- biohttp::body_or_null(clinvar_classification("rs113488022"))

  expect_identical(out$accession, "VCV000040389")
})

test_that("only error records is no_data, not a row of NA", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json('{"esearchresult":{"idlist":["999"]}}'))
    }
    mock_json(
      '{"result":{"999":{"error":"cannot get document summary","uid":"999"}}}'
    )
  })

  res <- clinvar_classification("rs113488022")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "no summary", fixed = TRUE)
})

test_that("retmax is sent, so the record limit is visible", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  clinvar_classification("rs113488022")

  expect_match(url, "retmax=20", fixed = TRUE)
})

test_that("protein changes are brought into the form ClinVar writes", {
  expect_identical(clinvar_protein_change("V600E"), "V600E")
  expect_identical(clinvar_protein_change("v600e"), "V600E")
  expect_identical(clinvar_protein_change("p.Val600Glu"), "V600E")
  expect_identical(clinvar_protein_change("p.(Val600Glu)"), "V600E")
  expect_identical(clinvar_protein_change("p.Arg213Ter"), "R213*")
  expect_identical(clinvar_protein_change("p.Gln1756fsTer74"), "Q1756fs")
  expect_identical(clinvar_protein_change("Q1756fs*74"), "Q1756fs")
})

test_that("a form that is not a protein change is passed through", {
  # It then matches no record, which is the honest answer.
  expect_identical(clinvar_protein_change("c.1799T>A"), "c.1799T>A")
})

test_that("a protein change is matched as a list, not as text", {
  record <- list(protein_change = "V600E, V512E, V578E")
  expect_true(clinvar_has_protein_change(record, "V512E"))
  expect_false(clinvar_has_protein_change(record, "V60"))
  expect_false(clinvar_has_protein_change(list(), "V600E"))
})

# --- The NCBI API key --------------------------------------------------------

test_that("the key is sent on both requests when one is configured", {
  # E-utilities takes it in the query string; there is no header form.
  reset_transport()
  withr::local_envvar(NCBI_API_KEY = "SECRET123")
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json('{"esearchresult":{"idlist":["40389"]}}'))
    }
    mock_json('{"result":{"40389":{"accession":"VCV000040389"}}}')
  })

  clinvar_classification("rs113488022")

  expect_length(urls, 2)
  expect_true(all(grepl("api_key=SECRET123", urls, fixed = TRUE)))
})

test_that("no key configured sends no api_key at all", {
  # An empty api_key= is not the same as omitting it.
  reset_transport()
  withr::local_envvar(NCBI_API_KEY = "")
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  clinvar_classification("rs113488022")

  expect_false(grepl("api_key", url, fixed = TRUE))
})

test_that("the key does not partition the cache", {
  # A rate-limit credential does not change the answer, so a call made with one
  # must be served from an entry warmed without one. Otherwise configuring or
  # rotating a key silently discards everything already fetched.
  reset_transport()
  calls <- 0L
  httr2::local_mocked_responses(function(req) {
    calls <<- calls + 1L
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  withr::with_envvar(c(NCBI_API_KEY = ""), clinvar_classification("rs1"))
  withr::with_envvar(
    c(NCBI_API_KEY = "SECRET123"),
    clinvar_classification("rs1")
  )

  expect_identical(calls, 1L)
})

test_that("a transport failure does not report the key back", {
  reset_transport()
  withr::local_envvar(NCBI_API_KEY = "SECRET123")
  httr2::local_mocked_responses(function(req) {
    stop("Could not resolve host: eutils.ncbi.nlm.nih.gov/?api_key=SECRET123")
  })

  res <- clinvar_classification("rs113488022")

  expect_false(isTRUE(res$ok))
  expect_false(grepl("SECRET123", res$detail, fixed = TRUE))
})

# --- The default throttle ----------------------------------------------------

test_that("the default throttle follows the documented rate for the key", {
  withr::with_envvar(c(NCBI_API_KEY = ""), {
    expect_identical(clinvar_throttle(), list(capacity = 3, fill_time_s = 1))
  })
  withr::with_envvar(c(NCBI_API_KEY = "SECRET123"), {
    expect_identical(clinvar_throttle(), list(capacity = 10, fill_time_s = 1))
  })
})

test_that("both requests are throttled by default", {
  # Without a throttle a loop over variants runs straight into the E-utilities
  # limit and reads as a flaky service rather than as a client sending too
  # fast. The default is supplied so a caller does not have to know that.
  reset_transport()
  withr::local_envvar(NCBI_API_KEY = "")
  realms <- character()
  httr2::local_mocked_responses(function(req) {
    realms <<- c(realms, req$policies$throttle_realm %||% NA_character_)
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json('{"esearchresult":{"idlist":["40389"]}}'))
    }
    mock_json('{"result":{"40389":{"accession":"VCV000040389"}}}')
  })

  clinvar_classification("rs113488022")

  expect_length(realms, 2)
  expect_false(anyNA(realms))
})

test_that("a caller can still supply its own throttle", {
  reset_transport()
  realm <- NULL
  httr2::local_mocked_responses(function(req) {
    realm <<- req$policies$throttle_realm
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  clinvar_classification(
    "rs113488022",
    throttle = list(capacity = 1, fill_time_s = 1, realm = "mine")
  )

  expect_identical(realm, "mine")
})

# --- Identifying the caller --------------------------------------------------

test_that("tool and email are sent when the identity is configured", {
  reset_transport()
  withr::local_envvar(
    BIOHTTP_CALLER_IDENTITY = "myapp",
    BIOHTTP_CONTACT_EMAIL = "dev@example.org"
  )
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    if (grepl("esearch", req$url, fixed = TRUE)) {
      return(mock_json('{"esearchresult":{"idlist":["40389"]}}'))
    }
    mock_json('{"result":{"40389":{"accession":"VCV000040389"}}}')
  })

  clinvar_classification("rs113488022")

  expect_length(urls, 2)
  expect_true(all(grepl("tool=myapp", urls, fixed = TRUE)))
  expect_true(all(grepl("email=dev%40example.org", urls, fixed = TRUE)))
})

test_that("a blank identity sends neither parameter", {
  # An empty tool= is not the same as omitting it.
  reset_transport()
  withr::local_envvar(BIOHTTP_CALLER_IDENTITY = "", BIOHTTP_CONTACT_EMAIL = "")
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"esearchresult":{"idlist":[]}}')
  })

  clinvar_classification("rs113488022")

  expect_false(grepl("tool=", url, fixed = TRUE))
  expect_false(grepl("email=", url, fixed = TRUE))
})

test_that("one half of the identity is sent without the other", {
  withr::with_envvar(
    c(BIOHTTP_CALLER_IDENTITY = "myapp", BIOHTTP_CONTACT_EMAIL = ""),
    expect_identical(clinvar_identity_query(), list(tool = "myapp"))
  )
  withr::with_envvar(
    c(BIOHTTP_CALLER_IDENTITY = "", BIOHTTP_CONTACT_EMAIL = "a@b.org"),
    expect_identical(clinvar_identity_query(), list(email = "a@b.org"))
  )
})
