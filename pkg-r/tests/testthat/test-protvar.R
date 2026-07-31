# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported ProtVar function fixture parses to clean prose", {
  body <- read_fixture("protvar_function_p04637_175.json")
  out <- protvar_parse_function(body)

  expect_type(out, "character")
  expect_match(out, "Multifunctional transcription factor")
  # The citations that dominated the raw comment are gone.
  expect_false(grepl("PubMed:", out, fixed = TRUE))
})

test_that("the ported ProtVar population fixture parses to variants", {
  body <- read_fixture("protvar_population_p04637_175.json")
  out <- protvar_parse_population(body)

  expect_s3_class(out, "tbl_df")
  expect_true("Cys" %in% out$change)
  expect_identical(out$sources[out$change == "Cys"], "NCI-TCGA")
})

# --- Sources are deduplicated ------------------------------------------------

test_that("a source repeated per supporting record is reported once", {
  # A variant carries one cross-reference per record, so the same database
  # recurs many times. The fixture lists NCI-TCGA four times for one change.
  # Reporting them raw makes a single source look like corroboration.
  body <- read_fixture("protvar_population_p04637_175.json")
  raw <- body$variants[[1]]$xrefs

  expect_true(length(raw) > 1)
  expect_identical(protvar_parse_population(body)$sources[1], "NCI-TCGA")
})

test_that("several distinct sources are joined", {
  body <- list(
    variants = list(list(
      alternativeSequence = "Cys",
      xrefs = list(list(name = "NCI-TCGA"), list(name = "ClinVar"))
    ))
  )
  expect_identical(
    protvar_parse_population(body)$sources,
    "NCI-TCGA, ClinVar"
  )
})

# --- Citation stripping ------------------------------------------------------

test_that("citation-only parentheticals are removed whole", {
  expect_identical(
    protvar_strip_citations(
      "Induces arrest (PubMed:11025664, PubMed:12524540)."
    ),
    "Induces arrest."
  )
})

test_that("a non-citation parenthetical survives", {
  # "(By similarity)" is a real qualifier about the evidence, not noise.
  expect_identical(
    protvar_strip_citations("Binds DNA (By similarity)."),
    "Binds DNA (By similarity)."
  )
})

test_that("citations mixed into a real parenthetical are stripped in place", {
  expect_identical(
    protvar_strip_citations("Binds DNA (By similarity, PubMed:123)."),
    "Binds DNA (By similarity)."
  )
})

test_that("blank text passes through untouched", {
  expect_identical(protvar_strip_citations(""), "")
  expect_null(protvar_strip_citations(NULL))
})

# --- Position parsing --------------------------------------------------------

test_that("a position is read from every form a change arrives in", {
  expect_identical(protvar_position("p.Arg175His"), 175L)
  expect_identical(protvar_position("R175H"), 175L)
  expect_identical(protvar_position("Arg175His"), 175L)
  expect_identical(protvar_position("175"), 175L)
})

test_that("a string with no number yields NULL", {
  expect_null(protvar_position("no digits"))
  expect_null(protvar_position(""))
  expect_null(protvar_position(NULL))
})

# --- Empty results -----------------------------------------------------------

test_that("no function comment is NA and no variants is NULL", {
  expect_true(is.na(protvar_parse_function(list(comments = list()))))
  expect_true(is.na(protvar_parse_function(list())))
  expect_null(protvar_parse_population(list(variants = list())))
})

# --- The client half ---------------------------------------------------------

test_that("the two endpoints are separate calls, not one combined one", {
  # A caller usually wants one of them. Combining them spent two requests
  # whichever was asked for.
  reset_transport()
  urls <- character()
  httr2::local_mocked_responses(function(req) {
    urls <<- c(urls, req$url)
    mock_json('{"comments":[{"type":"FUNCTION","text":[{"value":"x"}]}]}')
  })

  protvar_function("P04637", 175)

  expect_length(urls, 1)
  expect_match(urls[1], "function/P04637/175")
})

test_that("a bad accession or position never reaches the network", {
  reset_transport()
  expect_identical(protvar_function("", 175)$status, "no_data")
  expect_identical(protvar_function("P04637", NA)$status, "no_data")
  expect_identical(protvar_population("P04637", "nope")$status, "no_data")
})
