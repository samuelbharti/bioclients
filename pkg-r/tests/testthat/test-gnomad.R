# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported gnomAD constraint fixture parses to the expected values", {
  body <- read_fixture("gnomad_constraint_braf.json")
  out <- gnomad_parse_constraint(body, "BRAF")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$symbol, "BRAF")
  expect_equal(out$loeuf, 0.23234095358685583)
  expect_equal(out$oe_lof, 0.1529742462665005)
  expect_equal(out$mis_z, 5.520263159451872)
  expect_equal(out$syn_z, 1.0194537642329518)
  expect_equal(out$lof_z, 7.348724683473751)
  expect_equal(out$oe_mis, 0.5772130672664564)
  expect_true(out$pli > 0.999)
})

test_that("loeuf is gnomAD's oe_lof_upper", {
  # A consuming app's entire ranking model reads this one number, under the name
  # LOEUF, from the API field oe_lof_upper. Renaming it silently would break
  # ranking in a way no test over there would catch.
  body <- read_fixture("gnomad_constraint_braf.json")
  raw <- body$data$gene$gnomad_constraint$oe_lof_upper
  expect_equal(gnomad_parse_constraint(body, "BRAF")$loeuf, raw)
})

test_that("a gene with no constraint block parses to NULL", {
  body <- list(data = list(gene = list(gnomad_constraint = NULL)))
  expect_null(gnomad_parse_constraint(body, "NOSUCH"))
})

# --- Both query types exist, which is the reconciliation ---------------------

test_that("frequency and constraint are separate entry points", {
  # The three copies this replaces disagreed because they answered different
  # questions. Neither may be folded into the other, and constraint may not be
  # dropped just because two of the three callers wanted frequency.
  expect_true(is.function(gnomad_constraint))
  expect_true(is.function(gnomad_frequency))
  expect_false(identical(gnomad_constraint, gnomad_frequency))
})

# --- Populations -------------------------------------------------------------

test_that("exome and genome counts are summed per ancestry", {
  out <- gnomad_parse_populations(
    list(list(id = "nfe", ac = 3, an = 1000)),
    list(list(id = "nfe", ac = 1, an = 500))
  )
  expect_identical(nrow(out), 1L)
  expect_equal(out$ac, 4)
  expect_equal(out$an, 1500)
  expect_equal(out$af, 4 / 1500)
})

test_that("sex splits are dropped so each ancestry appears once", {
  # Keeping nfe_XX alongside nfe would double-count the same samples and put two
  # bars where there should be one.
  out <- gnomad_parse_populations(
    list(
      list(id = "nfe", ac = 3, an = 1000),
      list(id = "nfe_XX", ac = 2, an = 500),
      list(id = "XY", ac = 1, an = 500)
    ),
    NULL
  )
  expect_identical(out$pop, "nfe")
  expect_equal(out$ac, 3)
})

test_that("populations are sorted by frequency, most common first", {
  out <- gnomad_parse_populations(
    list(
      list(id = "nfe", ac = 1, an = 1000),
      list(id = "eas", ac = 50, an = 1000)
    ),
    NULL
  )
  expect_identical(out$pop, c("eas", "nfe"))
})

test_that("an ancestry with no called alleles is dropped, not divided by zero", {
  out <- gnomad_parse_populations(list(list(id = "nfe", ac = 0, an = 0)), NULL)
  expect_null(out)
})

test_that("codes get their display label", {
  out <- gnomad_parse_populations(list(list(id = "asj", ac = 1, an = 10)), NULL)
  expect_identical(out$label, "Ashkenazi Jewish")
})

# --- Frequency ---------------------------------------------------------------

test_that("a variant response parses into the frequency record", {
  body <- list(
    data = list(
      variant = list(
        variant_id = "7-140753336-A-T",
        exome = list(af = 0.001, ac = 2, an = 2000),
        genome = NULL
      )
    )
  )
  out <- gnomad_parse_frequency(body)

  expect_identical(out$variant_id, "7-140753336-A-T")
  expect_equal(out$exome$ac, 2)
  expect_null(out$genome)
})

test_that("a missing variant parses to NULL", {
  expect_null(gnomad_parse_frequency(list(data = list(variant = NULL))))
})

# --- The aliased batch -------------------------------------------------------

test_that("aliased results map back by index, in input order", {
  body <- list(
    data = list(
      g1 = list(symbol = "BRAF", gnomad_constraint = list(oe_lof_upper = 0.23)),
      g2 = list(symbol = "TP53", gnomad_constraint = list(oe_lof_upper = 0.55))
    )
  )
  out <- gnomad_parse_constraints(body, c("BRAF", "TP53"))

  expect_identical(out$symbol, c("BRAF", "TP53"))
  expect_equal(out$loeuf, c(0.23, 0.55))
})

test_that("a gene with no constraint keeps its row rather than shifting", {
  body <- list(
    data = list(
      g1 = list(symbol = "BRAF", gnomad_constraint = list(oe_lof_upper = 0.23)),
      g2 = NULL,
      g3 = list(symbol = "EGFR", gnomad_constraint = list(oe_lof_upper = 0.11))
    )
  )
  out <- gnomad_parse_constraints(body, c("BRAF", "MISSING", "EGFR"))

  expect_identical(nrow(out), 3L)
  expect_identical(out$symbol, c("BRAF", "MISSING", "EGFR"))
  expect_true(is.na(out$loeuf[2]))
  expect_equal(out$loeuf[3], 0.11)
})

test_that("the chunk size stays under gnomAD's query cost cap", {
  # gnomAD's real limit is a query COST cap of 25, not a request count, and each
  # alias is charged. A chunk over the cap fails in a way that reads like rate
  # limiting and is not.
  expect_lte(GNOMAD_CHUNK, 25)
})

# --- The client half ---------------------------------------------------------

test_that("gnomad_constraint returns an ok envelope carrying the table", {
  reset_transport()
  fixture <- readLines(
    testthat::test_path("fixtures", "gnomad_constraint_braf.json"),
    warn = FALSE
  )
  httr2::local_mocked_responses(function(req) {
    mock_json(paste(fixture, collapse = ""))
  })

  res <- gnomad_constraint("BRAF")

  expect_true(res$ok)
  expect_identical(res$source, "gnomAD")
  expect_equal(biohttp::body_or_null(res)$loeuf, 0.23234095358685583)
})

test_that("a GraphQL errors array inside a 200 is a failure", {
  # This is the trap GraphQL sets: the HTTP status says 200 and the body says
  # the query failed. Treating the 200 as success returns an empty result that
  # looks like "no data" rather than "your query was wrong".
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"errors":[{"message":"Cannot query field"}]}')
  })

  res <- gnomad_constraint("BRAF")

  expect_false(res$ok)
  expect_identical(res$status, "error")
})

test_that("a blank symbol is no_data and never reaches the network", {
  reset_transport()
  expect_identical(gnomad_constraint("")$status, "no_data")
  expect_identical(gnomad_frequency("")$status, "no_data")
})
