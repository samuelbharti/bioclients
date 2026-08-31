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

# --- Frequency by variant id -------------------------------------------------

test_that("a variant id is chrom-pos-ref-alt without a chr prefix", {
  expect_identical(
    gnomad_variant_id("1", 55516888, "G", "GA"),
    "1-55516888-G-GA"
  )
  expect_identical(
    gnomad_variant_id("chr1", 55516888, "g", "ga"),
    "1-55516888-G-GA"
  )
  expect_identical(
    gnomad_variant_id(
      c("17", "chrX"),
      c(7676154, 10),
      c("G", "A"),
      c("C", "T")
    ),
    c("17-7676154-G-C", "X-10-A-T")
  )
})

test_that("the recorded present variant parses to the expected row", {
  # Recorded live from gnomad.broadinstitute.org on 2026-08-30 for TP53
  # p.Pro72Arg, with the HGDP and 1000 Genomes subsets trimmed out of the
  # population arrays.
  body <- read_fixture("gnomad_variant_present.json")
  out <- gnomad_parse_variant(body, "17-7676154-G-C")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 1L)
  expect_identical(out$variant_id, "17-7676154-G-C")
  expect_identical(out$rsid, "rs1042522")
  expect_equal(out$exome_af, 0.7163184765845761)
  expect_equal(out$exome_ac, 1046941)
  expect_equal(out$exome_an, 1461558)
  expect_equal(out$exome_nhomalt, 380188)
  expect_equal(out$genome_af, 0.626776035)
  expect_equal(out$genome_ac, 95285)
  expect_equal(out$genome_an, 152024)
  expect_equal(out$genome_nhomalt, 31776)
  expect_equal(out$faf95, 0.74639068)
  expect_identical(out$faf95_pop, "nfe")
  # An empty filters array means the site passed.
  expect_true(is.na(out$filters))
})

test_that("nhomalt is gnomAD's homozygote_count", {
  body <- read_fixture("gnomad_variant_present.json")
  expect_equal(
    gnomad_parse_variant(body)$exome_nhomalt,
    body$data$variant$exome$homozygote_count
  )
})

test_that("grpmax is the top eligible group over summed exome and genome counts", {
  body <- read_fixture("gnomad_variant_present.json")
  out <- gnomad_parse_variant(body, "17-7676154-G-C")
  pops <- gnomad_parse_populations(
    body$data$variant$exome$populations,
    body$data$variant$genome$populations
  )
  eligible <- pops[!(pops$pop %in% GNOMAD_GRPMAX_EXCLUDED), ]

  expect_identical(out$grpmax_id, eligible$pop[1])
  expect_equal(out$grpmax_af, eligible$af[1])
  expect_equal(out$grpmax_an, eligible$an[1])
  expect_identical(out$grpmax_id, "nfe")
})

test_that("a bottlenecked group never becomes grpmax", {
  # Finnish and Ashkenazi founder effects are exactly what grpmax exists to
  # leave out. A rarity filter reading a Finnish-only frequency as the group
  # maximum would call a globally rare variant common.
  variant <- list(
    exome = list(
      populations = list(
        list(id = "fin", ac = 50, an = 100),
        list(id = "nfe", ac = 1, an = 100)
      )
    )
  )
  out <- gnomad_variant_row(variant, "x")

  expect_identical(out$grpmax_id, "nfe")
  expect_equal(out$grpmax_af, 0.01)
})

test_that("faf95 is the higher of the two sample sets, with its group", {
  variant <- list(
    exome = list(faf95 = list(popmax = 0.01, popmax_population = "afr")),
    genome = list(faf95 = list(popmax = 0.02, popmax_population = "eas"))
  )
  out <- gnomad_variant_row(variant, "x")

  expect_equal(out$faf95, 0.02)
  expect_identical(out$faf95_pop, "eas")
})

test_that("filters from both sample sets are joined", {
  variant <- list(
    exome = list(filters = list("AC0")),
    genome = list(filters = list("AS_VQSR", "AC0"))
  )
  expect_identical(gnomad_variant_row(variant, "x")$filters, "AC0;AS_VQSR")
})

test_that("the recorded absent variant parses to NULL", {
  # gnomAD answers a variant it has never seen with a 200, a null variant and
  # a "Variant not found" entry in errors.
  body <- read_fixture("gnomad_variant_absent.json")

  expect_null(gnomad_parse_variant(body, "17-7676154-G-GTTTTT"))
  expect_identical(body$errors[[1]]$message, "Variant not found")
})

test_that("a variant with only one sample set fills the other with NA", {
  body <- list(
    data = list(
      variant = list(
        variant_id = "x",
        genome = list(af = 0.5, ac = 1, an = 2, homozygote_count = 0)
      )
    )
  )
  out <- gnomad_parse_variant(body, "x")

  expect_true(is.na(out$exome_af))
  expect_equal(out$genome_af, 0.5)
  expect_true(is.na(out$rsid))
  expect_true(is.na(out$grpmax_af))
})

# --- The aliased variant batch -----------------------------------------------

test_that("aliased variants map back by index, with NA for an absent one", {
  present <- read_fixture("gnomad_variant_present.json")$data$variant
  body <- list(
    errors = list(list(message = "Variant not found")),
    data = list(v1 = NULL, v2 = present)
  )
  out <- gnomad_parse_variants(body, c("17-7676154-G-GTTTTT", "17-7676154-G-C"))

  expect_identical(nrow(out), 2L)
  expect_identical(out$variant_id, c("17-7676154-G-GTTTTT", "17-7676154-G-C"))
  expect_true(is.na(out$exome_af[1]))
  expect_identical(out$rsid[2], "rs1042522")
})

test_that("the variant alias query names one alias per id", {
  query <- gnomad_variant_alias_query(c("1-1-A-T", "2-2-C-G"), "gnomad_r4")

  expect_match(
    query,
    'v1: variant(variantId: "1-1-A-T", dataset: gnomad_r4)',
    fixed = TRUE
  )
  expect_match(
    query,
    'v2: variant(variantId: "2-2-C-G", dataset: gnomad_r4)',
    fixed = TRUE
  )
  expect_match(query, "homozygote_count", fixed = TRUE)
})

test_that("Variant not found is an answer, anything else in errors is a failure", {
  not_found <- biohttp::status_ok(
    data = list(errors = list(list(message = "Variant not found"))),
    source = "gnomAD"
  )
  expect_null(gnomad_variant_error(not_found))

  bad <- biohttp::status_ok(
    data = list(
      errors = list(
        list(message = "Variant not found"),
        list(
          message = "Query is too expensive (26). Maximum allowed cost is 25."
        )
      )
    ),
    source = "gnomAD"
  )
  expect_identical(gnomad_variant_error(bad)$status, "error")
})

# --- The client half ---------------------------------------------------------

test_that("gnomad_frequency_by_id returns an ok envelope carrying the row", {
  reset_transport()
  fixture <- readLines(
    testthat::test_path("fixtures", "gnomad_variant_present.json"),
    warn = FALSE
  )
  sent <- NULL
  httr2::local_mocked_responses(function(req) {
    sent <<- req$body$data
    mock_json(paste(fixture, collapse = ""))
  })

  res <- gnomad_frequency_by_id("17-7676154-G-C")

  expect_true(res$ok)
  expect_identical(res$source, "gnomAD")
  expect_identical(biohttp::body_or_null(res)$rsid, "rs1042522")
  expect_match(sent$query, "variantId: $id", fixed = TRUE)
  expect_match(sent$query, "dataset: gnomad_r4", fixed = TRUE)
  expect_identical(sent$variables$id, "17-7676154-G-C")
})

test_that("an absent variant is no_data, not an error", {
  reset_transport()
  fixture <- readLines(
    testthat::test_path("fixtures", "gnomad_variant_absent.json"),
    warn = FALSE
  )
  httr2::local_mocked_responses(function(req) {
    mock_json(paste(fixture, collapse = ""))
  })

  res <- gnomad_frequency_by_id("17-7676154-G-GTTTTT")

  expect_identical(res$status, "no_data")
  expect_match(res$detail, "no record", fixed = TRUE)
})

test_that("a dataset on the other assembly is refused before any request", {
  reset_transport()
  res <- gnomad_frequency_by_id("17-7676154-G-C", reference_genome = "GRCh37")
  expect_identical(res$status, "no_data")
  expect_identical(
    gnomad_frequencies("17-7676154-G-C", dataset = "gnomad_r2_1")$status,
    "no_data"
  )
})

test_that("a blank id is no_data and never reaches the network", {
  reset_transport()
  expect_identical(gnomad_frequency_by_id("")$status, "no_data")
  expect_identical(gnomad_frequencies(c("", NA))$status, "no_data")
})

test_that("gnomad_frequencies chunks, batches by alias and keeps input order", {
  skip_if_not_installed("jsonlite")
  reset_transport()
  present <- read_fixture("gnomad_variant_present.json")$data$variant
  queries <- character()
  httr2::local_mocked_responses(function(req) {
    query <- req$body$data$query
    queries <<- c(queries, query)
    ids <- regmatches(query, gregexpr('variantId: "[^"]+"', query))[[1]]
    ids <- sub('variantId: "([^"]+)"', "\\1", ids)
    data <- list()
    errors <- list()
    for (i in seq_along(ids)) {
      alias <- paste0("v", i)
      if (identical(ids[i], "17-7676154-G-C")) {
        data[[alias]] <- present
      } else {
        data[alias] <- list(NULL)
        errors <- c(errors, list(list(message = "Variant not found")))
      }
    }
    body <- list(data = data)
    if (length(errors) > 0) {
      body$errors <- errors
    }
    mock_json(jsonlite::toJSON(body, auto_unbox = TRUE, null = "null"))
  })

  ids <- c("1-1-A-T", "17-7676154-G-C", "2-2-C-G")
  res <- gnomad_frequencies(ids, chunk_size = 2)

  expect_true(res$ok)
  expect_length(queries, 2)
  out <- biohttp::body_or_null(res)
  expect_identical(out$variant_id, ids)
  expect_identical(out$rsid, c(NA, "rs1042522", NA))
})

test_that("a failed chunk yields NA rows rather than failing the call", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"errors":[{"message":"Cannot query field"}]}')
  })

  res <- gnomad_frequencies(c("1-1-A-T", "2-2-C-G"))

  expect_true(res$ok)
  out <- biohttp::body_or_null(res)
  expect_identical(nrow(out), 2L)
  expect_true(all(is.na(out$exome_af)))
})
