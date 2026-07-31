# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported Allele Registry fixture parses as expected", {
  body <- read_fixture("clingen_alleles_batch.json")
  out <- clingen_parse_batch(body)

  expect_s3_class(out, "tbl_df")
  expect_true(out$resolved[1])
  expect_identical(out$caid[1], "CA000072")
  expect_identical(out$clinvar_allele_id[1], 27390L)
  expect_match(out$title[1], "TP53")
})

# --- An unresolved input keeps its row ---------------------------------------

test_that("an input ClinGen could not resolve becomes a row with a reason", {
  # ClinGen answers with an object carrying no @id and an errorType instead.
  # Dropping it would shift every later row onto the wrong variant and lose the
  # reason the input was rejected.
  out <- clingen_parse_allele(list(
    errorType = "IncorrectHgvsPosition",
    message = "position does not match reference"
  ))

  expect_false(out$resolved)
  expect_identical(out$reason, "position does not match reference")
  expect_true(is.na(out$caid))
})

test_that("errorType is used when there is no message", {
  out <- clingen_parse_allele(list(errorType = "NoConsistentAlignment"))
  expect_identical(out$reason, "NoConsistentAlignment")
})

test_that("resolved and unresolved rows sit in one table", {
  body <- list(
    list(`@id` = "http://reg.genome.network/allele/CA000072"),
    list(errorType = "IncorrectHgvsPosition")
  )
  out <- clingen_parse_batch(body)

  expect_identical(nrow(out), 2L)
  expect_identical(out$resolved, c(TRUE, FALSE))
})

# --- The GRCh38 locus --------------------------------------------------------

test_that("only a chromosomal GRCh38 allele is read as the locus", {
  # Without the NC_ check this could pick up a scaffold or an older build and
  # report coordinates that do not mean what the caller assumes.
  element <- list(
    `@id` = "http://x/allele/CA1",
    genomicAlleles = list(
      list(
        referenceGenome = "GRCh37",
        chromosome = "17",
        hgvs = list("NC_000017.10:g.1A>T"),
        coordinates = list(list(start = 1, referenceAllele = "A", allele = "T"))
      ),
      list(
        referenceGenome = "GRCh38",
        chromosome = "17",
        hgvs = list("NC_000017.11:g.7676154G>C"),
        coordinates = list(list(
          start = 7676154,
          referenceAllele = "G",
          allele = "C"
        ))
      )
    )
  )
  out <- clingen_parse_allele(element)

  expect_identical(out$chrom, "17")
  expect_identical(out$pos, 7676154L)
  expect_identical(out$ref, "G")
})

test_that("a non-chromosomal GRCh38 allele is not used", {
  element <- list(
    `@id` = "http://x/allele/CA1",
    genomicAlleles = list(list(
      referenceGenome = "GRCh38",
      chromosome = "17",
      hgvs = list("NW_012345.1:g.1A>T"),
      coordinates = list(list(start = 1, referenceAllele = "A", allele = "T"))
    ))
  )
  expect_true(is.na(clingen_parse_allele(element)$pos))
})

# --- The client half ---------------------------------------------------------

test_that("the batch POST needs no key and sends newline-delimited text", {
  # Most identity services of this kind require registration. This one does not,
  # which is worth not forgetting.
  reset_transport()
  sent <- NULL
  url <- NULL
  headers <- NULL
  httr2::local_mocked_responses(function(req) {
    body <- req$body$data
    sent <<- if (is.raw(body)) rawToChar(body) else as.character(body)
    url <<- req$url
    headers <<- names(req$headers)
    mock_json('[{"@id":"http://x/allele/CA1"}]')
  })

  clingen_alleles(c("NM_1:c.1A>T", "NM_2:c.2C>G"))

  expect_identical(sent, "NM_1:c.1A>T\nNM_2:c.2C>G")
  expect_match(url, "file=hgvs", fixed = TRUE)
  # No credential of any kind goes out.
  expect_false(any(grepl("Authorization", headers, fixed = TRUE)))
  expect_false(any(grepl("api", tolower(headers), fixed = TRUE)))
})

test_that("the request still goes through biohttp's transport", {
  # This client assembles its own request because biohttp has no text/plain POST
  # wrapper. It must still perform through biohttp, or it loses the retries,
  # the breaker, and the envelope. An open breaker proves it does.
  reset_transport()
  for (i in seq_len(5)) {
    biohttp::breaker_record("reg.clinicalgenome.org", reachable = FALSE)
  }
  # No mock installed, so a dispatched request would hit the network.
  res <- clingen_alleles("NM_1:c.1A>T")

  expect_identical(res$status, "skipped")
  reset_transport()
})

test_that("no HGVS is no_data and never reaches the network", {
  reset_transport()
  expect_identical(clingen_alleles(character())$status, "no_data")
  expect_identical(clingen_alleles(c(NA, ""))$status, "no_data")
})
