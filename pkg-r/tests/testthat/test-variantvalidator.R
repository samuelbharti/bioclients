# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported VariantValidator fixture parses as expected", {
  body <- read_fixture("variantvalidator_tp53.json")
  out <- variantvalidator_parse(body, "NM_000546.6:c.215C>G")

  expect_s3_class(out, "tbl_df")
  expect_true(out$resolved)
  expect_identical(out$gene, "TP53")
  expect_type(out$warnings, "list")
})

# --- The response is keyed by the answer -------------------------------------

test_that("the record is found by excluding the fixed keys", {
  # VariantValidator keys the object by the NORMALIZED variant, which is not
  # known before the call. flag and metadata are the only fixed keys. A parser
  # expecting a fixed name finds nothing.
  body <- list(
    flag = "gene_variant",
    metadata = list(version = "2.2"),
    `NM_000546.6:c.215C>G` = list(gene_symbol = "TP53")
  )
  out <- variantvalidator_parse(body, "submitted")

  expect_identical(out$normalized, "NM_000546.6:c.215C>G")
  expect_identical(out$gene, "TP53")
  expect_identical(out$submitted, "submitted")
})

test_that("a response with only the fixed keys parses to NULL", {
  expect_null(variantvalidator_parse(list(flag = "warning", metadata = list())))
  expect_null(variantvalidator_parse(list()))
})

# --- Warnings are the rejection reasons --------------------------------------

test_that("validation warnings are kept, not collapsed to a pass or fail", {
  # They are VariantValidator's real rejection reasons, and a bare FALSE throws
  # away the only thing that tells a user what to fix.
  body <- list(
    flag = "warning",
    `NM_1:c.1A>T` = list(
      gene_symbol = "X",
      validation_warnings = list("position does not agree", "using RefSeq")
    )
  )
  out <- variantvalidator_parse(body, "NM_1:c.1A>T")

  expect_length(out$warnings[[1]], 2)
  expect_match(out$warnings[[1]][1], "position does not agree")
})

test_that("no warnings is an empty vector, not NULL", {
  body <- list(flag = "ok", `NM_1:c.1A>T` = list(gene_symbol = "X"))
  expect_identical(variantvalidator_parse(body)$warnings[[1]], character())
})

# --- The GRCh38 locus --------------------------------------------------------

test_that("the GRCh38 VCF locus is pulled out", {
  body <- list(
    flag = "ok",
    `NM_1:c.1A>T` = list(
      gene_symbol = "X",
      primary_assembly_loci = list(
        grch38 = list(
          vcf = list(
            chr = "chr17",
            pos = "7676154",
            ref = "G",
            alt = "C"
          )
        )
      )
    )
  )
  out <- variantvalidator_parse(body)

  expect_identical(out$chrom, "17")
  expect_identical(out$pos, 7676154L)
  expect_identical(out$ref, "G")
  expect_identical(out$alt, "C")
})

# --- The client half ---------------------------------------------------------

test_that("the documented rate limit is the default", {
  expect_identical(VARIANTVALIDATOR_THROTTLE$capacity, 4)
  expect_identical(VARIANTVALIDATOR_THROTTLE$fill_time_s, 1)
})

test_that("a blank HGVS is no_data and never reaches the network", {
  reset_transport()
  expect_identical(variantvalidator_normalize("")$status, "no_data")
  expect_identical(variantvalidator_normalize(NULL)$status, "no_data")
})
