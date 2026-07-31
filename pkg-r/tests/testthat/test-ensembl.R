# --- The ported fixture, which must pass unchanged ---------------------------

test_that("the ported Ensembl VEP fixture parses as expected", {
  body <- read_fixture("ensembl_vep_rs113488022.json")
  out <- ensembl_parse_vep(body)

  expect_identical(out$most_severe, "missense_variant")
  expect_identical(out$assembly, "GRCh38")
  expect_s3_class(out$consequences, "tbl_df")
  expect_identical(out$consequences$gene[1], "BRAF")
  expect_identical(out$consequences$transcript[1], "ENST00000288602")
  expect_identical(out$consequences$impact[1], "MODERATE")
})

test_that("the SIFT and PolyPhen calls come through", {
  body <- read_fixture("ensembl_vep_rs113488022.json")
  out <- ensembl_parse_vep(body)$consequences

  expect_identical(out$sift[1], "deleterious_low_confidence")
  expect_identical(out$polyphen[1], "possibly_damaging")
})

# --- Every biotype is returned -----------------------------------------------

test_that("non-coding transcripts are kept, with their biotype", {
  # variant-reviewer filters to protein_coding as "the rows worth showing".
  # Which biotypes are worth showing is the caller's call, and filtering here
  # would hide them from a caller that wanted them.
  consequences <- list(
    list(transcript_id = "ENST1", biotype = "protein_coding"),
    list(transcript_id = "ENST2", biotype = "retained_intron")
  )
  out <- ensembl_parse_consequences(consequences)

  expect_identical(nrow(out), 2L)
  expect_identical(out$biotype, c("protein_coding", "retained_intron"))
})

test_that("several consequence terms on one transcript are all kept", {
  # Taking the first drops the rest with no sign it happened.
  consequences <- list(list(
    transcript_id = "ENST1",
    consequence_terms = list("splice_region_variant", "intron_variant")
  ))
  out <- ensembl_parse_consequences(consequences)

  expect_identical(out$consequence, "splice_region_variant, intron_variant")
})

test_that("no consequences is NULL", {
  expect_null(ensembl_parse_consequences(list()))
  expect_null(ensembl_parse_consequences(NULL))
  expect_null(ensembl_parse_vep(list()))
})

# --- Trap: exon numbering follows the strand ---------------------------------

test_that("minus-strand exons are numbered from the highest coordinate", {
  # Exon 1 is the first one transcribed. On the minus strand that is the
  # highest genomic coordinate. Numbering by coordinate alone puts exon 1 at
  # the wrong end of every minus-strand gene, and the result looks fine.
  record <- list(
    seq_region_name = "17",
    Transcript = list(list(
      id = "ENST1",
      is_canonical = 1,
      strand = -1,
      Exon = list(
        list(start = 800, end = 900),
        list(start = 400, end = 500),
        list(start = 1, end = 100)
      )
    ))
  )
  out <- ensembl_parse_gene_model(record)

  expect_identical(out$exons$start, c(1, 400, 800))
  expect_identical(out$exons$number, c(3L, 2L, 1L))
})

test_that("plus-strand exons are numbered from the lowest coordinate", {
  record <- list(
    Transcript = list(list(
      id = "ENST1",
      is_canonical = 1,
      strand = 1,
      Exon = list(list(start = 800, end = 900), list(start = 1, end = 100))
    ))
  )
  out <- ensembl_parse_gene_model(record)

  expect_identical(out$exons$start, c(1, 800))
  expect_identical(out$exons$number, c(1L, 2L))
})

# --- Transcript choice -------------------------------------------------------

test_that("the canonical transcript wins", {
  record <- list(
    Transcript = list(
      list(
        id = "ENST_other",
        strand = 1,
        Exon = list(list(start = 1, end = 2))
      ),
      list(
        id = "ENST_canonical",
        is_canonical = 1,
        strand = 1,
        Exon = list(list(start = 5, end = 6))
      )
    )
  )
  expect_identical(
    ensembl_parse_gene_model(record)$transcript,
    "ENST_canonical"
  )
})

test_that("with no canonical flag the first transcript is used", {
  record <- list(
    Transcript = list(
      list(id = "ENST_first", strand = 1, Exon = list(list(start = 1, end = 2)))
    )
  )
  expect_identical(ensembl_parse_gene_model(record)$transcript, "ENST_first")
})

test_that("a record with no transcripts or no exons is NULL", {
  expect_null(ensembl_parse_gene_model(list()))
  expect_null(ensembl_parse_gene_model(list(Transcript = list())))
  expect_null(ensembl_parse_gene_model(list(
    Transcript = list(
      list(id = "ENST1", Exon = list())
    )
  )))
})

test_that("an exon with no coordinates is dropped", {
  record <- list(
    Transcript = list(list(
      id = "ENST1",
      strand = 1,
      Exon = list(list(start = 1, end = 100), list(end = 200))
    ))
  )
  expect_identical(nrow(ensembl_parse_gene_model(record)$exons), 1L)
})

# --- The response shape ------------------------------------------------------

test_that("an array body and a single-record body both work", {
  # The endpoint answers with an array, but the stored response is a bare
  # record. Taking [[1]] blindly on a record returns its first FIELD, which
  # parses to a row of NA rather than to an error.
  record <- read_fixture("ensembl_vep_rs113488022.json")

  expect_identical(ensembl_first_record(list(record)), record)
  expect_identical(ensembl_first_record(record), record)
  expect_null(ensembl_first_record(list()))
})

test_that("the client reads a single-record body correctly", {
  reset_transport()
  body <- paste(
    '{"most_severe_consequence":"missense_variant","assembly_name":"GRCh38",',
    '"transcript_consequences":[{"gene_symbol":"BRAF"}]}',
    sep = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(body))

  out <- biohttp::body_or_null(ensembl_vep_id("rs113488022"))

  expect_identical(out$most_severe, "missense_variant")
  expect_identical(out$consequences$gene[1], "BRAF")
})

# --- The client half ---------------------------------------------------------

test_that("expand=1 is sent, or the model comes back empty", {
  # Without it Ensembl returns a bare gene record with no Transcript array.
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json('{"Transcript":[]}')
  })

  ensembl_gene_model("ENSG00000157764")

  expect_match(url, "lookup/id/ENSG00000157764", fixed = TRUE)
  expect_match(url, "expand=1", fixed = TRUE)
})

test_that("the variant id goes into the VEP path", {
  reset_transport()
  url <- NULL
  httr2::local_mocked_responses(function(req) {
    url <<- req$url
    mock_json("[]")
  })

  ensembl_vep_id("rs113488022")

  expect_match(url, "vep/human/id/rs113488022", fixed = TRUE)
})

test_that("a malformed id never reaches the network", {
  reset_transport()
  expect_identical(ensembl_vep_id("rs1/../../x")$status, "no_data")
  expect_identical(ensembl_vep_id("")$status, "no_data")
  expect_identical(ensembl_gene_model("BRAF")$status, "no_data")
  expect_identical(ensembl_gene_model(NULL)$status, "no_data")
})
