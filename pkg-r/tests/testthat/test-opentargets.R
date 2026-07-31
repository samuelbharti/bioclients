# --- The ported fixtures, which must pass unchanged --------------------------

test_that("the ported target-to-disease fixture parses as expected", {
  body <- read_fixture("opentargets_tp53.json")
  out <- opentargets_parse_diseases(body, "ENSG00000141510")

  expect_s3_class(out, "tbl_df")
  expect_identical(nrow(out), 3L)
  expect_identical(out$disease[1], "Li-Fraumeni syndrome")
  expect_identical(out$disease_id[1], "MONDO_0018875")
  expect_equal(out$score[1], 0.8763216350824885)
  expect_identical(
    out$source_url[1],
    "https://platform.opentargets.org/evidence/ENSG00000141510/MONDO_0018875"
  )
})

test_that("association rows keep Open Targets' own score order", {
  # Open Targets returns them sorted, best first. Re-sorting here would be this
  # package deciding what "best" means, which is the app's call, not ours.
  body <- read_fixture("opentargets_tp53.json")
  out <- opentargets_parse_diseases(body, "ENSG00000141510")
  expect_true(all(diff(out$score) <= 0))
})

test_that("the ported disease-search fixture parses as expected", {
  body <- read_fixture("ot_disease_search_nf1.json")
  out <- opentargets_parse_matches(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$id[1], "MONDO_0018975")
  expect_identical(out$name[1], "neurofibromatosis type 1")
  expect_equal(out$score[1], 24.871593)
  expect_match(out$description[1], "tumor predisposition")
  expect_identical(
    out$source_url[1],
    "https://platform.opentargets.org/disease/MONDO_0018975"
  )
})

# --- One parser reads both response shapes -----------------------------------

test_that("the disease parser reads the single-lookup shape too", {
  # A caller resolving a term does not know whether a search or a direct lookup
  # ran, so both shapes have to land in the same table or the difference leaks
  # out into every caller.
  body <- list(
    data = list(
      disease = list(
        id = "EFO_0000508",
        name = "neurofibroma",
        description = "A benign nerve-sheath tumor."
      )
    )
  )
  out <- opentargets_parse_matches(body)

  expect_identical(nrow(out), 1L)
  expect_identical(out$id, "EFO_0000508")
  # A direct lookup has no relevance score to report.
  expect_true(is.na(out$score))
})

test_that("an empty search and an absent disease both parse to NULL", {
  expect_null(opentargets_parse_matches(list(
    data = list(search = list(hits = list()))
  )))
  expect_null(opentargets_parse_matches(list(data = list(disease = NULL))))
})

# --- Id handling -------------------------------------------------------------

test_that("ontology ids are told apart from free text", {
  expect_true(opentargets_is_id("MONDO:0018975"))
  expect_true(opentargets_is_id("EFO_0000508"))
  expect_true(opentargets_is_id("Orphanet_636"))
  expect_false(opentargets_is_id("neurofibromatosis type 1"))
  expect_false(opentargets_is_id(""))
  expect_false(opentargets_is_id(NULL))
})

test_that("only the colon separator is swapped, not the case", {
  # Open Targets ids are underscore-separated, but some ontologies keep mixed
  # case in the body of the id. Upcasing "Orphanet_636" stops it resolving.
  expect_identical(opentargets_normalize_id("MONDO:0018975"), "MONDO_0018975")
  expect_identical(opentargets_normalize_id("Orphanet:636"), "Orphanet_636")
  expect_identical(opentargets_normalize_id("Orphanet_636"), "Orphanet_636")
})

# --- The disease-to-target direction -----------------------------------------

test_that("disease-to-target rows parse into a gene table", {
  body <- list(
    data = list(
      disease = list(
        name = "neurofibromatosis type 1",
        associatedTargets = list(
          count = 2,
          rows = list(
            list(
              score = 0.9,
              target = list(
                id = "ENSG00000196712",
                approvedSymbol = "NF1"
              )
            ),
            list(
              score = 0.4,
              target = list(
                id = "ENSG00000133703",
                approvedSymbol = "KRAS"
              )
            )
          )
        )
      )
    )
  )
  out <- opentargets_parse_targets(body, "MONDO_0018975")

  expect_identical(out$symbol, c("NF1", "KRAS"))
  expect_identical(out$ensembl_id[1], "ENSG00000196712")
  expect_equal(out$score, c(0.9, 0.4))
  expect_match(out$source_url[1], "ENSG00000196712/MONDO_0018975")
})

test_that("an absent target or disease parses to NULL", {
  expect_null(opentargets_parse_diseases(list(data = list(target = NULL))))
  expect_null(opentargets_parse_targets(list(data = list(disease = NULL))))
})

# --- The client half ---------------------------------------------------------

test_that("opentargets_gene_diseases returns an ok envelope with the table", {
  reset_transport()
  fixture <- paste(
    readLines(
      testthat::test_path("fixtures", "opentargets_tp53.json"),
      warn = FALSE
    ),
    collapse = ""
  )
  httr2::local_mocked_responses(function(req) mock_json(fixture))

  res <- opentargets_gene_diseases("ENSG00000141510")

  expect_true(res$ok)
  expect_identical(res$source, "Open Targets")
  expect_identical(
    biohttp::body_or_null(res)$disease[1],
    "Li-Fraumeni syndrome"
  )
})

test_that("a GraphQL errors array inside a 200 is a failure", {
  reset_transport()
  httr2::local_mocked_responses(function(req) {
    mock_json('{"errors":[{"message":"Cannot query field"}]}')
  })
  res <- opentargets_gene_diseases("ENSG00000141510")

  expect_false(res$ok)
  expect_identical(res$status, "error")
})

test_that("resolve_disease looks an id up directly instead of searching", {
  # A search for an exact id wastes a relevance ranking on something already
  # unambiguous, and can rank the exact record below a fuzzy neighbour.
  reset_transport()
  seen <- NULL
  httr2::local_mocked_responses(function(req) {
    seen <<- req$body$data
    mock_json('{"data":{"disease":{"id":"MONDO_0018975","name":"NF1"}}}')
  })

  res <- opentargets_resolve_disease("MONDO:0018975")

  expect_true(res$ok)
  expect_true("efoId" %in% names(seen$variables))
  expect_false("queryString" %in% names(seen$variables))
  # The colon form was normalized before it went out.
  expect_identical(seen$variables$efoId, "MONDO_0018975")
})

test_that("resolve_disease searches free text", {
  reset_transport()
  seen <- NULL
  httr2::local_mocked_responses(function(req) {
    seen <<- req$body$data
    mock_json(
      '{"data":{"search":{"hits":[{"id":"MONDO_0018975","name":"NF1"}]}}}'
    )
  })

  res <- opentargets_resolve_disease("neurofibromatosis type 1")

  expect_true(res$ok)
  expect_true("q" %in% names(seen$variables))
  expect_false("efoId" %in% names(seen$variables))
  expect_match(seen$query, "queryString", fixed = TRUE)
})

test_that("blank inputs are no_data and never reach the network", {
  reset_transport()
  expect_identical(opentargets_gene_diseases("")$status, "no_data")
  expect_identical(opentargets_disease_targets("")$status, "no_data")
  expect_identical(opentargets_resolve_disease("  ")$status, "no_data")
})

# --- Known drugs -------------------------------------------------------------

test_that("the ported known-drugs fixture parses as expected", {
  body <- read_fixture("opentargets_drugs_braf.json")
  out <- opentargets_parse_drugs(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$drug[1], "BELVARAFENIB")
  expect_identical(out$drug_id[1], "CHEMBL3977543")
  expect_identical(out$drug_type[1], "Small molecule")
  expect_identical(out$max_phase[1], "PHASE_2")
})

test_that("every disease a drug was tried against is kept", {
  # Taking the first would report one of five here, and the first is the one
  # Open Targets could not map at all.
  body <- read_fixture("opentargets_drugs_braf.json")
  out <- opentargets_parse_drugs(body)

  expect_type(out$diseases, "list")
  expect_length(out$diseases[[1]], 5)
  expect_true("melanoma" %in% out$diseases[[1]])
})

test_that("an unmapped disease falls back to the source label", {
  # disease is null on the first entry; diseaseFromSource is all there is.
  body <- read_fixture("opentargets_drugs_braf.json")
  out <- opentargets_parse_drugs(body)

  expect_identical(out$diseases[[1]][1], "extracranial solid tumours")
  expect_true(is.na(out$disease_ids[[1]][1]))
})

test_that("the clinical stage is not reformatted for display", {
  # "PHASE_2" is what Open Targets sends. Turning it into "Phase 2" is
  # presentation, and presentation belongs to whatever is presenting.
  body <- read_fixture("opentargets_drugs_braf.json")

  expect_true(all(grepl(
    "^[A-Z_0-9]+$",
    opentargets_parse_drugs(body)$max_phase
  )))
})

test_that("a target with no known drugs is NULL", {
  expect_null(opentargets_parse_drugs(list(
    data = list(
      target = list(
        drugAndClinicalCandidates = list(count = 0, rows = list())
      )
    )
  )))
  expect_null(opentargets_parse_drugs(list()))
})

# --- Pharmacogenomics --------------------------------------------------------

test_that("the ported pharmacogenomics fixture parses as expected", {
  body <- read_fixture("opentargets_pgx_cyp2c19.json")
  out <- opentargets_parse_pgx(body)

  expect_s3_class(out, "tbl_df")
  expect_identical(out$phenotype[1], "decreased metabolism of venlafaxine")
  expect_identical(out$drugs[[1]], "venlafaxine")
})

test_that("a genotype-keyed annotation has no rsid, and that is not an error", {
  # variantRsId is null whenever the annotation is keyed on a genotype instead.
  body <- read_fixture("opentargets_pgx_cyp2c19.json")

  expect_true(is.na(opentargets_parse_pgx(body)$rsid[1]))
})

test_that("several drugs on one annotation stay separate", {
  body <- list(
    data = list(
      target = list(
        pharmacogenomics = list(
          list(
            drugs = list(
              list(drugFromSource = "venlafaxine"),
              list(drugFromSource = "citalopram"),
              list(drugFromSource = "venlafaxine")
            )
          )
        )
      )
    )
  )
  out <- opentargets_parse_pgx(body)

  expect_identical(out$drugs[[1]], c("venlafaxine", "citalopram"))
})

test_that("a gene with no pharmacogenomics is NULL", {
  # The normal case for most genes, not a failure.
  expect_null(opentargets_parse_pgx(list(
    data = list(
      target = list(
        pharmacogenomics = list()
      )
    )
  )))
})

# --- The client halves -------------------------------------------------------

test_that("both new queries post the Ensembl id as a variable", {
  reset_transport()
  seen <- list()
  httr2::local_mocked_responses(function(req) {
    seen[[length(seen) + 1]] <<- req$body$data
    mock_json('{"data":{"target":null}}')
  })

  opentargets_drugs("ENSG00000157764")
  opentargets_pgx("ENSG00000165841")

  expect_identical(seen[[1]]$variables$id, "ENSG00000157764")
  expect_match(seen[[1]]$query, "drugAndClinicalCandidates", fixed = TRUE)
  expect_identical(seen[[2]]$variables$id, "ENSG00000165841")
  expect_match(seen[[2]]$query, "pharmacogenomics", fixed = TRUE)
})

test_that("a blank id never reaches the network for either", {
  reset_transport()
  expect_identical(opentargets_drugs("")$status, "no_data")
  expect_identical(opentargets_pgx(NULL)$status, "no_data")
})
