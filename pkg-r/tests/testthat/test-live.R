# Live service checks.
#
# Every parser in this package is pinned by a stored response, but those
# responses were ported out of the apps rather than fetched. A fixture cannot
# tell you that a service changed its mind. This is the only file that asks a
# real server, and it is the only file that can answer that question.
#
# It never runs on its own:
#
#   BIOCLIENTS_LIVE=true Rscript -e 'devtools::test("pkg-r", filter = "live")'
#
# Three separate gates, because they fail for different reasons and a reader
# should be able to tell which one stopped the run. CRAN forbids network access
# outright, a laptop with no network is not a regression, and even with a
# network this should never run by accident. No CI job sets the variable, so the
# file skips on every runner regardless of what else changes.
#
# READING A FAILURE.
#
# A 5xx means the service is down. That says nothing about whether the port is
# correct, so treat it as an outage and run again later. What matters here is a
# column that comes back entirely NA, because that is what a field renamed
# upstream looks like from inside a parser that is still doing its job.

skip_if_not(
  nzchar(Sys.getenv("BIOCLIENTS_LIVE")),
  "live service checks: set BIOCLIENTS_LIVE=true to run"
)
skip_on_cran()
skip_if_offline()

# --- Helpers -----------------------------------------------------------------

# A live call is allowed to find nothing. It is not allowed to fault.
expect_answered <- function(res) {
  if (res$status %in% c("ok", "no_data")) {
    return(succeed())
  }
  fail(sprintf(
    "%s answered %s (HTTP %s): %s",
    res$source,
    res$status,
    format(res$http),
    res$detail %||% "no detail"
  ))
}

# The parsers build their tibbles field by field, so the column set cannot drift
# no matter what a service returns. Asserting on column names would only be
# testing this package against itself. What actually drifts is a field renamed
# at the source, which leaves the column in place and fills it with NA, so that
# is what this looks for.
expect_populated <- function(out, cols, source) {
  # Checked before anything touches nrow(), so a service that answered with
  # nothing reports that in one line instead of an unreadable pile of errors.
  if (!is.data.frame(out) || nrow(out) == 0L) {
    return(fail(sprintf(
      "%s: expected rows, got %s",
      source,
      if (is.null(out)) {
        "NULL, so the call came back no_data"
      } else {
        paste0("a ", class(out)[1], " of ", NROW(out), " rows")
      }
    )))
  }

  absent <- setdiff(cols, names(out))
  if (length(absent)) {
    return(fail(sprintf(
      "%s: the parser produced no column named %s",
      source,
      paste(absent, collapse = ", ")
    )))
  }

  blank <- cols[vapply(cols, function(x) all(is.na(out[[x]])), logical(1))]
  if (length(blank)) {
    return(fail(sprintf(
      "%s answered, but %s came back entirely NA, which is what a renamed field upstream looks like",
      source,
      paste(blank, collapse = ", ")
    )))
  }
  succeed()
}

# --- One real call per service -----------------------------------------------
#
# The subject of each probe is the subject of that service's stored fixture, so
# the live call and the offline test ask the same question of the same gene. A
# disagreement between them is then worth something.

test_that("live: MyGene answers with a populated gene record", {
  res <- mygene_gene("TP53")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("symbol", "entrez", "ensembl_gene", "uniprot", "hgnc"),
    "MyGene"
  )
})

test_that("live: gnomAD answers with gene constraint", {
  res <- gnomad_constraint("BRAF")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("symbol", "pli", "loeuf"),
    "gnomAD"
  )
})

test_that("live: ClinVar answers with a classification", {
  res <- clinvar_classification("rs113488022")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("uid", "accession", "title", "significance"),
    "ClinVar"
  )
})

test_that("live: MyVariant answers with a variant record", {
  res <- myvariant_variants("chr17:g.7676154G>C")
  expect_answered(res)
  expect_populated(biohttp::body_or_null(res), c("id", "gene"), "MyVariant")
})

test_that("live: Ensembl VEP answers for a region", {
  res <- vep_variants("7", 140753336, "A", "T")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("key", "gene", "consequence"),
    "Ensembl VEP"
  )
})

test_that("live: the Ensembl VEP id route answers", {
  res <- ensembl_vep_id("rs113488022")
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  # Not a tibble. A list of the two scalars and the consequences table.
  expect_type(out, "list")
  expect_false(is.na(out$most_severe))
  expect_populated(out$consequences, c("gene", "consequence"), "Ensembl VEP id")
})

test_that("live: VariantValidator normalises an HGVS description", {
  res <- variantvalidator_normalize("NM_000546.6:c.215C>G")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("resolved", "gene", "chrom", "pos", "ref", "alt"),
    "VariantValidator"
  )
})

test_that("live: the ClinGen Allele Registry resolves an HGVS description", {
  res <- clingen_alleles("NM_000546.6:c.215C>G")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("resolved", "caid"),
    "ClinGen Allele Registry"
  )
})

test_that("live: Open Targets answers with gene to disease associations", {
  res <- opentargets_gene_diseases("ENSG00000141510")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("disease", "disease_id", "score"),
    "Open Targets"
  )
})

test_that("live: DGIdb answers with an interaction count", {
  res <- dgidb_gene("NF1")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("symbol", "concept_id", "interaction_count"),
    "DGIdb"
  )
})

test_that("live: Pharos answers with a target development level", {
  res <- pharos_target("NF1")
  expect_answered(res)
  expect_populated(biohttp::body_or_null(res), c("symbol", "tdl"), "Pharos")
})

test_that("live: CIViC answers with a gene record", {
  res <- civic_gene("NF1")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("id", "symbol", "entrez_id"),
    "CIViC"
  )
})

test_that("live: JensenLab DISEASES answers on the Knowledge channel", {
  # DOID:1612, breast cancer. The fixture's own subject, DOID:0060293, now comes
  # back empty on both channels: the service dropped those annotations, and a
  # probe that cannot tell an empty answer from a broken one is not worth
  # running.
  res <- diseases_channel("DOID:1612", "Knowledge")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("symbol", "protein", "score", "channel"),
    "DISEASES"
  )
})

test_that("live: PanelApp answers with a page of the panel index", {
  res <- panelapp_panels()
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  expect_type(out, "list")
  expect_populated(out$panels, c("id", "name", "version"), "PanelApp")
})

test_that("live: ClinGen gene validity answers with the classification table", {
  res <- clingen_gene_validity()
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("gene", "disease", "classification"),
    "ClinGen gene validity"
  )
})

test_that("live: UniProt answers with disease involvement", {
  res <- uniprot_diseases("P21359")
  expect_answered(res)
  expect_populated(biohttp::body_or_null(res), c("id", "name"), "UniProt")
})

test_that("live: ProtVar answers with a function description", {
  res <- protvar_function("P04637", 175)
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  # A single description string, not a table.
  expect_type(out, "character")
  expect_length(out, 1L)
  expect_true(nzchar(out))
})

test_that("live: AlphaFold answers with a model", {
  res <- alphafold_model("P15056")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("accession", "pdb_url", "cif_url"),
    "AlphaFold"
  )
})

test_that("live: PDBe answers with experimental structures", {
  res <- pdbe_structures("P21359")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("pdb_id", "method"),
    "PDBe"
  )
})

test_that("live: STRING answers with interaction partners", {
  res <- string_partners("TP53")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("partner", "score"),
    "STRING"
  )
})

test_that("live: GTEx answers with its own versioned GENCODE id", {
  res <- gtex_gene_reference("TP53")
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  expect_populated(out, c("symbol", "gencode_id"), "GTEx")
  # The versioned form is the whole point of this lookup. An unversioned id
  # would parse and then fail every expression call made with it.
  expect_match(out$gencode_id[[1]], "^ENSG[0-9]+\\.[0-9]+$")
})

test_that("live: the Human Protein Atlas answers for an Ensembl id", {
  res <- hpa_gene("ENSG00000141510")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("symbol", "ensembl"),
    "Human Protein Atlas"
  )
})

test_that("live: QuickGO answers with GO annotations", {
  res <- quickgo_annotations("P21359")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("go_id", "go_name", "aspect"),
    "QuickGO"
  )
})

test_that("live: Reactome answers with pathways", {
  res <- reactome_pathways("NF1")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("pathway_id", "name"),
    "Reactome"
  )
})

test_that("live: the HPO answers a term search", {
  res <- hpo_search("seizure")
  expect_answered(res)
  expect_populated(biohttp::body_or_null(res), c("id", "name"), "HPO")
})

test_that("live: Monarch answers a search", {
  res <- monarch_search("Marfan syndrome")
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  expect_type(out, "list")
  expect_populated(out$matches, c("id", "name", "category"), "Monarch")
  # Monarch's own count, so a caller can say what was left behind rather than
  # presenting the first page as everything.
  expect_gte(out$total, nrow(out$matches))
})

test_that("live: IMPC answers with a mouse ortholog", {
  res <- impc_mouse_ortholog("NF1")
  expect_answered(res)
  expect_populated(
    biohttp::body_or_null(res),
    c("mgi", "marker_symbol"),
    "IMPC"
  )
})

test_that("live: Europe PMC answers a literature search", {
  res <- europepmc_search("NF1")
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  expect_type(out, "list")
  expect_populated(out$results, c("title", "source_id"), "Europe PMC")
  expect_gte(out$count, nrow(out$results))
})

test_that("live: PubTator3 answers with a literature count and results", {
  res <- pubtator_gene_literature("TP53", 7157)
  expect_answered(res)
  out <- biohttp::body_or_null(res)
  expect_type(out, "list")
  expect_gt(out$count, 0L)
  expect_populated(out$results, c("pmid", "title"), "PubTator3")
})

# --- The five behaviours only a live call can confirm ------------------------
#
# Each of these is a claim written as a comment at a call site, carried across
# from an app during the port. Where the claim is about the response body, this
# asserts against the raw body rather than the parsed result, because a parser
# that had quietly adapted to a change would hide exactly what is being asked.

test_that("live: MyGene names the HGNC id in upper case", {
  # Asking for `hgnc` returns nothing at all, and nothing reads as "this gene
  # has no HGNC id" rather than as a mistake. See the note in R/mygene.R.
  res <- biohttp::get_json(
    MYGENE_BASE,
    path = "query",
    query = list(
      q = "TP53",
      species = "human",
      fields = "symbol,HGNC,entrezgene",
      size = 1
    ),
    source = "MyGene"
  )
  expect_true(isTRUE(res$ok))

  hit <- res$data$hits[[1]]
  expect_true("HGNC" %in% names(hit))
  expect_false("hgnc" %in% names(hit))
  expect_true(nzchar(as.character(hit$HGNC)))
})

test_that("live: Reactome answers 404 for a gene it maps with no pathways", {
  # The client special-cases the status code, so if Reactome ever switches to an
  # empty 200 that branch goes dead without anyone noticing. ZNF876P is a gene
  # Reactome knows and has never annotated.
  raw <- biohttp::get_json(
    REACTOME_URL,
    path = "data/mapping/HGNC/ZNF876P/pathways",
    source = "Reactome"
  )
  expect_identical(raw$http, 404L)

  # And the client turns that into an answer rather than a fault.
  res <- reactome_pathways("ZNF876P")
  expect_identical(res$status, "no_data")
  expect_identical(res$http, 404L)
})

test_that("live: the PanelApp index search parameter does not filter", {
  # `panels/` accepts `search`, answers 200, and returns the unfiltered index.
  # This is the whole reason panelapp_all_panels() exists. A day when this test
  # fails is a day the client gets simpler.
  res <- biohttp::get_json(
    PANELAPP_URL,
    path = "panels/",
    query = list(search = "Marfan", page_size = 10),
    source = "PanelApp"
  )
  expect_true(isTRUE(res$ok))

  names_returned <- vapply(res$data$results, function(p) p$name, character(1))
  expect_gt(length(names_returned), 0L)
  # Not one of them matches, and the reported total is the whole index.
  expect_false(any(grepl("Marfan", names_returned, ignore.case = TRUE)))
  expect_gt(res$data$count, length(names_returned))
})

test_that("live: the Ensembl VEP id route answers with an array, not a record", {
  # ensembl_first_record() exists because of this. Taking [[1]] of a body that
  # was already a single record would return its first field instead.
  res <- biohttp::get_json(
    ENSEMBL_URL,
    path = "vep/human/id/rs113488022",
    query = list(`content-type` = "application/json"),
    source = "Ensembl VEP"
  )
  expect_true(isTRUE(res$ok))

  expect_type(res$data, "list")
  expect_null(names(res$data))
  expect_gte(length(res$data), 1L)
  expect_true("most_severe_consequence" %in% names(res$data[[1]]))
})

test_that("live: every Monarch route answers on the one host", {
  # This is what let the second host go. See the note at the top of R/monarch.R.
  routes <- list(
    search = list(
      path = "search",
      query = list(q = "Marfan syndrome", limit = 2)
    ),
    association = list(
      path = "association",
      query = list(subject = "HGNC:11998", limit = 2)
    ),
    entity = list(
      path = paste(
        "entity",
        "HGNC:11998",
        MONARCH_GENE_PHENOTYPE,
        sep = "/"
      ),
      query = list(limit = 2)
    )
  )

  for (name in names(routes)) {
    res <- biohttp::get_json(
      MONARCH_URL,
      path = routes[[name]]$path,
      query = routes[[name]]$query,
      source = "Monarch"
    )
    expect_identical(res$http, 200L)
    expect_true(isTRUE(res$ok))
  }
})
