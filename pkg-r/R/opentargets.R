# Open Targets Platform: gene and disease association.
#
# Ported from three places that all query the same endpoint:
#
#   variant-reviewer/R/api_opentargets.R    target -> diseases, by Ensembl id
#   genescout/R/tools/opentargets.R         target -> diseases, plus disease -> targets
#   genescout/R/tools/disease_resolver.R    free-text and ontology-id disease lookup
#
# A sibling app queries the same API again, through its own source client and its
# own disease resolver. The disease resolver in particular is the last obvious
# duplication in the family after the clients themselves, which is why it is here
# rather than left in two apps.
#
# Endpoint: https://api.platform.opentargets.org/api/v4/graphql

OPENTARGETS_URL <- "https://api.platform.opentargets.org/api/v4/graphql"
OPENTARGETS_WEB <- "https://platform.opentargets.org"

# Ontology prefixes Open Targets understands. A term matching this is looked up
# directly rather than searched, which is both faster and exact.
OPENTARGETS_ID_PATTERN <- "^(EFO|MONDO|Orphanet|ORPHA|HP|DOID|NCIT|GO|OTAR)[:_]"

OPENTARGETS_DISEASES_QUERY <- paste(
  "query($id: String!, $size: Int!) {",
  "  target(ensemblId: $id) {",
  "    approvedSymbol",
  "    associatedDiseases(page: {index: 0, size: $size}) {",
  "      count",
  "      rows { score disease { id name } }",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

OPENTARGETS_TARGETS_QUERY <- paste(
  "query($id: String!, $size: Int!) {",
  "  disease(efoId: $id) {",
  "    name",
  "    associatedTargets(page: {index: 0, size: $size}) {",
  "      count",
  "      rows { score target { id approvedSymbol } }",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

OPENTARGETS_SEARCH_QUERY <- paste(
  "query($q: String!, $size: Int!) {",
  "  search(queryString: $q, entityNames: [\"disease\"],",
  "         page: {index: 0, size: $size}) {",
  "    hits { id name score description }",
  "  }",
  "}",
  sep = "\n"
)

OPENTARGETS_LOOKUP_QUERY <- paste(
  "query($efoId: String!) {",
  "  disease(efoId: $efoId) { id name description }",
  "}",
  sep = "\n"
)

# Known drugs and clinical candidates. Each row is one drug with its highest
# clinical stage and every disease it has been tried against.
OPENTARGETS_DRUGS_QUERY <- paste(
  "query($id: String!) {",
  "  target(ensemblId: $id) {",
  "    drugAndClinicalCandidates {",
  "      count",
  "      rows {",
  "        maxClinicalStage",
  "        drug { id name drugType }",
  "        diseases { diseaseFromSource disease { id name } }",
  "      }",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

# Pharmacogenomics: variant or genotype to drug-response effect, with an
# evidence level. Most genes carry none; a pharmacogene such as CYP2C19 carries
# many, so an empty result here is the normal case rather than a failure.
OPENTARGETS_PGX_QUERY <- paste(
  "query($id: String!) {",
  "  target(ensemblId: $id) {",
  "    pharmacogenomics {",
  "      variantRsId",
  "      genotypeId",
  "      drugs { drugFromSource }",
  "      phenotypeText",
  "      genotypeAnnotationText",
  "      evidenceLevel",
  "    }",
  "  }",
  "}",
  sep = "\n"
)

#' Is a term an ontology id rather than free text
#'
#' Open Targets takes an exact id for a direct lookup and free text for a search,
#' and the two are different queries. This is what tells them apart.
#'
#' @param term A disease term.
#'
#' @return A single logical.
#'
#' @examples
#' opentargets_is_id("MONDO:0018975")
#' opentargets_is_id("neurofibromatosis type 1")
#'
#' @export
opentargets_is_id <- function(term) {
  grepl(OPENTARGETS_ID_PATTERN, trimws(term %||% ""), ignore.case = TRUE)
}

# Open Targets ids are underscore-separated. Only the ":" separator is swapped,
# because some ontologies keep mixed case in the body of the id (Orphanet_636),
# and upcasing the whole thing would stop it resolving.
opentargets_normalize_id <- function(id) {
  gsub(":", "_", trimws(as.character(id)))
}

#' Turn target-to-disease rows into a table
#'
#' Pure. Rows arrive already sorted by score, best first, and that order is kept.
#'
#' @param body A parsed Open Targets GraphQL response body.
#' @param ensembl_id The Ensembl gene id that was queried, used to build
#'   `source_url`.
#'
#' @return A tibble of `disease`, `disease_id`, `score`, and `source_url`, or
#'   `NULL` when the target is absent or has no associations.
#'
#' @examples
#' body <- list(data = list(target = list(
#'   approvedSymbol = "TP53",
#'   associatedDiseases = list(count = 1, rows = list(
#'     list(score = 0.88, disease = list(id = "MONDO_0018875", name = "LFS"))
#'   ))
#' )))
#' opentargets_parse_diseases(body, "ENSG00000141510")
#'
#' @export
opentargets_parse_diseases <- function(body, ensembl_id = NA_character_) {
  rows <- biohttp::pluck_at(
    body,
    "data",
    "target",
    "associatedDiseases",
    "rows"
  )
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  disease_id <- col_chr(rows, "disease", "id")
  tibble::tibble(
    disease = col_chr(rows, "disease", "name"),
    disease_id = disease_id,
    score = col_num(rows, "score"),
    source_url = paste0(
      OPENTARGETS_WEB,
      "/evidence/",
      ensembl_id,
      "/",
      disease_id
    )
  )
}

#' Turn disease-to-target rows into a table
#'
#' Pure. The discovery direction: the genes Open Targets associates with a
#' disease, each with its overall association score.
#'
#' @param body A parsed Open Targets GraphQL response body.
#' @param disease_id The disease id that was queried, used to build `source_url`.
#'
#' @return A tibble of `symbol`, `ensembl_id`, `score`, and `source_url`, or
#'   `NULL` when the disease is absent or has no associations.
#'
#' @examples
#' body <- list(data = list(disease = list(
#'   name = "NF1",
#'   associatedTargets = list(count = 1, rows = list(
#'     list(score = 0.9, target = list(id = "ENSG00000196712", approvedSymbol = "NF1"))
#'   ))
#' )))
#' opentargets_parse_targets(body, "MONDO_0018975")
#'
#' @export
opentargets_parse_targets <- function(body, disease_id = NA_character_) {
  rows <- biohttp::pluck_at(
    body,
    "data",
    "disease",
    "associatedTargets",
    "rows"
  )
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  ensembl_id <- col_chr(rows, "target", "id")
  tibble::tibble(
    symbol = col_chr(rows, "target", "approvedSymbol"),
    ensembl_id = ensembl_id,
    score = col_num(rows, "score"),
    source_url = paste0(
      OPENTARGETS_WEB,
      "/evidence/",
      ensembl_id,
      "/",
      disease_id
    )
  )
}

#' Turn a disease search or lookup into a table
#'
#' Pure, and it reads **both** response shapes on purpose: the `search` query
#' returns a `hits` array, the direct `disease` lookup returns a single record.
#' A caller resolving a term does not know in advance which query ran, so one
#' parser reads either and the difference stops mattering downstream.
#'
#' @param body A parsed Open Targets GraphQL response body.
#'
#' @return A tibble of `id`, `name`, `score`, `description`, and `source_url`,
#'   best match first, or `NULL` when nothing matched. `score` is `NA` for a
#'   direct lookup, which has no relevance score to report.
#'
#' @examples
#' hits <- list(data = list(search = list(hits = list(
#'   list(id = "MONDO_0018975", name = "neurofibromatosis type 1", score = 24.9)
#' ))))
#' opentargets_parse_matches(hits)
#'
#' single <- list(data = list(disease = list(id = "EFO_0000508", name = "neurofibroma")))
#' opentargets_parse_matches(single)
#'
#' @export
opentargets_parse_matches <- function(body) {
  hits <- biohttp::pluck_at(body, "data", "search", "hits")
  if (is.null(hits) || length(hits) == 0) {
    # The single-disease shape, from a direct id lookup.
    disease <- biohttp::pluck_at(body, "data", "disease")
    if (is.null(disease)) {
      return(NULL)
    }
    hits <- list(disease)
  }
  ids <- col_chr(hits, "id")
  tibble::tibble(
    id = ids,
    name = col_chr(hits, "name"),
    score = col_num(hits, "score"),
    description = col_chr(hits, "description"),
    source_url = paste0(OPENTARGETS_WEB, "/disease/", ids)
  )
}

#' Turn known-drug rows into a table
#'
#' Pure.
#'
#' @section Every disease, not just the first:
#' A drug row carries the whole list of diseases it has been tried against, and
#' the first entry is often the least useful one. In the stored BRAF response,
#' BELVARAFENIB lists five, and the first has no mapped `disease` at all. So the
#' diseases arrive as list columns rather than as one picked value.
#'
#' `diseases` prefers the mapped disease name and falls back to
#' `diseaseFromSource`, which is the label the trial registry used. `disease_ids`
#' is `NA` in the positions Open Targets could not map.
#'
#' `max_phase` is the raw `maxClinicalStage`, for example `"PHASE_2"`. Turning
#' that into "Phase 2" is presentation and belongs to whatever is presenting it.
#'
#' @param body A parsed Open Targets GraphQL response body.
#'
#' @return A tibble of `drug`, `drug_id`, `drug_type`, `max_phase`, `diseases`,
#'   `disease_ids`, and `source_url`. `NULL` when the target has no known drugs.
#'
#' @examples
#' body <- list(data = list(target = list(
#'   drugAndClinicalCandidates = list(count = 1, rows = list(
#'     list(
#'       maxClinicalStage = "PHASE_2",
#'       drug = list(id = "CHEMBL1", name = "DRUGX", drugType = "Small molecule"),
#'       diseases = list(list(
#'         diseaseFromSource = "melanoma",
#'         disease = list(id = "MONDO_0005105", name = "melanoma")
#'       ))
#'     )
#'   ))
#' )))
#' opentargets_parse_drugs(body)
#'
#' @export
opentargets_parse_drugs <- function(body) {
  rows <- biohttp::pluck_at(
    body,
    "data",
    "target",
    "drugAndClinicalCandidates",
    "rows"
  )
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  drug_id <- col_chr(rows, "drug", "id")
  diseases <- lapply(rows, function(row) {
    entries <- biohttp::pluck_at(row, "diseases", default = list())
    names <- col_chr(entries, "disease", "name")
    from_source <- col_chr(entries, "diseaseFromSource")
    names[is.na(names)] <- from_source[is.na(names)]
    unname(names)
  })
  tibble::tibble(
    drug = col_chr(rows, "drug", "name"),
    drug_id = drug_id,
    drug_type = col_chr(rows, "drug", "drugType"),
    max_phase = col_chr(rows, "maxClinicalStage"),
    diseases = diseases,
    disease_ids = lapply(rows, function(row) {
      entries <- biohttp::pluck_at(row, "diseases", default = list())
      unname(col_chr(entries, "disease", "id"))
    }),
    source_url = paste0(OPENTARGETS_WEB, "/drug/", drug_id)
  )
}

#' Turn pharmacogenomics rows into a table
#'
#' Pure.
#'
#' `drugs` is a list column because one annotation can name several. `rsid` is
#' `NA` where the annotation is keyed on a genotype rather than a variant, which
#' is common.
#'
#' @param body A parsed Open Targets GraphQL response body.
#'
#' @return A tibble of `rsid`, `genotype_id`, `drugs`, `phenotype`,
#'   `annotation`, and `evidence_level`. `NULL` when the target has none, which
#'   is the normal case for most genes.
#'
#' @examples
#' body <- list(data = list(target = list(pharmacogenomics = list(
#'   list(
#'     variantRsId = "rs4244285",
#'     drugs = list(list(drugFromSource = "venlafaxine")),
#'     phenotypeText = "decreased metabolism of venlafaxine",
#'     evidenceLevel = "3"
#'   )
#' ))))
#' opentargets_parse_pgx(body)
#'
#' @export
opentargets_parse_pgx <- function(body) {
  rows <- biohttp::pluck_at(body, "data", "target", "pharmacogenomics")
  if (is.null(rows) || length(rows) == 0) {
    return(NULL)
  }
  tibble::tibble(
    rsid = col_chr(rows, "variantRsId"),
    genotype_id = col_chr(rows, "genotypeId"),
    drugs = lapply(rows, function(row) {
      entries <- biohttp::pluck_at(row, "drugs", default = list())
      names <- unname(col_chr(entries, "drugFromSource"))
      unique(names[!is.na(names) & nzchar(names)])
    }),
    phenotype = col_chr(rows, "phenotypeText"),
    annotation = col_chr(rows, "genotypeAnnotationText"),
    evidence_level = col_chr(rows, "evidenceLevel")
  )
}

# Every entry point here is the same POST with a different query and variables,
# so the shared half lives once.
opentargets_post <- function(query, variables, ...) {
  res <- biohttp::post_json(
    OPENTARGETS_URL,
    body = list(query = query, variables = variables),
    source = "Open Targets",
    ...
  )
  biohttp::graphql_error(res, "Open Targets") %||% res
}

#' Diseases associated with a gene
#'
#' Takes an **Ensembl gene id**, not a symbol, because that is what the Open
#' Targets `target` query accepts. Resolve a symbol first with [mygene_gene()]
#' and read its `ensembl_gene` column.
#'
#' Keeping the resolution out of this function is deliberate. A client that
#' silently made a second call to a different service would make one failure look
#' like the other's, and a MyGene outage would read as an Open Targets outage.
#'
#' @param ensembl_id An Ensembl gene id, for example `"ENSG00000141510"`.
#' @param size How many associations to ask for.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [opentargets_parse_diseases()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(opentargets_gene_diseases("ENSG00000141510"))
#' }
#'
#' @export
opentargets_gene_diseases <- function(ensembl_id, size = 20, ...) {
  if (biohttp::is_blank(ensembl_id)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      detail = "no Ensembl gene id was supplied"
    ))
  }
  id <- trimws(as.character(ensembl_id))
  res <- opentargets_post(
    OPENTARGETS_DISEASES_QUERY,
    list(id = id, size = as.integer(size)),
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- opentargets_parse_diseases(res$data, id)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      http = res$http,
      detail = paste0("no Open Targets associations for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Open Targets", http = res$http)
}

#' Genes associated with a disease
#'
#' The discovery direction. `disease_id` is an Open Targets id in either
#' separator form; `MONDO:0018975` and `MONDO_0018975` both work.
#'
#' @param disease_id An EFO, MONDO, HP, or Orphanet id.
#' @param size How many associations to ask for.
#' @inheritParams opentargets_gene_diseases
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [opentargets_parse_targets()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(opentargets_disease_targets("MONDO_0018975"))
#' }
#'
#' @export
opentargets_disease_targets <- function(disease_id, size = 1000, ...) {
  if (biohttp::is_blank(disease_id)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      detail = "no disease id was supplied"
    ))
  }
  id <- opentargets_normalize_id(disease_id)
  res <- opentargets_post(
    OPENTARGETS_TARGETS_QUERY,
    list(id = id, size = as.integer(size)),
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- opentargets_parse_targets(res$data, id)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      http = res$http,
      detail = paste0("no Open Targets targets for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Open Targets", http = res$http)
}

#' Resolve a disease term to ontology records
#'
#' Free text is searched; something that already looks like an ontology id is
#' looked up directly, which is both exact and cheaper. [opentargets_is_id()] is
#' what decides.
#'
#' This is the `resolve_disease` that exists twice in the family, once in
#' `genescout/R/tools/disease_resolver.R` and once in a sibling app.
#'
#' @param term Free text, or an EFO/MONDO/HP/Orphanet id.
#' @param limit How many candidates to return for a free-text search.
#' @inheritParams opentargets_gene_diseases
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [opentargets_parse_matches()], best match first.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(opentargets_resolve_disease("neurofibromatosis type 1"))
#' biohttp::body_or_null(opentargets_resolve_disease("MONDO:0018975"))
#' }
#'
#' @export
opentargets_resolve_disease <- function(term, limit = 5, ...) {
  term <- trimws(as.character(term %||% ""))
  if (biohttp::is_blank(term)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      detail = "no disease term was supplied"
    ))
  }
  if (opentargets_is_id(term)) {
    res <- opentargets_post(
      OPENTARGETS_LOOKUP_QUERY,
      list(efoId = opentargets_normalize_id(term)),
      ...
    )
  } else {
    res <- opentargets_post(
      OPENTARGETS_SEARCH_QUERY,
      list(q = term, size = as.integer(limit)),
      ...
    )
  }
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- opentargets_parse_matches(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      http = res$http,
      detail = paste0("no Open Targets disease matched ", term)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Open Targets", http = res$http)
}

#' Known drugs and clinical candidates for a gene
#'
#' Takes an Ensembl gene id, the same as [opentargets_gene_diseases()].
#'
#' @inheritParams opentargets_gene_diseases
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [opentargets_parse_drugs()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(opentargets_drugs("ENSG00000157764"))
#' }
#'
#' @export
opentargets_drugs <- function(ensembl_id, ...) {
  if (biohttp::is_blank(ensembl_id)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      detail = "no Ensembl gene id was supplied"
    ))
  }
  id <- trimws(as.character(ensembl_id))
  res <- opentargets_post(OPENTARGETS_DRUGS_QUERY, list(id = id), ...)
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- opentargets_parse_drugs(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      http = res$http,
      detail = paste0("no Open Targets known drugs for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Open Targets", http = res$http)
}

#' Pharmacogenomics annotations for a gene
#'
#' Takes an Ensembl gene id, the same as [opentargets_gene_diseases()].
#'
#' Most genes have none, so `no_data` here is the ordinary answer rather than a
#' sign anything went wrong.
#'
#' @inheritParams opentargets_gene_diseases
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [opentargets_parse_pgx()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(opentargets_pgx("ENSG00000165841"))
#' }
#'
#' @export
opentargets_pgx <- function(ensembl_id, ...) {
  if (biohttp::is_blank(ensembl_id)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      detail = "no Ensembl gene id was supplied"
    ))
  }
  id <- trimws(as.character(ensembl_id))
  res <- opentargets_post(OPENTARGETS_PGX_QUERY, list(id = id), ...)
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- opentargets_parse_pgx(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Open Targets",
      http = res$http,
      detail = paste0("no Open Targets pharmacogenomics for ", id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Open Targets", http = res$http)
}
