# UniProt: curated disease involvement, and protein sequence features.
#
# Ported from genescout/R/tools/uniprot.R (curated diseases) and the two copies
# of the features client in variant-reviewer/R/api_proteins.R and
# a sibling app's api_proteins.R.
#
# TWO HOSTS, ON PURPOSE.
#
# The curated-disease data comes from rest.uniprot.org and the sequence features
# from the EBI Proteins API at www.ebi.ac.uk. They are different services with
# different availability, so they carry different `source` labels and land on
# different circuit breakers. One being down must not take the other out.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# genescout also matches disease names against a review's context, with a
# stop-word list so two unrelated diseases do not match on a shared word like
# "syndrome". That is relevance judgement for a particular review, and it stays
# with the app doing the reviewing.

UNIPROT_URL <- "https://rest.uniprot.org/uniprotkb"
UNIPROT_DISEASE_WEB <- "https://www.uniprot.org/diseases"
PROTEINS_URL <- "https://www.ebi.ac.uk/proteins/api"

# The feature types worth asking for when placing a variant in a protein's
# architecture. UniProt has many more, most of which are noise for that question.
UNIPROT_FEATURE_TYPES <- c(
  "DOMAIN",
  "REGION",
  "MOTIF",
  "REPEAT",
  "ZN_FING",
  "DNA_BIND",
  "CA_BIND",
  "NP_BIND",
  "ACT_SITE",
  "BINDING",
  "SITE"
)

# Readable labels for UniProt's terse type codes.
UNIPROT_FEATURE_LABELS <- c(
  DOMAIN = "Domain",
  REGION = "Region",
  MOTIF = "Motif",
  REPEAT = "Repeat",
  ZN_FING = "Zinc finger",
  DNA_BIND = "DNA binding",
  CA_BIND = "Calcium binding",
  NP_BIND = "Nucleotide binding",
  ACT_SITE = "Active site",
  BINDING = "Binding site",
  SITE = "Site"
)

#' Turn a UniProtKB entry into a curated-disease table
#'
#' Pure. Reads the `DISEASE` comments of an entry.
#'
#' @section The causal distinction:
#' UniProt's disease notes separate two very different claims. "The disease is
#' caused by variants affecting the gene" is a Mendelian statement; "may be
#' involved in the pathogenesis" is far weaker. Both arrive as `DISEASE`
#' comments, so treating them alike would promote a speculative association to a
#' causal one. `causal` carries the distinction, read from the note text.
#'
#' @param body A parsed UniProtKB entry.
#'
#' @return A tibble of `id`, `name`, `acronym`, `mim`, `causal`, and
#'   `source_url`. `NULL` when the entry curates no disease involvement.
#'
#' @examples
#' body <- list(comments = list(list(
#'   commentType = "DISEASE",
#'   disease = list(
#'     diseaseAccession = "DI-00218", diseaseId = "Li-Fraumeni syndrome",
#'     acronym = "LFS", diseaseCrossReference = list(database = "MIM", id = "151623")
#'   ),
#'   note = list(texts = list(list(value = "The disease is caused by variants")))
#' )))
#' uniprot_parse_diseases(body)
#'
#' @export
uniprot_parse_diseases <- function(body) {
  comments <- biohttp::pluck_at(body, "comments")
  if (is.null(comments) || length(comments) == 0) {
    return(NULL)
  }
  rows <- list()
  for (comment in comments) {
    if (!identical(biohttp::pluck_at(comment, "commentType"), "DISEASE")) {
      next
    }
    disease <- biohttp::pluck_at(comment, "disease")
    id <- as.character(
      biohttp::pluck_at(disease, "diseaseAccession", default = "")
    )
    name <- as.character(biohttp::pluck_at(disease, "diseaseId", default = ""))
    # A comment naming no disease carries nothing a caller can ground on.
    if (!nzchar(id) || !nzchar(name)) {
      next
    }
    xref <- biohttp::pluck_at(disease, "diseaseCrossReference")
    mim <- if (identical(biohttp::pluck_at(xref, "database"), "MIM")) {
      chr_at(xref, "id")
    } else {
      NA_character_
    }
    note <- tolower(paste(
      unlist(
        biohttp::pluck_at(comment, "note", "texts") %||% list(),
        use.names = FALSE
      ),
      collapse = " "
    ))
    rows[[length(rows) + 1L]] <- list(
      id = id,
      name = name,
      acronym = as.character(
        biohttp::pluck_at(disease, "acronym", default = "")
      ),
      mim = mim,
      causal = grepl("caused by", note, fixed = TRUE)
    )
  }
  if (length(rows) == 0) {
    return(NULL)
  }
  ids <- vapply(rows, function(row) row$id, character(1))
  tibble::tibble(
    id = ids,
    name = vapply(rows, function(row) row$name, character(1)),
    acronym = vapply(rows, function(row) row$acronym, character(1)),
    mim = vapply(rows, function(row) row$mim, character(1)),
    causal = vapply(rows, function(row) isTRUE(row$causal), logical(1)),
    source_url = paste0(UNIPROT_DISEASE_WEB, "/", ids)
  )
}

#' Turn an EBI Proteins features response into a table
#'
#' Pure. Sorted by start position.
#'
#' @param body A parsed EBI Proteins `features` response.
#'
#' @return A tibble of `type`, `label`, `description`, `begin`, and `end`. A
#'   zero-row tibble when the entry has no features of the requested types,
#'   because "this protein has no annotated domains" is a real answer.
#'
#' @examples
#' body <- list(features = list(
#'   list(type = "DOMAIN", description = "Protein kinase", begin = "457", end = "717")
#' ))
#' uniprot_parse_features(body)
#'
#' @export
uniprot_parse_features <- function(body) {
  empty <- tibble::tibble(
    type = character(),
    label = character(),
    description = character(),
    begin = integer(),
    end = integer()
  )
  features <- biohttp::pluck_at(body, "features")
  if (is.null(features) || length(features) == 0) {
    return(empty)
  }
  type <- col_chr(features, "type")
  labelled <- unname(UNIPROT_FEATURE_LABELS[type])
  out <- tibble::tibble(
    type = type,
    # An unmapped code falls back to itself rather than becoming NA, so a new
    # UniProt feature type still reads sensibly.
    label = ifelse(is.na(labelled), type, labelled),
    description = col_chr(features, "description"),
    # begin and end arrive as strings.
    begin = suppressWarnings(as.integer(col_chr(features, "begin"))),
    end = suppressWarnings(as.integer(col_chr(features, "end")))
  )
  out[order(out$begin), , drop = FALSE]
}

#' The features spanning a residue position
#'
#' Pure. Which annotated domains, sites, or regions contain a given residue,
#' which is the question a variant reviewer actually asks.
#'
#' @param features A feature tibble from [uniprot_parse_features()].
#' @param position A residue position.
#'
#' @return The rows of `features` whose `begin`/`end` span `position`.
#'
#' @examples
#' features <- uniprot_parse_features(list(features = list(
#'   list(type = "DOMAIN", description = "Kinase", begin = "457", end = "717")
#' )))
#' uniprot_features_at(features, 600)
#' uniprot_features_at(features, 100)
#'
#' @export
uniprot_features_at <- function(features, position) {
  if (!is.data.frame(features) || nrow(features) == 0) {
    return(features)
  }
  pos <- suppressWarnings(as.integer(position))
  if (length(pos) != 1L || is.na(pos)) {
    return(features[0, , drop = FALSE])
  }
  features[
    !is.na(features$begin) &
      !is.na(features$end) &
      features$begin <= pos &
      features$end >= pos,
    ,
    drop = FALSE
  ]
}

#' Diseases UniProt curates for an accession
#'
#' @param accession A UniProt accession, for example `"P04637"`.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [uniprot_parse_diseases()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(uniprot_diseases("P04637"))
#' }
#'
#' @export
uniprot_diseases <- function(accession, ...) {
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "UniProt",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  res <- biohttp::get_json(
    UNIPROT_URL,
    path = paste0(cleaned, ".json"),
    query = list(fields = "cc_disease"),
    source = "UniProt",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- uniprot_parse_diseases(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "UniProt",
      http = res$http,
      detail = paste0("UniProt curates no disease involvement for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "UniProt", http = res$http)
}

#' Sequence features for an accession
#'
#' From the EBI Proteins API, which is a different host to [uniprot_diseases()].
#'
#' @param accession A UniProt accession.
#' @param types Feature types to request. Defaults to the set worth showing when
#'   placing a variant in a protein's architecture.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [uniprot_parse_features()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(uniprot_features("P15056"))
#' }
#'
#' @export
uniprot_features <- function(
  accession,
  types = UNIPROT_FEATURE_TYPES,
  ...
) {
  cleaned <- clean_accession(accession)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "EBI Proteins",
      detail = "no usable UniProt accession was supplied"
    ))
  }
  res <- biohttp::get_json(
    PROTEINS_URL,
    path = paste0("features/", cleaned),
    query = list(types = paste(types, collapse = ",")),
    source = "EBI Proteins",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  biohttp::status_ok(
    data = uniprot_parse_features(res$data),
    source = "EBI Proteins",
    http = res$http
  )
}
