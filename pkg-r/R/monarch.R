# Monarch Initiative: entity search and the associations between entities.
#
# Ported from a sibling app's monarch.R (search and association) and
# variant-reviewer/R/api_monarch.R (gene to phenotype).
#
# A SEARCH RESULT DOES NOT DISAMBIGUATE ITSELF.
#
# Searching a gene symbol returns the human gene and then the orthologs, under
# the identical symbol. `FBN1` returns HGNC:3603 in Homo sapiens, then
# NCBIGene:373992 in Gallus gallus and NCBIGene:478293 in Canis lupus familiaris,
# all three named "FBN1". Taking the first result happens to be right here and is
# not right in general, so `taxon` is carried on every row and a caller that
# wants human results has to say so.
#
# THE CURIE COLON GOES THROUGH VERBATIM.
#
# Both segments of the entity path are CURIEs, `HGNC:11998` and
# `biolink:GeneToPhenotypicFeatureAssociation`, and Monarch wants the colons as
# they are rather than percent-encoded. biohttp's path builder leaves them
# alone, which is what makes this work, so a test pins it.
#
# ONE HOST, AFTER A LIVE CHECK.
#
# This used to be two. The apps this was ported from called search and
# association on `api-v3.monarchinitiative.org` and the entity route on
# `api.monarchinitiative.org`, and the ported client kept them apart rather than
# assume they were interchangeable. The cost was that biohttp saw two hosts and
# gave Monarch two circuit breakers, so one host going down only opened half of
# them.
#
# A live check settled it. All three routes answer on both hosts, and the entity
# payload is identical between them, the same total and the same ids in the same
# order. So they collapse to one, and `test-live.R` pins that.
#
# WHAT IS DELIBERATELY NOT HERE.
#
# One consuming app turns these associations into a graph model: it drops edges
# with no primary_knowledge_source, merges duplicates, and prunes nodes left
# stranded. That is the app's model, not Monarch's response, so it stays there.
# See the note on monarch_associations() about what that means for a caller
# fetching both directions.
#
# Endpoint: https://api-v3.monarchinitiative.org/v3/api

MONARCH_URL <- "https://api-v3.monarchinitiative.org/v3/api"
MONARCH_WEB <- "https://monarchinitiative.org"

# The association category for gene to phenotype, which the entity route takes
# as a path segment rather than as a parameter.
MONARCH_GENE_PHENOTYPE <- "biolink:GeneToPhenotypicFeatureAssociation"

#' Turn a Monarch search response into a table
#'
#' Pure.
#'
#' @section Read the taxon:
#' A gene symbol is not unique across species in Monarch, and a search returns
#' the orthologs alongside the human gene under the same name. `taxon` is the
#' only column that tells them apart. It is `NA` for a disease or a phenotype,
#' which are not species-scoped.
#'
#' @param body A parsed Monarch `search` response.
#'
#' @return A tibble of `id`, `name`, `category`, `description`, `taxon`, and
#'   `source_url`, one row per match. `NULL` when there are none.
#'
#' @examples
#' body <- list(items = list(list(
#'   id = "HGNC:3603",
#'   name = "FBN1",
#'   category = "biolink:Gene",
#'   full_name = "fibrillin 1",
#'   in_taxon_label = "Homo sapiens"
#' )))
#' monarch_parse_search(body)
#'
#' @export
monarch_parse_search <- function(body) {
  items <- biohttp::pluck_at(body, "items")
  if (is.null(items) || length(items) == 0) {
    return(NULL)
  }
  id <- col_chr(items, "id")
  # A gene carries full_name as its readable gloss and a disease carries
  # description. Whichever is present says more than neither.
  description <- col_chr(items, "description")
  full_name <- col_chr(items, "full_name")
  description[is.na(description)] <- full_name[is.na(description)]
  tibble::tibble(
    id = id,
    name = col_chr(items, "name"),
    category = col_chr(items, "category"),
    description = description,
    taxon = col_chr(items, "in_taxon_label"),
    source_url = paste0(MONARCH_WEB, "/", id)
  )
}

#' Turn Monarch association records into a table
#'
#' Pure.
#'
#' One row per association, exactly as Monarch sent it. Nothing is dropped and
#' nothing is merged. `publications` is a list column because an association
#' carries any number of them, and `null` and an array both occur.
#'
#' @param body A parsed Monarch association response, or the `items` list from
#'   one.
#'
#' @return A tibble of `subject`, `subject_label`, `subject_category`,
#'   `predicate`, `object`, `object_label`, `object_category`,
#'   `primary_knowledge_source`, `knowledge_level`, and `publications`. `NULL`
#'   when there are no associations.
#'
#' @examples
#' body <- list(items = list(list(
#'   subject = "MONDO:0017309",
#'   subject_label = "neonatal Marfan syndrome",
#'   predicate = "biolink:has_phenotype",
#'   object = "HP:0001653",
#'   object_label = "Mitral regurgitation",
#'   publications = NULL
#' )))
#' monarch_parse_associations(body)
#'
#' @export
monarch_parse_associations <- function(body) {
  items <- biohttp::pluck_at(body, "items") %||% body
  if (is.null(items) || length(items) == 0) {
    return(NULL)
  }
  tibble::tibble(
    subject = col_chr(items, "subject"),
    subject_label = col_chr(items, "subject_label"),
    subject_category = col_chr(items, "subject_category"),
    predicate = col_chr(items, "predicate"),
    object = col_chr(items, "object"),
    object_label = col_chr(items, "object_label"),
    object_category = col_chr(items, "object_category"),
    primary_knowledge_source = col_chr(items, "primary_knowledge_source"),
    knowledge_level = col_chr(items, "knowledge_level"),
    publications = lapply(items, function(item) {
      as.character(unlist(
        biohttp::pluck_at(item, "publications", default = list()),
        use.names = FALSE
      ))
    })
  )
}

#' Normalise an HGNC id to the CURIE form Monarch expects
#'
#' Monarch's entity route takes `HGNC:11998`. MyGene returns the bare digits.
#' Exported because a caller assembling its own entity id needs the same rule.
#'
#' @param hgnc An HGNC id, bare (`"11998"`) or prefixed (`"HGNC:11998"`).
#'
#' @return A single string, or `NULL` when there is nothing usable.
#'
#' @examples
#' monarch_hgnc_id("11998")
#' monarch_hgnc_id("hgnc:11998")
#'
#' @export
monarch_hgnc_id <- function(hgnc) {
  id <- toupper(trimws(as.character(hgnc %||% "")))
  id <- sub("^HGNC:", "", id)
  # The id goes into a URL path, so anything that is not one is refused rather
  # than turned into a request for some other entity.
  if (!grepl("^[0-9]+$", id)) {
    return(NULL)
  }
  paste0("HGNC:", id)
}

#' Search Monarch for an entity
#'
#' Turns typed text into the entity id the other calls need.
#'
#' @param text Free text, for example a gene symbol or a disease name.
#' @param limit Maximum matches to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `matches` (the tibble
#'   from [monarch_parse_search()]) and `total`. `total` is Monarch's own count,
#'   not the page size, so a caller can say how much was left behind rather than
#'   presenting the first few as everything.
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(monarch_search("Marfan syndrome"))$matches
#' }
#'
#' @export
monarch_search <- function(text, limit = 10, ...) {
  if (biohttp::is_blank(text)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      detail = "no search text was supplied"
    ))
  }
  res <- biohttp::get_json(
    MONARCH_URL,
    path = "search",
    query = list(q = as.character(text), limit = as.integer(limit)),
    source = "Monarch",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- monarch_parse_search(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      http = res$http,
      detail = paste0("Monarch found nothing for ", text)
    ))
  }
  biohttp::status_ok(
    data = list(
      matches = parsed,
      total = as.integer(
        biohttp::pluck_at(res$data, "total", default = nrow(parsed))
      )
    ),
    source = "Monarch",
    http = res$http
  )
}

#' Associations with an entity on one end
#'
#' @section One direction per call:
#' `end` picks which end of the association the entity sits on. Both directions
#' are two calls, and concatenating them returns any association that has the
#' entity on both ends **twice**. Collapsing those, and deciding what to do with
#' the publications on each copy, is the caller's model rather than Monarch's
#' answer, so this function does neither.
#'
#' @param entity An entity CURIE, for example `"MONDO:0007947"`.
#' @param end `"subject"` or `"object"`.
#' @param limit Maximum associations to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is the tibble described in
#'   [monarch_parse_associations()].
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(monarch_associations("MONDO:0007947"))
#' }
#'
#' @export
monarch_associations <- function(
  entity,
  end = c("subject", "object"),
  limit = 200,
  ...
) {
  end <- match.arg(end)
  if (biohttp::is_blank(entity)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      detail = "no entity id was supplied"
    ))
  }
  query <- list(limit = as.integer(limit))
  query[[end]] <- as.character(entity)
  res <- biohttp::get_json(
    MONARCH_URL,
    path = "association",
    query = query,
    source = "Monarch",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- monarch_parse_associations(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      http = res$http,
      detail = paste0("Monarch has no associations for ", entity)
    ))
  }
  biohttp::status_ok(data = parsed, source = "Monarch", http = res$http)
}

#' HPO phenotypes Monarch associates with a gene
#'
#' Keyed on the HGNC id, which [mygene_gene()] returns as `hgnc`. The rows come
#' back in the shape [monarch_parse_associations()] describes, so the phenotype
#' term is `object` and its label is `object_label`.
#'
#' @param hgnc An HGNC id, bare or prefixed. See [monarch_hgnc_id()].
#' @param limit Maximum associations to return.
#' @param ... Passed to [biohttp::get_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a list of `associations` (the
#'   tibble) and `total`.
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(monarch_gene_phenotypes("11998"))$associations
#' }
#'
#' @export
monarch_gene_phenotypes <- function(hgnc, limit = 50, ...) {
  id <- monarch_hgnc_id(hgnc)
  if (is.null(id)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      detail = "an HGNC id is required for a Monarch gene lookup"
    ))
  }
  res <- biohttp::get_json(
    MONARCH_URL,
    path = paste("entity", id, MONARCH_GENE_PHENOTYPE, sep = "/"),
    query = list(limit = as.integer(limit)),
    source = "Monarch",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  parsed <- monarch_parse_associations(res$data)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "Monarch",
      http = res$http,
      detail = paste0("Monarch has no phenotypes for ", id)
    ))
  }
  biohttp::status_ok(
    data = list(
      associations = parsed,
      total = as.integer(
        biohttp::pluck_at(res$data, "total", default = nrow(parsed))
      )
    ),
    source = "Monarch",
    http = res$http
  )
}
