# Ensembl VEP: consequence prediction.
#
# Ported from a sibling app's api_vep.R.
#
# FOUR TRAPS, ALL ENCODED AS CODE RATHER THAN LEFT AS COMMENTS TO REDISCOVER.
#
# 1. The POST cap is exactly 200. A 201st returns an error body. Verified
#    against the live API, not a tuning guess.
#
# 2. AlphaMissense lives at transcript_consequences[].alphamissense, never at the
#    top level. Reading it from the top level finds nothing, which reads as "VEP
#    does not serve AlphaMissense" rather than as a bug.
#
# 3. The `variants` body must serialise as a JSON ARRAY. A named R list becomes a
#    JSON object and VEP answers 500. unname() is load-bearing.
#
# 4. dbNSFP=<cols> is deliberately never requested. It returns comma-joined
#    multi-transcript strings in dbNSFP's own order, which is NOT aligned to the
#    VEP transcript being reported, so the values silently belong to a different
#    transcript than the rest of the row.
#
# 5. VEP left-trims and renumbers indels. An insertion sent as
#    `1 55516888 . T TA . . .` comes back with start 55516889 and allele_string
#    "-/A", so a key rebuilt from those fields never matches the key the caller
#    built from its own input, and every indel reads as "VEP returned nothing".
#    VEP echoes the exact input line in `input`, and with vcf_string=1 also
#    reports a normalised `vcf_string`, so results are matched on those first
#    and on the rebuilt key only as a fallback for a response that carries
#    neither.
#
# AND: VEP does not promise response order, so results are matched back by
# variant identity rather than by array position.
#
# Endpoint: https://rest.ensembl.org/vep/homo_sapiens/region

VEP_URL <- "https://rest.ensembl.org/vep/homo_sapiens/region"

# The verified hard limit. Recorded in a sibling app's source registry as
# "verified hard limits, not tuning".
VEP_BATCH <- 200L

#' Build a VEP region string from variant components
#'
#' Pure. VEP's region endpoint takes a VCF-like string.
#'
#' @param chrom A chromosome, with or without a `chr` prefix.
#' @param pos A 1-based position.
#' @param ref,alt Reference and alternate alleles.
#'
#' @return A single string.
#'
#' @examples
#' vep_region("7", 140753336, "A", "T")
#'
#' @export
vep_region <- function(chrom, pos, ref, alt) {
  chrom <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  paste(chrom, as.integer(pos), ".", toupper(ref), toupper(alt), ".", ".", ".")
}

#' Build the variant key VEP results are matched on
#'
#' Pure. `vep_variants()` returns a `key` column built this way, so this is what
#' a caller uses to join the result back onto its own variant table.
#'
#' @inheritParams vep_region
#'
#' @return A single string, `chrom-pos-ref-alt`.
#'
#' @examples
#' vep_key("7", 140753336, "A", "T")
#' vep_key("chr7", 140753336, "a", "t")
#'
#' @export
vep_key <- function(chrom, pos, ref, alt) {
  paste(
    sub("^chr", "", as.character(chrom), ignore.case = TRUE),
    as.integer(pos),
    toupper(ref),
    toupper(alt),
    sep = "-"
  )
}

# The key carried by the echoed input line. VEP hands back exactly the
# `chrom pos . ref alt . . .` string it was sent, so this is the one identity
# that survives left-trimming and renumbering of an indel.
vep_input_key <- function(input) {
  if (biohttp::is_blank(input)) {
    return(NA_character_)
  }
  parts <- strsplit(trimws(as.character(input)), "[[:space:]]+")[[1]]
  if (length(parts) < 5) {
    return(NA_character_)
  }
  vep_key(parts[1], parts[2], parts[4], parts[5])
}

# Every key a VEP element can be matched on, most reliable first: the echoed
# input, then the normalised vcf_string, then the key rebuilt from the
# reported coordinates. The last is wrong for an indel (trap 5 in the file
# header) and is kept only for a response that carries neither of the others.
vep_element_keys <- function(element) {
  vcf <- biohttp::pluck_at(element, "vcf_string", default = NA_character_)
  if (!is.na(vcf)) {
    vcf_parts <- strsplit(as.character(vcf), "-", fixed = TRUE)[[1]]
    vcf <- if (length(vcf_parts) == 4) {
      do.call(vep_key, as.list(vcf_parts))
    } else {
      NA_character_
    }
  }
  keys <- c(
    vep_input_key(biohttp::pluck_at(element, "input")),
    vcf,
    vep_element_key(element)
  )
  unique(keys[!is.na(keys)])
}

# The key rebuilt from what a VEP element reports about itself. Right for an
# SNV, wrong for an indel, so it is the last resort in vep_element_keys().
vep_element_key <- function(element) {
  alleles <- strsplit(
    biohttp::pluck_at(element, "allele_string", default = "/"),
    "/",
    fixed = TRUE
  )[[1]]
  vep_key(
    biohttp::pluck_at(element, "seq_region_name", default = NA_character_),
    biohttp::pluck_at(element, "start", default = NA_integer_),
    alleles[1] %||% NA_character_,
    alleles[length(alleles)] %||% NA_character_
  )
}

#' Choose which transcript to report for a variant
#'
#' Pure. A variant hits many transcripts and only one row can be reported.
#'
#' The order is MANE Select first, because that is the community's designated
#' representative transcript; then whichever transcript carries the variant's
#' most severe consequence; then the first VEP listed. Taking the first outright
#' reports an arbitrary transcript, often a non-coding one, which makes the
#' consequence look milder than it is.
#'
#' @param element One parsed VEP result element.
#'
#' @return One transcript consequence, or `NULL` when there are none.
#'
#' @examples
#' element <- list(
#'   most_severe_consequence = "missense_variant",
#'   transcript_consequences = list(
#'     list(consequence_terms = list("intron_variant")),
#'     list(consequence_terms = list("missense_variant"), mane_select = "NM_1")
#'   )
#' )
#' vep_pick_transcript(element)$mane_select
#'
#' @export
vep_pick_transcript <- function(element) {
  consequences <- biohttp::pluck_at(element, "transcript_consequences")
  if (is.null(consequences) || length(consequences) == 0) {
    return(NULL)
  }
  mane <- Filter(
    function(tc) !is.null(tc$mane_select),
    consequences
  )
  if (length(mane) > 0) {
    return(mane[[1]])
  }
  severe <- Filter(
    function(tc) {
      isTRUE(
        biohttp::pluck_at(element, "most_severe_consequence") %in%
          (tc$consequence_terms %||% "")
      )
    },
    consequences
  )
  if (length(severe) > 0) {
    return(severe[[1]])
  }
  consequences[[1]]
}

#' Turn one VEP element into a table row
#'
#' Pure.
#'
#' @section AlphaMissense is nested:
#' It is at `transcript_consequences[].alphamissense$am_pathogenicity`, per
#' transcript. There is nothing at the top level. Hoisting the read out of the
#' transcript is the single easiest way to get `NA` everywhere and conclude the
#' API does not serve it.
#'
#' @param element One parsed VEP result element.
#'
#' @return A one-row tibble.
#'
#' @examples
#' element <- list(
#'   most_severe_consequence = "missense_variant",
#'   transcript_consequences = list(list(
#'     gene_symbol = "BRAF", consequence_terms = list("missense_variant"),
#'     alphamissense = list(am_pathogenicity = 0.99, am_class = "pathogenic")
#'   ))
#' )
#' vep_parse_element(element)
#'
#' @export
vep_parse_element <- function(element) {
  transcript <- vep_pick_transcript(element)
  if (is.null(transcript)) {
    out <- vep_empty_row()
    out$consequence <- chr_at(element, "most_severe_consequence")
    return(out)
  }
  tibble::tibble(
    gene = chr_at(transcript, "gene_symbol"),
    consequence = as.character(
      (transcript$consequence_terms %||%
        list(biohttp::pluck_at(element, "most_severe_consequence")))[[1]]
    ),
    mane = chr_at(transcript, "mane_select"),
    impact = chr_at(transcript, "impact"),
    exon = chr_at(transcript, "exon"),
    # protein_start can be a range for an indel; the start is what positions a
    # residue marker.
    protein_pos = suppressWarnings(
      as.integer(biohttp::pluck_at(transcript, "protein_start", default = NA))
    ),
    sift = chr_at(transcript, "sift_prediction"),
    polyphen = chr_at(transcript, "polyphen_prediction"),
    # Nested, per transcript. Do not hoist this.
    alphamissense = num_at(
      biohttp::pluck_at(transcript, "alphamissense"),
      "am_pathogenicity"
    ),
    alphamissense_class = chr_at(
      biohttp::pluck_at(transcript, "alphamissense"),
      "am_class"
    )
  )
}

vep_empty_row <- function() {
  tibble::tibble(
    gene = NA_character_,
    consequence = NA_character_,
    mane = NA_character_,
    impact = NA_character_,
    exon = NA_character_,
    protein_pos = NA_integer_,
    sift = NA_character_,
    polyphen = NA_character_,
    alphamissense = NA_real_,
    alphamissense_class = NA_character_
  )
}

#' Turn a VEP batch response into a table
#'
#' Pure. One row per requested variant, in the order asked.
#'
#' @section Matched by identity, not by position:
#' VEP does not promise that results come back in the order they were sent.
#' Zipping the response onto the input by index therefore assigns consequences
#' to the wrong variants, silently. Elements are matched on the key built from
#' the echoed `input` line, which is the exact string that was sent.
#'
#' @section Indels are renumbered:
#' VEP left-trims and renumbers an indel, so the `start` and `allele_string` it
#' reports do not rebuild the key the caller asked with. An insertion sent as
#' `1 55516888 . T TA . . .` comes back as start 55516889 and `-/A`, and a
#' delins can come back re-anchored in `vcf_string` too. The echoed `input` is
#' the identity that survives, so it is matched first, then `vcf_string`, and
#' the rebuilt key only for a response carrying neither.
#'
#' @param body A parsed VEP response, an array of elements.
#' @param keys Variant keys from [vep_key()], in the order asked.
#'
#' @return A tibble with one row per entry in `keys`, plus a `key` column.
#'
#' @export
vep_parse_batch <- function(body, keys) {
  by_key <- list()
  for (element in body) {
    row <- vep_parse_element(element)
    for (key in vep_element_keys(element)) {
      # First identity wins. The rebuilt key of an indel is a different string
      # from any real input, so it never shadows a correct match.
      if (is.null(by_key[[key]])) {
        by_key[[key]] <- row
      }
    }
  }
  rows <- lapply(keys, function(key) {
    row <- by_key[[key]] %||% vep_empty_row()
    cbind(tibble::tibble(key = as.character(key)), row)
  })
  tibble::as_tibble(do.call(rbind, rows))
}

#' Consequence predictions for many variants
#'
#' @param chrom,pos,ref,alt Variant components, all the same length.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per variant,
#'   in the order asked. See [vep_parse_batch()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(vep_variants("7", 140753336, "A", "T"))
#' }
#'
#' @export
vep_variants <- function(chrom, pos, ref, alt, ...) {
  n <- length(chrom)
  if (n == 0) {
    return(biohttp::status_no_data(
      source = "VEP",
      detail = "no variants were supplied"
    ))
  }
  if (n > VEP_BATCH) {
    stop(
      "VEP accepts at most ",
      VEP_BATCH,
      " variants per request, not ",
      n,
      call. = FALSE
    )
  }
  regions <- vep_region(chrom, pos, ref, alt)
  keys <- vep_key(chrom, pos, ref, alt)
  res <- biohttp::post_json(
    # vcf_string=1 makes VEP report the normalised variant alongside the echoed
    # input, which is the second of the identities results are matched on.
    paste0(VEP_URL, "?AlphaMissense=1&mane=1&numbers=1&vcf_string=1"),
    # unname() is load-bearing: a named list serialises as a JSON object and VEP
    # answers 500. The value must be a JSON array.
    body = list(variants = unname(as.list(regions))),
    source = "VEP",
    ...
  )
  if (!isTRUE(res$ok)) {
    return(res)
  }
  biohttp::status_ok(
    data = vep_parse_batch(res$data, keys),
    source = "VEP",
    http = res$http
  )
}
