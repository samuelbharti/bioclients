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

#' The request flags VEP is asked for by default
#'
#' Pure. The default for the `options` argument of [vep_variants()] and
#' [vep_variants_all()]. Start from this and add to it rather than replacing
#' it, because `vcf_string` is one of the identities results are matched on
#' and `mane` is what [vep_pick_transcript()] chooses by.
#'
#' @section Supported flags:
#' Every entry is a VEP request parameter, sent as `name=value` on the query
#' string. `TRUE` and `1` send `name=1`; `FALSE`, `0` and `NULL` omit the flag.
#' Names are case sensitive, exactly as VEP spells them. The flags this package
#' parses are:
#'
#' * `AlphaMissense`, `mane`, `numbers`, `vcf_string`: the defaults.
#' * `af`, `af_gnomade`, `af_gnomadg`: colocated variant frequencies, read by
#'   [vep_parse_colocated()].
#' * `CADD`, `SpliceAI`, `REVEL`: per-transcript predictor scores, read by
#'   [vep_parse_element()].
#' * `hgvs`: `hgvsc` and `hgvsp` notation per transcript.
#' * `canonical`: marks the canonical transcript.
#' * `pick`, `pick_allele_gene`: ask VEP to return one transcript per variant
#'   or per allele and gene, rather than all of them.
#' * `protein`, `domains`, `variant_class`: extra transcript annotation, carried
#'   through untouched in the response for a caller parsing it directly.
#' * `LoF`: LOFTEE, read into the `lof` column.
#'
#' `dbNSFP` is refused. It returns comma-joined multi-transcript strings in
#' dbNSFP's own order, not aligned to the transcript being reported, so the
#' values silently belong to a different transcript than the rest of the row.
#'
#' @return A named list of flags.
#'
#' @examples
#' vep_default_options()
#' c(vep_default_options(), list(CADD = 1, REVEL = 1))
#'
#' @export
vep_default_options <- function() {
  list(AlphaMissense = 1, mane = 1, numbers = 1, vcf_string = 1)
}

# Turn the options list into the query string VEP takes. Pure.
vep_query_string <- function(options) {
  if (is.null(options) || length(options) == 0) {
    return("")
  }
  labels <- names(options)
  if (is.null(labels) || any(is.na(labels) | !nzchar(labels))) {
    stop("VEP options must be a named list of flags", call. = FALSE)
  }
  if ("dbNSFP" %in% labels) {
    stop(
      "dbNSFP is not supported: it returns values in dbNSFP's own transcript ",
      "order, not aligned to the transcript VEP reports. See the note at the ",
      "top of R/vep.R.",
      call. = FALSE
    )
  }
  parts <- character()
  for (label in labels) {
    value <- options[[label]]
    off <- is.null(value) ||
      isFALSE(value) ||
      identical(as.character(value), "0")
    if (off) {
      next
    }
    if (isTRUE(value)) {
      value <- 1
    }
    parts <- c(
      parts,
      paste0(label, "=", utils::URLencode(as.character(value), reserved = TRUE))
    )
  }
  paste(parts, collapse = "&")
}

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
#' @section Optional columns:
#' `hgvsc`, `hgvsp`, `canonical`, `cadd_phred`, `cadd_raw`, `revel`, the
#' `spliceai_*` scores and `lof` are only populated when the matching flag was
#' requested, see [vep_default_options()]. Otherwise they are `NA`. VEP marks
#' only the canonical transcript, so `canonical` is `TRUE` on that transcript
#' and `NA` everywhere else, including when the flag was never asked.
#' `spliceai_max` is the largest of the four SpliceAI delta scores.
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
  # Named `picked` rather than `transcript` because tibble() exposes each
  # column to the arguments after it, and a `transcript` column is defined
  # below. A local of the same name would be shadowed by the column halfway
  # through the call and every later field would read off a string.
  picked <- vep_pick_transcript(element)
  if (is.null(picked)) {
    out <- vep_empty_row()
    out$consequence <- chr_at(element, "most_severe_consequence")
    return(out)
  }
  spliceai <- biohttp::pluck_at(picked, "spliceai")
  spliceai_scores <- c(
    num_at(spliceai, "DS_AG"),
    num_at(spliceai, "DS_AL"),
    num_at(spliceai, "DS_DG"),
    num_at(spliceai, "DS_DL")
  )
  tibble::tibble(
    gene = chr_at(picked, "gene_symbol"),
    consequence = as.character(
      (picked$consequence_terms %||%
        list(biohttp::pluck_at(element, "most_severe_consequence")))[[1]]
    ),
    mane = chr_at(picked, "mane_select"),
    impact = chr_at(picked, "impact"),
    exon = chr_at(picked, "exon"),
    # protein_start can be a range for an indel; the start is what positions a
    # residue marker.
    protein_pos = suppressWarnings(
      as.integer(biohttp::pluck_at(picked, "protein_start", default = NA))
    ),
    sift = chr_at(picked, "sift_prediction"),
    polyphen = chr_at(picked, "polyphen_prediction"),
    # Nested, per transcript. Do not hoist this.
    alphamissense = num_at(
      biohttp::pluck_at(picked, "alphamissense"),
      "am_pathogenicity"
    ),
    alphamissense_class = chr_at(
      biohttp::pluck_at(picked, "alphamissense"),
      "am_class"
    ),
    transcript = chr_at(picked, "transcript_id"),
    gene_id = chr_at(picked, "gene_id"),
    biotype = chr_at(picked, "biotype"),
    hgvsc = chr_at(picked, "hgvsc"),
    hgvsp = chr_at(picked, "hgvsp"),
    # VEP sends canonical: 1 on the canonical transcript and nothing on the
    # others, so absence cannot tell "not canonical" from "not asked".
    canonical = if (identical(num_at(picked, "canonical"), 1)) TRUE else NA,
    codons = chr_at(picked, "codons"),
    amino_acids = chr_at(picked, "amino_acids"),
    cadd_phred = num_at(picked, "cadd_phred"),
    cadd_raw = num_at(picked, "cadd_raw"),
    revel = num_at(picked, "revel"),
    spliceai_ds_ag = spliceai_scores[[1]],
    spliceai_ds_al = spliceai_scores[[2]],
    spliceai_ds_dg = spliceai_scores[[3]],
    spliceai_ds_dl = spliceai_scores[[4]],
    spliceai_max = if (all(is.na(spliceai_scores))) {
      NA_real_
    } else {
      max(spliceai_scores, na.rm = TRUE)
    },
    lof = chr_at(picked, "lof")
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
    alphamissense_class = NA_character_,
    transcript = NA_character_,
    gene_id = NA_character_,
    biotype = NA_character_,
    hgvsc = NA_character_,
    hgvsp = NA_character_,
    canonical = NA,
    codons = NA_character_,
    amino_acids = NA_character_,
    cadd_phred = NA_real_,
    cadd_raw = NA_real_,
    revel = NA_real_,
    spliceai_ds_ag = NA_real_,
    spliceai_ds_al = NA_real_,
    spliceai_ds_dg = NA_real_,
    spliceai_ds_dl = NA_real_,
    spliceai_max = NA_real_,
    lof = NA_character_
  )
}

# The alternate allele as VEP reports it after trimming: "C" for an SNV, "A"
# for an insertion sent as T>TA, "-" for a deletion. This is the key VEP uses
# for the per-allele frequencies and clinical significance on a colocated
# variant, so it is read from the element rather than rebuilt from the input.
vep_element_allele <- function(element) {
  alleles <- strsplit(
    biohttp::pluck_at(element, "allele_string", default = ""),
    "/",
    fixed = TRUE
  )[[1]]
  if (length(alleles) < 2) NA_character_ else alleles[[length(alleles)]]
}

# The frequency block for the element's allele. VEP keys `frequencies` by
# allele, and a multi-allelic dbSNP record carries one block per allele it
# has frequencies for.
vep_allele_frequencies <- function(record, element) {
  frequencies <- biohttp::pluck_at(record, "frequencies")
  if (is.null(frequencies) || length(frequencies) == 0) {
    return(NULL)
  }
  allele <- vep_element_allele(element)
  if (!is.na(allele) && !is.null(frequencies[[allele]])) {
    return(frequencies[[allele]])
  }
  if (length(frequencies) == 1) {
    return(frequencies[[1]])
  }
  NULL
}

# The largest frequency among the per-population entries of one source, which
# are the keys carrying the source prefix: gnomadg_afr, gnomadg_nfe and so on.
# The bare gnomadg key is the overall frequency and is left out.
vep_max_population_af <- function(frequencies, prefix) {
  if (is.null(frequencies)) {
    return(NA_real_)
  }
  keys <- names(frequencies)
  keys <- keys[startsWith(keys, prefix)]
  values <- suppressWarnings(as.numeric(unlist(frequencies[keys])))
  values <- values[!is.na(values)]
  if (length(values) == 0) NA_real_ else max(values)
}

# Clinical significance for the element's allele. `clin_sig_allele` is the
# per-allele form, "C:pathogenic;T:benign", and is what tells the alleles of a
# multi-allelic site apart. A deletion's allele is "-" and its prefix in that
# string is empty. `clin_sig` is the union over all alleles and is the fallback
# when the per-allele form is absent.
vep_clin_sig <- function(record, element) {
  per_allele <- chr_at(record, "clin_sig_allele")
  if (!is.na(per_allele) && nzchar(per_allele)) {
    allele <- vep_element_allele(element)
    if (!is.na(allele)) {
      prefix <- if (identical(allele, "-")) "" else allele
      entries <- strsplit(per_allele, ";", fixed = TRUE)[[1]]
      mine <- vapply(
        entries,
        function(entry) {
          pair <- regmatches(
            entry,
            regexpr(":", entry, fixed = TRUE),
            invert = TRUE
          )[[1]]
          if (length(pair) == 2 && identical(pair[[1]], prefix)) {
            pair[[2]]
          } else {
            NA_character_
          }
        },
        character(1),
        USE.NAMES = FALSE
      )
      mine <- unique(mine[!is.na(mine)])
      if (length(mine) > 0) {
        return(paste(mine, collapse = ";"))
      }
    }
  }
  all <- unique(as.character(unlist(
    biohttp::pluck_at(record, "clin_sig", default = list()),
    use.names = FALSE
  )))
  if (length(all) == 0) NA_character_ else paste(all, collapse = ";")
}

#' Turn the colocated variants of a VEP element into a table row
#'
#' Pure. A VEP element lists the known variants at the same site under
#' `colocated_variants`: the dbSNP record, and COSMIC and HGMD entries beside
#' it. The dbSNP record is the one carrying population frequencies and
#' clinical significance, so that is the record read. When several dbSNP
#' records are listed, the one with frequencies wins.
#'
#' Frequencies are only present when `af`, `af_gnomade` or `af_gnomadg` was
#' requested, see [vep_default_options()]. `gnomadg_af_max` and
#' `gnomade_af_max` are the largest per-population frequency of that source,
#' which is what a rarity filter wants rather than the overall frequency.
#' `clin_sig` is the significance of the element's own allele where VEP
#' reports it per allele, so the alleles of a multi-allelic site do not share
#' one answer.
#'
#' @param element One parsed VEP result element.
#'
#' @return A one-row tibble of `rsid`, `gnomadg_af`, `gnomade_af`,
#'   `gnomadg_af_max`, `gnomade_af_max`, and `clin_sig`, all `NA` when the
#'   element has no dbSNP record.
#'
#' @examples
#' element <- list(
#'   allele_string = "G/C",
#'   colocated_variants = list(list(
#'     id = "rs1042522",
#'     frequencies = list(C = list(gnomadg = 0.62, gnomadg_afr = 0.38)),
#'     clin_sig = list("benign")
#'   ))
#' )
#' vep_parse_colocated(element)
#'
#' @export
vep_parse_colocated <- function(element) {
  colocated <- biohttp::pluck_at(element, "colocated_variants")
  if (is.null(colocated) || length(colocated) == 0) {
    return(vep_empty_colocated())
  }
  ids <- vapply(colocated, function(record) chr_at(record, "id"), character(1))
  dbsnp <- colocated[!is.na(ids) & grepl("^rs[0-9]+$", ids)]
  if (length(dbsnp) == 0) {
    return(vep_empty_colocated())
  }
  with_frequencies <- Filter(
    function(record) !is.null(record$frequencies),
    dbsnp
  )
  record <- if (length(with_frequencies) > 0) {
    with_frequencies[[1]]
  } else {
    dbsnp[[1]]
  }
  frequencies <- vep_allele_frequencies(record, element)
  tibble::tibble(
    rsid = chr_at(record, "id"),
    gnomadg_af = num_at(frequencies, "gnomadg"),
    gnomade_af = num_at(frequencies, "gnomade"),
    gnomadg_af_max = vep_max_population_af(frequencies, "gnomadg_"),
    gnomade_af_max = vep_max_population_af(frequencies, "gnomade_"),
    clin_sig = vep_clin_sig(record, element)
  )
}

vep_empty_colocated <- function() {
  tibble::tibble(
    rsid = NA_character_,
    gnomadg_af = NA_real_,
    gnomade_af = NA_real_,
    gnomadg_af_max = NA_real_,
    gnomade_af_max = NA_real_,
    clin_sig = NA_character_
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
#' @return A tibble with one row per entry in `keys`: a `key` column, the
#'   columns of [vep_parse_element()], then those of [vep_parse_colocated()].
#'
#' @export
vep_parse_batch <- function(body, keys) {
  by_key <- list()
  for (element in body) {
    row <- cbind(vep_parse_element(element), vep_parse_colocated(element))
    for (key in vep_element_keys(element)) {
      # First identity wins. The rebuilt key of an indel is a different string
      # from any real input, so it never shadows a correct match.
      if (is.null(by_key[[key]])) {
        by_key[[key]] <- row
      }
    }
  }
  rows <- lapply(keys, function(key) {
    row <- by_key[[key]] %||% cbind(vep_empty_row(), vep_empty_colocated())
    cbind(tibble::tibble(key = as.character(key)), row)
  })
  tibble::as_tibble(do.call(rbind, rows))
}

#' Consequence predictions for many variants
#'
#' @param chrom,pos,ref,alt Variant components, all the same length.
#' @param options A named list of VEP request flags. See
#'   [vep_default_options()] for the supported flags and what each one adds to
#'   the result.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per variant,
#'   in the order asked. See [vep_parse_batch()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(vep_variants("7", 140753336, "A", "T"))
#' biohttp::body_or_null(vep_variants(
#'   "7", 140753336, "A", "T",
#'   options = c(vep_default_options(), list(af_gnomadg = 1, CADD = 1))
#' ))
#' }
#'
#' @export
vep_variants <- function(
  chrom,
  pos,
  ref,
  alt,
  options = vep_default_options(),
  ...
) {
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
    # The defaults include vcf_string=1, which makes VEP report the normalised
    # variant alongside the echoed input, the second identity results are
    # matched on.
    paste0(VEP_URL, "?", vep_query_string(options)),
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

#' Consequence predictions for any number of variants
#'
#' [vep_variants()] refuses more than the 200 VEP accepts per POST. This
#' chunks the input at `chunk_size` and dispatches the chunks through
#' [biohttp::post_json_many()], so a variant table of any length is a handful
#' of requests rather than a loop written at every call site.
#'
#' @section A failed chunk is rows of NA, not a failed call:
#' Following [gnomad_constraints()]. Each chunk's envelope is inspected on its
#' own; a chunk that failed yields a row of `NA` per variant carrying the
#' envelope status in `status`, while the chunks that succeeded are parsed as
#' usual. Partial failure is the normal case across many requests, and taking
#' the whole call down would throw away the answers that did arrive. Only
#' when every chunk failed is the first failing envelope returned, so an
#' outage still reads as one.
#'
#' @inheritParams vep_variants
#' @param chunk_size Variants per request, at most `VEP_BATCH` (200).
#' @param ... Passed to [biohttp::post_json_many()], for example `throttle` or
#'   `max_active`.
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per
#'   variant, in the order asked: the columns of [vep_parse_batch()] plus
#'   `status`, which is `"ok"` for a row whose chunk was answered and the
#'   failing envelope's status otherwise.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(vep_variants_all(
#'   c("7", "17"),
#'   c(140753336, 7676154),
#'   c("A", "G"),
#'   c("T", "C")
#' ))
#' }
#'
#' @export
vep_variants_all <- function(
  chrom,
  pos,
  ref,
  alt,
  chunk_size = VEP_BATCH,
  options = vep_default_options(),
  ...
) {
  n <- length(chrom)
  if (n == 0) {
    return(biohttp::status_no_data(
      source = "VEP",
      detail = "no variants were supplied"
    ))
  }
  if (chunk_size < 1 || chunk_size > VEP_BATCH) {
    stop(
      "chunk_size must be between 1 and ",
      VEP_BATCH,
      ", the verified VEP limit per request",
      call. = FALSE
    )
  }
  # The query string is built once, before any request, so a refused option
  # such as dbNSFP stops the call before anything goes over the wire.
  url <- paste0(VEP_URL, "?", vep_query_string(options))
  regions <- vep_region(chrom, pos, ref, alt)
  keys <- vep_key(chrom, pos, ref, alt)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  bodies <- lapply(chunks, function(rows) {
    # unname() is load-bearing here too. See vep_variants().
    list(variants = unname(as.list(regions[rows])))
  })
  results <- biohttp::post_json_many(
    url,
    bodies = bodies,
    source = "VEP",
    ...
  )
  failed <- Filter(function(res) !isTRUE(res$ok), results)
  if (length(failed) == length(results)) {
    return(failed[[1]])
  }
  parsed <- lapply(seq_along(chunks), function(i) {
    res <- results[[i]]
    chunk_keys <- keys[chunks[[i]]]
    if (!isTRUE(res$ok)) {
      rows <- vep_parse_batch(list(), chunk_keys)
      rows$status <- as.character(res$status %||% "error")
      return(rows)
    }
    rows <- vep_parse_batch(res$data, chunk_keys)
    rows$status <- "ok"
    rows
  })
  biohttp::status_ok(
    data = tibble::as_tibble(do.call(rbind, parsed)),
    source = "VEP"
  )
}
