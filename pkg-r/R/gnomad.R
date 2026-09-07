# gnomAD: population frequency and gene constraint.
#
# Reconciled from the family's three divergent copies, which disagreed because
# they answer different questions:
#
#   variant-reviewer/R/api_gnomad.R   variant frequency, by rsID
#   a second app                      gene constraint (LOEUF), batched
#   a third app                       frequency in batch
#
# BOTH query types are exposed, as separate entry points. One function cannot
# serve both, and constraint must not be dropped because two of the three
# callers happen not to use it: the second app's whole ranking model is built on
# LOEUF.
#
# Endpoint: https://gnomad.broadinstitute.org/api (GraphQL)

GNOMAD_URL <- "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET <- "gnomad_r4"

# gnomAD's real limit is a query COST cap, not a request count, and the verified
# ceiling is exactly 25 aliases: a 26th returns "Query is too expensive". That
# was established against the live API and is recorded in the source registry of
# the app that verified it, as a hard limit rather than tuning.
#
# 20 rather than 25 on purpose. Sitting on the exact ceiling is how a batch
# starts failing for a reason that reads like rate limiting and is not, and
# the headroom costs one request in twenty-five. Raising this needs a live
# re-check, not a guess.
#
# Re-verified for the aliased variant query on 2026-08-30: 25 aliases of
# `variant()` with the whole GNOMAD_VARIANT_FIELDS selection set answered, and
# 26 came back HTTP 400 with "Query is too expensive (26). Maximum allowed
# cost is 25." So the cost is one per alias regardless of how many fields the
# alias selects, and the same chunk size serves both batched queries.
GNOMAD_CHUNK <- 20

# The groups gnomAD leaves out when it reports a group maximum frequency:
# the bottlenecked ancestries, whose founder effects inflate the frequency of
# their own variants, and the remaining/other bucket, which is not a group.
# The API does not serve grpmax itself, so gnomad_variant_row() derives it
# from the per-group counts and applies the same exclusion.
GNOMAD_GRPMAX_EXCLUDED <- c("ami", "asj", "fin", "mid", "remaining")

# Display labels for gnomAD's genetic-ancestry group codes.
GNOMAD_POP_LABELS <- c(
  afr = "African / African-American",
  ami = "Amish",
  amr = "Admixed American",
  asj = "Ashkenazi Jewish",
  eas = "East Asian",
  fin = "European (Finnish)",
  mid = "Middle Eastern",
  nfe = "European (non-Finnish)",
  sas = "South Asian",
  remaining = "Remaining"
)

# --- Gene constraint ---------------------------------------------------------

#' Turn a gnomAD constraint response into a table
#'
#' Pure. Takes an already-parsed response body and never touches the network.
#'
#' `loeuf` is gnomAD's `oe_lof_upper`. The two names are the same number, and
#' the field is called `oe_lof_upper` in the API but LOEUF everywhere else,
#' including in the ranking models that consume it. Both names appear here so a
#' reader of either can find it.
#'
#' @param body A parsed gnomAD GraphQL response body.
#' @param symbol The gene symbol that was queried.
#'
#' @return A one-row tibble with `symbol`, `pli`, `loeuf`, `oe_lof`, `oe_mis`,
#'   `mis_z`, `syn_z`, and `lof_z`. `NULL` when the gene has no constraint
#'   block, which is common and is an answer rather than a fault.
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' body <- list(data = list(gene = list(
#'   gnomad_constraint = list(pli = 1, oe_lof_upper = 0.23)
#' )))
#' gnomad_parse_constraint(body, "BRAF")
#'
#' @export
gnomad_parse_constraint <- function(body, symbol = NA_character_) {
  constraint <- biohttp::pluck_at(body, "data", "gene", "gnomad_constraint")
  gnomad_constraint_row(constraint, symbol)
}

# One constraint block to one tibble row.
gnomad_constraint_row <- function(constraint, symbol) {
  if (is.null(constraint)) {
    return(NULL)
  }
  tibble::tibble(
    symbol = as.character(symbol),
    pli = num_at(constraint, "pli"),
    loeuf = num_at(constraint, "oe_lof_upper"),
    oe_lof = num_at(constraint, "oe_lof"),
    oe_mis = num_at(constraint, "oe_mis"),
    mis_z = num_at(constraint, "mis_z"),
    syn_z = num_at(constraint, "syn_z"),
    lof_z = num_at(constraint, "lof_z")
  )
}

#' Gene constraint for one gene
#'
#' How intolerant a gene is to variation. `pli` and `loeuf` summarize
#' loss-of-function intolerance; the Z-scores summarize missense and synonymous
#' depletion.
#'
#' @param symbol A gene symbol.
#' @param reference_genome `"GRCh38"` or `"GRCh37"`.
#' @param ... Passed to [biohttp::post_json()], for example `throttle`.
#'
#' @return A biohttp envelope whose `data` is a one-row tibble. See
#'   [gnomad_parse_constraint()].
#'
#' @references
#' Chen et al. (2024). A genomic mutational constraint map using variation
#' in 76,156 human genomes. Nature 625(7993), 92-100.
#' \doi{10.1038/s41586-023-06045-0}
#'
#' Service documentation: <https://gnomad.broadinstitute.org/>
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(gnomad_constraint("BRAF"))
#' }
#'
#' @export
gnomad_constraint <- function(
  symbol,
  reference_genome = "GRCh38",
  ...
) {
  cleaned <- clean_symbol(symbol)
  if (is.null(cleaned)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = "no usable gene symbol was supplied"
    ))
  }
  query <- sprintf(
    paste(
      "query($sym: String!) {",
      "  gene(gene_symbol: $sym, reference_genome: %s) {",
      "    gnomad_constraint {",
      "      pli oe_lof oe_lof_upper mis_z syn_z oe_mis lof_z",
      "    }",
      "  }",
      "}",
      sep = "\n"
    ),
    reference_genome
  )
  res <- biohttp::post_json(
    GNOMAD_URL,
    body = list(query = query, variables = list(sym = cleaned)),
    source = "gnomAD",
    ...
  )
  # A GraphQL 200 carrying an errors array is a failure, and biohttp already
  # knows how to say so.
  bad <- biohttp::graphql_error(res, "gnomAD")
  if (!is.null(bad)) {
    return(bad)
  }
  parsed <- gnomad_parse_constraint(res$data, cleaned)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      http = res$http,
      detail = paste0("no gnomAD constraint for ", cleaned)
    ))
  }
  biohttp::status_ok(data = parsed, source = "gnomAD", http = res$http)
}

# Build one aliased GraphQL query covering a chunk of symbols. gnomAD takes one
# gene per `gene()` call, so many genes in one request means aliases.
gnomad_alias_query <- function(symbols, reference_genome) {
  aliases <- vapply(
    seq_along(symbols),
    function(i) {
      sprintf(
        paste0(
          'g%d: gene(gene_symbol: "%s", reference_genome: %s) ',
          "{ symbol gnomad_constraint ",
          "{ pli oe_lof oe_lof_upper mis_z syn_z oe_mis lof_z } }"
        ),
        i,
        symbols[[i]],
        reference_genome
      )
    },
    character(1)
  )
  paste0("{\n", paste(aliases, collapse = "\n"), "\n}")
}

#' Turn an aliased gnomAD constraint response into a table
#'
#' Pure. The batched query uses GraphQL aliases `g1`, `g2`, and so on, so the
#' response is keyed by alias rather than by symbol. Rows are mapped back by
#' alias index rather than by the echoed symbol, so the result stays aligned
#' even for a gene gnomAD does not return a symbol for.
#'
#' @param body A parsed gnomAD GraphQL response body.
#' @param symbols The gene symbols that were queried, in the order asked.
#'
#' @return A tibble with one row per entry in `symbols`, same order. A gene with
#'   no constraint block gets a row of `NA` rather than being dropped.
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' body <- list(data = list(
#'   g1 = list(symbol = "BRAF", gnomad_constraint = list(oe_lof_upper = 0.23)),
#'   g2 = list(symbol = "TP53", gnomad_constraint = NULL)
#' ))
#' gnomad_parse_constraints(body, c("BRAF", "TP53"))
#'
#' @export
gnomad_parse_constraints <- function(body, symbols) {
  data <- biohttp::pluck_at(body, "data")
  rows <- lapply(seq_along(symbols), function(i) {
    entry <- biohttp::pluck_at(data, paste0("g", i))
    row <- gnomad_constraint_row(
      biohttp::pluck_at(entry, "gnomad_constraint"),
      symbols[[i]]
    )
    row %||% gnomad_empty_constraint_row(symbols[[i]])
  })
  do.call(rbind, rows)
}

gnomad_empty_constraint_row <- function(symbol) {
  tibble::tibble(
    symbol = as.character(symbol),
    pli = NA_real_,
    loeuf = NA_real_,
    oe_lof = NA_real_,
    oe_mis = NA_real_,
    mis_z = NA_real_,
    syn_z = NA_real_,
    lof_z = NA_real_
  )
}

#' Gene constraint for many genes
#'
#' Batched through GraphQL aliases, so a gene list costs a handful of requests
#' rather than one per gene. Chunked at `chunk_size` to stay under gnomAD's
#' query cost cap of 25. See `GNOMAD_CHUNK`.
#'
#' @param symbols Gene symbols.
#' @param chunk_size Genes per request.
#' @inheritParams gnomad_constraint
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry in
#'   `symbols`, in the same order.
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(gnomad_constraints(c("BRAF", "TP53", "EGFR")))
#' }
#'
#' @export
gnomad_constraints <- function(
  symbols,
  reference_genome = "GRCh38",
  chunk_size = GNOMAD_CHUNK,
  ...
) {
  cleaned <- vapply(
    symbols,
    function(symbol) clean_symbol(symbol) %||% NA_character_,
    character(1),
    USE.NAMES = FALSE
  )
  if (all(is.na(cleaned))) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = "no usable gene symbols were supplied"
    ))
  }
  usable <- cleaned[!is.na(cleaned)]
  chunks <- split(usable, ceiling(seq_along(usable) / chunk_size))
  bodies <- lapply(
    chunks,
    function(chunk) {
      list(query = gnomad_alias_query(chunk, reference_genome))
    }
  )
  # One batched call for the chunks. Only the chunks the cache is missing go
  # over the wire.
  results <- biohttp::post_json_many(
    GNOMAD_URL,
    bodies = bodies,
    source = "gnomAD",
    ...
  )
  parsed <- lapply(seq_along(chunks), function(i) {
    res <- results[[i]]
    bad <- biohttp::graphql_error(res, "gnomAD")
    if (!is.null(bad)) {
      # A failed chunk yields NA rows rather than taking the whole call down.
      # Partial failure is the normal case across many genes.
      return(do.call(rbind, lapply(chunks[[i]], gnomad_empty_constraint_row)))
    }
    gnomad_parse_constraints(res$data, chunks[[i]])
  })
  table <- do.call(rbind, parsed)
  # Map back onto the caller's input, including the tokens that were unusable.
  out <- do.call(
    rbind,
    lapply(cleaned, function(symbol) {
      if (is.na(symbol)) {
        return(gnomad_empty_constraint_row(NA_character_))
      }
      table[match(symbol, table$symbol), , drop = FALSE]
    })
  )
  biohttp::status_ok(data = out, source = "gnomAD")
}

# --- Variant frequency -------------------------------------------------------

#' Combine gnomAD per-ancestry counts into one frequency table
#'
#' Pure. Sums the exome and genome sample sets per ancestry group.
#'
#' Sex-split ids such as `nfe_XX` and the bare `XX`/`XY` breakdowns are dropped,
#' by keeping only the known ancestry codes, so the result is one row per
#' ancestry rather than a mix of ancestries and sexes.
#'
#' @param exome_pops,genome_pops The `populations` lists from each sample set.
#'
#' @return A tibble of `pop`, `label`, `ac`, `an`, `af`, sorted by frequency
#'   descending. `NULL` when there is nothing to report.
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' gnomad_parse_populations(
#'   list(list(id = "nfe", ac = 3, an = 1000)),
#'   list(list(id = "nfe", ac = 1, an = 500))
#' )
#'
#' @export
gnomad_parse_populations <- function(exome_pops, genome_pops) {
  # Exome and genome counts for the same ancestry group add together, so both
  # lists fold into one accumulator. The accumulator is passed in and returned
  # rather than reached for with `<<-`, because a closure quietly reassigning a
  # variable in its parent is the kind of thing that reads fine and then
  # surprises whoever moves the function.
  add <- function(acc, pops) {
    for (pop in pops) {
      id <- biohttp::pluck_at(pop, "id")
      if (biohttp::is_blank(id) || !(id %in% names(GNOMAD_POP_LABELS))) {
        next
      }
      previous <- acc[[id]] %||% c(0, 0)
      acc[[id]] <- c(
        previous[[1]] + num_at(pop, "ac", default = 0),
        previous[[2]] + num_at(pop, "an", default = 0)
      )
    }
    acc
  }
  acc <- add(list(), exome_pops)
  acc <- add(acc, genome_pops)
  if (length(acc) == 0) {
    return(NULL)
  }
  ids <- names(acc)
  # Unnamed up front. acc is a named list, so vapply carries the ancestry code
  # through as a names attribute, and a named column reads back as a named
  # vector for every downstream caller.
  counts <- unname(vapply(acc, function(v) v[[1]], numeric(1)))
  totals <- unname(vapply(acc, function(v) v[[2]], numeric(1)))
  out <- tibble::tibble(
    pop = ids,
    label = unname(GNOMAD_POP_LABELS[ids]),
    ac = counts,
    an = totals,
    af = ifelse(totals > 0, counts / totals, NA_real_)
  )
  out$label[is.na(out$label)] <- out$pop[is.na(out$label)]
  out <- out[out$an > 0, , drop = FALSE]
  if (nrow(out) == 0) {
    return(NULL)
  }
  out[order(-out$af), , drop = FALSE]
}

# Normalize one frequency block (exome or genome), or NULL when gnomAD has no
# data for that sample set.
gnomad_freq_part <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  list(af = num_at(x, "af"), ac = num_at(x, "ac"), an = num_at(x, "an"))
}

#' Turn a gnomAD variant response into a frequency record
#'
#' Pure.
#'
#' @param body A parsed gnomAD GraphQL response body.
#' @param dataset The dataset that was queried, carried through to the result.
#'
#' @return A list of `variant_id`, `dataset`, `exome`, `genome`, and
#'   `populations`, or `NULL` when the body carries no variant.
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' body <- list(data = list(variant = list(
#'   variant_id = "7-140753336-A-T",
#'   exome = list(af = 0.001, ac = 2, an = 2000)
#' )))
#' gnomad_parse_frequency(body)$variant_id
#'
#' @export
gnomad_parse_frequency <- function(body, dataset = GNOMAD_DATASET) {
  variant <- biohttp::pluck_at(body, "data", "variant")
  if (is.null(variant)) {
    return(NULL)
  }
  list(
    variant_id = chr_at(variant, "variant_id"),
    dataset = dataset,
    exome = gnomad_freq_part(biohttp::pluck_at(variant, "exome")),
    genome = gnomad_freq_part(biohttp::pluck_at(variant, "genome")),
    populations = gnomad_parse_populations(
      biohttp::pluck_at(variant, "exome", "populations"),
      biohttp::pluck_at(variant, "genome", "populations")
    )
  )
}

#' Population allele frequency for a variant
#'
#' Looked up by rsID, which is what gnomAD's variant query takes.
#'
#' @section An rsID does not identify an allele:
#' `7:g.140753336A>T` and `7:g.140753336A>C` both map to `rs113488022`. Any
#' lookup routed through an rsID is therefore lossy, and it fails **silently**
#' while looking entirely plausible. The canonical key is
#' `(assembly, chromosome, position, ref, alt)`.
#'
#' This function is faithful to what gnomAD's variant query accepts, so the
#' limitation is gnomAD's rather than this package's. Know about it before you
#' rely on the answer for a multi-allelic site.
#'
#' @param rsid A dbSNP rsID.
#' @param dataset The gnomAD dataset.
#' @param ... Passed to [biohttp::post_json()].
#'
#' @return A biohttp envelope whose `data` is the list described in
#'   [gnomad_parse_frequency()].
#'
#' @inherit gnomad_constraint references
#'
#' @examples
#' \donttest{
#' biohttp::body_or_null(gnomad_frequency("rs113488022"))
#' }
#'
#' @export
gnomad_frequency <- function(rsid, dataset = GNOMAD_DATASET, ...) {
  if (biohttp::is_blank(rsid)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = "no rsID was supplied"
    ))
  }
  query <- sprintf(
    paste(
      "query($rsid: String!) {",
      "  variant(rsid: $rsid, dataset: %s) {",
      "    variant_id",
      "    exome { af ac an populations { id ac an } }",
      "    genome { af ac an populations { id ac an } }",
      "  }",
      "}",
      sep = "\n"
    ),
    dataset
  )
  res <- biohttp::post_json(
    GNOMAD_URL,
    body = list(query = query, variables = list(rsid = as.character(rsid))),
    source = "gnomAD",
    ...
  )
  bad <- biohttp::graphql_error(res, "gnomAD")
  if (!is.null(bad)) {
    return(bad)
  }
  parsed <- gnomad_parse_frequency(res$data, dataset)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      http = res$http,
      detail = paste0("gnomAD has no record for ", rsid)
    ))
  }
  biohttp::status_ok(data = parsed, source = "gnomAD", http = res$http)
}

# --- Variant frequency by id -------------------------------------------------

# The selection set for one variant, shared by the single and the batched
# query so the two cannot drift. `homozygote_count` is what gnomAD calls
# nhomalt. The faf95 block is the filtering allele frequency, and `filters` is
# the list of VCF filters the site failed, empty when it passed.
GNOMAD_VARIANT_FIELDS <- paste(
  "variant_id rsids",
  "exome { af ac an homozygote_count filters",
  "  faf95 { popmax popmax_population } populations { id ac an } }",
  "genome { af ac an homozygote_count filters",
  "  faf95 { popmax popmax_population } populations { id ac an } }"
)

#' Build a gnomAD variant id from variant components
#'
#' Pure. gnomAD's variant query takes `chrom-pos-ref-alt`, with no `chr`
#' prefix and the alleles in upper case. Vectorised over its arguments.
#'
#' @inheritParams vep_region
#'
#' @return A character vector of ids such as `"1-55516888-G-GA"`.
#'
#' @examples
#' gnomad_variant_id("chr1", 55516888, "g", "ga")
#'
#' @export
gnomad_variant_id <- function(chrom, pos, ref, alt) {
  paste(
    sub("^chr", "", as.character(chrom), ignore.case = TRUE),
    as.integer(pos),
    toupper(as.character(ref)),
    toupper(as.character(alt)),
    sep = "-"
  )
}

# The assembly a dataset's ids are on. The variant query is keyed by dataset
# alone, so this is what lets a caller's reference_genome be checked against
# it rather than silently ignored.
gnomad_dataset_genome <- function(dataset) {
  if (grepl("^gnomad_r2", dataset)) "GRCh37" else "GRCh38"
}

# One variant object to one tibble row.
gnomad_variant_row <- function(variant, variant_id) {
  if (is.null(variant)) {
    return(NULL)
  }
  exome <- biohttp::pluck_at(variant, "exome")
  genome <- biohttp::pluck_at(variant, "genome")
  rsids <- unlist(biohttp::pluck_at(variant, "rsids"), use.names = FALSE)
  # The group maximum, derived: exome and genome counts fold together per
  # ancestry group, the bottlenecked groups are left out, and the group with
  # the highest frequency is reported. gnomad_parse_populations() returns the
  # table sorted by frequency, so the first surviving row is the answer.
  populations <- gnomad_parse_populations(
    biohttp::pluck_at(exome, "populations"),
    biohttp::pluck_at(genome, "populations")
  )
  grpmax <- NULL
  if (!is.null(populations)) {
    eligible <- populations[!(populations$pop %in% GNOMAD_GRPMAX_EXCLUDED), ]
    if (nrow(eligible) > 0) {
      grpmax <- eligible[1, ]
    }
  }
  faf <- gnomad_faf95(exome, genome)
  filters <- unique(c(
    unlist(biohttp::pluck_at(exome, "filters"), use.names = FALSE),
    unlist(biohttp::pluck_at(genome, "filters"), use.names = FALSE)
  ))
  tibble::tibble(
    variant_id = as.character(variant_id),
    rsid = if (length(rsids) == 0) NA_character_ else as.character(rsids[[1]]),
    exome_af = num_at(exome, "af"),
    exome_ac = num_at(exome, "ac"),
    exome_an = num_at(exome, "an"),
    exome_nhomalt = num_at(exome, "homozygote_count"),
    genome_af = num_at(genome, "af"),
    genome_ac = num_at(genome, "ac"),
    genome_an = num_at(genome, "an"),
    genome_nhomalt = num_at(genome, "homozygote_count"),
    grpmax_af = if (is.null(grpmax)) NA_real_ else grpmax$af,
    grpmax_an = if (is.null(grpmax)) NA_real_ else grpmax$an,
    grpmax_id = if (is.null(grpmax)) NA_character_ else grpmax$pop,
    faf95 = faf$value,
    faf95_pop = faf$pop,
    filters = if (length(filters) == 0) {
      NA_character_
    } else {
      paste(filters, collapse = ";")
    }
  )
}

# The higher of the exome and genome filtering allele frequencies, with the
# group it belongs to.
gnomad_faf95 <- function(exome, genome) {
  candidates <- list(
    biohttp::pluck_at(exome, "faf95"),
    biohttp::pluck_at(genome, "faf95")
  )
  best <- list(value = NA_real_, pop = NA_character_)
  for (candidate in candidates) {
    value <- num_at(candidate, "popmax")
    if (!is.na(value) && (is.na(best$value) || value > best$value)) {
      best <- list(value = value, pop = chr_at(candidate, "popmax_population"))
    }
  }
  best
}

gnomad_empty_variant_row <- function(variant_id) {
  tibble::tibble(
    variant_id = as.character(variant_id),
    rsid = NA_character_,
    exome_af = NA_real_,
    exome_ac = NA_real_,
    exome_an = NA_real_,
    exome_nhomalt = NA_real_,
    genome_af = NA_real_,
    genome_ac = NA_real_,
    genome_an = NA_real_,
    genome_nhomalt = NA_real_,
    grpmax_af = NA_real_,
    grpmax_an = NA_real_,
    grpmax_id = NA_character_,
    faf95 = NA_real_,
    faf95_pop = NA_character_,
    filters = NA_character_
  )
}

#' Turn a gnomAD variant response into a frequency row
#'
#' Pure. The flat, one-row form of a variant record, for a variant table that
#' wants a frequency per row. [gnomad_parse_frequency()] is the nested form
#' with the per-ancestry table.
#'
#' @section grpmax is derived:
#' The API does not serve a group maximum. It is computed here from the
#' per-group counts, exome and genome summed, leaving out the bottlenecked
#' groups and the remaining bucket the same way gnomAD does. See
#' `GNOMAD_GRPMAX_EXCLUDED`.
#'
#' @param body A parsed gnomAD GraphQL response body.
#' @param variant_id The id that was queried, carried through to the row.
#'
#' @return A one-row tibble of `variant_id`, `rsid`, `exome_af`, `exome_ac`,
#'   `exome_an`, `exome_nhomalt`, `genome_af`, `genome_ac`, `genome_an`,
#'   `genome_nhomalt`, `grpmax_af`, `grpmax_an`, `grpmax_id`, `faf95`,
#'   `faf95_pop`, and `filters`. `NULL` when the body carries no variant,
#'   which is how gnomAD answers for a variant it has never seen.
#'
#' @examples
#' body <- list(data = list(variant = list(
#'   variant_id = "17-7676154-G-C",
#'   rsids = list("rs1042522"),
#'   exome = list(af = 0.72, ac = 1046941, an = 1461558, homozygote_count = 380188)
#' )))
#' gnomad_parse_variant(body, "17-7676154-G-C")
#'
#' @export
gnomad_parse_variant <- function(body, variant_id = NA_character_) {
  variant <- biohttp::pluck_at(body, "data", "variant")
  gnomad_variant_row(variant, variant_id)
}

# An errors array that only says "Variant not found" is not a failed query.
# gnomAD answers an absent variant with a 200, `variant: null` in data, and
# that message in errors, in both the single and the aliased form. Treating
# it as biohttp::graphql_error() does would turn every rare variant into an
# error envelope. Anything else in the array is a real failure.
gnomad_variant_error <- function(res) {
  if (!isTRUE(res$ok)) {
    return(res)
  }
  errors <- res$data$errors
  if (is.null(errors) || length(errors) == 0) {
    return(NULL)
  }
  messages <- vapply(errors, function(e) chr_at(e, "message"), character(1))
  if (all(!is.na(messages) & messages == "Variant not found")) {
    return(NULL)
  }
  biohttp::graphql_error(res, "gnomAD")
}

#' Population allele frequency for a variant, by id
#'
#' Looked up by `chrom-pos-ref-alt`, which names one allele exactly. This is
#' the lookup [gnomad_frequency()] cannot do, because an rsID is shared by
#' every allele at a site.
#'
#' @param variant_id A gnomAD variant id, see [gnomad_variant_id()].
#' @param dataset The gnomAD dataset. `gnomad_r4` is GRCh38; the `gnomad_r2`
#'   datasets are GRCh37.
#' @param reference_genome The assembly the id is on. The variant query is
#'   keyed by dataset alone, so this is checked against `dataset` and a
#'   mismatch is refused rather than sent, because gnomAD would answer with
#'   whatever sits at those coordinates on the other assembly.
#' @param ... Passed to [biohttp::post_json()].
#'
#' @return A biohttp envelope whose `data` is the one-row tibble described in
#'   [gnomad_parse_variant()]. `no_data` when gnomAD has no record of the
#'   variant.
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(gnomad_frequency_by_id("17-7676154-G-C"))
#' }
#'
#' @export
gnomad_frequency_by_id <- function(
  variant_id,
  dataset = GNOMAD_DATASET,
  reference_genome = "GRCh38",
  ...
) {
  if (biohttp::is_blank(variant_id)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = "no variant id was supplied"
    ))
  }
  if (!identical(gnomad_dataset_genome(dataset), reference_genome)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = paste0(dataset, " is not on ", reference_genome)
    ))
  }
  query <- sprintf(
    paste(
      "query($id: String!) {",
      "  variant(variantId: $id, dataset: %s) {",
      "    %s",
      "  }",
      "}",
      sep = "\n"
    ),
    dataset,
    GNOMAD_VARIANT_FIELDS
  )
  res <- biohttp::post_json(
    GNOMAD_URL,
    body = list(query = query, variables = list(id = as.character(variant_id))),
    source = "gnomAD",
    ...
  )
  bad <- gnomad_variant_error(res)
  if (!is.null(bad)) {
    return(bad)
  }
  parsed <- gnomad_parse_variant(res$data, variant_id)
  if (is.null(parsed)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      http = res$http,
      detail = paste0("gnomAD has no record for ", variant_id)
    ))
  }
  biohttp::status_ok(data = parsed, source = "gnomAD", http = res$http)
}

# One aliased query covering a chunk of variant ids, v1, v2, and so on.
gnomad_variant_alias_query <- function(variant_ids, dataset) {
  aliases <- vapply(
    seq_along(variant_ids),
    function(i) {
      sprintf(
        'v%d: variant(variantId: "%s", dataset: %s) { %s }',
        i,
        variant_ids[[i]],
        dataset,
        GNOMAD_VARIANT_FIELDS
      )
    },
    character(1)
  )
  paste0("{\n", paste(aliases, collapse = "\n"), "\n}")
}

#' Turn an aliased gnomAD variant response into a table
#'
#' Pure. The batched query uses aliases `v1`, `v2`, and so on, so rows are
#' mapped back by alias index rather than by the echoed id. A variant gnomAD
#' has not seen comes back as `null` under its alias, with a
#' "Variant not found" entry in `errors`, and becomes a row of `NA`.
#'
#' @param body A parsed gnomAD GraphQL response body.
#' @param variant_ids The ids that were queried, in the order asked.
#'
#' @return A tibble with one row per entry in `variant_ids`, same order. See
#'   [gnomad_parse_variant()] for the columns.
#'
#' @examples
#' body <- list(data = list(
#'   v1 = list(variant_id = "17-7676154-G-C", exome = list(af = 0.72)),
#'   v2 = NULL
#' ))
#' gnomad_parse_variants(body, c("17-7676154-G-C", "17-7676154-G-GTTTTT"))
#'
#' @export
gnomad_parse_variants <- function(body, variant_ids) {
  data <- biohttp::pluck_at(body, "data")
  rows <- lapply(seq_along(variant_ids), function(i) {
    row <- gnomad_variant_row(
      biohttp::pluck_at(data, paste0("v", i)),
      variant_ids[[i]]
    )
    row %||% gnomad_empty_variant_row(variant_ids[[i]])
  })
  do.call(rbind, rows)
}

#' Population allele frequency for many variants
#'
#' Batched through GraphQL aliases and chunked at `chunk_size` to stay under
#' gnomAD's query cost cap of 25, which was verified for this query. See
#' `GNOMAD_CHUNK`. Dispatched through [biohttp::post_json_many()], so only
#' the chunks the cache is missing go over the wire.
#'
#' A failed chunk yields a row of `NA` per variant rather than taking the
#' whole call down, following [gnomad_constraints()]. A variant gnomAD has no
#' record of is a row of `NA` too, because that is an answer.
#'
#' @param variant_ids gnomAD variant ids, see [gnomad_variant_id()].
#' @param chunk_size Variants per request.
#' @inheritParams gnomad_frequency_by_id
#' @param ... Passed to [biohttp::post_json_many()].
#'
#' @return A biohttp envelope whose `data` is a tibble with one row per entry
#'   in `variant_ids`, in the same order. See [gnomad_parse_variant()].
#'
#' @examples
#' \dontrun{
#' biohttp::body_or_null(gnomad_frequencies(
#'   c("17-7676154-G-C", "7-117559590-ATCT-A")
#' ))
#' }
#'
#' @export
gnomad_frequencies <- function(
  variant_ids,
  dataset = GNOMAD_DATASET,
  reference_genome = "GRCh38",
  chunk_size = GNOMAD_CHUNK,
  ...
) {
  ids <- as.character(variant_ids %||% character())
  usable <- unique(ids[!is.na(ids) & nzchar(trimws(ids))])
  if (length(usable) == 0) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = "no variant ids were supplied"
    ))
  }
  if (!identical(gnomad_dataset_genome(dataset), reference_genome)) {
    return(biohttp::status_no_data(
      source = "gnomAD",
      detail = paste0(dataset, " is not on ", reference_genome)
    ))
  }
  chunks <- split(usable, ceiling(seq_along(usable) / chunk_size))
  bodies <- lapply(chunks, function(chunk) {
    list(query = gnomad_variant_alias_query(chunk, dataset))
  })
  results <- biohttp::post_json_many(
    GNOMAD_URL,
    bodies = bodies,
    source = "gnomAD",
    ...
  )
  parsed <- lapply(seq_along(chunks), function(i) {
    res <- results[[i]]
    bad <- gnomad_variant_error(res)
    if (!is.null(bad)) {
      return(do.call(rbind, lapply(chunks[[i]], gnomad_empty_variant_row)))
    }
    gnomad_parse_variants(res$data, chunks[[i]])
  })
  table <- do.call(rbind, parsed)
  # Map back onto the caller's input, including blanks and repeats.
  out <- do.call(
    rbind,
    lapply(ids, function(id) {
      hit <- match(id, table$variant_id)
      if (is.na(hit)) {
        return(gnomad_empty_variant_row(id))
      }
      table[hit, , drop = FALSE]
    })
  )
  biohttp::status_ok(data = out, source = "gnomAD")
}
