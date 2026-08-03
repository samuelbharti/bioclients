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
# 20 rather than 25 on purpose. Cost is charged per field, not only per alias,
# and this query asks for eight fields per gene where the client that verified
# the 25 asked for two. Sitting on the exact ceiling with a wider selection set
# is how a batch starts failing for a reason that reads like rate limiting and
# is not. Raising this needs a live re-check, not a guess.
GNOMAD_CHUNK <- 20

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
#' @examples
#' \dontrun{
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
#' @examples
#' \dontrun{
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
#' @examples
#' \dontrun{
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
