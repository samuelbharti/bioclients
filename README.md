# bioclients <img src="pkg-r/man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml/badge.svg)](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml)
<!-- badges: end -->

One client per biological database, each with a pure parser that runs offline.

> **Status:** 29 clients, 125 exported functions, working. Installs from local
> source today; not published anywhere yet.

## Why

The same clients keep getting written across the app family, and the copies have
drifted in ways that are hard to see from inside any one app:

| Service | State today |
| --- | --- |
| gnomAD | Three divergent copies. They do not return the same thing and they do not query the same fields. |
| MyGene | Two copies, in `variant-reviewer` and `genescout`. |
| Disease resolution | Two copies, in `genescout` and `gene-list-builder`. |

`genescout/R/tools/` holds 24 client files with one file per service and the
client kept separate from the parser. That is the best layout in the family and
it is the reference this package follows.

## The shape of a client

Every service module ships two halves:

1. **The client** builds a request, calls into `biohttp`, and returns the
   envelope. It knows URLs, parameters, and rate limits. It touches the network.
2. **The parser** takes an already-parsed body and returns a canonical structure.
   It is pure and never touches the network.

A caller that already has a response body can use the parser alone. A caller that
wants both gets a thin wrapper. Fixture tests exercise the parser directly, with
no network and no mock server.

## What it does not do

- **No ranking or scoring.** `gene-list-builder`'s source-weighted model stays
  there.
- **No curation.** LLM curation and evidence review stay in `genescout` and
  `gene-list-builder`.
- **No disease-to-gene assembly.** That is `gene-list-builder`'s product.
- **No Shiny.** Same rule as `biohttp`.
- **No transport.** Retries, breakers, caching, and error normalization belong to
  `biohttp`.

## Repository layout

The R package is in `pkg-r/`, not at the repository root:

```text
bioclients/
├── .github/workflows/   R CMD check, imports-only, lint, pkgdown, secret scan
├── air.toml .lintr .pre-commit-config.yaml
└── pkg-r/
    ├── DESCRIPTION NAMESPACE
    ├── R/                one file per service, client and parser together
    └── tests/testthat/
        └── fixtures/     stored response bodies, ported from the apps
```

## What is in it

Every service the app family calls now has a client, 29 of them in all. A few of
them, to show the shape every client follows:

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| MyGene | `mygene_gene()` | `mygene_genes()` | `mygene_parse_hits()`, `mygene_parse_batch()` |
| gnomAD constraint | `gnomad_constraint()` | `gnomad_constraints()` | `gnomad_parse_constraint()`, `gnomad_parse_constraints()` |
| gnomAD frequency | `gnomad_frequency()` | | `gnomad_parse_frequency()`, `gnomad_parse_populations()` |
| ClinVar | `clinvar_classification()` | | `clinvar_parse_record()` |
| Open Targets | `opentargets_gene_diseases()`, `opentargets_disease_targets()`, `opentargets_resolve_disease()` | | `opentargets_parse_diseases()`, `opentargets_parse_targets()`, `opentargets_parse_matches()` |
| DGIdb | `dgidb_gene()` | `dgidb_genes()` | `dgidb_parse_genes()` |
| Pharos | `pharos_target()` | `pharos_targets()` | `pharos_parse_targets()` |
| CIViC | `civic_gene()` | | `civic_parse_gene()` |

```r
res <- mygene_genes(c("TP53", "BRCA1", "EGFR"))
biohttp::body_or_null(res)
#> # A tibble: 3 x 7
#>   symbol name              entrez ensembl_gene    uniprot ...
#>   TP53   tumor protein p53 7157   ENSG00000141510 P04637
#>   ...
```

Every client returns a `biohttp` envelope rather than raising, so you branch on
`res$status` and write no `tryCatch()`. Reaching a source that has nothing is
`no_data`, which is an answer rather than a fault.

The batch entry points return one row per input, **in input order**, so you can
zip results onto your inputs by position. An unmatched gene gets a row of `NA`
rather than being dropped, because a shorter table silently shifts every row
after it onto the wrong gene.

### gnomAD has two query types, not one

The three copies of this client in the family disagreed because they answer
different questions: variant frequency, and gene constraint. Both are here, as
separate entry points. Neither is folded into the other, and constraint is not
dropped because two of the three callers wanted frequency.
`gene-list-builder`'s whole ranking model is built on LOEUF, which is
`gnomad_constraint()`.

### Where the line sits on scoring

Two of these clients had a scoring step in their original app copy and it did not
come across. Pharos's source clients mapped its TDL category onto a 0 to 1
weight; this client returns the category. DGIdb's returned a count the app then
weighted; this one returns the count.

Turning a value into a weight is ranking, and ranking is the consuming app's
product. Copying the weights down here would put `gene-list-builder`'s model in
two places, and a change to it would then need a release of this package.

### Absence of evidence is not evidence of absence

`dgidb_parse_genes()` reports `NA` for a gene DGIdb has never heard of and `0`
for a gene it knows with no recorded interactions. `civic_parse_gene()` draws the
same line. The difference matters: collapsing them tells a caller that an unknown
gene is known not to be druggable, which is a much stronger claim than the data
supports.

## Dependencies

`Imports` is `biohttp` and `tibble`.

Everything a single service needs goes in `Suggests`, guarded at the call site
with `requireNamespace()`. This is an architectural requirement, not a cleanup
task: a caller that wants MyGene must not install the dependency closure of the
other services. The `imports-only` CI job is what enforces it, by installing hard
dependencies only and then running the suite.

## No compiled code

The package is pure R. There is no `src/`, and there will not be.

That was evaluated rather than assumed. Measured in `biohttp` against the live
MyGene API, a call spends about 210 ms on the network and 0.2 ms parsing the
response, so parsing is roughly one tenth of one percent of the work. bioclients
adds field extraction on top of that parse, which is smaller still. The largest
stored response in the ported fixture set is 8 KB.

Where the performance actually is, in order:

1. **Concurrent fan-out.** Twelve services queried one after another at ~300 ms
   each is 3.6 seconds; queried together it is closer to 400 ms.
2. **Batch endpoints.** MyGene takes a batch POST, gnomAD's GraphQL API takes
   aliases so many genes fit in one request, and Ensembl VEP takes exactly 200
   per POST. Each turns N round trips into a small fraction of N.
3. **The cache `biohttp` already has**, which a client gets for free by using the
   entry points rather than building requests by hand.

All three are I/O, not compute. If a parsing bottleneck ever does appear, reach
for an existing C or C++ implementation: `yyjsonr` or `RcppSimdJson` for JSON,
`data.table::fread` or `vroom` for the bulk flat-file sources.

## Installation

Not published anywhere. Both packages install from local source:

```r
install.packages("path/to/biohttp", repos = NULL, type = "source")
install.packages("path/to/bioclients/pkg-r", repos = NULL, type = "source")
```

`bioclients` needs `biohttp` installed first. Until `biohttp` is published, the
`R-CMD-check` and `imports-only` CI jobs cannot resolve it and will fail on a
runner. Local `R CMD check` is clean.

## Roadmap

| Phase | State |
| --- | --- |
| 0. Scaffold | done |
| 1. Three pilot clients: MyGene, gnomAD, ClinVar | done |
| 1b. Batch A, the GraphQL services: Open Targets, DGIdb, Pharos, CIViC | done |
| 2. Migrate `variant-reviewer` onto `biohttp` and `bioclients` | not started |
| 3. Expand to the remaining services, three at a time | not started |

## Testing

Offline. No test touches a real host. The parsers run against stored response
bodies ported unchanged from `variant-reviewer`, and a fixture that needed
editing would mean the parser changed behavior during the port.

## License

MIT. See `LICENSE`.
