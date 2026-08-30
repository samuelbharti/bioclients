# bioclients <img src="pkg-r/man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml/badge.svg)](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml)
[![r-universe](https://samuelbharti.r-universe.dev/badges/bioclients)](https://samuelbharti.r-universe.dev/bioclients)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21770870.svg)](https://doi.org/10.5281/zenodo.21770870)
<!-- badges: end -->

One client per biological database, each with a pure parser that runs offline.

> **Status:** 0.1.0, released. Read the docs at
> <https://www.samuelbharti.com/bioclients/>. The parser output shape is what can
> still move; a change to an existing column is a breaking change.

## Installation

Neither package is on CRAN. Both are on r-universe, which pulls `biohttp` in as a
dependency, so one call is enough:

```r
install.packages("bioclients", repos = "https://samuelbharti.r-universe.dev")
```

From GitHub instead, note the `subdir`. The package sits in `pkg-r/` rather than
at the repository root, and an install that leaves this out fails without saying
why:

```r
pak::pak("samuelbharti/bioclients/pkg-r")
# or
remotes::install_github("samuelbharti/bioclients", subdir = "pkg-r")
```

`DESCRIPTION` carries a `Remotes:` line pointing at
[`biohttp`](https://github.com/samuelbharti/biohttp) on GitHub, which is what
lets a clean CI runner resolve it without any credentials. That line has to come
back out before any CRAN submission.

## A first call

```r
library(bioclients)

res <- mygene_gene("TP53")
res$status
#> [1] "ok"

biohttp::body_or_null(res)
#> # A tibble: 1 x 8
#>   symbol name              summary entrez ensembl_gene    uniprot hgnc  type_of_gene
#>   TP53   tumor protein p53 This g. 7157   ENSG00000141510 P04637  11998 protein-cod.
```

Every client returns an envelope rather than raising, so you branch on
`res$status` and write no `tryCatch()` of your own. A source that has nothing to
say answers `no_data`, which is an answer rather than a fault.

## Why

The same clients keep getting written across the app family, and the copies have
drifted in ways that are hard to see from inside any one app:

| Service | State today |
| --- | --- |
| gnomAD | Three divergent copies. They do not return the same thing and they do not query the same fields. |
| MyGene | Two copies, in `variant-reviewer` and `genescout`. |
| Disease resolution | Two copies, in `genescout` and one other app. |

[`genescout`](https://github.com/samuelbharti/genescout)`/R/tools/` holds 24
client files with one file per service and the client kept separate from the
parser. That is the best layout in the family and it is the reference this
package follows.

## The shape of a client

Every service module ships two halves:

1. **The client** builds a request, calls into
   [`biohttp`](https://github.com/samuelbharti/biohttp), and returns the
   envelope. It knows URLs, parameters, and rate limits. It touches the network.
2. **The parser** takes an already-parsed body and returns a canonical structure.
   It is pure and never touches the network.

A caller that already has a response body can use the parser alone. A caller that
wants both gets a thin wrapper. Fixture tests exercise the parser directly, with
no network and no mock server.

## What it does not do

- **No ranking or scoring.** A source-weighted model belongs to the app that
  holds it.
- **No curation.** LLM curation and evidence review stay in the consuming apps.
- **No disease-to-gene assembly.** That is a consuming app's product.
- **No Shiny.** Same rule as `biohttp`.
- **No transport.** Retries, breakers, caching, and error normalization belong to
  `biohttp`.

## Repository layout

The R package is in `pkg-r/`, not at the repository root:

```text
bioclients/
├── .github/workflows/   R CMD check, imports-only, lint, pkgdown, secret scan
├── air.toml .lintr .pre-commit-config.yaml
├── README.md            this file: why it exists, scope, roadmap
└── pkg-r/
    ├── DESCRIPTION NAMESPACE NEWS.md
    ├── README.md         the package itself: install and usage
    ├── R/                one file per service, client and parser together
    ├── man/figures/      the hex logo
    └── tests/testthat/
        └── fixtures/     stored response bodies, ported from the apps
```

The workflows only run on a pull request into `main`, and on the push to `main`
that merges one. Everything else is checked locally through the prek hooks. See
`CONTRIBUTING.md`.

## What is in it

Every service the app family calls now has a client, 29 of them in all, grouped
below by the question being asked rather than alphabetically, because a caller
arrives knowing what they want to look up and not which service answers it.

### Which gene is this

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| MyGene | `mygene_gene()` | `mygene_genes()` | `mygene_parse_hits()`, `mygene_parse_batch()` |

Nearly every other lookup starts here, because the rest are keyed on an Entrez,
Ensembl, UniProt or HGNC id rather than a symbol.

### What is known about this variant

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| gnomAD constraint | `gnomad_constraint()` | `gnomad_constraints()` | `gnomad_parse_constraint()`, `gnomad_parse_constraints()` |
| gnomAD frequency | `gnomad_frequency()` | | `gnomad_parse_frequency()`, `gnomad_parse_populations()` |
| ClinVar | `clinvar_classification()` | | `clinvar_parse_record()`, `clinvar_category()` |
| MyVariant | | `myvariant_variants()` | `myvariant_parse_record()`, `myvariant_parse_batch()` |
| Ensembl VEP | | `vep_variants()` | `vep_parse_element()`, `vep_parse_colocated()`, `vep_parse_batch()` |
| Ensembl REST | `ensembl_vep_id()`, `ensembl_gene_model()` | | `ensembl_parse_vep()`, `ensembl_parse_gene_model()` |
| VariantValidator | `variantvalidator_normalize()` | | `variantvalidator_parse()` |
| ClinGen Allele Registry | | `clingen_alleles()` | `clingen_parse_allele()`, `clingen_parse_batch()` |

### Is this gene linked to disease, and is it druggable

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| Open Targets | `opentargets_gene_diseases()`, `opentargets_disease_targets()`, `opentargets_resolve_disease()`, `opentargets_drugs()`, `opentargets_pgx()` | | `opentargets_parse_diseases()`, `opentargets_parse_targets()`, `opentargets_parse_matches()` |
| DGIdb | `dgidb_gene()` | `dgidb_genes()` | `dgidb_parse_genes()` |
| Pharos | `pharos_target()` | `pharos_targets()` | `pharos_parse_targets()` |
| CIViC | `civic_gene()` | | `civic_parse_gene()` |
| JensenLab DISEASES | `diseases_channel()`, `diseases_gene_associations()` | | `diseases_parse_channel()`, `diseases_merge_channels()` |
| PanelApp | `panelapp_panel()`, `panelapp_panels()` | `panelapp_all_panels()` | `panelapp_parse_index()`, `panelapp_parse_panel()` |
| ClinGen gene validity | `clingen_validity_for()` | `clingen_gene_validity()` | `clingen_parse_validity()` |

### What does the protein look like

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| UniProt | `uniprot_diseases()`, `uniprot_features()` | | `uniprot_parse_diseases()`, `uniprot_parse_features()` |
| ProtVar | `protvar_function()`, `protvar_population()` | | `protvar_parse_function()`, `protvar_parse_population()` |
| AlphaFold | `alphafold_model()` | | `alphafold_parse_model()` |
| PDBe | `pdbe_structures()` | | `pdbe_parse_structures()` |
| STRING | `string_partners()`, `string_network()` | `string_map_ids()` | `string_parse_partners()`, `string_parse_network()` |

### Where is it expressed, and what does it do

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| GTEx | `gtex_median_expression()`, `gtex_gene_reference()` | | `gtex_parse_expression()`, `gtex_parse_reference()` |
| Human Protein Atlas | `hpa_gene()` | | `hpa_parse_gene()` |
| QuickGO | `quickgo_annotations()` | | `quickgo_parse_annotations()` |
| Reactome | `reactome_pathways()` | | `reactome_parse_pathways()` |

### What phenotype does it cause

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| HPO | `hpo_search()`, `hpo_term()`, `hpo_gene_annotation()` | | `hpo_parse_search()`, `hpo_parse_term()`, `hpo_parse_diseases()` |
| Monarch | `monarch_search()`, `monarch_associations()`, `monarch_gene_phenotypes()` | | `monarch_parse_search()`, `monarch_parse_associations()` |
| IMPC | `impc_mouse_ortholog()`, `impc_gene_phenotypes()` | | `impc_parse_ortholog()`, `impc_parse_phenotypes()` |

### Who has written about it

| Service | Ask about one | Ask about many | Pure parser |
| --- | --- | --- | --- |
| Europe PMC | `europepmc_search()`, `europepmc_count()` | | `europepmc_parse_results()`, `europepmc_parse_count()` |
| PubTator3 | `pubtator_gene_literature()` | | `pubtator_parse_results()`, `pubtator_parse_count()` |

`pkg-r/_pkgdown.yml` holds the same grouping for the reference index, and
`NEWS.md` lists the verified behaviour each client pins.

```r
res <- mygene_genes(c("TP53", "BRCA1", "EGFR"))
biohttp::body_or_null(res)
#> # A tibble: 3 x 7
#>   symbol name              entrez ensembl_gene    uniprot ...
#>   TP53   tumor protein p53 7157   ENSG00000141510 P04637
#>   ...
```

The batch entry points return one row per input, **in input order**, so you can
zip results onto your inputs by position. An unmatched gene gets a row of `NA`
rather than being dropped, because a shorter table silently shifts every row
after it onto the wrong gene.

### gnomAD has two query types, not one

The three copies of this client in the family disagreed because they answer
different questions: variant frequency, and gene constraint. Both are here, as
separate entry points. Neither is folded into the other, and constraint is not
dropped because two of the three callers wanted frequency. One caller's whole
ranking model is built on LOEUF, which is `gnomad_constraint()`.

### Where the line sits on scoring

Two of these clients had a scoring step in their original app copy and it did not
come across. Pharos's source clients mapped its TDL category onto a 0 to 1
weight; this client returns the category. DGIdb's returned a count the app then
weighted; this one returns the count.

Turning a value into a weight is ranking, and ranking is the consuming app's
product. Copying the weights down here would put that model in two places, and a
change to it would then need a release of this package.

### Absence of evidence is not evidence of absence

`dgidb_parse_genes()` reports `NA` for a gene DGIdb has never heard of and `0`
for a gene it knows with no recorded interactions. `civic_parse_gene()` draws the
same line. The difference matters: collapsing them tells a caller that an unknown
gene is known not to be druggable, which is a much stronger claim than the data
supports.

## Dependencies

`Imports` is `biohttp` for transport, `httr2` for the one client that has to
assemble its own request, and `tibble`.

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

## Roadmap

| Phase | State |
| --- | --- |
| 0. Scaffold | done |
| 1. Three pilot clients: MyGene, gnomAD, ClinVar | done |
| 1b. Batch A, the GraphQL services: Open Targets, DGIdb, Pharos, CIViC | done |
| 2. Expand to every remaining service | done, 29 clients and 125 exports |
| 3. Confirm the ported behaviour against live services | done, 28 of 29 confirmed |
| 4. Migrate `variant-reviewer` onto `biohttp` and `bioclients` | demonstrated, not landed |

Phase 3 called every service once and then made the specific assertions the five
unchecked claims depended on. Four held: MyGene's upper case `HGNC`, Reactome's
404, PanelApp's `search` that does not filter, and Ensembl's array wrapping a
single record. The fifth did not, and usefully so. Monarch's two hosts turned out
to serve all three routes with identical payloads, so they collapsed to one, and
Monarch now has one circuit breaker rather than two.

Pharos is the one service the run could not confirm, because it was answering
HTTP 502 from its own gateway at the time. That is an outage rather than a
finding.

The checks live in `pkg-r/tests/testthat/test-live.R` and never run on their own.
No CI job sets the variable that turns them on.

```sh
BIOCLIENTS_LIVE=true Rscript -e 'devtools::test("pkg-r", filter = "live")'
```

Phase 4 has a working demonstration. Swapping `variant-reviewer`'s HTTP layer for
`biohttp` touched one file, left all eleven of its API clients unchanged, and
passed its full suite at 340 of 340 with live calls confirmed against four
services. It also surfaced three real bugs in `biohttp`, since fixed in 0.1.1.

## Testing

Offline by default. The parsers run against stored response bodies ported
unchanged from the apps this package replaces, mostly `genescout` and
`variant-reviewer` with a smaller number from two other apps in the family. A
fixture that needed editing would mean the parser changed behaviour during the
port.

One file is the exception. `pkg-r/tests/testthat/test-live.R` calls real
services, which is the only way to answer whether a ported claim is still true.
It is gated three ways, on `BIOCLIENTS_LIVE`, on not being CRAN, and on having a
network, so it skips unless it is asked for by name. No CI job sets the variable.

## Citing bioclients

Each release is archived on Zenodo. Use the concept DOI, which always resolves to
the newest release:

> Bharti, S. (2026). *bioclients: Clients for Biological Database Web Services*.
> Zenodo. <https://doi.org/10.5281/zenodo.21770870>

To pin the exact version you used, cite its own DOI instead. Version 0.1.0 is
[10.5281/zenodo.21770871](https://doi.org/10.5281/zenodo.21770871).

`CITATION.cff` carries the same metadata, so `citation("bioclients")` in R and
the "Cite this repository" button on GitHub both work.

## License

MIT. See `LICENSE`.
