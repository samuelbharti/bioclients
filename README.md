# bioclients <img src="pkg-r/man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN status](https://www.r-pkg.org/badges/version/bioclients)](https://CRAN.R-project.org/package=bioclients)
[![r-universe](https://samuelbharti.r-universe.dev/badges/bioclients)](https://samuelbharti.r-universe.dev/bioclients)
[![R-CMD-check](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml/badge.svg)](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21770870-1682D4)](https://doi.org/10.5281/zenodo.21770870)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/samuelbharti/bioclients/blob/main/LICENSE)
<!-- badges: end -->

Look up genes, variants and proteins from R. One consistent way to call gnomAD,
ClinVar, UniProt, Ensembl and the other databases you already use, instead of
writing a client for each one.

Documentation is at <https://www.samuelbharti.com/bioclients/>.

## Installation

```r
install.packages("bioclients")
```

For the development version:

```r
pak::pak("samuelbharti/bioclients/pkg-r")
```

r-universe serves prebuilt binaries of the latest release:

```r
install.packages("bioclients", repos = "https://samuelbharti.r-universe.dev")
```

## Usage

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

Ask about many genes at once and the answer comes back one row per input, in the
order you asked:

```r
res <- mygene_genes(c("TP53", "NOT_A_GENE", "BRAF"))
biohttp::body_or_null(res)
#> # A tibble: 3 x 8
#>   symbol     entrez ...
#>   TP53       7157
#>   NOT_A_GENE NA        <- a row of NA, not a dropped row
#>   BRAF       673
```

A miss keeps its row, because a shorter table silently shifts every row after it
onto the wrong gene.

## Motivation

The same clients keep getting written across the app family, and the copies
drift in ways that are hard to see from inside any one app. This package is the
one copy, laid out the way `genescout` already lays them out: one file per
service, with the client kept separate from the parser.

## Scope

bioclients fetches and parses, and stops there. Ranking, scoring, curation and
disease-to-gene assembly belong to the app that calls it. Retries, caching and
error handling belong to [`biohttp`](https://github.com/samuelbharti/biohttp)
underneath. There is no Shiny code here.

## Services

Each service below has a client. The tables group them by what you are looking
up rather than alphabetically.

### Genes

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| MyGene | `mygene_gene()` | `mygene_genes()` | `mygene_parse_hits()`, `mygene_parse_batch()` |

Nearly every other lookup starts here, because the rest are keyed on an Entrez,
Ensembl, UniProt or HGNC id rather than a symbol.

### Variants

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| gnomAD constraint | `gnomad_constraint()` | `gnomad_constraints()` | `gnomad_parse_constraint()`, `gnomad_parse_constraints()` |
| gnomAD frequency | `gnomad_frequency()`, `gnomad_frequency_by_id()` | `gnomad_frequencies()` | `gnomad_parse_frequency()`, `gnomad_parse_variant()`, `gnomad_parse_variants()`, `gnomad_parse_populations()` |
| ClinVar | `clinvar_classification()` | | `clinvar_parse_record()`, `clinvar_category()` |
| MyVariant | | `myvariant_variants()` | `myvariant_parse_record()`, `myvariant_parse_batch()` |
| Ensembl VEP | | `vep_variants()`, `vep_variants_all()` | `vep_parse_element()`, `vep_parse_colocated()`, `vep_parse_batch()` |
| Ensembl REST | `ensembl_vep_id()`, `ensembl_gene_model()` | | `ensembl_parse_vep()`, `ensembl_parse_gene_model()` |
| VariantValidator | `variantvalidator_normalize()` | | `variantvalidator_parse()` |
| ClinGen Allele Registry | | `clingen_alleles()` | `clingen_parse_allele()`, `clingen_parse_batch()` |

### Disease and drugs

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| Open Targets | `opentargets_gene_diseases()`, `opentargets_disease_targets()`, `opentargets_resolve_disease()`, `opentargets_drugs()`, `opentargets_pgx()` | | `opentargets_parse_diseases()`, `opentargets_parse_targets()`, `opentargets_parse_matches()` |
| DGIdb | `dgidb_gene()` | `dgidb_genes()` | `dgidb_parse_genes()` |
| Pharos | `pharos_target()` | `pharos_targets()` | `pharos_parse_targets()` |
| CIViC | `civic_gene()` | | `civic_parse_gene()` |
| JensenLab DISEASES | `diseases_channel()`, `diseases_gene_associations()` | | `diseases_parse_channel()`, `diseases_merge_channels()` |
| PanelApp | `panelapp_panel()`, `panelapp_panels()` | `panelapp_all_panels()` | `panelapp_parse_index()`, `panelapp_parse_panel()` |
| ClinGen gene validity | `clingen_validity_for()` | `clingen_gene_validity()` | `clingen_parse_validity()` |

### Proteins

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| UniProt | `uniprot_diseases()`, `uniprot_features()` | | `uniprot_parse_diseases()`, `uniprot_parse_features()` |
| ProtVar | `protvar_function()`, `protvar_population()` | | `protvar_parse_function()`, `protvar_parse_population()` |
| AlphaFold | `alphafold_model()` | | `alphafold_parse_model()` |
| PDBe | `pdbe_structures()` | | `pdbe_parse_structures()` |
| STRING | `string_partners()`, `string_network()` | `string_map_ids()` | `string_parse_partners()`, `string_parse_network()` |

### Expression and function

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| GTEx | `gtex_median_expression()`, `gtex_gene_reference()` | | `gtex_parse_expression()`, `gtex_parse_reference()` |
| Human Protein Atlas | `hpa_gene()` | | `hpa_parse_gene()` |
| QuickGO | `quickgo_annotations()` | | `quickgo_parse_annotations()` |
| Reactome | `reactome_pathways()` | | `reactome_parse_pathways()` |

### Phenotypes

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| HPO | `hpo_search()`, `hpo_term()`, `hpo_gene_annotation()` | | `hpo_parse_search()`, `hpo_parse_term()`, `hpo_parse_diseases()` |
| Monarch | `monarch_search()`, `monarch_associations()`, `monarch_gene_phenotypes()` | | `monarch_parse_search()`, `monarch_parse_associations()` |
| IMPC | `impc_mouse_ortholog()`, `impc_gene_phenotypes()` | | `impc_parse_ortholog()`, `impc_parse_phenotypes()` |

### Literature

| Service | Single call | Batch call | Parser |
| --- | --- | --- | --- |
| Europe PMC | `europepmc_search()`, `europepmc_count()` | | `europepmc_parse_results()`, `europepmc_parse_count()` |
| PubTator3 | `pubtator_gene_literature()` | | `pubtator_parse_results()`, `pubtator_parse_count()` |

`pkg-r/_pkgdown.yml` holds the same grouping for the reference index, and
`NEWS.md` lists the verified behaviour each client pins.

## Citation

The package is archived on Zenodo. Use the concept DOI, which always resolves to
the newest archived release:

> Bharti, S. (2026). *bioclients: Clients for Biological Database Web Services*.
> Zenodo. <https://doi.org/10.5281/zenodo.21770870>

`CITATION.cff` carries the same metadata and one identifier per archived
version, so `citation("bioclients")` in R and the "Cite this repository" button
on GitHub both work.

## Acknowledgements

Barret Schloerke and Carson Sievert advise this work as thesis advisors.
Posit Software, PBC funded early work on this package and holds copyright
together with the author.

## License

MIT. See `LICENSE`.
