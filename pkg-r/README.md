# bioclients (R package) <img src="man/figures/logo.png" align="right" height="139" alt="bioclients logo" />

Look up genes, variants and proteins from R. One consistent way to call gnomAD,
ClinVar, UniProt, Ensembl and the other databases you already use, instead of
writing a client for each one.

Each database gets one client, split in two. One half makes the request. The
other turns the answer into a table, and that half needs no network at all, so
it is tested against a saved response.

The repository root `README.md` covers scope and how the repository is laid
out.

## Installation

```r
pak::pak("samuelbharti/bioclients/pkg-r")
```

The `pkg-r` on the end matters. The package sits in a subdirectory, not at the
root of the repository, and an install that leaves it off fails without saying
why.

r-universe works too, and pulls in `biohttp` for you:

```r
install.packages("bioclients", repos = "https://samuelbharti.r-universe.dev")
```

bioclients is not on CRAN yet.
[`biohttp`](https://github.com/samuelbharti/biohttp), the transport underneath
it, is.

## Usage

Every call returns a `biohttp` envelope rather than raising, so you branch on
`res$status` and write no `tryCatch()` of your own.

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

### Nothing found is an answer, not a fault

```r
mygene_gene("NOT_A_REAL_GENE")$status
#> [1] "no_data"
```

`no_data` means the source answered and had nothing. `error`, `timeout` and
`rate_limited` mean the source did not answer. Collapsing the two loses the
difference between "this gene has no published literature" and "Europe PMC is
down", which is the distinction the envelope exists to keep.

### A batch returns one row per input, in input order

```r
res <- mygene_genes(c("TP53", "NOT_A_GENE", "BRAF"))
biohttp::body_or_null(res)
#> # A tibble: 3 x 8
#>   symbol     entrez ...
#>   TP53       7157
#>   NOT_A_GENE NA        <- a row of NA, not a dropped row
#>   BRAF       673
```

A miss is a row of `NA` rather than a missing row, because a shorter table
silently shifts every row after it onto the wrong gene.

### The parser runs without a network

```r
body <- list(hits = list(list(
  symbol = "TP53", name = "tumor protein p53", entrezgene = 7157, HGNC = "11998"
)))
mygene_parse_hits(body, symbol = "TP53")
#> # A tibble: 1 x 8
#>   symbol entrez hgnc  ...
#>   TP53   7157   11998
```

Every client exports its parser separately, so a caller that already holds a
response body never has to make the request again, and every parser is tested
against a stored body rather than a live service.

## Where the traps are written down

Behaviour that cost real debugging time to establish is recorded as a comment at
the call site and pinned by a test. `NEWS.md` lists them together: MyVariant's
silent `notfound` without `assembly=hg38`, IMPC's missing `human_gene_symbol`
field, GTEx's own versioned GENCODE id, MyGene's upper case `HGNC`, and the
rest.

## Contributing

See `CONTRIBUTING.md` at the repository root. Tests are offline, fixtures are
ported rather than regenerated, and a ported fixture that needs editing is a
signal that a parser changed behaviour during the port.
