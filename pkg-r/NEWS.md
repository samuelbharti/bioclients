# bioclients 0.0.0.9000

First working version, still in development. Nothing is published, so anything
here can still change.

The heading carries the real version rather than the usual
"(development version)" wording, because `R CMD check` parses `NEWS.md` for a
package name followed by a version and reports "no news entries found" when it
finds prose there instead.

## What it covers

* One client per service, 29 of them, covering every web service the app family
  calls. 125 exported functions in all.
* Every client is split in two. The parser takes an already parsed response body
  and returns a tibble, or `NULL` when there is genuinely nothing there. It
  touches no network, so it is tested directly against a stored response.
* The request half calls `biohttp`, passes a failing envelope straight through,
  turns a `NULL` parse into `no_data`, and wraps a success in `ok`. No client
  writes its own `tryCatch()`.
* Where a service supports it, a batch entry point returns one row per input
  **in input order**, with a row of `NA` for a miss. A shorter table would
  silently shift every row after it onto the wrong gene.

## Fixtures

* 51 stored response bodies, ported byte for byte out of the apps this package
  replaces. A test asserts they still match.
* A ported fixture that needs editing is treated as a signal rather than a
  chore: it means the parser changed behaviour during the port, and that has to
  be justified rather than fixed by editing the recording.

## Verified service behaviour

These cost real debugging time to establish in the apps and are easy to lose in
a port, so each one has a comment at the call site and a test where it is
testable.

* **MyVariant silently returns `notfound`** for GRCh38 coordinates when the
  request omits `assembly=hg38`. Silent, not an error, which is exactly why it
  needs a test.
* **Ensembl VEP takes exactly 200 per POST**, gnomAD exactly 25 GraphQL aliases,
  MyVariant 1000 per batch. Verified ceilings, not tuning.
* **IMPC's `genotype-phenotype` core has no `human_gene_symbol` field**, so
  querying it by symbol returns HTTP 200 with zero documents. A miss is
  `no_data` and never `0`, because a gene never phenotyped and a gene phenotyped
  with no significant result are not the same answer.
* **Europe PMC and PubTator count zero as an answer**, so both come back `ok`. A
  gene with no literature is not an outage.
* **GTEx needs its own versioned GENCODE id**, from GTEx's own reference lookup.
  An Ensembl gene id from anywhere else will not work.
* **MyGene returns the HGNC id under `HGNC`, upper case.** Asking for `hgnc`
  returns nothing at all.
* **Monarch wants the CURIE colon unencoded**, so `HGNC:11998` must reach it
  as-is rather than as `HGNC%3A11998`.
* **PanelApp's index `search` parameter does not filter**, and Reactome answers
  a gene it does not know with a 404 rather than an empty array.
* **STRING replies with `text/json`**, so the content type check has to be off.
* **On the minus strand, Ensembl numbers exons from the highest coordinate**,
  because they are numbered in transcription order rather than by position.

## Credentials

* `clinvar_classification()` picks up `NCBI_API_KEY` from the environment at
  call time and passes it through `biohttp`'s `secret_query`, which attaches it
  at dispatch. The key never reaches the cache key, a printed request, or an
  error message. It is optional; without it NCBI allows 3 requests a second
  instead of 10.

## Known limits

* No live call has been made against Monarch, MyGene's `HGNC` field, Ensembl's
  array versus record shape, Reactome's 404, or PanelApp's non-filtering
  `search`. All five are pinned by tests built from stored responses, but none
  has been seen against a real server since the port.
* `biohttp` is not on CRAN, so `DESCRIPTION` carries a `Remotes:` line.
