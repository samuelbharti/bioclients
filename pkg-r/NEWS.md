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

Everything below has now been confirmed against the live service, not only
against a stored response. `tests/testthat/test-live.R` is where that happens.
It calls all 29 services once and then makes the specific assertions the claims
below depend on. It never runs on its own, and no CI job turns it on:

```sh
BIOCLIENTS_LIVE=true Rscript -e 'devtools::test("pkg-r", filter = "live")'
```

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
* **Ensembl takes `content-type` as a query parameter, not a header.** Leave it
  out and Ensembl serves its HTML browser page with HTTP 200, so the call looks
  like a success until the JSON parser reaches the first tag.
* **On the minus strand, Ensembl numbers exons from the highest coordinate**,
  because they are numbered in transcription order rather than by position.

## Credentials

* `clinvar_classification()` picks up `NCBI_API_KEY` from the environment at
  call time and passes it through `biohttp`'s `secret_query`, which attaches it
  at dispatch. The key never reaches the cache key, a printed request, or an
  error message. It is optional; without it NCBI allows 3 requests a second
  instead of 10.

## What the live run changed

Four of the five behaviours that had never been checked against a real server
held exactly as the port described them: MyGene's upper case `HGNC`, Reactome's
404, PanelApp's `search` that does not filter, and Ensembl's array wrapping a
single record.

The fifth did not, and it was the one with an action already written into it.
Monarch was being called on two hosts, `api-v3.monarchinitiative.org` for search
and association and `api.monarchinitiative.org` for the entity route, because
that was how the apps did it and nobody had checked whether they were
interchangeable. They are. All three routes answer on both, and the entity
payload is identical between them down to the order of the ids. So Monarch is
one host now, which also means one circuit breaker instead of two. A host going
down used to open only half of them.

## Known limits

* `biohttp` is not on CRAN, so `DESCRIPTION` carries a `Remotes:` line.
* Pharos was returning HTTP 502 from its own gateway throughout the live run, so
  its probe is the one service the run could not confirm. That is an outage
  rather than a finding about the port, and the check is left in place to say so
  the next time it is run.
