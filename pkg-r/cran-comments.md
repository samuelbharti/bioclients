# cran-comments

## Submission

This is a new submission.

**This package is not ready to be sent yet, and this file records why.**
bioclients imports `biohttp`, which is not on CRAN. Until it is, `DESCRIPTION`
carries a `Remotes:` line pointing at its GitHub repository, which is what lets
a clean CI runner and r-universe resolve the dependency. That line comes out,
and this paragraph with it, at the point biohttp is accepted. Everything else
below is the submission as it stands.

## Notes for the reviewer

Three points came back on the first submission of a sibling package,
biobouncer. They are answered here in advance.

* **References.** The package implements no published method of its own. It is
  a client layer over 29 third-party web services, and the references that
  matter are theirs. Citing all 29 in the `Description` field would make it
  unreadable, so five anchors are there in the requested `authors (year)
  <doi:...>` form, with no space after `doi:`: Ensembl, UniProt, gnomAD, Open
  Targets and the AlphaFold Protein Structure Database. The other 24 are on the
  help pages. Every one of the 125 exported topics has a `References` section
  naming the service it calls, its canonical publication with a DOI, and the
  service's own documentation URL. Each DOI was resolved through the CrossRef
  API and checked against the title, first author, year and pagination.

* **Commented-out example code.** There is none. No example in the package
  contains a commented line.

* **Examples needing a package from `Suggests`.** There are none, so no
  `requireNamespace()` guard is needed. Examples call only `biohttp` and
  `tibble`, both in `Imports`. `Suggests` lists jsonlite, knitr, rmarkdown,
  testthat and withr, and none of them is called from an example.

Two further points about how the examples are structured.

* **53 examples call a live service, and all of them run.** They are wrapped in
  `\donttest{}` rather than `\dontrun{}` because they are executable: every
  service is public and none needs a key. They are safe to run during a check
  for the reason the package exists. A client returns a value describing the
  outcome instead of raising, so with no network reachable the example gets a
  failed envelope, `biohttp::body_or_null()` returns `NULL`, and the example
  prints `NULL` with no error and no warning. This was verified by running the
  examples behind an unroutable proxy.

* **Nothing is written to the user's filespace.** biohttp's disk cache is
  opt-in and off by default, so an example writes nothing anywhere. The check
  directory is clean afterwards.

## R CMD check results

`R CMD check --as-cran --run-donttest`, with incoming and remote checks
enabled: 0 errors | 1 warning | 1 note.

The warning is the `Remotes:` field and the unavailable `biohttp` dependency
described at the top. It goes away when that line does.

The note is two examples over 5 seconds:

```
                      user system elapsed
ensembl_gene_model   0.036  0.003  15.031
vep_variants         0.019  0.002  20.017
```

Both are single requests to Ensembl, and the elapsed time is entirely Ensembl's
own latency; the CPU time is a rounding error. `ensembl_gene_model()` fetches
one gene model and `vep_variants()` posts one variant, so there is nothing left
to make smaller. Repeated timing puts them at 10 to 20 seconds depending on how
Ensembl is feeling, and a run that exceeds biohttp's timeout returns a `timeout`
envelope rather than failing the example. Both are in `\donttest{}`, which is
the wrapper the CRAN Cookbook prescribes for an example over 5 seconds, so the
note is reporting that they are correctly wrapped.

All five DOIs in the `Description` field resolve, as do the 29 on the help
pages.

## Test environments

* Local: macOS 26.5.2, R 4.6.0, aarch64
* GitHub Actions, on the pull request into `main`: ubuntu-latest (R-devel,
  R-release, R-oldrel-1), windows-latest (R-release), macos-latest (R-release),
  plus a job that installs hard dependencies only and runs the suite with no
  `Suggests` present

## Tests

The default suite is offline. No test reaches a real host: parsers run against
stored response bodies and clients run against `httr2::local_mocked_responses()`.
One file, `tests/testthat/test-live.R`, calls the real services, and it is gated
three ways, on `skip_on_cran()`, on `skip_if_offline()`, and on the
`BIOCLIENTS_LIVE` environment variable that no CI job sets.
