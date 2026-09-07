# cran-comments

## Submission

This is a new submission.

bioclients imports `biohttp`, which is on CRAN at 0.1.2. `DESCRIPTION` asks for
`biohttp (>= 0.1.2)` and has no `Remotes:` field. I ran the tests against the
CRAN build of biohttp, not the newer one on GitHub, so I know nothing here needs
a version you do not have.

## Notes for the reviewer

Three points came back on the first submission of a sibling package,
biobouncer. They are answered here in advance.

* **References.** The package implements no published method of its own. It is
  a client layer over 29 third-party web services, and the references that
  matter are theirs. Citing all 29 in the `Description` field would make it
  unreadable, so five anchors are there in the requested `authors (year)
  <doi:...>` form, with no space after `doi:`: Ensembl, UniProt, gnomAD, Open
  Targets and the AlphaFold Protein Structure Database. The other 24 are on the
  help pages. Every one of the 134 exported topics has a `References` section
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

* **56 examples call a live service, and all of them run.** They are wrapped in
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
enabled: 0 errors | 0 warnings | 1 or 2 notes depending on the run. One is
always the new submission. The other, when it appears, is examples over 5
seconds. From the most recent run:

```
                    user system elapsed
vep_variants       0.096  0.008  20.463
uniprot_features   0.025  0.000  15.026
diseases_channel   0.082  0.003   5.475
```

Which examples show up changes from run to run. Look at the two columns: 0.096
seconds of CPU against 20 seconds elapsed. That time is the service answering,
not this package working. Each one is a single request, so there is nothing
left to make smaller.

If a service is slow enough to pass biohttp's timeout, the example gets a
`timeout` envelope back and still prints without an error. All of these sit in
`\donttest{}`, which the CRAN Cookbook asks for when an example takes more
than 5 seconds. So the note is telling us they are wrapped the right way.

The whole check fits well inside the 10 minute limit. Examples take about two
minutes, tests about one.

Every DOI resolves: the five in the `Description` field and the 29 distinct
ones across the help pages. MyGene and MyVariant share a citation, and
`uniprot_features()` carries two, because it calls the EBI Proteins API rather
than UniProt's own REST service and both deserve naming.

## Test environments

* Local: Windows 11, R 4.6.1, x86_64
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
