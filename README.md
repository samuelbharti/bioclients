# bioclients

<!-- badges: start -->
[![R-CMD-check](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml/badge.svg)](https://github.com/samuelbharti/bioclients/actions/workflows/r.yml)
<!-- badges: end -->

One client per biological database, each with a pure parser that runs offline.

> **Status:** scaffold only. No clients yet. Phase 1 is gated on `biohttp`
> becoming installable.

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

## Dependencies

`Imports` is empty right now and that is deliberate. `biohttp`, `jsonlite`,
`rlang`, and `tibble` all belong there, but `biohttp` is not installable from any
repository until it reaches r-universe, and declaring it now would fail
`R CMD check` for everyone including CI. All four go in together with the first
client.

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

Not published yet, and it will not install until `biohttp` is on r-universe.

## Roadmap

| Phase | State |
| --- | --- |
| 0. Scaffold | in progress |
| 1. Three pilot clients: MyGene, gnomAD, ClinVar | blocked on `biohttp` r-universe |
| 2. Migrate `variant-reviewer` onto `biohttp` and `bioclients` | not started |
| 3. Expand to the remaining services, three at a time | not started |

## License

MIT. See `LICENSE`.
