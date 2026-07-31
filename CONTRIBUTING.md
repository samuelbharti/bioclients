# Contributing to bioclients

Thanks for helping. This guide covers the workflow and the local tooling.

## Repository layout

The R package lives in `pkg-r/`, not at the repository root. Tooling and
community files sit at the root and apply to the whole repository. Run R
commands from `pkg-r/`:

```sh
cd pkg-r
Rscript -e 'devtools::check()'
```

## What belongs here

bioclients is one client per service and nothing else. Before opening a pull
request, check that the change fits inside these lines:

- **No ranking or scoring.** Anything that weights, scores, or orders genes
  belongs to `gene-list-builder`. Its source-weighted model is the intellectual
  core of that app.
- **No curation.** LLM curation, citation gating, and evidence review stay in
  `genescout` and `gene-list-builder`.
- **No disease-to-gene assembly.** That is `gene-list-builder`'s product.
- **No Shiny**, in Imports, in Suggests, or in tests.
- **No transport.** Retries, breakers, caching, and error normalization belong
  to `biohttp`. A client reaching for `httr2` directly is a sign something is
  wrong.

## The client and parser split

Every service module ships two halves, and the split is the point:

1. **The client** builds a request, calls into `biohttp`, and returns the
   envelope. It knows URLs, parameters, and rate limits. It touches the network.
2. **The parser** takes an already-parsed body and returns a canonical
   structure. It is pure and never touches the network. This is what the fixture
   tests exercise.

A caller that already has a response body can use the parser alone. Keep them
separable, and do not fold one into the other for convenience.

## Dependencies

`Imports` stays small: `biohttp` for transport, plus the handful of packages
every client needs.

Everything a single service needs goes in `Suggests`, guarded at the call site
with `requireNamespace()` and a message naming the missing package. This is an
architectural requirement, not a tidiness rule. A caller that wants MyGene must
not install the dependency closure of the other services.

Two things enforce it: the `imports-only` CI job, which installs hard
dependencies only and then runs the suite, and `tests/testthat/test-dependencies.R`,
which is the fast local version of the same check. A test that needs a suggested
package must call `skip_if_not_installed()` rather than assume it is there.

## No compiled code

The package is pure R and stays that way. There is no `src/`, no C, no C++, and
no Rust.

This was measured rather than assumed, in `biohttp`. Against the live MyGene API
a call spends about 210 ms on the network and 0.2 ms parsing the response, so
parsing is roughly one tenth of one percent of the work. bioclients adds field
extraction on top of that parse, which is smaller still. Every entry point here
waits on somebody else's server, so a faster native parser optimizes the wrong
end of the call, and native code would put a toolchain requirement into every
consuming app's Docker build.

If a parsing bottleneck ever does appear, reach for an existing C or C++
implementation first: `yyjsonr` or `RcppSimdJson` for JSON, `data.table::fread`
or `vroom` for the bulk flat-file sources. All of them install everywhere today.

## Branches and commits

- `dev` is the integration branch. **Every pull request targets `dev`**, not
  `main`. `dev` is merged into `main` at a release.
- Both `main` and `dev` are protected. Do not commit to either directly. The
  `no-commit-to-branch` hook blocks it locally.
- Name branches with a type prefix: `feat/<slug>`, `fix/<slug>`, or
  `chore/<slug>`.
- Use Conventional Commit messages, for example `feat: add the MyGene client`.
  Keep commits small and focused. The commit-msg hook checks the format.
- The PR title also follows Conventional Commits. A CI check enforces it.

## Local setup

Install the git hooks once:

```sh
prek install --install-hooks
prek install --hook-type commit-msg
```

Then before every push:

```sh
prek run --all-files
```

The hooks run air for R formatting plus a set of general checks. CI additionally
runs `R CMD check` on five platforms, the imports-only job, lintr, and gitleaks.

## Tests

Tests are offline. No test hits a real host, and CI runs with no network. Parsers
are tested directly against a stored response body; clients are tested against a
`webfakes` server.

Two rules about fixtures:

1. **Fixtures are ported, not regenerated.** They come from the apps this package
   replaces and they should pass unchanged.
2. **A ported fixture that needs editing is a signal, not a chore.** It means the
   parser changed behavior during the port, and that needs justifying in the pull
   request rather than fixing by editing the fixture.

Verified API behavior gets a comment at the call site and a test where it is
testable. These cost real debugging time to establish and they are easy to lose
in a port. The clearest example: MyVariant silently returns `notfound` for GRCh38
coordinates if the request omits `assembly=hg38`. Silent, not an error, which is
exactly why it needs a test.

## Secrets

Never commit secrets. Put local values in `.env`, which is gitignored. For a file
whose name is itself sensitive, add it to `.git/info/exclude`, which is never
committed, rather than to `.gitignore`. See `SECURITY.md`.
