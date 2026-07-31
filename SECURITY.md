# Security policy

## Reporting a vulnerability

Please report suspected vulnerabilities privately. Open a GitHub security
advisory on this repository, or contact the maintainer. Do not open a public
issue for a security problem.

## Credentials in requests

This package sends HTTP requests on behalf of a caller, so tokens pass through
it. Headers marked sensitive are passed to `httr2::req_headers(.redact = )`, so
the value is not printed when a request is inspected, logged, or included in an
error message. If you find a path where a redacted header leaks into visible
output, treat it as a security bug and report it as one.

## Keeping secrets out of the repository

- The package itself needs no credentials. Any token is supplied by the caller
  at request time and is never stored or written to the cache.
- Put local values in `.env`, which is gitignored.
- The `.gitignore` ignores whole categories and directories (for example
  `secrets/`, `*.key`, `*.pem`, `.env`) so that no individually sensitive
  filename has to be listed in a tracked file.
- For a file whose name itself would reveal something sensitive, add it to
  `.git/info/exclude`. That file is local and never committed.
- Two backstops run automatically. The `detect-private-key` and
  `detect-aws-credentials` hooks run locally on every commit, and gitleaks scans
  the full history in CI.
