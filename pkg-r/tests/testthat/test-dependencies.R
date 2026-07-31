# The Suggests split is an architectural requirement, not a tidiness rule: a
# caller that wants one service must not install the dependency closure of every
# other one. See section 5.2 of the build plan.
#
# The Imports-only CI job is the primary enforcement, because it actually
# installs nothing but Imports and then runs the suite. This test is the fast
# local version of the same rule, so a guard that was forgotten shows up during
# development rather than twenty minutes later in CI.

# The package source is present under devtools::test() but not when the tests
# run from a built tarball, so locate it and skip rather than fail.
r_source_dir <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "R"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", "."), "R")
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) NA_character_ else normalizePath(hit[1])
}

description_field <- function(field) {
  path <- testthat::test_path("..", "..", "DESCRIPTION")
  if (!file.exists(path)) {
    return(character())
  }
  raw <- read.dcf(path, fields = field)[1, 1]
  if (is.na(raw)) {
    return(character())
  }
  parts <- trimws(strsplit(raw, ",")[[1]])
  # Strip a version constraint such as "testthat (>= 3.0.0)".
  parts <- sub("\\s*\\(.*\\)$", "", parts)
  parts[nzchar(parts)]
}

test_that("every namespaced call is either an Import or guarded", {
  dir <- r_source_dir()
  skip_if(is.na(dir), "package source not available")

  files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)

  allowed <- c(
    description_field("Imports"),
    description_field("Depends"),
    rownames(installed.packages(priority = "base"))
  )

  violations <- character()

  for (file in files) {
    text <- readLines(file, warn = FALSE)
    # Ignore comment lines so a package named in prose is not a false positive.
    code <- text[!grepl("^\\s*#", text)]
    used <- unique(unlist(regmatches(
      code,
      gregexpr("\\b[a-zA-Z][a-zA-Z0-9.]*(?=:::?)", code, perl = TRUE)
    )))

    for (pkg in setdiff(used, allowed)) {
      guarded <- any(grepl(
        paste0("requireNamespace\\(\\s*[\"']", pkg, "[\"']"),
        code
      ))
      if (!guarded) {
        violations <- c(
          violations,
          paste0(basename(file), " calls ", pkg, "::")
        )
      }
    }
  }

  # One assertion either way, so a run with nothing to flag still counts as a
  # test rather than reporting itself as empty.
  expect_equal(violations, character())
})
