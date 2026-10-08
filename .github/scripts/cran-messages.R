# Fail the CRAN check job on any message R CMD check gives that is not on the
# known list below.
#
# Every warning and every note is read, and each one is split into its
# separate messages first. R puts several findings into one warning: the
# incoming feasibility check, for one, can report the package version, a link
# and a ROR id all at once. Judging the warning as a whole would let one known
# message hide the others, which is what happened in October 2026 when a
# "possibly invalid ROR IDs" message passed under the known version message.
#
# Notes are read too, because a dead link or DOI is a note, not a warning.
#
# Run: Rscript .github/scripts/cran-messages.R path/to/00check.log

# Messages expected for reasons that are not defects. Each is matched against
# the start of one message.
known <- c(
  # Between releases the repository holds the version CRAN already has.
  "Insufficient package version",
  # The examples call live services, so their elapsed time is the services'
  # time. cran-comments.md explains this to CRAN.
  "Examples with CPU (user + system) or elapsed time > 5s"
)

# One warning or note, as rcmdcheck gives it, to its separate messages. The
# first line names the check. After it, findings are paragraphs split by blank
# lines. The incoming check always opens with a "Maintainer:" paragraph, which
# is not a finding. A note with no detail is judged by its first line.
messages_of <- function(block) {
  lines <- strsplit(block, "\n", fixed = TRUE)[[1]]
  body <- lines[-1]
  blank <- !nzchar(trimws(body))
  paragraphs <- vapply(
    split(body[!blank], cumsum(blank)[!blank]),
    paste,
    character(1),
    collapse = "\n"
  )
  paragraphs <- paragraphs[!startsWith(paragraphs, "Maintainer:")]
  if (length(paragraphs) == 0) {
    paragraphs <- lines[[1]]
  }
  data.frame(check = lines[[1]], message = unname(paragraphs))
}

log <- commandArgs(trailingOnly = TRUE)[[1]]
res <- rcmdcheck::parse_check(log)

# rcmdcheck drops a warning or note that has no detail under it, so a count of
# what was read can fall short of what the check reported. The last line of
# the log, "Status: 1 WARNING, 2 NOTEs", is the check's own count. If it
# reports more than was read, fail: something was not looked at.
status <- grep("^Status:", readLines(log, warn = FALSE), value = TRUE)
count_in_status <- function(kind) {
  hit <- regmatches(status, regexpr(paste0("[0-9]+ ", kind), status))
  if (length(hit) == 0) 0L else as.integer(sub(" .*", "", hit))
}
reported <- c(
  warnings = count_in_status("WARNING"),
  notes = count_in_status("NOTE")
)
read <- c(warnings = length(res$warnings), notes = length(res$notes))
if (any(reported > read)) {
  cat(
    "The check reported ",
    reported[["warnings"]],
    " warning(s) and ",
    reported[["notes"]],
    " note(s), but only ",
    read[["warnings"]],
    " and ",
    read[["notes"]],
    " could be read. Look at ",
    log,
    " directly.\n",
    sep = ""
  )
  quit(status = 1)
}

found <- do.call(rbind, lapply(c(res$warnings, res$notes), messages_of))

if (is.null(found)) {
  cat("No warnings or notes.\n")
  quit(status = 0)
}

found$known <- vapply(
  found$message,
  function(message) any(startsWith(message, known)),
  logical(1)
)

for (i in seq_len(nrow(found))) {
  cat(
    if (found$known[[i]]) "known:   " else "UNKNOWN: ",
    found$check[[i]],
    "\n  ",
    gsub("\n", "\n  ", found$message[[i]], fixed = TRUE),
    "\n\n",
    sep = ""
  )
}

if (any(!found$known)) {
  cat(sum(!found$known), "message(s) not on the known list. See above.\n")
  quit(status = 1)
}
cat("Every message is on the known list.\n")
