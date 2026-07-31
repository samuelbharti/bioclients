# Internal helpers.
#
# Deliberately thin. `pluck_at()` and `is_blank()` are NOT redefined here:
# biohttp exports both and this package uses those. Every helper duplicated here
# is one more thing that can drift away from the transport layer.

`%||%` <- function(x, y) if (is.null(x)) y else x # nolint: object_name_linter.

# Strip anything that is not a plausible gene-identifier character.
#
# Ported from genescout/R/tools/mygene.R. It is here rather than in mygene.R
# because gnomAD and ClinVar take the same symbols and had their own ad hoc
# versions of this check.
clean_symbol <- function(symbol) {
  if (biohttp::is_blank(symbol)) {
    return(NULL)
  }
  cleaned <- gsub("[^A-Za-z0-9._-]", "", trimws(as.character(symbol)))
  if (cleaned == "") NULL else cleaned
}

# A numeric field from a parsed JSON list. NA by default when absent, but
# counts that are being summed want 0 instead, so the default is settable.
num_at <- function(x, key, default = NA) {
  as.numeric(biohttp::pluck_at(x, key, default = default))
}

# A character field from a parsed JSON list, NA when absent.
chr_at <- function(x, key) {
  as.character(biohttp::pluck_at(x, key, default = NA_character_))
}
