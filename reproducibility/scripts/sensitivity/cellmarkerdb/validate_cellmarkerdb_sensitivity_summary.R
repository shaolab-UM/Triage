#!/usr/bin/env Rscript

# Validate the publication-level CellMarkerDB sensitivity summary.
# No API or LLM call is made.
#
# This script checks the Fig. 3A numeric input from the released
# publication-derived summary data. It does not regenerate the historical
# fourth-reviewer assignments from DEG files.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1L]])
  default
}

repo_root <- get_arg("--repo-root", Sys.getenv("TRIAGE_HOME", unset = getwd()))
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)

summary_file <- file.path(
  repo_root, "reproducibility", "sensitivity",
  "cellmarkerdb_additional_reviewer_summary.tsv"
)
by_dataset_file <- file.path(
  repo_root, "reproducibility", "sensitivity",
  "cellmarkerdb_additional_reviewer_by_dataset.tsv"
)

if (!file.exists(summary_file)) stop("Missing: ", summary_file)
if (!file.exists(by_dataset_file)) stop("Missing: ", by_dataset_file)

s <- read.delim(summary_file, sep = "\t", stringsAsFactors = FALSE)
d <- read.delim(by_dataset_file, sep = "\t", stringsAsFactors = FALSE)

expected <- c(
  "CellMarkerDB added" = 78.7712340295342,
  "CellMarkerDB added with reliability filtering" = 80.1512474030248
)

for (nm in names(expected)) {
  obs <- s$mean_cl_similarity[s$condition == nm]
  if (length(obs) != 1L) stop("Expected exactly one row for: ", nm)
  if (abs(obs - expected[[nm]]) > 1e-8) {
    stop("Publication mean mismatch for ", nm, ": ", obs)
  }

  ds <- d[d$condition == nm, , drop = FALSE]
  if (nrow(ds) != 6L) stop("Expected 5 datasets + Overall for: ", nm)
  calc <- ds$full_triage_mean + ds$delta_vs_full_triage
  if (max(abs(calc - ds$condition_mean)) > 1e-10) {
    stop("Dataset reconstruction mismatch for: ", nm)
  }
}

cat("CellMarkerDB additional-reviewer sensitivity summary\n\n")
print(s, row.names = FALSE, digits = 8)
cat("\nPASS: publication means and Fig. 3A dataset-level deltas are internally consistent.\n")
cat("See docs/CELLMARKERDB_SENSITIVITY.md for the reproducibility boundary.\n")
