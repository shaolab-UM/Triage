#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(jsonlite))

root <- Sys.getenv("TRIAGE_HOME", unset = getwd())

expected <- c(
  "reproducibility/primary/Census_immune/final" = 16L,
  "reproducibility/primary/Sikkema_lung/final" = 11L,
  "reproducibility/primary/TS_kidney/final" = 4L,
  "reproducibility/primary/TS_pancreas/final" = 10L,
  "reproducibility/primary/Zheng_blood/final" = 9L,
  "reproducibility/external/Anderson_DLPFC/final" = 18L,
  "reproducibility/external/Zha_AD_mouse/final" = 9L,
  "reproducibility/external/S6K1_organoid/final" = 15L
)

required_top <- c(
  "schema_version", "dataset", "cluster_id", "final_decision",
  "reviewers", "evidence", "qc", "evaluation", "execution", "provenance"
)
forbidden_top <- c(
  "decision_trace", "reviewer_support", "adjudication_flags", "manual_review_plan",
  "shadow_pool", "fallback_diagnostics"
)

forbidden_reviewer <- c(
  "final_cl_id", "support_class", "is_correct", "matches_final_label"
)

failed <- FALSE
all_ids <- character(0)

for (rel in names(expected)) {
  d <- file.path(root, rel)
  files <- if (dir.exists(d)) sort(list.files(d, pattern = "^cluster_.*\\.json$", full.names = TRUE)) else character(0)
  if (length(files) != expected[[rel]]) {
    cat("COUNT MISMATCH:", rel, length(files), "/", expected[[rel]], "\n")
    failed <- TRUE
  }
  for (fp in files) {
    x <- jsonlite::fromJSON(fp, simplifyVector = FALSE)
    miss <- setdiff(required_top, names(x))
    bad <- intersect(forbidden_top, names(x))
    if (length(miss)) {
      cat("MISSING FIELDS:", fp, paste(miss, collapse = ","), "\n")
      failed <- TRUE
    }
    if (length(bad)) {
      cat("FORBIDDEN PUBLIC FIELDS:", fp, paste(bad, collapse = ","), "\n")
      failed <- TRUE
    }
    for (reviewer_name in names(x$reviewers)) {
      bad_reviewer <- intersect(names(x$reviewers[[reviewer_name]]), forbidden_reviewer)
      if (length(bad_reviewer)) {
        cat("FORBIDDEN REVIEWER-DERIVED FIELDS:", fp, reviewer_name,
            paste(bad_reviewer, collapse = ","), "\n")
        failed <- TRUE
      }
    }
    if (is.null(x$evaluation$reference_cl_id)) {
      cat("MISSING REFERENCE CL:", fp, "\n")
      failed <- TRUE
    }
    all_ids <- c(all_ids, paste(x$dataset, x$cluster_id, sep = "::"))
  }
}

if (anyDuplicated(all_ids)) {
  cat("DUPLICATE DATASET/CLUSTER IDs detected\n")
  failed <- TRUE
}

if (failed) quit(status = 1L)
cat("PASS: public release structure, counts and schema checks completed.\n")
