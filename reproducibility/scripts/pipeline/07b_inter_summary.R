#!/usr/bin/env Rscript
# ============================================================
# 07b_inter_summary.R — summarize enrichment results
# Input: 06b_inter Output
# Output: outputs/<run>/07b_inter_summary/summary.csv
# ============================================================
rm(list = ls())
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1) {
    if (is.na(a)) return(b)
    if (is.character(a) && !nzchar(a)) return(b)
  }
  a
}
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
suppressMessages(library(Triage))

option_list <- list(
  make_option("--in_dir", type = "character", default = ""),
  make_option("--out_dir", type = "character", default = ""),
  make_option("--dataset_name", type = "character", default = "",
              help = "Dataset name for default paths (default: basename(getwd())).")
)

opt <- parse_args(OptionParser(option_list = option_list))
dataset_name <- if (nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- get_dataset_config(dataset_name, project_root)
if (!nzchar(opt$in_dir)) opt$in_dir <- file.path(cfg$llm_outputs_root, "inter", "final_passed")
if (!nzchar(opt$out_dir)) opt$out_dir <- file.path(cfg$llm_outputs_root, "inter", "post_summary")
ensure_dir(opt$out_dir)

fs <- list.files(opt$in_dir, pattern = "\\.json$", full.names = TRUE)
if (length(fs) == 0) {
  stop("No inter outputs found in: ", opt$in_dir)
}

rows <- lapply(fs, function(fp) {
  x <- tryCatch(jsonlite::fromJSON(fp, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(x) || !is.list(x)) return(NULL)
  tibble::tibble(
    cluster_id = as.character(x$cluster_id %||% NA_character_),
    predicted_label = as.character(x$predicted_label %||% NA_character_),
    cl_id = as.character(x$cl_id %||% NA_character_),
    confidence = {
      # 06b LLM Return (High/Medium/Low) or a numeric value
      conf_raw <- tolower(trimws(as.character(x$confidence %||% NA_character_)))
      conf_num <- suppressWarnings(as.numeric(conf_raw))
      if (!is.na(conf_num)) {
        conf_num
      } else {
        switch(conf_raw,
          high = 0.9, medium = 0.7, low = 0.5,
          `very high` = 0.95, `very low` = 0.3,
          NA_real_
        )
      }
    },
    evidence_markers = if (length(x$evidence_markers %||% character(0)) > 0) paste(x$evidence_markers, collapse = "; ") else "",
    evidence_terms = if (length(x$evidence_terms %||% character(0)) > 0) paste(x$evidence_terms, collapse = "; ") else "",
    rationale = as.character(x$rationale %||% NA_character_)
  )
})

out <- dplyr::bind_rows(rows)
readr::write_csv(out, file.path(opt$out_dir, "summary.csv"))
cat("[OK] wrote inter summary: ", file.path(opt$out_dir, "summary.csv"), "\n", sep = "")
