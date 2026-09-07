#!/usr/bin/env Rscript
# ============================================================
# Build post_summary/summary.csv from judge final outputs
# Input:  llm_judge_outputs/final/
# Output: llm_judge_outputs/judge_post_summary/summary.csv
#
# - PASSED clusters: parse *_LLM_JUDGE_FINAL.json
# - FAILED clusters: parse *_LLM_JUDGE_FAILED_DEBUG.json (if present),
#   and fill fields from last head/chief if possible.
# ============================================================

rm(list=ls())

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(optparse)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
suppressMessages(library(Triage))
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive=TRUE)

read_json_safe <- function(p) {
  tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error=function(e) NULL)
}

# Always return character(1)
as_chr1 <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.atomic(x) && length(x)==1 && !is.na(x)) return(as.character(x))
  tryCatch(
    jsonlite::toJSON(x, auto_unbox=TRUE, null="null"),
    error=function(e) {
      y <- as.character(x)[1]
      if (is.na(y) || !nzchar(y)) NA_character_ else y
    }
  )
}

# Turn vector/list to "a | b | c"
as_pipe_txt <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.character(x)) {
    x <- x[!is.na(x) & nzchar(x)]
    return(if (length(x)==0) NA_character_ else paste(x, collapse=" | "))
  }
  if (is.atomic(x)) {
    y <- as.character(x)
    y <- y[!is.na(y) & nzchar(y)]
    return(if (length(y)==0) NA_character_ else paste(y, collapse=" | "))
  }
  if (is.list(x)) {
    y <- unlist(x, recursive=TRUE, use.names=FALSE)
    y <- as.character(y)
    y <- y[!is.na(y) & nzchar(y)]
    return(if (length(y)==0) NA_character_ else paste(y, collapse=" | "))
  }
  as_chr1(x)
}

# Extract a safe scalar numeric
as_num1 <- function(x) {
  if (is.null(x)) return(NA_real_)
  suppressWarnings(as.numeric(x)[1])
}

# ============================================================
# CLI
# ============================================================
option_list <- list(
  make_option("--out_root", type="character", default="",
              help="Root output dir containing final/ (default: llm_judge_outputs/<dataset_name>)"),
  make_option("--out_dir", type="character", default="judge_post_summary",
              help="Subdir under out_root to write summary.csv"),
  make_option("--dataset_name", type="character", default="",
              help="Dataset name for default paths (default: basename(getwd()))."),
  make_option("--include_failed", action="store_true", default=TRUE,
              help="Include FAILED_DEBUG clusters in summary.csv [default TRUE]")
)
opt <- parse_args(OptionParser(option_list=option_list))

dataset_name <- if (nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- Triage:::get_dataset_config(dataset_name, project_root)
root <- opt$out_root
if (!nzchar(root)) root <- cfg$judge_outputs_root
final_dir <- file.path(root, "final")
if (!dir.exists(final_dir)) stop("Missing final dir: ", final_dir)

out_dir <- file.path(root, opt$out_dir)
ensure_dir(out_dir)

# ============================================================
# Discover clusters from final/
# ============================================================
final_ok_files   <- list.files(final_dir, pattern="_LLM_JUDGE_FINAL\\.json$", full.names=TRUE)
final_fail_files <- list.files(final_dir, pattern="_LLM_JUDGE_FAILED_DEBUG\\.json$", full.names=TRUE)

if (length(final_ok_files)==0 && length(final_fail_files)==0) {
  stop("No judge final/failed files found in: ", final_dir)
}

cluster_from_name <- function(p) {
  bn <- basename(p)
  sub("_LLM_JUDGE_.*$", "", bn)
}

ok_cids <- unique(vapply(final_ok_files, cluster_from_name, character(1)))
fail_cids <- unique(vapply(final_fail_files, cluster_from_name, character(1)))

# build union list (ok + failed)
clusters <- sort(unique(c(ok_cids, fail_cids)))

# ============================================================
# Parse one cluster
# ============================================================
parse_final_obj <- function(cid) {
  ok_path   <- file.path(final_dir, paste0(cid, "_LLM_JUDGE_FINAL.json"))
  fail_path <- file.path(final_dir, paste0(cid, "_LLM_JUDGE_FAILED_DEBUG.json"))
  
  status <- NA_character_
  obj <- NULL
  
  if (file.exists(ok_path)) {
    status <- "PASSED"
    obj <- read_json_safe(ok_path)
  } else if (isTRUE(opt$include_failed) && file.exists(fail_path)) {
    status <- "FAILED"
    obj <- read_json_safe(fail_path)
  } else {
    return(NULL)
  }
  
  # If FAILED debug bundle, prefer its head object if present
  if (status=="FAILED" && is.list(obj)) {
    if (is.list(obj$head)) obj <- obj$head
  }
  
  fd <- obj$final_decision %||% list()
  mv <- obj$method_verdict %||% list()
  pi <- obj$post_issues %||% list()
  ev <- obj$evidence %||% list()
  
  data.frame(
    cluster_id = as_chr1(obj$cluster_id %||% cid),
    judge_status = status,
    
    # final decision
    final_cell_type = as_chr1(fd$primary_cell_type %||% fd$final_cell_type %||% NULL),
    final_cell_ontology_id = as_chr1(fd$final_cell_ontology_id %||% NULL),
    decision_category = as_chr1(fd$decision_category %||% NULL),
    confidence = as_num1(fd$confidence %||% NA_real_),
    
    # Cassia verdict
    cassia_predicted_cell_type = as_chr1(mv$cassia$predicted_cell_type %||% NULL),
    cassia_cell_ontology_id = as_chr1(mv$cassia$cell_ontology_id %||% NULL),
    cassia_is_correct = as_chr1(mv$cassia$is_correct %||% NULL),
    cassia_error_type = as_chr1(mv$cassia$error_type %||% NULL),
    
    # Our verdict
    our_predicted_cell_type = as_chr1(mv$our_method$predicted_cell_type %||% NULL),
    our_cell_ontology_id = as_chr1(mv$our_method$cell_ontology_id %||% NULL),
    our_is_correct = as_chr1(mv$our_method$is_correct %||% NULL),
    our_error_type = as_chr1(mv$our_method$error_type %||% NULL),
    
    # post issues
    needs_manual_review = as_chr1(pi$needs_manual_review %||% NULL),
    flags = as_pipe_txt(pi$flags %||% NULL),
    notes = as_chr1(pi$notes %||% NULL),
    
    # evidence (optional but helpful for summary browsing)
    supporting_evidence = as_pipe_txt(ev$supporting %||% NULL),
    conflicting_evidence = as_pipe_txt(ev$conflicting %||% NULL),
    
    stringsAsFactors = FALSE
  )
}

rows <- purrr::map(clusters, parse_final_obj)
rows <- rows[!vapply(rows, is.null, logical(1))]
summary_df <- dplyr::bind_rows(rows)

# Order columns similar to “post_summary/summary.csv”
summary_df <- summary_df %>%
  arrange(cluster_id)

out_path <- file.path(out_dir, "summary.csv")
readr::write_csv(summary_df, out_path)

cat("[OK] wrote judge post summary: ", out_path, "\n", sep="")
cat(sprintf(
  "Rows: %d (PASSED=%d, FAILED=%d)\n",
  nrow(summary_df),
  sum(summary_df$judge_status == "PASSED", na.rm = TRUE),
  sum(summary_df$judge_status == "FAILED", na.rm = TRUE)
))
