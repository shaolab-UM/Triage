#!/usr/bin/env Rscript
# =============================================================
# 03b_run_cassia.R — CASSIA annotation（LLM API）
#
# Input：filtered_deg.csv（Stage-03a output）
# Output：CASSIA_<tissue>_<species>_<timestamp>/
#         01_annotation_results/annotation_cassia_FINAL_RESULTS.csv
#
# Usage：
#   Rscript 03b_run_cassia.R --deg <filtered_deg.csv> --out_dir <dir> \
#     --tissue blood --species human [--workers N] [--model ...]
# =============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(optparse)
  library(CASSIA)
  library(dplyr)
  library(stringr)
})
option_list <- list(
  make_option("--deg", type = "character", default = "filtered_deg.csv", help = "Stage-03a output"),
  make_option("--out_dir", type = "character", default = ".", help = "Output directory"),
  make_option("--tissue", type = "character", default = "blood", help = "Tissue"),
  make_option("--species", type = "character", default = "human", help = "Species"),
  make_option("--study_context", type = "character", default = "reference", help = "Study context"),
  make_option("--workers", type = "integer", default = 0, help = "Number of workers，0=automatic(CPU*0.2)"),
  make_option("--score_threshold", type = "integer", default = 75, help = "Minimum score"),
  make_option("--api_key_env", type = "character", default = "DEEPSEEK_API_KEY", help = "Environment-variable name containing the API key"),
  make_option("--api_base_url", type = "character",
              default = Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = "XXXXX")),
              help = "API base URL; set LLM_API_BASE_URL or CASSIA_API_BASE_URL"),
  make_option("--model", type = "character", default = "deepseek-v4-flash", help = "Model identifier")
)
opt <- parse_args(OptionParser(option_list = option_list))
suppressMessages(library(Triage))

if (!file.exists(opt$deg)) stop("DEG file does not exist: ", opt$deg)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# API（read from an environment variable; no hard-coded credential）
api_key <- Sys.getenv(opt$api_key_env, unset = "")
if (!nzchar(api_key)) stop("Missing API key; export ", opt$api_key_env)
api_provider <- sub("/chat/completions$", "", opt$api_base_url)
api_provider <- sub("/v1/chat/completions$", "", api_provider)
setLLMApiKey(api_key, provider = api_provider, persist = FALSE)

workers <- if (opt$workers > 0) opt$workers else max(1L, floor(parallel::detectCores() * 0.2))

deg <- Triage:::read_deg(opt$deg)
cat("DEG clusters:", length(unique(deg$cluster)), "| rows:", nrow(deg), "\n")

additional_info <- paste(
  "species:", opt$species,
  "tissue:", opt$tissue,
  "study_context:", opt$study_context,
  "Goal: marker-based candidate cell type annotation at cluster/state level.",
  sep = "\n"
)

old_wd <- getwd()
setwd(opt$out_dir)
on.exit(setwd(old_wd), add = TRUE)

runCASSIA_pipeline(
  output_file_name = "annotation_cassia",
  tissue = opt$tissue,
  species = opt$species,
  marker = deg,
  max_workers = workers,
  annotation_model = opt$model,
  annotation_provider = api_provider,
  score_model = opt$model,
  score_provider = api_provider,
  annotationboost_model = opt$model,
  annotationboost_provider = api_provider,
  merge_model = opt$model,
  merge_provider = api_provider,
  score_threshold = opt$score_threshold,
  additional_info = additional_info
)
cat("[DONE] 03b_run_cassia | Output: ", opt$out_dir, "\n")
