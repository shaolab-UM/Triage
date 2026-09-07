#!/usr/bin/env Rscript
# =============================================================
# 03a_filter_deg.R — Noise-gene filtering（03/04 
#
# Input：maskdeg.csv（02b Output）
# Output：filtered_deg.csv（/Mitochondrial genes//
#
# Usage：
#   Rscript 03a_filter_deg.R --deg <maskdeg.csv> --out_dir <dir>
# =============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})
option_list <- list(
  make_option("--deg", type = "character", default = "maskdeg.csv", help = "02b Output"),
  make_option("--out_dir", type = "character", default = ".", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
suppressMessages(library(Triage))

if (!file.exists(opt$deg)) stop("deg file does not exist: ", opt$deg)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

deg <- Triage:::read_deg(opt$deg)
filtered <- Triage:::filter_deg(deg)

readr::write_csv(filtered, file.path(opt$out_dir, "filtered_deg.csv"))
cat("[OK] filtered_deg.csv -> ", file.path(opt$out_dir, "filtered_deg.csv"),
    "| rows:", nrow(filtered), "（original", nrow(deg), "）",
    "| clusters:", length(unique(filtered$cluster)), "\n")
