#!/usr/bin/env Rscript
# =============================================================
# 02b_build_maskdeg.R — anonymize DEG cluster identifiers
#
# Input：deg.csv（Stage-01a output）、cluster_map.csv（Stage-02a output）
# Output：maskdeg.csv（cluster columns = cluster_N 
#
# Usage：
#   Rscript 02b_build_maskdeg.R --deg <deg.csv> --map <cluster_map.csv> --out_dir <dir>
# =============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
})
option_list <- list(
  make_option("--deg", type = "character", default = "deg.csv", help = "Stage-01a output"),
  make_option("--map", type = "character", default = "cluster_map.csv", help = "Stage-02a output"),
  make_option("--out_dir", type = "character", default = ".", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
suppressMessages(library(Triage))

if (!file.exists(opt$deg)) stop("deg file does not exist: ", opt$deg)
if (!file.exists(opt$map)) stop("cluster_map file does not exist: ", opt$map)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

deg <- read_deg(opt$deg)
cluster_map <- readr::read_csv(opt$map, show_col_types = FALSE) %>%
  dplyr::mutate(cluster = stringr::str_trim(as.character(cluster)))

mask_deg <- deg %>%
  dplyr::mutate(cluster = stringr::str_trim(as.character(cluster))) %>%
  dplyr::left_join(cluster_map, by = "cluster") %>%
  dplyr::mutate(cluster = cluster_anon)

if (any(is.na(mask_deg$cluster))) {
  stop("cluster_map does not cover all clusters; missing: ",
       paste(unique(mask_deg$cluster[is.na(mask_deg$cluster)]), collapse = ", "))
}

mask_deg <- mask_deg %>%
  dplyr::select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene)

readr::write_csv(mask_deg, file.path(opt$out_dir, "maskdeg.csv"))
cat("[OK] maskdeg.csv -> ", file.path(opt$out_dir, "maskdeg.csv"),
    "| clusters:", length(unique(mask_deg$cluster)), "| rows:", nrow(mask_deg), "\n")
