#!/usr/bin/env Rscript
# =============================================================
# 02a_build_true_label.R — extract reference labels and create anonymous cluster mapping
#
# Input：deg.csv（Stage-01a output；cluster columns = 
# Output：true_label.csv（cluster_id + cluster_label）
# cluster_map.csv（ <-> cluster_N 
#
# Usage：
#   Rscript 02a_build_true_label.R --input <deg.csv> --out_dir <dir>
# =============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

option_list <- list(
  make_option("--input", type = "character", default = "deg.csv", help = "Stage-01a output deg.csv"),
  make_option("--out_dir", type = "character", default = ".", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
source(file.path(dirname(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/", mustWork = FALSE)))), "lib", "00_utils.R"))

deg <- read_deg(opt$input)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Stable anonymization：sort clusters before assigning anonymous identifiers。
# numeric cluster identifiers ( "0","1","10") sort numerically; otherwise sort alphabetically
cluster_map <- tibble::tibble(cluster = stringr::str_trim(as.character(deg$cluster))) %>%
  dplyr::filter(nzchar(cluster)) %>%
  dplyr::distinct(cluster)

# ("0","1") "cluster_N" → sort numerically;
cluster_numeric <- suppressWarnings(as.numeric(cluster_map$cluster))
cluster_n_from_prefix <- suppressWarnings(as.numeric(sub("^cluster_", "", cluster_map$cluster)))
if (all(!is.na(cluster_numeric))) {
  cluster_map <- cluster_map %>% dplyr::arrange(cluster_numeric)
} else if (all(!is.na(cluster_n_from_prefix))) {
  cluster_map <- cluster_map %>% dplyr::arrange(cluster_n_from_prefix)
} else {
  cluster_map <- cluster_map %>% dplyr::arrange(cluster)
}
cluster_map <- cluster_map %>%
  dplyr::mutate(cluster_anon = paste0("cluster_", dplyr::row_number()))

true_label <- cluster_map %>%
  dplyr::mutate(
    cluster_label = standardize_true_label(cluster),
    cluster_id = sub("^cluster_", "", cluster_anon)
  ) %>%
  dplyr::select(cluster_id, cluster_label) %>%
  dplyr::arrange(as.integer(cluster_id))

readr::write_csv(cluster_map, file.path(opt$out_dir, "cluster_map.csv"))
readr::write_csv(true_label, file.path(opt$out_dir, "true_label.csv"))
cat("[OK] true_label.csv (", nrow(true_label), " clusters)\n")
cat("     cluster_map.csv (anonymous mapping)\n")
print(true_label)
