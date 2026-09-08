#!/usr/bin/env Rscript
# 01b_prepare_user_deg.R — deterministic user-input preparation
#
# Takes a cluster-level DEG/marker table from the user, validates the required
# columns, and produces the anonymized inputs consumed by pipeline stage 03a:
#   <out_dir>/cluster_map.csv  (cluster_id,cluster_label)
#   <out_dir>/maskdeg.csv      (user DEG rows with the anonymous cluster id in the `cluster` column)
#
# No true labels are generated. Original cluster names are preserved ONLY in
# the local cluster_map.csv for traceability and are never interpreted
# biologically. The generated maskdeg.csv contains only anonymous IDs, so
# stages 03a-10 never receive the original cluster identifiers.

suppressMessages({
  library(optparse); library(dplyr); library(readr); library(stringr)
})

option_list <- list(
  make_option("--deg", type = "character", help = "User DEG/marker table (csv or tsv)"),
  make_option("--out-dir", type = "character", help = "Output directory"),
  make_option("--cluster-id", type = "character", default = NULL,
              help = "Optional: restrict preparation to this anonymous cluster id (for example cluster_1; single-cluster preview)"),
  make_option("--n-clusters", type = "numeric", default = NULL,
              help = "Optional: restrict preparation to the first N clusters after deterministic mapping (single-cluster preview)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$deg) || is.null(opt$`out-dir`)) {
  stop("01b_prepare_user_deg.R: both --deg and --out-dir are required.")
}

REQUIRED_COLS <- c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")

deg <- tryCatch(
  read_delim(opt$deg, delim = NULL, show_col_types = FALSE, guess_max = 100000),
  error = function(e) stop("01b_prepare_user_deg.R: cannot read DEG file '", opt$deg, "': ", conditionMessage(e))
)
missing <- setdiff(REQUIRED_COLS, names(deg))
if (length(missing) > 0) {
  stop("01b_prepare_user_deg.R: user DEG table is missing required column(s): ",
       paste(missing, collapse = ", "),
       "\nRequired columns: ", paste(REQUIRED_COLS, collapse = ", "))
}

deg <- as.data.frame(deg)
deg$cluster <- as.character(deg$cluster)
deg$gene <- as.character(deg$gene)

if (any(is.na(deg$cluster)) || any(!nzchar(deg$cluster))) {
  stop("01b_prepare_user_deg.R: 'cluster' column contains missing or empty values.")
}

# Deterministic cluster mapping: numeric sort when all ids are numeric,
# lexicographic otherwise. cluster_1..N; original label kept only in the map.
orig <- unique(deg$cluster)
numeric_all <- !any(is.na(suppressWarnings(as.numeric(orig))))
orig <- if (numeric_all) {
  orig[order(as.numeric(orig))]
} else {
  sort(orig)
}
map <- tibble::tibble(cluster_label = orig,
                      cluster_id = paste0("cluster_", seq_along(orig)))

maskdeg <- deg %>%
  inner_join(map, by = c("cluster" = "cluster_label")) %>%
  select(cluster = cluster_id, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2)

if (!is.null(opt$`cluster-id`)) {
  keep_id <- as.character(opt$`cluster-id`)
  if (!keep_id %in% map$cluster_id) {
    stop("01b_prepare_user_deg.R: --cluster-id '", keep_id,
         "' is not an anonymous cluster id generated from this DEG table.")
  }
  map <- map %>% filter(cluster_id == keep_id)
  maskdeg <- maskdeg %>% filter(cluster == keep_id)
}

if (!is.null(opt$`n-clusters`)) {
  n <- as.integer(opt$`n-clusters`)
  if (is.na(n) || n < 1L) stop("01b_prepare_user_deg.R: --n-clusters must be a positive integer.")
  keep <- paste0("cluster_", seq_len(n))
  map <- map %>% filter(cluster_id %in% keep)
  maskdeg <- maskdeg %>% filter(cluster %in% keep)
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
write_csv(map, file.path(opt$`out-dir`, "cluster_map.csv"))
write_csv(maskdeg, file.path(opt$`out-dir`, "maskdeg.csv"))

cat(sprintf("01b_prepare_user_deg: %d clusters -> cluster_1..cluster_%d; %d DEG rows\n",
            nrow(map), nrow(map), nrow(maskdeg)))
cat(sprintf("01b_prepare_user_deg: wrote cluster_map.csv and maskdeg.csv to %s\n",
            opt$`out-dir`))
