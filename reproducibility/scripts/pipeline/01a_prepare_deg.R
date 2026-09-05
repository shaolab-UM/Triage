#!/usr/bin/env Rscript
# =============================================================
# 01a_prepare_deg.R — generate a standardized deg.csv（one task per stage）
#
# Input：
# Mode top100: findmarkersTop100.txt（cluster columns = 
#   Mode seurat: Seurat object in RDS format
# Output：
# deg.csv — columns: p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster,gene
#
# Usage：
#   Rscript 01a_prepare_deg.R --mode top100 --input <file> --out_dir <dir> [--n_clusters N]
#   Rscript 01a_prepare_deg.R --mode seurat --input <rds> --out_dir <dir> [--n_clusters N]
# =============================================================

rm(list = ls())
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--mode", type = "character", default = "top100",
              help = "seurat | top100 [default %default]"),
  make_option("--input", type = "character", default = "",
              help = "Input file path"),
  make_option("--out_dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--n_clusters", type = "integer", default = NULL,
              help = "Optionally retain the first N clusters（for validation；default is all clusters）"),
  make_option("--resolution", type = "double", default = 0.05,
              help = "Clustering resolution（mode=seurat）[default %default]"),
  make_option("--n_workers", type = "integer", default = 30,
              help = "FindAllMarkers Number of workers（mode=seurat）[default %default]"),
  make_option("--seed", type = "integer", default = 20251230,
              help = "Random seed [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opt$input) || !file.exists(opt$input)) stop("--input is required and must exist: ", opt$input)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(dplyr))

# ---- Mode B: top100 use directly ----
if (identical(opt$mode, "top100")) {
  deg <- readr::read_delim(opt$input, delim = if (grepl("\t", readLines(opt$input, 1))) "\t" else ",",
                           show_col_types = FALSE)
  std <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  miss <- setdiff(std, names(deg))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))
  deg <- deg[, std]
  if (!is.null(opt$n_clusters)) {
    keep <- head(sort(unique(as.character(deg$cluster))), opt$n_clusters)
    deg <- deg[as.character(deg$cluster) %in% keep, ]
    cat("  Subset: retaining", length(keep), "clusters\n")
  }
  readr::write_csv(deg, file.path(opt$out_dir, "deg.csv"))
  cat("[OK] deg.csv -> ", file.path(opt$out_dir, "deg.csv"),
      "| clusters:", length(unique(deg$cluster)), "| rows:", nrow(deg), "\n")
  quit(save = "no", status = 0)
}

# ---- Mode A: Seurat ----
if (!identical(opt$mode, "seurat")) stop("--mode must be seurat or top100")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(future)
})
set.seed(opt$seed)
rna <- readRDS(opt$input)
wsnn_graphs <- grep("^wsnn", Seurat::Graphs(rna), value = TRUE)
if (length(wsnn_graphs) == 0) stop("No wsnn graph was found in the Seurat object")
rna <- FindClusters(rna, graph.name = wsnn_graphs[1], algorithm = 3,
                    resolution = opt$resolution, random.seed = opt$seed)

DefaultAssay(rna) <- "SCT"
Idents(rna) <- rna$seurat_clusters
plan("multisession", workers = opt$n_workers)
all_markers <- FindAllMarkers(rna, assay = "SCT", slot = "data",
                              test.use = "wilcox", min.pct = 0.25,
                              logfc.threshold = 0.25, only.pos = TRUE)
plan("sequential")

deg <- all_markers %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::select(dplyr::any_of(c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")))

if (!is.null(opt$n_clusters)) {
  keep <- head(sort(unique(deg$cluster)), opt$n_clusters)
  deg <- deg[deg$cluster %in% keep, ]
}
readr::write_csv(deg, file.path(opt$out_dir, "deg.csv"))
cat("[OK] deg.csv -> ", file.path(opt$out_dir, "deg.csv"),
    "| clusters:", length(unique(deg$cluster)), "| rows:", nrow(deg), "\n")
