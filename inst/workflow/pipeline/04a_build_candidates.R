#!/usr/bin/env Rscript
# =============================================================
# 04a_build_candidates.R — candidates（OUR Core
#
# Input：filtered_deg.csv（Stage-03a output）
# Output：cellanno/structured_candidates.csv
#       cellanno/detailed_annotations.xlsx
#
# Usage：
#   Rscript 04a_build_candidates.R --deg <filtered_deg.csv> --out_dir <dir> \
#     --ontology <CL-ontology.json> [--species human] [--n_degs 40] [--top_n 30]
# =============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})
option_list <- list(
  make_option("--deg", type = "character", default = "filtered_deg.csv", help = "Stage-03a output"),
  make_option("--out_dir", type = "character", default = ".", help = "Output directory"),
  make_option("--ontology", type = "character",
              default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"),
              help = "CL ontology JSON"),
  make_option("--species", type = "character", default = "human", help = "human|mouse"),
  make_option("--n_degs", type = "integer", default = 40, help = "Top N DEGs per cluster"),
  make_option("--top_n", type = "integer", default = 30, help = "Number of candidates to retain")
)
opt <- parse_args(OptionParser(option_list = option_list))
suppressMessages(library(Triage)); suppressMessages(library(Triage))

if (!file.exists(opt$deg)) stop("DEG file does not exist: ", opt$deg)
if (!file.exists(opt$ontology)) stop("Ontology file does not exist: ", opt$ontology)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

species_label <- if (tolower(opt$species) == "mouse") "Mouse" else "Human"

# with caching
cache_dir <- file.path(opt$out_dir, ".kb_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
ontology <- Triage:::load_cl_ontology(opt$ontology, file.path(cache_dir, "cl_ontology.rds"))
kb <- Triage:::load_accordion_kb(species_label, file.path(cache_dir, "accordion_kb.rds"))

deg <- Triage:::read_deg(opt$deg)
# 03a 
cat("DEG clusters:", length(unique(deg$cluster)), "| rows:", nrow(deg), "\n")

query_degs_list <- deg %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(genes = list(gene), .groups = "drop") %>%
  { stats::setNames(.$genes, .$cluster) }

all_results <- list()
for (query_id in names(query_degs_list)) {
  degs_to_use <- head(query_degs_list[[query_id]], opt$n_degs)
  rankings_v4 <- Triage:::annotate_specificity_weighted_V4(degs_to_use, kb = kb, penalty_factor = 2)
  if (nrow(rankings_v4) == 0 || max(rankings_v4$final_score) < 0.1) {
    final_res <- Triage:::annotate_simple_ratio_V1(degs_to_use, kb = kb, penalty_factor = 2)
    rescue_info <- "Rescued_by_V1"
  } else {
    final_res <- rankings_v4
    rescue_info <- "V4_Success"
  }
  if (nrow(final_res) > 0) {
    final_res$Query_ID <- query_id
    final_res$rescue_status <- rescue_info
    all_results[[length(all_results) + 1]] <- final_res
  }
}

combined <- dplyr::bind_rows(all_results)
if (nrow(combined) == 0) stop("No candidate results were generated")

out_cellanno <- file.path(opt$out_dir, "cellanno")
dir.create(out_cellanno, showWarnings = FALSE, recursive = TRUE)

curated <- combined %>%
  dplyr::filter(rank <= opt$top_n) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(top_evidence_genes = paste(head(unlist(evidence_genes), 15), collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::select(Query_ID, cell_type, final_score, rank, top_evidence_genes, rescue_status)

structured_path <- file.path(out_cellanno, "structured_candidates.csv")
data.table::fwrite(curated, structured_path)
cat("[OK] structured_candidates.csv -> ", structured_path, "\n")

if (requireNamespace("writexl", quietly = TRUE)) {
  detailed <- combined %>%
    dplyr::rowwise() %>%
    dplyr::mutate(evidence_genes = paste(unlist(evidence_genes), collapse = ", ")) %>%
    dplyr::ungroup()
  writexl::write_xlsx(detailed, file.path(out_cellanno, "detailed_annotations.xlsx"))
  cat("[OK] detailed_annotations.xlsx\n")
}
cat("[DONE] 04a_build_candidates | clusters:", length(unique(combined$Query_ID)), "\n")
