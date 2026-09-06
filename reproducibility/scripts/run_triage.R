#!/usr/bin/env Rscript
# =========================================================================
# run_triage.R — single-user-entry orchestrator for the Triage workflow.
#
# Runs the released pipeline stages (no scientific reimplementation):
#   user DEG -> deterministic anonymization -> 03a -> 03b (CASSIA)
#   -> 04a candidates -> 05 evidence -> 06 in-house reviewer
#   -> 06b enrichment reviewer -> 07/07b summaries -> 07.5 CL mapping
#   -> 08 adjudication inputs -> 09 adjudication (Handling Editor / Chief QC)
#   -> 10 final summary [-> 11 optional evaluation when --reference-labels]
#
# Generic mode:
#   Rscript reproducibility/scripts/run_triage.R \
#     --deg path/to/markers.csv --species human --tissue pancreas --out results/
#   optional: --study-context normal_adult --dataset-name my_dataset
#             --reference-labels reference.csv --workers 4
#             --cluster-id cluster_3 --n-clusters 1   (single-cluster preview)
#             --preflight-only
#
# Benchmark mode (bundled release data):
#   Rscript reproducibility/scripts/run_triage.R --benchmark Census_immune --out results/
#
# Reference labels are NEVER used for adjudication; they are only read when
# --reference-labels is explicitly supplied (post-adjudication evaluation).
# =========================================================================

suppressMessages({
  library(optparse)
})

option_list <- list(
  make_option("--deg", type = "character", default = NULL,
              help = "Cluster-level DEG/marker table (csv/tsv; required columns: cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2). Generic mode."),
  make_option("--species", type = "character", default = "human",
              help = "human or mouse (generic mode; default human)."),
  make_option("--tissue", type = "character", default = NULL,
              help = "Tissue context for the CASSIA stage (generic mode)."),
  make_option("--study-context", type = "character", default = "reference",
              help = "Optional study context for the CASSIA stage (default 'reference')."),
  make_option("--dataset-name", type = "character", default = "user_dataset",
              help = "Name for this run (default 'user_dataset')."),
  make_option("--out", type = "character", default = "results",
              help = "Output root directory (default 'results')."),
  make_option("--benchmark", type = "character", default = NULL,
              help = "Run a bundled benchmark dataset (e.g. Census_immune). Uses the released masked DEG input and dataset configuration."),
  make_option("--reference-labels", type = "character", default = NULL,
              help = "Optional reference-label CSV for post-adjudication evaluation (stage 11). Never used in adjudication."),
  make_option("--workers", type = "numeric", default = 4,
              help = "Worker count for parallel stages (default 4)."),
  make_option("--run-tag", type = "character", default = NULL,
              help = "Run tag; defaults to a timestamp."),
  make_option("--cluster-id", type = "character", default = NULL,
              help = "Restrict a generic run to one anonymous cluster id (for example cluster_1; single-cluster preview, applied after anonymization)."),
  make_option("--n-clusters", type = "numeric", default = NULL,
              help = "Restrict a generic run to the first N anonymous clusters after deterministic anonymization (single-cluster preview)."),
  make_option("--preflight-only", action = "store_true", default = FALSE,
              help = "Run the full-workflow preflight check and exit without analysis.")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               usage = "run_triage.R [options]"))

script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
triage_home <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")
pipeline_dir <- file.path(triage_home, "reproducibility", "scripts", "pipeline")

Sys.setenv(TRIAGE_HOME = triage_home)
Sys.setenv(PROJECT_ROOT = triage_home)
Sys.setenv(CL_LOCAL_JSON = Sys.getenv("CL_LOCAL_JSON",
  unset = file.path(triage_home, "inputs", "raw", "ontology",
                    "CL-ontology-v2025-07-30.json")))

# ------------------------------ helpers ---------------------------------
fail <- function(msg) {
  message("run_triage: ERROR: ", msg)
  quit(status = 1)
}
run_stage <- function(script, args, stage_label) {
  message(">>> [", stage_label, "] Rscript ", script)
  status <- system2("Rscript",
                    c(file.path(pipeline_dir, script), args),
                    stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) {
    fail(paste0("stage ", stage_label, " (", script, ") failed with status ", status))
  }
}

# Benchmark context resolution (bundled release data only)
benchmark_maskdeg <- function(ds) {
  p <- file.path(triage_home, "reproducibility", "primary", ds,
                 "model_inputs", "maskdeg.csv")
  if (!file.exists(p)) {
    fail(paste0("benchmark dataset '", ds, "' not found (missing ", p, ")"))
  }
  p
}

suppressMessages(library(Triage))
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

# Resolve dataset context (species/tissue/study_context) for benchmark mode;
# generic mode values come from CLI args.
if (!is.null(opt$benchmark)) {
  cfg <- Triage:::get_dataset_config(opt$benchmark, project_root = triage_home)
  species <- tolower(as.character(cfg$species %||% "human"))
  tissue <- as.character(cfg$tissue[[1]] %||% "")
  study_context <- as.character(cfg$study_context[[1]] %||% "reference")
  dataset_name <- opt$benchmark
  deg_source <- benchmark_maskdeg(opt$benchmark)
  if (!is.null(opt$species) && !identical(opt$species, "human")) {
    message("run_triage: --species ignored in benchmark mode (dataset config species = ", species, ").")
  }
} else {
  if (is.null(opt$deg)) fail("--deg is required in generic mode (or use --benchmark).")
  if (is.null(opt$tissue)) fail("--tissue is required in generic mode.")
  if (!opt$species %in% c("human", "mouse")) {
    fail("--species must be 'human' or 'mouse'.")
  }
  species <- opt$species
  tissue <- opt$tissue
  study_context <- opt$`study-context`
  dataset_name <- opt$`dataset-name`
  deg_source <- opt$deg
}

run_tag <- opt$`run-tag` %||% format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(opt$out, dataset_name, run_tag)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
message("run_triage: dataset = ", dataset_name,
        " | species = ", species, " | tissue = ", tissue,
        " | study context = ", study_context)
message("run_triage: DEG source = ", deg_source)
message("run_triage: output root = ", out_root)

# ------------------------------ preflight -------------------------------
preflight_args <- c("--species", species)
if (!is.null(opt$deg)) preflight_args <- c(preflight_args, "--deg", opt$deg)
message(">>> [preflight] full-workflow dependency check")
status <- system2("Rscript", c(file.path(script_dir, "preflight_check.R"),
                               preflight_args), stdout = "", stderr = "")
if (!identical(as.integer(status), 0L)) {
  fail("preflight check failed: resolve the reported prerequisites and rerun.")
}
if (opt$`preflight-only`) {
  message("run_triage: --preflight-only requested; stopping before analysis.")
  quit(status = 0)
}

# ------------------------- input preparation ----------------------------
map_dir <- file.path(out_root, "01_input_prep")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
if (!is.null(opt$benchmark)) {
  # Bundled benchmark: the released masked DEG is already anonymized; reuse it
  # directly (no true labels involved).
  benchmark_deg <- readr::read_csv(deg_source, show_col_types = FALSE)
  benchmark_cluster_col <- if ("cluster_anon" %in% names(benchmark_deg)) {
    "cluster_anon"
  } else if ("cluster" %in% names(benchmark_deg)) {
    "cluster"
  } else {
    fail(paste0("benchmark masked DEG lacks an anonymous cluster column: ", deg_source))
  }
  benchmark_ids <- unique(as.character(benchmark_deg[[benchmark_cluster_col]]))
  cluster_num <- suppressWarnings(as.integer(sub("^cluster_", "", benchmark_ids)))
  benchmark_ids <- if (all(!is.na(cluster_num))) {
    benchmark_ids[order(cluster_num)]
  } else {
    sort(benchmark_ids)
  }
  if (!is.null(opt$`cluster-id`)) {
    if (!opt$`cluster-id` %in% benchmark_ids) {
      fail(paste0("--cluster-id must name an anonymous benchmark ID; available IDs include ",
                  paste(utils::head(benchmark_ids, 5), collapse = ", ")))
    }
    benchmark_ids <- opt$`cluster-id`
  }
  if (!is.null(opt$`n-clusters`)) {
    n <- as.integer(opt$`n-clusters`)
    if (is.na(n) || n < 1L) fail("--n-clusters must be a positive integer.")
    benchmark_ids <- utils::head(benchmark_ids, n)
  }
  benchmark_deg <- dplyr::filter(benchmark_deg,
                                 .data[[benchmark_cluster_col]] %in% benchmark_ids)
  if (nrow(benchmark_deg) == 0L) fail("benchmark cluster filtering produced no DEG rows.")
  readr::write_csv(benchmark_deg, file.path(map_dir, "maskdeg.csv"))
  maskdeg <- file.path(map_dir, "maskdeg.csv")
  message("run_triage: benchmark mode uses ", length(benchmark_ids),
          " anonymous cluster(s) from released masked DEG: ", maskdeg)
} else {
  prep_args <- c("--deg", opt$deg, "--out-dir", map_dir)
  if (!is.null(opt$`cluster-id`)) {
    prep_args <- c(prep_args, "--cluster-id", opt$`cluster-id`)
  }
  if (!is.null(opt$`n-clusters`)) {
    # deterministic first-N restriction happens after mapping in 01b
    prep_args <- c(prep_args, "--n-clusters", as.character(opt$`n-clusters`))
  }
  run_stage("01b_prepare_user_deg.R", prep_args, "01b user-DEG preparation")
  maskdeg <- file.path(map_dir, "maskdeg.csv")
}
deg_for_stages <- maskdeg

# ------------------------------- stages ---------------------------------
run_dir <- out_root

# 03a: deterministic DEG filtering
deg_filtered <- file.path(run_dir, "filtered_deg.csv")
run_stage("03a_filter_deg.R", c("--deg", deg_for_stages,
                                "--out_dir", run_dir), "03a")

# 03b: CASSIA reviewer (LLM API; CASSIA strips the /chat/completions suffix)
cassia_out <- file.path(run_dir, "03b_cassia")
run_stage("03b_run_cassia.R", c(
  "--deg", deg_filtered, "--out_dir", cassia_out,
  "--tissue", tissue, "--species", species,
  "--study_context", study_context,
  "--workers", as.character(opt$workers),
  "--api_base_url", Sys.getenv("LLM_API_BASE_URL")
), "03b")

# 04a: structured candidates
candidates_out <- file.path(run_dir, "04a_candidates")
run_stage("04a_build_candidates.R", c(
  "--deg", deg_filtered, "--out_dir", candidates_out,
  "--species", species
), "04a")

candidates_csv <- file.path(candidates_out, "cellanno", "structured_candidates.csv")

# 05: evidence / dossier construction (writes relative to cwd = run_dir)
withr::with_dir(run_dir, {
  run_stage("05_prepare_llm_inputs.R", c(
    "--mode", "full",
    "--deg_file", deg_filtered,
    "--candidates_file", candidates_csv,
    "--out_root", file.path(run_dir, "05c_llm_queries"),
    "--dataset_name", dataset_name,
    "--project_root", triage_home
  ), "05")
})

# 06: in-house reviewer (LLM API)
run_stage("06_run_llm_pipeline.R", c(
  "--mode", "full",
  "--input_root", file.path(run_dir, "05c_llm_queries"),
  "--output_root", file.path(run_dir, "06_llm_outputs"),
  "--dataset_name", dataset_name,
  "--project_root", triage_home,
  "--workers", as.character(opt$workers),
  "--max_rounds", "3",
  "--cl_local_json", Sys.getenv("CL_LOCAL_JSON")
), "06")

# 06b: enrichment reviewer (clusterProfiler + CellMarkerDB)
run_stage("06b_run_inter.R", c(
  "--marker_csv", deg_filtered,
  "--step1_dir", file.path(run_dir, "05c_llm_queries", "step1_report_queries"),
  "--bioinfo_dir", file.path(run_dir, "intermediate_outputs",
                             paste0(dataset_name, "_LLM_Input_Run"),
                             "bioinformatics_tsv"),
  "--out_dir", file.path(run_dir, "06b_inter"),
  "--dataset_name", dataset_name,
  "--species", species,
  "--workers", as.character(opt$workers),
  "--cl_local_json", Sys.getenv("CL_LOCAL_JSON")
), "06b")

# 07 / 07b: reviewer summaries
run_stage("07_our_llm_summary.R", c(
  "--out_root", file.path(run_dir, "06_llm_outputs"),
  "--dataset_name", dataset_name,
  "--out_dir", file.path(run_dir, "07_our_summary")
), "07")

run_stage("07b_inter_summary.R", c(
  "--in_dir", file.path(run_dir, "06b_inter", "final_passed"),
  "--out_dir", file.path(run_dir, "07b_inter_summary"),
  "--dataset_name", dataset_name
), "07b")

# 07.5: CL-Linker reviewer mapping (LLM-backed)
manifest <- file.path(run_dir, "evidence_mapping_run.tsv")
writeLines(c("dataset\trun_dir",
             paste0(dataset_name, "\t", run_dir)),
           manifest, sep = "\n")
run_stage("07.5_build_reviewer_mapping.R", c(
  "--cl_json", Sys.getenv("CL_LOCAL_JSON"),
  "--out_dir", file.path(run_dir, "07.5_mapping"),
  "--manifest", manifest,
  "--dataset", dataset_name
), "07.5")

# 08: adjudication input assembly
registry_csv <- list.files(file.path(run_dir, "07.5_mapping"),
                           pattern = "reviewer_mapping_registry.*[.]tsv$",
                           full.names = TRUE)
if (length(registry_csv) != 1L) {
  fail("expected exactly one reviewer mapping registry under 07.5_mapping")
}
run_stage("08_build_judge_inputs.R", c(
  "--cassia_csv", file.path(cassia_out,
    list.files(cassia_out, pattern = "annotation_cassia_FINAL_RESULTS.csv$",
               recursive = TRUE)[1]),
  "--in_house_summary_csv", file.path(run_dir, "07_our_summary", "summary.csv"),
  "--enrichment_summary_csv", file.path(run_dir, "07b_inter_summary", "summary.csv"),
  "--mapping_registry_csv", registry_csv,
  "--intermediate_outputs_dir", file.path(run_dir, "intermediate_outputs",
                                          paste0(dataset_name, "_LLM_Input_Run")),
  "--step1_dir", file.path(run_dir, "05c_llm_queries", "step1_report_queries"),
  "--out_dir", file.path(run_dir, "08_judge_inputs"),
  "--dataset_name", dataset_name
), "08")

# 09: adjudication (Handling Editor; Chief QC / release-state in batch workflow)
run_stage("09_run_judge.R", c(
  "--judge_input_dir", file.path(run_dir, "08_judge_inputs"),
  "--out_root", file.path(run_dir, "09_judge_outputs"),
  "--dataset_name", dataset_name,
  "--temperature", "0",
  "--workers", as.character(opt$workers),
  "--max_rounds", "3",
  "--always_run_chief",
  "--release_policy", "auto",
  "--ols_first", "FALSE",
  "--cl_local_json", Sys.getenv("CL_LOCAL_JSON")
), "09")

# 10: final summary
run_stage("10_judge_post_summary.R", c(
  "--out_root", file.path(run_dir, "09_judge_outputs"),
  "--dataset_name", dataset_name
), "10")

# 11: optional evaluation (reference labels enter ONLY here)
if (!is.null(opt$`reference-labels`)) {
  run_stage("11_eval_accuracy.R", c(
    "--dataset_name", dataset_name,
    "--true_label_csv", opt$`reference-labels`,
    "--cassia_csv", file.path(cassia_out,
      list.files(cassia_out, pattern = "annotation_cassia_FINAL_RESULTS.csv$",
                 recursive = TRUE)[1]),
    "--our_csv", file.path(run_dir, "07_our_summary", "summary.csv"),
    "--inter_csv", file.path(run_dir, "07b_inter_summary", "summary.csv"),
    "--judge_csv", file.path(run_dir, "09_judge_outputs", "summary_final.csv"),
    "--judge_final_dir", file.path(run_dir, "09_judge_outputs", "final"),
    "--out_dir", file.path(run_dir, "11_eval"),
    "--cl_json", Sys.getenv("CL_LOCAL_JSON")
  ), "11 evaluation")
} else {
  message("run_triage: no --reference-labels supplied; evaluation (stage 11) skipped. ",
          "Reference labels are never used for adjudication.")
}

message("run_triage: complete. Outputs under ", out_root)
