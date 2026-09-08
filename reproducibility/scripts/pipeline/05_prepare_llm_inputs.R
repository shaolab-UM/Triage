#!/usr/bin/env Rscript
# ============================================================
# 05_prepare_llm_inputs.R — LLM Input (step1 report queries)
# Input: 03a filtered_deg.csv + 04a candidates
# Output: outputs/<run>/05c_llm_queries/step1_report_queries/<cluster>_step1_report_query.json
# ============================================================

rm(list = ls())

# --- 1. Environment setup and package loading ---
suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(glue); library(stringr); library(readr)
  library(dplyr); library(purrr); library(memoise); library(cachem)
  library(clusterProfiler); library(org.Hs.eg.db); library(DOSE); library(ReactomePA)
  library(enrichR); library(data.table); library(tidyr)
  library(knitr)
  library(future); library(future.apply); library(tictoc)
  library(decoupleR); library(tibble)
  library(readxl)    # retained as part of the released workflow structure
  library(optparse)  # retained as part of the released workflow structure
  library(rlang)
})

# source configuration and llm_run（locate from TRIAGE_HOME or the script location）
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
suppressMessages(library(Triage))

# --- 2. Load custom functions ---
# Required build_validator_instructions() / generate_expert_report_query() /
# generate_citation_fix_query() are implemented in llm_run.R; no separate llm_simple.R is required.
# TRIAGE_WORKFLOW_DIR lets the installed-package runtime (inst/workflow)
# point at its bundled copy of llm_run.R; repository runs use the default.
llm_run_dir <- Sys.getenv("TRIAGE_WORKFLOW_DIR",
                          unset = file.path(triageHome, "reproducibility", "scripts"))
source(file.path(llm_run_dir, "llm_run.R"))

# --- 2b. CLI args (adds dryrun mode) ---
option_list <- list(
  make_option(c("--mode"), type="character", default="full",
              help="Run mode: full | dryrun [default %default]"),
  make_option(c("--cluster_id"), type="character", default=NULL,
              help="If set, only process this cluster id (Query_ID)."),
  make_option(c("--n_clusters"), type="integer", default=NULL,
              help="In dryrun, optionally process first N clusters (after filtering)."),
  make_option(c("--seed"), type="integer", default=20251230,
              help="Random seed for any sampling [default %default]"),
  make_option(c("--deg_file"), type="character", default=NULL,
              help="Path to DEG CSV (overrides deg_file_path)."),
  make_option(c("--candidates_file"), type="character", default=NULL,
              help="Path to candidates CSV (overrides candidates_file_path)."),
  make_option(c("--out_root"), type="character", default=NULL,
              help="Output root directory for runner inputs (default: llm_inputs/<dataset_name>)."),
  make_option(c("--dataset_name"), type="character", default=NULL,
              help="Dataset name for default paths (default: basename(getwd()))."),
  make_option(c("--project_root"), type="character", default=getwd(),
              help="Project root directory (used to resolve relative resource files like collectri_*.rds)."),
  make_option(c("--top_k_candidates"), type="integer", default=10,
              help="Keep only top K candidates per cluster in LLM dossier [default %default]"),
  make_option(c("--disable_smart_context"), action="store_true", default=FALSE,
              help="(Recommended) Do not inject smart_context into the prompt context."),
  make_option(c("--cache_only"), action="store_true", default=FALSE,
              help="Use cache only; skip new API/bioinfo calls if cache miss."),
  make_option(c("--intermediate_root"), type="character", default="",
              help="Absolute root for intermediate enrichment TSVs (wired through TRIAGE_INTERMEDIATE_ROOT). Default keeps the historical working-directory-relative location.")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Ensure resource files are resolvable in future workers
Sys.setenv(PROJECT_ROOT = opt$project_root)
try(setwd(opt$project_root), silent = TRUE)
# Separate resource/project root from the intermediate-output root: when an
# explicit --intermediate_root is supplied it is wired through
# TRIAGE_INTERMEDIATE_ROOT so enrichment TSVs land under the run directory
# regardless of working directories.
if (nzchar(opt$intermediate_root)) {
  dir.create(opt$intermediate_root, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(TRIAGE_INTERMEDIATE_ROOT = normalizePath(opt$intermediate_root, mustWork = TRUE))
}

set.seed(opt$seed)

# ===============================================================
# --- 3. Global configuration (configure parameters here) ---
# ===============================================================

# ##--- 3a. 【Core】Input ---
current_dataset_name <- if (!is.null(opt$dataset_name) && nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
cfg <- Triage:::get_dataset_config(current_dataset_name, opt$project_root)
# Release structure：does not switch into final/<dataset>（paths are supplied explicitly via CLI）
deg_file_path <- cfg$deg_file
candidates_file_path <- cfg$candidates_file

# --- CLI overrides ---
if (!is.null(opt$deg_file)) deg_file_path <- opt$deg_file
if (!is.null(opt$candidates_file)) candidates_file_path <- opt$candidates_file
out_root <- opt$out_root %||% cfg$llm_inputs_root
mode <- tolower(opt$mode)

# (4) 【Important】Specify the actual column names used by the DEG CSV file
deg_original_colnames <- list(
  cluster_col = "cluster",      # <-- column containing the cell-cluster name or ID
  gene_col    = "gene",         # <-- column containing the gene symbol
  logfc_col   = "avg_log2FC"    # <-- column containing the log fold-change value
)

# (5) 【Important】Specify the actual column names used by the candidate-cell CSV file
candidates_original_colnames <- list(
  cluster_col = "Query_ID"      # <-- column containing the cell-cluster name or ID (must correspond to the DEG cluster identifiers)
)

###--- 3b. 【Core】analysis-context configuration ---
prompt_context_config <- list(
  species = cfg$species,
  tissue = cfg$tissue,
  data_type = cfg$data_type,
  study_context = cfg$study_context,
  user_notes = cfg$user_notes,
  dataset_scope = cfg$dataset_scope,
  scope_profile = cfg$scope_profile,
  allowed_lineages = cfg$allowed_lineages %||% NULL,
  gate_mode = cfg$gate_mode %||% "flag_only"
)

get_primary_context_value <- function(x, fallback = NULL) {
  if (is.null(x) || length(x) == 0) return(fallback)
  x[[1]]
}

###--- 3c. Other fixed configuration ---
# DisGeNET credential is read from the environment when available
# (the optional disgenet2r dependency is handled when llm_run.R is sourced).

# --- Optional: Step 2 validation (judge) ---
# Set to TRUE only if you explicitly want to generate Step 2 validation queries.
enable_validation <- TRUE

n_deg_per_cluster <- 40 # Number of top-ranked genes supplied to the LLM

# SpeciesPPI。 '9606', mouse is '10090'
primary_species <- get_primary_context_value(prompt_context_config$species) %||% "human"
species_ppi_code <- ifelse(primary_species == "human", "9606", "10090")
ppi_files <- list(
  aliases = paste0(Sys.getenv("TRIAGE_PPI_ROOT", unset = file.path(triageHome, "inputs", "raw", "ppi")), "/", species_ppi_code, ".protein.aliases.v12.0.txt"),
  interactions = paste0(Sys.getenv("TRIAGE_PPI_ROOT", unset = file.path(triageHome, "inputs", "raw", "ppi")), "/", species_ppi_code, ".protein.physical.links.v12.0.txt")
)

###--- 3d. Analysis & RAG Parameters (normally unchanged) ---
analysis_config <- list(
  n_bp = 15, n_cc = 3, n_mf = 3, n_pathway_each = 10, n_tf = 20,
  n_disease_each = 5, n_hubs = 5, litsense_primary_threshold = 0.6,
  litsense_fallback_threshold = 0.5,
  genes_per_sample = 8, n_deg_for_other_analyses = 40,
  top_k_candidates = opt$top_k_candidates,
  disable_smart_context = isTRUE(opt$disable_smart_context),
  eupmc_require_gene_hit = TRUE,
  eupmc_download_workers = 2,
  eupmc_parse_workers = 2,
  cache_only = isTRUE(opt$cache_only)
)

###--- 3e. Dimension keyword library (unchanged) ---
master_keyword_library <- list(
  "base" = list(
    "Biological Processes" = c("biological process", "cellular process", "metabolic process", "GO:BP"),
    "Signaling Pathways" = c("signaling pathway", "signal transduction", "cascade", "KEGG", "Reactome"),
    "Regulatory Mechanisms" = c("transcription factor", "gene regulation", "promoter", "enhancer"),
    "Protein Interaction Networks" = c("protein-protein interaction", "protein complex", "STRING-db", "BioGRID"),
    "Disease-Associated Functions" = c("disease", "syndrome", "pathology", "phenotype")
  ),
  "healthy_immune" = list(
    "Biological Processes" = c("immune response", "lymphocyte activation", "cell differentiation", "cytokine production", "leukocyte migration"),
    "Signaling Pathways" = c("TCR signaling", "BCR signaling", "JAK-STAT pathway", "NF-kappa B signaling"),
    "Regulatory Mechanisms" = c("immune cell differentiation", "tolerance", "somatic recombination"),
    "Disease-Associated Functions" = c("homeostasis", "inflammation", "host defense", "innate immunity", "adaptive immunity")
  ),
  "cancer" = list(
    "Biological Processes" = c("cell proliferation", "apoptosis", "cell migration", "angiogenesis", "DNA damage response"),
    "Signaling Pathways" = c("cell cycle", "PI3K-Akt signaling", "MAPK signaling"),
    "Regulatory Mechanisms" = c("tumor microenvironment", "oncogene activation", "tumor suppressor"),
    "Disease-Associated Functions" = c("cancer", "carcinoma", "neoplasm", "metastasis", "prognosis", "tumorigenesis")
  ),
  "cancer_immunotherapy" = list(
    "Biological Processes" = c("T cell activation", "immune evasion", "antigen processing and presentation"),
    "Signaling Pathways" = c("T cell exhaustion", "interferon-gamma signaling"),
    "Regulatory Mechanisms" = c("tumor infiltrating lymphocytes", "immune checkpoint"),
    "Disease-Associated Functions" = c("immunotherapy", "checkpoint inhibitor", "CAR-T", "treatment response", "resistance")
  ),
  "aml" = list(
    "Biological Processes" = c("myeloid cell differentiation", "hematopoiesis", "regulation of cell cycle"),
    "Regulatory Mechanisms" = c("hematopoietic stem cell differentiation", "leukemic stem cell"),
    "Disease-Associated Functions" = c("leukemia", "AML", "blast cells", "chemotherapy", "relapse", "remission")
  ),
  "development" = list(
    "Biological Processes" = c("embryonic morphogenesis", "pattern specification process", "cell fate commitment", "gastrulation"),
    "Signaling Pathways" = c("Wnt pathway", "Notch signaling", "TGF-beta signaling"),
    "Regulatory Mechanisms" = c("lineage commitment", "organogenesis", "stem cell maintenance"),
    "Disease-Associated Functions" = c("developmental disorder", "congenital abnormality")
  ),
  "neuro_healthy" = list(
    "Biological Processes" = c("synaptic signaling", "neurotransmitter secretion", "axonogenesis", "learning or memory"),
    "Signaling Pathways" = c("synaptic transmission", "calcium signaling", "axon guidance"),
    "Regulatory Mechanisms" = c("neurogenesis", "synaptic plasticity", "long-term potentiation"),
    "Disease-Associated Functions" = c("cognition", "neuronal activity", "behavior")
  ),
  "neuro_degenerative" = list(
    "Biological Processes" = c("protein misfolding", "neuronal apoptosis", "response to oxidative stress", "neuroinflammation"),
    "Signaling Pathways" = c("unfolded protein response", "autophagy"),
    "Regulatory Mechanisms" = c("neuronal death", "gliosis", "microglial activation"),
    "Disease-Associated Functions" = c("neurodegeneration", "Alzheimer's disease", "Parkinson's disease", "amyloid")
  ),
  "fibrosis" = list(
    "Biological Processes" = c("extracellular matrix organization", "collagen fibril organization", "epithelial to mesenchymal transition"),
    "Signaling Pathways" = c("TGF-beta signaling", "Wnt signaling pathway"),
    "Regulatory Mechanisms" = c("myofibroblast differentiation", "tissue remodeling", "wound healing"),
    "Disease-Associated Functions" = c("fibrosis", "cirrhosis", "sclerosis", "scarring")
  )
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ##--- 3f. Dimension keyword library (unchanged) ---
build_dynamic_dimensions <- function(study_context, keyword_library) {
  if (is.null(study_context) || !study_context %in% names(keyword_library)) {
    warning(paste("Study context '", study_context, "' not found. Using 'base' keywords only."))
    return(keyword_library[["base"]])
  }
  base_dims <- keyword_library[["base"]]
  context_dims <- keyword_library[[study_context]]
  final_dims <- purrr::map(names(base_dims), function(dim_name) {
    base_keywords <- base_dims[[dim_name]]
    context_keywords <- context_dims[[dim_name]] %||% c()
    unique_keywords <- unique(c(base_keywords, context_keywords))
    return(unique_keywords)
  })
  names(final_dims) <- names(base_dims)
  return(final_dims)
}

cat("--- Configuration loaded successfully. ---\n")
tictoc::tic(glue::glue("Total processing for Dataset {current_dataset_name}"))

# ===============================================================
# --- 4. Main processing flow (single analysis flow) ---
# ===============================================================
cat(glue::glue("\n\n=========================================================\n"))
cat(glue::glue("### Starting Single Analysis for Dataset: {current_dataset_name} ###\n"))

# --- 4a. Input ---
cat("--- Loading and preprocessing input files... ---\n")
if (!file.exists(deg_file_path)) stop(glue("FATAL: DEG file not found at: {deg_file_path}"))
if (!file.exists(candidates_file_path)) stop(glue("FATAL: Candidates file not found at: {candidates_file_path}"))

# ---- Print column names for sanity check (helps you verify mapping) ----
deg_cols_preview <- data.table::fread(deg_file_path, nrows = 2)
cat("DEG columns detected: ", paste(names(deg_cols_preview), collapse = ", "), "\n")

cand_ext <- tolower(tools::file_ext(candidates_file_path))
cand_cols_preview <- if (cand_ext %in% c("xlsx", "xls")) {
  readxl::read_excel(candidates_file_path, n_max = 2)
} else {
  data.table::fread(candidates_file_path, nrows = 2)
}
cat("Candidates columns detected: ", paste(names(cand_cols_preview), collapse = ", "), "\n")

# ---------------------------
# ✅ MODE A CHANGE:
# Do NOT anonymize clusters here.
# Just standardize columns; llm_run.R will anonymize consistently.
# ---------------------------

# Read the DEG CSV and【rename columns explicitly】columns
processed_deg_data <- data.table::fread(deg_file_path) %>%
  dplyr::rename(
    cluster_name = !!rlang::sym(deg_original_colnames$cluster_col),
    geneSymbol   = !!rlang::sym(deg_original_colnames$gene_col),
    avg_log2FC   = !!rlang::sym(deg_original_colnames$logfc_col)
  ) %>%
  dplyr::mutate(
    cluster_name = as.character(cluster_name),
    geneSymbol   = as.character(geneSymbol),
    avg_log2FC   = as.numeric(avg_log2FC),
    
    # ---- compatibility aliases (critical) ----
    cluster = cluster_name,
    gene    = geneSymbol,
    logfc   = avg_log2FC,
    logFC   = avg_log2FC
  ) %>%
  tibble::as_tibble()

processed_deg_data <- processed_deg_data %>%
  dplyr::mutate(cluster_name_true = cluster_name)

cat(glue::glue("  -> Successfully loaded and standardized DEG file. Rows: {nrow(processed_deg_data)}\n"))

# Read the candidate-cell CSV and【rename columns explicitly】columns
all_candidates_df <- if (cand_ext %in% c("xlsx", "xls")) {
  readxl::read_excel(candidates_file_path)
} else {
  data.table::fread(candidates_file_path)
} %>%
  dplyr::rename(Query_ID = !!rlang::sym(candidates_original_colnames$cluster_col)) %>%
  dplyr::mutate(
    Query_ID = as.character(Query_ID),
    
    # compatibility aliases
    cluster = Query_ID,
    cluster_name = Query_ID
  ) %>%
  tibble::as_tibble()

cat(glue::glue("  -> Successfully loaded and standardized Candidates file. Rows: {nrow(all_candidates_df)}\n"))

# ✅ Ensure stable ordering so cluster_1 mapping is deterministic
clusters_to_process <- all_candidates_df %>%
  dplyr::distinct(Query_ID) %>%
  dplyr::pull(Query_ID)

# --- Mode control: dryrun can limit clusters ---
if (identical(mode, "dryrun")) {
  if (!is.null(opt$cluster_id)) {
    clusters_to_process <- clusters_to_process[clusters_to_process == opt$cluster_id]
  }
  if (!is.null(opt$n_clusters)) {
    clusters_to_process <- head(clusters_to_process, opt$n_clusters)
  } else if (is.null(opt$cluster_id)) {
    clusters_to_process <- head(clusters_to_process, 1)
  }
} else {
  if (!is.null(opt$cluster_id)) clusters_to_process <- clusters_to_process[clusters_to_process == opt$cluster_id]
}

clusters_to_process <- sort(clusters_to_process)

if (!is.null(opt$cluster_id)) {
  if (!(opt$cluster_id %in% clusters_to_process)) {
    stop("cluster_id not found in candidates: ", opt$cluster_id)
  }
}

if (length(clusters_to_process) == 0) {
  stop("FATAL: No clusters found in the 'Query_ID' column of the candidates file. Please check column names.")
}
cat(glue::glue("  -> Found {length(clusters_to_process)} clusters to process.\n"))

# --- 4b. Prepare run parameters ---
dimension_definitions <- build_dynamic_dimensions(prompt_context_config$study_context, master_keyword_library)
run_name <- glue::glue("{current_dataset_name}_LLM_Input_Run")

# ✅ Keep original cluster IDs for filenames
cluster_id_map <- setNames(clusters_to_process, clusters_to_process)

validator_instructions <- if (enable_validation) build_validator_instructions() else NULL

# --- 4c. [STEP 1] Generate expert-evaluation query files ---
cat(glue::glue("\n--- Generating STEP 1 (Expert Report) queries for {current_dataset_name} ---\n"))

expert_report_queries <- generate_expert_report_query(
  processed_deg_data = processed_deg_data,
  all_candidates_df = all_candidates_df,
  clusters_to_process = clusters_to_process,
  n_deg_per_cluster = n_deg_per_cluster,
  ppi_files = ppi_files,
  run_name = run_name,
  analysis_config = analysis_config,
  dimension_definitions = dimension_definitions,
  prompt_context_config = prompt_context_config,
  top_k_candidates = analysis_config$top_k_candidates,
  include_preanalysis = FALSE
)

# --- 4d. Output ---
if (length(expert_report_queries) > 0) {
  
  # ---------- Step 1 ----------
  output_dir_step1 <- file.path(out_root, "step1_report_queries")
  if (!dir.exists(output_dir_step1)) dir.create(output_dir_step1, recursive = TRUE)
  
  purrr::iwalk(expert_report_queries, function(json_object, name) {
    original_cluster_name <- cluster_id_map[[name]] %||% name
    safe_cluster_name <- stringr::str_replace_all(original_cluster_name, "[^a-zA-Z0-9_.-]", "_") %>%
      stringr::str_replace_all("_+", "_")
    file_path <- file.path(output_dir_step1, glue("{safe_cluster_name}_step1_report_query.json"))
    jsonlite::write_json(json_object, file_path, auto_unbox = TRUE, pretty = TRUE, force = TRUE)
  })
  cat(glue::glue("  -> Successfully saved {length(expert_report_queries)} STEP 1 query file(s) to {output_dir_step1}\n"))
  
  # ---------- Step 1.5 ----------
  citation_fix_queries <- purrr::map(expert_report_queries, generate_citation_fix_query)
  
  output_dir_step1_5 <- file.path(out_root, "step1.5_citation_fix_queries")
  if (!dir.exists(output_dir_step1_5)) dir.create(output_dir_step1_5, recursive = TRUE)
  
  purrr::iwalk(citation_fix_queries, function(json_object, name) {
    original_cluster_name <- cluster_id_map[[name]] %||% name
    safe_cluster_name <- stringr::str_replace_all(original_cluster_name, "[^a-zA-Z0-9_.-]", "_") %>%
      stringr::str_replace_all("_+", "_")
    file_path <- file.path(output_dir_step1_5, glue("{safe_cluster_name}_step1.5_citation_fix_query.json"))
    jsonlite::write_json(json_object, file_path, auto_unbox = TRUE, pretty = TRUE, force = TRUE)
  })
  cat(glue::glue("  -> Successfully saved {length(citation_fix_queries)} STEP 1.5 query file(s) to {output_dir_step1_5}\n"))
  
  # ---------- Step 2 (template queries; runner will inject runtime inputs) ----------
  if (isTRUE(enable_validation)) {
    
    if (is.null(validator_instructions)) {
      validator_instructions <- build_validator_instructions()
    }
    
    output_dir_step2 <- file.path(out_root, "step2_validation_queries")
    if (!dir.exists(output_dir_step2)) dir.create(output_dir_step2, recursive = TRUE)
    
    purrr::iwalk(expert_report_queries, function(json_object, name) {
      original_cluster_name <- cluster_id_map[[name]] %||% name
      safe_cluster_name <- stringr::str_replace_all(original_cluster_name, "[^a-zA-Z0-9_.-]", "_") %>%
        stringr::str_replace_all("_+", "_")
      
      # IMPORTANT: runner only needs instructions_for_llm + ids; it will build runtime input_data
      validation_query_object <- list(
        query_id = glue("{name}_validation"),
        analysis_type = "Step 2: Annotation Report Validation",
        instructions_for_llm = validator_instructions
      )
      
      file_path <- file.path(output_dir_step2, glue("{safe_cluster_name}_step2_validation_query.json"))
      jsonlite::write_json(validation_query_object, file_path, auto_unbox = TRUE, pretty = TRUE, force = TRUE)
    })
    
    cat(glue::glue("  -> Successfully saved {length(expert_report_queries)} STEP 2 validation template file(s) to {output_dir_step2}\n"))
    
  } else {
    cat("  -> Step 2 validation generation disabled (enable_validation=FALSE).\n")
  }
  
} else {
  cat("  -> No query files were generated. Check the input files and configurations.\n")
}

tictoc::toc()
