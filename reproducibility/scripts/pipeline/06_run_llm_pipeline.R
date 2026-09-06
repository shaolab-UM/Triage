#!/usr/bin/env Rscript
# =========================================================
# Parallel CellAnnotation Pipeline (cluster-level parallel)
# Step1 -> Step1.5 -> Gate -> (optional Step2 auto-gen if missing)
# + timing + token/usage logging
# + per-stage model routing (step1/step1.5/step2/repair)
#
# Step-2 behavior:
# - enable_step2 defaults to TRUE.
# - When enabled and the Step-2 query file is absent, validator instructions are generated automatically.
# =========================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(glue)
  library(stringr)
  library(purrr)
  library(readr)
  library(future)
  library(furrr)
  library(dplyr)
  library(optparse)
})

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1) {
    if (is.na(a)) return(b)
    if (is.character(a) && !nzchar(a)) return(b)
  }
  a
}
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
if (!file.exists(file.path(triageHome, "config", "cl_normalizer.R"))) {
  triageHome <- Sys.getenv("PROJECT_ROOT", unset = triageHome)
}
suppressMessages(library(Triage))
suppressMessages(library(Triage))

# =========================
# 0) CASSIA-style Step2 validator instructions (JSON output)
# =========================
build_validator_instructions_cassia <- function() {
  list(
    role = "You are a strict but practical QC validator for LLM-generated cell type annotations.",
    task = paste(
      "You will be given an input_data object with:",
      "1) report_to_be_validated (the generated annotation JSON),",
      "2) original_dossier (DEGs + bioinfo + literature evidence).",
      "Your job is to decide PASS/FAIL using a CASSIA-like philosophy:",
      "PASS if the CORE LINEAGE / CORE IDENTITY is supported and there are no hard contradictions.",
      "DO NOT fail solely because of weak/incorrect citations or over-strong wording in phenotypic details—log those as AUDIT WARNINGS instead.",
      "IMPORTANT: validation_status MUST reflect ONLY core identity (main/sub1/sub2.core_identity).",
      "Note: You may adjust identity fields (candidate_cell_type / cell_ontology_id / phenotypic_label) to fix core identity issues; do NOT rewrite narrative text beyond those fields."
    ),
    
    pass_fail_policy = list(
      pass_definition = c(
        "PASS if main_type and subtype_level_1 (and subtype_level_2.core_identity) are consistent with DEG markers AND not contradicted by bioinfo.",
        "PASS even if phenotype/state claims are weakly supported, as long as they are not directly contradictory.",
        "validation_status = VALIDATION FAILED ONLY when core identity fails."
      ),
      fail_conditions = c(
        "FAIL if the proposed CORE IDENTITY (main/sub1/sub2.core_identity) contradicts the DEG evidence.",
        "FAIL if the report asserts key markers that are NOT present in the provided DEG list as if they were observed.",
        "FAIL if the report contains direct contradictions against the dossier (e.g., claims pathway enrichment that does not appear in bioinfo lists).",
        "FAIL if there is SCHEMA CONTRADICTION: main_type/subtype_level_1/subtype_level_2.core_identity and open_world_summary.best_cell_type imply different broad lineages.",
        "FAIL if main_type/subtype_level_1 relies on low-rank markers (>20 in degs) as PRIMARY evidence while in_scope_top_genes support a different lineage/program.",
        "FAIL if main_type/subtype_level_1 is NOT supported by any lineage-consistent markers within in_scope_top_genes or ranks_16_50_in_scope. If in_scope_top_genes are empty or lack coherent anchors, allow ranks_16_50_in_scope to pick the best-supported in-scope lineage and do NOT FAIL unless no plausible lineage exists across the evidence summary.",
        "FAIL if core_identity implies a coarse lineage incompatible with dataset_scope/tissue constraints when the context indicates a restricted lineage scope."
      ),
      
      non_fail_audit_items = c(
        "Over-strong marker language for non-canonical genes -> AUDIT WARNING.",
        "Weak/indirect citations for phenotype/state claims -> AUDIT WARNING.",
        "Wrong or weak citation-to-claim match (PMID exists but does not support the specific claim) -> AUDIT WARNING.",
        "Uncited statements that are clearly traceable to bioinfo fields -> AUDIT WARNING + suggest adding explicit source pointer.",
        "Confidence calibration/hierarchy issues that do not change identity -> AUDIT WARNING (recommend fix).",
        "If mixed signals are present, choose the single best-supported lineage for core identity and record ambiguity only in phenotypic_label. Mixed/ambiguous core identities are FINAL ROUND last resort only after in_scope_top_genes + ranks_16_50_in_scope review.",
        "If out_of_scope_top_genes are present, treat them as contamination and do NOT use them to fail core identity.",
        "In restricted scopes, do not treat out_of_scope_top_genes as coherent lineage evidence if in_scope markers exist in ranks_16_50_in_scope; log as contamination instead.",
        "unknown_top_genes are contamination-only; do NOT use them to override core identity.",
        "If in_scope_top_genes are weak but not empty, prefer a conservative in-scope parent over unknown/other and record uncertainty in phenotypic_label.",
        "Global summary is for consistency checks only; do NOT override in_scope_top_genes with global summary signals."
      )
    ),
    
    rules = c(
      "**0) Layered Workflow Rule:** (1) Use Top50 + recovery tier to score main lineage candidates (do NOT hard-lock from a single tier). (2) Refine subtype only if evidence is sufficient. (3) If subtype evidence is insufficient, stop at the parent lineage (do NOT fall back to generic 'cell').",
      "**0a) Evidence Coverage Rule:** Use all available ranked in-scope genes (Top50 + recovery tier) to score candidates; do NOT rely on Top15 alone.",
      "**0b) Weighted Candidate Ranking:** Compute a weighted score for plausible lineages/subtypes (Top50 canonical = 2 points, recovery tier canonical = 1 point). Unknown markers may provide weak support (weight <=0.5) but cannot define lineage alone. Choose the highest-scoring in-scope candidate.",
      "**0c) Subtype Threshold Rule:** Do NOT finalize a subtype unless there are >=2 subtype-specific markers supporting it. A single marker is only a hint.",
      "**0d) Evidence Routing:** Top50 = strong evidence, recovery tier = weak evidence, unknown = weak support only (cannot define lineage alone), out-of-scope = contradiction only.",
      "**0f) Confidence Calibration:** Set confidence_score_breakdown strictly from evidence. Do NOT force main_type confidence to exceed subtype confidence; subtype may be higher if evidence is stronger.",
      "**0e) Candidate Consistency Rule:** Candidate list is consistency-only; do not override evidence with candidates.",
      "**0c) Contamination Override:** If out_of_scope_top_genes are present but ranks_16_50_in_scope contain coherent in-scope lineage markers with supporting bioinfo/literature, prefer the in-scope lineage and treat out_of_scope_top_genes as contamination.",
      "**0d) Dominant Lineage Preservation:** If multiple canonical markers support a coherent lineage, keep that lineage unless an alternative lineage has stronger, coherent support across in_scope_top_genes + ranks_16_50_in_scope.",
      "**0.5) Schema Consistency Check:** Ensure all core identity fields and open_world_summary agree on the same broad lineage.",
      "**1) Core Identity Check (Primary):** Verify main_type, subtype_level_1, and subtype_level_2.core_identity against DEGs. Decide whether lineage is correct.",
      "**2) Contradiction Check (Hard):** Identify any direct contradictions against DEGs/bioinfo. Contradictions can trigger FAIL.",
      "**3) Citation Integrity Check:** Do NOT require every claim to be cited. Do NOT FAIL for citation issues. Log mismatches as audit warnings.",
      "**4) Bioinfo Traceability:** If the report mentions KEGG/Reactome/TF results, check they exist in original_dossier.bioinfo. Missing traceability -> audit warning, not fail.",
      "**5) Confidence Logic:** If confidence scores violate hierarchy constraints or look uncalibrated, flag as audit warning unless it fundamentally misrepresents identity.",
      "**6) Output:** Always return JSON in the schema below."
    ),
    
    output_format = list(
      description = "Return a single valid JSON object. PASS/FAIL is for core identity only; audit is separate.",
      schema = list(
        validation_status = "string (Either 'VALIDATION PASSED' or 'VALIDATION FAILED')",
        failure_reasons = "list[string] (Only include if FAILED; concise, evidence-based)",
        core_identity_verdict = list(
          main_type_ok = "boolean",
          subtype_level_1_ok = "boolean",
          subtype_level_2_core_ok = "boolean",
          core_identity_summary = "string (1-3 sentences explaining why identity is or isn't supported)"
        ),
        verdicts = list(
          core_identity = "PASS|FAIL",
          marker_rank_support = "PASS|FAIL",
          schema_consistency = "PASS|FAIL",
          bioinfo_traceability = "PASS|WARN",
          confidence_calibration = "PASS|WARN"
        ),
        audit_warnings = "list[string] (Non-fatal issues: weak citations, overstated marker claims, missing bioinfo pointers, etc.)",
        suggested_fixes = "list[string] (Concrete edits to improve the report without changing core identity)",
        final_verdict_summary = "string (One sentence conclusion; PASS/FAIL reflects core identity only)",
        revision_directive = list(
          description = "REQUIRED when VALIDATION FAILED. Classify the failure type to enable automated triage.",
          schema = list(
            validation_outcome = "string (pass | revise | hard_reject)",
            failure_type = "string (A_conservative_ok | B_overspecific | C_wrong_lineage | null if passed)",
            action = "string (accept_as_is | downgrade_to_parent | stop)",
            suggested_parent = "string | null (CL label of suggested parent if downgrade needed)",
            suggested_parent_clid = "string | null (CL:XXXXXXX format if known)",
            notes = "string (brief explanation of the classification)"
          ),
          classification_rules = c(
            "A_conservative_ok: Prediction is correct DIRECTION but not specific enough (e.g., 'epithelial cell' when truth is 'RPE'). Prediction is an ANCESTOR of the true type. Action: accept_as_is.",
            "B_overspecific: Prediction is correct DIRECTION but TOO specific (e.g., 'IT-projecting glutamatergic neuron' when evidence only supports 'glutamatergic neuron'). Prediction is a DESCENDANT of what evidence supports. Action: downgrade_to_parent.",
            "C_wrong_lineage: Prediction is WRONG DIRECTION entirely (e.g., 'melanocyte' when markers support 'epithelial'). Neither ancestor nor descendant relationship. Action: stop."
          )
        )
      )
    )
  )
}

# =========================
# CLI
# =========================
option_list <- list(
  make_option(c("--mode"), type="character", default="full", help="full | dryrun [default %default]"),
  make_option(c("--input_root"), type="character", default=NULL, help="Root dir containing step1/step1.5/step2 queries."),
  make_option(c("--output_root"), type="character", default=NULL, help="Root dir to write outputs."),
  make_option(c("--dataset_name"), type="character", default=NULL, help="Dataset name for default paths (default: basename(getwd()))."),
  make_option(c("--project_root"), type="character", default=getwd(), help="Project root (resource resolution)."),
  make_option(c("--clusters"), type="integer", default=1, help="In dryrun: process first N clusters [default %default]"),
  make_option(c("--cluster_id"), type="character", default=NULL, help="If set, only process this cluster id prefix."),
  make_option(c("--workers"), type="integer", default=12, help="Parallel workers [default %default]"),
  make_option(c("--disable_step2"), action="store_true", default=FALSE, help="Disable Step2 validation [default FALSE]"),
  make_option(c("--enable_repair"), action="store_true", default=FALSE, help="Enable optional repair stage [default %default]"),
  make_option(c("--max_rounds"), type="integer", default=3, help="Max validate/rewrite rounds [default %default]"),
  make_option(c("--cl_local_json"), type="character", default=file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"),
              help="Local CL ontology JSON path [default %default]"),
  make_option(c("--ols_first"), type="logical", default=TRUE,
              help="Use OLS first for CL normalization [default TRUE]"),
  make_option(c("--ols_cache_dir"), type="character", default="",
              help="OLS cache dir (default: <output_root>/.ols_cache)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Ensure resource files are resolvable in future workers
Sys.setenv(PROJECT_ROOT = opt$project_root)
if (!is.null(opt$project_root) && dir.exists(opt$project_root)) {
  try(setwd(opt$project_root), silent = TRUE)
}
mode <- tolower(opt$mode)
dataset_name <- opt$dataset_name %||% Sys.getenv("CASSIA_DATASET", unset = "")
if (!nzchar(dataset_name)) dataset_name <- basename(getwd())
cfg <- Triage:::get_dataset_config(dataset_name, opt$project_root)

# =========================
# 0) CONFIG
# =========================
config <- list(
  api_key_env = Sys.getenv("LLM_API_KEY_ENV", unset = Sys.getenv("CASSIA_API_KEY_ENV", unset = "DEEPSEEK_API_KEY")),
  api_base_url = Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = "XXXXX")),
  deepseek_api_key = Sys.getenv(Sys.getenv("LLM_API_KEY_ENV", unset = Sys.getenv("CASSIA_API_KEY_ENV", unset = "DEEPSEEK_API_KEY"))),
  deepseek_model   = "deepseek-v4-flash",
  
  # ---- model routing ----
  model_default = Sys.getenv("LLM_MODEL_STEP1", unset = Sys.getenv("CASSIA_MODEL_STEP1", unset = "deepseek-v4-flash")),
  model_step1   = Sys.getenv("LLM_MODEL_STEP1", unset = Sys.getenv("CASSIA_MODEL_STEP1", unset = "deepseek-v4-flash")),
  model_step15  = Sys.getenv("LLM_MODEL_STEP15", unset = Sys.getenv("CASSIA_MODEL_STEP15", unset = "deepseek-v4-flash")),
  model_step2   = Sys.getenv("LLM_MODEL_STEP2", unset = Sys.getenv("CASSIA_MODEL_STEP2", unset = "deepseek-v4-flash")),
  model_repair  = Sys.getenv("LLM_MODEL_REPAIR", unset = Sys.getenv("CASSIA_MODEL_REPAIR", unset = "deepseek-v4-flash")),
  
  temperature = 0.0,
  max_retries = 5,
  timeout_seconds = 1200,
  retry_delay_seconds = 5,
  stream = FALSE,
  
  input_root = cfg$llm_inputs_root,
  in_step1   = "step1_report_queries",
  in_step15  = "step1.5_citation_fix_queries",
  in_step2   = "step2_validation_queries",
  
  output_root = cfg$llm_outputs_root,
  out_step1   = "step1_report_outputs",
  out_step15  = "step1.5_citation_fix_outputs",
  out_step2   = "step2_validation_outputs",
  out_final   = "final_passed",
  out_debug   = "debug_failed",
  out_logs    = "run_logs",
  
  # CL normalization
  cl_local_json = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"),
  ols_first = TRUE,
  ols_cache_dir = "",
  
  # loops
  max_rounds = 3,
  
  # evidence pruning
  prune_step1  = TRUE,
  prune_step15 = TRUE,
  prune_step2  = TRUE,
  keep_per_dimension_step1  = 3,
  keep_per_dimension_step15 = 3,
  keep_per_dimension_step2  = 3,
  keep_chars_step1  = 1200,
  keep_chars_step15 = 1200,
  keep_chars_step2  = 1200,
  
  # Gate v1
  gate_pmid_overuse_threshold = 8,
  
  # Parallelism
  parallel_workers = 12,
  parallel_seed = TRUE,
  
  # flags (DEFAULT: Step2 enabled)
  enable_step2  = TRUE,
  enable_repair = FALSE,
  
  # Single-cluster run
  test_one_cluster = FALSE,
  test_cluster_prefix = NULL
)

# --- CLI overrides on config ---
if (!is.null(opt$input_root)) config$input_root <- opt$input_root
if (!is.null(opt$output_root)) config$output_root <- opt$output_root
if (isTRUE(opt$disable_step2)) {
  config$enable_step2 <- FALSE
}
config$enable_repair <- isTRUE(opt$enable_repair)
config$max_rounds    <- as.integer(opt$max_rounds)
config$cl_local_json  <- opt$cl_local_json
config$ols_first      <- isTRUE(opt$ols_first)
config$ols_cache_dir  <- opt$ols_cache_dir
if (!nzchar(config$ols_cache_dir)) config$ols_cache_dir <- file.path(config$output_root, ".ols_cache")
workers <- max(1L, as.integer(opt$workers))
config$parallel_workers <- workers

# --- Dryrun / cluster selection ---
if (!is.null(opt$cluster_id)) {
  config$test_one_cluster <- TRUE
  config$test_cluster_prefix <- opt$cluster_id
}

cat(glue("[CONFIG] enable_step2={config$enable_step2}\n"))
cat(glue("[CONFIG] workers={config$parallel_workers}\n"))
if (is.null(config$deepseek_api_key) || !nzchar(config$deepseek_api_key)) {
  cat("[WARN] API key is empty. Set env ", config$api_key_env, " before running.\n", sep = "")
}

cl_cfg <- Triage:::make_cl_cfg(config$cl_local_json, prefer_ols = config$ols_first, cache_dir = config$ols_cache_dir)
prompt_context_config <- list(
  species = cfg$species,
  tissue = cfg$tissue,
  dataset_scope = cfg$dataset_scope,
  allowed_lineages = cfg$allowed_lineages %||% NULL,
  gate_mode = cfg$gate_mode %||% "flag_only"
)

get_primary_context_value <- function(x, fallback = "") {
  if (is.null(x) || length(x) == 0) return(fallback)
  if (is.list(x)) {
    v <- x[[1]]
    if (is.null(v)) return(fallback)
    return(as.character(v))
  }
  if (is.character(x)) return(x[[1]])
  fallback
}

enforce_scope_coarse_gate <- function(report_json, prompt_context_list) {
  if (is.null(report_json) || !is.list(report_json)) return(report_json)

  allowed <- prompt_context_list$allowed_lineages %||% NULL
  if (is.null(allowed) || length(allowed) == 0) return(report_json)

  gate_mode <- tolower(prompt_context_list$gate_mode %||% "flag_only")

  normalize_lineage_for_scope <- function(lineage, allowed) {
    allowed <- allowed %||% character(0)
    if ("neural_glial" %in% allowed && lineage %in% c("neural", "glial")) return("neural_glial")
    lineage
  }

  normalize_conservative_parent <- function(lineage) {
    if (!nzchar(lineage) || lineage == "unknown") return("cell (ambiguous lineage)")
    if (grepl("cell", lineage, ignore.case = TRUE)) return(lineage)
    paste(lineage, "cell")
  }

  get_lin <- function(lbl, clid = "") {
    if (exists("infer_coarse_lineage", mode = "function")) {
      out <- tryCatch(infer_coarse_lineage(lbl %||% "", clid %||% ""), error = function(e) "unknown")
    } else {
      out <- infer_coarse_lineage_simple(lbl %||% "", clid %||% "", cl_cfg)
    }
    out <- normalize_lineage_for_scope(out, allowed)
    if (!nzchar(out)) "unknown" else out
  }

  main_lin <- get_lin(report_json$main_type$candidate_cell_type, report_json$main_type$cell_ontology_id)
  sub1_lin <- get_lin(report_json$subtype_level_1$candidate_cell_type, report_json$subtype_level_1$cell_ontology_id)
  sub2_lin <- get_lin(report_json$subtype_level_2$core_identity$candidate_cell_type, report_json$subtype_level_2$core_identity$cell_ontology_id)

  bad <- !(main_lin %in% allowed && sub1_lin %in% allowed && sub2_lin %in% allowed)

  if (isTRUE(bad)) {
    if (identical(gate_mode, "flag_only")) {
      report_json$subtype_level_2$phenotypic_label <- paste(
        report_json$subtype_level_2$phenotypic_label %||% "",
        "| NOTE: scope gate flag (out-of-scope core identity detected).",
        sep = " "
      )
      return(report_json)
    }

    pick_in_scope <- function(...) {
      vals <- c(...)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      for (v in vals) {
        lin <- get_lin(v, "")
        if (lin %in% allowed) return(v)
      }
      ""
    }

    candidate_parent <- pick_in_scope(
      report_json$subtype_level_2$core_identity$candidate_cell_type,
      report_json$subtype_level_1$candidate_cell_type,
      report_json$main_type$candidate_cell_type,
      report_json$open_world_summary$best_cell_type
    )
    fallback_lineage <- if (length(allowed) > 0) allowed[[1]] else "unknown"
    conservative <- if (nzchar(candidate_parent)) candidate_parent else normalize_conservative_parent(fallback_lineage)

    main_oos <- !(main_lin %in% allowed)
    sub1_oos <- !(sub1_lin %in% allowed)
    sub2_oos <- !(sub2_lin %in% allowed)

    if (main_oos) {
      report_json$main_type$candidate_cell_type <- conservative
      report_json$main_type$cell_ontology_id <- ""
    }
    if (sub1_oos) {
      report_json$subtype_level_1$candidate_cell_type <- conservative
      report_json$subtype_level_1$cell_ontology_id <- ""
    }
    if (sub2_oos) {
      report_json$subtype_level_2$core_identity$candidate_cell_type <- conservative
      report_json$subtype_level_2$core_identity$cell_ontology_id <- ""
    }

    report_json$subtype_level_2$phenotypic_label <- paste(
      report_json$subtype_level_2$phenotypic_label %||% "",
      "| NOTE: scope-coarse gate triggered; core identity downgraded to a conservative parent.",
      sep = " "
    )
  }

  report_json
}

pick_conservative_parent_from_summary <- function(step1_query_obj, prompt_context_config) {
  allowed <- prompt_context_config$allowed_lineages %||% character(0)
  allowed_effective <- setdiff(allowed, "unknown")
  if (length(allowed_effective) == 0) return("cell (ambiguous lineage)")

  summary <- step1_query_obj$input_data$cluster_dossier$global_summary$evidence_lineage_summary %||% list()
  counts <- summary$lineage_counts %||% list()
  if (length(counts) == 0) return("cell (ambiguous lineage)")

  best_lineage <- NULL
  best_count <- -Inf
  for (lin in names(counts)) {
    if (!(lin %in% allowed_effective)) next
    val <- suppressWarnings(as.numeric(counts[[lin]]))
    if (is.na(val)) next
    if (val > best_count) {
      best_count <- val
      best_lineage <- lin
    }
  }

  if (is.null(best_lineage) || !nzchar(best_lineage)) {
    if ("immune" %in% allowed_effective) return("immune cell")
    return(paste(allowed_effective[[1]], "cell"))
  }

  if (grepl("cell", best_lineage, ignore.case = TRUE)) return(best_lineage)
  paste(best_lineage, "cell")
}

# =========================
# 1) IO helpers
# =========================
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

save_json_pretty <- function(obj, path) {
  ensure_dir(dirname(path))
  jsonlite::write_json(obj, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

read_json_safely <- function(path) {
  tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) stop("Failed JSON: ", path, " | ", conditionMessage(e))
  )
}

apply_core_downgrade <- function(step_obj, parent_label, parent_clid) {
  if (is.null(step_obj) || !nzchar(parent_label)) return(step_obj)
  if (!is.null(step_obj$subtype_level_1_schema)) {
    step_obj$subtype_level_1_schema$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$subtype_level_1_schema$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$subtype_level_1)) {
    step_obj$subtype_level_1$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$subtype_level_1$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$subtype_level_2_schema$core_identity)) {
    step_obj$subtype_level_2_schema$core_identity$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$subtype_level_2_schema$core_identity$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$subtype_level_2$core_identity)) {
    step_obj$subtype_level_2$core_identity$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$subtype_level_2$core_identity$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$main_type_schema) &&
      !is.null(step_obj$subtype_level_1_schema$candidate_cell_type) &&
      identical(step_obj$main_type_schema$candidate_cell_type, step_obj$subtype_level_1_schema$candidate_cell_type)) {
    step_obj$main_type_schema$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$main_type_schema$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$main_type)) {
    step_obj$main_type$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$main_type$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$open_world_summary)) {
    step_obj$open_world_summary$best_cell_type <- parent_label
    if (nzchar(parent_clid)) step_obj$open_world_summary$cell_ontology_id <- parent_clid
  }
  if (!is.null(step_obj$subtype_level_2$phenotypic_label)) {
    step_obj$subtype_level_2$phenotypic_label <- paste(
      step_obj$subtype_level_2$phenotypic_label,
      "| NOTE: downgraded to conservative parent.",
      sep = " "
    )
  }
  step_obj
}

# =========================
# 2) DeepSeek API call
# =========================
invoke_deepseek_api <- function(prompt_json_string, config, model = NULL) {
  base_url <- config$api_base_url %||% "XXXXX"
  retries <- 0
  delay <- config$retry_delay_seconds
  model_to_use <- model %||% config$model_default %||% config$deepseek_model
  
  while (retries < config$max_retries) {
    req_body <- list(
      model = model_to_use,
      messages = list(
        list(
          role = "system",
          content = paste(
            "Return ONLY valid JSON (no markdown fences, no extra text).",
            "Preserve the input schema as much as possible; fill missing fields with null.",
            "Do not add new keys unless necessary."
          )
        ),
        list(role = "user", content = prompt_json_string)
      ),
      temperature = config$temperature,
      top_p = 1,
      presence_penalty = 0,
      frequency_penalty = 0,
      stream = config$stream
    )
    
    response <- tryCatch({
      httr::POST(
        url = base_url,
        httr::add_headers(
          `Content-Type`  = "application/json",
          `Authorization` = paste("Bearer", config$deepseek_api_key)
        ),
        body = jsonlite::toJSON(req_body, auto_unbox = TRUE, null = "null"),
        encode = "raw",
        httr::timeout(config$timeout_seconds)
      )
    }, error = function(e) e)
    
    if (inherits(response, "error")) {
      retries <- retries + 1
      Sys.sleep(delay)
      next
    }
    
    status <- httr::status_code(response)
    if (status == 200) {
      content <- httr::content(response, as = "parsed")
      out <- tryCatch(content$choices[[1]]$message$content, error = function(e) NULL)
      usage <- content$usage %||% NULL
      if (!is.null(out) && nchar(out) > 0) {
        return(list(ok = TRUE, status = status, text = out, usage = usage, model = model_to_use))
      }
    }
    
    retries <- retries + 1
    Sys.sleep(delay)
  }
  
  list(ok = FALSE, status = NA_integer_, text = NULL, usage = NULL, model = model_to_use)
}

# =========================
# 3) Robust JSON extraction
# =========================
extract_first_json_object_stack <- function(text) {
  if (is.null(text) || text == "") return(NA_character_)
  ok <- tryCatch({ jsonlite::fromJSON(text, simplifyVector = FALSE); TRUE }, error = function(e) FALSE)
  if (ok) return(text)
  
  text2 <- stringr::str_replace_all(text, "^```json\\s*|\\s*```$", "")
  ok2 <- tryCatch({ jsonlite::fromJSON(text2, simplifyVector = FALSE); TRUE }, error = function(e) FALSE)
  if (ok2) return(text2)

  score_obj <- function(obj) {
    if (is.null(obj) || !is.list(obj)) return(0)
    nms <- names(obj) %||% character(0)
    s <- 0
    # Step2 validator outputs
    if ("validation_status" %in% nms) s <- s + 100
    # Judge outputs (not used here, but safe)
    if ("final_decision" %in% nms) s <- s + 80
    # Step1/Step1.5 outputs
    if ("main_type" %in% nms) s <- s + 40
    if ("subtype_level_1" %in% nms) s <- s + 10
    if ("subtype_level_2" %in% nms) s <- s + 10
    # Penalize schema-only blobs
    if (("main_type_schema" %in% nms) && !("main_type" %in% nms)) s <- s - 30
    if (("subtype_level_1_schema" %in% nms) && !("subtype_level_1" %in% nms)) s <- s - 10
    if (("subtype_level_2_schema" %in% nms) && !("subtype_level_2" %in% nms)) s <- s - 10

    mt <- obj$main_type %||% NULL
    if (is.list(mt)) {
      cand <- mt$candidate_cell_type %||% NA_character_
      cand <- as.character(cand)[1]
      if (!is.na(cand) && nzchar(stringr::str_trim(cand))) s <- s + 50
    }

    s
  }
  
  chars <- strsplit(text2, "", fixed = TRUE)[[1]]
  start <- NA_integer_
  depth <- 0L
  in_str <- FALSE
  esc <- FALSE

  found <- list()
  
  for (i in seq_along(chars)) {
    ch <- chars[[i]]
    if (in_str) {
      if (esc) esc <- FALSE
      else if (ch == "\\") esc <- TRUE
      else if (ch == "\"") in_str <- FALSE
      next
    } else {
      if (ch == "\"") { in_str <- TRUE; next }
    }
    
    if (!in_str) {
      if (ch == "{") {
        if (depth == 0L && is.na(start)) start <- i
        depth <- depth + 1L
      } else if (ch == "}") {
        if (depth > 0L) depth <- depth - 1L
        if (depth == 0L && !is.na(start)) {
          candidate <- paste0(chars[start:i], collapse = "")
          obj3 <- tryCatch(jsonlite::fromJSON(candidate, simplifyVector = FALSE), error = function(e) NULL)
          if (!is.null(obj3)) {
            found[[length(found) + 1L]] <- list(text = candidate, obj = obj3, end_i = i)
          }
          start <- NA_integer_
        }
      }
    }
  }

  if (length(found) == 0) return(NA_character_)
  scores <- vapply(found, function(x) score_obj(x$obj), numeric(1))
  ends <- vapply(found, function(x) as.integer(x$end_i %||% 0L), integer(1))
  pick <- order(scores, ends, decreasing = TRUE)[1]
  found[[pick]]$text
}

# =========================
# 4) Evidence pruning
# =========================
prune_step1_evidence <- function(x, keep_per_dim = 3, keep_chars = 1200) {
  ra <- x$input_data$cluster_dossier$evidence$retrieved_articles
  if (is.null(ra) || is.null(ra$articles_db) || is.null(ra$articles_db$pmid) || is.null(ra$articles_db$text)) {
    return(x)
  }
  
  pmid_vec  <- as.character(ra$articles_db$pmid)
  text_vec  <- ra$articles_db$text
  
  corr_vec <- ra$articles_db$correlation %||%
    ra$articles_db$corr %||%
    ra$articles_db$corr_score %||%
    ra$articles_db$correlation_score %||% NULL
  
  score_vec <- ra$articles_db$score %||% rep(NA_real_, length(pmid_vec))
  rank_primary <- if (!is.null(corr_vec)) as.numeric(corr_vec) else as.numeric(score_vec)
  
  keep_pmids <- NULL
  if (!is.null(ra$relevance_map) && is.list(ra$relevance_map) && length(ra$relevance_map) > 0) {
    keep_pmids <- unique(unlist(lapply(ra$relevance_map, function(v) head(as.character(v), keep_per_dim))))
  }
  
  if (is.null(keep_pmids) || length(keep_pmids) == 0) {
    ord <- order(rank_primary, decreasing = TRUE, na.last = TRUE)
    keep_pmids <- head(pmid_vec[ord], keep_per_dim * 5)
  }
  
  keep_idx <- which(pmid_vec %in% keep_pmids)
  pmid_vec2  <- pmid_vec[keep_idx]
  text_vec2  <- text_vec[keep_idx]
  score_vec2 <- score_vec[keep_idx]
  
  text_vec2 <- vapply(text_vec2, function(t) {
    if (is.null(t) || is.na(t)) return(NA_character_)
    t <- as.character(t)
    if (nchar(t) > keep_chars) substr(t, 1, keep_chars) else t
  }, character(1))
  
  x$input_data$cluster_dossier$evidence$retrieved_articles$articles_db <- list(
    pmid  = pmid_vec2,
    text  = text_vec2,
    score = score_vec2
  )
  
  if (!is.null(ra$relevance_map) && is.list(ra$relevance_map)) {
    x$input_data$cluster_dossier$evidence$retrieved_articles$relevance_map <- lapply(
      ra$relevance_map,
      function(v) {
        vv <- as.character(v)
        vv <- vv[vv %in% pmid_vec2]
        head(vv, keep_per_dim)
      }
    )
  }
  x
}

# =========================
# 5) filename helpers
# =========================
file_prefix_from_step1 <- function(step1_file) {
  bn <- basename(step1_file)
  sub("_step1_report_query\\.json$", "", bn)
}

# =========================
# 6) Step1->Step1.5 injection
# =========================
extract_allowed_pmids <- function(step15_query_obj, step1_query_obj) {
  pmids <- step15_query_obj$input_data$literature_evidence$retrieved_articles$articles_db$pmid %||%
    step1_query_obj$input_data$cluster_dossier$evidence$retrieved_articles$articles_db$pmid
  unique(as.character(pmids %||% character(0)))
}

inject_step15_report_and_guardrails <- function(step15_query_obj, step1_output_obj, allowed_pmids) {
  if (is.null(step15_query_obj$input_data) || !is.list(step15_query_obj$input_data)) {
    step15_query_obj$input_data <- list()
  }
  step15_query_obj$input_data$report_to_be_corrected <- step1_output_obj
  step15_query_obj$input_data$allowed_pmids <- allowed_pmids
  
  extra_rules <- c(
    "HALLUCINATION_HUNTER_RULES:",
    "1) You MUST ONLY use PMIDs present in input_data.allowed_pmids.",
    "2) Do NOT invent any PMID. If you cannot find support, leave the sentence uncited (do not add any tag).",
    "3) Do NOT change any text content; only append citation tags like [pmid:12345].",
    "4) Prefer the most directly supporting PMID(s); do not overuse one PMID everywhere."
  )
  
  if (is.character(step15_query_obj$instructions_for_llm)) {
    step15_query_obj$instructions_for_llm <- paste(step15_query_obj$instructions_for_llm,
                                                   paste(extra_rules, collapse = "\n"),
                                                   sep = "\n\n")
  } else if (is.list(step15_query_obj$instructions_for_llm)) {
    step15_query_obj$instructions_for_llm$hallucination_hunter_rules <- extra_rules
  } else {
    step15_query_obj$instructions_for_llm <- extra_rules
  }
  step15_query_obj
}

append_instruction_text <- function(prompt_obj, extra_text) {
  if (is.null(extra_text) || !nzchar(extra_text)) return(prompt_obj)
  if (is.list(prompt_obj$instructions_for_llm)) {
    prompt_obj$instructions_for_llm$round_guardrail <- extra_text
  } else if (is.character(prompt_obj$instructions_for_llm)) {
    prompt_obj$instructions_for_llm <- paste(prompt_obj$instructions_for_llm, extra_text, sep = "\n\n")
  } else {
    prompt_obj$instructions_for_llm <- extra_text
  }
  prompt_obj
}

# =========================
# 7) Step2 runtime prompt (auto-gen template supported)
# =========================
build_step2_runtime_prompt <- function(step2_instructions, step15_output_obj, step1_query_obj, query_id = "validation") {
  list(
    query_id = query_id,
    analysis_type = "Step 2: Annotation Report Validation",
    instructions_for_llm = step2_instructions,
    input_data = list(
      report_to_be_validated = step15_output_obj,
        original_dossier = step1_query_obj$input_data$cluster_dossier %||% step1_query_obj$cluster_dossier
    )
  )
}

get_validation_status <- function(step2_out) {
  if (is.null(step2_out$validation_status)) return(NA_character_)
  as.character(step2_out$validation_status)
}

# =========================
# 8) Gate
# =========================
extract_pmids_from_report <- function(report_obj) {
  txt <- jsonlite::toJSON(report_obj, auto_unbox = TRUE, null = "null")
  txt <- stringr::str_replace_all(txt, "[\\u3010\\uFF3B]", "[")
  txt <- stringr::str_replace_all(txt, "[\\u3011\\uFF3D]", "]")
  m1 <- stringr::str_match_all(txt, "\\[pmid:([0-9]+)\\]")
  m2 <- stringr::str_match_all(txt, "\\(pmid:([0-9]+)\\)")
  pmids <- unique(c(m1[[1]][,2] %||% character(0), m2[[1]][,2] %||% character(0)))
  pmids[!is.na(pmids) & pmids != ""]
}

citation_gate <- function(step15_report_obj, allowed_pmids, overuse_threshold = 8) {
  allowed_pmids <- unique(as.character(allowed_pmids %||% character(0)))
  used_pmids <- extract_pmids_from_report(step15_report_obj)
  invalid_pmids <- setdiff(used_pmids, allowed_pmids)
  
  txt <- jsonlite::toJSON(step15_report_obj, auto_unbox = TRUE, null = "null")
  pmid_counts <- purrr::map_int(used_pmids, function(x) {
    stringr::str_count(txt, paste0("\\[pmid:", x, "\\]")) +
      stringr::str_count(txt, paste0("\\(pmid:", x, "\\)"))
  })
  names(pmid_counts) <- used_pmids
  overused <- names(pmid_counts[pmid_counts >= overuse_threshold])
  
  list(
  gate_passed = (length(invalid_pmids) == 0 && length(overused) == 0),
    used_pmids = used_pmids,
    invalid_pmids = invalid_pmids,
    pmid_overuse = overused,
    pmid_counts = pmid_counts
  )
}

# =========================
# 9) Logging
# =========================
append_log <- function(log_file, event) {
  line <- jsonlite::toJSON(event, auto_unbox = TRUE, null = "null")
  write(line, file = log_file, append = TRUE)
}

log_api_usage <- function(log_file, cluster, round, stage, prompt_str, resp) {
  usage <- resp$usage %||% NULL
  event <- list(
    ts = as.character(Sys.time()),
    cluster = cluster,
    round = round,
    stage = stage,
    ok = isTRUE(resp$ok),
    status = resp$status %||% NA_integer_,
    model = resp$model %||% NA_character_,
    prompt_chars = nchar(prompt_str %||% ""),
    resp_chars = nchar(resp$text %||% "")
  )
  if (!is.null(usage)) {
    event$prompt_tokens <- usage$prompt_tokens %||% NA_integer_
    event$completion_tokens <- usage$completion_tokens %||% NA_integer_
    event$total_tokens <- usage$total_tokens %||% NA_integer_
  }
  append_log(log_file, event)
}

# =========================
# 9.5) Open x TopK fusion (Step1 -> Step1.5)
# =========================
cl_graph_cache <- NULL
cl_graph_cache_path <- NULL

get_section <- function(obj, base) {
  if (is.null(obj) || !is.list(obj)) return(NULL)
  sec <- obj[[base]] %||% obj[[paste0(base, "_schema")]] %||% NULL
  if (!is.list(sec)) return(NULL)
  sec
}

get_cl_graph_cache <- function(cl_cfg) {
  path <- cl_cfg$local_json_path %||% ""
  if (!is.null(cl_graph_cache) && identical(cl_graph_cache_path, path)) return(cl_graph_cache)
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  cl <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  onehop_set <- new.env(parent = emptyenv())
  add_pair <- function(a, b) {
    if (is.null(a) || is.null(b)) return()
    key <- paste(sort(c(a, b)), collapse = "||")
    assign(key, TRUE, envir = onehop_set)
  }
  for (id in names(cl)) {
    if (!startsWith(id, "CL:")) next
    anc <- cl[[id]]$ancestors %||% NULL
    if (is.null(anc)) next
    d <- unlist(anc, use.names = TRUE)
    if (length(d) == 0) next
    for (p in names(d[d == 1])) add_pair(id, p)
  }
  cl_graph_cache <<- list(cl = cl, onehop_set = onehop_set)
  cl_graph_cache_path <<- path
  cl_graph_cache
}

is_valid_clid <- function(x, cl) {
  !is.na(x) && nzchar(x) && startsWith(x, "CL:") && !is.null(cl[[x]])
}

is_onehop <- function(a, b, onehop_set) {
  key <- paste(sort(c(a, b)), collapse = "||")
  exists(key, envir = onehop_set, inherits = FALSE)
}

dist_between <- function(a, b, cl) {
  if (!is_valid_clid(a, cl) || !is_valid_clid(b, cl)) return(Inf)
  if (a == b) return(0)
  anc_a <- cl[[a]]$ancestors %||% NULL
  anc_b <- cl[[b]]$ancestors %||% NULL
  if (!is.null(anc_a)) {
    d <- unlist(anc_a, use.names = TRUE)
    if (b %in% names(d)) return(as.numeric(d[[b]]))
  }
  if (!is.null(anc_b)) {
    d <- unlist(anc_b, use.names = TRUE)
    if (a %in% names(d)) return(as.numeric(d[[a]]))
  }
  Inf
}

normalize_core_labels <- function(step1_obj, cl_cfg) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(step1_obj)
  update_label <- function(sec, label_key, clid_key) {
    if (is.null(sec) || !is.list(sec)) return(sec)
    lbl <- as.character(sec[[label_key]] %||% NA_character_)[1]
    if (is.na(lbl) || !nzchar(lbl)) return(sec)
    norm <- Triage:::normalize_cl_three_state(lbl, "", cl_cfg)
    sec[[label_key]] <- norm$final_name %||% lbl
    sec[[clid_key]] <- norm$final_clid %||% sec[[clid_key]] %||% ""
    sec
  }

  main_key <- if (!is.null(step1_obj$main_type)) "main_type" else if (!is.null(step1_obj$main_type_schema)) "main_type_schema" else NULL
  sub1_key <- if (!is.null(step1_obj$subtype_level_1)) "subtype_level_1" else if (!is.null(step1_obj$subtype_level_1_schema)) "subtype_level_1_schema" else NULL
  sub2_key <- if (!is.null(step1_obj$subtype_level_2)) "subtype_level_2" else if (!is.null(step1_obj$subtype_level_2_schema)) "subtype_level_2_schema" else NULL

  if (!is.null(main_key)) {
    step1_obj[[main_key]] <- update_label(step1_obj[[main_key]], "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(sub1_key)) {
    step1_obj[[sub1_key]] <- update_label(step1_obj[[sub1_key]], "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(sub2_key) && is.list(step1_obj[[sub2_key]]$core_identity)) {
    step1_obj[[sub2_key]]$core_identity <- update_label(step1_obj[[sub2_key]]$core_identity, "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(step1_obj$open_world_summary)) {
    step1_obj$open_world_summary <- update_label(step1_obj$open_world_summary, "best_cell_type", "cell_ontology_id")
  }
  step1_obj
}

patch_fill_only_clid <- function(report, cl_cfg) {
  if (is.null(report) || !is.list(report) || is.null(cl_cfg)) return(report)
  fill_clid <- function(sec, label_key, clid_key) {
    if (is.null(sec) || !is.list(sec)) return(sec)
    if (nzchar(sec[[clid_key]] %||% "")) return(sec)
    lbl <- as.character(sec[[label_key]] %||% NA_character_)[1]
    if (!nzchar(lbl)) return(sec)
    m <- Triage:::normalize_cl_three_state(lbl, "", cl_cfg)
    sec[[clid_key]] <- as.character(m$final_clid %||% "")
    sec
  }

  main_key <- if (!is.null(report$main_type)) "main_type" else if (!is.null(report$main_type_schema)) "main_type_schema" else NULL
  sub1_key <- if (!is.null(report$subtype_level_1)) "subtype_level_1" else if (!is.null(report$subtype_level_1_schema)) "subtype_level_1_schema" else NULL
  sub2_key <- if (!is.null(report$subtype_level_2)) "subtype_level_2" else if (!is.null(report$subtype_level_2_schema)) "subtype_level_2_schema" else NULL

  if (!is.null(main_key)) {
    report[[main_key]] <- fill_clid(report[[main_key]], "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(sub1_key)) {
    report[[sub1_key]] <- fill_clid(report[[sub1_key]], "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(sub2_key) && is.list(report[[sub2_key]]$core_identity)) {
    report[[sub2_key]]$core_identity <- fill_clid(report[[sub2_key]]$core_identity, "candidate_cell_type", "cell_ontology_id")
  }
  if (!is.null(report$open_world_summary)) {
    report$open_world_summary <- fill_clid(report$open_world_summary, "best_cell_type", "cell_ontology_id")
  }
  report
}

patch_repair_schema_contradiction <- function(report, cl_cfg, top15_genes = NULL) {
  if (is.null(report) || !is.list(report)) return(report)
  ow_genes <- report$open_world_summary$evidence_genes %||% character(0)
  if (is.list(ow_genes)) ow_genes <- unlist(ow_genes, recursive = TRUE, use.names = FALSE)
  ow_genes <- as.character(ow_genes)
  if (!is.null(top15_genes)) {
    top15_genes <- as.character(top15_genes)
    if (sum(ow_genes %in% top15_genes, na.rm = TRUE) < 2) return(report)
  }
  ow <- tolower(as.character(report$open_world_summary$best_cell_type %||% ""))
  if (!nzchar(ow)) return(report)
  get_lbl <- function(x) tolower(as.character(x %||% ""))
  main_lbl <- get_lbl(report$main_type_schema$candidate_cell_type %||% report$main_type$candidate_cell_type)
  sub1_lbl <- get_lbl(report$subtype_level_1_schema$candidate_cell_type %||% report$subtype_level_1$candidate_cell_type)
  sub2_lbl <- get_lbl(report$subtype_level_2_schema$core_identity$candidate_cell_type %||% report$subtype_level_2$core_identity$candidate_cell_type)
  conflict <- (nzchar(main_lbl) && main_lbl != ow) ||
    (nzchar(sub1_lbl) && sub1_lbl != ow) ||
    (nzchar(sub2_lbl) && sub2_lbl != ow)
  ca <- report$main_type_schema$confidence_score_breakdown$candidate_agreement %||%
    report$main_type$confidence_score_breakdown$candidate_agreement %||% NA_integer_
  if (!isTRUE(conflict) || !identical(as.integer(ca), 0L)) return(report)

  report$main_type_schema$candidate_cell_type <- report$open_world_summary$best_cell_type
  report$subtype_level_1_schema$candidate_cell_type <- report$open_world_summary$best_cell_type
  report$subtype_level_2_schema$core_identity$candidate_cell_type <- report$open_world_summary$best_cell_type

  ow_clid <- report$open_world_summary$cell_ontology_id %||% ""
  if (nzchar(ow_clid)) {
    report$main_type_schema$cell_ontology_id <- ow_clid
    report$subtype_level_1_schema$cell_ontology_id <- ow_clid
    report$subtype_level_2_schema$core_identity$cell_ontology_id <- ow_clid
  } else if (!is.null(cl_cfg)) {
    m <- Triage:::normalize_cl_three_state(report$open_world_summary$best_cell_type, "", cl_cfg)
    clid <- as.character(m$final_clid %||% "")
    if (nzchar(clid)) {
      report$main_type_schema$cell_ontology_id <- clid
      report$subtype_level_1_schema$cell_ontology_id <- clid
      report$subtype_level_2_schema$core_identity$cell_ontology_id <- clid
    }
  }
  report
}

patch_program_triggered_refine <- function(report, top15_genes = NULL, cl_cfg = NULL, allowed_lineages = NULL) {
  if (is.null(report) || !is.list(report)) return(report)
  if (is.null(top15_genes) || length(top15_genes) == 0) return(report)
  allowed_effective <- setdiff(allowed_lineages %||% character(0), "unknown")
  if (length(allowed_effective) > 0 && !("epithelial" %in% allowed_effective)) return(report)
  program <- report$program_evidence %||% report$cluster_dossier$program_evidence %||% NULL
  if (is.null(program)) return(report)
  terms <- program$programs_from_terms %||% list()
  term_text <- tolower(paste(unlist(terms, recursive = TRUE, use.names = FALSE), collapse = " "))
  pigment_program <- stringr::str_detect(term_text, "pigment|melanin|melanogenesis|melanosome")
  if (!isTRUE(pigment_program)) return(report)
  support <- program$support_level %||% list()
  epi <- tolower(as.character(support$epithelial %||% ""))
  transport <- tolower(as.character(support$transport %||% ""))
  secretory <- tolower(as.character(support$secretory %||% ""))

  has_epi <- epi %in% c("moderate", "strong")
  has_transport <- transport %in% c("moderate", "strong")
  has_secretory <- secretory %in% c("moderate", "strong")
  if (!(has_epi && (has_transport || has_secretory))) return(report)

  # enforce pigmented epithelial lineage for subtype levels
  target_lbl <- "pigmented epithelial cell"
  report$subtype_level_1_schema$candidate_cell_type <- target_lbl
  report$subtype_level_2_schema$core_identity$candidate_cell_type <- target_lbl
  if (!is.null(cl_cfg)) {
    m <- Triage:::normalize_cl_three_state(target_lbl, "", cl_cfg)
    clid <- as.character(m$final_clid %||% "")
    if (nzchar(clid)) {
      report$subtype_level_1_schema$cell_ontology_id <- clid
      report$subtype_level_2_schema$core_identity$cell_ontology_id <- clid
    }
  }
  report
}

infer_coarse_lineage_simple <- function(label, clid = NULL, cl_cfg = NULL) {
  lineage_from_name <- function(nm) {
    x <- tolower(as.character(nm %||% ""))
    if (!nzchar(x)) return(NA_character_)
    if (stringr::str_detect(x, "epithel")) return("epithelial")
    if (stringr::str_detect(x, "endothel")) return("endothelial")
    if (stringr::str_detect(x, "fibroblast|stromal|mesenchym")) return("stromal_mesenchymal")
    if (stringr::str_detect(x, "smooth muscle|myocyte|muscle")) return("muscle")
  if (stringr::str_detect(x, "glia|astro|olig|microglia|ependym|neur|neuron|glutamatergic|gaba|interneuron")) return("neural_glial")
    if (stringr::str_detect(x, "immune|lymph|t cell|b cell|nk|myeloid|mono|macro|dendritic")) return("immune")
    if (stringr::str_detect(x, "eryth|platelet|megakaryo")) return("erythroid_megakaryocytic")
    NA_character_
  }

  if (!is.null(cl_cfg) && !is.null(clid) && nzchar(clid)) {
    cache <- get_cl_graph_cache(cl_cfg)
    cl <- if (!is.null(cache)) cache$cl else NULL
    if (!is.null(cl) && !is.null(cl[[clid]])) {
      ids <- c(clid)
      anc <- cl[[clid]]$ancestors %||% NULL
      if (!is.null(anc)) {
        ids <- unique(c(ids, names(unlist(anc, use.names = TRUE))))
      }
      for (id in ids) {
        nm <- cl[[id]]$name %||% ""
        lin <- lineage_from_name(nm)
        if (!is.na(lin)) return(lin)
      }
    }
  }

  lin <- lineage_from_name(label)
  if (!is.na(lin)) return(lin)
  "unknown"
}

patch_repair_cross_lineage <- function(report) {
  if (is.null(report) || !is.list(report)) return(report)
  main_lbl <- report$main_type_schema$candidate_cell_type %||% report$main_type$candidate_cell_type %||% ""
  sub1_lbl <- report$subtype_level_1_schema$candidate_cell_type %||% report$subtype_level_1$candidate_cell_type %||% ""
  sub2_lbl <- report$subtype_level_2_schema$core_identity$candidate_cell_type %||% report$subtype_level_2$core_identity$candidate_cell_type %||% ""
  ow_lbl <- report$open_world_summary$best_cell_type %||% ""

  main_lin <- infer_coarse_lineage_simple(main_lbl, report$main_type_schema$cell_ontology_id %||% report$main_type$cell_ontology_id, cl_cfg)
  sub1_lin <- infer_coarse_lineage_simple(sub1_lbl, report$subtype_level_1_schema$cell_ontology_id %||% report$subtype_level_1$cell_ontology_id, cl_cfg)
  sub2_lin <- infer_coarse_lineage_simple(sub2_lbl, report$subtype_level_2_schema$core_identity$cell_ontology_id %||% report$subtype_level_2$core_identity$cell_ontology_id, cl_cfg)
  ow_lin <- infer_coarse_lineage_simple(ow_lbl, report$open_world_summary$cell_ontology_id, cl_cfg)

  if (main_lin == "unknown") return(report)

  if (ow_lin != "unknown" && ow_lin != main_lin && nzchar(main_lbl)) {
    report$open_world_summary$best_cell_type <- main_lbl
    if (!is.null(cl_cfg)) {
      m <- Triage:::normalize_cl_three_state(main_lbl, "", cl_cfg)
      clid <- as.character(m$final_clid %||% "")
      if (nzchar(clid)) report$open_world_summary$cell_ontology_id <- clid
    }
  }

  is_stage_only_label <- function(x) {
    grepl("(progenitor|precursor|intermediate|blast|cycling|proliferat)",
          tolower(as.character(x %||% "")))
  }

  effective_sub1_lin <- sub1_lin
  effective_sub2_lin <- sub2_lin
  if (sub1_lin == "unknown" && is_stage_only_label(sub1_lbl)) effective_sub1_lin <- main_lin
  if (sub2_lin == "unknown" && is_stage_only_label(sub2_lbl)) effective_sub2_lin <- main_lin

  conflict <- (effective_sub1_lin != "unknown" && effective_sub1_lin != main_lin) ||
    (effective_sub2_lin != "unknown" && effective_sub2_lin != main_lin)
  if (!isTRUE(conflict)) return(report)

  target_lbl <- if (ow_lin == main_lin && nzchar(ow_lbl)) ow_lbl else main_lbl

  report$subtype_level_1_schema$candidate_cell_type <- target_lbl
  report$subtype_level_2_schema$core_identity$candidate_cell_type <- target_lbl

  if (nzchar(target_lbl) && !is.null(cl_cfg)) {
    m <- Triage:::normalize_cl_three_state(target_lbl, "", cl_cfg)
    clid <- as.character(m$final_clid %||% "")
    if (nzchar(clid)) {
      report$subtype_level_1_schema$cell_ontology_id <- clid
      report$subtype_level_2_schema$core_identity$cell_ontology_id <- clid
    }
  }
  report
}


sync_schema_to_core <- function(report) {
  if (is.null(report) || !is.list(report)) return(report)

  if (!is.null(report$main_type_schema)) {
    if (is.null(report$main_type)) report$main_type <- list()
    if (nzchar(report$main_type_schema$candidate_cell_type %||% "")) {
      report$main_type$candidate_cell_type <- report$main_type_schema$candidate_cell_type
    }
    if (nzchar(report$main_type_schema$cell_ontology_id %||% "")) {
      report$main_type$cell_ontology_id <- report$main_type_schema$cell_ontology_id
    }
  }

  if (!is.null(report$subtype_level_1_schema)) {
    if (is.null(report$subtype_level_1)) report$subtype_level_1 <- list()
    if (nzchar(report$subtype_level_1_schema$candidate_cell_type %||% "")) {
      report$subtype_level_1$candidate_cell_type <- report$subtype_level_1_schema$candidate_cell_type
    }
    if (nzchar(report$subtype_level_1_schema$cell_ontology_id %||% "")) {
      report$subtype_level_1$cell_ontology_id <- report$subtype_level_1_schema$cell_ontology_id
    }
  }

  if (!is.null(report$subtype_level_2_schema$core_identity)) {
    if (is.null(report$subtype_level_2)) report$subtype_level_2 <- list()
    if (is.null(report$subtype_level_2$core_identity)) report$subtype_level_2$core_identity <- list()
    if (nzchar(report$subtype_level_2_schema$core_identity$candidate_cell_type %||% "")) {
      report$subtype_level_2$core_identity$candidate_cell_type <- report$subtype_level_2_schema$core_identity$candidate_cell_type
    }
    if (nzchar(report$subtype_level_2_schema$core_identity$cell_ontology_id %||% "")) {
      report$subtype_level_2$core_identity$cell_ontology_id <- report$subtype_level_2_schema$core_identity$cell_ontology_id
    }
  }

  report
}

apply_openxtopk_fusion <- function(step1_obj, cl_cfg, log_file = NULL, prefix = NULL, allowed_lineages = NULL) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(step1_obj)
  st <- as.character(step1_obj$evaluation_status %||% NA_character_)
  if (!is.na(st) && stringr::str_detect(st, "^Rejected")) return(step1_obj)

  main_sec <- get_section(step1_obj, "main_type")
  sub1_sec <- get_section(step1_obj, "subtype_level_1")
  sub2_sec <- get_section(step1_obj, "subtype_level_2")
  if (is.null(main_sec)) return(step1_obj)

  open_label <- as.character(step1_obj$open_world_summary$best_cell_type %||% NA_character_)[1]
  open_label <- stringr::str_trim(open_label)
  if (!nzchar(open_label)) return(step1_obj)

  topk_labels <- c(
    as.character(main_sec$candidate_cell_type %||% NA_character_)[1],
    as.character(sub1_sec$candidate_cell_type %||% NA_character_)[1],
    as.character(sub2_sec$core_identity$candidate_cell_type %||% NA_character_)[1]
  )
  topk_labels <- stringr::str_trim(topk_labels)
  topk_labels[topk_labels == ""] <- NA_character_

  topk_clids <- c(
    as.character(main_sec$cell_ontology_id %||% NA_character_)[1],
    as.character(sub1_sec$cell_ontology_id %||% NA_character_)[1],
    as.character(sub2_sec$core_identity$cell_ontology_id %||% NA_character_)[1]
  )

  open_norm <- Triage:::norm_name(open_label)
  topk_norms <- vapply(topk_labels, Triage:::norm_name, character(1))

  open_map <- Triage:::normalize_cl_three_state(open_label, "", cl_cfg)
  open_clid <- as.character(open_map$final_clid %||% NA_character_)[1]

  allowed_effective <- setdiff(allowed_lineages %||% character(0), "unknown")
  if (length(allowed_effective) > 0) {
    open_lin <- infer_coarse_lineage_simple(open_label, open_clid, cl_cfg)
    if (open_lin == "unknown" || !(open_lin %in% allowed_effective)) return(step1_obj)
  }

  cache <- get_cl_graph_cache(cl_cfg)
  cl <- if (!is.null(cache)) cache$cl else NULL
  onehop_set <- if (!is.null(cache)) cache$onehop_set else NULL

  valid_clid_format <- function(x) !is.na(x) && nzchar(x) && startsWith(x, "CL:")
  if (!is.null(cl) && !is_valid_clid(open_clid, cl)) open_clid <- NA_character_
  if (is.null(cl) && !valid_clid_format(open_clid)) open_clid <- NA_character_

  for (i in seq_along(topk_clids)) {
    if (is.na(topk_clids[i]) || !nzchar(topk_clids[i])) {
      if (!is.na(topk_labels[i]) && nzchar(topk_labels[i])) {
        m <- Triage:::normalize_cl_three_state(topk_labels[i], "", cl_cfg)
        topk_clids[i] <- as.character(m$final_clid %||% NA_character_)[1]
      }
    }
  }

  match_idx <- which(!is.na(topk_norms) & nzchar(topk_norms) & (open_norm == topk_norms))
  if (length(match_idx) == 0 && nzchar(open_clid)) {
    match_idx <- which(topk_clids == open_clid)
  }

  chosen_label <- NULL
  chosen_clid <- NA_character_
  parent_label <- NULL
  child_label <- NULL
  parent_clid <- NA_character_
  child_clid <- NA_character_
  reason <- NULL

  if (length(match_idx) > 0) {
    chosen_label <- open_label
    chosen_clid <- if (nzchar(open_clid)) open_clid else topk_clids[match_idx[[1]]]
    reason <- "exact_match"
  } else if (!is.null(cl) && nzchar(open_clid)) {
    for (i in seq_along(topk_clids)) {
      tc <- topk_clids[i]
      if (!is_valid_clid(tc, cl)) next
      if (!is_onehop(open_clid, tc, onehop_set)) next
      d_open_to_tc <- dist_between(open_clid, tc, cl)
      d_tc_to_open <- dist_between(tc, open_clid, cl)
      if (is.finite(d_open_to_tc) && d_open_to_tc == 1) {
        parent_clid <- tc
        child_clid <- open_clid
        parent_label <- Triage:::local_lookup_by_clid(parent_clid, cl_cfg) %||% topk_labels[i]
        child_label <- open_label
        chosen_label <- parent_label
        chosen_clid <- parent_clid
        reason <- "onehop_open_child"
      } else if (is.finite(d_tc_to_open) && d_tc_to_open == 1) {
        parent_clid <- open_clid
        child_clid <- tc
        parent_label <- open_label
        child_label <- Triage:::local_lookup_by_clid(child_clid, cl_cfg) %||% topk_labels[i]
        chosen_label <- parent_label
        chosen_clid <- parent_clid
        reason <- "onehop_topk_child"
      }
      if (!is.null(chosen_label)) break
    }
  }

  if (!is.null(chosen_label) && nzchar(chosen_label)) {
    main_key <- if (!is.null(step1_obj$main_type)) "main_type" else if (!is.null(step1_obj$main_type_schema)) "main_type_schema" else NULL
    sub1_key <- if (!is.null(step1_obj$subtype_level_1)) "subtype_level_1" else if (!is.null(step1_obj$subtype_level_1_schema)) "subtype_level_1_schema" else NULL
    sub2_key <- if (!is.null(step1_obj$subtype_level_2)) "subtype_level_2" else if (!is.null(step1_obj$subtype_level_2_schema)) "subtype_level_2_schema" else NULL

    if (!is.null(main_key)) {
      norm_main <- Triage:::normalize_cl_three_state(chosen_label, "", cl_cfg)
      step1_obj[[main_key]]$candidate_cell_type <- norm_main$final_name %||% chosen_label
      step1_obj[[main_key]]$cell_ontology_id <- norm_main$final_clid %||% ""
    }

    if (!is.null(child_label) && nzchar(child_label)) {
      if (!is.null(sub1_key)) {
        norm_sub1 <- Triage:::normalize_cl_three_state(child_label, "", cl_cfg)
        step1_obj[[sub1_key]]$candidate_cell_type <- norm_sub1$final_name %||% child_label
        step1_obj[[sub1_key]]$cell_ontology_id <- norm_sub1$final_clid %||% ""
      }
      if (!is.null(sub2_key) && is.list(step1_obj[[sub2_key]]$core_identity)) {
        norm_sub2 <- Triage:::normalize_cl_three_state(child_label, "", cl_cfg)
        step1_obj[[sub2_key]]$core_identity$candidate_cell_type <- norm_sub2$final_name %||% child_label
        step1_obj[[sub2_key]]$core_identity$cell_ontology_id <- norm_sub2$final_clid %||% ""
      }
    }

    if (!is.null(log_file)) {
      append_log(log_file, list(
        ts = as.character(Sys.time()),
        cluster = prefix %||% NA_character_,
        stage = "step1_openxtopk_fusion",
        ok = TRUE,
        reason = reason,
        open_label = open_label,
        chosen_label = step1_obj[[main_key]]$candidate_cell_type %||% NA_character_,
        chosen_clid = step1_obj[[main_key]]$cell_ontology_id %||% NA_character_,
        subtype_label = if (!is.null(sub1_key)) step1_obj[[sub1_key]]$candidate_cell_type %||% NA_character_ else NA_character_
      ))
    }
  }

  step1_obj
}

# =========================
# 10) Worker
# =========================
process_cluster <- function(step1_file, config) {
  input_step1_dir  <- file.path(config$input_root, config$in_step1)
  input_step15_dir <- file.path(config$input_root, config$in_step15)
  input_step2_dir  <- file.path(config$input_root, config$in_step2)
  
  output_step1_dir  <- file.path(config$output_root, config$out_step1)
  output_step15_dir <- file.path(config$output_root, config$out_step15)
  output_step2_dir  <- file.path(config$output_root, config$out_step2)
  output_final_dir  <- file.path(config$output_root, config$out_final)
  output_debug_dir  <- file.path(config$output_root, config$out_debug)
  output_logs_dir   <- file.path(config$output_root, config$out_logs)
  
  ensure_dir(output_step1_dir); ensure_dir(output_step15_dir); ensure_dir(output_step2_dir)
  ensure_dir(output_final_dir); ensure_dir(output_debug_dir); ensure_dir(output_logs_dir)
  
  prefix <- file_prefix_from_step1(step1_file)
  cat(glue("[CLUSTER] {prefix}\n"))
  
  log_file <- file.path(output_logs_dir, paste0(prefix, "_runlog.jsonl"))
  append_log(log_file, list(ts=as.character(Sys.time()), cluster=prefix, stage="cluster_start", ok=TRUE))
  validator_feedback <- NULL
  
  step15_file <- file.path(input_step15_dir, paste0(prefix, "_step1.5_citation_fix_query.json"))
  step2_file  <- file.path(input_step2_dir,  paste0(prefix, "_step2_validation_query.json"))
  
  if (!file.exists(step15_file)) {
    append_log(log_file, list(ts=as.character(Sys.time()), cluster=prefix, stage="missing_step15", ok=FALSE))
    return(list(cluster=prefix, final_pass=FALSE, fail_stage="missing_step15"))
  }
  
  step1_query_obj  <- read_json_safely(step1_file)
  step15_query_tpl <- read_json_safely(step15_file)
  
  # --- Step2 instructions source:
  # If step2_file exists, prefer it; else auto-generate CASSIA-style instructions.
  step2_instructions <- NULL
  if (isTRUE(config$enable_step2)) {
    if (file.exists(step2_file)) {
      step2_query_tpl <- read_json_safely(step2_file)
      step2_instructions <- step2_query_tpl$instructions_for_llm %||% build_validator_instructions_cassia()
      cat(glue("[INFO] Step2 query file found: {basename(step2_file)}\n"))
    } else {
      step2_instructions <- build_validator_instructions_cassia()
      cat(glue("[INFO] Step2 query file missing; auto-generating Step2 instructions for {prefix}\n"))
    }
  }
  
  allowed_pmids <- extract_allowed_pmids(step15_query_tpl, step1_query_obj)
  for (round in seq_len(config$max_rounds)) {
    is_final_round <- round == config$max_rounds
    round_guardrail <- if (is_final_round) {
      "FINAL ROUND: Mixed/ambiguous/unknown/other labels or rejection schemas are allowed only if no plausible lineage remains after in_scope_top_genes + ranks_16_50_in_scope review. Otherwise choose the best-supported lineage and express uncertainty only in phenotypic_label."
    } else {
      "NOT FINAL ROUND: Do NOT output mixed/ambiguous/unknown/other labels or rejection schemas. Choose the best-supported lineage and express uncertainty only in phenotypic_label."
    }
    
    # ---------- Step1 prompt (inject feedback if any) ----------
    step1_prompt_obj <- step1_query_obj
    if (!is.null(validator_feedback)) {
      # Minimal: append feedback to instructions (works whether instructions_for_llm is list or string)
      if (is.list(step1_prompt_obj$instructions_for_llm)) {
        step1_prompt_obj$instructions_for_llm$validator_feedback <- validator_feedback
      } else if (is.character(step1_prompt_obj$instructions_for_llm)) {
        step1_prompt_obj$instructions_for_llm <- paste(step1_prompt_obj$instructions_for_llm,
                                                       "\n\n[VALIDATOR_FEEDBACK_FROM_PREV_ROUND]\n",
                                                       validator_feedback,
                                                       sep = "")
      }
      if (is.null(step1_prompt_obj$input_data)) step1_prompt_obj$input_data <- list()
      step1_prompt_obj$input_data$validator_feedback <- validator_feedback
    }
    step1_prompt_obj <- append_instruction_text(step1_prompt_obj, round_guardrail)
    
    if (isTRUE(config$prune_step1)) {
      step1_prompt_obj <- prune_step1_evidence(step1_prompt_obj, config$keep_per_dimension_step1, config$keep_chars_step1)
    }
    
    step1_prompt_str <- jsonlite::toJSON(step1_prompt_obj, auto_unbox = TRUE, null = "null")
    resp1 <- invoke_deepseek_api(step1_prompt_str, config, model = config$model_step1)
    log_api_usage(log_file, prefix, round, "Step1", step1_prompt_str, resp1)
    if (!isTRUE(resp1$ok)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step1_api_round", round)))
    
    step1_clean <- extract_first_json_object_stack(resp1$text)
    if (is.na(step1_clean)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step1_nonjson_round", round)))
    step1_obj <- jsonlite::fromJSON(step1_clean, simplifyVector = FALSE)
    step1_obj <- normalize_core_labels(step1_obj, cl_cfg)
    step1_obj <- patch_fill_only_clid(step1_obj, cl_cfg)
    
    save_json_pretty(step1_obj, file.path(output_step1_dir, paste0(prefix, "_round", round, "_raw_step1_report_output.json")))
    
    degs_data <- step1_query_obj$input_data$cluster_dossier$degs$data %||% list()
    top15_genes <- degs_data$geneSymbol %||% degs_data$gene %||% character(0)
    top15_genes <- as.character(top15_genes)
    if (length(top15_genes) > 15) top15_genes <- top15_genes[1:15]
    
  step1_obj <- apply_openxtopk_fusion(step1_obj, cl_cfg, log_file = log_file, prefix = prefix, allowed_lineages = prompt_context_config$allowed_lineages)
    step1_obj <- patch_repair_schema_contradiction(step1_obj, cl_cfg, top15_genes = top15_genes)
    step1_obj <- patch_repair_cross_lineage(step1_obj)
    step1_obj <- sync_schema_to_core(step1_obj)
    step1_obj$program_evidence <- step1_query_obj$input_data$cluster_dossier$program_evidence %||% NULL
  step1_obj <- patch_program_triggered_refine(step1_obj, top15_genes = top15_genes, cl_cfg = cl_cfg, allowed_lineages = prompt_context_config$allowed_lineages)
    step1_obj <- enforce_scope_coarse_gate(step1_obj, prompt_context_config)
    
    save_json_pretty(step1_obj, file.path(output_step1_dir, paste0(prefix, "_round", round, "_step1_report_output.json")))
    
    # ---------- Step1.5 ----------
    step15_prompt_obj <- inject_step15_report_and_guardrails(step15_query_tpl, step1_obj, allowed_pmids)
    step15_prompt_obj <- append_instruction_text(step15_prompt_obj, round_guardrail)
    if (isTRUE(config$prune_step15)) {
      tmp <- list(input_data = list(cluster_dossier = list(evidence = step15_prompt_obj$input_data$literature_evidence)))
      tmp2 <- prune_step1_evidence(tmp, keep_per_dim = config$keep_per_dimension_step15, keep_chars = config$keep_chars_step15)
      step15_prompt_obj$input_data$literature_evidence <- tmp2$input_data$cluster_dossier$evidence
    }
    
    step15_prompt_str <- jsonlite::toJSON(step15_prompt_obj, auto_unbox = TRUE, null = "null")
    resp15 <- invoke_deepseek_api(step15_prompt_str, config, model = config$model_step15)
    log_api_usage(log_file, prefix, round, "Step1.5", step15_prompt_str, resp15)
    if (!isTRUE(resp15$ok)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step15_api_round", round)))
    
    step15_clean <- extract_first_json_object_stack(resp15$text)
    if (is.na(step15_clean)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step15_nonjson_round", round)))
    step15_obj <- jsonlite::fromJSON(step15_clean, simplifyVector = FALSE)
    step15_obj <- Triage:::normalize_step15_cl(step15_obj, cl_cfg)
    step15_obj <- normalize_core_labels(step15_obj, cl_cfg)
    step15_obj <- patch_fill_only_clid(step15_obj, cl_cfg)
  step15_obj <- apply_openxtopk_fusion(step15_obj, cl_cfg, log_file = log_file, prefix = prefix, allowed_lineages = prompt_context_config$allowed_lineages)
    step15_obj <- patch_repair_schema_contradiction(step15_obj, cl_cfg, top15_genes = top15_genes)
    step15_obj <- patch_repair_cross_lineage(step15_obj)
    step15_obj <- sync_schema_to_core(step15_obj)
    step15_obj$program_evidence <- step1_obj$program_evidence %||% NULL
  step15_obj <- patch_program_triggered_refine(step15_obj, top15_genes = top15_genes, cl_cfg = cl_cfg, allowed_lineages = prompt_context_config$allowed_lineages)
    step15_obj <- enforce_scope_coarse_gate(step15_obj, prompt_context_config)
    
    save_json_pretty(step15_obj, file.path(output_step15_dir, paste0(prefix, "_round", round, "_step1.5_citation_fix_output.json")))
    
    # ---------- Gate ----------
    gate <- citation_gate(step15_obj, allowed_pmids, overuse_threshold = config$gate_pmid_overuse_threshold)
    save_json_pretty(gate, file.path(output_step15_dir, paste0(prefix, "_round", round, "_GATE.json")))
    
    # If Step2 disabled: final = gate passed
    if (!isTRUE(config$enable_step2)) {
      if (isTRUE(gate$gate_passed)) {
        save_json_pretty(step15_obj, file.path(output_final_dir, paste0(prefix, "_FINAL_passed.json")))
        return(list(cluster=prefix, final_pass=TRUE, fail_stage=NA_character_))
      } else {
        save_json_pretty(list(gate=gate, step15=step15_obj),
                         file.path(output_debug_dir, paste0(prefix, "_FINAL_failed_debug.json")))
        return(list(cluster=prefix, final_pass=FALSE, fail_stage="gate_failed"))
      }
    }
    
    # ---------- Step2 ----------
    step2_prompt_obj <- build_step2_runtime_prompt(
      step2_instructions = step2_instructions,
      step15_output_obj = step15_obj,
      step1_query_obj = step1_query_obj,
      query_id = paste0(prefix, "_validation_round", round)
    )
    step2_prompt_obj <- append_instruction_text(step2_prompt_obj, round_guardrail)
    
    if (isTRUE(config$prune_step2)) {
      od <- step2_prompt_obj$input_data$original_dossier %||% list()
      ev <- od$evidence %||% NULL
      tmp <- list(input_data = list(cluster_dossier = list(evidence = ev)))
      tmp2 <- prune_step1_evidence(tmp, keep_per_dim = config$keep_per_dimension_step2, keep_chars = config$keep_chars_step2)
      if (!is.null(tmp2$input_data$cluster_dossier$evidence)) {
        step2_prompt_obj$input_data$original_dossier$evidence <- tmp2$input_data$cluster_dossier$evidence
      }
    }
    
    step2_prompt_str <- jsonlite::toJSON(step2_prompt_obj, auto_unbox = TRUE, null = "null")
    resp2 <- invoke_deepseek_api(step2_prompt_str, config, model = config$model_step2)
    log_api_usage(log_file, prefix, round, "Step2", step2_prompt_str, resp2)
    if (!isTRUE(resp2$ok)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step2_api_round", round)))
    
    step2_clean <- extract_first_json_object_stack(resp2$text)
    if (is.na(step2_clean)) return(list(cluster=prefix, final_pass=FALSE, fail_stage=paste0("step2_nonjson_round", round)))
    step2_obj <- jsonlite::fromJSON(step2_clean, simplifyVector = FALSE)
    save_json_pretty(step2_obj, file.path(output_step2_dir, paste0(prefix, "_round", round, "_step2_validation_output.json")))
    
    status <- get_validation_status(step2_obj)
    llm_pass <- (!is.na(status) && grepl("PASSED", toupper(status)))
    final_pass <- (llm_pass && isTRUE(gate$gate_passed))

    revision <- step2_obj$revision_directive %||% list()
    action <- tolower(revision$action %||% "")
    suggested_parent <- as.character(revision$suggested_parent %||% "")
    suggested_parent_clid <- as.character(revision$suggested_parent_clid %||% "")
    allowed_effective <- setdiff(prompt_context_config$allowed_lineages %||% character(0), "unknown")
    if (length(allowed_effective) > 0 && nzchar(suggested_parent)) {
      lin <- infer_coarse_lineage_simple(suggested_parent, suggested_parent_clid, cl_cfg)
      if (lin == "unknown" || !(lin %in% allowed_effective)) {
        suggested_parent <- pick_conservative_parent_from_summary(step1_query_obj, prompt_context_config)
        suggested_parent_clid <- ""
      }
    }

    if (!isTRUE(final_pass) && isTRUE(gate$gate_passed) && action == "downgrade_to_parent" && nzchar(suggested_parent)) {
      step15_obj <- apply_core_downgrade(step15_obj, suggested_parent, suggested_parent_clid)
      step15_obj$auto_downgrade_note <- list(
        reason = "B_overspecific",
        suggested_parent = suggested_parent,
        suggested_parent_clid = suggested_parent_clid,
        source = "step2"
      )
      save_json_pretty(step15_obj, file.path(output_final_dir, paste0(prefix, "_FINAL_passed.json")))
      save_json_pretty(list(gate = gate, step2 = step2_obj, step15 = step15_obj, auto_downgraded = TRUE),
                       file.path(output_debug_dir, paste0(prefix, "_FINAL_downgraded_debug.json")))
      return(list(cluster = prefix, final_pass = TRUE, fail_stage = NA_character_))
    }

    in_scope_primary <- step1_query_obj$input_data$cluster_dossier$degs$in_scope_top_genes %||% character(0)
    in_scope_recovery <- step1_query_obj$input_data$cluster_dossier$degs$ranks_16_50_in_scope %||% character(0)
    has_in_scope <- length(in_scope_primary) > 0 || length(in_scope_recovery) > 0

    if (!isTRUE(final_pass) && isTRUE(gate$gate_passed) && has_in_scope) {
      conservative_parent <- pick_conservative_parent_from_summary(step1_query_obj, prompt_context_config)
      step15_obj <- apply_core_downgrade(step15_obj, conservative_parent, "")
      step15_obj$auto_downgrade_note <- list(
        reason = "weak_in_scope",
        suggested_parent = conservative_parent,
        suggested_parent_clid = "",
        source = "step2"
      )
      step15_obj$subtype_level_2$phenotypic_label <- paste(
        step15_obj$subtype_level_2$phenotypic_label %||% "",
        "| NOTE: weak in-scope evidence; conservative parent applied.",
        sep = " "
      )
      save_json_pretty(step15_obj, file.path(output_final_dir, paste0(prefix, "_FINAL_passed.json")))
      save_json_pretty(list(gate = gate, step2 = step2_obj, step15 = step15_obj, auto_downgraded = TRUE),
                       file.path(output_debug_dir, paste0(prefix, "_FINAL_downgraded_debug.json")))
      return(list(cluster = prefix, final_pass = TRUE, fail_stage = NA_character_))
    }

    if (!isTRUE(final_pass) && isTRUE(gate$gate_passed) && !has_in_scope) {
      conservative_parent <- pick_conservative_parent_from_summary(step1_query_obj, prompt_context_config)
      step15_obj <- apply_core_downgrade(step15_obj, conservative_parent, "")
      step15_obj$auto_downgrade_note <- list(
        reason = "abstain_no_in_scope",
        suggested_parent = conservative_parent,
        suggested_parent_clid = "",
        source = "step2"
      )
      step15_obj$subtype_level_2$phenotypic_label <- paste(
        step15_obj$subtype_level_2$phenotypic_label %||% "",
        "| NOTE: abstain due to no in-scope evidence.",
        sep = " "
      )
      save_json_pretty(step15_obj, file.path(output_final_dir, paste0(prefix, "_FINAL_passed.json")))
      save_json_pretty(list(gate = gate, step2 = step2_obj, step15 = step15_obj, auto_downgraded = TRUE),
                       file.path(output_debug_dir, paste0(prefix, "_FINAL_downgraded_debug.json")))
      return(list(cluster = prefix, final_pass = TRUE, fail_stage = NA_character_))
    }
    
    if (isTRUE(final_pass)) {
      save_json_pretty(step15_obj, file.path(output_final_dir, paste0(prefix, "_FINAL_passed.json")))
      return(list(cluster=prefix, final_pass=TRUE, fail_stage=NA_character_))
    }
    
    # ---------- prepare feedback for next round ----------
    failure_reasons <- step2_obj$failure_reasons %||% character(0)
    suggested_fixes <- step2_obj$suggested_fixes %||% character(0)
    core_sum <- step2_obj$core_identity_verdict$core_identity_summary %||% ""
    
    validator_feedback <- paste(
      "Prev round failed Step2 validation. You MUST address the following:",
      if (length(failure_reasons)) paste0("- failure_reasons: ", paste(failure_reasons, collapse=" | ")) else "- failure_reasons: (none provided)",
      if (nzchar(core_sum)) paste0("- core_identity_summary: ", core_sum) else "- core_identity_summary: (empty)",
      if (length(suggested_fixes)) paste0("- suggested_fixes: ", paste(suggested_fixes, collapse=" | ")) else "- suggested_fixes: (none provided)",
      "Do NOT invent markers not present in DEGs/program_evidence. If forced to downgrade, keep organ-specific guess only in phenotypic_label with '-like' wording.",
      sep = "\n"
    )
    
    # last round -> dump debug and return
    if (round == config$max_rounds) {
      save_json_pretty(list(gate=gate, step2=step2_obj, step15=step15_obj, last_validator_feedback=validator_feedback),
                       file.path(output_debug_dir, paste0(prefix, "_FINAL_failed_debug.json")))
      return(list(cluster=prefix, final_pass=FALSE, fail_stage="failed"))
    }
  }
}

# =========================
# 11) MAIN
# =========================
input_step1_dir  <- file.path(config$input_root, config$in_step1)
if (!dir.exists(input_step1_dir)) stop(glue("Missing Step1 input dir: {input_step1_dir}"))

step1_files <- list.files(input_step1_dir, pattern = "_step1_report_query\\.json$", full.names = TRUE, recursive = TRUE)
if (length(step1_files) == 0) stop(glue("No Step1 query files under: {input_step1_dir}"))
cat(glue("Found {length(step1_files)} Step1 query files.\n"))

# Dryrun: cap to first N
if (identical(mode, "dryrun")) {
  step1_files <- sort(step1_files)
  step1_files <- head(step1_files, max(1L, as.integer(opt$clusters)))
  cat(glue("Dryrun: using {length(step1_files)} Step1 files.\n"))
}

# Single cluster test mode
if (isTRUE(config$test_one_cluster)) {
  step1_files <- step1_files[grepl(paste0("/", config$test_cluster_prefix, "_step1_report_query\\.json$"), step1_files)]
  if (length(step1_files) == 0) stop("TEST MODE: no matching cluster file found.")
  cat(glue(">>> TEST MODE enabled. Running {length(step1_files)} cluster(s).\n"))
}

# Ensure output dirs
ensure_dir(file.path(config$output_root, config$out_step1))
ensure_dir(file.path(config$output_root, config$out_step15))
ensure_dir(file.path(config$output_root, config$out_step2))
ensure_dir(file.path(config$output_root, config$out_final))
ensure_dir(file.path(config$output_root, config$out_debug))
ensure_dir(file.path(config$output_root, config$out_logs))

  options(future.globals.maxSize = 8 * 1024^3)
  plan(multisession, workers = config$parallel_workers)

res <- furrr::future_map(step1_files, ~ process_cluster(.x, config), .options = furrr::furrr_options(seed = config$parallel_seed))
res_df <- dplyr::bind_rows(res)
save_json_pretty(list(results = res), file.path(config$output_root, "run_summary.json"))
readr::write_csv(res_df, file.path(config$output_root, "run_summary.csv"))
cat(glue::glue("Summary saved: {file.path(config$output_root, 'run_summary.csv')}\n"))
