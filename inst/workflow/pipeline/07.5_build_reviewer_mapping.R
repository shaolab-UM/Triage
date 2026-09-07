#!/usr/bin/env Rscript
# =============================================================
# 07.5_build_reviewer_mapping.R — CL-Linker: mapping for three reviewers with a fixed registry
#
# Workflow (single entry point):
#   Input: 03b CASSIA + 07 in-house summary + 07b enrichment summary (150 reviewer instances)
#         + step-1 queries (markers) + frozen CL ontology
#
#   For each reviewer instance:
# Step 1: (CL-Linker: canonical/synonym/normalized/alias/parenthetical/provided_id)
#       -> on success: retain directly (mapping_status=mapped, without an LLM call)
#       -> on failure: continue
#     Step 2: lexical top-k candidates (pure R, stringdist)
# Step 3: LLM Mapper A/B (with marker/context evidence, candidates)
#     Step 4: A/B stability + specificity consistency checks
# Step 5: risk-based Verifier (high-risk cases only, markers/context/MapperOutput)
#     Step 6: deterministic ontology validation (frozen CL graph)
#
#   Output (a single mapping registry):
#     data/primary/reviewer_mapping_registry.tsv
#
# Public release version: v1.0.0
# =============================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(jsonlite)
})
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required. Install it before running this script.")
}

option_list <- list(
  make_option("--cl_json", type = "character",
              default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json")),
  make_option("--out_dir", type = "character", default = "reports/cl_mapping"),
  make_option("--manifest", type = "character", default = "",
              help = "Optional run manifest TSV; default reproducibility/config/evidence_mapping_runs.tsv"),
  make_option("--model", type = "character", default = "deepseek-v4-flash"),
  make_option("--api_key_env", type = "character", default = "DEEPSEEK_API_KEY"),
  make_option("--temperature", type = "double", default = 0),
  make_option("--no_cache", action = "store_true", default = FALSE),
  make_option("--max_tokens", type = "integer", default = 300L),
  make_option("--dataset", type = "character", default = "", help = "Process only the specified dataset (subset run)"),
  make_option("--cluster_id", type = "character", default = "", help = "Process only the specified cluster (subset run)"),
  make_option("--reviewer", type = "character", default = "", help = "Process only the specified reviewer: cassia/in_house/enrichment (subset run)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Project root (all relative paths are resolved from this directory) ----
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0L) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
triageHome <- normalizePath(triageHome, winslash = "/", mustWork = TRUE)

# ---- Convert out_dir to an absolute path ----
is_absolute_path <- function(x) grepl("^/", x) || grepl("^[A-Za-z]:[/\\\\]", x)
out_dir <- if (is_absolute_path(opt$out_dir)) opt$out_dir else file.path(triageHome, opt$out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Version (update when prompts, decision rules or verifier schema change) ----
PROMPT_VERSION <- "v1.0.0"
PIPELINE_VERSION <- "v1.0.0"
TOP_K_CANDIDATES <- 5L
MARKER_LIMIT <- 20L

# ---- Use a fixed run manifest (do not select the latest run) ----
manifest_path <- if (nzchar(opt$manifest)) opt$manifest else file.path(triageHome, "reproducibility", "config", "evidence_mapping_runs.tsv")
if (!file.exists(manifest_path)) stop("Fixed run manifest is missing: ", manifest_path)
manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
get_run <- function(ds) {
  m <- manifest %>% filter(dataset == ds)
  if (nrow(m) == 0L) stop("Manifest is missing dataset: ", ds)
  if (nrow(m) > 1L) stop("Dataset appears more than once in manifest: ", ds)
  run_dir <- Triage:::resolve_manifest_run_dir(triageHome, ds, m$run_dir[1])
  if (!dir.exists(run_dir)) stop("Run directory referenced by the manifest does not exist: ", run_dir)
  run_dir
}

# ---- Source shared modules ----
suppressMessages(library(Triage))
suppressMessages(library(Triage))
suppressMessages(library(Triage))
# Package-internal null coalescing is used throughout this script.
`%||%` <- Triage:::`%||%`

# ---- Load packaged mapper/verifier prompts (single source of truth) ----
prompts_path <- system.file("prompts", "cl_mapper_prompts.R", package = "Triage")
if (!nzchar(prompts_path) || !file.exists(prompts_path)) {
  stop("Installed Triage package is missing inst/prompts/cl_mapper_prompts.R; reinstall the package.")
}
source(prompts_path, local = TRUE)
idx <- Triage:::build_cl_index(opt$cl_json)
ontology_checksum <- unname(tools::md5sum(opt$cl_json))

# ---- Deterministic methods ( CL-Linker Return) ----
DETERMINISTIC_METHODS <- c(
  "provided_valid_clid", "canonical_exact", "synonym_exact",
  "normalized", "alias"
)

# ============ LLM calls ============
call_llm_json <- function(prompt, system_prompt, api_key, model = opt$model,
                          temperature = opt$temperature, max_tokens = opt$max_tokens,
                          timeout_sec = 180L, retries = 3L) {
  if (!nzchar(api_key)) return(list(ok = FALSE, data = NULL, error_type = "missing_api_key"))
  body <- list(model = model,
               messages = list(list(role = "system", content = system_prompt),
                               list(role = "user", content = prompt)),
               thinking = list(type = "disabled"), temperature = temperature,
               top_p = 1, max_tokens = max_tokens,
               response_format = list(type = "json_object"))
  for (attempt in seq_len(retries)) {
    response <- tryCatch(
      httr::POST(Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = "XXXXX")),
                 httr::add_headers(Authorization = paste("Bearer", api_key)),
                 body = body, encode = "json",
                 httr::config(connecttimeout = 15), httr::timeout(timeout_sec)),
      error = function(e) e)
    if (inherits(response, "error")) {
      if (attempt < retries) { Sys.sleep(2^(attempt-1L)); next }
      return(list(ok = FALSE, data = NULL, error_type = "request_error", error_message = conditionMessage(response)))
    }
    status <- httr::status_code(response)
    raw_text <- httr::content(response, as = "text", encoding = "UTF-8")
    if (status == 429L || status >= 500L) { if (attempt < retries) { Sys.sleep(2^(attempt-1L)); next } }
    if (status != 200L) return(list(ok = FALSE, data = NULL, error_type = "http_error", http_status = status))
    payload <- tryCatch(jsonlite::fromJSON(raw_text, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(payload) || is.null(payload$choices)) {
      if (attempt < retries) { Sys.sleep(2^(attempt-1L)); next }
      return(list(ok = FALSE, data = NULL, error_type = "invalid_payload"))
    }
    finish <- tryCatch(payload$choices[[1]]$finish_reason, error = function(e) NULL)
    if (identical(finish, "length")) { if (attempt < retries) { Sys.sleep(2^(attempt-1L)); next }
      return(list(ok = FALSE, data = NULL, error_type = "truncated_length")) }
    txt <- tryCatch(payload$choices[[1]]$message$content, error = function(e) NULL)
    if (is.null(txt) || !nzchar(trimws(txt))) { if (attempt < retries) { Sys.sleep(2^(attempt-1L)); next }
      return(list(ok = FALSE, data = NULL, error_type = "empty_content")) }
    parsed <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
    if (is.null(parsed)) {
      if (attempt < retries) {
        Sys.sleep(2^(attempt - 1L))
        next
      }
      return(list(
        ok = FALSE,
        data = NULL,
        error_type = "invalid_json",
        raw_content = substr(txt, 1L, 1000L)
      ))
    }
    return(list(ok = TRUE, data = parsed, finish_reason = finish))
  }
}
api_key <- Sys.getenv(opt$api_key_env, unset = "")

# ============ Normalize LLM scalar outputs ============
normalise_llm_scalar <- function(x, default) {
  value <- as.character(x %||% default)
  if (length(value) == 0L || is.na(value[1])) value <- default else value <- value[1]
  tolower(trimws(value))
}
ALLOWED_MAPPER_STATUS <- c("mapped", "parent_only", "mixed_identity", "mixed_signal",
                           "ambiguous", "malformed", "unmapped", "not_a_cell_type")
STABLE_NON_SINGLE_STATUS <- c("mixed_identity", "mixed_signal", "ambiguous",
                              "malformed", "unmapped", "not_a_cell_type")
ALLOWED_VERIFIER_VERDICT <- c(
  "accept", "accept_parent", "ambiguous", "composite_multi_identity", "review"
)

normalise_verifier_verdict <- function(x) {
  value <- normalise_llm_scalar(x, "review")
  value <- gsub("[ -]+", "_", value)
  aliases <- c(
    "parent_only" = "accept_parent",
    "composite" = "composite_multi_identity",
    "mark_composite" = "composite_multi_identity"
  )
  if (value %in% names(aliases)) value <- unname(aliases[[value]])
  value
}

# Append the strict registry schema here so
# this script is self-contained and the cache hash changes with the schema.
MAPPER_SCHEMA_ADDENDUM <- paste(
  "REGISTRY MAPPING OUTPUT CONTRACT:",
  "Return one JSON object only.",
  "Allowed status values: mapped, parent_only, mixed_identity, mixed_signal, ambiguous, malformed, unmapped, not_a_cell_type.",
  "For mapped or parent_only: selected_candidate_rank must be one integer referring to the displayed candidate list.",
  "For mixed_identity: component_candidate_ranks must contain at least two distinct candidate ranks.",
  "mixed_identity is reserved for multiple core cell identities (e.g. platelet/megakaryocyte, macrophage/neutrophil).",
  "Do not call one core identity with state/modifier terms mixed_identity (e.g. TREM2+ lipid-associated macrophage, activated B cell).",
  "For mixed_signal: use when one primary identity is accompanied by contamination/doublet/ambient signal.",
  "For malformed: use when the label is truncated or syntactically incomplete.",
  "Do not generate a CL ID outside the candidate list.",
  "Required keys: status, selected_candidate_rank, component_candidate_ranks, specificity_relation, confidence_category, reason.",
  sep = "\n"
)

VERIFIER_SCHEMA_ADDENDUM <- paste(
  "REGISTRY VERIFIER OUTPUT CONTRACT:",
  "Allowed verdict values: accept, accept_parent, ambiguous, composite_multi_identity, review.",
  "If verdict=accept_parent, selected_parent_rank is required and must refer to the displayed forward candidate list.",
  "If verdict=composite_multi_identity, component_candidate_ranks must contain at least two distinct ranks.",
  "Do not treat multiple modifiers of one core identity as multiple identities.",
  "Return one JSON object only.",
  sep = "\n"
)

# ============ candidates ============
safe_text <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(NA_character_)
  value <- trimws(as.character(x[1]))
  if (!nzchar(value) || toupper(value) == "NA" || tolower(value) == "null") return(NA_character_)
  value
}

split_pipe <- function(x) {
  value <- safe_text(x)
  if (is.na(value)) return(character(0))
  out <- trimws(unlist(strsplit(value, "\\|", perl = TRUE), use.names = FALSE))
  unique(out[!is.na(out) & nzchar(out)])
}

extract_rank <- function(d, cands_df) {
  rk <- suppressWarnings(as.integer(d$selected_candidate_rank %||% NA))
  if (length(rk) == 1L && !is.na(rk) && rk >= 1L && rk <= nrow(cands_df)) return(rk)
  NA_integer_
}

extract_component_set <- function(d, cands_df) {
  ranks <- suppressWarnings(as.integer(unlist(d$component_candidate_ranks %||% integer(0), use.names = FALSE)))
  ranks <- sort(unique(ranks[!is.na(ranks) & ranks >= 1L & ranks <= nrow(cands_df)]))
  list(
    ranks = ranks,
    ids = if (length(ranks) > 0L) as.character(cands_df$cl_id[ranks]) else character(0),
    labels = if (length(ranks) > 0L) as.character(cands_df$canonical_label[ranks]) else character(0)
  )
}

# ============ Cache ============
cache_get <- function(hash) {
  fp <- file.path(cache_dir, paste0(hash, ".rds"))
  if (file.exists(fp)) readRDS(fp) else NULL
}
cache_set <- function(hash, obj) saveRDS(obj, file.path(cache_dir, paste0(hash, ".rds")))
prompt_hash <- function(system_prompt, user_prompt, model, temperature) {
  digest::digest(list(system_prompt = system_prompt, user_prompt = user_prompt,
                      model = model, temperature = temperature, top_p = 1,
                      max_tokens = opt$max_tokens, ontology_checksum = ontology_checksum,
                      prompt_version = PROMPT_VERSION),
                 algo = "sha256")
}

# ============ Deterministic ontology checks ============
is_ancestor_cl <- function(a, b) {
  if (is.na(a) || is.na(b)) return(FALSE)
  node <- idx$cl[[b]]
  if (is.null(node)) return(FALSE)
  anc <- node$ancestors %||% NULL
  if (is.null(anc)) return(FALSE)
  ancestor_ids <- unique(c(names(anc), if (is.atomic(anc)) as.character(unname(anc)) else character(0)))
  a %in% ancestor_ids
}
resolve_mismatch <- function(clid_a, clid_b) {
  if (is_ancestor_cl(clid_a, clid_b)) return(list(clid = clid_a, status = "parent_only", note = "A_ancestor_of_B"))
  if (is_ancestor_cl(clid_b, clid_a)) return(list(clid = clid_b, status = "parent_only", note = "B_ancestor_of_A"))
  list(clid = NA_character_, status = "review", note = "sibling_or_cross_lineage")
}

# ============ Extract marker context (top-ranked genes) ============
get_markers <- function(ds, cid) {
  run <- get_run(ds)
  fp <- file.path(run, "05c_llm_queries", "step1_report_queries", paste0(cid, "_step1_report_query.json"))
  if (!file.exists(fp)) return(character(0))
  d <- tryCatch(jsonlite::fromJSON(fp, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(d)) return(character(0))
  genes <- d$input_data$cluster_dossier$degs$in_scope_top_genes
  if (is.null(genes) || length(genes) == 0L) return(character(0))
  if (is.data.frame(genes)) {
    gene_col <- intersect(c("gene", "gene_symbol", "symbol"), names(genes))[1]
    if (is.na(gene_col)) stop("Cannot identify gene column in marker table: ", fp)
    genes <- genes[[gene_col]]
  } else {
    genes <- unlist(genes, use.names = FALSE)
  }
  genes <- trimws(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes) & toupper(genes) != "NA"]
  head(unique(genes), MARKER_LIMIT)
}

# ============ context automatic (from dataset_config) ============
safe_meta_value <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return("not provided")
  value <- trimws(as.character(x[1]))
  if (!nzchar(value) || toupper(value) == "NA") return("not provided")
  value
}

get_ds_meta <- function(ds) {
  cfg <- tryCatch(Triage:::get_dataset_config(ds, triageHome), error = function(e) NULL)
  if (is.null(cfg)) {
    return(list(
      tissue = "not provided",
      species = "not provided",
      condition = "not provided"
    ))
  }
  list(
    tissue = safe_meta_value(cfg$tissue),
    species = safe_meta_value(cfg$species),
    condition = safe_meta_value(cfg$study_context)
  )
}

# ============ Collect reviewer instances (deduplicate) ============
REVIEWER_MAP <- c(cassia = "cassia", in_house = "in_house", enrichment = "enrichment")
REVIEWER_SOURCE <- list(
  cassia = list(label_col = "Predicted Main Cell Type", cluster_col = "True Cell Type",
                clid_col = NA_character_, conf_col = "Score"),
  in_house = list(label_col = "main_pred", cluster_col = "cluster_id",
                  clid_col = "main_cell_ontology_id", conf_col = "main_confidence"),
  enrichment = list(label_col = "predicted_label", cluster_col = "cluster_id",
                    clid_col = "cl_id", conf_col = "confidence")
)

datasets <- if (nzchar(opt$dataset)) opt$dataset else c("Census_immune", "Sikkema_lung", "TS_kidney", "TS_pancreas", "Zheng_blood")

get_step1_cluster_ids <- function(run_dir) {
  qdir <- file.path(run_dir, "05c_llm_queries", "step1_report_queries")
  if (!dir.exists(qdir)) stop("Missing step1 query directory: ", qdir)
  fs <- list.files(qdir, pattern = "_step1_report_query\\.json$", full.names = FALSE)
  ids <- sub("_step1_report_query\\.json$", "", fs)
  ids <- safe_text(ids)
  # safe_text is scalar; process vector explicitly
  ids <- trimws(as.character(sub("_step1_report_query\\.json$", "", fs)))
  unique(ids[!is.na(ids) & nzchar(ids)])
}

pick_cassia_source <- function(run_dir, expected_ids, dataset_name) {
  candidates <- Sys.glob(file.path(
    run_dir, "03b_cassia", "CASSIA_*", "01_annotation_results",
    "annotation_cassia_FINAL_RESULTS.csv"
  ))
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0L) stop("No CASSIA result found for dataset: ", dataset_name)

  score_one <- function(fp) {
    d <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(d) || !"True Cell Type" %in% names(d)) return(-Inf)
    ids <- trimws(as.character(d[["True Cell Type"]]))
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])
    length(intersect(ids, expected_ids))
  }
  scores <- vapply(candidates, score_one, numeric(1))
  if (!any(is.finite(scores))) stop("No readable CASSIA result with 'True Cell Type' for: ", dataset_name)
  best <- which(scores == max(scores, na.rm = TRUE))
  if (length(best) > 1L) {
    mt <- suppressWarnings(as.numeric(file.info(candidates[best])$mtime))
    best <- best[which.max(mt)]
  }
  fp <- candidates[best[1]]
  d <- readr::read_csv(fp, show_col_types = FALSE)
  ids <- trimws(as.character(d[["True Cell Type"]]))
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!setequal(ids, expected_ids)) {
    stop(
      "CASSIA cluster set does not match step1 clusters for ", dataset_name,
      "\nMissing in CASSIA: ", paste(setdiff(expected_ids, ids), collapse = ", "),
      "\nExtra in CASSIA: ", paste(setdiff(ids, expected_ids), collapse = ", ")
    )
  }
  fp
}

source_path_for <- function(run_dir, reviewer, expected_ids, dataset_name) {
  if (reviewer == "cassia") return(pick_cassia_source(run_dir, expected_ids, dataset_name))
  if (reviewer == "in_house") return(file.path(run_dir, "07_our_summary", "summary.csv"))
  if (reviewer == "enrichment") return(file.path(run_dir, "07b_inter_summary", "summary.csv"))
  stop("Unknown reviewer: ", reviewer)
}

instances <- list()
for (ds in datasets) {
  run <- get_run(ds)
  expected_ids <- get_step1_cluster_ids(run)
  if (length(expected_ids) == 0L) stop("No step1 clusters found for dataset: ", ds)

  reviewer_cluster_sets <- list()
  for (rv in names(REVIEWER_MAP)) {
    src <- REVIEWER_SOURCE[[rv]]
    fp <- source_path_for(run, rv, expected_ids, ds)
    if (length(fp) != 1L || is.na(fp) || !file.exists(fp)) {
      stop("Expected one source file for ", ds, " / ", rv, ": ", paste(fp, collapse = ", "))
    }
    d <- readr::read_csv(fp, show_col_types = FALSE)
    required_cols <- c(src$label_col, src$cluster_col)
    missing_cols <- setdiff(required_cols, names(d))
    if (length(missing_cols) > 0L) {
      stop("Missing columns in ", fp, ": ", paste(missing_cols, collapse = ", "))
    }

    ids_seen <- character(0)
    for (j in seq_len(nrow(d))) {
      lbl <- safe_text(d[[src$label_col]][j])
      cid <- safe_text(d[[src$cluster_col]][j])
      if (is.na(cid)) stop("Missing cluster_id in ", fp, " row ", j)
      ids_seen <- c(ids_seen, cid)
      if (is.na(lbl)) stop("Missing reviewer label in ", fp, " row ", j, " (cluster ", cid, ")")

      provided_clid <- if (!is.na(src$clid_col) && src$clid_col %in% names(d)) {
        safe_text(d[[src$clid_col]][j])
      } else NA_character_
      conf <- if (!is.na(src$conf_col) && src$conf_col %in% names(d)) {
        suppressWarnings(as.numeric(d[[src$conf_col]][j]))
      } else NA_real_

      instances[[length(instances) + 1L]] <- list(
        dataset = ds,
        cluster_id = cid,
        reviewer = unname(REVIEWER_MAP[[rv]]),
        raw_label = lbl,
        provided_clid = provided_clid,
        conf = conf,
        source_file = normalizePath(fp, winslash = "/", mustWork = TRUE),
        source_run_id = basename(run)
      )
    }
    reviewer_cluster_sets[[rv]] <- unique(ids_seen)
    if (!setequal(reviewer_cluster_sets[[rv]], expected_ids)) {
      stop(
        "Reviewer cluster set mismatch for ", ds, " / ", rv,
        "\nMissing: ", paste(setdiff(expected_ids, reviewer_cluster_sets[[rv]]), collapse = ", "),
        "\nExtra: ", paste(setdiff(reviewer_cluster_sets[[rv]], expected_ids), collapse = ", ")
      )
    }
  }
}

instance_df_full <- dplyr::bind_rows(lapply(instances, as.data.frame, stringsAsFactors = FALSE))
dup <- instance_df_full %>% count(dataset, cluster_id, reviewer) %>% filter(n != 1L)
if (nrow(dup) > 0L) stop("Duplicate/missing reviewer records:\n", paste(capture.output(print(dup)), collapse = "\n"))
reviewer_counts <- instance_df_full %>% count(reviewer)
if (!setequal(reviewer_counts$reviewer, names(REVIEWER_MAP))) stop("Reviewer set is incomplete")
if (nzchar(opt$dataset)) {
  expected_per_reviewer <- length(unique(instance_df_full$cluster_id))
  if (nrow(instance_df_full) != 3L * expected_per_reviewer ||
      any(reviewer_counts$n != expected_per_reviewer)) {
    stop("Expected one record per cluster for each of three reviewers; observed ",
         nrow(instance_df_full), "\n",
         paste(capture.output(print(reviewer_counts)), collapse = "\n"))
  }
} else if (nrow(instance_df_full) != 150L || any(reviewer_counts$n != 50L)) {
  stop("Expected 150 instances and 50 per reviewer; observed ", nrow(instance_df_full),
       "\n", paste(capture.output(print(reviewer_counts)), collapse = "\n"))
}

is_filtered_run <- nzchar(opt$dataset) || nzchar(opt$cluster_id) || nzchar(opt$reviewer)
scope_parts <- c(if (nzchar(opt$dataset)) opt$dataset else NULL,
                 if (nzchar(opt$cluster_id)) opt$cluster_id else NULL,
                 if (nzchar(opt$reviewer)) opt$reviewer else NULL)
scope_tag <- if (is_filtered_run) {
  paste0("_subset_", gsub("[^A-Za-z0-9_]+", "_", paste(scope_parts, collapse = "_")))
} else ""

if (nzchar(opt$dataset)) instances <- Filter(function(x) x$dataset == opt$dataset, instances)
if (nzchar(opt$cluster_id)) instances <- Filter(function(x) x$cluster_id == opt$cluster_id, instances)
if (nzchar(opt$reviewer)) instances <- Filter(function(x) x$reviewer == opt$reviewer, instances)
if (length(instances) == 0L) stop("No instances remain after subset filters")
cat("Total instances:", length(instances), if (is_filtered_run) "(subset run)" else "(full registry)", "\n")

registry_filename <- paste0("reviewer_mapping_registry", scope_tag, ".tsv")
# ============ Resume from partial output (Version) ============
partial_file <- file.path(out_dir, paste0("cl_link_", PIPELINE_VERSION, scope_tag, "_", substr(ontology_checksum, 1, 8), ".partial.tsv"))
if (file.exists(partial_file)) {
  existing <- readr::read_tsv(partial_file, show_col_types = FALSE)
  if (nrow(existing) > 0) {
    req_fields <- c("prompt_version", "pipeline_version", "ontology_checksum", "model",
                    "temperature", "max_tokens", "top_k", "marker_limit")
    if (!all(req_fields %in% names(existing))) stop("Partial output is missing version fields: ", partial_file)
    if (!all(existing$prompt_version == PROMPT_VERSION) ||
        !all(existing$pipeline_version == PIPELINE_VERSION) ||
        !all(existing$ontology_checksum == ontology_checksum) ||
        !all(existing$model == opt$model) ||
        !all(as.numeric(existing$temperature) == opt$temperature) ||
        !all(as.integer(existing$max_tokens) == opt$max_tokens) ||
        !all(as.integer(existing$top_k) == TOP_K_CANDIDATES) ||
        !all(as.integer(existing$marker_limit) == MARKER_LIMIT)) {
      stop("Partial output version mismatch; refusing to mix versions: ", partial_file)
    }
  }
  done_keys <- paste(existing$dataset, existing$cluster_id, existing$reviewer, existing$raw_label, sep = "||")
  rows <- split(existing, seq_len(nrow(existing)))
} else {
  done_keys <- character(0); rows <- list()
}

# ============ Main loop ============
i <- 0
for (inst in instances) {
  key <- paste(inst$dataset, inst$cluster_id, inst$reviewer, inst$raw_label, sep = "||")
  if (key %in% done_keys) next
  i <- i + 1
  cat(sprintf("[%d/%d] %s %s %s: %s\n", i, length(instances), inst$dataset, inst$cluster_id,
              inst$reviewer, substr(inst$raw_label, 1, 40))); flush.console()

  # ---- Step 1: ----
  m <- Triage:::cl_link(inst$raw_label, inst$provided_clid %||% NA_character_, idx)
  det_clid <- as.character(m$cl_id %||% NA_character_)[1]
  det_method <- as.character(m$mapping_method %||% NA_character_)[1]
  det_label <- as.character(m$canonical_label %||% NA_character_)[1]
  if (!is.na(det_clid) && nzchar(det_clid) &&
      !is.na(det_method) && det_method %in% DETERMINISTIC_METHODS) {
    meta_det <- get_ds_meta(inst$dataset)
    single_ready <- TRUE
    row <- data.frame(dataset = inst$dataset, cluster_id = inst$cluster_id, reviewer = inst$reviewer,
                      raw_label = inst$raw_label, normalised_label = idx$norm_form(inst$raw_label),
                      provided_cl_id = inst$provided_clid %||% NA_character_,
                      reviewer_confidence = inst$conf,
                      markers = NA_character_, tissue = meta_det$tissue, species = meta_det$species,
                      condition = meta_det$condition,
                      candidate_ids = det_clid, candidate_labels = det_label,
                      candidate_scores = NA_character_,
                      component_cl_ids = NA_character_, component_labels = NA_character_, component_agreement = NA,
                      final_status = "mapped", decision = "deterministic",
                      cl_id = det_clid, canonical_label = det_label,
                      mapping_method = det_method, specificity_relation = "exact",
                      conf_cat = NA, reason = "deterministic",
                      proposed_cl_id = det_clid, proposed_canonical_label = det_label,
                      single_cl_ready = single_ready, adjudication_ready = TRUE,
                      evaluation_ready = single_ready,
                      editor_input_type = "single_cl_id",
                      mapper_a_status = NA, mapper_a_rank = NA_integer_, mapper_a_cl_id = NA,
                      mapper_a_specificity = NA, mapper_a_confidence = NA,
                      mapper_b_status = NA, mapper_b_rank = NA_integer_, mapper_b_cl_id = NA,
                      mapper_b_specificity = NA, mapper_b_confidence = NA,
                      ab_clid_agreement = NA, ab_specificity_agreement = NA,
                      verifier_called = FALSE, verifier_verdict = NA, verifier_reason = NA, verifier_action = NA,
                      hash_a = NA, hash_b = NA, hash_v = NA,
                      source_file = inst$source_file, source_run_id = inst$source_run_id,
                      model = opt$model, prompt_version = PROMPT_VERSION, pipeline_version = PIPELINE_VERSION,
                      ontology_version = "2025-07-30", ontology_checksum = ontology_checksum,
                      temperature = opt$temperature, max_tokens = opt$max_tokens,
                      top_k = TOP_K_CANDIDATES, marker_limit = MARKER_LIMIT,
                      stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- row
    write_tsv(bind_rows(rows), partial_file)
    next
  }

  # ---- Step 2: lexical top-k ----
  cands <- Triage:::lexical_candidates(inst$raw_label, idx, top_k = TOP_K_CANDIDATES)
  markers <- get_markers(inst$dataset, inst$cluster_id)
  meta <- get_ds_meta(inst$dataset)
  context_lines <- paste(sprintf("Species: %s\nTissue: %s\nDisease/condition: %s\nDataset: %s\nCluster: %s",
                                 meta$species, meta$tissue, meta$condition, inst$dataset, inst$cluster_id))
  marker_lines <- if (length(markers) > 0) paste("Top-ranked marker genes (fixed order):\n",
                                                 paste(sprintf("  %s", markers), collapse = "\n")) else "Top-ranked marker genes: not provided"

  # Initialize (reset each iteration to prevent state carry-over)
  final_status <- "review"; decision <- ""; final_clid <- NA; final_label <- NA; conf_cat <- NA; reason <- ""
  ab_clid_agree <- NA; ab_spec_agree <- NA
  proposed_clid <- NA_character_; proposed_label <- NA_character_
  ma_status <- NA; ma_rank <- NA_integer_; ma_clid <- NA; ma_spec <- NA; ma_conf <- NA
  mb_status <- NA; mb_rank <- NA_integer_; mb_clid <- NA; mb_spec <- NA; mb_conf <- NA
  verifier_called <- FALSE; verifier_verdict <- NA; verifier_reason <- NA; verifier_action <- NA
  component_clids <- character(0); component_labels <- character(0); component_agree <- NA
  hash_a <- NA; hash_b <- NA; hash_v <- NA

  if (nrow(cands) == 0) {
    final_status <- "unmapped"; decision <- "no_candidates"
    reason <- "no lexical candidates"
  } else {
    # ---- Step 3: Mapper A (forward order) ----
    prompt_a <- paste(
      MAPPER_USER_TEMPLATE_EVIDENCE(inst$raw_label, idx$norm_form(inst$raw_label), cands, marker_lines, context_lines),
      "\n\n", MAPPER_SCHEMA_ADDENDUM
    )
    hash_a <- prompt_hash(MAPPER_SYSTEM_PROMPT, prompt_a, opt$model, opt$temperature)
    cached_a <- if (isTRUE(opt$no_cache)) NULL else cache_get(hash_a)
    ma <- if (!is.null(cached_a)) cached_a else {
      r <- call_llm_json(prompt_a, MAPPER_SYSTEM_PROMPT, api_key, model = opt$model)
      if (r$ok) cache_set(hash_a, r)
      r
    }
    # Mapper B (candidates)
    cands_rev <- cands[rev(seq_len(nrow(cands))), ]
    prompt_b <- paste(
      MAPPER_USER_TEMPLATE_EVIDENCE(inst$raw_label, idx$norm_form(inst$raw_label), cands_rev, marker_lines, context_lines),
      "\n\n", MAPPER_SCHEMA_ADDENDUM
    )
    hash_b <- prompt_hash(MAPPER_SYSTEM_PROMPT, prompt_b, opt$model, opt$temperature)
    cached_b <- if (isTRUE(opt$no_cache)) NULL else cache_get(hash_b)
    mb <- if (!is.null(cached_b)) cached_b else {
      r <- call_llm_json(prompt_b, MAPPER_SYSTEM_PROMPT, api_key, model = opt$model)
      if (r$ok) cache_set(hash_b, r)
      r
    }

    if (!ma$ok || !mb$ok) {
      decision <- paste0("llm_error_", ma$error_type %||% "", "_", mb$error_type %||% "")
      reason <- decision
    } else {
      da <- ma$data; db <- mb$data
      ma_status <- normalise_llm_scalar(da$status, "unmapped")
      mb_status <- normalise_llm_scalar(db$status, "unmapped")
      ma_spec <- normalise_llm_scalar(da$specificity_relation, "exact")
      mb_spec <- normalise_llm_scalar(db$specificity_relation, "exact")
      ma_conf <- as.character(da$confidence_category %||% NA)
      mb_conf <- as.character(db$confidence_category %||% NA)
      if (!ma_status %in% ALLOWED_MAPPER_STATUS) ma_status <- "invalid_output"
      if (!mb_status %in% ALLOWED_MAPPER_STATUS) mb_status <- "invalid_output"

      ma_rank <- extract_rank(da, cands)
      mb_rank <- extract_rank(db, cands_rev)
      ma_clid <- if (!is.na(ma_rank)) cands$cl_id[ma_rank] else NA_character_
      mb_clid <- if (!is.na(mb_rank)) cands_rev$cl_id[mb_rank] else NA_character_
      ab_clid_agree <- identical(ma_clid, mb_clid) && !is.na(ma_clid)
      ab_spec_agree <- identical(tolower(ma_spec), tolower(mb_spec))
      proposed_clid <- if (!is.na(ma_clid)) ma_clid else if (!is.na(mb_clid)) mb_clid else NA_character_
      proposed_label <- if (!is.na(proposed_clid)) cands$canonical_label[which(cands$cl_id == proposed_clid)][1] else NA_character_

      # ---- Step 4/5: decision logic ----
      if (identical(ma_status, mb_status) && ma_status == "mixed_identity") {
        comp_a <- extract_component_set(da, cands)
        comp_b <- extract_component_set(db, cands_rev)
        component_agree <- setequal(comp_a$ids, comp_b$ids)
        if (isTRUE(component_agree) && length(comp_a$ids) >= 2L) {
          final_status <- "mixed_identity"
          decision <- "ab_consistent_mixed_identity"
          component_clids <- sort(unique(comp_a$ids))
          component_labels <- vapply(component_clids, function(id) {
            hit <- match(id, cands$cl_id)
            if (!is.na(hit)) as.character(cands$canonical_label[hit]) else id
          }, character(1))
          reason <- "Mapper A/B agreed on the mixed-identity component set"
        } else {
          final_status <- "review"
          decision <- "mixed_component_mismatch"
          reason <- paste0("A components=", paste(comp_a$ids, collapse = "|"),
                           "; B components=", paste(comp_b$ids, collapse = "|"))
        }
      } else if (identical(ma_status, mb_status) && ma_status %in% STABLE_NON_SINGLE_STATUS) {
        final_status <- ma_status
        decision <- paste0("ab_consistent_", ma_status)
        reason <- paste0("Mapper A/B consistently returned ", ma_status)
      } else if (identical(ma_status, mb_status) && ma_status == "parent_only") {
        if (ab_clid_agree && ab_spec_agree) {
          final_status <- "parent_only"
          decision <- "ab_consistent_parent_only"
          proposed_clid <- ma_clid
          proposed_label <- cands$canonical_label[match(ma_clid, cands$cl_id)]
          reason <- "Mapper A/B agreed on the same broader CL term"
        } else {
          final_status <- "review"
          decision <- "parent_only_mismatch"
          reason <- paste0("A=", ma_clid, "; B=", mb_clid,
                           "; specA=", ma_spec, "; specB=", mb_spec)
        }
      } else if (!identical(ma_status, mb_status)) {
        final_status <- "review"
        decision <- paste0("status_mismatch_", ma_status, "_", mb_status)
        reason <- decision
      } else if (ma_status == "mapped" && !is.na(ma_clid) && !is.na(mb_clid)) {
        if (ab_clid_agree) {
          if (!ab_spec_agree) {
            final_status <- "review"
            decision <- "specificity_mismatch"
            reason <- paste0("A spec=", ma_spec, "; B spec=", mb_spec)
          } else {
            margin <- if (nrow(cands) > 1L) cands$lexical_score[1] - cands$lexical_score[2] else NA_real_
            rk_a <- ma_rank
            high_risk <- grepl("\\b(or|and|mixed|doublet|like)\\b|/", tolower(inst$raw_label)) ||
              nchar(inst$raw_label) <= 6L || is.na(margin) || margin < 0.10 ||
              (!is.na(rk_a) && rk_a != 1L) || identical(ma_spec, "potentially_over_specific")

            if (high_risk) {
              verifier_called <- TRUE
              other_idx <- setdiff(seq_len(nrow(cands)), rk_a)
              runner_up_idx <- if (length(other_idx) > 0L) other_idx[which.max(cands$lexical_score[other_idx])] else NA_integer_
              alt_label <- if (!is.na(runner_up_idx)) {
                sprintf("%s (%s)", cands$canonical_label[runner_up_idx], cands$cl_id[runner_up_idx])
              } else "none"
              sel_label <- cands$canonical_label[match(ma_clid, cands$cl_id)]
              v_prompt <- paste(
                VERIFIER_USER_TEMPLATE(
                  inst$raw_label,
                  sprintf("%s (%s)", sel_label, ma_clid),
                  alt_label,
                  paste(sprintf("%d. %s (%s)", seq_len(nrow(cands)), cands$canonical_label, cands$cl_id), collapse = "\n")
                ),
                "\n\nCLUSTER CONTEXT:\n", context_lines,
                "\n\nMARKER EVIDENCE:\n", marker_lines,
                "\n\nMAPPER A JSON:\n", jsonlite::toJSON(da, auto_unbox = TRUE),
                "\n\nMAPPER B JSON:\n", jsonlite::toJSON(db, auto_unbox = TRUE),
                "\n\n", VERIFIER_SCHEMA_ADDENDUM
              )
              hash_v <- prompt_hash(VERIFIER_SYSTEM_PROMPT, v_prompt, opt$model, opt$temperature)
              cached_v <- if (isTRUE(opt$no_cache)) NULL else cache_get(hash_v)
              vv <- if (!is.null(cached_v)) cached_v else {
                r <- call_llm_json(v_prompt, VERIFIER_SYSTEM_PROMPT, api_key, model = opt$model)
                if (r$ok) cache_set(hash_v, r)
                r
              }

              if (vv$ok) {
                verifier_verdict <- normalise_verifier_verdict(vv$data$verdict)
                verifier_reason <- safe_text(vv$data$warning %||% vv$data$reason %||% NA_character_)
                verifier_action <- safe_text(vv$data$suggested_action %||% NA_character_)
                if (!verifier_verdict %in% ALLOWED_VERIFIER_VERDICT) {
                  verifier_reason <- paste0("Invalid verifier verdict: ", verifier_verdict,
                                            if (!is.na(verifier_reason)) paste0(" | ", verifier_reason) else "")
                  verifier_verdict <- "invalid_output"
                  final_status <- "review"
                  decision <- "verifier_invalid_output"
                } else if (verifier_verdict == "accept") {
                  final_status <- "mapped_evidence"
                  decision <- "accept_after_verifier"
                } else if (verifier_verdict == "accept_parent") {
                  parent_rank <- suppressWarnings(as.integer(vv$data$selected_parent_rank %||% NA))
                  valid_parent_rank <- length(parent_rank) == 1L && !is.na(parent_rank) &&
                    parent_rank >= 1L && parent_rank <= nrow(cands)
                  if (valid_parent_rank) {
                    parent_candidate <- cands$cl_id[parent_rank]
                    if (identical(parent_candidate, proposed_clid) || is_ancestor_cl(parent_candidate, proposed_clid)) {
                      final_status <- "parent_only"
                      decision <- "verifier_accept_parent"
                      proposed_clid <- parent_candidate
                      proposed_label <- idx$clid_info[[parent_candidate]]$canonical_label %||% parent_candidate
                    } else {
                      final_status <- "review"
                      decision <- "verifier_invalid_parent"
                    }
                  } else {
                    final_status <- "review"
                    decision <- "verifier_accept_parent_missing_rank"
                  }
                } else if (verifier_verdict == "ambiguous") {
                  final_status <- "ambiguous"
                  decision <- "verifier_ambiguous"
                } else if (verifier_verdict == "composite_multi_identity") {
                  comp_v <- extract_component_set(vv$data, cands)
                  if (length(comp_v$ids) >= 2L) {
                    final_status <- "mixed_identity"
                    decision <- "verifier_mixed_identity"
                    component_clids <- sort(unique(comp_v$ids))
                    component_labels <- vapply(component_clids, function(id) {
                      hit <- match(id, cands$cl_id)
                      if (!is.na(hit)) as.character(cands$canonical_label[hit]) else id
                    }, character(1))
                  } else {
                    final_status <- "review"
                    decision <- "verifier_mixed_missing_components"
                  }
                } else {
                  final_status <- "review"
                  decision <- "verifier_review"
                }
              } else {
                final_status <- "review"
                decision <- paste0("verifier_error_", vv$error_type)
              }
            } else {
              final_status <- "mapped_evidence"
              decision <- "accept_low_risk"
            }
            conf_cat <- ma_conf
            reason <- if (isTRUE(verifier_called) && !is.na(verifier_reason) && nzchar(verifier_reason)) {
              paste0("A/B consistent; spec=", ma_spec, " | verifier: ", verifier_reason)
            } else {
              paste0("A/B consistent; spec=", ma_spec)
            }
          }
        } else {
          mm <- resolve_mismatch(ma_clid, mb_clid)
          if (mm$status == "parent_only") {
            final_status <- "parent_only"
            decision <- "ab_parent_child"
            proposed_clid <- mm$clid
            proposed_label <- idx$clid_info[[mm$clid]]$canonical_label %||% mm$clid
            reason <- mm$note
          } else {
            final_status <- "review"
            decision <- "ab_mismatch"
            reason <- paste0("A=", ma_clid, " B=", mb_clid, " | ", mm$note)
          }
        }
      } else {
        final_status <- "review"
        decision <- "mapped_but_incomplete"
        reason <- paste0("statusA=", ma_status, " statusB=", mb_status)
      }
    }
  }

  # ---- Step 6: write cl_id only for formally accepted mappings ----
  accepted_statuses <- c("mapped_evidence", "parent_only")
  if (final_status %in% accepted_statuses && !is.na(proposed_clid)) {
    final_clid <- proposed_clid
    final_label <- proposed_label
  } else if (!(final_status %in% accepted_statuses)) {
    final_clid <- NA_character_
    final_label <- NA_character_
  }

  # registry status flags
  single_ready <- final_status %in% accepted_statuses && !is.na(final_clid)
  editor_type <- dplyr::case_when(
    final_status %in% accepted_statuses ~ "single_cl_id",
    final_status %in% c("mixed_identity") ~ "component_set",
    final_status %in% c("ambiguous", "mixed_signal", "review") ~ "candidate_set",
    TRUE ~ "raw_report"
  )
  row <- data.frame(dataset = inst$dataset, cluster_id = inst$cluster_id, reviewer = inst$reviewer,
                    raw_label = inst$raw_label, normalised_label = idx$norm_form(inst$raw_label),
                    provided_cl_id = inst$provided_clid %||% NA_character_,
                    reviewer_confidence = inst$conf,
                    markers = paste(markers, collapse = "|"), tissue = meta$tissue,
                    species = meta$species, condition = meta$condition,
                    candidate_ids = paste(cands$cl_id, collapse = "|"),
                    candidate_labels = paste(cands$canonical_label, collapse = "|"),
                    candidate_scores = paste(round(cands$lexical_score, 4), collapse = "|"),
                    component_cl_ids = if (length(component_clids) > 0L) paste(component_clids, collapse = "|") else NA_character_,
                    component_labels = if (length(component_labels) > 0L) paste(component_labels, collapse = "|") else NA_character_,
                    component_agreement = component_agree,
                    final_status = final_status, decision = decision,
                    cl_id = final_clid, canonical_label = final_label,
                    mapping_method = if (final_status %in% accepted_statuses) "evidence_llm" else "evidence_llm_abstain",
                    specificity_relation = if (!is.na(ma_spec)) ma_spec else NA_character_,
                    conf_cat = conf_cat, reason = reason,
                    proposed_cl_id = proposed_clid, proposed_canonical_label = proposed_label,
                    single_cl_ready = single_ready, adjudication_ready = TRUE,
                    evaluation_ready = single_ready,
                    editor_input_type = editor_type,
                    mapper_a_status = ma_status, mapper_a_rank = ma_rank, mapper_a_cl_id = ma_clid,
                    mapper_a_specificity = ma_spec, mapper_a_confidence = ma_conf,
                    mapper_b_status = mb_status, mapper_b_rank = mb_rank, mapper_b_cl_id = mb_clid,
                    mapper_b_specificity = mb_spec, mapper_b_confidence = mb_conf,
                    ab_clid_agreement = ab_clid_agree, ab_specificity_agreement = ab_spec_agree,
                    verifier_called = verifier_called, verifier_verdict = verifier_verdict,
                    verifier_reason = verifier_reason, verifier_action = verifier_action,
                    hash_a = hash_a, hash_b = hash_b, hash_v = hash_v,
                    source_file = inst$source_file, source_run_id = inst$source_run_id,
                    model = opt$model, prompt_version = PROMPT_VERSION, pipeline_version = PIPELINE_VERSION,
                    ontology_version = "2025-07-30", ontology_checksum = ontology_checksum,
                    temperature = opt$temperature, max_tokens = opt$max_tokens,
                    top_k = TOP_K_CANDIDATES, marker_limit = MARKER_LIMIT,
                    stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <- row
  write_tsv(bind_rows(rows), partial_file)
}

# ============ Complete ============
res <- bind_rows(rows)
write_tsv(res, file.path(out_dir, registry_filename))
if (file.exists(partial_file)) file.remove(partial_file)
cat("\n=== CL-Linker Final Summary ===\n")
print(res %>% count(final_status) %>% mutate(pct = n / sum(n) * 100))
cat("\nOutput: ", file.path(out_dir, "reviewer_mapping_registry.tsv"), "\n")
