#!/usr/bin/env Rscript
# ============================================================
# 08_build_judge_inputs_v4_3_7_1_mapping_safe.R — Judge Input
# Input: 07.5 registry + 07/07b summary + 03b cassia + step1 queries
# Output: final/<ds>/llm_outputs/llm_judge_inputs_v3/*_LLM_JUDGE_INPUT.json
# ============================================================
rm(list = ls())
triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
})

suppressMessages(library(Triage))

`%||%` <- function(a,b) if (!is.null(a)) a else b
suppressMessages(library(Triage))
suppressMessages(library(Triage))
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

MAX_DEGS <- 80L
MAX_ENRICH_TERMS <- 50L
MAX_TOP_GENES <- 15L

truncate_text <- function(x, max_chars = 2500L) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return(NULL)
  x <- stringr::str_squish(x)
  if (nchar(x) <= max_chars) return(x)
  paste0(substr(x, 1, max_chars), "...")
}

as_char_vec <- function(x, k = 5L) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) {
    out <- x
  } else if (is.list(x)) {
    out <- unlist(x, use.names = FALSE)
  } else {
    out <- as.character(x)
  }
  out <- out[!is.na(out) & nzchar(out)]
  out <- unique(out)
  if (length(out) == 0) return(NULL)
  if (length(out) > k) out <- out[1:k]
  as.character(out)
}

as_list_vec <- function(x, k = 5L) {
  v <- as_char_vec(x, k = k)
  if (is.null(v) || length(v) == 0) return(list())
  as.list(v)
}

cl_json <- Sys.getenv("CL_LOCAL_JSON", unset = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"))
cl_cfg <- if (file.exists(cl_json)) make_cl_cfg(cl_json, prefer_ols = FALSE, cache_dir = "") else NULL
cl_graph <- if (!is.null(cl_cfg)) load_cl_graph(cl_cfg) else NULL
cl_idx <- if (file.exists(cl_json)) build_cl_index(cl_json) else NULL

map_label_to_clid <- function(label, cl_cfg) {
  if (is.null(label) || !nzchar(label) || is.null(cl_cfg)) return(NA_character_)
  res <- normalize_cl_three_state(label, "", cl_cfg)
  res$final_clid %||% NA_character_
}

get_onehop_parent_labels <- function(label, cl_graph, cl_cfg) {
  if (is.null(label) || is.null(cl_graph) || is.null(cl_graph$cl)) return(character(0))
  res <- normalize_cl_three_state(label, "", cl_cfg)
  clid <- res$final_clid %||% NA_character_
  if (is.na(clid) || !nzchar(clid)) return(character(0))
  term <- cl_graph$cl[[clid]]
  if (is.null(term)) return(character(0))
  anc <- term$ancestors %||% NULL
  if (is.null(anc)) return(character(0))
  d <- unlist(anc, use.names = TRUE)
  if (length(d) == 0) return(character(0))
  onehop_ids <- names(d[d == 1])
  if (length(onehop_ids) == 0) return(character(0))
  labels <- vapply(onehop_ids, function(id) {
    lab <- cl_graph$cl[[id]]$label %||% NA_character_
    as.character(lab)
  }, character(1))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  unique(labels)
}

build_backoff_candidates <- function(labels, cl_graph, cl_cfg) {
  if (length(labels) == 0) return(list())
  out <- character(0)
  for (lab in labels) {
    out <- c(out, get_onehop_parent_labels(lab, cl_graph, cl_cfg))
  }
  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

extract_report_labels <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(character(0))
  out <- character(0)
  get_lab <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x) && !is.null(x$candidate_cell_type)) return(as.character(x$candidate_cell_type))
    if (is.character(x)) return(as.character(x))
    NULL
  }

  for (key in c("main_type_schema", "main_type", "subtype_level_1_schema", "subtype_level_1")) {
    lab <- get_lab(obj[[key]] %||% NULL)
    if (!is.null(lab)) out <- c(out, lab)
  }

  sub2 <- obj$subtype_level_2_schema %||% obj$subtype_level_2 %||% NULL
  if (is.list(sub2) && is.list(sub2$core_identity)) {
    lab <- get_lab(sub2$core_identity)
    if (!is.null(lab)) out <- c(out, lab)
  }

  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

get_hop_parent_labels <- function(label, cl_graph, cl_cfg, hop = 1L) {
  if (is.null(label) || is.null(cl_graph) || is.null(cl_graph$cl)) return(character(0))
  res <- normalize_cl_three_state(label, "", cl_cfg)
  clid <- res$final_clid %||% NA_character_
  if (is.na(clid) || !nzchar(clid)) return(character(0))
  term <- cl_graph$cl[[clid]]
  if (is.null(term)) return(character(0))
  anc <- term$ancestors %||% NULL
  if (is.null(anc)) return(character(0))
  d <- unlist(anc, use.names = TRUE)
  if (length(d) == 0) return(character(0))
  hop_ids <- names(d[d == hop])
  if (length(hop_ids) == 0) return(character(0))
  labels <- vapply(hop_ids, function(id) {
    lab <- cl_graph$cl[[id]]$label %||% NA_character_
    as.character(lab)
  }, character(1))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  labels <- setdiff(labels, STOP_BACKOFF)
  unique(labels)
}

build_backoff_candidates <- function(labels, cl_graph, cl_cfg, max_hops = 2L) {
  if (length(labels) == 0) return(character(0))
  out <- character(0)
  for (lab in labels) {
    for (h in seq_len(max_hops)) {
      out <- c(out, get_hop_parent_labels(lab, cl_graph, cl_cfg, hop = h))
    }
  }
  out <- out[!is.na(out) & nzchar(out)]
  out <- setdiff(unique(out), STOP_BACKOFF)
  unique(out)
}

STOP_BACKOFF <- c(
  "cell", "native cell", "animal cell",
  "precursor cell", "progenitor cell",
  "neuron", "neural cell"
)

# ============================================================
# CLI (paths follow your structure; do NOT change defaults)
# - NC-ready paradigm:
#   * NO marker_panels in judge inputs (remove dataset-specific heuristics)
#   * Keep low-dimensional lineage scores as a SOFT, reproducible cue
#   * Do NOT neuron-filter our markers (avoid bias / mixed inflation)
# ============================================================
option_list <- list(
  make_option("--cassia_csv", type="character",
              default="",
              help="Path to Cassia annotation_cassia_FINAL_RESULTS.csv. If empty or missing, auto-detect from CASSIA_*_*_*/01_annotation_results/"),
  make_option("--in_house_summary_csv", type="character",
              default=""),
  # ---- INTER ----
  make_option("--enrichment_summary_csv", type="character",
              default=""),
  make_option("--mapping_registry_csv", type="character", default="",
              help="Frozen reviewer mapping registry from Step 07.5. Valid registry CLIDs are preserved; missing/ambiguous free-text candidates may use deterministic local CL-Linker fallback."),
  make_option("--enrichment_mode", type="character", default="csv",
              help="csv | intermediate_outputs | step1_json"),
  make_option("--intermediate_outputs_dir", type="character",
              default=""),
  make_option("--step1_dir", type="character",
              default=""),
  make_option("--exclude_reviewer", type="character", default="",
              help="Ablation: reviewer to exclude (cassia|in_house|enrichment). Empty = full 3-reviewer run."),
  make_option("--out_dir", type="character",
              default=""),
  make_option("--dataset_name", type="character", default="",
              help="Dataset name for default paths (default: basename(getwd()))."),
  make_option("--our_level", type="character", default="main",
              help="main|sub1|sub2 [default %default]"),
  make_option("--cluster_id", type="character", default="",
              help="Optional: only build one cluster_id (must match step1 query filename prefix)"),
  # Keep these flags for backward compatibility, but they are UNUSED in this optimized version
  make_option("--deg_csv", type="character", default="maskdeg.csv",
              help="(Unused in optimized version; kept for compatibility)"),
  make_option("--seurat_qs", type="character", default="",
              help="(Unused in optimized version; kept for compatibility)"),
  make_option("--seurat_cluster_field", type="character", default="",
              help="(Unused in optimized version; kept for compatibility)"),
  make_option("--seurat_assay", type="character", default="",
              help="(Unused in optimized version; kept for compatibility)")
)
opt <- parse_args(OptionParser(option_list=option_list))
dataset_name <- if (nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- get_dataset_config(dataset_name, project_root)
if (!nzchar(opt$in_house_summary_csv)) opt$in_house_summary_csv <- file.path(cfg$llm_outputs_root, "post_summary", "summary.csv")
if (!nzchar(opt$enrichment_summary_csv)) opt$enrichment_summary_csv <- file.path(cfg$llm_outputs_root, "inter", "post_summary", "summary.csv")
if (!nzchar(opt$mapping_registry_csv)) opt$mapping_registry_csv <- file.path(triageHome, "reports", "cl_mapping", "reviewer_mapping_registry.tsv")
if (!nzchar(opt$intermediate_outputs_dir)) opt$intermediate_outputs_dir <- cfg$intermediate_outputs_root
if (!nzchar(opt$step1_dir)) opt$step1_dir <- file.path(cfg$llm_inputs_root, "step1_report_queries")
if (!nzchar(opt$out_dir)) opt$out_dir <- file.path(cfg$llm_outputs_root, "llm_judge_inputs_v3")
ensure_dir(opt$out_dir)

get_final_dir <- function(our_summary_csv) {
  root <- normalizePath(dirname(dirname(our_summary_csv)), mustWork = FALSE)
  file.path(root, "final_passed")
}

get_step1_cluster_ids <- function(step1_dir) {
  if (!dir.exists(step1_dir)) return(character(0))
  fs <- list.files(step1_dir, pattern = "_step1_report_query\\.json$", full.names = FALSE)
  if (length(fs) == 0) return(character(0))
  sub("_step1_report_query\\.json$", "", fs)
}

# ---- INTER ----
read_enrichment_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  cols <- names(df)
  pick <- function(cands) {
    hit <- intersect(cands, cols)
    if (length(hit) == 0) NA_character_ else hit[[1]]
  }
  col_cluster <- pick(c("cluster_id", "cluster", "clusterID"))
  col_label <- pick(c("predicted_label", "cell_type", "label", "predicted_cell_type"))
  col_clid <- pick(c("cl_id", "clid", "cell_ontology_id"))
  col_conf <- pick(c("confidence", "score"))
  col_markers <- pick(c("evidence_markers", "markers"))
  col_terms <- pick(c("evidence_terms", "terms"))
  col_rationale <- pick(c("rationale", "reasoning", "explanation"))
  if (is.na(col_cluster)) return(NULL)
  df %>%
    transmute(
      cluster_id = as.character(.data[[col_cluster]]),
      predicted_label = if (!is.na(col_label)) as.character(.data[[col_label]]) else NA_character_,
      cl_id = if (!is.na(col_clid)) as.character(.data[[col_clid]]) else NA_character_,
      confidence = if (!is.na(col_conf)) suppressWarnings(as.numeric(.data[[col_conf]])) else NA_real_,
      evidence_markers = if (!is.na(col_markers)) as.character(.data[[col_markers]]) else NA_character_,
      evidence_terms = if (!is.na(col_terms)) as.character(.data[[col_terms]]) else NA_character_,
      rationale = if (!is.na(col_rationale)) as.character(.data[[col_rationale]]) else NA_character_
    )
}

extract_inter_terms_from_step1 <- function(step1_obj) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(character(0))
  bio <- step1_obj$bioinfo$analysis_results %||% NULL
  if (is.null(bio) || !is.list(bio)) return(character(0))
  terms <- character(0)
  for (k in names(bio)) {
    v <- bio[[k]]
    if (is.character(v)) terms <- c(terms, v)
    if (is.list(v)) {
      for (item in v) {
        if (is.character(item)) terms <- c(terms, item)
        if (is.list(item)) {
          cand <- item$term %||% item$description %||% item$name %||% item$pathway %||% item$id
          if (!is.null(cand)) terms <- c(terms, as.character(cand))
        }
      }
    }
  }
  terms <- terms[!is.na(terms) & nzchar(terms)]
  unique(terms)
}

extract_inter_markers_from_step1 <- function(step1_obj) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(character(0))
  te <- step1_obj$candidates_for_evaluation$data$Top_Evidence_Genes %||% NULL
  if (is.null(te)) return(character(0))
  if (is.list(te)) te <- unlist(te, recursive = TRUE, use.names = FALSE)
  if (!is.character(te)) return(character(0))
  te <- trimws(te)
  te <- te[nzchar(te)]
  unique(te)
}

load_enrichment_for_cluster <- function(cid, inter_mode, inter_summary_df, step1_obj) {
  if (identical(inter_mode, "csv") && !is.null(inter_summary_df)) {
    hit <- inter_summary_df %>% filter(.data$cluster_id == cid)
    if (nrow(hit) > 0) return(hit[1, ])
  }
  if (identical(inter_mode, "step1_json")) {
    return(tibble::tibble(
      cluster_id = cid,
      predicted_label = NA_character_,
      cl_id = NA_character_,
      confidence = NA_real_,
      evidence_markers = if (length(extract_inter_markers_from_step1(step1_obj)) > 0) paste(extract_inter_markers_from_step1(step1_obj), collapse = "; ") else "",
      evidence_terms = if (length(extract_inter_terms_from_step1(step1_obj)) > 0) paste(extract_inter_terms_from_step1(step1_obj), collapse = "; ") else "",
      rationale = NA_character_
    ))
  }
  if (identical(inter_mode, "intermediate_outputs")) {
    terms <- extract_inter_terms_from_step1(step1_obj)
    return(tibble::tibble(
      cluster_id = cid,
      predicted_label = NA_character_,
      cl_id = NA_character_,
      confidence = NA_real_,
      evidence_markers = "",
      evidence_terms = if (length(terms) > 0) paste(terms, collapse = "; ") else "",
      rationale = NA_character_
    ))
  }
  NULL
}

auto_pick_cassia_csv <- function(step1_dir, dataset_name) {
  candidates <- Sys.glob(file.path(project_root, "final", dataset_name, "CASSIA_*_*_*", "01_annotation_results", "annotation_cassia_FINAL_RESULTS.csv"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return(NA_character_)
  
  expected_ids <- get_step1_cluster_ids(step1_dir)
  expected_ids <- as.character(expected_ids)
  expected_ids <- expected_ids[nzchar(expected_ids)]
  
  score_one <- function(p) {
    df <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || !"True Cell Type" %in% names(df)) return(list(score = -1, n = 0))
    ids <- as.character(df[["True Cell Type"]])
    ids <- ids[!is.na(ids)]
    ids <- stringr::str_trim(ids)
    ids <- ids[nzchar(ids)]
    if (length(expected_ids) == 0) return(list(score = 0, n = length(unique(ids))))
    overlap <- length(intersect(unique(ids), unique(expected_ids)))
    list(score = overlap, n = length(unique(ids)))
  }
  
  scored <- lapply(candidates, score_one)
  scores <- vapply(scored, function(x) x$score, numeric(1))
  mt <- suppressWarnings(as.numeric(file.info(candidates)$mtime))
  pick <- order(scores, mt, decreasing = TRUE)[1]
  candidates[pick]
}

if (!nzchar(opt$cassia_csv) || !file.exists(opt$cassia_csv)) {
  picked <- auto_pick_cassia_csv(opt$step1_dir, dataset_name)
  if (is.na(picked)) {
    stop(
      "cassia_csv not provided/found, and auto-detection failed.\n",
      "Provide --cassia_csv explicitly, or ensure CASSIA_*_*_*/01_annotation_results/annotation_cassia_FINAL_RESULTS.csv exists."
    )
  }
  opt$cassia_csv <- picked
  message("[INFO] auto-picked cassia_csv: ", opt$cassia_csv)
}

SYSTEM_PROMPT <- "You are preparing judge inputs. Keep content concise; do not add instructions."
INSTRUCTION_PROMPT <- "Adjudicate the final cell type using the provided inputs. Return ONLY the final adjudication JSON."

# ============================================================
# Helpers
# ============================================================
split_first_token <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x)) return(NA_character_)
  x <- stringr::str_trim(x)
  if (!nzchar(x)) return(NA_character_)
  x
}

is_bad_main <- function(x) {
  x2 <- as.character(x)
  x2 <- ifelse(is.na(x2), "", x2)
  x2 <- tolower(stringr::str_trim(x2))
  !nzchar(x2) | stringr::str_detect(x2, "mixed|doublet|triplet|multiplet|artifact|unresolved|technical|contamination|sequencing")
}

is_too_generic_for_judge <- function(x) {
  x2 <- as.character(x)
  x2 <- ifelse(is.na(x2), "", x2)
  x2 <- tolower(stringr::str_trim(x2))
  !nzchar(x2) | (x2 %in% c("cell", "unknown", "unassigned"))
}

slice_atomic <- function(x, n) {
  if (is.null(x)) return(x)
  if (!is.atomic(x)) return(x)
  if (length(x) <= n) return(x)
  head(x, n)
}

slice_list_columns <- function(d, n) {
  if (is.null(d) || !is.list(d)) return(d)
  lapply(d, function(v) slice_atomic(v, n))
}

extract_pmids <- function(ev) {
  if (is.null(ev) || !is.list(ev)) return(character(0))
  pmids <- ev$retrieved_articles$articles_db$pmid %||% ev$retrieved_articles$pmid %||% character(0)
  pmids <- as.character(pmids)
  pmids <- pmids[nzchar(pmids)]
  unique(pmids)
}

slim_evidence <- function(cd) {
  ev <- cd$evidence %||% NULL
  if (is.null(ev) || !is.list(ev)) return(NULL)
  pmids <- extract_pmids(ev)
  list(
    description = ev$description %||% NULL,
    retrieved_articles = list(
      articles_db = list(pmid = pmids)
    )
  )
}

extract_term_fields <- function(item) {
  if (!is.list(item)) return(character(0))
  cand <- item$term %||% item$description %||% item$name %||% item$pathway %||% item$id %||% NULL
  if (is.null(cand)) return(character(0))
  as.character(cand)
}

slim_bioinfo <- function(cd, max_terms = MAX_ENRICH_TERMS) {
  bio <- cd$bioinfo %||% NULL
  if (is.null(bio) || !is.list(bio)) return(NULL)
  ar <- bio$analysis_results %||% NULL
  if (is.null(ar) || !is.list(ar)) return(NULL)
  out <- list()
  for (k in names(ar)) {
    v <- ar[[k]]
    if (is.character(v)) {
      out[[k]] <- head(unique(v), max_terms)
      next
    }
    if (is.list(v)) {
      if (length(v) > 0 && is.list(v[[1]]) && !is.atomic(v[[1]])) {
        terms <- unlist(lapply(v, extract_term_fields), use.names = FALSE)
        terms <- terms[nzchar(terms)]
        out[[k]] <- head(unique(terms), max_terms)
        next
      }
      out[[k]] <- lapply(v, function(x) if (is.character(x)) head(unique(x), max_terms) else x)
      next
    }
    out[[k]] <- v
  }
  list(description = bio$description %||% NULL, analysis_results = out)
}

slim_degs <- function(cd, max_rows = MAX_DEGS) {
  degs <- cd$degs %||% NULL
  if (is.null(degs) || !is.list(degs)) return(NULL)
  d <- degs$data %||% NULL
  d <- slice_list_columns(d, max_rows)
  list(description = degs$description %||% NULL, data = d)
}

slim_degs_state <- function(cd, max_rows = MAX_DEGS) {
  degs <- cd$degs_state_support %||% NULL
  if (is.null(degs) || !is.list(degs)) return(NULL)
  d <- degs$data %||% NULL
  d <- slice_list_columns(d, max_rows)
  list(description = degs$description %||% NULL, data = d)
}

extract_top_evidence_genes <- function(top_evidence_genes) {
  if (is.null(top_evidence_genes)) return(character(0))
  genes <- top_evidence_genes
  if (is.list(genes)) genes <- unlist(genes, recursive = TRUE, use.names = FALSE)
  if (!is.character(genes)) return(character(0))
  parts <- unlist(strsplit(genes, "[,;]"), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

build_marker_texts_from_top_genes <- function(top_evidence_genes, max_genes = MAX_TOP_GENES) {
  genes <- extract_top_evidence_genes(top_evidence_genes)
  if (length(genes) == 0) {
    return(list(markers_raw = NULL, markers_no_mito = NULL, genes_raw = character(0), genes_no_mito = character(0)))
  }
  genes_u <- unique(genes)
  raw <- head(genes_u, max_genes)
  no_mito <- raw[!stringr::str_detect(toupper(raw), "^MT-")]
  list(
    markers_raw = if (length(raw) > 0) paste(raw, collapse = ", ") else NULL,
    markers_no_mito = if (length(no_mito) > 0) paste(no_mito, collapse = ", ") else NULL,
    genes_raw = raw,
    genes_no_mito = no_mito
  )
}

get_top_degs <- function(step1_obj, n = 15L) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(NULL)
  gs <- step1_obj$degs$data$geneSymbol %||% character(0)
  gs <- as.character(gs)
  gs <- gs[!is.na(gs) & nzchar(gs)]
  if (length(gs) == 0) return(NULL)
  paste(head(gs, n), collapse = ", ")
}

get_top_terms <- function(step1_obj, n = 12L) {
  if (is.null(step1_obj) || !is.list(step1_obj)) return(NULL)
  bp <- step1_obj$bioinfo$analysis_results$biological_processes %||% character(0)
  cc <- step1_obj$bioinfo$analysis_results$cellular_components %||% character(0)
  kegg <- step1_obj$bioinfo$analysis_results$kegg_pathways %||% character(0)
  re <- step1_obj$bioinfo$analysis_results$reactome_pathways %||% character(0)
  if (is.list(re)) re <- unlist(re, use.names = FALSE) else re <- as.character(re)
  terms <- unique(c(bp, cc, kegg, re))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0) return(NULL)
  paste(head(terms, n), collapse = "; ")
}

extract_fallback_terms <- function(inter_summary) {
  if (is.null(inter_summary) || !is.list(inter_summary)) return(NULL)
  raw <- c(
    inter_summary$evidence_terms_method %||% NULL,
    inter_summary$evidence_terms %||% NULL
  )
  if (is.null(raw)) return(NULL)
  if (is.list(raw)) raw <- unlist(raw, recursive = TRUE, use.names = FALSE)
  raw <- as.character(raw)
  raw <- unlist(strsplit(raw, "[;|,]"), use.names = FALSE)
  raw <- trimws(raw)
  raw <- raw[!is.na(raw) & nzchar(raw)]
  if (length(raw) == 0) return(NULL)
  paste(head(unique(raw), 12L), collapse = "; ")
}

extract_step15_pmids <- function(step15_obj) {
  pm <- step15_obj$pmids %||% character(0)
  pm <- as.character(pm)
  pm <- pm[!is.na(pm) & nzchar(pm)]
  unique(pm)
}

load_final_for_cluster <- function(cid) {
  dir_final <- get_final_dir(opt$in_house_summary_csv)
  if (!dir.exists(dir_final)) return(NULL)
  f <- file.path(dir_final, paste0(cid, "_FINAL_passed.json"))
  if (!file.exists(f)) return(NULL)
  tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

extract_overall_assessment <- function(sec) {
  if (is.null(sec) || !is.list(sec)) return(NULL)
  sr <- sec$structured_reasoning %||% NULL
  if (is.null(sr) || !is.list(sr)) return(NULL)
  oa <- sr$overall_assessment %||% NULL
  if (is.null(oa)) return(NULL)
  oa <- as.character(oa)
  oa <- oa[nzchar(oa)]
  if (length(oa) == 0) return(NULL)
  oa[[1]]
}

build_final_reasoning_text <- function(final_obj, max_chars = 800L) {
  if (is.null(final_obj) || !is.list(final_obj)) return(NULL)
  parts <- c(
    extract_overall_assessment(final_obj$main_type_schema %||% NULL),
    extract_overall_assessment(final_obj$subtype_level_1_schema %||% NULL),
    extract_overall_assessment(final_obj$subtype_level_2_schema %||% NULL)
  )
  parts <- parts[!is.na(parts)]
  parts <- parts[nzchar(parts)]
  parts <- unique(parts)
  if (length(parts) == 0) return(NULL)
  out <- paste(parts, collapse = " | ")
  out <- stringr::str_trim(out)
  if (!nzchar(out)) return(NULL)
  if (nchar(out) > max_chars) out <- substr(out, 1, max_chars)
  out
}

compute_prefix_counts <- function(genes) {
  if (is.null(genes) || length(genes) == 0) return(list())
  genes_u <- unique(toupper(as.character(genes)))
  prefixes <- c("MT-", "RPL", "RPS", "HSP", "IGH", "IGK", "IGL", "ALB", "HBB", "HBA")
  out <- list()
  for (p in prefixes) {
    out[[p]] <- sum(stringr::str_detect(genes_u, paste0("^", p)))
  }
  out
}

count_enrichment_terms <- function(step1_obj) {
  bio <- step1_obj$bioinfo$analysis_results %||% NULL
  if (is.null(bio) || !is.list(bio)) return(0L)
  terms <- character(0)
  for (k in names(bio)) {
    v <- bio[[k]]
    if (is.character(v)) {
      terms <- c(terms, v)
    } else if (is.list(v)) {
      if (length(v) > 0 && is.list(v[[1]]) && !is.atomic(v[[1]])) {
        terms <- c(terms, unlist(lapply(v, extract_term_fields), use.names = FALSE))
      } else {
        for (item in v) if (is.character(item)) terms <- c(terms, item)
      }
    }
  }
  terms <- terms[nzchar(as.character(terms))]
  length(unique(terms))
}

count_degs <- function(step1_obj) {
  g <- step1_obj$degs$data$geneSymbol %||% character(0)
  g <- as.character(g)
  g <- g[nzchar(g)]
  length(unique(g))
}

# ============================================================
# Load Step1 dossier (minimize tokens, keep audit pointers)
# ============================================================
load_step1 <- function(cid) {
  f <- file.path(opt$step1_dir, paste0(cid, "_step1_report_query.json"))
  if (!file.exists(f)) return(NULL)
  x <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  cd <- x$input_data$cluster_dossier %||% NULL
  if (is.null(cd) || !is.list(cd)) return(NULL)
  
  # keep only minimal candidate fields
  if (is.list(cd$candidates_for_evaluation) && is.list(cd$candidates_for_evaluation$data)) {
    d <- cd$candidates_for_evaluation$data
    keep <- c("Candidate_Cell_Type", "Cell_Ontology_ID", "Final_Score", "Top_Evidence_Genes")
    keep <- keep[keep %in% names(d)]
    cd$candidates_for_evaluation$data <- if (length(keep) > 0) d[keep] else d
  }
  
  cd$degs <- slim_degs(cd)
  if (!is.null(cd$degs_state_support)) cd$degs_state_support <- slim_degs_state(cd)
  cd$evidence <- slim_evidence(cd)
  cd$bioinfo <- slim_bioinfo(cd)
  cd
}

normalize_ct_label <- function(x) {
  if (is.null(x)) return(NA_character_)
  x2 <- as.character(x)[1]
  if (is.na(x2)) return(NA_character_)
  x2 <- stringr::str_trim(x2)
  x2 <- stringr::str_replace_all(x2, "\\s+", " ")
  x2 <- tolower(x2)
  if (!nzchar(x2)) NA_character_ else x2
}

# v2: prioritize CL-ID matching, avoid permissive token matching, Return score + provenance
# Match order: 1) fixed CL ID ∈ step1 Cell_Ontology_ID  2) normalized exact label
# 3) synonym/alias exact 4) unavailable (without imputation)
candidate_score_match_v2 <- function(step1_obj, label, clid = NA_character_) {
  if (is.null(step1_obj) || !is.list(step1_obj)) {
    return(list(score = NA_real_, source = "unavailable", match_type = "no_step1"))
  }
  d <- step1_obj$candidates_for_evaluation$data %||% NULL
  if (is.null(d) || !is.list(d)) {
    return(list(score = NA_real_, source = "unavailable", match_type = "no_candidate_table"))
  }
  # data may be a column-oriented dictionary (columns) or a row-oriented list ( list)
  if (!is.null(d$Candidate_Cell_Type)) {
    cand <- as.character(d$Candidate_Cell_Type %||% NULL)
    sc <- suppressWarnings(as.numeric(d$Final_Score %||% NULL))
    cand_clids <- as.character(d$Cell_Ontology_ID %||% NULL)
  } else if (length(d) > 0 && is.list(d[[1]])) {
    cand <- vapply(d, function(row) as.character(row$Candidate_Cell_Type %||% NA_character_), character(1))
    sc <- suppressWarnings(vapply(d, function(row) as.numeric(row$Final_Score %||% NA_real_), numeric(1)))
    cand_clids <- vapply(d, function(row) as.character(row$Cell_Ontology_ID %||% NA_character_), character(1))
  } else {
    return(list(score = NA_real_, source = "unavailable", match_type = "unrecognized_data_shape"))
  }
  if (length(cand) == 0 || length(sc) == 0) {
    return(list(score = NA_real_, source = "unavailable", match_type = "no_score_column"))
  }

  # 1) CL ID exact match (prioritized, consistent with the fixed registry)
  if (!is.na(clid) && nzchar(clid) && length(cand_clids) > 0) {
    clid_clean <- trimws(clid)
    hit <- which(!is.na(cand_clids) & trimws(cand_clids) == clid_clean)
    if (length(hit) > 0) {
      return(list(score = max(sc[hit], na.rm = TRUE), source = "step1_exact_clid",
                  match_type = "exact_clid"))
    }
  }

  # 2) normalized exact label
  lab <- normalize_ct_label(label)
  if (!is.na(lab)) {
    cand_n <- vapply(cand, normalize_ct_label, character(1))
    hit <- which(!is.na(cand_n) & cand_n == lab)
    if (length(hit) > 0) {
      return(list(score = max(sc[hit], na.rm = TRUE), source = "step1_exact_label",
                  match_type = "exact_label"))
    }
  }

  # 3) unavailable: without imputation, do not apply permissive matching
  list(score = NA_real_, source = "unavailable", match_type = "no_match")
}

candidate_score_for_label <- function(step1_obj, label) {
  r <- candidate_score_match_v2(step1_obj, label, NA_character_)
  r$score
}

# Validation helper: check label-CLID consistency
validate_label_clid <- function(label, clid, cl_cfg, method_name, cluster_id) {
  if (is.null(label) || is.na(label) || !nzchar(label)) {
    warning(sprintf("[VALIDATE] %s cluster %s: missing label", method_name, cluster_id))
    return(list(valid = FALSE, corrected_clid = NA_character_))
  }
  
  # Map label to expected CLID
  expected_clid <- map_label_to_clid(label, cl_cfg)
  
  # If input CLID is missing, use the mapped one
  if (is.null(clid) || is.na(clid) || !nzchar(clid)) {
    if (!is.na(expected_clid)) {
      cat(sprintf("[VALIDATE] %s cluster %s: missing CLID, mapped from label '%s' -> %s\n", 
                  method_name, cluster_id, label, expected_clid))
      return(list(valid = TRUE, corrected_clid = expected_clid))
    } else {
      warning(sprintf("[VALIDATE] %s cluster %s: cannot map label '%s' to CLID", 
                      method_name, cluster_id, label))
      return(list(valid = FALSE, corrected_clid = NA_character_))
    }
  }
  
  # If input CLID doesn't match expected, warn and use label-mapped CLID
  if (!is.na(expected_clid) && clid != expected_clid) {
    cat(sprintf("[VALIDATE] %s cluster %s: CLID mismatch! label='%s', input_CLID=%s, expected_CLID=%s. Using expected.\n",
                method_name, cluster_id, label, clid, expected_clid))
    return(list(valid = FALSE, corrected_clid = expected_clid))
  }
  
  return(list(valid = TRUE, corrected_clid = clid))
}

pick_our_top1_by_candidate_score <- function(step1_obj, our_main_raw, our_sub1_raw, our_sub2_raw, preferred_level = c("main","sub1","sub2")) {
  preferred_level <- match.arg(preferred_level)
  main_first <- split_first_token(our_main_raw)
  sub1_first <- split_first_token(our_sub1_raw)
  sub2_first <- split_first_token(our_sub2_raw)
  
  s_main <- candidate_score_for_label(step1_obj, main_first)
  s_sub1 <- candidate_score_for_label(step1_obj, sub1_first)
  s_sub2 <- candidate_score_for_label(step1_obj, sub2_first)
  
  if (all(is.na(c(s_main, s_sub1, s_sub2)))) {
    if (preferred_level == "sub2" && !is.na(sub2_first) && nzchar(sub2_first)) return(list(label=sub2_first, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
    if (preferred_level == "sub1" && !is.na(sub1_first) && nzchar(sub1_first)) return(list(label=sub1_first, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
    if (!is.na(main_first) && nzchar(main_first)) return(list(label=main_first, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
    if (!is.na(sub1_first) && nzchar(sub1_first)) return(list(label=sub1_first, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
    if (!is.na(sub2_first) && nzchar(sub2_first)) return(list(label=sub2_first, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
    return(list(label=NA_character_, scores=list(main=s_main, sub1=s_sub1, sub2=s_sub2)))
  }
  
  scores <- c(main = s_main, sub1 = s_sub1, sub2 = s_sub2)
  best <- names(scores)[which.max(replace(scores, is.na(scores), -Inf))][1]
  tied <- names(scores)[which(replace(scores, is.na(scores), -Inf) == max(replace(scores, is.na(scores), -Inf)))]
  if (length(tied) > 1) {
    if (preferred_level %in% tied) best <- preferred_level
    else if ("sub2" %in% tied) best <- "sub2"
    else if ("sub1" %in% tied) best <- "sub1"
    else best <- tied[[1]]
  }
  
  lab <- if (best == "sub2") sub2_first else if (best == "sub1") sub1_first else main_first
  list(label = lab, scores = as.list(scores))
}

build_topk_labels <- function(x, k = 3L) {
  if (is.null(x)) return(character(0))
  if (is.list(x) && !is.character(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- stringr::str_trim(x)
  x <- x[nzchar(x)]
  x <- unique(x)
  if (length(x) == 0) return(character(0))
  head(x, as.integer(k))
}

build_cassia_top3 <- function(main_raw, sub_raw, mix_raw) {
  cand <- c(split_first_token(main_raw), split_first_token(sub_raw), split_first_token(mix_raw))
  cand <- cand[!is.na(cand)]
  cand <- cand[nzchar(cand)]
  cand <- cand[!vapply(cand, is_bad_main, logical(1))]
  build_topk_labels(unique(cand), 3L)
}

build_our_top3_by_candidate_score <- function(step1_obj, our_main_raw, our_sub1_raw, our_sub2_raw) {
  main_first <- split_first_token(our_main_raw)
  sub1_first <- split_first_token(our_sub1_raw)
  sub2_first <- split_first_token(our_sub2_raw)
  labs <- c(main = main_first, sub1 = sub1_first, sub2 = sub2_first)
  labs <- labs[!is.na(labs)]
  labs <- labs[nzchar(labs)]
  if (length(labs) == 0) return(list(labels = character(0), scores = list()))
  
  sc <- lapply(labs, function(l) candidate_score_for_label(step1_obj, l))
  scv <- suppressWarnings(as.numeric(unlist(sc)))
  names(scv) <- names(labs)
  ord <- order(replace(scv, is.na(scv), -Inf), decreasing = TRUE)
  labs_ord <- unname(labs[ord])
  labs_ord <- build_topk_labels(labs_ord, 3L)
  
  out_scores <- list(
    main = candidate_score_for_label(step1_obj, main_first),
    sub1 = candidate_score_for_label(step1_obj, sub1_first),
    sub2 = candidate_score_for_label(step1_obj, sub2_first)
  )
  list(labels = labs_ord, scores = out_scores)
}

build_candidate_details <- function(step1_obj, cassia_top1, cassia_topk, our_top1, our_topk,
                                 enrich_top1, enrich_topk,
                                 clid_map = list()) {
  entries <- list()
  # clid_map: named list, method -> CL ID (frozen registry). Falls back to label match if absent.
  add_entries <- function(method, top1, topk) {
    labels <- c(as_char_vec(top1, k = 1L) %||% character(0), as_char_vec(topk, k = 3L) %||% character(0))
    labels <- labels[!is.na(labels) & nzchar(labels)]
    if (length(labels) == 0) return(NULL)
    labels <- unique(labels)
    method_clid <- clid_map[[method]] %||% NA_character_
    for (i in seq_along(labels)) {
      lab <- labels[[i]]
      m2 <- candidate_score_match_v2(step1_obj, lab, method_clid)
      entries[[length(entries) + 1]] <<- list(
        label = lab,
        method = method,
        rank = i,
        candidate_score = m2$score,
        candidate_score_available = !is.na(m2$score),
        candidate_score_source = m2$source
      )
    }
    NULL
  }
  add_entries("cassia", cassia_top1, cassia_topk)
  add_entries("our", our_top1, our_topk)
  add_entries("enrich", enrich_top1, enrich_topk)
  entries
}

# ============================================================
# Step1.5 citation payload (unchanged; optional)
# ============================================================
infer_out_root <- function() {
  normalizePath(dirname(dirname(opt$in_house_summary_csv)), mustWork = FALSE)
}

load_step15_pmid_payload <- function(cid, round_i = 1L) {
  out_root <- infer_out_root()
  step15_dir <- file.path(out_root, "step1.5_citation_fix_outputs")
  p <- file.path(step15_dir, paste0(cid, sprintf("_round%d_step1.5_citation_fix_output.json", round_i)))
  if (!file.exists(p)) return(NULL)
  o <- tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(o) || !is.list(o)) return(NULL)
  
  txt <- tryCatch(paste0(readLines(p, warn = FALSE), collapse = "\n"), error = function(e) "")
  pmids <- unique(stringr::str_match_all(txt, regex("pmid\\s*:?\\s*(\\d{6,10})", ignore_case = TRUE))[[1]][,2])
  pmids <- pmids[!is.na(pmids)]
  pmids <- pmids[nzchar(pmids)]
  pmids <- sort(unique(pmids))
  
  sr <- o$main_type$structured_reasoning %||% NULL
  sr_keep <- NULL
  if (is.list(sr)) {
    sr_keep <- list(
      overall_assessment = sr$overall_assessment %||% NULL,
      lineage_identity = sr$lineage_identity %||% NULL,
      functional_capability = sr$functional_capability %||% NULL,
      cellular_state = sr$cellular_state %||% NULL
    )
  }
  
  list(
    pmids = if (length(pmids) > 0) pmids else NULL,
    main_type_structured_reasoning = sr_keep
  )
}

build_reasoning_text_from_step15 <- function(step15_obj, max_items = 5L) {
  if (is.null(step15_obj) || is.null(step15_obj$main_type)) return(NULL)
  sr <- step15_obj$main_type$structured_reasoning %||% NULL
  if (is.null(sr)) return(NULL)
  parts <- c(sr$lineage_identity %||% character(0), sr$functional_capability %||% character(0))
  if (!is.null(sr$overall_assessment)) parts <- c(sr$overall_assessment, parts)
  if (!is.character(parts)) return(NULL)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NULL)
  parts <- head(unique(parts), max_items)
  paste(parts, collapse = " | ")
}

# ============================================================
# Frozen registry helpers
# - Reviewer TOP1 CL IDs remain authoritative from the frozen registry.
# - Reviewer TOP-K labels are mapped ONCE here against the same local frozen CL graph
#   and serialized as `frozen_reviewer_candidates`.
# - Step 09 must consume these frozen candidate CLIDs and must not remap TOP-K text.
# ============================================================
safe_text_registry <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(NA_character_)
  value <- trimws(as.character(x[1]))
  if (!nzchar(value) || toupper(value) == "NA" || tolower(value) == "null") return(NA_character_)
  value
}

split_registry_field <- function(x) {
  value <- safe_text_registry(x)
  if (is.na(value)) return(character(0))
  out <- trimws(unlist(strsplit(value, "\\|", perl = TRUE), use.names = FALSE))
  unique(out[!is.na(out) & nzchar(out)])
}

registry_bool <- function(x) {
  if (is.logical(x) && length(x) > 0L && !is.na(x[1])) return(isTRUE(x[1]))
  value <- tolower(safe_text_registry(x))
  !is.na(value) && value %in% c("true", "t", "1", "yes")
}

registry_record <- function(registry_df, cid, reviewer_id) {
  # Ablation: excluded reviewer has no record -> return NULL (caller must handle)
  if (nzchar(exclude_reviewer) && identical(reviewer_id, exclude_reviewer)) return(NULL)
  hit <- registry_df %>% filter(.data$cluster_id == cid, .data$reviewer == reviewer_id)
  if (nrow(hit) != 1L) stop("Registry must contain exactly one record for ", dataset_name, " / ", cid, " / ", reviewer_id)
  as.list(hit[1, , drop = FALSE])
}

registry_mapping_bundle <- function(rec) {
  if (is.null(rec)) return(NULL)
  single_ready <- registry_bool(rec$single_cl_ready)
  final_clid <- if (single_ready) safe_text_registry(rec$cl_id) else NA_character_
  list(
    reviewer = safe_text_registry(rec$reviewer),
    raw_label = safe_text_registry(rec$raw_label),
    canonical_label = safe_text_registry(rec$canonical_label),
    mapping_status = safe_text_registry(rec$final_status),
    mapping_decision = safe_text_registry(rec$decision),
    mapping_reason = safe_text_registry(rec$reason),
    mapping_method = safe_text_registry(rec$mapping_method),
    specificity_relation = safe_text_registry(rec$specificity_relation),
    editor_input_type = safe_text_registry(rec$editor_input_type),
    single_cl_ready = single_ready,
    adjudication_ready = registry_bool(rec$adjudication_ready),
    evaluation_ready = registry_bool(rec$evaluation_ready),
    cell_ontology_id = final_clid,
    candidate_cl_ids = as.list(split_registry_field(rec$candidate_ids)),
    candidate_labels = as.list(split_registry_field(rec$candidate_labels)),
    component_cl_ids = as.list(split_registry_field(rec$component_cl_ids)),
    component_labels = as.list(split_registry_field(rec$component_labels)),
    proposed_cl_id = safe_text_registry(rec$proposed_cl_id),
    ontology_version = safe_text_registry(rec$ontology_version),
    ontology_checksum = safe_text_registry(rec$ontology_checksum),
    registry_pipeline_version = safe_text_registry(rec$pipeline_version)
  )
}


# Build a frozen ontology record for each reviewer candidate that can enter Head/Chief
# ontology context. This is deliberately performed in Step 08, never inside Step 09.
valid_local_clid_08 <- function(clid) {
  id <- as.character(clid %||% NA_character_)
  !is.na(id) && nzchar(id) && !is.null(cl_graph) && !is.null(cl_graph$cl[[id]])
}

match_registry_candidate_clid <- function(label, mapping_bundle) {
  if (is.null(mapping_bundle) || !is.list(mapping_bundle)) return(NA_character_)
  labs <- as.character(unlist(mapping_bundle$candidate_labels %||% list(), recursive = TRUE, use.names = FALSE))
  ids  <- as.character(unlist(mapping_bundle$candidate_cl_ids %||% list(), recursive = TRUE, use.names = FALSE))
  comp_labs <- as.character(unlist(mapping_bundle$component_labels %||% list(), recursive = TRUE, use.names = FALSE))
  comp_ids  <- as.character(unlist(mapping_bundle$component_cl_ids %||% list(), recursive = TRUE, use.names = FALSE))
  labs <- c(labs, comp_labs)
  ids <- c(ids, comp_ids)
  if (length(labs) == 0L || length(ids) == 0L || length(labs) != length(ids)) return(NA_character_)
  q <- normalize_ct_label(label)
  if (is.na(q)) return(NA_character_)
  nl <- vapply(labs, normalize_ct_label, character(1))
  hit <- which(!is.na(nl) & nl == q & vapply(ids, valid_local_clid_08, logical(1)))
  if (length(hit) != 1L) return(NA_character_)
  ids[[hit[[1]]]]
}

map_topk_label_once_08 <- function(label, mapping_bundle) {
  # 1) Preserve a valid frozen registry candidate mapping when available.
  rid <- match_registry_candidate_clid(label, mapping_bundle)
  if (!is.na(rid) && nzchar(rid)) {
    return(list(
      clid = rid,
      mapping_source = "frozen_registry_candidate",
      mapping_status = "registry_candidate",
      map_quality = 3L
    ))
  }

  # 2) Otherwise use the shared CL-Linker exactly once in Step 08.
  #    Step 09 must never remap this free-text candidate.
  if (is.null(cl_idx)) {
    return(list(
      clid = NA_character_,
      mapping_source = "step08_unmapped",
      mapping_status = "unmapped",
      map_quality = 0L
    ))
  }

  m <- tryCatch(
    cl_link(label, NA_character_, cl_idx),
    error = function(e) NULL
  )
  cid <- as.character(m$cl_id %||% NA_character_)
  if (!valid_local_clid_08(cid)) {
    return(list(
      clid = NA_character_,
      mapping_source = "step08_unmapped",
      mapping_status = as.character(m$mapping_status %||% "unmapped"),
      map_quality = 0L
    ))
  }

  list(
    clid = cid,
    mapping_source = paste0("step08_cl_link:", as.character(m$mapping_method %||% "mapped")),
    mapping_status = as.character(m$mapping_status %||% "mapped"),
    map_quality = 2L
  )
}

build_frozen_reviewer_candidates <- function(method, top1_label, topk_labels,
                                             top1_clid, mapping_bundle,
                                             step1_obj, max_k = 3L) {
  labs <- c(
    as_char_vec(top1_label, k = 1L) %||% character(0),
    as_char_vec(topk_labels, k = max_k) %||% character(0)
  )
  labs <- labs[!is.na(labs) & nzchar(trimws(labs))]
  if (length(labs) == 0L) return(list())

  # Preserve reviewer order/rank. Ontology mappability must never promote a lower
  # reviewer candidate into TOP1.
  key <- vapply(labs, normalize_ct_label, character(1))
  labs <- labs[!duplicated(key)]
  if (length(labs) > max_k) labs <- labs[seq_len(max_k)]

  out <- list()
  for (i in seq_along(labs)) {
    lab <- labs[[i]]

    if (i == 1L && valid_local_clid_08(top1_clid)) {
      # Valid registry TOP1 remains authoritative.
      mp <- list(
        clid = as.character(top1_clid),
        mapping_source = "frozen_registry_top1",
        mapping_status = as.character(mapping_bundle$mapping_status %||% "registry_top1"),
        map_quality = 3L
      )
    } else {
      # If registry TOP1 is unmapped/ambiguous, try the SAME submitted TOP1 text
      # through the shared deterministic CL-Linker. Lower ranks remain lower ranks.
      mp <- map_topk_label_once_08(lab, mapping_bundle)
      if (i == 1L && valid_local_clid_08(mp$clid)) {
        mp$mapping_source <- paste0("top1_", mp$mapping_source)
      }
    }

    sc <- candidate_score_for_label(step1_obj, lab)
    out[[length(out) + 1L]] <- list(
      method = as.character(method),
      rank = as.integer(i),
      label = as.character(lab),
      clid = if (valid_local_clid_08(mp$clid)) as.character(mp$clid) else NULL,
      mapping_source = as.character(mp$mapping_source),
      mapping_status = as.character(mp$mapping_status),
      map_quality = as.integer(mp$map_quality),
      candidate_score = if (is.na(sc)) NULL else as.numeric(sc),
      candidate_score_available = !is.na(sc),
      ontology_version = as.character(mapping_bundle$ontology_version %||% ""),
      ontology_checksum = as.character(mapping_bundle$ontology_checksum %||% "")
    )
  }
  out
}

assert_frozen_candidate_top1 <- function(records, method, expected_clid, cluster_id) {
  hit <- Filter(function(x) {
    identical(as.character(x$method %||% ""), method) &&
      as.integer(x$rank %||% 99L) == 1L
  }, records)
  if (length(hit) != 1L) {
    stop("Frozen candidate integrity: expected exactly one TOP1 for ", cluster_id, " / ", method)
  }

  # A valid registry CLID must never be changed. If registry TOP1 was ambiguous/
  # unmapped, a NULL CLID is allowed and a deterministic CL-Linker fallback may map
  # the SAME submitted TOP1 text.
  if (valid_local_clid_08(expected_clid)) {
    got <- as.character(hit[[1]]$clid %||% NA_character_)
    if (!identical(got, as.character(expected_clid))) {
      stop("Frozen candidate integrity: TOP1 CLID changed for ", cluster_id, " / ", method,
           ": registry=", expected_clid, " frozen_candidate=", got)
    }
  }
  invisible(TRUE)
}

# ============================================================
# Load source summaries for evidence/provenance, and frozen registry for CL IDs
# ============================================================
if (!file.exists(opt$mapping_registry_csv)) stop("mapping_registry_csv not found: ", opt$mapping_registry_csv)
registry_all <- readr::read_tsv(opt$mapping_registry_csv, show_col_types = FALSE)
registry_required <- c(
  "dataset", "cluster_id", "reviewer", "raw_label", "final_status", "cl_id",
  "single_cl_ready", "adjudication_ready", "evaluation_ready", "editor_input_type",
  "candidate_ids", "candidate_labels", "component_cl_ids", "component_labels",
  "canonical_label", "decision", "reason", "mapping_method", "specificity_relation",
  "proposed_cl_id", "ontology_version", "ontology_checksum", "pipeline_version"
)
missing_registry_cols <- setdiff(registry_required, names(registry_all))
if (length(missing_registry_cols) > 0L) stop("Registry missing columns: ", paste(missing_registry_cols, collapse = ", "))
registry_df <- registry_all %>%
  mutate(dataset = as.character(.data$dataset), cluster_id = as.character(.data$cluster_id), reviewer = as.character(.data$reviewer)) %>%
  filter(.data$dataset == dataset_name)
if (nrow(registry_df) == 0L) stop("No registry records for dataset: ", dataset_name)
# ---- Ablation: exclude reviewer ----
exclude_reviewer <- trimws(opt$exclude_reviewer)
if (nzchar(exclude_reviewer)) {
  if (!exclude_reviewer %in% c("cassia", "in_house", "enrichment")) {
    stop("Invalid --exclude_reviewer: ", exclude_reviewer, " (must be cassia|in_house|enrichment)")
  }
  cat("[ABLATION] Excluding reviewer: ", exclude_reviewer, "\n", sep = "")
  registry_df <- registry_df %>% filter(.data$reviewer != exclude_reviewer)
}

if (!setequal(unique(registry_df$reviewer), c("cassia", "in_house", "enrichment"))) {
  if (nzchar(exclude_reviewer)) {
    expected <- setdiff(c("cassia", "in_house", "enrichment"), exclude_reviewer)
    if (!setequal(unique(registry_df$reviewer), expected)) {
      stop("After excluding ", exclude_reviewer, ", registry reviewer set = ",
           paste(unique(registry_df$reviewer), collapse = ","), " expected ", paste(expected, collapse = ","))
    }
  } else {
    stop("Registry reviewer set must be cassia/in_house/enrichment for dataset: ", dataset_name)
  }
}
registry_dup <- registry_df %>% count(cluster_id, reviewer) %>% filter(n != 1L)
if (nrow(registry_dup) > 0L) stop("Duplicate/missing registry keys:\n", paste(capture.output(print(registry_dup)), collapse = "\n"))
registry_counts <- registry_df %>% count(cluster_id)
expected_n_reviewers <- if (nzchar(exclude_reviewer)) 2L else 3L
if (any(registry_counts$n != expected_n_reviewers)) {
  stop("Each cluster must have exactly ", expected_n_reviewers,
       " registry reviewer records (exclude=", exclude_reviewer, ")")
}
checksum_values <- unique(registry_df$ontology_checksum[!is.na(registry_df$ontology_checksum)])
if (length(checksum_values) != 1L) stop("Registry contains mixed ontology checksums")
current_ontology_checksum <- unname(tools::md5sum(cl_json))
if (!identical(as.character(checksum_values[1]), as.character(current_ontology_checksum))) {
  stop("Registry ontology checksum does not match Step 08 ontology: registry=", checksum_values[1],
       " current=", current_ontology_checksum)
}

# Integrity checks: Step 08 must never silently coerce a registry abstention into a single CL ID.
for (j in seq_len(nrow(registry_df))) {
  rec <- registry_df[j, , drop = FALSE]
  single_ready <- registry_bool(rec$single_cl_ready)
  clid <- safe_text_registry(rec$cl_id)
  status <- safe_text_registry(rec$final_status)
  components <- split_registry_field(rec$component_cl_ids)
  if (single_ready && is.na(clid)) stop("single_cl_ready row has no cl_id: ", rec$cluster_id, " / ", rec$reviewer)
  if (!single_ready && !is.na(clid)) stop("non-single registry row carries released cl_id: ", rec$cluster_id, " / ", rec$reviewer)
  if (identical(status, "mixed_identity") && length(components) < 2L) {
    stop("mixed_identity row lacks >=2 component CL IDs: ", rec$cluster_id, " / ", rec$reviewer)
  }
}

cassia_raw <- readr::read_csv(opt$cassia_csv, show_col_types = FALSE)
in_house_sum <- readr::read_csv(opt$in_house_summary_csv, show_col_types = FALSE)

# ---- ENRICH option aliases retained for input compatibility ----
enrich_summary_df <- NULL
if (identical(tolower(opt$enrichment_mode), "csv")) {
  enrich_summary_df <- read_enrichment_summary(opt$enrichment_summary_csv)
}

req_cassia <- c("True Cell Type", "Predicted Main Cell Type")
miss_c <- setdiff(req_cassia, names(cassia_raw))
if (length(miss_c) > 0L) stop("Cassia csv missing columns: ", paste(miss_c, collapse = ", "))
if (!"cluster_id" %in% names(in_house_sum)) stop("In-house summary.csv missing column: cluster_id")

our_level <- tolower(opt$our_level)  # legacy CLI name retained
our_col <- if (our_level == "sub1") "subtype1_pred" else if (our_level == "sub2") "subtype2_core_pred" else "main_pred"
if (!our_col %in% names(in_house_sum)) stop("In-house summary.csv missing column: ", our_col)

# Source tables are used only for rationale/top-k evidence. Registry is authoritative for raw label + CL mapping.
cassia_source <- cassia_raw %>%
  transmute(
    cluster_id = as.character(.data[["True Cell Type"]]),
    cassia_source_main = as.character(.data[["Predicted Main Cell Type"]] %||% ""),
    cassia_source_sub = as.character(.data[["Predicted Sub Cell Types"]] %||% ""),
    cassia_source_mixed = as.character(.data[["Possible Mixed Cell Types"]] %||% "")
  )

cassia_full <- cassia_raw %>% mutate(cluster_id = as.character(.data[["True Cell Type"]]))
get_cassia_raw_inputs <- function(cid) {
  rows <- cassia_full %>% filter(.data$cluster_id == cid)
  if (nrow(rows) == 0L) return(NULL)
  r1 <- rows[1, , drop = FALSE]
  keep_cols <- intersect(c("Predicted Main Cell Type", "Predicted Sub Cell Types", "Possible Mixed Cell Types",
                           "Marker Number", "Marker List", "Score", "Scoring_Reasoning"), names(r1))
  as.list(if (length(keep_cols) > 0L) r1[, keep_cols, drop = FALSE] else r1)
}

in_house_sum2 <- in_house_sum %>%
  mutate(
    main_first = if ("main_pred" %in% names(.)) vapply(.data$main_pred, split_first_token, character(1)) else NA_character_,
    sub1_first = if ("subtype1_pred" %in% names(.)) vapply(.data$subtype1_pred, split_first_token, character(1)) else NA_character_,
    sub2_first = if ("subtype2_core_pred" %in% names(.)) vapply(.data$subtype2_core_pred, split_first_token, character(1)) else NA_character_
  )

in_house_source <- in_house_sum2 %>%
  transmute(
    cluster_id = as.character(.data$cluster_id),
    our_main = if ("main_pred" %in% names(in_house_sum2)) as.character(.data$main_pred) else "",
    our_sub1 = if ("subtype1_pred" %in% names(in_house_sum2)) as.character(.data$subtype1_pred) else "",
    our_sub2 = if ("subtype2_core_pred" %in% names(in_house_sum2)) as.character(.data$subtype2_core_pred) else ""
  )

all_clusters <- sort(unique(registry_df$cluster_id))
if (nzchar(opt$cluster_id)) {
  all_clusters <- all_clusters[all_clusters == opt$cluster_id]
  if (length(all_clusters) == 0L) stop("cluster_id not found in registry: ", opt$cluster_id)
}
idx <- tibble::tibble(cluster_id = all_clusters) %>%
  left_join(cassia_source, by = "cluster_id") %>%
  left_join(in_house_source, by = "cluster_id")

# Source cluster sets must cover registry clusters; they may contain no extra rows.
for (nm in c("cassia", "in_house")) {
  source_ids <- if (nm == "cassia") unique(cassia_source$cluster_id) else unique(in_house_source$cluster_id)
  missing_ids <- setdiff(all_clusters, source_ids)
  if (length(missing_ids) > 0L) stop("Source summary missing registry clusters for ", nm, ": ", paste(missing_ids, collapse = ", "))
}

# ============================================================
# Write per-cluster judge input JSON
# - marker_panels: NULL
# - lineage: soft numeric + hits
# - reviewer_reports: structured reviewer bundle
# ============================================================
for (i in seq_len(nrow(idx))) {
  cid <- idx$cluster_id[i]
  
  step1_obj <- load_step1(cid)
  step15_obj <- load_step15_pmid_payload(cid, round_i = 1L)
  final_obj <- load_final_for_cluster(cid)
  
  reg_cassia <- registry_record(registry_df, cid, "cassia")
  reg_in_house <- registry_record(registry_df, cid, "in_house")
  reg_enrichment <- registry_record(registry_df, cid, "enrichment")
  map_cassia <- registry_mapping_bundle(reg_cassia)
  map_in_house <- registry_mapping_bundle(reg_in_house)
  map_enrichment <- registry_mapping_bundle(reg_enrichment)

  preferred_level <- if (tolower(opt$our_level) == "sub2") "sub2" else if (tolower(opt$our_level) == "sub1") "sub1" else "main"
  our_pick_source <- pick_our_top1_by_candidate_score(step1_obj, idx$our_main[i], idx$our_sub1[i], idx$our_sub2[i], preferred_level = preferred_level)
  our_top3 <- build_our_top3_by_candidate_score(step1_obj, idx$our_main[i], idx$our_sub1[i], idx$our_sub2[i])
  cassia_top3 <- build_cassia_top3(idx$cassia_source_main[i] %||% "", idx$cassia_source_sub[i] %||% "", idx$cassia_source_mixed[i] %||% "")
  cassia_label <- map_cassia$raw_label
  in_house_label <- map_in_house$raw_label
  enrichment_label <- map_enrichment$raw_label
  
  top_evidence_genes <- NULL
  if (is.list(step1_obj$candidates_for_evaluation) && is.list(step1_obj$candidates_for_evaluation$data)) {
    top_evidence_genes <- step1_obj$candidates_for_evaluation$data$Top_Evidence_Genes %||% NULL
  }
  marker_texts <- build_marker_texts_from_top_genes(top_evidence_genes)
  final_reasoning_text <- build_final_reasoning_text(final_obj)
  reasoning_text <- final_reasoning_text %||% build_reasoning_text_from_step15(step15_obj)

  std_markers <- get_top_degs(step1_obj, n = 50L)
  std_terms <- get_top_terms(step1_obj, n = 50L)
  pmids <- extract_pmids(step1_obj$evidence %||% NULL)
  step15_pmids <- extract_step15_pmids(step15_obj)
  if (length(pmids) == 0 && length(step15_pmids) > 0) pmids <- step15_pmids
  pmids_short <- if (length(pmids) > 0) paste(head(pmids, 8), collapse = ", ") else NULL

  cassia_raw <- get_cassia_raw_inputs(cid) %||% list()
  cassia_long <- cassia_raw$Scoring_Reasoning %||% cassia_raw$scoring_reasoning %||% NULL
  cassia_reasoning_short <- truncate_text(cassia_long %||%
    paste0("Cassia top1=", cassia_label %||% "NA",
           "; top markers=", std_markers %||% "NA"), 2500L)

  our_reasoning_short <- truncate_text(reasoning_text %||%
    paste0("In-house top1=", (in_house_label %||% our_pick_source$label %||% "NA"),
           "; top markers=", std_markers %||% "NA"), 2500L)

  # ---- ENRICH ----
  enrich_row <- load_enrichment_for_cluster(cid, tolower(opt$enrichment_mode), enrich_summary_df, step1_obj)
  enrich_summary <- if (!is.null(enrich_row) && nrow(enrich_row) > 0) {
    top1 <- enrich_row$predicted_label[[1]] %||% NULL
    enrich_reasoning_short <- truncate_text(
      (enrich_row$rationale[[1]] %||% NULL) %||%
        paste0("Enrich top1=", (top1 %||% "NA"), "; markers=", std_markers %||% "NA"),
      2500L
    )
    list(
      top1_cell_type = top1,
      topk_cell_types = as_list_vec(top1, k = 3L),
      cell_ontology_id = enrich_row$cl_id[[1]] %||% NULL,
      score = enrich_row$confidence[[1]] %||% NA_real_,
      confidence = enrich_row$confidence[[1]] %||% NA_real_,
      reasoning_short = enrich_reasoning_short,
      evidence_markers_std = std_markers,
      evidence_terms_std = std_terms,
      pmids_top = pmids_short,
      evidence_markers_method = enrich_row$evidence_markers[[1]] %||% NULL,
      evidence_terms_method = enrich_row$evidence_terms[[1]] %||% NULL,
      predicted_label = top1,
      cl_id = enrich_row$cl_id[[1]] %||% NULL,
      rationale = enrich_row$rationale[[1]] %||% NULL
    )
  } else {
    NULL
  }

  # Evidence fallback: some datasets have sparse bioinfo terms in step1.
  # Keep std terms when available, otherwise backfill from inter evidence terms.
  if (is.null(std_terms) || !nzchar(std_terms)) {
    std_terms <- extract_fallback_terms(enrich_summary)
  }
  if (is.null(std_terms) || !nzchar(std_terms)) {
    std_terms <- "NO_ENRICHMENT_TERMS_AVAILABLE"
  }
  if (is.null(pmids_short) || !nzchar(pmids_short)) {
    pmids_short <- "NO_PMID_EVIDENCE_AVAILABLE"
  }

  # Registry is the sole source of reviewer CL mappings. Step 08 performs integrity checks only.
  cassia_clid <- map_cassia$cell_ontology_id %||% NA_character_
  in_house_clid <- map_in_house$cell_ontology_id %||% NA_character_
  enrichment_clid <- map_enrichment$cell_ontology_id %||% NA_character_

  if (is.null(enrich_summary)) enrich_summary <- list()
  enrich_summary$top1_cell_type <- enrichment_label
  enrich_summary$predicted_label <- enrichment_label
  enrich_summary$cell_ontology_id <- enrichment_clid
  enrich_summary$cl_id <- enrichment_clid
  enrich_summary$mapping <- map_enrichment
  enrich_summary$topk_cell_types <- as.list(unique(c(
    enrichment_label,
    unlist(map_enrichment$candidate_labels, use.names = FALSE),
    unlist(map_enrichment$component_labels, use.names = FALSE)
  )))

  candidate_labels <- unique(c(
    as_char_vec(cassia_top3, k = 20L) %||% character(0),
    as_char_vec(cassia_label, k = 1L) %||% character(0),
    as_char_vec(our_top3$labels, k = 20L) %||% character(0),
    as_char_vec(in_house_label, k = 1L) %||% character(0),
    as_char_vec(enrichment_label, k = 1L) %||% character(0),
    as_char_vec(enrich_summary$topk_cell_types %||% NULL, k = 20L) %||% character(0),
    unlist(map_cassia$candidate_labels, use.names = FALSE),
    unlist(map_cassia$component_labels, use.names = FALSE),
    unlist(map_in_house$candidate_labels, use.names = FALSE),
    unlist(map_in_house$component_labels, use.names = FALSE),
    unlist(map_enrichment$candidate_labels, use.names = FALSE),
    unlist(map_enrichment$component_labels, use.names = FALSE)
  ))
  fine_labels <- unique(c(
    extract_report_labels(step15_obj),
    extract_report_labels(step1_obj)
  ))
  backoff_labels <- build_backoff_candidates(candidate_labels, cl_graph, cl_cfg, max_hops = 2L)
  allowed_greedy_candidates <- as_list_vec(unique(c(candidate_labels, fine_labels)), k = 30L)
  allowed_backoff_candidates <- as_list_vec(backoff_labels, k = 20L)

  candidate_details <- build_candidate_details(
    step1_obj,
    cassia_label,
    cassia_top3,
    in_house_label,
    our_top3$labels,
    enrichment_label,
    enrich_summary$topk_cell_types %||% NULL,
    clid_map = list(
      cassia = cassia_clid,
      our = in_house_clid,
      enrich = enrichment_clid
    )
  )

  frozen_reviewer_candidates <- c(
    build_frozen_reviewer_candidates(
      "cassia", cassia_label, cassia_top3, cassia_clid, map_cassia, step1_obj, max_k = 3L
    ),
    build_frozen_reviewer_candidates(
      "our", in_house_label, our_top3$labels, in_house_clid, map_in_house, step1_obj, max_k = 3L
    ),
    build_frozen_reviewer_candidates(
      "enrich", enrichment_label, enrich_summary$topk_cell_types %||% NULL,
      enrichment_clid, map_enrichment, step1_obj, max_k = 3L
    )
  )
  assert_frozen_candidate_top1(frozen_reviewer_candidates, "cassia", cassia_clid, cid)
  assert_frozen_candidate_top1(frozen_reviewer_candidates, "our", in_house_clid, cid)
  assert_frozen_candidate_top1(frozen_reviewer_candidates, "enrich", enrichment_clid, cid)
  
  evidence_provenance <- list(
    n_pmids = if (length(pmids) > 0) length(pmids) else 0L,
    n_enrichment_terms = count_enrichment_terms(step1_obj),
    n_degs_available = count_degs(step1_obj),
    top_evidence_genes_count = length(extract_top_evidence_genes(top_evidence_genes))
  )
  
  signals <- list(
    top_gene_prefix_counts = compute_prefix_counts(marker_texts$genes_raw)
  )
  
  packed <- list(
    cluster_id = cid,
    llm_prompt = list(system = SYSTEM_PROMPT, instruction = INSTRUCTION_PROMPT),
    inputs = list(
      cassia_summary = list(
        top1_cell_type = cassia_label,
        topk_cell_types = as_list_vec(unique(c(cassia_label, cassia_top3, unlist(map_cassia$candidate_labels, use.names = FALSE))), k = 8L),
        cell_ontology_id = cassia_clid,
        mapping = map_cassia,
        evidence_markers = std_markers,
        reasoning_short = cassia_reasoning_short,
        evidence_markers_std = std_markers,
        evidence_terms_std = std_terms,
        pmids_top = pmids_short,
        evidence_markers_method = truncate_text(cassia_raw$`Marker List` %||% NULL, 600L),
        cassia_limitations = "Automated annotation; mapping uncertainty is preserved from the frozen registry.",
        raw_inputs = list(
          predicted_main = idx$cassia_source_main[i] %||% "",
          predicted_sub = idx$cassia_source_sub[i] %||% "",
          possible_mixed = idx$cassia_source_mixed[i] %||% "",
          score = cassia_raw$Score %||% NULL
        )
      ),
      in_house_summary = list(
        top1_cell_type = in_house_label,
        topk_cell_types = as_list_vec(unique(c(in_house_label, our_top3$labels, unlist(map_in_house$candidate_labels, use.names = FALSE))), k = 8L),
        cell_ontology_id = in_house_clid,
        mapping = map_in_house,
        evidence_markers = std_markers,
        reasoning_short = our_reasoning_short,
        evidence_markers_std = std_markers,
        evidence_terms_std = std_terms,
        pmids_top = pmids_short,
        evidence_markers_method = marker_texts$markers_raw %||% NULL,
        raw_inputs = list(
          main_pred = idx$our_main[i] %||% "",
          subtype1_pred = idx$our_sub1[i] %||% "",
          subtype2_core_pred = idx$our_sub2[i] %||% "",
          candidate_score_main = our_top3$scores$main %||% our_pick_source$scores$main %||% NA_real_,
          candidate_score_sub1 = our_top3$scores$sub1 %||% our_pick_source$scores$sub1 %||% NA_real_,
          candidate_score_sub2 = our_top3$scores$sub2 %||% our_pick_source$scores$sub2 %||% NA_real_
        )
      ),
      # Mirror the canonical field because Step 09 also reads `our_summary` in several helpers.
      our_summary = list(
        top1_cell_type = in_house_label,
        topk_cell_types = as_list_vec(unique(c(in_house_label, our_top3$labels, unlist(map_in_house$candidate_labels, use.names = FALSE))), k = 8L),
        cell_ontology_id = in_house_clid,
        mapping = map_in_house,
        evidence_markers = std_markers,
        reasoning_short = our_reasoning_short,
        evidence_markers_std = std_markers,
        evidence_terms_std = std_terms,
        pmids_top = pmids_short,
        evidence_markers_method = marker_texts$markers_raw %||% NULL
      ),
      enrich_summary = enrich_summary,
      inter_summary = enrich_summary,
      candidate_details = list(entries = candidate_details),
      frozen_reviewer_candidates = frozen_reviewer_candidates,
      allowed_greedy_candidates = allowed_greedy_candidates,
      allowed_backoff_candidates = allowed_backoff_candidates,
      reviewer_reports = list(
        list(
          reviewer_role = "cassia",
          top1_cell_type = cassia_label,
          mapping = map_cassia,
          topk_cell_types = as_list_vec(cassia_top3, k = 3L),
          evidence_markers_std = std_markers,
          evidence_terms_std = std_terms,
          pmids_top = pmids_short,
          evidence_markers_method = truncate_text(cassia_raw$`Marker List` %||% NULL, 600L),
          reasoning_short = cassia_reasoning_short
        ),
        list(
          reviewer_role = "in_house",
          top1_cell_type = in_house_label,
          mapping = map_in_house,
          topk_cell_types = as_list_vec(our_top3$labels, k = 3L),
          evidence_markers_std = std_markers,
          evidence_terms_std = std_terms,
          pmids_top = pmids_short,
          evidence_markers_method = marker_texts$markers_raw %||% NULL,
          reasoning_short = our_reasoning_short
        ),
        # ---- ENRICH ----
        list(
          reviewer_role = "enrichment",
          top1_cell_type = enrichment_label,
          mapping = map_enrichment,
          topk_cell_types = enrich_summary$topk_cell_types %||% list(),
          evidence_markers_std = std_markers,
          evidence_terms_std = std_terms,
          pmids_top = pmids_short,
          evidence_markers_method = enrich_summary$evidence_markers_method %||% NULL,
          evidence_terms_method = enrich_summary$evidence_terms_method %||% NULL,
          reasoning_short = enrich_summary$reasoning_short %||% NULL
        )
      ),
      dossier = list(
        marker_genes = std_markers,
        enrichment_terms = std_terms,
        pmids_top = pmids_short
      ),
      original_evidence = list(step1_report_query = step1_obj),
      original_evidence_summary = list(
        evidence_provenance = evidence_provenance,
        signals = signals
      ),
      evidence_summary = list(
        marker_panels = NULL,
        marker_genes = std_markers,
        enrichment_terms = std_terms,
        pmids_top = pmids_short
      )
    )
  )
  
  out_file <- file.path(opt$out_dir, paste0(cid, "_LLM_JUDGE_INPUT.json"))
  jsonlite::write_json(packed, out_file, pretty = TRUE, auto_unbox = TRUE, null = "null")
}

# ============================================================
# Index CSV
# ============================================================
index_out <- file.path(opt$out_dir, "judge_input_index.csv")
readr::write_csv(
  registry_df %>%
    dplyr::select(dataset, cluster_id, reviewer, raw_label, final_status, cl_id, editor_input_type, single_cl_ready),
  index_out
)

cat("[OK] wrote judge inputs to: ", opt$out_dir, "\n", sep="")
cat("[OK] wrote index: ", index_out, "\n", sep="")
