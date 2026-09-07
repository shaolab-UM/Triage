#!/usr/bin/env Rscript
# ============================================================
# Step 12: 3-method evaluation using CL ontology similarity scoring
#
# Scoring per cluster:
# - strict: 1.0 if exact same CL term, 0.0 otherwise
# - similarity: ontology-aware percentage similarity (0-100)
#
# Inputs:
# - true_label.csv (wsnn_res.0.05, cluster_label)
# - Cassia final CSV
# - Our Step7 summary.csv (has *_cell_ontology_id)
# - Judge Step9 summary_final.csv (text only)
# - CL ontology JSON: ../SC/CL-ontology-v2025-07-30.json
#
# Outputs (under --out_dir):
# - per_cluster_detail.csv
# - method_metrics.csv
# - per_label_metrics.csv
# - confusion_<method>.csv
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
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  if (requireNamespace('ontologyIndex', quietly = TRUE)) library(ontologyIndex)
  if (requireNamespace('ontologySimilarity', quietly = TRUE)) library(ontologySimilarity)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
suppressMessages(library(Triage))
suppressMessages(library(Triage))

option_list <- list(
  make_option("--cl_json", type = "character", default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"),
              help = "CL ontology JSON [default %default]"),
  make_option("--true_label_csv", type = "character", default = "",
              help = "Ground truth CSV with wsnn_res.0.05 + cluster_label (default: final/<dataset_name>/true_label.csv)"),
  make_option("--cassia_csv", type = "character",
              default = "",
              help = "Cassia final results CSV. If empty or missing, auto-detect from CASSIA_*_*_*/01_annotation_results/"),
  make_option("--our_csv", type = "character",
              default = "",
              help = "Our Step7 post summary CSV (default: llm_outputs/<dataset_name>/post_summary/summary.csv)"),
  make_option("--judge_csv", type = "character",
              default = "",
              help = "Judge Step9 summary_final.csv (default: llm_judge_outputs/<dataset_name>/summary_final.csv)"),
  make_option("--judge_final_dir", type = "character",
              default = "",
              help = "Judge Step9 final JSON dir (default: llm_judge_outputs/<dataset_name>/final)"),
  # ---- INTER ----
  make_option("--inter_csv", type = "character", default = "",
              help = "Inter method summary CSV (default: llm_outputs/<dataset_name>/inter/post_summary/summary.csv)"),
  make_option("--out_dir", type = "character", default = "",
              help = "Output directory (default: acc_out_3methods_cl/<dataset_name>)"),
  make_option("--dataset_name", type = "character", default = "",
              help = "Dataset name for default paths (default: basename(getwd()))."),
  make_option("--exclude_clusters", type = "character", default = "cluster_99",
              help = "Comma/space-separated cluster_ids to exclude [default %default]"),
  make_option("--ols_first", type = "logical", default = TRUE,
              help = "Use OLS first for CL normalization [default TRUE]"),
  make_option("--ols_cache_dir", type = "character", default = "",
              help = "OLS cache dir (default: <out_dir>/.ols_cache)")
)

opt <- parse_args(OptionParser(option_list = option_list))
dataset_name <- if (nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- get_dataset_config(dataset_name, project_root)
if (!nzchar(opt$true_label_csv)) opt$true_label_csv <- file.path(project_root, "final", dataset_name, "true_label.csv")
if (!nzchar(opt$our_csv)) opt$our_csv <- file.path(cfg$llm_outputs_root, "post_summary", "summary.csv")
if (!nzchar(opt$judge_csv)) opt$judge_csv <- file.path(cfg$judge_outputs_root, "summary_final.csv")
if (!nzchar(opt$judge_final_dir)) opt$judge_final_dir <- file.path(cfg$judge_outputs_root, "final")
if (!nzchar(opt$inter_csv)) opt$inter_csv <- file.path(cfg$llm_outputs_root, "inter", "post_summary", "summary.csv")
if (!nzchar(opt$out_dir)) opt$out_dir <- file.path(project_root, "final", dataset_name, "acc_out_3methods_cl")
ensure_dir(opt$out_dir)

ols_cache_dir <- opt$ols_cache_dir
if (!nzchar(ols_cache_dir)) ols_cache_dir <- file.path(opt$out_dir, ".ols_cache")
cl_cfg <- make_cl_cfg(opt$cl_json, prefer_ols = isTRUE(opt$ols_first), cache_dir = ols_cache_dir)

auto_pick_cassia_csv <- function(dataset_name) {
  candidates <- Sys.glob(file.path(project_root, "final", dataset_name, "CASSIA_*_*_*", "01_annotation_results", "annotation_cassia_FINAL_RESULTS.csv"))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return(NA_character_)
  mt <- suppressWarnings(as.numeric(file.info(candidates)$mtime))
  candidates[order(mt, decreasing = TRUE)[1]]
}

if (!nzchar(opt$cassia_csv) || !file.exists(opt$cassia_csv)) {
  picked <- auto_pick_cassia_csv(dataset_name)
  if (!is.na(picked)) {
    opt$cassia_csv <- picked
    message("[INFO] auto-picked cassia_csv: ", opt$cassia_csv)
  }
}

if (!file.exists(opt$cl_json)) stop("cl_json not found: ", opt$cl_json)
if (!file.exists(opt$true_label_csv)) stop("true_label_csv not found: ", opt$true_label_csv)
if (!file.exists(opt$cassia_csv)) stop("cassia_csv not found: ", opt$cassia_csv)
if (!file.exists(opt$our_csv)) stop("our_csv not found: ", opt$our_csv)
if (!file.exists(opt$judge_csv) && !dir.exists(opt$judge_final_dir)) {
  stop("Need either judge_csv or judge_final_dir. Missing: ", opt$judge_csv, " and ", opt$judge_final_dir)
}

# ---------------------------
# CL helpers
# ---------------------------
normalize_txt <- function(x) {
  x <- as.character(x)[1] %||% ""
  if (is.na(x)) x <- ""
  x <- stringr::str_to_lower(stringr::str_trim(x))
  x <- stringr::str_replace_all(x, "[_\\-]+", " ")
  x <- stringr::str_replace_all(x, "[,;/]+", " ")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

normalize_cluster_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  ifelse(
    is.na(x) | !nzchar(x),
    NA_character_,
    ifelse(stringr::str_detect(x, "^cluster_"), x, paste0("cluster_", x))
  )
}

split_candidates <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- stringr::str_trim(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return(character(0))
  # NOTE: Do not split on "and"; it breaks valid labels like
  # "Hematopoietic Stem and Progenitor Cell (HSPC)".
  parts <- unlist(stringr::str_split(x, "\\s*/\\s*|\\s*,\\s*|\\s*;\\s*|\\s*\\|\\s*|\\s*\\+\\s*"))
  parts <- stringr::str_trim(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

split_subtypes <- function(x) {
  split_candidates(x)
}

# Cassia reports multiple fields (main/sub/mixed). For evaluation and reporting, we must keep
# Cassia's *reported/top1* label (no cherry-picking a different substring just because it maps better).
cassia_is_bad_label <- function(x) {
  x2 <- as.character(x)[1] %||% ""
  if (is.na(x2)) x2 <- ""
  x2 <- stringr::str_to_lower(stringr::str_trim(x2))
  if (!nzchar(x2)) return(TRUE)
  stringr::str_detect(x2, "mixed|doublet|triplet|multiplet|artifact|unresolved|technical|contamination|sequencing")
}

cassia_pick_top1_label <- function(main, sub, mix) {
  m1 <- split_candidates(main)
  s1 <- split_candidates(sub)
  x1 <- split_candidates(mix)
  m_first <- if (length(m1) > 0) m1[[1]] else NA_character_
  s_first <- if (length(s1) > 0) s1[[1]] else NA_character_
  x_first <- if (length(x1) > 0) x1[[1]] else NA_character_

  if (!is.na(m_first) && nzchar(m_first) && !cassia_is_bad_label(m_first)) return(m_first)
  if (!is.na(s_first) && nzchar(s_first) && !cassia_is_bad_label(s_first)) return(s_first)
  if (!is.na(x_first) && nzchar(x_first) && !cassia_is_bad_label(x_first)) return(x_first)
  NA_character_
}

strip_parenthetical <- function(x) {
  x <- as.character(x)[1] %||% ""
  x <- stringr::str_trim(x)
  x <- stringr::str_replace(x, "\\s*\\(.*\\)\\s*$", "")
  x <- stringr::str_replace(x, "\\s*\\[.*\\]\\s*$", "")
  stringr::str_trim(x)
}

cl <- jsonlite::fromJSON(opt$cl_json, simplifyVector = FALSE)

# ---- OntologySimilarity setup (percentage similarity 0-100) ----
# We infer direct parents from the precomputed 'ancestors' map in the CL JSON (distance==1).
cl_ont <- NULL
if (requireNamespace('ontologyIndex', quietly = TRUE)) {
  ids_all <- names(cl)
  parents_all <- lapply(ids_all, function(id) {
    anc <- cl[[id]]$ancestors %||% NULL
    if (is.null(anc)) return(character(0))
    d <- unlist(anc, use.names = TRUE)
    if (length(d) == 0) return(character(0))
    names(d[d == 1])
  })
  names(parents_all) <- ids_all
  nm_all <- ids_all
  exports <- getNamespaceExports('ontologyIndex')
  if ('ontology_index' %in% exports) {
    cl_ont <- ontologyIndex::ontology_index(id = ids_all, name = nm_all, parents = parents_all)
  } else if ('ontologyIndex' %in% exports) {
    cl_ont <- ontologyIndex::ontologyIndex(id = ids_all, name = nm_all, parents = parents_all)
  }
}

sim_cache <- new.env(parent = emptyenv())
get_sim_pair_cached <- function(a, b) {
  # returns NA if cannot compute
  key <- paste(a, b, sep = '||')
  if (exists(key, envir = sim_cache, inherits = FALSE)) return(get(key, envir = sim_cache))
  if (is.null(cl_ont) || !requireNamespace('ontologySimilarity', quietly = TRUE)) {
    assign(key, NA_real_, envir = sim_cache)
    assign(paste(b, a, sep = '||'), NA_real_, envir = sim_cache)
    return(NA_real_)
  }
  # Correct usage: ontology first, term_sets as list of character vectors
  g <- tryCatch(
    ontologySimilarity::get_sim_grid(ontology = cl_ont, term_sets = list(c(a), c(b))),
    error = function(e) NULL
  )
  if (is.null(g)) {
    assign(key, NA_real_, envir = sim_cache)
    assign(paste(b, a, sep = '||'), NA_real_, envir = sim_cache)
    return(NA_real_)
  }
  # Result is a 2x2 matrix; similarity is at [1,2] or [2,1]
  s <- suppressWarnings(as.numeric(g[1, 2]))
  if (length(s) == 0 || is.na(s)) s <- NA_real_
  # normalize to 0-100 if the package returns 0-1
  if (!is.na(s) && s <= 1) s <- s * 100
  assign(key, s, envir = sim_cache)
  assign(paste(b, a, sep = '||'), s, envir = sim_cache)
  s
}

label_to_clid <- new.env(parent = emptyenv())
add_label <- function(label, clid) {
  if (is.null(label)) return()
  label <- as.character(label)[1]
  if (is.na(label) || !nzchar(label)) return()
  n <- normalize_txt(label)
  if (!nzchar(n)) return()
  if (!exists(n, envir = label_to_clid, inherits = FALSE)) {
    assign(n, clid, envir = label_to_clid)
  }
  # plural variant
  if (stringr::str_ends(n, "s") && !exists(stringr::str_replace(n, "s$", ""), envir = label_to_clid, inherits = FALSE)) {
    assign(stringr::str_replace(n, "s$", ""), clid, envir = label_to_clid)
  }
}

for (id in names(cl)) {
  if (!startsWith(id, "CL:")) next
  term <- cl[[id]]
  add_label(term$label %||% NULL, id)
  syn <- term$synonyms %||% NULL
  if (!is.null(syn)) {
    syns <- if (is.list(syn)) unlist(syn, recursive = TRUE, use.names = FALSE) else syn
    for (s in syns) add_label(s, id)
  }
}

map_to_clid <- function(x) {
  n <- normalize_txt(x)
  if (is.na(n) || !nzchar(n)) return(NA_character_)
  if (exists(n, envir = label_to_clid, inherits = FALSE)) return(get(n, envir = label_to_clid, inherits = FALSE))
  # try strip parenthetical
  n2 <- normalize_txt(strip_parenthetical(x))
  if (nzchar(n2) && exists(n2, envir = label_to_clid, inherits = FALSE)) return(get(n2, envir = label_to_clid, inherits = FALSE))
  # try plural->singular
  n3 <- stringr::str_replace(n, "s$", "")
  if (nzchar(n3) && exists(n3, envir = label_to_clid, inherits = FALSE)) return(get(n3, envir = label_to_clid, inherits = FALSE))

  # Heuristic fallbacks for common pipeline labels (Cassia/Our/Judge)
  # Use canonical CL IDs present in ../SC/CL-ontology-v2025-07-30.json
  if (stringr::str_detect(n, "oligodendrocyte precursor|\\bopc\\b")) return("CL:0002453")
  if (stringr::str_detect(n, "oligodendrocyte")) return("CL:0000128")
  if (stringr::str_detect(n, "astrocyte")) return("CL:0000127")
  if (stringr::str_detect(n, "microglia")) return("CL:0000129")
  if (stringr::str_detect(n, "endothelial")) return("CL:0000115")
  if (stringr::str_detect(n, "pericyte")) return("CL:0000669")
  if (stringr::str_detect(n, "fibroblast")) return("CL:0000057")
  if (stringr::str_detect(n, "vascular leptomeningeal")) return("CL:4023051")
  if (stringr::str_detect(n, "vlmc")) return("CL:4023051")
  if (stringr::str_detect(n, "perivascular stromal")) return("CL:0000669")
  if (stringr::str_detect(n, "vascular")) return("CL:0000115")
  if (stringr::str_detect(n, "gaba|gabaergic|interneuron")) return("CL:0000617")
  if (stringr::str_detect(n, "glutamate|glutamatergic|excitatory")) return("CL:0000679")
  if (stringr::str_detect(n, "myeloid")) return("CL:0000763")
  if (stringr::str_detect(n, "\\bneuron\\b")) return("CL:0000540")

  NA_character_
}

map_label_to_clid <- function(label) {
  res <- normalize_cl_three_state(label, "", cl_cfg)
  res$final_clid %||% NA_character_
}

map_labels_to_clids <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0) return(character(0))
  clids <- vapply(labels, map_label_to_clid, character(1))
  clids <- clids[!is.na(clids) & nzchar(clids)]
  unique(clids)
}

is_valid_clid <- function(x) {
  !is.na(x) && nzchar(x) && startsWith(x, "CL:") && !is.null(cl[[x]])
}




# Manual equivalence mapping for evaluation-only canonicalization.
# Use this to merge biologically near-equivalent labels into one eval target.
EVAL_CANONICAL_CLID <- c(
  "CL:0000826" = "CL:0000817"  # pro-B cell -> precursor B cell
)

canonicalize_eval_clid <- function(x) {
  x <- as.character(x)[1] %||% ""
  x <- stringr::str_trim(x)
  if (!nzchar(x)) return(NA_character_)
  if (x %in% names(EVAL_CANONICAL_CLID)) return(as.character(EVAL_CANONICAL_CLID[[x]]))
  x
}

# ontology-aware percentage similarity (0-100), computed via ontologySimilarity::get_sim_grid
score_similarity <- function(pred, truth) {
  pred <- canonicalize_eval_clid(pred)
  truth <- canonicalize_eval_clid(truth)
  if (!is_valid_clid(pred) || !is_valid_clid(truth)) return(0)
  if (pred == truth) return(100)
  s <- get_sim_pair_cached(pred, truth)
  if (!is.na(s)) return(s)
  0
}

dist_between <- function(a, b) {
  if (!is_valid_clid(a) || !is_valid_clid(b)) return(Inf)
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

# ---------------------------
# GT loading + GT CLIDs (direct label -> CLID)
# ---------------------------
true_df <- readr::read_csv(opt$true_label_csv, show_col_types = FALSE)
find_col <- function(df, choices) {
  for (c in choices) if (c %in% names(df)) return(c)
  NA_character_
}
gt_label_col <- find_col(true_df, c("cluster_label", "gt_label", "truth", "True Cell Type"))
if (is.na(gt_label_col)) {
  stop("true_label_csv missing label column. Expected one of: cluster_label, gt_label, truth, True Cell Type.")
}

gt <- true_df %>%
  transmute(
    cluster_id = if ("cluster_id" %in% names(true_df)) normalize_cluster_id(.data[["cluster_id"]]) else paste0("cluster_", as.character(.data[["wsnn_res.0.05"]])),
    gt_label = as.character(.data[[gt_label_col]])
  )

gt <- gt %>%
  mutate(gt_clid = vapply(gt_label, function(lbl) {
    res <- normalize_cl_three_state(lbl, "", cl_cfg)
    res$final_clid %||% NA_character_
  }, character(1)))

exclude_ids <- opt$exclude_clusters %||% ""
exclude_ids <- as.character(exclude_ids)
exclude_ids <- stringr::str_split(exclude_ids, "[,;\\s]+", simplify = TRUE)
exclude_ids <- as.character(exclude_ids)
exclude_ids <- exclude_ids[nzchar(exclude_ids)]
if (length(exclude_ids) > 0) {
  gt <- gt %>% filter(!(.data$cluster_id %in% exclude_ids))
}

if (any(is.na(gt$gt_clid) | !nzchar(gt$gt_clid))) {
  stop("Failed to map some GT labels to CL IDs: ", paste(gt$gt_label[is.na(gt$gt_clid) | !nzchar(gt$gt_clid)], collapse = ", "))
}

# ---------------------------
# Load predictions
# ---------------------------
cassia_raw <- readr::read_csv(opt$cassia_csv, show_col_types = FALSE)
cassia_id_col <- find_col(cassia_raw, c("True Cell Type", "cluster_id", "cell_id", "ID"))
cassia_main_col <- find_col(cassia_raw, c("Predicted Main Cell Type", "pred_main", "cassia_top1_label"))
cassia_sub_col <- find_col(cassia_raw, c("Predicted Sub Cell Types", "pred_subtypes", "cassia_subtypes"))
cassia_mix_col <- find_col(cassia_raw, c("Possible Mixed Cell Types", "possible_mixed", "cassia_mixed"))

if (is.na(cassia_id_col) || is.na(cassia_main_col) || is.na(cassia_sub_col) || is.na(cassia_mix_col)) {
  stop("Cassia csv missing columns. Expected one of: ",
       "id={True Cell Type, cluster_id, cell_id, ID}; ",
       "pred_main={Predicted Main Cell Type, pred_main, cassia_top1_label}; ",
       "pred_subtypes={Predicted Sub Cell Types, pred_subtypes, cassia_subtypes}; ",
       "mixed={Possible Mixed Cell Types, possible_mixed, cassia_mixed}.")
}

our_raw <- readr::read_csv(opt$our_csv, show_col_types = FALSE)
req_our <- c("cluster_id", "main_pred", "main_cell_ontology_id")
miss_o <- setdiff(req_our, names(our_raw))
if (length(miss_o) > 0) stop("Our csv missing columns: ", paste(miss_o, collapse = ", "))

load_judge_from_final_json <- function(final_dir) {
  fs <- list.files(final_dir, pattern = "_LLM_JUDGE_FINAL\\.json$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  rows <- lapply(fs, function(fp) {
    x <- tryCatch(jsonlite::fromJSON(fp, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(x) || !is.list(x)) return(NULL)
    fd <- x$final_decision %||% list()
    mv <- x$method_verdict %||% list()
    cass <- mv$cassia %||% list()
    our <- mv$our_method %||% list()
    # Accept both supported final-label field names: `final_cell_type` and `primary_cell_type`.
    final_ct <- if (!is.null(fd$primary_cell_type) && nzchar(fd$primary_cell_type)) {
      fd$primary_cell_type
    } else if (!is.null(fd$final_cell_type) && nzchar(fd$final_cell_type)) {
      fd$final_cell_type
    } else {
      ""
    }
    tibble::tibble(
      cluster_id = as.character(x$cluster_id %||% NA_character_),
      final_cell_type = final_ct,
      final_cell_ontology_id = as.character(fd$final_cell_ontology_id %||% ""),
      decision_category = as.character(fd$decision_category %||% ""),
      cassia_predicted_cell_type = as.character(cass$predicted_cell_type %||% ""),
      cassia_cell_ontology_id = as.character(cass$cell_ontology_id %||% ""),
      our_predicted_cell_type = as.character(our$predicted_cell_type %||% ""),
      our_cell_ontology_id = as.character(our$cell_ontology_id %||% "")
    )
  })
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(NULL)
  out
}

judge_raw <- NULL
if (dir.exists(opt$judge_final_dir)) {
  judge_raw <- load_judge_from_final_json(opt$judge_final_dir)
}
if (is.null(judge_raw)) {
  judge_raw <- readr::read_csv(opt$judge_csv, show_col_types = FALSE)
}
judge_raw$cluster_id <- normalize_cluster_id(judge_raw$cluster_id)

# ---- INTER ----
inter_raw <- NULL
if (file.exists(opt$inter_csv)) {
  inter_raw <- readr::read_csv(opt$inter_csv, show_col_types = FALSE)
}
if (is.null(inter_raw)) {
  inter_raw <- tibble::tibble(cluster_id = character(), predicted_label = character(), cl_id = character())
}
inter_raw$cluster_id <- normalize_cluster_id(inter_raw$cluster_id)

if (!"cluster_id" %in% names(judge_raw)) stop("Judge data missing column: cluster_id")
if (!"final_cell_type" %in% names(judge_raw) && !"primary_cell_type" %in% names(judge_raw)) {
  stop("Judge data missing column: final_cell_type or primary_cell_type")
}
if (!"decision_category" %in% names(judge_raw)) judge_raw$decision_category <- ""
if (!"cassia_predicted_cell_type" %in% names(judge_raw)) judge_raw$cassia_predicted_cell_type <- ""
if (!"our_predicted_cell_type" %in% names(judge_raw)) judge_raw$our_predicted_cell_type <- ""
if (!"cassia_cell_ontology_id" %in% names(judge_raw)) judge_raw$cassia_cell_ontology_id <- ""
if (!"our_cell_ontology_id" %in% names(judge_raw)) judge_raw$our_cell_ontology_id <- ""
if (!"final_cell_ontology_id" %in% names(judge_raw)) judge_raw$final_cell_ontology_id <- ""

# ---------------------------
# Per-cluster scoring helpers
# ---------------------------
best_from_candidates <- function(cand_labels, cand_clids, gt_clid) {
  if (length(cand_clids) == 0) {
    return(list(best_score = 0, best_label = NA_character_, best_clid = NA_character_))
  }
  scores <- vapply(cand_clids, score_similarity, numeric(1), truth = gt_clid)
  i <- which.max(scores)
  list(best_score = scores[i], best_label = cand_labels[i] %||% NA_character_, best_clid = cand_clids[i] %||% NA_character_)
}

# ---------------------------
# Assemble per-cluster detail
# ---------------------------
detail <- gt %>%
  mutate(cluster_id_int = suppressWarnings(as.integer(cluster_id))) %>%
  arrange(cluster_id_int) %>%
  select(-cluster_id_int) %>%
  # Cassia
  left_join(
    cassia_raw %>%
      transmute(
        cluster_id = normalize_cluster_id(.data[[cassia_id_col]]),
        cassia_main = as.character(.data[[cassia_main_col]] %||% ""),
        cassia_sub  = as.character(.data[[cassia_sub_col]] %||% ""),
        cassia_mix  = as.character(.data[[cassia_mix_col]] %||% "")
      ),
    by = "cluster_id"
  ) %>%
  # Our
  left_join(
    our_raw %>%
      transmute(
        cluster_id = normalize_cluster_id(.data[["cluster_id"]]),
        our_main = as.character(.data[["main_pred"]] %||% ""),
        our_main_clid = as.character(.data[["main_cell_ontology_id"]] %||% ""),
        our_main_conf = suppressWarnings(as.numeric(.data[["main_confidence"]] %||% NA_real_)),
        our_sub1 = as.character(.data[["subtype1_pred"]] %||% ""),
        our_sub1_clid = as.character(.data[["subtype1_cell_ontology_id"]] %||% ""),
        our_sub1_conf = suppressWarnings(as.numeric(.data[["subtype1_confidence"]] %||% NA_real_)),
        our_sub2 = as.character(.data[["subtype2_core_pred"]] %||% ""),
        our_sub2_clid = as.character(.data[["subtype2_core_cell_ontology_id"]] %||% ""),
        our_sub2_conf = suppressWarnings(as.numeric(.data[["subtype2_confidence"]] %||% NA_real_))
      ),
    by = "cluster_id"
  ) %>%
  # Judge
  left_join(
    judge_raw %>%
      transmute(
        cluster_id = normalize_cluster_id(.data[["cluster_id"]]),
        judge_final = as.character(.data[["final_cell_type"]] %||% ""),
        judge_final_clid_raw = as.character(.data[["final_cell_ontology_id"]] %||% ""),
        judge_decision_category = as.character(.data[["decision_category"]] %||% ""),
        judge_cassia_label = as.character(.data[["cassia_predicted_cell_type"]] %||% ""),
        judge_cassia_clid_raw = as.character(.data[["cassia_cell_ontology_id"]] %||% ""),
        judge_our_label = as.character(.data[["our_predicted_cell_type"]] %||% ""),
        judge_our_clid_raw = as.character(.data[["our_cell_ontology_id"]] %||% "")
      ),
    by = "cluster_id"
  ) %>%
  # ---- INTER ----
  left_join(
    inter_raw %>%
      transmute(
        cluster_id = normalize_cluster_id(.data[["cluster_id"]]),
        inter_label = as.character(.data[["predicted_label"]] %||% ""),
        inter_clid_raw = as.character(.data[["cl_id"]] %||% "")
      ),
    by = "cluster_id"
  ) %>%
  rowwise() %>%
  mutate(
    # Cassia reported/top1 label (main > sub > mixed; avoid technical tokens)
    cassia_top1_label = cassia_pick_top1_label(cassia_main, cassia_sub, cassia_mix),
    cassia_top1_clid = map_label_to_clid(cassia_top1_label),
    cassia_top1_score = score_similarity(cassia_top1_clid, gt_clid),
    cassia_sub_clids = list(map_labels_to_clids(split_subtypes(cassia_sub))),
    score_strict = if (cassia_top1_clid == gt_clid && is_valid_clid(cassia_top1_clid)) 1 else 0,

    # Cassia candidates
    cassia_candidates = list(unique(c(
      split_candidates(cassia_main),
      split_candidates(cassia_sub),
      split_candidates(cassia_mix)
    ))),
    cassia_candidate_clids = list(vapply(cassia_candidates, map_label_to_clid, character(1))),
    cassia_best = list(best_from_candidates(cassia_candidates, cassia_candidate_clids, gt_clid)),
    cassia_best_score = cassia_best$best_score,
    cassia_best_label = cassia_best$best_label,
    cassia_best_clid  = cassia_best$best_clid,

    # Our candidates (prefer ontology ids; fallback to label mapping)
    our_candidate_labels = list(unique(c(our_main, our_sub1, our_sub2))),
    our_candidate_clids_raw = list(c(our_main_clid, our_sub1_clid, our_sub2_clid)),
    our_candidate_clids = list({
      clids <- unlist(our_candidate_clids_raw)
      clids <- ifelse(nzchar(clids), clids, NA_character_)
      clids <- ifelse(vapply(clids, is_valid_clid, logical(1)), clids, NA_character_)
      # fill missing using label mapping
      labs <- unlist(our_candidate_labels)
      for (i in seq_along(clids)) {
        if (is.na(clids[i]) && nzchar(labs[i])) clids[i] <- map_label_to_clid(labs[i])
      }
      clids
    }),
    our_best = list(best_from_candidates(our_candidate_labels, our_candidate_clids, gt_clid)),
    our_best_score = our_best$best_score,
    our_best_label = our_best$best_label,
    our_best_clid  = our_best$best_clid,

    # Our top1 (deterministic main prediction; no best-of-candidates)
    our_top1_label = {
      x <- as.character(our_main)[1] %||% ""
      x <- stringr::str_trim(x)
      if (!nzchar(x)) NA_character_ else x
    },
    our_top1_clid = {
      x <- as.character(our_main_clid)[1] %||% ""
      x <- stringr::str_trim(x)
      if (nzchar(x) && is_valid_clid(x)) x else map_label_to_clid(our_top1_label)
    },
    our_top1_score = score_similarity(our_top1_clid, gt_clid),

    # Judge candidates (text -> clid)
    judge_candidate_labels = list(split_candidates(judge_final)),
    judge_candidate_clids = list(vapply(judge_candidate_labels, map_label_to_clid, character(1))),
    judge_best = list(best_from_candidates(judge_candidate_labels, judge_candidate_clids, gt_clid)),
    judge_best_score = judge_best$best_score,
    judge_best_label = judge_best$best_label,
    judge_best_clid  = judge_best$best_clid,

    # Judge top1 (final label only; prefer raw CLID if provided)
    judge_top1_label = {
      x <- as.character(judge_final)[1] %||% ""
      x <- stringr::str_trim(x)
      if (!nzchar(x)) NA_character_ else x
    },
    judge_top1_clid = {
      raw <- as.character(judge_final_clid_raw)[1] %||% ""
      raw <- stringr::str_trim(raw)
       if (nzchar(raw) && is_valid_clid(raw)) raw else map_label_to_clid(judge_top1_label)
    },
    judge_top1_score = score_similarity(judge_top1_clid, gt_clid),

    # Judge as selector (preferred): use decision_category to choose Cassia vs Our prediction;
    # only allow third label when both_incorrect.
    judge_choice_label = {
      dc <- normalize_txt(judge_decision_category)
      if (dc == "cassia_better") judge_cassia_label
      else if (dc == "our_better") judge_our_label
      else judge_final
    },
    judge_choice_clid = {
      dc <- normalize_txt(judge_decision_category)
      clid <- NA_character_
      pick_raw <- function(x) {
        x2 <- as.character(x)[1] %||% ""
        x2 <- stringr::str_trim(x2)
        if (!nzchar(x2)) return(NA_character_)
        if (is_valid_clid(x2)) return(x2)
        NA_character_
      }
      if (dc == "cassia_better") clid <- pick_raw(judge_cassia_clid_raw)
      if (dc == "our_better") clid <- pick_raw(judge_our_clid_raw)
      if (is.na(clid)) clid <- map_label_to_clid(judge_choice_label)
      clid
    },
    judge_choice_score = score_similarity(judge_choice_clid, gt_clid)
    ,
    # ---- INTER ----
    inter_clid = {
      raw <- as.character(inter_clid_raw %||% "")
      raw <- stringr::str_trim(raw)
      if (nzchar(raw) && is_valid_clid(raw)) raw else map_label_to_clid(inter_label)
    },
    inter_score = score_similarity(inter_clid, gt_clid),
    inter_sub_clids = list(character(0))
  ) %>%
  ungroup() %>%
  select(
    cluster_id, gt_label, gt_clid,
    cassia_top1_score, cassia_top1_label, cassia_top1_clid,
    cassia_sub,
    score_strict,
    cassia_best_score, cassia_best_label, cassia_best_clid,
    our_best_score, our_best_label, our_best_clid,
    judge_best_score, judge_best_label, judge_best_clid,
    our_top1_score, our_top1_label, our_top1_clid,
    judge_top1_score, judge_top1_label, judge_top1_clid,
    judge_choice_score, judge_choice_label, judge_choice_clid,
    inter_score, inter_label, inter_clid,
    judge_decision_category
  )

readr::write_csv(detail, file.path(opt$out_dir, "per_cluster_detail.csv"))

simple_table <- detail %>%
  transmute(
    cluster_id,
    gt_label,
    gt_clid,
    cassia_label = cassia_top1_label,
    cassia_clid = cassia_top1_clid,
    cassia_score = cassia_top1_score,
    our_label = our_top1_label,
    our_clid = our_top1_clid,
    our_score = our_top1_score,
    judge_label = judge_choice_label,
    judge_clid = judge_choice_clid,
    judge_score = judge_choice_score,
    inter_label = inter_label,
    inter_clid = inter_clid,
    inter_score = inter_score
  )
readr::write_csv(simple_table, file.path(opt$out_dir, "simple_table.csv"))

results_eval_with_2scores <- detail %>%
  transmute(
    cluster_id,
    gt_label,
    gt_clid,
    pred_main_label = cassia_top1_label,
    pred_main_clid = cassia_top1_clid,
    pred_subtypes_raw = cassia_sub,
    score_strict,
    score_similarity = cassia_top1_score
  )
readr::write_csv(results_eval_with_2scores, file.path(opt$out_dir, "results_eval_with_2scores.csv"))

# ---------------------------
# Metrics (CLID-only)
# ---------------------------
metrics_one <- function(df, best_score_col, method_name) {
  s <- df[[best_score_col]]
  s[is.na(s)] <- 0
  cl_weighted_top1 <- mean(s)
  cl_strict_top1 <- mean(s == 100)
  list(
    method = tibble::tibble(
      method = method_name,
      cl_weighted_top1 = cl_weighted_top1,
      cl_strict_top1 = cl_strict_top1
    )
  )
}

metrics_2scores <- function(df, pred_main_clid_col, method_name) {
  pred_main <- df[[pred_main_clid_col]]
  gt <- df[["gt_clid"]]
  score_strict <- mapply(function(p, g) {
    if (!is_valid_clid(p) || !is_valid_clid(g)) return(0)
    if (p == g) return(1) else return(0)
  }, pred_main, gt)
  score_sim <- mapply(score_similarity, pred_main, gt)
  tibble::tibble(
    method = method_name,
    n = length(score_sim),
    acc_strict = mean(score_strict, na.rm = TRUE),
    cl_avg_similarity = mean(score_sim, na.rm = TRUE)
  )
}

# Main table: deterministic top-1 for all methods (no oracle, no best-of)
m_cassia <- metrics_one(detail, "cassia_top1_score", "cassia")
m_our    <- metrics_one(detail, "our_top1_score", "our")
m_judge  <- metrics_one(detail, "judge_top1_score", "judge")
m_enrich <- metrics_one(detail, "inter_score", "enrich")

# Best-of-candidate metrics (reference only; Methods declares best-of evaluation)
m_cassia_t1 <- metrics_one(detail, "cassia_best_score", "cassia_best")
m_our_t1    <- metrics_one(detail, "our_best_score", "our_best")
m_judge_t1  <- metrics_one(detail, "judge_best_score", "judge_best")
m_enrich_t1 <- metrics_one(detail, "inter_score", "enrich")

method_metrics <- bind_rows(m_cassia$method, m_our$method, m_judge$method, m_enrich$method)
readr::write_csv(method_metrics, file.path(opt$out_dir, "method_metrics.csv"))

method_metrics_top1 <- bind_rows(m_cassia_t1$method, m_our_t1$method, m_judge_t1$method, m_enrich_t1$method)
readr::write_csv(method_metrics_top1, file.path(opt$out_dir, "method_metrics_top1.csv"))

metrics_2 <- bind_rows(
  metrics_2scores(detail, "cassia_top1_clid", "CASSIA"),
  metrics_2scores(detail, "our_top1_clid", "OUR"),
  metrics_2scores(detail, "judge_top1_clid", "JUDGE"),
  metrics_2scores(detail, "inter_clid", "INTER")
)
readr::write_csv(metrics_2, file.path(opt$out_dir, "method_metrics_2scores.csv"))

cat("\n[OK] wrote outputs under: ", opt$out_dir, "\n", sep = "")
print(method_metrics)

cat("\n[INFO] top1-only metrics:\n")
print(method_metrics_top1)

cat("\n[INFO] 2-score metrics (strict + similarity):\n")
print(metrics_2)
