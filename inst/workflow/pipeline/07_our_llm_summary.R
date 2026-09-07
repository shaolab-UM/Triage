#!/usr/bin/env Rscript
# ============================================================
# Step 07 (post)
# Build llm_outputs/<dataset>/post_summary/summary.csv
# from Step6 outputs ("post" artifacts):
#   - step1.5_citation_fix_outputs/*_step1.5_citation_fix_output.json
#   - step2_validation_outputs/*_step2_validation_output.json (optional)
#   - step1.5_citation_fix_outputs/*_GATE.json (optional)
#   - run_summary.csv (optional)
#
# Note:
# Some clusters may have schema-only Step1.5 outputs (main_type/subtype null).
# In that case we fall back to Step1 report outputs (and then final_passed) to
# populate predictions, while still attaching post (gate/step2/run_summary).
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
  library(purrr)
  library(readr)
  library(stringr)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
suppressMessages(library(Triage))

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

option_list <- list(
  make_option("--out_root", type = "character",
              default = "",
              help = "Dataset output root produced by Step6 (contains step1.5_citation_fix_outputs/, step2_validation_outputs/, etc.). Default: llm_outputs/<dataset_name>"),
  make_option("--dataset_name", type = "character", default = "",
              help = "Dataset name for default paths (default: basename(getwd()))."),
  make_option("--out_dir", type = "character",
              default = "",
              help = "Output dir for post_summary/summary.csv (default: <out_root>/post_summary)"),
  make_option("--round", type = "integer",
              default = 1,
              help = "Round number to summarize (default 1)"),
  make_option("--include_failed", action = "store_true", default = TRUE,
              help = "Include clusters even if Step2 failed/missing (default TRUE)")
)
opt <- parse_args(OptionParser(option_list = option_list))
dataset_name <- if (nzchar(opt$dataset_name)) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- Triage:::get_dataset_config(dataset_name, project_root)

root <- opt$out_root
if (!nzchar(root)) root <- cfg$llm_outputs_root
if (!dir.exists(root)) stop("out_root not found: ", root)
out_dir <- opt$out_dir
if (!nzchar(out_dir)) out_dir <- file.path(root, "post_summary")
ensure_dir(out_dir)

dir_step15 <- file.path(root, "step1.5_citation_fix_outputs")
dir_step1  <- file.path(root, "step1_report_outputs")
dir_step2  <- file.path(root, "step2_validation_outputs")
dir_final  <- file.path(root, "final_passed")

if (!dir.exists(dir_step15) && !dir.exists(dir_final)) {
  stop(
    "Neither step1.5_citation_fix_outputs nor final_passed exist under out_root: ",
    root,
    "\nExpected at least one of: ",
    dir_step15,
    " or ",
    dir_final
  )
}

round_i <- as.integer(opt$round)

# Prefer final_passed as canonical output; fall back to Step1.5/Step1 if missing.
pattern_step15 <- sprintf("_round%d_step1\\.5_citation_fix_output\\.json$", round_i)
pattern_final  <- "_FINAL_passed\\.json$"

files_step15 <- if (dir.exists(dir_step15)) list.files(dir_step15, pattern = pattern_step15, full.names = TRUE) else character(0)
files_final  <- if (dir.exists(dir_final))  list.files(dir_final,  pattern = pattern_final,  full.names = TRUE) else character(0)

cluster_id_from_step15 <- function(p) sub(pattern_step15, "", basename(p))
cluster_id_from_final  <- function(p) sub(pattern_final,  "", basename(p))

clusters <- sort(unique(c(vapply(files_step15, cluster_id_from_step15, character(1)), vapply(files_final, cluster_id_from_final, character(1)))))
if (length(clusters) == 0) stop("No clusters found under step1.5_citation_fix_outputs or final_passed.")

get1 <- function(x, ...) {
  # safe nested getter; returns NA if missing
  keys <- list(...)
  cur <- x
  for (k in keys) {
    if (is.null(cur) || !is.list(cur) || is.null(cur[[k]])) return(NA)
    cur <- cur[[k]]
  }
  if (is.null(cur)) return(NA)
  if (is.atomic(cur) && length(cur) >= 1) return(cur[[1]])
  # if list/object, stringify
  tryCatch(jsonlite::toJSON(cur, auto_unbox = TRUE, null = "null"), error = function(e) NA)
}

get_conf_score <- function(sec) {
  if (is.null(sec) || !is.list(sec)) return(NA_real_)
  direct <- suppressWarnings(as.numeric(get1(sec, "confidence")))
  if (!is.na(direct)) return(direct)
  bd <- sec$confidence_score_breakdown %||% NULL
  if (is.null(bd) || !is.list(bd)) return(NA_real_)
  vals <- suppressWarnings(as.numeric(unlist(bd)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(0)
  sum(vals)
}

get_section <- function(obj, base) {
  if (is.null(obj) || !is.list(obj)) return(NULL)
  sec <- obj[[base]] %||% obj[[paste0(base, "_schema")]] %||% NULL
  if (!is.list(sec)) return(NULL)
  sec
}

has_any_prediction <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(FALSE)
  main_sec <- get_section(obj, "main_type")
  sub1_sec <- get_section(obj, "subtype_level_1")
  sub2_sec <- get_section(obj, "subtype_level_2")
  cand <- c(
    get1(main_sec, "candidate_cell_type"),
    get1(sub1_sec, "candidate_cell_type"),
    get1(sub2_sec, "core_identity", "candidate_cell_type")
  )
  cand <- as.character(cand)
  cand <- cand[!is.na(cand)]
  cand <- stringr::str_trim(cand)
  any(nzchar(cand))
}

is_rejection_output <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(FALSE)
  st <- as.character(obj$evaluation_status %||% NA_character_)
  if (is.na(st)) return(FALSE)
  stringr::str_detect(st, "^Rejected")
}

rejection_label <- function(obj) {
  art <- as.character(obj$probable_artifact_type %||% NA_character_)
  st  <- as.character(obj$evaluation_status %||% NA_character_)
  core <- dplyr::coalesce(art, st, "Unknown")
  paste0("Rejected: ", core)
}

parse_report_for_cluster <- function(cid) {
  p15 <- file.path(dir_step15, paste0(cid, sprintf("_round%d_step1.5_citation_fix_output.json", round_i)))
  p1  <- file.path(dir_step1,  paste0(cid, sprintf("_round%d_step1_report_output.json", round_i)))
  pfinal <- file.path(dir_final, paste0(cid, "_FINAL_passed.json"))
  src <- NA_character_
  obj <- NULL

  # Prefer final_passed if it contains predictions; otherwise fall back.
  if (file.exists(pfinal)) {
    of <- tryCatch(jsonlite::fromJSON(pfinal, simplifyVector = FALSE), error = function(e) NULL)
    if (has_any_prediction(of)) {
      obj <- of
      src <- "final_passed"
    } else if (is_rejection_output(of)) {
      obj <- of
      src <- "final_passed_rejected"
    }
  }

  if (is.null(obj) && file.exists(p15)) {
    o15 <- tryCatch(jsonlite::fromJSON(p15, simplifyVector = FALSE), error = function(e) NULL)
    if (has_any_prediction(o15)) {
      obj <- o15
      src <- "step1.5"
    } else if (is_rejection_output(o15)) {
      obj <- o15
      src <- "step1.5_rejected"
    }
  }

  if (is.null(obj) && file.exists(p1)) {
    o1 <- tryCatch(jsonlite::fromJSON(p1, simplifyVector = FALSE), error = function(e) NULL)
    if (has_any_prediction(o1)) {
      obj <- o1
      src <- "step1"
    } else if (is_rejection_output(o1)) {
      obj <- o1
      src <- "step1_rejected"
    }
  }

  if (is.null(obj) || !is.list(obj)) return(NULL)
  if (is.null(obj) || !is.list(obj)) return(NULL)

  main_sec <- get_section(obj, "main_type")
  sub1_sec <- get_section(obj, "subtype_level_1")
  sub2_sec <- get_section(obj, "subtype_level_2")

  main_pred_val <- if (is_rejection_output(obj)) rejection_label(obj) else as.character(get1(main_sec, "candidate_cell_type"))
  subtype1_pred_val <- if (is_rejection_output(obj)) NA_character_ else as.character(get1(sub1_sec, "candidate_cell_type"))
  subtype2_pred_val <- if (is_rejection_output(obj)) NA_character_ else as.character(get1(sub2_sec, "core_identity", "candidate_cell_type"))

  out <- tibble::tibble(
    cluster_id = as.character(cid),

    # provenance
    report_source = as.character(src),

    main_pred = main_pred_val,
    main_cell_ontology_id = as.character(get1(main_sec, "cell_ontology_id")),
    main_confidence = get_conf_score(main_sec),

    subtype1_pred = subtype1_pred_val,
    subtype1_cell_ontology_id = as.character(get1(sub1_sec, "cell_ontology_id")),
    subtype1_confidence = get_conf_score(sub1_sec),

    subtype2_core_pred = subtype2_pred_val,
    subtype2_core_cell_ontology_id = as.character(get1(sub2_sec, "core_identity", "cell_ontology_id")),
    subtype2_confidence = get_conf_score(sub2_sec)
  )

  # allow subtype confidence to exceed main when evidence is stronger
  if (is.na(out$main_confidence)) out$main_confidence <- 0
  if (is.na(out$subtype1_confidence)) out$subtype1_confidence <- 0
  if (is.na(out$subtype2_confidence)) out$subtype2_confidence <- 0

  out
}

read_step2_status <- function(cid) {
  if (!dir.exists(dir_step2)) return(tibble::tibble(
    cluster_id = as.character(cid),
    step2_validation_status = NA_character_
  ))
  p2 <- file.path(dir_step2, paste0(cid, sprintf("_round%d_step2_validation_output.json", round_i)))
  if (!file.exists(p2)) return(tibble::tibble(
    cluster_id = as.character(cid),
    step2_validation_status = NA_character_
  ))
  obj <- tryCatch(jsonlite::fromJSON(p2, simplifyVector = FALSE), error = function(e) NULL)
  tibble::tibble(
    cluster_id = as.character(cid),
    step2_validation_status = as.character(obj$validation_status %||% NA_character_)
  )
}

read_gate <- function(cid) {
  # Gate JSON is written into step1.5_citation_fix_outputs
  pG <- file.path(dir_step15, paste0(cid, sprintf("_round%d_GATE.json", round_i)))
  if (!file.exists(pG)) return(tibble::tibble(
    cluster_id = as.character(cid),
    gate_passed = NA,
    gate_invalid_pmids_n = NA_integer_,
    gate_pmid_overuse_n = NA_integer_
  ))
  g <- tryCatch(jsonlite::fromJSON(pG, simplifyVector = FALSE), error = function(e) NULL)
  tibble::tibble(
    cluster_id = as.character(cid),
    gate_passed = as.logical(g$gate_passed %||% NA),
    gate_invalid_pmids_n = length(g$invalid_pmids %||% character(0)),
    gate_pmid_overuse_n = length(g$pmid_overuse %||% character(0))
  )
}

rows <- purrr::map(clusters, parse_report_for_cluster)
rows <- rows[!vapply(rows, is.null, logical(1))]
if (length(rows) == 0) stop("Failed to parse any Step1.5/final reports.")

df_rep <- dplyr::bind_rows(rows)
df_step2 <- dplyr::bind_rows(purrr::map(clusters, read_step2_status))
df_gate  <- dplyr::bind_rows(purrr::map(clusters, read_gate))

summary_df <- df_rep %>%
  dplyr::left_join(df_step2, by = "cluster_id") %>%
  dplyr::left_join(df_gate, by = "cluster_id")

# Optional: include Step6 run summary if present
run_sum_path <- file.path(root, "run_summary.csv")
if (file.exists(run_sum_path)) {
  rs <- tryCatch(readr::read_csv(run_sum_path, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(rs) && ("cluster" %in% names(rs))) {
    rs2 <- rs %>%
      dplyr::transmute(
        cluster_id = as.character(.data$cluster),
        step6_final_pass = as.logical(.data$final_pass %||% NA),
        step6_fail_stage = as.character(.data$fail_stage %||% NA_character_)
      )
    summary_df <- summary_df %>% dplyr::left_join(rs2, by = "cluster_id")
  }
}

summary_df <- summary_df %>% dplyr::arrange(cluster_id)

out_path <- file.path(out_dir, "summary.csv")
readr::write_csv(summary_df, out_path)

cat("[OK] wrote post summary: ", out_path, "\n", sep = "")
cat("Rows: ", nrow(summary_df), "\n", sep = "")
