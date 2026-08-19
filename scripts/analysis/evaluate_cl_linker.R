#!/usr/bin/env Rscript

# CL-Linker publication evaluation.
#
# This script performs two deterministic checks:
#   1) Recomputes the 140-annotation operational evaluation from the released
#      reviewer mapping registry and manual reference file.
#   2) Recomputes the OLS/exact-synonym comparator outcomes on the 97-label
#      direct-mapping set from archived per-label predictions.
#
# The final 97-label CL-Linker direct-mapping row is retained as publication
# source data because the separate per-label output file for that direct
# evaluation was not present in the supplied project snapshot. No values are
# silently reconstructed or imputed.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(jsonlite)
})

option_list <- list(
  make_option("--repo-root", type="character",
              default=Sys.getenv("TRIAGE_HOME", unset=getwd())),
  make_option("--cl-json", type="character", default=""),
  make_option("--out-dir", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))

root <- normalizePath(opt$`repo-root`, winslash="/", mustWork=TRUE)
cl_json <- if (nzchar(opt$`cl-json`)) opt$`cl-json` else
  Sys.getenv("CL_LOCAL_JSON", unset=file.path(root, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"))
if (!file.exists(cl_json)) stop("Cell Ontology JSON not found: ", cl_json)

out_dir <- if (nzchar(opt$`out-dir`)) opt$`out-dir` else file.path(root, "results", "cl_linker_evaluation")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

gold_file <- file.path(root, "data", "cl_linker", "manual_mapping_gold.tsv")
registry_file <- file.path(root, "data", "primary", "reviewer_mapping_registry.tsv")
comp_file <- file.path(root, "data", "cl_linker", "mapping_comparator_predictions.tsv")
published_direct_file <- file.path(root, "data", "cl_linker", "publication_direct_mapping_summary.tsv")
published_oper_file <- file.path(root, "data", "cl_linker", "publication_operational_summary.tsv")

for (fp in c(gold_file, registry_file, comp_file, published_direct_file, published_oper_file)) {
  if (!file.exists(fp)) stop("Missing evaluation input: ", fp)
}

cl <- jsonlite::fromJSON(cl_json, simplifyVector=FALSE)
gold <- readr::read_tsv(gold_file, show_col_types=FALSE)
registry <- readr::read_tsv(registry_file, show_col_types=FALSE)
comp <- readr::read_tsv(comp_file, show_col_types=FALSE)
published_direct <- readr::read_tsv(published_direct_file, show_col_types=FALSE)
published_oper <- readr::read_tsv(published_oper_file, show_col_types=FALSE)

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- NA_character_
  x
}

ancestor_distance <- function(descendant, ancestor) {
  if (is.na(descendant) || is.na(ancestor) || is.null(cl[[descendant]])) return(Inf)
  anc <- cl[[descendant]]$ancestors
  if (is.null(anc) || is.null(anc[[ancestor]])) return(Inf)
  suppressWarnings(as.numeric(anc[[ancestor]]))
}

classify_mapping <- function(pred, ref) {
  pred <- safe_chr(pred)
  ref <- safe_chr(ref)
  if (is.na(pred) || is.na(ref)) return("unmapped")
  if (identical(pred, ref)) return("exact")

  # Predicted term is an ancestor of the reference.
  d_parent <- ancestor_distance(ref, pred)
  if (is.finite(d_parent) && d_parent == 1) return("direct_parent")

  # Predicted term is a descendant of the reference.
  d_desc <- ancestor_distance(pred, ref)
  if (is.finite(d_desc)) return("overspecific")

  "incorrect"
}

summarize_outcomes <- function(method, outcomes) {
  tibble::tibble(
    mapping_approach = method,
    n = length(outcomes),
    exact_n = sum(outcomes == "exact"),
    exact_or_direct_parent_n = sum(outcomes %in% c("exact", "direct_parent")),
    overspecific_n = sum(outcomes == "overspecific"),
    incorrect_n = sum(outcomes == "incorrect"),
    unmapped_n = sum(outcomes == "unmapped")
  )
}

# ------------------------------------------------------------------
# 97 unique single-identity labels: deterministic comparator check
# ------------------------------------------------------------------
single97 <- gold %>%
  filter(gt_status == "single_identity", !is.na(gt_cl_id), gt_cl_id != "") %>%
  distinct(raw_label, .keep_all=TRUE) %>%
  select(raw_label, gt_cl_id)

if (nrow(single97) != 97L) stop("Expected 97 unique single-identity labels; found ", nrow(single97))

cmp <- single97 %>%
  left_join(comp, by=c("raw_label","gt_cl_id"))

if (nrow(cmp) != 97L) stop("Comparator merge did not preserve 97 labels")

exact_syn_out <- mapply(classify_mapping, cmp$exact_syn_clid, cmp$gt_cl_id, USE.NAMES=FALSE)
ols_exact_out <- mapply(classify_mapping, cmp$ols_exact_clid, cmp$gt_cl_id, USE.NAMES=FALSE)
ols_search_out <- mapply(classify_mapping, cmp$ols_search_clid, cmp$gt_cl_id, USE.NAMES=FALSE)

computed_comparators <- bind_rows(
  summarize_outcomes("OLS top-result search", ols_search_out),
  summarize_outcomes("OLS exact lookup", ols_exact_out),
  summarize_outcomes("Exact/synonym lookup", exact_syn_out)
)

expected_comparators <- published_direct %>%
  filter(mapping_approach != "CL-Linker")

check_cols <- c("n","exact_n","exact_or_direct_parent_n","overspecific_n","incorrect_n","unmapped_n")
cc <- computed_comparators %>%
  arrange(mapping_approach)
ee <- expected_comparators %>%
  arrange(mapping_approach)

if (!identical(cc$mapping_approach, ee$mapping_approach) ||
    any(as.matrix(cc[check_cols]) != as.matrix(ee[check_cols]))) {
  stop("Direct-mapping comparator results do not match the publication source data.")
}

# ------------------------------------------------------------------
# 140 single-identity reviewer annotations: operational CL-Linker
# ------------------------------------------------------------------
join_keys <- c("dataset","cluster_id","reviewer","raw_label")
op <- gold %>%
  left_join(
    registry %>% select(all_of(join_keys), cl_id, single_cl_ready, final_status, mapping_method),
    by=join_keys
  )

single140 <- op %>% filter(gt_status == "single_identity")
if (nrow(single140) != 140L) stop("Expected 140 single-identity reviewer annotations; found ", nrow(single140))

op_outcomes <- mapply(classify_mapping, single140$cl_id, single140$gt_cl_id, USE.NAMES=FALSE)

operational <- tibble::tibble(
  metric = c(
    "single_identity_reviewer_annotations",
    "exact_mapping",
    "direct_parent_mapping",
    "exact_or_direct_parent_mapping",
    "overspecific_mapping",
    "incorrect_mapping",
    "unmapped_or_abstained",
    "fraction_mapped_numerator",
    "fraction_mapped_denominator",
    "non_single_reference_withheld"
  ),
  value = c(
    length(op_outcomes),
    sum(op_outcomes == "exact"),
    sum(op_outcomes == "direct_parent"),
    sum(op_outcomes %in% c("exact","direct_parent")),
    sum(op_outcomes == "overspecific"),
    sum(op_outcomes == "incorrect"),
    sum(op_outcomes == "unmapped"),
    sum(op_outcomes != "unmapped"),
    length(op_outcomes),
    sum(op$gt_status != "single_identity" & (is.na(op$cl_id) | op$cl_id == ""))
  )
)

expected_oper <- published_oper %>% mutate(value=as.numeric(value))
if (!identical(operational$metric, expected_oper$metric) ||
    any(operational$value != expected_oper$value)) {
  stop("Operational CL-Linker evaluation does not match the publication source data.")
}

# ------------------------------------------------------------------
# Write release-facing results
# ------------------------------------------------------------------
readr::write_tsv(computed_comparators, file.path(out_dir, "direct_mapping_comparators_recomputed.tsv"))
readr::write_tsv(operational, file.path(out_dir, "operational_cl_linker_recomputed.tsv"))
readr::write_tsv(published_direct, file.path(out_dir, "publication_direct_mapping_summary.tsv"))

cat("\nCL-Linker evaluation checks\n")
cat("---------------------------\n")
print(computed_comparators)
cat("\nOperational evaluation:\n")
print(operational)

cat("\nPASS: OLS/exact-synonym comparator outcomes and the 140-annotation operational CL-Linker evaluation match the publication source data.\n")
cat("NOTE: the 97-label CL-Linker direct-mapping row is distributed as publication source data; its separate per-label direct-evaluation output was not available in the supplied project snapshot.\n")
