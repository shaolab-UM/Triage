#!/usr/bin/env Rscript

# Deterministic reproduction of the selector controls reported in Table S7 / Fig. 3C-D.
# No LLM/API calls are made.
#
# Inputs:
#   data/primary/selector_inputs.tsv
#   Cell Ontology JSON used by the manuscript (release 2025-07-30)
#
# Outputs:
#   selector_cluster_results.tsv
#   selector_summary.tsv
#
# Important:
#   The retrospective oracle is the cluster-wise maximum across CASSIA,
#   In-house, clusterProfiler and Triage CL similarity to the reference.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1L]])
  default
}

repo_root <- get_arg("--repo-root", Sys.getenv("TRIAGE_HOME", unset = getwd()))
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)

input_file <- get_arg(
  "--input",
  file.path(repo_root, "data", "primary", "selector_inputs.tsv")
)
cl_json <- get_arg(
  "--cl-json",
  Sys.getenv(
    "CL_LOCAL_JSON",
    unset = file.path(repo_root, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json")
  )
)
out_dir <- get_arg(
  "--out-dir",
  file.path(repo_root, "results", "selector_controls")
)

if (!file.exists(input_file)) stop("Missing selector input: ", input_file)
if (!file.exists(cl_json)) stop(
  "Missing Cell Ontology JSON: ", cl_json,
  "\nProvide --cl-json or set CL_LOCAL_JSON."
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "lib", "00_cl_similarity.R"))
clsim <- build_cl_similarity(cl_json)
cl <- clsim$cl

dat <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)

if (nrow(dat) != 50L) stop("Expected 50 primary benchmark clusters; found ", nrow(dat))

pct_rank <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  out[ok] <- rank(x[ok], ties.method = "average") / sum(ok) * 100
  out
}

dat$cassia_percentile <- pct_rank(dat$cassia_native_score)
dat$in_house_percentile <- pct_rank(dat$in_house_native_score)
dat$clusterprofiler_percentile <- pct_rank(dat$clusterprofiler_native_score)

valid_cl <- function(id) {
  !is.na(id) && nzchar(id) && !is.null(cl[[id]])
}

score_cl <- function(pred, ref) {
  clsim$score_similarity(pred, ref)
}

majority_id <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) < 2L) return(NA_character_)
  tab <- table(ids)
  mx <- max(tab)
  winners <- names(tab)[tab == mx]
  if (mx >= 2L && length(winners) == 1L) winners[[1L]] else NA_character_
}

# Ontology helpers used by the deterministic ontology-only control.
depth_to_root <- function(id) {
  if (!valid_cl(id)) return(0)
  anc <- cl[[id]]$ancestors
  if (is.null(anc) || length(anc) == 0L) return(0)
  vals <- suppressWarnings(as.numeric(unlist(anc, use.names = FALSE)))
  if (length(vals) == 0L || all(!is.finite(vals))) return(0)
  max(vals[is.finite(vals)])
}

ancestor_ids_including_self <- function(id) {
  if (!valid_cl(id)) return(character(0))
  anc <- names(unlist(cl[[id]]$ancestors, use.names = TRUE))
  unique(c(id, anc))
}

is_descendant_or_same <- function(child, ancestor) {
  if (!valid_cl(child) || !valid_cl(ancestor)) return(FALSE)
  if (identical(child, ancestor)) return(TRUE)
  ancestor %in% names(unlist(cl[[child]]$ancestors, use.names = TRUE))
}

lca_two <- function(a, b) {
  if (!valid_cl(a) || !valid_cl(b)) return(NA_character_)
  common <- intersect(ancestor_ids_including_self(a), ancestor_ids_including_self(b))
  if (length(common) == 0L) return(NA_character_)
  d <- vapply(common, depth_to_root, numeric(1))
  # Deterministic tie break by CL ID; ties do not affect the reported dataset.
  ord <- order(-d, common)
  common[[ord[[1L]]]]
}

ontology_only <- function(ids, min_depth = 3L, k = 2L) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ids <- ids[vapply(ids, valid_cl, logical(1))]
  eligible <- ids[vapply(ids, depth_to_root, numeric(1)) >= min_depth]

  if (length(eligible) == 0L) {
    return(list(cl_id = NA_character_, decision = "Deferred"))
  }

  tab <- table(eligible)
  exact_anchor <- names(tab)[tab >= k]
  if (length(exact_anchor) > 0L) {
    d <- vapply(exact_anchor, depth_to_root, numeric(1))
    exact_anchor <- exact_anchor[order(-d, exact_anchor)]
    return(list(cl_id = exact_anchor[[1L]], decision = "Majority agreement"))
  }

  if (length(eligible) == 1L) {
    return(list(cl_id = eligible[[1L]], decision = "Single-reviewer anchor"))
  }

  pairs <- combn(seq_along(eligible), 2L)
  lcas <- character(0)
  for (j in seq_len(ncol(pairs))) {
    z <- lca_two(eligible[[pairs[1L, j]]], eligible[[pairs[2L, j]]])
    if (!is.na(z) && nzchar(z)) lcas <- c(lcas, z)
  }
  lcas <- unique(lcas)

  supported <- lcas[vapply(
    lcas,
    function(a) sum(vapply(eligible, is_descendant_or_same, logical(1), ancestor = a)) >= k,
    logical(1)
  )]

  if (length(supported) == 0L) {
    return(list(cl_id = NA_character_, decision = "Deferred"))
  }

  d <- vapply(supported, depth_to_root, numeric(1))
  supported <- supported[order(-d, supported)]
  list(cl_id = supported[[1L]], decision = "Supported common ancestor")
}

n <- nrow(dat)
dat$majority_vote_cl_id <- NA_character_
dat$majority_vote_similarity <- 0
dat$top_scoring_reviewer <- NA_character_
dat$top_scoring_cl_id <- NA_character_
dat$top_scoring_similarity <- 0
dat$ontology_only_cl_id <- NA_character_
dat$ontology_only_decision <- NA_character_
dat$ontology_only_similarity <- 0

dat$cassia_similarity <- 0
dat$in_house_similarity <- 0
dat$clusterprofiler_similarity <- 0
dat$triage_similarity <- 0
dat$published_retrospective_oracle_similarity <- 0

method_names <- c("CASSIA", "In-house", "clusterProfiler")

for (i in seq_len(n)) {
  ids <- c(dat$cassia_cl_id[[i]], dat$in_house_cl_id[[i]], dat$clusterprofiler_cl_id[[i]])
  pcts <- c(
    dat$cassia_percentile[[i]],
    dat$in_house_percentile[[i]],
    dat$clusterprofiler_percentile[[i]]
  )

  dat$cassia_similarity[[i]] <- score_cl(dat$cassia_cl_id[[i]], dat$reference_cl_id[[i]])
  dat$in_house_similarity[[i]] <- score_cl(dat$in_house_cl_id[[i]], dat$reference_cl_id[[i]])
  dat$clusterprofiler_similarity[[i]] <- score_cl(dat$clusterprofiler_cl_id[[i]], dat$reference_cl_id[[i]])
  dat$triage_similarity[[i]] <- score_cl(dat$triage_cl_id[[i]], dat$reference_cl_id[[i]])

  mid <- majority_id(ids)
  dat$majority_vote_cl_id[[i]] <- mid
  dat$majority_vote_similarity[[i]] <- score_cl(mid, dat$reference_cl_id[[i]])

  ok <- which(is.finite(pcts))
  if (length(ok) > 0L) {
    best <- ok[[which.max(pcts[ok])]]
    dat$top_scoring_reviewer[[i]] <- method_names[[best]]
    dat$top_scoring_cl_id[[i]] <- ids[[best]]
    dat$top_scoring_similarity[[i]] <- score_cl(ids[[best]], dat$reference_cl_id[[i]])
  }

  ont <- ontology_only(ids, min_depth = 3L, k = 2L)
  dat$ontology_only_cl_id[[i]] <- ont$cl_id
  dat$ontology_only_decision[[i]] <- ont$decision
  dat$ontology_only_similarity[[i]] <- score_cl(ont$cl_id, dat$reference_cl_id[[i]])
  # Retrospective best-available oracle:
  # maximum CL similarity across the three reviewer outputs plus Triage.
  # This is reference-using and is reported only as a retrospective ceiling.
  oracle_candidates <- c(
    dat$cassia_similarity[[i]],
    dat$in_house_similarity[[i]],
    dat$clusterprofiler_similarity[[i]],
    dat$triage_similarity[[i]]
  )
  dat$published_retrospective_oracle_similarity[[i]] <- max(oracle_candidates, na.rm = TRUE)
}

summarize_subset <- function(x, label) {
  data.frame(
    Relationship = label,
    n = nrow(x),
    Majority_vote = mean(x$majority_vote_similarity),
    Top_reviewer_percentile = mean(x$top_scoring_similarity),
    Ontology_only = mean(x$ontology_only_similarity),
    Triage = mean(x$triage_similarity),
    Published_retrospective_oracle = mean(x$published_retrospective_oracle_similarity),
    stringsAsFactors = FALSE
  )
}

rels <- c(
  "Incomplete reviewer resolution",
  "Reviewer consensus",
  "Supported refinement",
  "Unresolved disagreement"
)

summary_rows <- lapply(rels, function(r) summarize_subset(dat[dat$reviewer_relationship == r, , drop = FALSE], r))
summary_rows[[length(summary_rows) + 1L]] <- summarize_subset(dat, "Overall")
summary_tbl <- do.call(rbind, summary_rows)

write.table(
  dat,
  file = file.path(out_dir, "selector_cluster_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write.table(
  summary_tbl,
  file = file.path(out_dir, "selector_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("\nSelector-control reproduction\n")
print(summary_tbl, row.names = FALSE, digits = 8)

cat("\nDeferred / unavailable outputs\n")
cat("Majority vote:", sum(is.na(dat$majority_vote_cl_id)), "/ 50\n")
cat("Top-scoring reviewer with no mapped CL:", sum(is.na(dat$top_scoring_cl_id)), "/ 50\n")
cat("Ontology-only:", sum(is.na(dat$ontology_only_cl_id)), "/ 50\n")

expected <- c(
  Majority_vote = 41.3207,
  Top_reviewer_percentile = 61.6773,
  Ontology_only = 62.9211,
  Triage = 79.9790,
  Published_retrospective_oracle = 83.7732
)
overall <- summary_tbl[summary_tbl$Relationship == "Overall", ]
observed <- c(
  Majority_vote = overall$Majority_vote,
  Top_reviewer_percentile = overall$Top_reviewer_percentile,
  Ontology_only = overall$Ontology_only,
  Triage = overall$Triage,
  Published_retrospective_oracle = overall$Published_retrospective_oracle
)

if (any(abs(observed - expected) > 0.02)) {
  stop(
    "Selector reproduction differs from the publication values by >0.02.\n",
    paste(names(observed), sprintf("%.6f", observed), collapse = "\n")
  )
}

cat("\nPASS: published selector values reproduced within tolerance.\n")
