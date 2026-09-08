#!/usr/bin/env Rscript
# ============================================================
# 06b_run_inter.R — run enrichment-based annotation (third reviewer)
# Input: 05c_llm_queries (step1 queries)
# Output: outputs/<run>/06b_inter/ (enrichment results)
# ============================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(rio)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(future)
  library(furrr)
})

triageHome <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome)) {
  scriptArgV <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome <- if (length(scriptArgV) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
suppressMessages(library(Triage))

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1) {
    if (is.na(a)) return(b)
    if (is.character(a) && !nzchar(a)) return(b)
  }
  a
}

ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(Triage))

# The enrichment reviewer (stage 06b) requires the interpret() API that only
# exists from clusterProfiler 4.19.4 onward, and the tested execution route
# is the fanyi-backed interpret() at the exact revision below. Never install
# a moving HEAD: modern revisions rerouted interpret() through aisdk.
TRIAGE_TESTED_CLUSTERPROFILER_SHA <- "f9f0d502508cacd258ac1a1cba6d5d497b98fe6c"

clusterprofiler_sha_matches <- function() {
  desc <- tryCatch(utils::packageDescription("clusterProfiler"), error = function(e) NULL)
  if (is.null(desc)) return(FALSE)
  sha <- desc$RemoteSha %||% desc$GithubSHA1 %||% ""
  nzchar(sha) && identical(sha, TRIAGE_TESTED_CLUSTERPROFILER_SHA)
}

ensure_clusterprofiler_github <- function() {
  if (clusterprofiler_sha_matches()) return(invisible(TRUE))
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("YuLab-SMU/clusterProfiler",
                          ref = TRIAGE_TESTED_CLUSTERPROFILER_SHA,
                          upgrade = "never")
  if (!clusterprofiler_sha_matches()) {
    stop("06b_run_inter: clusterProfiler is not the tested revision ",
         "(GitHub sha ", TRIAGE_TESTED_CLUSTERPROFILER_SHA,
         "). Install it with: remotes::install_github(\"YuLab-SMU/clusterProfiler\", ref = \"",
         TRIAGE_TESTED_CLUSTERPROFILER_SHA, "\")")
  }
  invisible(TRUE)
}

option_list <- list(
  make_option("--marker_csv", type = "character", default = normalizePath("maskdeg.csv", mustWork = FALSE)),
  make_option("--step1_dir", type = "character", default = NULL),
  make_option("--bioinfo_dir", type = "character", default = NULL),
  make_option("--out_dir", type = "character", default = NULL),
  make_option("--dataset_name", type = "character", default = NULL),
  make_option("--model", type = "character", default = Sys.getenv("LLM_MODEL_INTER", unset = Sys.getenv("CASSIA_MODEL_INTER", unset = "deepseek-v4-flash"))),
  make_option("--context", type = "character", default = NULL),
  make_option("--n_markers", type = "integer", default = 50),
  make_option("--cluster_id", type = "character", default = ""),
  make_option("--workers", type = "integer", default = 12),
  make_option("--species", type = "character", default = "human"),
  make_option("--cellmarker_human_path", type = "character", default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "cellmarker", "Cell_marker_Human.xlsx")),
  make_option("--cellmarker_mouse_path", type = "character", default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "cellmarker", "Cell_marker_Mouse.xlsx")),
  make_option("--install_clusterprofiler_github", action = "store_true", default = FALSE),
  make_option("--cl_local_json", type = "character", default = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json")),
  make_option("--ols_first", type = "logical", default = TRUE),
  make_option("--ols_cache_dir", type = "character", default = "")
)

opt <- parse_args(OptionParser(option_list = option_list))
dataset_name <- if (nzchar(opt$dataset_name %||% "")) opt$dataset_name else basename(getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
cfg <- Triage:::get_dataset_config(dataset_name, project_root)
if (is.null(opt$step1_dir) || !nzchar(opt$step1_dir)) {
  opt$step1_dir <- file.path(cfg$llm_inputs_root, "step1_report_queries")
}
if (is.null(opt$bioinfo_dir) || !nzchar(opt$bioinfo_dir)) {
  opt$bioinfo_dir <- file.path(cfg$intermediate_outputs_root, "bioinformatics_tsv")
  if (!dir.exists(opt$bioinfo_dir)) {
    legacy <- file.path(project_root, "intermediate_outputs", paste0(dataset_name, "_LLM_Input_Run"), "bioinformatics_tsv")
    if (dir.exists(legacy)) opt$bioinfo_dir <- legacy
  }
}
if (is.null(opt$out_dir) || !nzchar(opt$out_dir)) {
  opt$out_dir <- file.path(cfg$llm_outputs_root, "inter")
}
if (is.null(opt$context) || !nzchar(opt$context)) {
  opt$context <- cfg$user_notes
}

default_marker_csv <- normalizePath("maskdeg.csv", mustWork = FALSE)
if (is.null(opt$marker_csv) || !nzchar(opt$marker_csv) || identical(opt$marker_csv, default_marker_csv)) {
  if (file.exists(cfg$deg_file)) opt$marker_csv <- cfg$deg_file
}
if (isTRUE(opt$install_clusterprofiler_github)) {
  ensure_clusterprofiler_github()
} else if (!clusterprofiler_sha_matches()) {
  warning("clusterProfiler is not the tested revision (GitHub sha ",
          TRIAGE_TESTED_CLUSTERPROFILER_SHA, "). The interpret() API used by ",
          "the enrichment reviewer requires clusterProfiler >= 4.19.4 with the ",
          "fanyi route; run with --install_clusterprofiler_github when GitHub ",
          "access is available.")
}
ensure_dir(opt$out_dir)
ensure_dir(file.path(opt$out_dir, "final_passed"))
ensure_dir(file.path(opt$out_dir, "debug_failed"))

ols_cache_dir <- opt$ols_cache_dir
if (!nzchar(ols_cache_dir)) ols_cache_dir <- file.path(opt$out_dir, ".ols_cache")
cl_cfg <- Triage:::make_cl_cfg(opt$cl_local_json, prefer_ols = isTRUE(opt$ols_first), cache_dir = ols_cache_dir)

get_step1_cluster_ids <- function(step1_dir) {
  if (!dir.exists(step1_dir)) return(character(0))
  fs <- list.files(step1_dir, pattern = "_step1_report_query\\.json$", full.names = FALSE)
  if (length(fs) == 0) return(character(0))
  sub("_step1_report_query\\.json$", "", fs)
}

load_cellmarker_db <- function(species, opt) {
  cm_path <- if (tolower(species) == "mouse") opt$cellmarker_mouse_path else opt$cellmarker_human_path
  if (!file.exists(cm_path)) stop("CellMarker file not found: ", cm_path)
  cm <- rio::import(cm_path)
  req_cols <- c("cell_name", "marker")
  miss <- setdiff(req_cols, names(cm))
  if (length(miss) > 0) stop("CellMarker missing columns: ", paste(miss, collapse = ", "))
  cm[, req_cols]
}

load_top_markers <- function(marker_csv, cluster_id, n_markers) {
  if (!file.exists(marker_csv)) stop("marker_csv not found: ", marker_csv)
  deg <- readr::read_csv(marker_csv, show_col_types = FALSE)
  req <- c("cluster", "gene", "avg_log2FC")
  miss <- setdiff(req, names(deg))
  if (length(miss) > 0) stop("marker_csv missing columns: ", paste(miss, collapse = ", "))
  deg %>%
    dplyr::filter(cluster == cluster_id) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::pull(gene) %>%
    unique() %>%
    utils::head(n_markers)
}

read_bioinfo_cache <- function(cluster_id, bioinfo_dir, pattern) {
  fp <- file.path(bioinfo_dir, cluster_id, pattern)
  if (!file.exists(fp)) return(NULL)
  readr::read_tsv(fp, show_col_types = FALSE, col_names = FALSE)
}

make_compare_cluster_result_from_cache <- function(df, cluster_id, fun, organism, genes) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (ncol(df) < 12) return(NULL)
  colnames(df)[1:12] <- c(
    "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
    "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
    "geneID", "Count"
  )
  df$Cluster <- cluster_id
  df$cluster <- cluster_id
  df <- df[, c(
    "Cluster", "cluster", "ID", "Description", "GeneRatio", "BgRatio",
    "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust",
    "qvalue", "geneID", "Count"
  )]
  methods::new(
    "compareClusterResult",
    compareClusterResult = df,
    geneClusters = setNames(list(genes), cluster_id),
    fun = fun,
    gene2Symbol = character(0),
    keytype = "SYMBOL",
    readable = FALSE,
    .call = quote(clusterProfiler::compareCluster),
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = "ORA",
    dr = list(),
    organism = organism
  )
}

build_enrichment_objects <- function(genes, cluster_id, opt) {
  gene_df <- tibble::tibble(gene = genes, cluster = cluster_id)
  cm <- load_cellmarker_db(opt$species, opt)

  x_cell <- clusterProfiler::compareCluster(
    gene ~ cluster,
    data = gene_df,
    fun = clusterProfiler::enricher,
    TERM2GENE = cm
  )

  go_df <- read_bioinfo_cache(cluster_id, opt$bioinfo_dir, "go_biological_process_full_results.tsv")
  kegg_df <- read_bioinfo_cache(cluster_id, opt$bioinfo_dir, "kegg_pathways_full_results.tsv")

  organism <- if (tolower(opt$species) == "mouse") "mouse" else "human"
  x_go <- make_compare_cluster_result_from_cache(go_df, cluster_id, "enrichGO", organism, genes)
  x_kegg <- make_compare_cluster_result_from_cache(kegg_df, cluster_id, "enrichKEGG", organism, genes)

  list(x_cell, x_go, x_kegg)
}

normalize_key <- function(x) {
  stringr::str_replace_all(tolower(x), "[^a-z0-9]", "")
}

get_field <- function(res, keys) {
  if (is.null(res) || is.null(names(res))) return(NULL)
  nms <- names(res)
  norm <- normalize_key(nms)
  for (k in keys) {
    idx <- which(norm == normalize_key(k))
    if (length(idx) > 0) return(res[[nms[[idx[[1]]]]]])
  }
  NULL
}

parse_inter_output <- function(res) {
  out <- list()
  if (is.character(res)) {
    txt <- res[[1]]
    out$raw_text <- txt
    parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(parsed) && is.list(parsed)) res <- parsed
  }
  if (is.list(res)) {
    if (length(res) == 1 && is.list(res[[1]])) {
      res <- res[[1]]
    }
    out$cell_type <- get_field(res, c("cell_type", "celltype", "cell type", "cell-type", "celltype_label", "celltype_label", "prediction", "predicted_label", "label"))
    out$confidence <- get_field(res, c("confidence", "score", "probability"))
    out$reasoning <- get_field(res, c("reasoning", "rationale", "explanation", "overview"))
    out$markers <- get_field(res, c("markers", "marker_genes", "genes"))
    out$terms <- get_field(res, c("evidence_terms", "terms", "pathways"))
    if (is.null(out$raw_text)) {
      out$raw_text <- tryCatch(jsonlite::toJSON(res, auto_unbox = TRUE, pretty = TRUE), error = function(e) NULL)
    }
  }
  out
}

write_debug <- function(cluster_id, err_msg, raw_text, out_dir) {
  payload <- list(
    cluster_id = cluster_id,
    error = err_msg,
    raw_text = raw_text
  )
  jsonlite::write_json(payload, file.path(out_dir, "debug_failed", paste0(cluster_id, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
}

# ---- INTER ----
cluster_ids <- get_step1_cluster_ids(opt$step1_dir)
if (nzchar(opt$cluster_id)) cluster_ids <- cluster_ids[cluster_ids == opt$cluster_id]
if (length(cluster_ids) == 0) stop("No clusters found to process.")

plan(multisession, workers = opt$workers)

process_cluster <- function(cid, opt, cl_cfg) {
  genes <- load_top_markers(opt$marker_csv, cid, opt$n_markers)
  if (length(genes) == 0) {
    write_debug(cid, "No marker genes found", NA_character_, opt$out_dir)
    return(NULL)
  }

  api_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = NA_character_)
  # Canonical LLM transport for the enrichment reviewer: the same full
  # chat-completions endpoint (LLM_API_BASE_URL) and the same model
  # (opt$model) used by every other Triage LLM stage. The endpoint/model
  # reach clusterProfiler::interpret() through Triage:::.triage_interpret().
  base_url <- Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = ""))
  if (is.na(api_key) || !nzchar(api_key)) api_key <- NULL
  if (!is.na(api_key) && nzchar(api_key)) {
    options(yulab_translate = list(dsk = list(key = api_key, user_model = opt$model)))
    if (requireNamespace("fanyi", quietly = TRUE)) {
      tryCatch({
        fanyi::set_translate_option(source = "dsk", key = api_key, user_model = opt$model)
      }, error = function(e) {
        warning("Failed to set fanyi translate options: ", e$message)
      })
    }
  }

  # Localize: replace fanyi::search_gene with org.db（avoid an EBI API dependency）
  if (requireNamespace("fanyi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
      assignInNamespace("search_gene", function(x, organism = "Homo sapiens") {
        if (tolower(organism) %in% c("homo sapiens", "human")) {
          suppressMessages({
            map <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(x),
                                         columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
          })
          map <- map[!is.na(map$ENTREZID), ]
          res <- data.frame(SYMBOL = map$SYMBOL, ENTREZID = map$ENTREZID, stringsAsFactors = FALSE)
          return(res)
        }
        fanyi::search_gene(x, organism)
      }, ns = "fanyi")
      cat("[LOCAL] fanyi::search_gene replaced with org.Hs.eg.db; no EBI dependency\n")
    }, error = function(e) {
      warning("Local search_gene replacement failed: ", e$message)
    })
  }

  context_text <- paste(
    opt$context,
    paste0("Top markers: ", paste(genes, collapse = ", ")),
    sep = "\n"
  )

  res <- tryCatch({
    enrich_list <- build_enrichment_objects(genes, cid, opt)
    enrich_list <- Filter(Negate(is.null), enrich_list)
    if (length(enrich_list) == 0) stop("No enrichment results available")
    # LLM callsretry（up to three attempts，for transient API errors or timeouts）
    interpret_out <- NULL
    last_err <- NULL
    for (attempt in 1:3) {
      tryCatch({
        interpret_out <- Triage:::.triage_interpret(
          enrich_list, context = context_text, model = opt$model,
          api_key = api_key, base_url = base_url, task = "annotation")
        break
      }, error = function(e) {
        last_err <<- e
        cat(sprintf("  [WARN] interpret attempt %d failed: %s | retrying in 5 s\n", attempt, conditionMessage(e)))
        Sys.sleep(5)
      })
    }
    if (is.null(interpret_out) && !is.null(last_err)) {
      stop("interpret retry failed after three attempts: ", conditionMessage(last_err))
    }
    interpret_out
  }, error = function(e) e)

  if (inherits(res, "error")) {
    write_debug(cid, res$message, NA_character_, opt$out_dir)
    return(NULL)
  }

  parsed <- parse_inter_output(res)
  label <- as.character(parsed$cell_type %||% "")
  if (!nzchar(label)) {
    write_debug(cid, "Missing cell_type in interpret output", parsed$raw_text %||% NA_character_, opt$out_dir)
    return(NULL)
  }

  cl_norm <- Triage:::normalize_cl_three_state(label, "", cl_cfg)
  terms <- parsed$terms %||% character(0)
  markers <- parsed$markers %||% character(0)

  out <- list(
    cluster_id = cid,
    predicted_label = cl_norm$final_name %||% label,
    cl_id = cl_norm$final_clid %||% NULL,
    confidence = parsed$confidence %||% NULL,
    rationale = parsed$reasoning %||% NULL,
    evidence_markers = markers,
    evidence_terms = terms
  )
  jsonlite::write_json(out, file.path(opt$out_dir, "final_passed", paste0(cid, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
  out
}

furrr::future_map(cluster_ids, process_cluster, opt = opt, cl_cfg = cl_cfg, .options = furrr::furrr_options(seed = TRUE))
plan(sequential)
