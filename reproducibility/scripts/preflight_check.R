#!/usr/bin/env Rscript
# preflight_check.R — full-workflow preflight validation
#
# Reports ALL missing prerequisites for running pipeline stages 03a-10 in one
# pass, BEFORE any analysis starts. Exits non-zero when required items are
# missing. Checks only; performs no analysis and makes no LLM calls.
#
# Usage:
#   Rscript reproducibility/scripts/preflight_check.R --species human
#   Rscript reproducibility/scripts/preflight_check.R --species mouse [--deg markers.csv]

suppressMessages({
  library(optparse); library(jsonlite)
})

option_list <- list(
  make_option("--species", type = "character", default = "human", help = "human or mouse"),
  make_option("--deg", type = "character", default = NULL, help = "Optional user DEG table to column-validate"),
  make_option("--skip-internet", type = "logical", default = FALSE, action = "true",
              help = "Skip the endpoint reachability probe")
)
opt <- parse_args(OptionParser(option_list = option_list))

species <- tolower(opt$species)
if (!species %in% c("human", "mouse")) {
  cat("PREFLIGHT: unsupported species '", opt$species, "' (use human or mouse)\n", sep = "")
  quit(status = 2)
}

triage_home <- Sys.getenv("TRIAGE_HOME", unset = getwd())
project_root <- Sys.getenv("PROJECT_ROOT", unset = triage_home)

missing <- character(0)
warn <- character(0)
ok <- function(x) cat("  OK   ", x, "\n")
bad <- function(x) { cat("  MISS ", x, "\n"); missing <<- c(missing, x) }
soft <- function(x) { cat("  WARN ", x, "\n"); warn <<- c(warn, x) }

cat("=== Triage full-workflow preflight ===\n")
cat("species:", species, "\n")
cat("TRIAGE_HOME:", triage_home, "\n\n")

# --- 1. R packages required by stages 03a-10 -------------------------------
pkg_sets <- list(
  "core (03a-11)" = c("optparse", "jsonlite", "dplyr", "readr", "stringr",
                      "purrr", "tibble", "digest", "writexl", "rlang"),
  "evidence/LLM (05)" = c("httr", "glue", "data.table", "tidyr", "memoise",
                          "cachem", "knitr", "readxl", "tictoc", "future",
                          "future.apply", "xml2", "fs"),
  "in-house reviewer (06/06b)" = c("furrr", "rio", "AnnotationDbi",
                                   "org.Hs.eg.db"),
  "enrichment (05/06b)" = c("clusterProfiler", "DOSE", "ReactomePA", "enrichR",
                            "decoupleR")
)
cat("R packages:\n")
for (nm in names(pkg_sets)) {
  miss <- pkg_sets[[nm]][!vapply(pkg_sets[[nm]], function(p)
    requireNamespace(p, quietly = TRUE), logical(1))]
  if (length(miss) > 0) bad(paste0("R packages [", nm, "]: ", paste(miss, collapse = ", ")))
  else ok(paste0("R packages [", nm, "]"))
}
if (!requireNamespace("Triage", quietly = TRUE)) {
  bad("R package: Triage (install from this repository)")
} else {
  ok("R package: Triage")
}
if (!requireNamespace("CASSIA", quietly = TRUE)) {
  bad("R package: CASSIA (stage 03b; install separately — not distributed with Triage)")
} else {
  ok("R package: CASSIA")
}
if (requireNamespace("CASSIA", quietly = TRUE)) {
  # Version-compatible lookup: check_python_env is present in the tested
  # CASSIA revision (b008c0ac); getFromNamespace works whether or not the
  # historical revision exports it.
  py_ok <- tryCatch({
    check_fun <- utils::getFromNamespace("check_python_env", "CASSIA")
    suppressWarnings(isTRUE(check_fun()))
  }, error = function(e) FALSE)
  if (isTRUE(py_ok)) {
    ok("CASSIA Python backend")
  } else {
    bad(paste0("CASSIA Python backend is not usable; run CASSIA::setup_cassia_env() ",
               "once (one-time R command; creates the Python environment CASSIA needs)"))
  }
  cassia_sha <- tryCatch(packageDescription("CASSIA")$RemoteSha,
                         error = function(e) NA_character_)
  if (is.null(cassia_sha) || is.na(cassia_sha)) {
    cassia_sha <- tryCatch(packageDescription("CASSIA")$GithubSHA1,
                           error = function(e) NA_character_)
  }
  cassia_tested <- "b008c0ac3dd81b2c2dff131d20f5081a58aca027"
  if (!is.na(cassia_sha) && cassia_sha != cassia_tested) {
    soft(paste0("CASSIA revision ", substr(cassia_sha, 1, 12),
                " differs from the tested revision ",
                substr(cassia_tested, 1, 12),
                "; reinstall with install_triage_dependencies() for the canonical configuration"))
  } else {
    ok("CASSIA revision matches the tested revision")
  }
}
# DisGeNET is optional: disgenet2r is required only when its API key is set.
if (nzchar(Sys.getenv("DISGENET_API_KEY", unset = "")) &&
    !requireNamespace("disgenet2r", quietly = TRUE)) {
  bad(paste0("disgenet2r is required when DISGENET_API_KEY is set; install with: ",
             "remotes::install_gitlab(\"medbio/disgenet2r\")"))
}
# Canonical KEGG evidence (stage 05 enrichKEGG use_internal_data=TRUE) is
# backed by KEGG.db (removed from current Bioconductor release).
if (!requireNamespace("KEGG.db", quietly = TRUE)) {
  bad(paste0("R package: KEGG.db (canonical KEGG evidence stage 05); install with install_triage_dependencies() ",
             "or remotes::install_url(\"https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz\")"))
} else if (!tryCatch(length(AnnotationDbi::keys(get("KEGGPATHID2EXTID",
                                                    envir = asNamespace("KEGG.db")))) > 0,
                     error = function(e) FALSE)) {
  bad("R package: KEGG.db is installed but its KEGG data cannot be loaded")
} else {
  ok("R package: KEGG.db (canonical KEGG evidence)")
}
# Canonical Reactome evidence (stage 05 Reactome local enrichment).
if (!requireNamespace("reactome.db", quietly = TRUE)) {
  bad("R package: reactome.db (Reactome enrichment stage 05; install with install_triage_dependencies())")
} else {
  ok("R package: reactome.db")
}
soft("Seurat is required only for 01a --mode seurat (not for the CSV workflow)")
if (species == "mouse" && !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  bad("R package: org.Mm.eg.db (required for mouse evidence analysis; install with BiocManager::install(\"org.Mm.eg.db\"))")
}

# --- 2. API credentials -----------------------------------------------------
key_env <- Sys.getenv("LLM_API_KEY_ENV", unset = "DEEPSEEK_API_KEY")
api_key <- Sys.getenv(key_env, unset = "")
base_url <- Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = ""))
if (!nzchar(api_key) || identical(api_key, "XXXXX")) {
  bad(paste0("API key env ", key_env, " (empty or placeholder)"))
} else {
  ok("API key")
}
if (!nzchar(base_url) || identical(base_url, "XXXXX")) {
  bad("LLM_API_BASE_URL (full chat-completions endpoint)")
} else if (!grepl("/chat/completions$", base_url)) {
  bad(paste0("LLM_API_BASE_URL must be the FULL chat-completions endpoint (ends with /chat/completions); got: ", base_url))
} else {
  ok(paste0("LLM_API_BASE_URL: ", base_url))
}

# --- 3. Bundled Cell Ontology ------------------------------------------------
cl_json <- Sys.getenv("CL_LOCAL_JSON", unset = file.path(triage_home, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"))
if (file.exists(cl_json)) {
  ontology_ok <- tryCatch({
    parsed <- jsonlite::fromJSON(cl_json, simplifyVector = FALSE)
    is.list(parsed) || is.data.frame(parsed)
  }, error = function(e) FALSE)
  if (ontology_ok) ok(paste0("Cell Ontology JSON parses: ", cl_json))
  else bad(paste0("Cell Ontology JSON is unreadable or invalid: ", cl_json))
} else {
  bad(paste0("Cell Ontology JSON not found at ", cl_json))
}

# --- 4. External resources (species-specific) --------------------------------
check_string_file <- function(path, label) {
  if (!file.exists(path)) {
    bad(paste0("STRING ", label, " file: ", path,
               " (download per resources/README.md; or set TRIAGE_PPI_ROOT)"))
    return(invisible(NULL))
  }
  if (file.info(path)$size <= 0L || file.access(path, 4L) != 0L) {
    bad(paste0("STRING ", label, " file is unreadable or empty: ", path))
    return(invisible(NULL))
  }
  fields <- tryCatch({
    lines <- readLines(path, n = 20L, warn = FALSE)
    line <- lines[nzchar(trimws(lines))][1]
    strsplit(trimws(line), "\\s+")[[1]]
  }, error = function(e) character(0))
  if (length(fields) != 3L) {
    bad(paste0("STRING ", label,
               " file must have three whitespace-delimited columns for the pipeline reader: ", path))
  } else {
    ok(paste0("STRING ", label, " readable with three columns: ", path))
  }
}

ppi_code <- if (species == "human") "9606" else "10090"
ppi_root <- Sys.getenv("TRIAGE_PPI_ROOT", unset = file.path(triage_home, "inputs", "raw", "ppi"))
check_string_file(file.path(ppi_root, paste0(ppi_code, ".protein.aliases.v12.0.txt")), "aliases")
check_string_file(file.path(ppi_root, paste0(ppi_code, ".protein.physical.links.v12.0.txt")), "physical links")

collectri <- file.path(project_root, "inputs", "raw", "collectri",
                       paste0("collectri_", species, "_network.rds"))
if (file.exists(collectri)) {
  collectri_ok <- tryCatch({
    net <- readRDS(collectri)
    is.data.frame(net) && nrow(net) > 0L &&
      all(c("source", "target", "weight") %in% names(net))
  }, error = function(e) FALSE)
  if (collectri_ok) ok(paste0("CollecTRI RDS readable with source/target/weight: ", collectri))
  else bad(paste0("CollecTRI RDS is invalid for decoupleR (requires non-empty source, target, weight columns): ", collectri))
} else {
  bad(paste0("CollecTRI network: ", collectri))
}

cm_xlsx <- file.path(triage_home, "inputs", "raw", "cellmarker",
                     if (species == "human") "Cell_marker_Human.xlsx" else "Cell_marker_Mouse.xlsx")
if (file.exists(cm_xlsx)) {
  cellmarker_ok <- tryCatch({
    sheet <- readxl::read_excel(cm_xlsx, n_max = 1L)
    all(c("cell_name", "marker") %in% names(sheet))
  }, error = function(e) FALSE)
  if (cellmarker_ok) ok(paste0("CellMarkerDB readable with cell_name/marker: ", cm_xlsx))
  else bad(paste0("CellMarkerDB spreadsheet is invalid (requires cell_name and marker columns): ", cm_xlsx))
} else {
  bad(paste0("CellMarkerDB spreadsheet: ", cm_xlsx, " (download per resources/cellmarkerdb/README.md)"))
}

# --- 5. Optional user DEG validation -----------------------------------------
if (!is.null(opt$deg)) {
  REQUIRED_COLS <- c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
  cols <- tryCatch(names(readr::read_delim(opt$deg, delim = NULL, n_max = 1, show_col_types = FALSE)),
                   error = function(e) NULL)
  if (is.null(cols)) bad(paste0("cannot read DEG file: ", opt$deg))
  else {
    m <- setdiff(REQUIRED_COLS, cols)
    if (length(m) > 0) bad(paste0("user DEG missing column(s): ", paste(m, collapse = ", ")))
    else ok("user DEG columns")
  }
}

# --- 6. Endpoint reachability (informational) --------------------------------
if (!opt$`skip-internet` && !grepl("^XXXXX", base_url) && nzchar(base_url)) {
  reach <- tryCatch({
    r <- httr::HEAD(sub("/chat/completions$", "", base_url), httr::timeout(10))
    TRUE
  }, error = function(e) FALSE)
  if (reach) {
    ok("endpoint host reachable")
  } else {
    soft("endpoint host not reachable right now (network down or URL wrong)")
  }
}

# --- Summary ------------------------------------------------------------------
cat("\n=== preflight summary ===\n")
if (length(missing) > 0) {
  cat("MISSING (", length(missing), "):\n", paste0(" - ", missing, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("All required prerequisites present.\n")
if (length(warn) > 0) cat("Warnings: ", paste(warn, collapse = "; "), "\n", sep = "")
quit(status = 0)
