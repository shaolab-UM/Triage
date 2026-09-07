# Installed-package full-workflow runtime.
#
# run_triage() and triage_preflight() let a user who only installed the
# package run the complete DEG-input workflow without a repository
# checkout, terminal, or shell scripts. The stage implementations are NOT
# reimplemented here: run_triage() invokes the validated pipeline scripts
# bundled byte-identically under inst/workflow/ (see inst/workflow/README.md).
# The repository reproducibility/ tree remains the manuscript reproduction
# layer.

.triage_resource_root <- function(resource_root = NULL) {
  root <- resource_root %||% Sys.getenv("TRIAGE_DATA_ROOT", unset = "")
  if (!nzchar(root)) root <- tools::R_user_dir("Triage", "data")
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

.triage_workflow_dir <- function() {
  wd <- system.file("workflow", package = "Triage")
  if (!nzchar(wd) || !dir.exists(file.path(wd, "pipeline"))) {
    stop("run_triage: the bundled workflow runtime is missing from the ",
         "installed package (system.file(\"workflow\", package=\"Triage\") ",
         "does not contain pipeline/). Reinstall the Triage package.")
  }
  wd
}

.triage_bundled_ontology <- function() {
  p <- system.file("extdata", "ontology", "CL-ontology-v2025-07-30.json",
                   package = "Triage")
  if (!nzchar(p) || !file.exists(p)) {
    stop("run_triage: the bundled Cell Ontology snapshot is missing from ",
         "the installed package. Reinstall the Triage package.")
  }
  p
}

# Resource structure validators (mirror reproducibility/scripts/preflight_check.R)
.triage_check_string_file <- function(path, label, missing) {
  if (!file.exists(path)) {
    return(c(missing, paste0("STRING ", label, " file: ", path,
                             " (run setup_triage_resources())")))
  }
  if (file.info(path)$size <= 0L || file.access(path, 4L) != 0L) {
    return(c(missing, paste0("STRING ", label,
                             " file is unreadable or empty: ", path)))
  }
  fields <- tryCatch({
    lines <- readLines(path, n = 20L, warn = FALSE)
    line <- lines[nzchar(trimws(lines))][1]
    strsplit(trimws(line), "\\s+")[[1]]
  }, error = function(e) character(0))
  if (length(fields) != 3L) {
    missing <- c(missing, paste0("STRING ", label,
                                 " file must have three whitespace-delimited columns: ",
                                 path))
  }
  missing
}

.triage_check_collectri <- function(path, missing) {
  if (!file.exists(path)) {
    return(c(missing, paste0("CollecTRI network: ", path,
                             " (run setup_triage_resources())")))
  }
  ok <- tryCatch({
    net <- readRDS(path)
    is.data.frame(net) && nrow(net) > 0L &&
      all(c("source", "target", "weight") %in% names(net))
  }, error = function(e) FALSE)
  if (!ok) {
    missing <- c(missing, paste0("CollecTRI RDS is invalid for decoupleR ",
                                 "(requires non-empty source, target, weight columns): ",
                                 path))
  }
  missing
}

.triage_check_cellmarker <- function(path, missing) {
  if (!file.exists(path)) {
    return(c(missing, paste0("CellMarkerDB spreadsheet: ", path,
                             " (download manually, then register with ",
                             "setup_triage_resources(cellmarker_file = ...))")))
  }
  ok <- tryCatch({
    sheet <- readxl::read_excel(path, n_max = 1L)
    all(c("cell_name", "marker") %in% names(sheet))
  }, error = function(e) FALSE)
  if (!ok) {
    missing <- c(missing, paste0("CellMarkerDB spreadsheet is invalid ",
                                 "(requires cell_name and marker columns): ", path))
  }
  missing
}

# Resolve the effective LLM key and endpoint once, honoring the primary
# variables and the historical fallback variables (LLM_API_KEY_ENV-named
# key variable, CASSIA_API_BASE_URL endpoint). Callers map the result onto
# the canonical DEEPSEEK_API_KEY / LLM_API_BASE_URL names.
.triage_resolve_api <- function(api_key, api_base_url) {
  eff_key <- api_key
  if (is.null(eff_key) || !nzchar(eff_key)) {
    key_env <- Sys.getenv("LLM_API_KEY_ENV", unset = "DEEPSEEK_API_KEY")
    eff_key <- Sys.getenv(key_env, unset = "")
  }
  eff_url <- api_base_url
  if (is.null(eff_url) || !nzchar(eff_url)) {
    eff_url <- Sys.getenv("LLM_API_BASE_URL", unset = "")
    if (!nzchar(eff_url)) {
      eff_url <- Sys.getenv("CASSIA_API_BASE_URL", unset = "")
    }
  }
  list(key = eff_key, url = eff_url)
}

# Resolve the Rscript executable portably (Windows uses Rscript.exe);
# never assume a shell can find it on PATH.
.triage_rscript <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  rscript <- file.path(R.home("bin"), exe)
  if (!file.exists(rscript)) {
    stop("run_triage: Rscript executable not found at ", rscript,
         call. = FALSE)
  }
  rscript
}

#' Check the installed-package full-workflow prerequisites
#'
#' Validates everything `run_triage()` needs in one pass, before any
#' analysis: required R packages (pipeline stages 03a-10), CASSIA
#' availability (including its Python backend via
#' a version-compatible lookup of CASSIA's internal
#' `check_python_env()`), external resources under the Triage
#' user-data directory
#' (STRING, CollecTRI, CellMarkerDB), the bundled Cell Ontology, an API
#' key, and a full chat-completions endpoint. No analysis and no LLM calls
#' are performed, and no repository checkout is assumed.
#'
#' @param species "human" or "mouse".
#' @param api_key API key for the LLM-backed stages. When `NULL`, the
#'   `DEEPSEEK_API_KEY` environment variable (or the variable named by
#'   `LLM_API_KEY_ENV`) is used.
#' @param api_base_url Full chat-completions endpoint URL (ending with
#'   `/chat/completions`). When `NULL`, the `LLM_API_BASE_URL` environment
#'   variable (fallback `CASSIA_API_BASE_URL`) is used.
#' @param resource_root Resource directory; defaults to the Triage
#'   user-data directory populated by `setup_triage_resources()`.
#' @param deg optional path to a user DEG table to column-validate.
#'
#' @return Invisibly, a list with `ok`, `missing` and `notes`. Missing
#'   items are also printed. `run_triage()` calls this automatically and
#'   stops when `ok` is `FALSE`.
#' @export
triage_preflight <- function(species = "human",
                             api_key = NULL,
                             api_base_url = NULL,
                             resource_root = NULL,
                             deg = NULL) {
  species <- match.arg(tolower(species), c("human", "mouse"))
  resource_root <- .triage_resource_root(resource_root)

  missing <- character(0)
  notes <- character(0)
  note <- function(x) { message("  note ", x); notes <<- c(notes, x) }

  message("=== Triage installed-package workflow preflight ===")
  message("species: ", species)
  message("resources: ", resource_root)

  pkg_sets <- list(
    "core (03a-11)" = c("optparse", "jsonlite", "dplyr", "readr", "stringr",
                        "purrr", "tibble", "digest", "writexl", "rlang"),
    "evidence/LLM (05)" = c("httr", "glue", "data.table", "tidyr", "memoise",
                            "cachem", "knitr", "readxl", "tictoc", "future",
                            "future.apply", "xml2", "fs"),
    "in-house reviewer (06/06b)" = c("furrr", "rio", "AnnotationDbi",
                                     "org.Hs.eg.db"),
    "enrichment (05/06b)" = c("clusterProfiler", "DOSE", "ReactomePA",
                              "enrichR", "decoupleR")
  )
  for (nm in names(pkg_sets)) {
    miss <- pkg_sets[[nm]][!vapply(pkg_sets[[nm]], function(p)
      requireNamespace(p, quietly = TRUE), logical(1))]
    if (length(miss) > 0) {
      missing <- c(missing, paste0("R packages [", nm, "]: ",
                                   paste(miss, collapse = ", ")))
    }
  }
  if (!requireNamespace("CASSIA", quietly = TRUE)) {
    missing <- c(missing, "R package: CASSIA (stage 03b; install with install_triage_dependencies())")
  } else {
    # Version-compatible lookup: check_python_env is present in the tested
    # CASSIA revision (b008c0ac); getFromNamespace works whether or not the
    # historical revision exports it.
    py_ok <- tryCatch({
      check_fun <- utils::getFromNamespace("check_python_env", "CASSIA")
      isTRUE(check_fun())
    }, error = function(e) FALSE)
    if (!py_ok) {
      missing <- c(missing,
                   "CASSIA Python backend (run CASSIA::setup_cassia_env() once, then re-run preflight)")
    }
    cassia_sha <- tryCatch(utils::packageDescription("CASSIA")$RemoteSha,
                           error = function(e) NA_character_)
    if (is.null(cassia_sha) || is.na(cassia_sha)) {
      cassia_sha <- tryCatch(utils::packageDescription("CASSIA")$GithubSHA1,
                             error = function(e) NA_character_)
    }
    cassia_tested <- "b008c0ac3dd81b2c2dff131d20f5081a58aca027"
    if (!is.na(cassia_sha) && cassia_sha != cassia_tested) {
      note(paste0("CASSIA revision ", substr(cassia_sha, 1, 12),
                  " differs from the tested revision ",
                  substr(cassia_tested, 1, 12),
                  "; reinstall with install_triage_dependencies() for the canonical configuration"))
    }
  }
  # Canonical KEGG evidence (stage 05 enrichKEGG use_internal_data=TRUE) is
  # backed by KEGG.db. Missing or unloadable KEGG data must fail clearly.
  if (!requireNamespace("KEGG.db", quietly = TRUE)) {
    missing <- c(missing,
                 paste0("R package: KEGG.db (canonical KEGG evidence stage 05; ",
                        "install with install_triage_dependencies() or ",
                        "remotes::install_url(\"https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz\"))"))
  } else if (!tryCatch(length(AnnotationDbi::keys(get("KEGGPATHID2EXTID",
                                                       envir = asNamespace("KEGG.db")))) > 0,
                       error = function(e) FALSE)) {
    missing <- c(missing, "R package: KEGG.db (installed but KEGG data unloads)")
  }
  # Canonical Reactome evidence (stage 05 Reactome local enrichment).
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    missing <- c(missing,
                 "R package: reactome.db (Reactome enrichment stage 05; install with install_triage_dependencies())")
  }
  disgenet_key <- Sys.getenv("DISGENET_API_KEY", unset = "")
  if (nzchar(disgenet_key) && !requireNamespace("disgenet2r", quietly = TRUE)) {
    missing <- c(missing,
                 paste0("R package: disgenet2r (DISGENET_API_KEY is set; ",
                        "install with remotes::install_gitlab(\"medbio/disgenet2r\"))"))
  }
  if (species == "mouse" && !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    missing <- c(missing,
                 "R package: org.Mm.eg.db (required for mouse evidence analysis; install with BiocManager::install(\"org.Mm.eg.db\"))")
  }

  # API credentials
  eff <- .triage_resolve_api(api_key, api_base_url)
  api_key <- eff$key
  api_base_url <- eff$url
  if (!nzchar(api_key) || identical(api_key, "XXXXX")) {
    missing <- c(missing, "API key (pass api_key= or set DEEPSEEK_API_KEY)")
  }
  if (!nzchar(api_base_url) || identical(api_base_url, "XXXXX")) {
    missing <- c(missing, "LLM endpoint (pass api_base_url= as the FULL chat-completions URL, or set LLM_API_BASE_URL)")
  } else if (!grepl("/chat/completions$", api_base_url)) {
    missing <- c(missing, paste0("LLM endpoint must be the FULL chat-completions URL (ends with /chat/completions); got: ",
                                 api_base_url))
  }

  # Bundled Cell Ontology
  cl_json <- tryCatch(.triage_bundled_ontology(), error = function(e) NA_character_)
  if (is.na(cl_json)) {
    missing <- c(missing, "bundled Cell Ontology snapshot (reinstall the Triage package)")
  } else {
    ontology_ok <- tryCatch({
      parsed <- jsonlite::fromJSON(cl_json, simplifyVector = FALSE)
      is.list(parsed) || is.data.frame(parsed)
    }, error = function(e) FALSE)
    if (!ontology_ok) {
      missing <- c(missing, "bundled Cell Ontology JSON is unreadable or invalid")
    }
  }

  # External resources under the Triage user-data directory
  ppi_code <- if (species == "human") "9606" else "10090"
  ppi_root <- file.path(resource_root, "ppi")
  missing <- .triage_check_string_file(
    file.path(ppi_root, paste0(ppi_code, ".protein.aliases.v12.0.txt")),
    "aliases", missing)
  missing <- .triage_check_string_file(
    file.path(ppi_root, paste0(ppi_code, ".protein.physical.links.v12.0.txt")),
    "physical links", missing)
  missing <- .triage_check_collectri(
    file.path(resource_root, "collectri",
              paste0("collectri_", species, "_network.rds")), missing)
  missing <- .triage_check_cellmarker(
    file.path(resource_root, "cellmarker",
              if (species == "human") "Cell_marker_Human.xlsx" else "Cell_marker_Mouse.xlsx"),
    missing)

  # Optional DEG column validation
  if (!is.null(deg) && nzchar(deg)) {
    required_cols <- c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj",
                       "pct.1", "pct.2")
    cols <- tryCatch(
      names(readr::read_delim(deg, delim = NULL, n_max = 1, show_col_types = FALSE)),
      error = function(e) NULL)
    if (is.null(cols)) {
      missing <- c(missing, paste0("cannot read DEG file: ", deg))
    } else {
      m <- setdiff(required_cols, cols)
      if (length(m) > 0) {
        missing <- c(missing, paste0("user DEG missing column(s): ",
                                     paste(m, collapse = ", ")))
      }
    }
  }

  ok <- length(missing) == 0L
  if (!ok) {
    message("MISSING (", length(missing), "):")
    for (m in missing) message(" - ", m)
  } else {
    message("All required prerequisites present.")
    if (length(notes) > 0) message("Notes: ", paste(notes, collapse = "; "))
  }
  invisible(list(ok = ok, missing = missing, notes = notes))
}

#' Run the full Triage workflow on your own DEG input
#'
#' Installed-package entry point for the complete DEG-input workflow: it
#' anonymizes your cluster-level DEG table, runs the upstream
#' reviewer-generation stages, maps reviewer labels with CL-Linker, and
#' runs Triage adjudication through the final post-summary — invoking the
#' validated pipeline scripts bundled with the installed package
#' (`inst/workflow/`); no scientific logic is reimplemented and no
#' repository checkout, shell, or `Rscript` knowledge is required.
#'
#' Required user input is only the DEG table, a species and a tissue (plus
#' API credentials for the LLM-backed stages). You never provide candidate
#' tables, reviewer outputs, CL-Linker mappings, judge-input JSON, Handling
#' Editor output or final adjudication JSON; reference labels are never
#' used for adjudication.
#'
#' External resources are resolved automatically from the Triage user-data
#' directory populated by `setup_triage_resources()`: the bundled Cell
#' Ontology ships with the package, STRING is located via `TRIAGE_PPI_ROOT`,
#' and CollecTRI / CellMarkerDB are staged from
#' `tools::R_user_dir("Triage", "data")` into a per-run runtime home that
#' the validated stage path resolution reads. Users never copy files into
#' `inputs/raw/`.
#'
#' @param deg path to a cluster-level DEG/marker table (csv/tsv) with
#'   columns `cluster`, `gene`, `avg_log2FC`, `p_val`, `p_val_adj`,
#'   `pct.1`, `pct.2`.
#' @param species "human" or "mouse".
#' @param tissue tissue context for the CASSIA stage (e.g. "pancreas").
#' @param study_context optional study context for the CASSIA stage
#'   (default "reference").
#' @param api_key API key for the LLM-backed stages; overrides the
#'   `DEEPSEEK_API_KEY` environment variable. Never printed.
#' @param api_base_url Full chat-completions endpoint URL (ending with
#'   `/chat/completions`); overrides the `LLM_API_BASE_URL` environment
#'   variable.
#' @param out output root directory (default "results").
#' @param dataset_name name for this run (default "user_dataset").
#' @param workers worker count for parallel stages (default 4).
#' @param run_tag optional run tag; defaults to a timestamp.
#' @param reference_labels optional reference-label CSV for
#'   post-adjudication evaluation (stage 11) only; never used for
#'   adjudication.
#' @param preflight_only logical; stage the runtime, run
#'   `triage_preflight()` and return without running any analysis stage.
#' @param resource_root resource directory; defaults to the Triage
#'   user-data directory populated by `setup_triage_resources()`.
#' @param ... reserved for future options; currently unused.
#'
#' @return Invisibly, a list with `out_root`, `runtime_home`, `final_dir`,
#'   `summary_path`, `summary` (parsed stage-10 summary data frame),
#'   `dataset_name`, `species`, `tissue` and `run_tag`. A `preflight_only`
#'   run returns the same list without `final_dir`/`summary_path`/`summary`
#'   stage outputs. All environment-variable changes made for the run are
#'   restored when the function exits (success, preflight-only, or error).
#' @export
run_triage <- function(deg,
                       species,
                       tissue,
                       study_context = "reference",
                       api_key = NULL,
                       api_base_url = NULL,
                       out = "results",
                       dataset_name = "user_dataset",
                       workers = 4,
                       run_tag = NULL,
                       reference_labels = NULL,
                       preflight_only = FALSE,
                       resource_root = NULL,
                       ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    warning("run_triage: unused option(s) ignored: ",
            paste(names(dots), collapse = ", "))
  }
  species <- match.arg(tolower(species), c("human", "mouse"))
  if (is.null(deg) || !nzchar(deg) || !file.exists(deg)) {
    stop("run_triage: deg file not found: ", deg)
  }
  # Canonicalize every path the child stages will receive so that relative
  # arguments (and the default out = "results") resolve identically from any
  # working directory, and stage-to-stage handoffs stay absolute.
  deg <- normalizePath(deg, mustWork = TRUE)
  if (is.null(tissue) || !nzchar(tissue)) {
    stop("run_triage: tissue is required (e.g. tissue = \"pancreas\").")
  }

  # Canonicalize API configuration once. run_triage() accepts the primary
  # variables directly (api_key / api_base_url) and falls back to the
  # historical fallback variables (LLM_API_KEY_ENV-named key variable and
  # CASSIA_API_BASE_URL), but the child stages only consume the canonical
  # DEEPSEEK_API_KEY / LLM_API_BASE_URL names. The effective values are
  # resolved here once and mapped onto the canonical names for the duration
  # of the run; original environment values are restored on exit (success,
  # preflight_only, or stage error). The key is never printed.
  eff <- .triage_resolve_api(api_key, api_base_url)
  api_key <- eff$key
  api_base_url <- eff$url
  withr::local_envvar(DEEPSEEK_API_KEY = api_key,
                      LLM_API_BASE_URL = api_base_url)

  resource_root <- .triage_resource_root(resource_root)
  workflow_dir <- .triage_workflow_dir()
  pipeline_dir <- file.path(workflow_dir, "pipeline")
  ontology_src <- .triage_bundled_ontology()

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
  }
  run_tag <- run_tag %||% format(Sys.time(), "%Y%m%d_%H%M%S")
  out_root <- file.path(out, dataset_name, run_tag)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  # The run directory now exists, so it canonicalizes to an absolute path.
  out_root <- normalizePath(out_root, mustWork = TRUE)
  runtime_home <- file.path(out_root, "runtime_home")

  message("run_triage: dataset = ", dataset_name,
          " | species = ", species, " | tissue = ", tissue,
          " | study context = ", study_context)
  message("run_triage: DEG source = ", deg)
  message("run_triage: output root = ", out_root)

  # --- stage external resources into a per-run runtime home --------------
  # The validated stage path resolution reads TRIAGE_HOME / PROJECT_ROOT
  # defaults; run_triage() stages the user-data resources (and the bundled
  # ontology) so those defaults resolve without any manual copying.
  onto_dir <- file.path(runtime_home, "inputs", "raw", "ontology")
  coll_dir <- file.path(runtime_home, "inputs", "raw", "collectri")
  cm_dir <- file.path(runtime_home, "inputs", "raw", "cellmarker")
  for (d in c(onto_dir, coll_dir, cm_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  cl_json <- file.path(onto_dir, basename(ontology_src))
  file.copy(ontology_src, cl_json, overwrite = TRUE)
  collectri_src <- file.path(resource_root, "collectri",
                             paste0("collectri_", species, "_network.rds"))
  file.copy(collectri_src, coll_dir, overwrite = TRUE)
  cm_name <- if (species == "human") "Cell_marker_Human.xlsx" else
    "Cell_marker_Mouse.xlsx"
  cm_src <- file.path(resource_root, "cellmarker", cm_name)
  cm_dest <- file.path(cm_dir, cm_name)
  file.copy(cm_src, cm_dest, overwrite = TRUE)

  # --- per-run generic dataset context ------------------------------------
  # Generic runs must NOT inherit any manuscript dataset profile merely
  # because the dataset name collides with a benchmark alias. A neutral
  # context file is written into the run's runtime home and wired through
  # DATASET_CONTEXT_JSON_PATH (supported by get_dataset_config as an
  # explicit override applied on top of every other profile source).
  ctx_json <- file.path(runtime_home, "dataset_context.json")
  jsonlite::write_json(list(
    species = species,
    tissue = tissue,
    study_context = study_context,
    dataset_scope = "mixed_unknown",
    scope_profile = "mixed_unknown",
    gate_mode = "flag_only",
    allowed_lineages = character(0),
    user_notes = paste(
      "User-supplied dataset; generic workflow context.",
      "Goal: marker-based candidate cell type annotation at cluster/state level.",
      "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
    )
  ), ctx_json, auto_unbox = TRUE, pretty = TRUE)

  # Stage-specific environment variables are also restored on exit (withr
  # registers an exit handler on this frame); child Rscript processes
  # inherit the staged values while the function runs.
  withr::local_envvar(
    TRIAGE_HOME = runtime_home,
    PROJECT_ROOT = runtime_home,
    CL_LOCAL_JSON = cl_json,
    TRIAGE_PPI_ROOT = file.path(resource_root, "ppi"),
    TRIAGE_WORKFLOW_DIR = workflow_dir,
    DATASET_CONTEXT_JSON_PATH = ctx_json
  )

  # --- preflight (installed-package runtime) ------------------------------
  message(">>> [preflight] full-workflow dependency check")
  pf <- triage_preflight(species = species, resource_root = resource_root,
                         deg = deg)
  if (!isTRUE(pf$ok)) {
    stop("run_triage: preflight failed; resolve the reported prerequisites ",
         "and rerun.")
  }
  if (isTRUE(preflight_only)) {
    message("run_triage: preflight_only requested; stopping before analysis.")
    return(invisible(list(out_root = out_root, runtime_home = runtime_home,
                          dataset_name = dataset_name, species = species,
                          tissue = tissue, run_tag = run_tag,
                          preflight = pf)))
  }

  fail <- function(msg) stop("run_triage: ", msg, call. = FALSE)
  run_stage <- function(script, args, stage_label) {
    message(">>> [", stage_label, "] Rscript ", script)
    status <- system2(.triage_rscript(),
                      shQuote(c(file.path(pipeline_dir, script), args)),
                      stdout = "", stderr = "")
    if (!identical(as.integer(status), 0L)) {
      fail(paste0("stage ", stage_label, " (", script,
                  ") failed with status ", status))
    }
  }

  # --- input preparation ---------------------------------------------------
  map_dir <- file.path(out_root, "01_input_prep")
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  prep_args <- c("--deg", deg, "--out-dir", map_dir)
  run_stage("01b_prepare_user_deg.R", prep_args, "01b user-DEG preparation")
  maskdeg <- file.path(map_dir, "maskdeg.csv")
  deg_for_stages <- maskdeg
  run_dir <- out_root
  # ONE canonical per-run evidence directory. Stage 05 writes its enrichment
  # TSVs relative to its working directory (= run_dir) under
  # intermediate_outputs/<dataset>_LLM_Input_Run/bioinformatics_tsv; both
  # downstream consumers (06b and 08) are pointed at exactly that directory.
  intermediate_run_dir <- file.path(run_dir, "intermediate_outputs",
                                   paste0(dataset_name, "_LLM_Input_Run"))

  # 03a: deterministic DEG filtering
  deg_filtered <- file.path(run_dir, "filtered_deg.csv")
  run_stage("03a_filter_deg.R", c("--deg", deg_for_stages,
                                  "--out_dir", run_dir), "03a")

  # 03b: CASSIA reviewer (LLM API; CASSIA strips the /chat/completions suffix)
  cassia_out <- file.path(run_dir, "03b_cassia")
  run_stage("03b_run_cassia.R", c(
    "--deg", deg_filtered, "--out_dir", cassia_out,
    "--tissue", tissue, "--species", species,
    "--study_context", study_context,
    "--workers", as.character(workers),
    "--api_base_url", Sys.getenv("LLM_API_BASE_URL")
  ), "03b")

  # 04a: structured candidates
  candidates_out <- file.path(run_dir, "04a_candidates")
  run_stage("04a_build_candidates.R", c(
    "--deg", deg_filtered, "--out_dir", candidates_out,
    "--species", species, "--ontology", cl_json
  ), "04a")

  candidates_csv <- file.path(candidates_out, "cellanno", "structured_candidates.csv")

  # 05: evidence / dossier construction (writes relative to cwd = run_dir)
  withr::with_dir(run_dir, {
    run_stage("05_prepare_llm_inputs.R", c(
      "--mode", "full",
      "--deg_file", deg_filtered,
      "--candidates_file", candidates_csv,
      "--out_root", file.path(run_dir, "05c_llm_queries"),
      "--dataset_name", dataset_name,
      "--project_root", runtime_home
    ), "05")
  })

  # 06: in-house reviewer (LLM API)
  run_stage("06_run_llm_pipeline.R", c(
    "--mode", "full",
    "--input_root", file.path(run_dir, "05c_llm_queries"),
    "--output_root", file.path(run_dir, "06_llm_outputs"),
    "--dataset_name", dataset_name,
    "--project_root", runtime_home,
    "--workers", as.character(workers),
    "--max_rounds", "3",
    "--cl_local_json", cl_json
  ), "06")

  # 06b: enrichment reviewer (clusterProfiler + CellMarkerDB)
  run_stage("06b_run_inter.R", c(
    "--marker_csv", deg_filtered,
    "--step1_dir", file.path(run_dir, "05c_llm_queries", "step1_report_queries"),
    "--bioinfo_dir", file.path(intermediate_run_dir, "bioinformatics_tsv"),
    "--out_dir", file.path(run_dir, "06b_inter"),
    "--dataset_name", dataset_name,
    "--species", species,
    "--workers", as.character(workers),
    "--cl_local_json", cl_json,
    "--cellmarker_human_path", cm_dest,
    "--cellmarker_mouse_path", cm_dest
  ), "06b")

  # 07 / 07b: reviewer summaries
  run_stage("07_our_llm_summary.R", c(
    "--out_root", file.path(run_dir, "06_llm_outputs"),
    "--dataset_name", dataset_name,
    "--out_dir", file.path(run_dir, "07_our_summary")
  ), "07")

  run_stage("07b_inter_summary.R", c(
    "--in_dir", file.path(run_dir, "06b_inter", "final_passed"),
    "--out_dir", file.path(run_dir, "07b_inter_summary"),
    "--dataset_name", dataset_name
  ), "07b")

  # 07.5: CL-Linker reviewer mapping (LLM-backed)
  manifest <- file.path(run_dir, "evidence_mapping_run.tsv")
  writeLines(c("dataset\trun_dir",
               paste0(dataset_name, "\t", run_dir)),
             manifest, sep = "\n")
  run_stage("07.5_build_reviewer_mapping.R", c(
    "--cl_json", cl_json,
    "--out_dir", file.path(run_dir, "07.5_mapping"),
    "--manifest", manifest,
    "--dataset", dataset_name
  ), "07.5")

  # 08: adjudication input assembly
  registry_csv <- list.files(file.path(run_dir, "07.5_mapping"),
                             pattern = "reviewer_mapping_registry.*[.]tsv$",
                             full.names = TRUE)
  if (length(registry_csv) != 1L) {
    fail("expected exactly one reviewer mapping registry under 07.5_mapping")
  }
  run_stage("08_build_judge_inputs.R", c(
    "--cassia_csv", file.path(cassia_out,
      list.files(cassia_out, pattern = "annotation_cassia_FINAL_RESULTS.csv$",
                 recursive = TRUE)[1]),
    "--in_house_summary_csv", file.path(run_dir, "07_our_summary", "summary.csv"),
    "--enrichment_summary_csv", file.path(run_dir, "07b_inter_summary", "summary.csv"),
    "--mapping_registry_csv", registry_csv,
    "--intermediate_outputs_dir", intermediate_run_dir,
    "--step1_dir", file.path(run_dir, "05c_llm_queries", "step1_report_queries"),
    "--out_dir", file.path(run_dir, "08_judge_inputs"),
    "--dataset_name", dataset_name
  ), "08")

  # 09: adjudication (Handling Editor; Chief QC / release-state in batch workflow)
  run_stage("09_run_judge.R", c(
    "--judge_input_dir", file.path(run_dir, "08_judge_inputs"),
    "--out_root", file.path(run_dir, "09_judge_outputs"),
    "--dataset_name", dataset_name,
    "--temperature", "0",
    "--workers", as.character(workers),
    "--max_rounds", "3",
    "--always_run_chief",
    "--release_policy", "auto",
    "--ols_first", "FALSE",
    "--cl_local_json", cl_json
  ), "09")

  # 10: final summary
  run_stage("10_judge_post_summary.R", c(
    "--out_root", file.path(run_dir, "09_judge_outputs"),
    "--dataset_name", dataset_name
  ), "10")

  # 11: optional evaluation (reference labels enter ONLY here)
  if (!is.null(reference_labels) && nzchar(reference_labels)) {
    reference_labels <- normalizePath(reference_labels, mustWork = TRUE)
    run_stage("11_eval_accuracy.R", c(
      "--dataset_name", dataset_name,
      "--true_label_csv", reference_labels,
      "--cassia_csv", file.path(cassia_out,
        list.files(cassia_out, pattern = "annotation_cassia_FINAL_RESULTS.csv$",
                   recursive = TRUE)[1]),
      "--our_csv", file.path(run_dir, "07_our_summary", "summary.csv"),
      "--inter_csv", file.path(run_dir, "07b_inter_summary", "summary.csv"),
      "--judge_csv", file.path(run_dir, "09_judge_outputs", "summary_final.csv"),
      "--judge_final_dir", file.path(run_dir, "09_judge_outputs", "final"),
      "--out_dir", file.path(run_dir, "11_eval"),
      "--cl_json", cl_json
    ), "11 evaluation")
  } else {
    message("run_triage: no reference labels supplied; evaluation (stage 11) skipped. ",
            "Reference labels are never used for adjudication.")
  }

  message("run_triage: complete. Outputs under ", out_root)
  final_dir <- file.path(run_dir, "09_judge_outputs", "final")
  summary_path <- file.path(run_dir, "09_judge_outputs",
                            "judge_post_summary", "summary.csv")
  summary_df <- if (file.exists(summary_path)) {
    tryCatch(readr::read_csv(summary_path, show_col_types = FALSE),
             error = function(e) NULL)
  } else NULL
  invisible(list(out_root = out_root, runtime_home = runtime_home,
                 final_dir = final_dir, summary_path = summary_path,
                 summary = summary_df,
                 dataset_name = dataset_name, species = species,
                 tissue = tissue, run_tag = run_tag))
}
