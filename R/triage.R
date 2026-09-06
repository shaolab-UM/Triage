# =========================================================================
# triage.R — public API wrappers
# Thin, validated entry points over the internal deterministic core.
# Loading this package never triggers network access or API calls;
# API-backed paths require credentials only when explicitly invoked.
# =========================================================================

#' Default Cell Ontology JSON path (package convenience)
#'
#' Resolves the bundled Cell Ontology JSON (v2025-07-30) shipped with the
#' package. The environment variable `CL_LOCAL_JSON` overrides it.
#'
#' @return Character path to a Cell Ontology JSON file.
#' @noRd
default_ontology_path <- function() {
  env <- Sys.getenv("CL_LOCAL_JSON", unset = "")
  if (nzchar(env)) return(env)
  p <- system.file("extdata", "ontology", "CL-ontology-v2025-07-30.json",
                   package = "Triage")
  if (nzchar(p)) return(p)
  stop("Cell Ontology JSON not found. Set CL_LOCAL_JSON to a valid ontology JSON path.")
}

#' Load the Cell Ontology handle (package convenience)
#'
#' Builds the configuration and graph objects used by ontology-aware
#' adjudication and CL-Linker mapping. No network access: only the local
#' ontology JSON is used.
#'
#' @param ontology_path Path to a Cell Ontology JSON file. Defaults to
#'   \code{default_ontology_path()}.
#' @param cache_dir Directory for ontology graph caches.
#'
#' @return A list with elements `cfg` (ontology configuration) and `graph`
#'   (ontology graph).
#' @export
load_triage_ontology <- function(ontology_path = default_ontology_path(),
                                 cache_dir = ".ols_cache") {
  cfg <- make_cl_cfg(ontology_path, prefer_ols = FALSE, cache_dir = cache_dir)
  graph <- load_cl_graph(cfg)
  list(cfg = cfg, graph = graph)
}

#' Read a Triage adjudication input file (package convenience)
#'
#' Reads (or validates, if already parsed) a Triage reviewer/adjudication
#' input object: a named list containing at least `cluster_id` and an
#' `inputs` block with reviewer summaries. This is the JSON schema produced
#' by the manuscript input-construction pipeline (stage 08) and expected by
#' the deterministic adjudication core. The pipeline code that constructs
#' these inputs (reviewer parsing, dossier assembly) is in
#' `reproducibility/scripts/pipeline/08_build_judge_inputs.R`.
#'
#' @param x Path to a JSON file, a JSON string, or an already-parsed list.
#'
#' @return The parsed input as a list.
#' @export
read_triage_input <- function(x) {
  obj <- NULL
  if (is.character(x) && length(x) == 1) {
    if (file.exists(x)) {
      obj <- jsonlite::fromJSON(x, simplifyVector = FALSE)
    } else {
      obj <- tryCatch(jsonlite::fromJSON(x, simplifyVector = FALSE),
                      error = function(e) NULL)
      if (is.null(obj)) {
        stop("read_triage_input: '", x, "' is neither an existing file nor valid JSON text.")
      }
    }
  } else if (is.list(x)) {
    obj <- x
  } else {
    stop("read_triage_input: expects a JSON path, JSON string, or parsed list.")
  }
  if (is.null(obj$cluster_id) || !nzchar(as.character(obj$cluster_id))) {
    stop("read_triage_input: input is missing required field 'cluster_id'.")
  }
  if (is.null(obj$inputs) || !is.list(obj$inputs)) {
    stop("read_triage_input: input is missing required field 'inputs' (reviewer summaries).")
  }
  obj
}

#' Map a cell-type label onto the Cell Ontology (package convenience)
#'
#' Exposes the internal label-matching used by the adjudication core.
#' This is a matching utility, not an annotation step: candidate annotations
#' are generated upstream of Triage.
#'
#' @param label Free-text cell-type label.
#' @param provided_clid Optional reviewer-proposed CL identifier.
#' @param ontology Ontology handle from [load_triage_ontology()].
#' @param min_best,min_delta Matching thresholds (see internal `coerce_clid`).
#'
#' @return The mapped value (a `CL:` identifier) or NA when unmapped.
#' @noRd
map_cell_ontology <- function(label,
                              provided_clid = NULL,
                              ontology = NULL,
                              min_best = 0.6,
                              min_delta = 0.1) {
  if (is.null(ontology)) ontology <- load_triage_ontology()
  coerce_clid(label, provided_clid, ontology$cfg, ontology$graph,
              min_best = min_best, min_delta = min_delta)
}

#' Run CL-Linker cell-type label mapping
#'
#' Published lexical CL-Linker mapping (label + optional provided CL ID to
#' a Cell Ontology identifier). Fully deterministic and offline.
#'
#' @param labels Character vector of cell-type labels.
#' @param provided_clids Optional vector/list of provided CL identifiers.
#' @param ontology_path Path to a Cell Ontology JSON file.
#'
#' @return A data frame of mapping results (predicted CL ID, scores, status).
#' @export
run_cl_linker <- function(labels,
                          provided_clids = NULL,
                          ontology_path = default_ontology_path()) {
  idx <- build_cl_index(ontology_path)
  cl_link_batch(labels, provided_clids = provided_clids, idx = idx)
}

#' Build a normalized reviewer summary block (package convenience)
#'
#' Validates and normalizes reviewer records (CASSIA, in-house,
#' clusterProfiler/enrichment) into the summary schema consumed by the
#' adjudication core: `top1_cell_type`, `topk_cell_types`,
#' `cell_ontology_id`, `mapping`, `reasoning_short`, plus raw inputs.
#'
#' @param reviewers Named list of reviewer records, e.g.
#'   `list(cassia = list(...), in_house = list(...), enrich = list(...))`.
#'   Each record requires `top1_cell_type`; all other fields are optional.
#' @param cl Only for internal validation messages.
#'
#' @return A named list with normalized reviewer summaries.
#' @noRd
build_reviewer_summary <- function(reviewers, cl = NULL) {
  if (!is.list(reviewers) || length(reviewers) == 0) {
    stop("build_reviewer_summary: 'reviewers' must be a non-empty named list.")
  }
  required <- "top1_cell_type"
  out <- list()
  for (nm in names(reviewers)) {
    rec <- reviewers[[nm]]
    if (!is.list(rec)) {
      stop("build_reviewer_summary: reviewer '", nm, "' must be a list.")
    }
    missing <- setdiff(required, names(rec))
    if (length(missing) > 0) {
      stop("build_reviewer_summary: reviewer '", nm, "' missing required field(s): ",
           paste(missing, collapse = ", "))
    }
    rec$topk_cell_types <- rec$topk_cell_types %||% list(rec$top1_cell_type)
    rec$cell_ontology_id <- rec$cell_ontology_id %||% NA_character_
    rec$mapping <- rec$mapping %||% list()
    rec$reasoning_short <- rec$reasoning_short %||% ""
    out[[nm]] <- rec
  }
  out
}

#' Build a Triage adjudication input object
#'
#' Assembles the adjudication (judge) input consumed by
#' [run_triage_adjudication()] from reviewer summaries and optional DEG /
#' dossier evidence. Pure data assembly; no API access. This is a simplified
#' constructor: the full stage-08 dossier assembly lives in the pipeline
#' script `reproducibility/scripts/pipeline/08_build_judge_inputs.R`.
#'
#' @param cluster_id Cluster identifier string.
#' @param reviewers Named list of reviewer records (see
#'   \code{build_reviewer_summary()}).
#' @param deg Optional DEG data frame (read with the internal \code{read_deg}
#'   reader, columns: cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2).
#' @param dossier Optional list of dossier evidence (marker genes,
#'   enrichment terms, PMIDs).
#'
#' @return A list shaped like the stage-08 judge input JSON.
#' @export
build_adjudication_input <- function(cluster_id, reviewers,
                                     deg = NULL, dossier = NULL) {
  if (missing(cluster_id) || !nzchar(as.character(cluster_id))) {
    stop("build_adjudication_input: 'cluster_id' is required.")
  }
  summaries <- build_reviewer_summary(reviewers)
  inputs <- list()
  alias <- list(
    cassia = "cassia_summary",
    in_house = "our_summary",
    enrich = "enrich_summary",
    enrichment = "enrich_summary"
  )
  reviewer_reports <- list()
  for (nm in names(summaries)) {
    key <- alias[[nm]] %||% paste0(nm, "_summary")
    inputs[[key]] <- summaries[[nm]]
    rr <- summaries[[nm]]
    rr$reviewer_role <- nm
    reviewer_reports[[length(reviewer_reports) + 1L]] <- rr
  }
  inputs$reviewer_reports <- reviewer_reports
  inputs$dossier <- dossier %||% list()
  jin <- list(cluster_id = as.character(cluster_id), inputs = inputs)
  if (!is.null(deg) && is.data.frame(deg)) {
    top <- utils::head(deg[order(-abs(deg$avg_log2FC %||% rep(0, nrow(deg)))), , drop = FALSE], 50)
    jin$inputs$deg_top <- as.list(top$gene %||% list())
  }
  jin
}

#' Run Triage adjudication for one cluster
#'
#' Runs the deterministic adjudication core for one cluster. Two modes:
#' \itemize{
#'   \item `use_api = TRUE` invokes the Handling Editor stage and then runs
#'     the deterministic downstream processing (release policy, CL
#'     normalization, citation gates, local validation gate). Requires
#'     `DEEPSEEK_API_KEY` (or an explicit `api_key`) and `LLM_API_BASE_URL`;
#'     fails with a clear message when credentials are missing.
#'   \item `use_api = FALSE` requires a precomputed Handling Editor output
#'     (`head_output`: a path to JSON, a JSON string, or a parsed list) and
#'     executes only the deterministic downstream processing. No network
#'     access.
#' }
#'
#' Full manuscript workflow orchestration is available under
#' `reproducibility/scripts/`.
#
#' @param input Triage adjudication input (from [read_triage_input()] or
#'   [build_adjudication_input()]).
#' @param ontology Ontology handle from [load_triage_ontology()]; loaded
#'   automatically when NULL.
#' @param use_api Logical; run the API-backed path?
#' @param head_output Precomputed head-editor adjudication (path, JSON string
#'   or list); required for the deterministic path.
#' @param model Model identifier for the API path (primary pipeline default:
#'   `deepseek-v4-flash`).
#' @param api_key API key for the API path; defaults to the environment
#'   variable named by `TRIAGE_LLM_API_KEY_ENV` (default `DEEPSEEK_API_KEY`).
#' @param api_base_url Full chat-completions endpoint URL for the API path;
#'   defaults to the `LLM_API_BASE_URL` environment variable (with the
#'   documented fallback).
#' @param temperature Sampling temperature (primary pipeline: 0).
#' @param prompt_profile Head-editor prompt profile: "compact" (primary
#'   pipeline default) or "default".
#' @param species_value Species label for marker context.
#' @param dataset_name Dataset alias for the release dataset configuration
#'   (e.g. "Census_immune", "screview_clean", "validation_test").
#' @param dataset_cfg Optional prebuilt dataset configuration list; overrides
#'   `dataset_name`.
#' @param project_root Root used to resolve dataset configuration files.
#'
#' @return The final adjudication record (list), with attribute
#'   `gate` = local gate result.
#' @export
run_triage_adjudication <- function(input,
                                    ontology = NULL,
                                    use_api = FALSE,
                                    head_output = NULL,
                                    model = "deepseek-v4-flash",
                                    api_key = NULL,
                                    api_base_url = NULL,
                                    temperature = 0,
                                    prompt_profile = "compact",
                                    species_value = "human",
                                    dataset_name = NULL,
                                    dataset_cfg = NULL,
                                    project_root = getwd()) {
  jin <- if (inherits(input, "list")) input else read_triage_input(input)
  if (is.null(ontology)) ontology <- load_triage_ontology()
  if (is.null(dataset_cfg) && !is.null(dataset_name)) {
    dataset_cfg <- get_dataset_config(dataset_name, project_root = project_root)
  }
  if (is.null(dataset_cfg)) {
    dataset_cfg <- get_dataset_config("unknown", project_root = project_root)
    dataset_cfg$dataset_name <- as.character(jin$cluster_id)
  }

  if (isTRUE(use_api)) {
    sys_prompt <- if (identical(prompt_profile, "compact")) {
      build_head_editor_system_prompt_compact(species_value = species_value)
    } else {
      build_head_editor_system_prompt(species_value = species_value)
    }
    query <- build_head_editor_query(jin, ontology$graph, ontology$cfg)
    key_env <- Sys.getenv("TRIAGE_LLM_API_KEY_ENV", unset = "DEEPSEEK_API_KEY")
    key <- api_key %||% Sys.getenv(key_env, unset = "")
    if (!nzchar(key) || identical(key, "XXXXX")) {
      stop("run_triage_adjudication(use_api = TRUE): no API key found. Set ", key_env,
           " (or pass api_key =). The deterministic path (use_api = FALSE) needs no key.")
    }
    # Publication-pipeline pattern: serialize ONLY the Handling Editor query
    # object; the wrapper builds the chat-completions request itself and
    # receives the system prompt separately.
    query_str <- jsonlite::toJSON(query, auto_unbox = TRUE, null = "null")
    raw <- invoke_deepseek_api(query_str, api_key = key, model = model,
                               temperature = temperature, system_prompt = sys_prompt,
                               base_url = api_base_url)
    if (!isTRUE(raw$ok)) {
      stop("run_triage_adjudication: Handling Editor API request failed (status: ",
           raw$status, "; error: ", raw$error, ")")
    }
    head_txt <- extract_first_json_object_stack(raw$text,
                                                expected_cluster_id = as.character(jin$cluster_id))
    if (is.null(head_txt) || (is.character(head_txt) &&
        (length(head_txt) != 1 || is.na(head_txt) || !nzchar(head_txt)))) {
      stop("run_triage_adjudication: could not parse a JSON adjudication object from the model response.")
    }
    head_out <- if (is.list(head_txt)) head_txt else
      jsonlite::fromJSON(head_txt, simplifyVector = FALSE)
  } else {
    if (is.null(head_output)) {
      stop("run_triage_adjudication(use_api = FALSE) requires 'head_output': ",
           "a path to a precomputed head-editor JSON, a JSON string, or a parsed list.")
    }
    head_out <- if (is.list(head_output)) head_output
      else if (is.character(head_output) && length(head_output) == 1 && file.exists(head_output)) {
        read_json_safely(head_output)
      } else {
        extract_first_json_object_stack(paste(head_output, collapse = "\n"))
      }
    if (is.null(head_out) || !is.list(head_out)) {
      stop("run_triage_adjudication: could not read 'head_output' as an adjudication object.")
    }
  }

  cite_req <- extract_citation_requirements(jin)
  out <- postprocess_head_out(
    head_out, jin, cite_req, ontology$cfg, ontology$graph,
    dataset_cfg = dataset_cfg, species_value = species_value, release_policy = "auto"
  )
  out <- normalize_manual_review_plan(out)
  out <- force_array_fields(out)
  out <- ensure_enrich_verdict(out, jin)
  out <- normalize_judge_final_decision_cl_preserve(out, ontology$cfg)
  gate <- judge_local_gate(out, cite_req, jin)
  attr(out, "gate") <- gate
  out
}

#' Validate a Triage adjudication result (package convenience)
#'
#' Deterministic validation of a final adjudication record against the
#' release schema, citation policy and (optionally) the adjudication input.
#' Package convenience wrapper around the deterministic local gate used by
#' the release pipeline.
#'
#' @param result Adjudication record (list) or path/JSON string.
#' @param input Optional matching adjudication input for citation checks.
#'
#' @return Invisibly, the parsed result. Raises an error listing all gate
#'   violations when the record is invalid.
#' @export
validate_triage_result <- function(result, input = NULL) {
  j <- if (is.list(result)) result else {
    if (is.character(result) && length(result) == 1 && file.exists(result)) {
      read_json_safely(result)
    } else {
      extract_first_json_object_stack(paste(result, collapse = "\n"))
    }
  }
  if (is.null(j) || !is.list(j)) {
    stop("validate_triage_result: could not read the adjudication result.")
  }
  cite_req <- if (is.null(input)) NULL else extract_citation_requirements(read_triage_input(input))
  gate <- judge_local_gate(j, cite_req, if (is.null(input)) NULL else read_triage_input(input))
  if (!isTRUE(gate$ok)) {
    stop("validate_triage_result: adjudication record failed the local gate: ",
         paste(gate$reasons, collapse = "; "))
  }
  invisible(j)
}

#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom utils data head
#' @importFrom stringr str_trim str_squish str_detect str_replace str_replace_all str_to_lower str_extract
#' @importFrom rlang .data
NULL
