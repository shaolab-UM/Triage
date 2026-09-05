norm_name <- function(x) {
  x %||% "" |>
    stringr::str_trim() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_replace("\\s*\\(.*?\\)", "") |>
    stringr::str_trim()
}

canonicalize_name <- function(name) {
  n <- norm_name(name)
  alias <- list(
    "excitatory neuron" = "glutamatergic neuron",
    "cortical excitatory neuron" = "glutamatergic neuron",
    "glutamatergic neuron" = "glutamatergic neuron",
    "inhibitory neuron" = "GABAergic neuron",
    "gaba neuron" = "GABAergic neuron",
    "neuronal" = "neuron",
    "neuronal progenitor" = "neural progenitor cell",
    "neuronal progenitor cell" = "neural progenitor cell",
    "intermediate progenitor" = "neural progenitor cell",
    "intermediate progenitor cell" = "neural progenitor cell",
    "radial glia" = "radial glial cell",
    "radial glia-1" = "radial glial cell",
    "radial glia-2" = "radial glial cell",
    "cycling progenitor" = "neural progenitor cell",
    "choroid" = "choroid plexus epithelial cell",
    "mesenchyme" = "mesenchymal cell",
    "retinal pigment epithelium" = "retinal pigment epithelial cell",
    "retinal progenitor" = "retinal progenitor cell",
    "retinal ganglion cell" = "retinal ganglion cell",
    "ribosomal protein" = "cell",
    "ribosomal protein (rp)" = "cell",
    # Blood cell plurals -> singulars (Cell Ontology uses singular forms)
    "monocytes" = "monocyte",
    "classical monocytes" = "classical monocyte",
    "macrophages" = "macrophage",
    "neutrophils" = "neutrophil",
    "lymphocytes" = "lymphocyte",
    "t cells" = "T cell",
    "b cells" = "B cell",
    "nk cells" = "natural killer cell",
    "natural killer" = "natural killer cell",
    "nk cell" = "natural killer cell",
    "mesenchymal stromal cell" = "mesenchymal stem cell of the bone marrow",
    "msc" = "mesenchymal stem cell of the bone marrow",
    "mesenchymal stem cell" = "mesenchymal stem cell of the bone marrow",
    "basal epithelial cell" = "basal cell",
    "basal epithelial cells" = "basal cell",
    "naive t cell" = "naive T cell",
    "naive t cells" = "naive T cell",
    "plasmacytoid dendritic cell" = "plasmacytoid dendritic cell",
    "plasmacytoid dendritic cells" = "plasmacytoid dendritic cell",
    "pdc" = "plasmacytoid dendritic cell",
    "pdc cell" = "plasmacytoid dendritic cell",
    "pdcs" = "plasmacytoid dendritic cell",
    "dendritic cells" = "dendritic cell",
    "erythrocytes" = "erythrocyte",
    "platelets" = "platelet",
    "basophils" = "basophil",
    "eosinophils" = "eosinophil",
    "granulocytes" = "granulocyte",
    "mast cells" = "mast cell",
    "plasma cells" = "plasma cell",
    "stem cells" = "stem cell",
    "progenitor cells" = "progenitor cell",
    "fibroblasts" = "fibroblast",
    "endothelial cells" = "endothelial cell",
    "epithelial cells" = "epithelial cell",
    "neurons" = "neuron",
    "astrocytes" = "astrocyte",
    "microglia" = "microglial cell",
    "oligodendrocytes" = "oligodendrocyte"
  )
  if (n %in% names(alias)) return(alias[[n]])
  # Generic plural stripping for cell types ending in 's' (but not already ending in 'cell' or 'cells')
  if (grepl("s$", n) && !grepl("cell$|cells$", n)) {
    singular <- sub("s$", "", n)
    # Return singular form for lookup
    return(singular)
  }
  # Handle "cells" -> "cell" suffix
  if (grepl(" cells$", n)) {
    return(sub(" cells$", " cell", n))
  }
  name
}

LOCAL_INDEX_CACHE <- new.env(parent = emptyenv())

make_cl_cfg <- function(local_json_path, prefer_ols = TRUE, cache_dir = NULL) {
  if (is.null(cache_dir) || !nzchar(cache_dir)) cache_dir <- ".ols_cache"
  list(
    local_json_path = local_json_path,
    prefer_ols = isTRUE(prefer_ols),
    cache_dir = cache_dir,
    local_index = NULL
  )
}

build_local_index <- function(local_json_path) {
  if (is.null(local_json_path) || !file.exists(local_json_path)) return(NULL)
  cl <- jsonlite::fromJSON(local_json_path, simplifyVector = FALSE)
  label_to <- new.env(parent = emptyenv())
  clid_to <- new.env(parent = emptyenv())

  add_label <- function(lbl, clid) {
    if (is.null(lbl)) return()
    lbl <- as.character(lbl)[1]
    if (is.na(lbl) || !nzchar(lbl)) return()
    n <- norm_name(lbl)
    if (!nzchar(n)) return()
    if (!exists(n, envir = label_to, inherits = FALSE)) {
      assign(n, list(clid = clid, label = lbl), envir = label_to)
    }
  }

  for (id in names(cl)) {
    if (!startsWith(id, "CL:")) next
    term <- cl[[id]]
    if (isTRUE(term$deprecated)) next
    assign(id, term$label %||% NA_character_, envir = clid_to)
    add_label(term$label %||% NULL, id)
    syn <- term$synonyms %||% NULL
    if (!is.null(syn)) {
      syns <- if (is.list(syn)) unlist(syn, recursive = TRUE, use.names = FALSE) else syn
      for (s in syns) add_label(s, id)
    }
  }

  list(label_to = label_to, clid_to = clid_to)
}

get_local_index <- function(cl_cfg) {
  key <- as.character(cl_cfg$local_json_path %||% "")
  if (!nzchar(key)) return(NULL)
  if (exists(key, envir = LOCAL_INDEX_CACHE, inherits = FALSE)) {
    return(get(key, envir = LOCAL_INDEX_CACHE, inherits = FALSE))
  }
  idx <- build_local_index(key)
  assign(key, idx, envir = LOCAL_INDEX_CACHE)
  idx
}

local_lookup_by_clid <- function(clid, cl_cfg) {
  idx <- get_local_index(cl_cfg)
  if (is.null(idx) || is.null(idx$clid_to)) return(NA_character_)
  if (!exists(clid, envir = idx$clid_to, inherits = FALSE)) return(NA_character_)
  get(clid, envir = idx$clid_to, inherits = FALSE)
}

local_lookup_by_name <- function(name, cl_cfg) {
  idx <- get_local_index(cl_cfg)
  if (is.null(idx) || is.null(idx$label_to)) return(NULL)
  n <- norm_name(name)
  if (!exists(n, envir = idx$label_to, inherits = FALSE)) return(NULL)
  get(n, envir = idx$label_to, inherits = FALSE)
}

ols_search_cached <- function(q, ontology = "cl", exact = TRUE, rows = 10, include_obsolete = FALSE, queryFields = "label,synonym", cache_dir = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) return(NULL)
  if (is.null(cache_dir) || !nzchar(cache_dir)) cache_dir <- ".ols_cache"
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- paste(q, ontology, exact, rows, include_obsolete, queryFields, sep = "|")
  key <- gsub("[^A-Za-z0-9]+", "_", key)
  key <- substr(key, 1, 120)
  cache_path <- file.path(cache_dir, paste0("ols_", key, ".json"))
  if (file.exists(cache_path)) {
    return(jsonlite::fromJSON(cache_path, simplifyVector = FALSE))
  }
  req <- httr2::request("https://www.ebi.ac.uk/ols4/api/search") |>
    httr2::req_url_query(
      q = q,
      ontology = ontology,
      exact = if (isTRUE(exact)) "true" else "false",
      rows = rows,
      obsoletes = if (isTRUE(include_obsolete)) "true" else "false",
      queryFields = queryFields
    ) |>
    httr2::req_timeout(10)
  resp <- httr2::req_perform(req)
  txt <- httr2::resp_body_string(resp)
  writeLines(txt, cache_path)
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

extract_candidates <- function(js) {
  if (is.null(js$response) || is.null(js$response$docs)) return(tibble::tibble())
  docs <- js$response$docs
  if (length(docs) == 0) return(tibble::tibble())
  tibble::tibble(
    obo_id = purrr::map_chr(docs, ~(.x$obo_id %||% .x$short_form %||% NA_character_)),
    label  = purrr::map_chr(docs, ~(.x$label %||% .x$lbl %||% NA_character_)),
    is_obsolete = purrr::map_lgl(docs, ~{
      v <- .x$is_obsolete %||% .x$obsolete %||% .x$isObsolete
      if (is.null(v) || is.na(v)) FALSE else as.logical(v)
    })
  ) |>
    dplyr::filter(!is.na(.data$obo_id) & !is.na(.data$label)) |>
    dplyr::filter(stringr::str_detect(.data$obo_id, "^CL:"))
}

ols_lookup_by_clid <- function(clid, cl_cfg) {
  js <- ols_search_cached(clid, exact = TRUE, rows = 20, include_obsolete = TRUE, queryFields = "label,synonym,obo_id", cache_dir = cl_cfg$cache_dir)
  if (is.null(js)) return(tibble::tibble())
  cands <- extract_candidates(js)
  if (nrow(cands) == 0) return(cands)
  cands2 <- cands |> dplyr::filter(.data$obo_id == clid)
  if (nrow(cands2) == 0) cands2 <- cands
  cands2
}

ols_lookup_by_name_exact <- function(name, cl_cfg) {
  q <- norm_name(name)
  js <- ols_search_cached(q, exact = TRUE, rows = 20, include_obsolete = FALSE, queryFields = "label,synonym", cache_dir = cl_cfg$cache_dir)
  if (is.null(js)) return(tibble::tibble())
  extract_candidates(js)
}

ols_lookup_by_name_fuzzy <- function(name, cl_cfg) {
  q <- norm_name(name)
  js <- ols_search_cached(q, exact = FALSE, rows = 20, include_obsolete = FALSE, queryFields = "label,synonym", cache_dir = cl_cfg$cache_dir)
  if (is.null(js)) return(tibble::tibble())
  extract_candidates(js)
}

normalize_cl_three_state <- function(name, clid, cl_cfg) {
  name0 <- as.character(name %||% "")
  name_lookup <- canonicalize_name(name0)
  clid0 <- as.character(clid %||% "")
  prefer_ols <- isTRUE(cl_cfg$prefer_ols)
  
  if (nzchar(clid0)) {
    if (prefer_ols) {
      cands <- ols_lookup_by_clid(clid0, cl_cfg)
      if (nrow(cands) > 0 && !isTRUE(cands$is_obsolete[[1]])) {
        lbl <- cands$label[[1]]
        return(list(final_name = lbl %||% name0, final_clid = clid0, status = ifelse(norm_name(lbl) == norm_name(name0), "VALID", "AUTO_CORRECTED")))
      }
    }
    lbl_local <- local_lookup_by_clid(clid0, cl_cfg)
    if (!is.na(lbl_local)) {
      return(list(final_name = lbl_local %||% name0, final_clid = clid0, status = ifelse(norm_name(lbl_local) == norm_name(name0), "VALID", "AUTO_CORRECTED")))
    }
  }
  
  if (nzchar(name_lookup)) {
    if (prefer_ols) {
      cands <- ols_lookup_by_name_exact(name_lookup, cl_cfg)
      if (nrow(cands) > 0) {
        cands2 <- cands |> dplyr::filter(norm_name(.data$label) == norm_name(name_lookup))
        if (nrow(cands2) == 1) {
          return(list(final_name = cands2$label[[1]], final_clid = cands2$obo_id[[1]], status = "AUTO_CORRECTED"))
        }
      }
      cands_f <- ols_lookup_by_name_fuzzy(name_lookup, cl_cfg)
      if (nrow(cands_f) > 0) {
        cands_f <- cands_f |> dplyr::filter(!.data$is_obsolete)
        cands_e <- cands_f |> dplyr::filter(norm_name(.data$label) == norm_name(name_lookup))
        pick <- if (nrow(cands_e) >= 1) cands_e[1, ] else cands_f[1, ]
        return(list(final_name = pick$label[[1]], final_clid = pick$obo_id[[1]], status = "AUTO_CORRECTED"))
      }
    }
    loc <- local_lookup_by_name(name_lookup, cl_cfg)
    if (!is.null(loc)) {
      return(list(final_name = loc$label %||% name0, final_clid = loc$clid %||% clid0, status = "AUTO_CORRECTED"))
    }
  }
  
  list(final_name = name0 %||% NA_character_, final_clid = clid0 %||% NA_character_, status = "UNMAPPED")
}

normalize_section_candidate <- function(sec, cl_cfg) {
  if (is.null(sec) || !is.list(sec)) return(sec)
  res <- normalize_cl_three_state(sec$candidate_cell_type %||% "", sec$cell_ontology_id %||% "", cl_cfg)
  sec$candidate_cell_type <- res$final_name %||% sec$candidate_cell_type
  sec$cell_ontology_id <- res$final_clid %||% sec$cell_ontology_id
  sec
}

normalize_step15_cl <- function(step15_obj, cl_cfg) {
  if (is.null(step15_obj) || !is.list(step15_obj)) return(step15_obj)
  key_main <- if (!is.null(step15_obj$main_type_schema)) "main_type_schema" else "main_type"
  key_sub1 <- if (!is.null(step15_obj$subtype_level_1_schema)) "subtype_level_1_schema" else "subtype_level_1"
  key_sub2 <- if (!is.null(step15_obj$subtype_level_2_schema)) "subtype_level_2_schema" else "subtype_level_2"
  step15_obj[[key_main]] <- normalize_section_candidate(step15_obj[[key_main]], cl_cfg)
  step15_obj[[key_sub1]] <- normalize_section_candidate(step15_obj[[key_sub1]], cl_cfg)
  if (is.list(step15_obj[[key_sub2]])) {
    if (is.list(step15_obj[[key_sub2]]$core_identity)) {
      step15_obj[[key_sub2]]$core_identity <- normalize_section_candidate(step15_obj[[key_sub2]]$core_identity, cl_cfg)
    }
  }
  step15_obj
}

normalize_judge_final_decision_cl <- function(final_obj, cl_cfg) {
  if (is.null(final_obj) || !is.list(final_obj)) return(final_obj)
  if (!is.list(final_obj$final_decision)) return(final_obj)
  fd <- final_obj$final_decision
  res <- normalize_cl_three_state(fd$primary_cell_type %||% "", fd$final_cell_ontology_id %||% "", cl_cfg)
  fd$primary_cell_type <- res$final_name %||% fd$primary_cell_type
  fd$final_cell_ontology_id <- res$final_clid %||% fd$final_cell_ontology_id
  if (!is.null(fd$greedy_cell_type) || !is.null(fd$greedy_cell_ontology_id)) {
    resg <- normalize_cl_three_state(fd$greedy_cell_type %||% "", fd$greedy_cell_ontology_id %||% "", cl_cfg)
    fd$greedy_cell_type <- resg$final_name %||% fd$greedy_cell_type
    fd$greedy_cell_ontology_id <- resg$final_clid %||% fd$greedy_cell_ontology_id
  }
  final_obj$final_decision <- fd
  final_obj
}

# ==============================================================================
# CL ONTOLOGY HIERARCHY FUNCTIONS FOR TRIAGE
# ==============================================================================

#' Load the CL ontology graph with ancestor/descendant relationships
#' @param cl_cfg Configuration list with local_json_path
#' @return Environment with cl (full ontology) and ancestor lookup
#' @noRd
load_cl_graph <- function(cl_cfg) {
  path <- cl_cfg$local_json_path %||% ""
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  cl <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  list(cl = cl)
}

#' Check if clid_a is an ancestor of clid_b (a is higher/more general)
#' @param clid_a The potential ancestor CL ID
#' @param clid_b The potential descendant CL ID
#' @param cl_graph Result from load_cl_graph()
#' @return TRUE if clid_a is ancestor of clid_b, FALSE otherwise
#' @noRd
is_ancestor_of <- function(clid_a, clid_b, cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(FALSE)
  if (is.na(clid_a) || is.na(clid_b) || !nzchar(clid_a) || !nzchar(clid_b)) return(FALSE)
  if (clid_a == clid_b) return(FALSE)
  
  cl <- cl_graph$cl
  # Check if clid_a appears in ancestors of clid_b
  term_b <- cl[[clid_b]]
  if (is.null(term_b)) return(FALSE)
  
  anc <- term_b$ancestors %||% NULL
  if (is.null(anc)) return(FALSE)
  
  # ancestors is a named list/vector where names are ancestor CLIDs
  anc_ids <- names(unlist(anc, use.names = TRUE))
  clid_a %in% anc_ids
}

#' Check if clid_a is a descendant of clid_b (a is lower/more specific)
#' @param clid_a The potential descendant CL ID
#' @param clid_b The potential ancestor CL ID
#' @param cl_graph Result from load_cl_graph()
#' @return TRUE if clid_a is descendant of clid_b, FALSE otherwise
#' @noRd
is_descendant_of <- function(clid_a, clid_b, cl_graph) {
  # a is descendant of b means b is ancestor of a

  is_ancestor_of(clid_b, clid_a, cl_graph)
}

#' Get the distance between two CL IDs in the ontology
#' @param clid_a First CL ID
#' @param clid_b Second CL ID
#' @param cl_graph Result from load_cl_graph()
#' @return Integer distance (1 = direct parent/child), Inf if not related
#' @noRd
get_ontology_distance <- function(clid_a, clid_b, cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(Inf)
  if (is.na(clid_a) || is.na(clid_b) || !nzchar(clid_a) || !nzchar(clid_b)) return(Inf)
  if (clid_a == clid_b) return(0)
  
  cl <- cl_graph$cl
  
  # Check a -> b (a is descendant of b)
  term_a <- cl[[clid_a]]
  if (!is.null(term_a) && !is.null(term_a$ancestors)) {
    anc_a <- unlist(term_a$ancestors, use.names = TRUE)
    if (clid_b %in% names(anc_a)) {
      return(as.numeric(anc_a[[clid_b]]))
    }
  }
  
  # Check b -> a (b is descendant of a)
  term_b <- cl[[clid_b]]
  if (!is.null(term_b) && !is.null(term_b$ancestors)) {
    anc_b <- unlist(term_b$ancestors, use.names = TRUE)
    if (clid_a %in% names(anc_b)) {
      return(as.numeric(anc_b[[clid_a]]))
    }
  }
  
  Inf
}

#' Triage a prediction against ground truth using CL ontology
#' 
#' Classifies the prediction as:
#' - EXACT_MATCH: pred_clid == gt_clid
#' - A_CONSERVATIVE_OK: pred is ancestor of gt (correct direction, not specific enough)
#' - B_OVERSPECIFIC: pred is descendant of gt (too specific for evidence)
#' - C_WRONG_LINEAGE: neither ancestor nor descendant (wrong direction)
#' 
#' NOTE: We do NOT force sibling matches. If the ontology doesn't have a direct
#' ancestor/descendant relationship, we report C_WRONG_LINEAGE and let the
#' feedback loop handle revision.
#' 
#' @param pred_clid Predicted CL ID (e.g., "CL:0000127")
#' @param gt_clid Ground truth CL ID
#' @param cl_graph Result from load_cl_graph()
#' @return List with triage_type, distance, accept (boolean)
#' @noRd
triage_prediction <- function(pred_clid, gt_clid, cl_graph) {
  result <- list(
    triage_type = "C_WRONG_LINEAGE",
    distance = Inf,
    accept = FALSE,
    pred_clid = pred_clid,
    gt_clid = gt_clid
  )
  
  if (is.null(cl_graph) || is.null(cl_graph$cl)) {
    result$notes <- "No CL graph available"
    return(result)
  }
  
  if (is.na(pred_clid) || is.na(gt_clid) || !nzchar(pred_clid) || !nzchar(gt_clid)) {
    result$notes <- "Missing CLID"
    return(result)
  }
  
  # Exact match
  if (pred_clid == gt_clid) {
    return(list(
      triage_type = "EXACT_MATCH",
      distance = 0,
      accept = TRUE,
      pred_clid = pred_clid,
      gt_clid = gt_clid,
      notes = "Exact match"
    ))
  }
  
  # Check if pred is ancestor of gt (conservative - prediction is more general)
  if (is_ancestor_of(pred_clid, gt_clid, cl_graph)) {
    dist <- get_ontology_distance(pred_clid, gt_clid, cl_graph)
    return(list(
      triage_type = "A_CONSERVATIVE_OK",
      distance = dist,
      accept = TRUE,  # Accept conservative predictions
      pred_clid = pred_clid,
      gt_clid = gt_clid,
      notes = paste0("Prediction is ancestor of ground truth (", dist, " hops up)")
    ))
  }
  
  # Check if pred is descendant of gt (overspecific - prediction is more specific than evidence supports)
  if (is_descendant_of(pred_clid, gt_clid, cl_graph)) {
    dist <- get_ontology_distance(pred_clid, gt_clid, cl_graph)
    return(list(
      triage_type = "B_OVERSPECIFIC",
      distance = dist,
      accept = FALSE,  # Need to downgrade
      pred_clid = pred_clid,
      gt_clid = gt_clid,
      notes = paste0("Prediction is descendant of ground truth (", dist, " hops down) - consider downgrade")
    ))
  }
  
  # Neither ancestor nor descendant - wrong lineage
  # Do NOT force match - let the feedback loop handle revision
  list(
    triage_type = "C_WRONG_LINEAGE",
    distance = Inf,
    accept = FALSE,
    pred_clid = pred_clid,
    gt_clid = gt_clid,
    notes = "No ancestor/descendant relationship - different lineage"
  )
}

#' Apply triage logic to handle revision_directive from Step 2 validator
#' 
#' @param step2_obj The Step 2 validation output object
#' @param step15_obj The Step 1.5 report object
#' @param cl_cfg CL configuration
#' @return List with final_accept, final_report, triage_info
#' @noRd
apply_revision_directive_triage <- function(step2_obj, step15_obj, cl_cfg) {
  result <- list(
    final_accept = FALSE,
    final_report = step15_obj,
    triage_info = NULL,
    action_taken = "none"
  )
  
  # Check validation status
  status <- step2_obj$validation_status %||% ""
  if (grepl("PASSED", toupper(status))) {
    result$final_accept <- TRUE
    result$action_taken <- "passed"
    return(result)
  }
  
  # Get revision_directive if present
  rd <- step2_obj$revision_directive %||% NULL
  
  if (!is.null(rd)) {
    failure_type <- rd$failure_type %||% ""
    action <- rd$action %||% ""
    suggested_parent <- rd$suggested_parent %||% NULL
    suggested_parent_clid <- rd$suggested_parent_clid %||% NULL
    
    result$triage_info <- rd
    
    # Handle A_conservative_ok - accept as is
    if (failure_type == "A_conservative_ok" || action == "accept_as_is") {
      result$final_accept <- TRUE
      result$action_taken <- "accepted_conservative"
      return(result)
    }
    
    # Handle B_overspecific - downgrade to parent
    if (failure_type == "B_overspecific" && action == "downgrade_to_parent") {
      if (!is.null(suggested_parent) && nzchar(suggested_parent)) {
        # Apply the downgrade
        downgraded <- downgrade_report_to_parent(step15_obj, suggested_parent, suggested_parent_clid, cl_cfg)
        result$final_report <- downgraded
        result$final_accept <- TRUE
        result$action_taken <- "downgraded_to_parent"
        return(result)
      }
    }
    
    # Handle C_wrong_lineage - reject
    if (failure_type == "C_wrong_lineage" || action == "stop") {
      result$final_accept <- FALSE
      result$action_taken <- "rejected_wrong_lineage"
      return(result)
    }
  }
  
  # No revision_directive or unhandled case - use original pass/fail
  result$action_taken <- "no_directive"
  result
}

#' Downgrade a report's core identity to a parent term
#' 
#' @param report The Step 1.5 report object
#' @param parent_label The parent cell type label
#' @param parent_clid The parent CL ID (optional)
#' @param cl_cfg CL configuration
#' @return Modified report with downgraded identity
#' @noRd
downgrade_report_to_parent <- function(report, parent_label, parent_clid = NULL, cl_cfg = NULL) {
  if (is.null(report) || !is.list(report)) return(report)
  
  # Resolve CLID if not provided
  if (is.null(parent_clid) || !nzchar(parent_clid)) {
    if (!is.null(cl_cfg)) {
      m <- normalize_cl_three_state(parent_label, "", cl_cfg)
      parent_clid <- m$final_clid %||% ""
    }
  }
  
  # Update all core identity fields
  main_key <- if (!is.null(report$main_type)) "main_type" else if (!is.null(report$main_type_schema)) "main_type_schema" else NULL
  sub1_key <- if (!is.null(report$subtype_level_1)) "subtype_level_1" else if (!is.null(report$subtype_level_1_schema)) "subtype_level_1_schema" else NULL
  sub2_key <- if (!is.null(report$subtype_level_2)) "subtype_level_2" else if (!is.null(report$subtype_level_2_schema)) "subtype_level_2_schema" else NULL
  
  if (!is.null(main_key)) {
    report[[main_key]]$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) report[[main_key]]$cell_ontology_id <- parent_clid
  }
  
  if (!is.null(sub1_key)) {
    report[[sub1_key]]$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) report[[sub1_key]]$cell_ontology_id <- parent_clid
  }
  
  if (!is.null(sub2_key) && is.list(report[[sub2_key]]$core_identity)) {
    report[[sub2_key]]$core_identity$candidate_cell_type <- parent_label
    if (nzchar(parent_clid)) report[[sub2_key]]$core_identity$cell_ontology_id <- parent_clid
  }
  
  if (!is.null(report$open_world_summary)) {
    report$open_world_summary$best_cell_type <- parent_label
    if (nzchar(parent_clid)) report$open_world_summary$cell_ontology_id <- parent_clid
  }
  
  # Add a note about the downgrade
  report$downgrade_applied <- list(
    original_prediction = report$open_world_summary$best_cell_type %||% "unknown",
    downgraded_to = parent_label,
    downgraded_clid = parent_clid,
    reason = "B_overspecific - evidence insufficient for specific subtype"
  )
  
  report
}
