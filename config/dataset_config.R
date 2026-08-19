`%||%` <- function(a, b) if (!is.null(a)) a else b

# Developmental-stage ontology roots for rule policies.
# Keep this list as the single source of truth for stage filtering logic.
STAGE_ROOT_CLIDS <- c(
  "CL:0000034", # stem cell
  "CL:0000037"  # hematopoietic stem cell
)

# Strict mode: stage policy is CL-based by default.
# Token fallback is disabled unless explicitly enabled.
ALLOW_STAGE_TOKEN_FALLBACK <- FALSE

get_dataset_config <- function(dataset_name = basename(getwd()), project_root = getwd()) {
  as_char_vec <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NULL)
    unique(trimws(x))
  }

  load_context_override <- function() {
    p <- Sys.getenv("DATASET_CONTEXT_JSON_PATH", unset = Sys.getenv("CASSIA_CONTEXT_JSON_PATH", unset = ""))
    if (!nzchar(p) || !file.exists(p)) return(list())
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      warning("DATASET_CONTEXT_JSON_PATH set but jsonlite unavailable; ignoring context override.")
      return(list())
    }
    obj <- tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.list(obj)) {
      warning("Failed to parse context JSON override: ", p)
      return(list())
    }
    obj
  }

  resolve_dataset_dir_name <- function(name, root) {
    final_root <- file.path(root, "final")
    if (!dir.exists(final_root)) return(name)
    dirs <- list.dirs(final_root, full.names = FALSE, recursive = FALSE)
    if (name %in% dirs) return(name)
    match <- dirs[tolower(dirs) == tolower(trimws(name))]
    if (length(match) > 0) return(match[[1]])
    name
  }

  dataset_name <- resolve_dataset_dir_name(dataset_name, project_root)
  dataset_key <- tolower(trimws(dataset_name))

  default_context <- list(
    species = c("human"),
    tissue = c("tissue"),
    data_type = c("single-cell RNA"),
    study_context = c("reference"),
    allowed_lineages = NULL,
    user_notes = paste(
      "Human single-cell RNA dataset.",
      "Goal: marker-based candidate cell type annotation at cluster/state level.",
      "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
    )
  )

  scope_profiles <- list(
    immune_enriched = list(
      allowed_lineages = c("immune", "erythroid_megakaryocytic", "unknown"),
      gate_mode = "flag_only"
    ),
    whole_tissue = list(
      allowed_lineages = c("immune", "epithelial", "endothelial", "stromal_mesenchymal", "muscle", "neural_glial", "erythroid_megakaryocytic", "unknown"),
      gate_mode = "flag_only"
    ),
    bone_marrow = list(
      allowed_lineages = c("immune", "erythroid_megakaryocytic", "stromal_mesenchymal", "endothelial", "unknown"),
      gate_mode = "flag_only"
    ),
    tme = list(
      allowed_lineages = c("immune", "epithelial", "stromal_mesenchymal", "endothelial", "muscle", "unknown"),
      gate_mode = "flag_only"
    ),
    mixed_unknown = list(
      allowed_lineages = NULL,
      gate_mode = "flag_only"
    )
  )

  dataset_context_map <- list(
    census_immune = list(
      tissue = c("blood"),
      study_context = c("healthy_immune"),
      dataset_scope = "immune_enriched",
      scope_profile = "immune_enriched",
      user_notes = paste(
        "Human blood single-cell RNA dataset (Census immune).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    zheng_blood = list(
      tissue = c("blood"),
      study_context = c("healthy_immune"),
      dataset_scope = "immune_enriched",
      scope_profile = "immune_enriched",
      user_notes = paste(
        "Human blood single-cell RNA dataset (Zheng).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    ts_pancreas = list(
      tissue = c("pancreas"),
      study_context = c("normal_adult"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human pancreas single-cell RNA dataset (Tabula Sapiens).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    ts_kidney = list(
      tissue = c("kidney"),
      study_context = c("normal_adult"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human kidney single-cell RNA dataset (Tabula Sapiens).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    sikkema_lung = list(
      tissue = c("lung"),
      study_context = c("health_and_disease_reference_atlas"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human lung single-cell RNA dataset (Sikkema).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    sikkema_lung_rawdeg = list(
      tissue = c("lung"),
      study_context = c("health_and_disease_reference_atlas"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human lung single-cell RNA dataset (Sikkema, rawdeg compare).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    ts_pancreas_rawdeg = list(
      tissue = c("pancreas"),
      study_context = c("normal_adult"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human pancreas single-cell RNA dataset (Tabula Sapiens, rawdeg compare).",
        "Goal: marker-based candidate cell type annotation at cluster/state level.",
        "Context is a soft prior; scope gating is flag-only: out-of-scope evidence is contamination-only and must not change core identity."
      )
    ),
    validation_test = list(
      species = c("human"),
      tissue = c("dorsolateral prefrontal cortex", "DLPFC"),
      data_type = c("single nucleus multiome", "snRNA-seq", "snATAC-seq"),
      study_context = c("alzheimer_disease_vs_control", "brain_multiomics"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Human DLPFC single-nucleus multiomics validation dataset (AD vs Control).",
        "Source context: Anderson et al., Cell Genomics 2023 (ZEB1/MAFB regulatory analysis).",
        "Use brain whole-tissue scope; apply flag-only gating so out-of-scope evidence is contamination-only and must not override core identity."
      )
    ),
    screview_clean = list(
      species = c("human"),
      tissue = c("dorsolateral prefrontal cortex", "DLPFC", "brain"),
      data_type = c("single nucleus multiome", "snRNA-seq", "snATAC-seq"),
      study_context = c("alzheimer_disease_vs_control"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Held-out human DLPFC validation dataset (screview_clean).",
        "Same brain whole-tissue scope as validation_test; flag-only gating."
      )
    ),
    screview_mouse = list(
      species = c("mouse"),
      tissue = c("hippocampi", "brain"),
      data_type = c("single-cell RNA"),
      study_context = c("neuro_degenerative"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Held-out mouse hippocampus validation dataset (screview_mouse, 5xFAD AD model).",
        "Mouse species; gene symbols are title-case. Whole-tissue scope; flag-only gating."
      )
    ),
    screview_ture = list(
      species = c("human"),
      tissue = c("cortex organoid", "organoid", "cerebral cortex"),
      data_type = c("single-cell RNA"),
      study_context = c("organoid_development"),
      dataset_scope = "whole_tissue",
      scope_profile = "whole_tissue",
      user_notes = paste(
        "Held-out human cerebral cortex organoid validation dataset (screview_ture).",
        "Developing organoid; includes retinal/progenitor populations. Whole-tissue scope; flag-only gating."
      )
    )
  )

  context_override <- load_context_override()

  if (is.null(dataset_context_map[[dataset_key]]) && length(context_override) == 0) {
    warning("Unrecognized dataset_key; falling back to default context and mixed_unknown scope.")
  }

  env_species <- as_char_vec(Sys.getenv("DATASET_SPECIES", unset = Sys.getenv("CASSIA_SPECIES", unset = "")))
  env_tissue <- as_char_vec(Sys.getenv("DATASET_TISSUE", unset = Sys.getenv("CASSIA_TISSUE", unset = "")))
  env_data_type <- as_char_vec(Sys.getenv("DATASET_DATA_TYPE", unset = Sys.getenv("CASSIA_DATA_TYPE", unset = "")))
  env_study_context <- as_char_vec(Sys.getenv("DATASET_STUDY_CONTEXT", unset = Sys.getenv("CASSIA_STUDY_CONTEXT", unset = "")))
  env_user_notes <- Sys.getenv("DATASET_USER_NOTES", unset = Sys.getenv("CASSIA_USER_NOTES", unset = ""))
  env_dataset_scope <- Sys.getenv("DATASET_SCOPE", unset = Sys.getenv("CASSIA_DATASET_SCOPE", unset = ""))
  env_scope_profile <- Sys.getenv("DATASET_SCOPE_PROFILE", unset = Sys.getenv("CASSIA_SCOPE_PROFILE", unset = ""))
  env_gate_mode <- Sys.getenv("DATASET_GATE_MODE", unset = Sys.getenv("CASSIA_GATE_MODE", unset = ""))
  env_allowed_lineages <- as_char_vec(strsplit(Sys.getenv("DATASET_ALLOWED_LINEAGES", unset = Sys.getenv("CASSIA_ALLOWED_LINEAGES", unset = "")), "[,;]")[[1]])

  context <- modifyList(default_context, dataset_context_map[[dataset_key]] %||% list())
  context <- modifyList(context, context_override)

  if (!is.null(env_species)) context$species <- env_species
  if (!is.null(env_tissue)) context$tissue <- env_tissue
  if (!is.null(env_data_type)) context$data_type <- env_data_type
  if (!is.null(env_study_context)) context$study_context <- env_study_context
  if (nzchar(env_user_notes)) context$user_notes <- env_user_notes
  if (nzchar(env_dataset_scope)) context$dataset_scope <- env_dataset_scope
  if (nzchar(env_scope_profile)) context$scope_profile <- env_scope_profile
  if (!is.null(env_allowed_lineages)) context$allowed_lineages <- env_allowed_lineages
  if (nzchar(env_gate_mode)) context$gate_mode <- env_gate_mode

  dataset_scope <- context$dataset_scope %||% "mixed_unknown"
  scope_profile <- context$scope_profile %||% dataset_scope
  profile_cfg <- scope_profiles[[scope_profile]] %||% list()

  allowed_lineages <- context$allowed_lineages %||% profile_cfg$allowed_lineages %||% NULL
  gate_mode <- context$gate_mode %||% profile_cfg$gate_mode %||% "flag_only"
  allowed_lineages_source <- if (!is.null(context$allowed_lineages)) {
    "dataset"
  } else if (!is.null(profile_cfg$allowed_lineages)) {
    "profile"
  } else {
    "default"
  }
  gate_mode_source <- if (!is.null(context$gate_mode)) {
    "dataset"
  } else if (!is.null(profile_cfg$gate_mode)) {
    "profile"
  } else {
    "default"
  }

  list(
    dataset_name = dataset_name,
    dataset_key = dataset_key,
    species = context$species,
    tissue = context$tissue,
    data_type = context$data_type,
    study_context = context$study_context,
    user_notes = context$user_notes,
    dataset_scope = dataset_scope,
    scope_profile = scope_profile,
    allowed_lineages = allowed_lineages,
    allowed_lineages_source = allowed_lineages_source,
    gate_mode = gate_mode,
    gate_mode_source = gate_mode_source,
    deg_file = file.path(project_root, "final", dataset_name, "maskdeg.csv"),
    candidates_file = file.path(project_root, "final", dataset_name, "cellanno", "AD_structured_candidates.csv"),
    true_label_csv = file.path(project_root, "final", dataset_name, "true_label.csv"),
    llm_inputs_root = file.path(project_root, "final", dataset_name, "llm_inputs"),
    llm_outputs_root = file.path(project_root, "final", dataset_name, "llm_outputs"),
    intermediate_outputs_root = file.path(project_root, "final", dataset_name, "intermediate_outputs"),
    judge_outputs_root = file.path(project_root, "final", dataset_name, "llm_judge_outputs"),
    cassia_results_glob = file.path(project_root, "final", dataset_name, "CASSIA_*_*_*", "01_annotation_results", "annotation_cassia_FINAL_RESULTS.csv")
  )
}
