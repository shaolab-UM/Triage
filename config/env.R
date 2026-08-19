# Optional centralized environment/path helpers for public release.
# Pipeline scripts also accept the same environment variables directly.

`%||%` <- function(a, b) if (!is.null(a)) a else b

project_root <- normalizePath(
  Sys.getenv("TRIAGE_HOME", unset = Sys.getenv("PROJECT_ROOT", unset = getwd())),
  winslash = "/", mustWork = FALSE
)

resource_root <- Sys.getenv(
  "TRIAGE_RESOURCE_ROOT",
  unset = file.path(project_root, "inputs", "raw")
)

paths <- list(
  cl_ontology      = file.path(resource_root, "ontology", "CL-ontology-v2025-07-30.json"),
  ppi_aliases      = file.path(resource_root, "ppi", "{ppi_code}.protein.aliases.v12.0.txt"),
  ppi_interactions = file.path(resource_root, "ppi", "{ppi_code}.protein.physical.links.v12.0.txt"),
  cell_marker_h    = file.path(resource_root, "cellmarker", "Cell_marker_Human.xlsx"),
  cell_marker_m    = file.path(resource_root, "cellmarker", "Cell_marker_Mouse.xlsx"),
  collectri_h      = file.path(resource_root, "collectri", "collectri_human_network.rds"),
  collectri_m      = file.path(resource_root, "collectri", "collectri_mouse_network.rds")
)

api <- list(
  key_env  = Sys.getenv("LLM_API_KEY_ENV", unset = Sys.getenv("CASSIA_API_KEY_ENV", unset = "DEEPSEEK_API_KEY")),
  base_url = Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = "XXXXX")),
  disgenet = Sys.getenv("DISGENET_API_KEY", unset = "")
)

normalize_provider <- function(url) {
  if (grepl("/chat/completions$", url)) {
    url <- sub("/v1/chat/completions$", "", url)
    url <- sub("/chat/completions$", "", url)
  }
  url
}
api$provider <- normalize_provider(api$base_url)

ppi_path <- function(ppi_code) {
  list(
    aliases = sub("{ppi_code}", ppi_code, paths$ppi_aliases, fixed = TRUE),
    interactions = sub("{ppi_code}", ppi_code, paths$ppi_interactions, fixed = TRUE)
  )
}

outputs_root <- file.path(project_root, "outputs")
run_tag <- Sys.getenv("RUN_TAG", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
dataset_name <- Sys.getenv("DATASET_NAME", unset = "")

dataset_out_dir <- function(ds = dataset_name, tag = run_tag) {
  if (!nzchar(ds)) stop("DATASET_NAME is not set")
  file.path(outputs_root, ds, tag)
}
