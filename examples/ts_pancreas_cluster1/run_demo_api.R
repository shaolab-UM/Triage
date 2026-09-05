# ts_pancreas_cluster1: Handling Editor API call for one TS_pancreas
# cluster. Requires DEEPSEEK_API_KEY (and optionally LLM_API_BASE_URL).
# This script makes real LLM calls; it is NOT part of the test suite.
# Chief QC is not invoked by the package.

suppressMessages(library(Triage))

if (!nzchar(Sys.getenv("DEEPSEEK_API_KEY"))) {
  stop("ts_pancreas_cluster1: set DEEPSEEK_API_KEY (see examples/ts_pancreas_cluster1/README.md).")
}

ontology <- load_triage_ontology()

jin <- read_triage_input(
  file.path("examples", "ts_pancreas_cluster1", "cluster_1_round1.json"))

fin <- run_triage_adjudication(
  jin,
  ontology = ontology,
  use_api = TRUE,
  model = "deepseek-v4-flash",   # package default
  temperature = 0,
  prompt_profile = "compact",
  dataset_name = "ts_pancreas",
  project_root = file.path("reproducibility", "primary")
)

cat("Local gate:", isTRUE(attr(fin, "gate")$ok), "\n")
cat("Final decision:", fin$final_decision$primary_cell_type,
    "|", fin$final_decision$final_cell_ontology_id, "\n")
