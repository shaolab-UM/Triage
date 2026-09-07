# ts_pancreas_cluster1: Handling Editor API call for one TS_pancreas
# cluster. Supply your API key and chat-completions endpoint through the
# api_key / api_base_url variables below, or set them to NULL to use the
# DEEPSEEK_API_KEY / LLM_API_BASE_URL environment variables.
# This script makes real LLM calls; it is NOT part of the test suite.

suppressMessages(library(Triage))

api_key <- "YOUR_API_KEY"                                     # <- replace me
api_base_url <- "https://api.deepseek.com/chat/completions"   # <- or NULL for env var

if (identical(api_key, "YOUR_API_KEY")) {
  stop("ts_pancreas_cluster1: replace api_key at the top of ",
       "examples/ts_pancreas_cluster1/run_demo_api.R ",
       "(see examples/ts_pancreas_cluster1/README.md).")
}

ontology <- load_triage_ontology()

jin <- read_triage_input(
  file.path("examples", "ts_pancreas_cluster1", "cluster_1_round1.json"))

fin <- run_triage_adjudication(
  jin,
  ontology = ontology,
  use_api = TRUE,
  api_key = api_key,
  api_base_url = api_base_url,
  temperature = 0,
  prompt_profile = "compact",
  dataset_name = "ts_pancreas",
  project_root = file.path("reproducibility", "primary")
)

cat("Local gate:", isTRUE(attr(fin, "gate")$ok), "\n")
cat("Final decision:", fin$final_decision$primary_cell_type,
    "|", fin$final_decision$final_cell_ontology_id, "\n")
