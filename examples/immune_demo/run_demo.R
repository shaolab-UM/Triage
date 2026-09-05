# immune_demo: deterministic, no-API adjudication of one cluster.
# Run from the repository root: source("examples/immune_demo/run_demo.R")

suppressMessages(library(Triage))

ontology <- load_triage_ontology()

jin <- read_triage_input(
  file.path("examples", "immune_demo", "cluster_1_round1.json"))

fin <- run_triage_adjudication(
  jin,
  ontology = ontology,
  use_api = FALSE,
  head_output = file.path("examples", "immune_demo", "head_round1.json"),
  dataset_name = "census_immune",
  project_root = file.path("reproducibility", "primary")
)

cat("Local gate:", isTRUE(attr(fin, "gate")$ok), "\n")
cat("Final decision:", fin$final_decision$primary_cell_type,
    "|", fin$final_decision$final_cell_ontology_id,
    "|", fin$final_decision$decision_category, "\n")

validate_triage_result(fin, jin)
cat("validate_triage_result: OK\n")
