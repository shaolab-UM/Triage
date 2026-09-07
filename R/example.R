#' Run the bundled deterministic example
#'
#' Runs the packaged Census_immune cluster 1 example end to end inside the
#' installed package: loads the bundled Cell Ontology, reads the packaged
#' reviewer-input fixture, and executes the deterministic adjudication chain
#' (`use_api = FALSE`, precomputed Handling Editor fixture). No network
#' access and no API key are required, and no repository checkout is needed
#' — all inputs are located with [system.file()].
#'
#' The example reproduces the reported primary-benchmark identity for
#' Census_immune cluster 1: classical monocyte (CL:0000860).
#'
#' @return The final adjudication record (list) with attribute `gate`
#'   (local gate result), invisibly.
#' @export
run_triage_example <- function() {
  input_path <- system.file("extdata", "examples", "census_immune_cluster1",
                            "cluster_1_round1.json", package = "Triage")
  head_path <- system.file("extdata", "examples", "census_immune_cluster1",
                           "head_round1.json", package = "Triage")
  if (!nzchar(input_path) || !file.exists(input_path) ||
      !nzchar(head_path) || !file.exists(head_path)) {
    stop("run_triage_example: packaged example fixtures are missing from the installed package.")
  }

  message("Loading the bundled Cell Ontology ...")
  ontology <- load_triage_ontology()

  message("Reading the Census_immune cluster 1 reviewer-input fixture ...")
  jin <- read_triage_input(input_path)

  message("Running deterministic adjudication (no API) ...")
  fin <- run_triage_adjudication(
    jin,
    ontology = ontology,
    use_api = FALSE,
    head_output = head_path,
    dataset_name = "census_immune",
    project_root = tempdir()
  )

  if (!isTRUE(attr(fin, "gate")$ok)) {
    stop("run_triage_example: the deterministic local gate did not pass.")
  }
  validate_triage_result(fin, jin)

  fd <- fin$final_decision
  message("Local gate: OK")
  message("Result: ", fd$primary_cell_type,
          " | ", fd$final_cell_ontology_id)
  invisible(fin)
}
