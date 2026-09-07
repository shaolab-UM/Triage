# Triage package tests (self-attaching so test_dir works standalone)
library(Triage)

# Deterministic tests only: no network, no API keys, no LLM calls.

test_that("packaged example fixtures ship with the installed package", {
  input_path <- system.file("extdata", "examples", "census_immune_cluster1",
                            "cluster_1_round1.json", package = "Triage")
  head_path <- system.file("extdata", "examples", "census_immune_cluster1",
                           "head_round1.json", package = "Triage")
  expect_true(nzchar(input_path) && file.exists(input_path))
  expect_true(nzchar(head_path) && file.exists(head_path))
})

test_that("run_triage_example reproduces the reported identity without a repository or API", {
  fin <- withr::with_dir(tempdir(), run_triage_example())
  expect_true(is.list(fin))
  expect_true(isTRUE(attr(fin, "gate")$ok))
  expect_equal(fin$final_decision$primary_cell_type, "classical monocyte")
  expect_equal(fin$final_decision$final_cell_ontology_id, "CL:0000860")
})
