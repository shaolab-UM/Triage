# Triage package tests (self-attaching so test_dir works standalone)
library(Triage)

# Deterministic tests only: no network, no API keys, no LLM calls.

test_that("package loads and exports the public API", {
  expect_true("Triage" %in% rownames(utils::installed.packages()))
  for (fn in c("read_triage_input", "run_cl_linker",
               "build_adjudication_input", "run_triage_adjudication",
               "validate_triage_result", "load_triage_ontology",
               "run_triage_example", "install_triage_dependencies",
               "setup_triage_resources", "run_triage", "triage_preflight")) {
    expect_true(exists(fn, where = "package:Triage"), info = fn)
  }
  for (fn in c("map_cell_ontology", "build_reviewer_summary",
               "default_ontology_path")) {
    expect_true(exists(fn, where = asNamespace("Triage")), info = fn)
  }
})

test_that("ontology loads from the shipped fixture and maps labels", {
  o <- load_triage_ontology()
  expect_true(is.list(o))
  res <- Triage:::map_cell_ontology("monocyte", ontology = o)
  expect_equal(res$clid, "CL:0000576")
  expect_true(nzchar(res$status))
})

test_that("CL-Linker maps labels with and without provided ids", {
  out <- run_cl_linker(c("macrophage", "T cell"), c("", "CL:0000084"))
  expect_type(out, "list")
  expect_length(out, 2)
  expect_true(grepl("^CL:", out[[1]]$cl_id %||% out[[1]]$clid %||% ""))
})

test_that("reviewer summaries parse and validate", {
  rs <- Triage:::build_reviewer_summary(list(
    cassia = list(top1_cell_type = "classical monocyte",
                  topk_cell_types = c("classical monocyte", "intermediate monocyte"),
                  reasoning_short = "marker match"),
    in_house = list(top1_cell_type = "monocyte",
                    cell_ontology_id = "CL:0000576")
  ))
  expect_equal(rs$cassia$top1_cell_type, "classical monocyte")
  expect_error(Triage:::build_reviewer_summary(list(cassia = list())),
               regexp = "top1_cell_type|top1")
})

test_that("schema parsing extracts JSON from fenced LLM text", {
  fenced <- "Here is the output:\n```json\n{\"cluster_id\": \"c1\", \"final_decision\": {\"primary_cell_type\": \"monocyte\"}}\n```\nThanks."
  obj <- jsonlite::fromJSON(Triage:::extract_first_json_object_stack(fenced), simplifyVector = FALSE)
  expect_equal(obj$cluster_id, "c1")
  expect_true(is.na(Triage:::extract_first_json_object_stack("no json here at all")))
})

test_that("no-API adjudication reproduces the head decision on the fixture", {
  jin <- read_triage_input(test_path("fixtures", "cluster_1_round1.json"))
  fin <- run_triage_adjudication(jin, use_api = FALSE,
                                 head_output = test_path("fixtures", "head_round1.json"),
                                 dataset_name = "census_immune",
                                 project_root = test_path("fixtures"))
  expect_true(isTRUE(attr(fin, "gate")$ok))
  expect_equal(fin$final_decision$primary_cell_type, "classical monocyte")
  expect_equal(fin$final_decision$final_cell_ontology_id, "CL:0000860")
  expect_silent(validate_triage_result(fin, jin))
})

test_that("invalid inputs fail with informative errors", {
  expect_error(read_triage_input("does_not_exist.json"), regexp = "not found|exist")
  expect_error(read_triage_input('{"foo": 1}'), regexp = "cluster_id|inputs")
  expect_error(suppressWarnings(run_triage_adjudication(list(cluster_id = "x"), use_api = FALSE)),
               regexp = "head_output|head")
})
