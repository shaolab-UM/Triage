# Triage package tests (self-attaching so test_dir works standalone)
library(Triage)

# Deterministic tests only: no network, no API keys, no LLM calls.
#
# Regression guard for apply_consensus_subtype_policy(): the stage-policy
# switches STAGE_ROOT_CLIDS / ALLOW_STAGE_TOKEN_FALLBACK are package-internal
# constants (R/dataset_config.R). The function previously resolved them via
# get0(..., inherits = TRUE), a string-based dynamic lookup that silently
# succeeded only when a caller leaked the constants into the global
# environment; after the config refactor it found nothing and hit the strict
# stop "STAGE_ROOT_CLIDS is empty or unavailable; strict CL-based stage policy
# cannot run". These tests pin the fixed lexical-namespace lookup: the policy
# must work with the constants ABSENT from (and, in the adversarial case,
# CONFLICTING in) .GlobalEnv.

test_that("stage-policy constants stay package-internal with pinned values", {
  expect_false(exists("STAGE_ROOT_CLIDS", envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists("ALLOW_STAGE_TOKEN_FALLBACK", envir = .GlobalEnv, inherits = FALSE))
  expect_false("STAGE_ROOT_CLIDS" %in% getNamespaceExports("Triage"))
  expect_false("ALLOW_STAGE_TOKEN_FALLBACK" %in% getNamespaceExports("Triage"))
  expect_identical(get("STAGE_ROOT_CLIDS", envir = asNamespace("Triage")),
                   c("CL:0000034", "CL:0000037"))
  expect_false(isTRUE(get("ALLOW_STAGE_TOKEN_FALLBACK", envir = asNamespace("Triage"))))
})

test_that("consensus stage policy resolves configured CL roots without .GlobalEnv", {
  ontology <- load_triage_ontology()
  # Same filter apply_consensus_subtype_policy applies (adjudication.R: the
  # configured roots must survive intersection with the loaded CL graph), so
  # the strict-stop branch must be unreachable for the shipped ontology.
  resolved <- Triage:::STAGE_ROOT_CLIDS[Triage:::STAGE_ROOT_CLIDS %in% names(ontology$graph$cl)]
  expect_true(length(resolved) > 0)

  jin_path <- test_path("fixtures", "cluster_1_round1.json")
  head_path <- test_path("fixtures", "head_round1.json")

  # End-to-end deterministic postprocessing (use_api = FALSE): must not throw
  # "STAGE_ROOT_CLIDS is empty or unavailable".
  out <- run_triage_adjudication(input = jin_path, ontology = ontology,
                                 use_api = FALSE, head_output = head_path)
  expect_true(is.list(out))
  gate <- attr(out, "gate")
  expect_true(is.list(gate))

  # Direct policy invocation with hostile-free globals (belt and braces).
  jin <- Triage:::read_triage_input(jin_path)
  head_out <- Triage:::read_json_safely(head_path)
  pol <- Triage:::apply_consensus_subtype_policy(head_out, jin, ontology$cfg,
                                                 ontology$graph,
                                                   dataset_cfg = Triage:::get_dataset_config("unknown"))
  expect_true(is.list(pol))
  expect_false(identical(pol$final_decision, "ERROR"))
})

test_that("policy ignores hostile .GlobalEnv copies of the constants", {
  ontology <- load_triage_ontology()
  jin <- Triage:::read_triage_input(test_path("fixtures", "cluster_1_round1.json"))
  head_out <- Triage:::read_json_safely(test_path("fixtures", "head_round1.json"))

  had_root <- exists("STAGE_ROOT_CLIDS", envir = .GlobalEnv, inherits = FALSE)
  had_fallback <- exists("ALLOW_STAGE_TOKEN_FALLBACK", envir = .GlobalEnv, inherits = FALSE)
  prev_root <- if (had_root) get("STAGE_ROOT_CLIDS", envir = .GlobalEnv) else NULL
  prev_fallback <- if (had_fallback) get("ALLOW_STAGE_TOKEN_FALLBACK", envir = .GlobalEnv) else NULL
  assign("STAGE_ROOT_CLIDS", character(0), envir = .GlobalEnv)
  assign("ALLOW_STAGE_TOKEN_FALLBACK", TRUE, envir = .GlobalEnv)
  cleanup <- function() {
    if (had_root) assign("STAGE_ROOT_CLIDS", prev_root, envir = .GlobalEnv) else
      rm("STAGE_ROOT_CLIDS", envir = .GlobalEnv)
    if (had_fallback) assign("ALLOW_STAGE_TOKEN_FALLBACK", prev_fallback, envir = .GlobalEnv) else
      rm("ALLOW_STAGE_TOKEN_FALLBACK", envir = .GlobalEnv)
  }

  tryCatch({
    # Lexical package lookup must win: empty hostile roots would trigger the
    # strict stop; a hostile fallback=TRUE must not change the policy.
    pol <- Triage:::apply_consensus_subtype_policy(head_out, jin, ontology$cfg,
                                                   ontology$graph,
                                                   dataset_cfg = Triage:::get_dataset_config("unknown"))
    expect_true(is.list(pol))
  }, finally = cleanup())

  expect_false(exists("STAGE_ROOT_CLIDS", envir = .GlobalEnv, inherits = FALSE))
  expect_false(exists("ALLOW_STAGE_TOKEN_FALLBACK", envir = .GlobalEnv, inherits = FALSE))
})
