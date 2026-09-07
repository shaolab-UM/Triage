# Review-round regression tests: API fallback canonicalization,
# environment restoration, stage-10 confidence extraction, and
# child-Rscript path quoting. No LLM API calls are made.

library(Triage)

skip_if_runner_scripts_unavailable <- function() {
  testthat::skip_if_not(
    file.exists(test_path("..", "..", "reproducibility", "scripts",
                          "run_triage.R")),
    "Repository-level runner scripts are excluded from the R package bundle.")
}

# Build a staged resource root that passes the format-aware preflight.
triage_staged_resource_root <- function(species = "human") {
  root <- tempfile("triage_res_")
  dir.create(file.path(root, "ppi"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "collectri"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "cellmarker"), recursive = TRUE, showWarnings = FALSE)
  code <- if (species == "human") "9606" else "10090"
  writeLines(paste0("ensp", 1:3, "\tGeneA\t5"),
             file.path(root, "ppi",
                       paste0(code, ".protein.aliases.v12.0.txt")))
  writeLines("ensp1 ensp2 900",
             file.path(root, "ppi",
                       paste0(code, ".protein.physical.links.v12.0.txt")))
  net <- data.frame(source = c("A", "B"), target = c("C", "C"),
                    weight = c(1, 2))
  saveRDS(net, file.path(root, "collectri",
                         paste0("collectri_", species, "_network.rds")))
  cm <- data.frame(cell_name = c("monocyte", "T cell"),
                   marker = c("GENE1", "GENE2"))
  writexl::write_xlsx(cm,
                      file.path(root, "cellmarker",
                                if (species == "human") "Cell_marker_Human.xlsx"
                                else "Cell_marker_Mouse.xlsx"))
  root
}

triage_deg_fixture <- function() {
  deg <- tempfile("deg_", fileext = ".csv")
  writeLines(c("cluster,gene,avg_log2FC,p_val,p_val_adj,pct.1,pct.2",
               "Cluster A,Gene1,1.5,0.01,0.02,0.5,0.1",
               "Cluster B,Gene2,2.0,0.01,0.03,0.6,0.1"), deg)
  deg
}

test_that(".triage_resolve_api canonicalizes the primary variables directly", {
  eff <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "", LLM_API_BASE_URL = "", CASSIA_API_BASE_URL = "",
      LLM_API_KEY_ENV = ""),
    Triage:::.triage_resolve_api("k-direct", "https://e.example/chat/completions"))
  expect_identical(eff$key, "k-direct")
  expect_identical(eff$url, "https://e.example/chat/completions")
})

test_that(".triage_resolve_api honors a custom LLM_API_KEY_ENV key variable", {
  eff <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "", LLM_API_KEY_ENV = "MY_API_KEY",
      MY_API_KEY = "k-custom",
      LLM_API_BASE_URL = "https://e.example/chat/completions"),
    Triage:::.triage_resolve_api(NULL, NULL))
  expect_identical(eff$key, "k-custom")
  expect_identical(eff$url, "https://e.example/chat/completions")
})

test_that(".triage_resolve_api falls back to CASSIA_API_BASE_URL for the endpoint", {
  eff <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "k-b", LLM_API_BASE_URL = "",
      CASSIA_API_BASE_URL = "https://legacy.example/chat/completions"),
    Triage:::.triage_resolve_api(NULL, NULL))
  expect_identical(eff$key, "k-b")
  expect_identical(eff$url, "https://legacy.example/chat/completions")
})

test_that("run_triage restores environment variables on preflight failure", {
  skip_if_runner_scripts_unavailable()
  deg <- triage_deg_fixture()
  before_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  before_url <- Sys.getenv("LLM_API_BASE_URL", unset = "")
  before_home <- Sys.getenv("TRIAGE_HOME", unset = "")
  err <- tryCatch(
    run_triage(deg = deg, species = "human", tissue = "pancreas",
               api_key = "secret-key-value",
               api_base_url = "https://example.org/chat/completions",
               preflight_only = TRUE,
               out = tempfile("triage_out_"),
               resource_root = tempfile("triage_empty_")),
    error = function(e) e)
  expect_true(inherits(err, "error"))
  expect_identical(Sys.getenv("DEEPSEEK_API_KEY", unset = ""), before_key)
  expect_identical(Sys.getenv("LLM_API_BASE_URL", unset = ""), before_url)
  expect_identical(Sys.getenv("TRIAGE_HOME", unset = ""), before_home)
})

test_that("run_triage accepts a custom key variable (LLM_API_KEY_ENV) and restores the environment on success", {
  skip_if_runner_scripts_unavailable()
  root <- triage_staged_resource_root()
  deg <- triage_deg_fixture()
  before_key <- Sys.getenv("DEEPSEEK_API_KEY", unset = "")
  before_url <- Sys.getenv("LLM_API_BASE_URL", unset = "")
  res <- withr::with_envvar(
    c(LLM_API_KEY_ENV = "MY_API_KEY", MY_API_KEY = "k-case-a",
      DEEPSEEK_API_KEY = "", LLM_API_BASE_URL = "", CASSIA_API_BASE_URL = ""),
    run_triage(deg = deg, species = "human", tissue = "pancreas",
               api_base_url = "https://example.org/chat/completions",
               preflight_only = TRUE, resource_root = root,
               out = tempfile("triage_out_")))
  expect_true(isTRUE(res$preflight$ok))
  expect_false(file.exists(file.path(res$out_root, "09_judge_outputs")))
  expect_identical(Sys.getenv("DEEPSEEK_API_KEY", unset = ""), before_key)
  expect_identical(Sys.getenv("LLM_API_BASE_URL", unset = ""), before_url)
})

test_that("run_triage accepts CASSIA_API_BASE_URL as endpoint fallback", {
  skip_if_runner_scripts_unavailable()
  root <- triage_staged_resource_root()
  deg <- triage_deg_fixture()
  res <- withr::with_envvar(
    c(CASSIA_API_BASE_URL = "https://legacy.example/chat/completions",
      LLM_API_BASE_URL = "", DEEPSEEK_API_KEY = "k-case-b",
      LLM_API_KEY_ENV = "DEEPSEEK_API_KEY"),
    run_triage(deg = deg, species = "human", tissue = "pancreas",
               preflight_only = TRUE, resource_root = root,
               out = tempfile("triage_out_")))
  expect_true(isTRUE(res$preflight$ok))
})

test_that("stage 10 extracts confidence_primary when only that field exists", {
  skip_if_runner_scripts_unavailable()
  script <- test_path("..", "..", "reproducibility", "scripts", "pipeline",
                      "10_judge_post_summary.R")
  if (!file.exists(script)) {
    script <- system.file("workflow", "pipeline", "10_judge_post_summary.R",
                          package = "Triage")
  }
  skip_if_not(file.exists(script), "stage-10 script unavailable")
  root <- tempfile("judge_out_")
  dir.create(file.path(root, "final"), recursive = TRUE, showWarnings = FALSE)
  fin <- list(
    cluster_id = "cluster_1",
    final_decision = list(primary_cell_type = "classical monocyte",
                          final_cell_ontology_id = "CL:0000860",
                          confidence_primary = 0.73),
    method_verdict = list(),
    post_issues = list(),
    evidence = list())
  jsonlite::write_json(fin,
                       file.path(root, "final",
                                 "cluster_1_LLM_JUDGE_FINAL.json"),
                       auto_unbox = TRUE)
  run_dir <- tempfile("stage10_cwd_")
  dir.create(run_dir)
  status <- system2(
    Triage:::.triage_rscript(),
    shQuote(c(script, "--out_root", root, "--dataset_name", "census_immune")),
    stdout = "", stderr = "",
    env = c(paste0("TRIAGE_HOME=", tempdir()),
            paste0("PROJECT_ROOT=", tempdir())))
  expect_identical(as.integer(status), 0L)
  summary_path <- file.path(root, "judge_post_summary", "summary.csv")
  expect_true(file.exists(summary_path))
  s <- read.csv(summary_path)
  expect_true("confidence" %in% names(s))
  expect_equal(as.numeric(s$confidence[s$cluster_id == "cluster_1"]), 0.73)
})

test_that("the resolved Rscript runs scripts and arguments containing spaces", {
  wd <- file.path(tempdir(), "dir with spaces")
  dir.create(wd, showWarnings = FALSE)
  script <- file.path(wd, "child run me.R")
  writeLines("cat('GOT:', commandArgs(trailingOnly = TRUE)[1], '\\n')", script)
  out <- system2(
    Triage:::.triage_rscript(),
    shQuote(c(script, "value with spaces")),
    stdout = TRUE, stderr = "")
  expect_match(paste(out, collapse = "\n"), "GOT: value with spaces")
})

test_that("CASSIA Python readiness uses a version-compatible lookup", {
  skip_if_not(requireNamespace("CASSIA", quietly = TRUE))
  # fresh-install provenance of the tested revision must be intact
  pd <- utils::packageDescription("CASSIA")
  sha <- if (!is.null(pd$RemoteSha) && nzchar(pd$RemoteSha)) pd$RemoteSha else pd$GithubSHA1
  skip_if_not(identical(sha, "b008c0ac3dd81b2c2dff131d20f5081a58aca027"),
              "CASSIA not installed from the pinned tested revision")
  # getFromNamespace works regardless of export status at that revision
  fun <- tryCatch(utils::getFromNamespace("check_python_env", "CASSIA"),
                  error = function(e) NULL)
  expect_false(is.null(fun))
  res <- tryCatch(isTRUE(fun()), error = function(e) NA)
  expect_true(is.logical(res) && !is.na(res))
})
