# Triage package tests (self-attaching so test_dir works standalone)
library(Triage)

# Installed-package full-workflow runtime: preflight, resource contract and
# byte-identity of the bundled pipeline scripts. Deterministic; no API calls.

test_that("bundled workflow runtime ships and matches the repository scripts", {
  wf <- system.file("workflow", package = "Triage")
  expect_true(nzchar(wf) && dir.exists(file.path(wf, "pipeline")))
  required <- c("01b_prepare_user_deg.R", "03a_filter_deg.R",
                "03b_run_cassia.R", "04a_build_candidates.R",
                "05_prepare_llm_inputs.R", "06_run_llm_pipeline.R",
                "06b_run_inter.R", "07_our_llm_summary.R",
                "07b_inter_summary.R", "07.5_build_reviewer_mapping.R",
                "08_build_judge_inputs.R", "09_run_judge.R",
                "10_judge_post_summary.R", "11_eval_accuracy.R")
  for (f in required) {
    expect_true(file.exists(file.path(wf, "pipeline", f)), info = f)
  }
  expect_true(file.exists(file.path(wf, "llm_run.R")))

  repo <- test_path("..", "..", "reproducibility")
  testthat::skip_if_not(dir.exists(repo),
                        "Repository scripts are excluded from the R package bundle.")
  for (f in c(required, file.path("..", "..", "llm_run.R"))) {
    inst_path <- if (grepl("llm_run", f)) file.path(wf, "llm_run.R") else
      file.path(wf, "pipeline", basename(f))
    repo_path <- if (grepl("llm_run", f)) file.path(repo, "scripts", "llm_run.R") else
      file.path(repo, "scripts", "pipeline", basename(f))
    expect_identical(
      tools::md5sum(inst_path)[[1]], tools::md5sum(repo_path)[[1]],
      info = paste0("byte-identity: ", basename(f)))
  }
})

test_that("triage_preflight reports all missing prerequisites together", {
  rr <- tempfile("triage_res_")
  dir.create(file.path(rr, "ppi"), recursive = TRUE)
  pf <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "", LLM_API_BASE_URL = "", LLM_API_KEY_ENV = "DEEPSEEK_API_KEY"),
    triage_preflight(species = "human", resource_root = rr))
  expect_false(pf$ok)
  expect_true(any(grepl("API key", pf$missing)))
  expect_true(any(grepl("LLM endpoint", pf$missing)))
  expect_true(any(grepl("STRING aliases", pf$missing)))
  expect_true(any(grepl("STRING physical links", pf$missing)))
  expect_true(any(grepl("CollecTRI", pf$missing)))
  expect_true(any(grepl("CellMarkerDB", pf$missing)))
})

test_that("triage_preflight passes a structurally valid resource root", {
  rr <- tempfile("triage_res_ok_")
  dir.create(file.path(rr, "ppi"), recursive = TRUE)
  dir.create(file.path(rr, "collectri"), recursive = TRUE)
  dir.create(file.path(rr, "cellmarker"), recursive = TRUE)
  writeLines(c("ensp1\tsymbol1\t9606", "ensp2\tsymbol2\t9606"),
             file.path(rr, "ppi", "9606.protein.aliases.v12.0.txt"))
  writeLines(c("ensp1\tensp2\t900", "ensp2\tensp3\t700"),
             file.path(rr, "ppi", "9606.protein.physical.links.v12.0.txt"))
  net <- data.frame(source = c("TF1", "TF2"), target = c("G1", "G2"),
                    weight = c(2, 1))
  saveRDS(net, file.path(rr, "collectri", "collectri_human_network.rds"))
  df <- data.frame(cell_name = "monocyte", marker = "LST1")
  writexl::write_xlsx(df, file.path(rr, "cellmarker", "Cell_marker_Human.xlsx"))
  pf <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "test-key", LLM_API_KEY_ENV = "DEEPSEEK_API_KEY",
      LLM_API_BASE_URL = "https://example.test/chat/completions"),
    triage_preflight(species = "human", resource_root = rr))
  expect_true(pf$ok)
  expect_length(pf$missing, 0L)
})

test_that("run_triage validates arguments before any staging", {
  expect_error(run_triage(deg = "no_such_file.csv", species = "human",
                          tissue = "pancreas"),
               regexp = "deg file not found")
  deg <- file.path(tempdir(), "bad_deg.csv")
  writeLines("wrong,columns\n1,2", deg)
  expect_error(run_triage(deg = deg, species = "human", tissue = NULL),
               regexp = "tissue is required")
  expect_error(run_triage(deg = deg, species = "drosophila", tissue = "x"),
               regexp = "human")
})

test_that("run_triage preflight_only stages resources and resolves installed paths", {
  wf <- system.file("workflow", package = "Triage")
  testthat::skip_if_not(nzchar(wf), "workflow runtime not bundled")

  # fake but structurally valid resource root + DEG
  rr <- tempfile("triage_res_run_")
  dir.create(file.path(rr, "ppi"), recursive = TRUE)
  dir.create(file.path(rr, "collectri"), recursive = TRUE)
  dir.create(file.path(rr, "cellmarker"), recursive = TRUE)
  writeLines(c("ensp1\tsymbol1\t9606"),
             file.path(rr, "ppi", "9606.protein.aliases.v12.0.txt"))
  writeLines(c("ensp1\tensp2\t900"),
             file.path(rr, "ppi", "9606.protein.physical.links.v12.0.txt"))
  saveRDS(data.frame(source = "TF1", target = "G1", weight = 2),
          file.path(rr, "collectri", "collectri_human_network.rds"))
  writexl::write_xlsx(data.frame(cell_name = "monocyte", marker = "LST1"),
                      file.path(rr, "cellmarker", "Cell_marker_Human.xlsx"))
  out_dir <- tempfile("triage_run_")
  dir.create(out_dir)
  deg <- file.path(out_dir, "markers.csv")
  readr::write_csv(
    data.frame(cluster = c("Alpha", "Beta"), gene = c("G1", "G2"),
               avg_log2FC = c(1, 0.5), p_val = c(0.01, 0.02),
               p_val_adj = c(0.05, 0.1), pct.1 = c(0.5, 0.4),
               pct.2 = c(0.1, 0.2)),
    deg)

  # environment contract: variables are set for the child stages and
  # RESTORED to their pre-call values after run_triage() returns
  env_pre <- vapply(c("TRIAGE_HOME", "PROJECT_ROOT", "CL_LOCAL_JSON",
                      "TRIAGE_PPI_ROOT", "TRIAGE_WORKFLOW_DIR",
                      "DATASET_CONTEXT_JSON_PATH"),
                    function(v) Sys.getenv(v, unset = ""), character(1))
  res <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "test-key", LLM_API_KEY_ENV = "DEEPSEEK_API_KEY",
      LLM_API_BASE_URL = "https://example.test/chat/completions"),
    run_triage(deg = deg, species = "human", tissue = "pancreas",
               api_key = NULL, api_base_url = NULL, out = out_dir,
               dataset_name = "user_dataset", preflight_only = TRUE,
               resource_root = rr))
  expect_true(res$preflight$ok)
  expect_false(dir.exists(file.path(res$out_root, "01_input_prep")))

  # setup output paths are exactly the paths consumed by run_triage:
  # collectri + cellmarker staged from the resource root into runtime home
  staged_coll <- file.path(res$runtime_home, "inputs", "raw", "collectri",
                           "collectri_human_network.rds")
  staged_cm <- file.path(res$runtime_home, "inputs", "raw", "cellmarker",
                         "Cell_marker_Human.xlsx")
  expect_true(file.exists(staged_coll))
  expect_true(file.exists(staged_cm))
  expect_identical(tools::md5sum(staged_coll)[[1]],
                   tools::md5sum(file.path(rr, "collectri",
                                           "collectri_human_network.rds"))[[1]])
  expect_identical(tools::md5sum(staged_cm)[[1]],
                   tools::md5sum(file.path(rr, "cellmarker",
                                           "Cell_marker_Human.xlsx"))[[1]])
  # bundled ontology staged; workflow dir is the installed package location
  expect_true(file.exists(file.path(res$runtime_home, "inputs", "raw",
                                    "ontology", "CL-ontology-v2025-07-30.json")))
  expect_identical(normalizePath(res$runtime_home, mustWork = FALSE),
                   normalizePath(file.path(res$out_root, "runtime_home"),
                                 mustWork = FALSE))
  # environment contract: restored to pre-call values after the call
  env_post <- vapply(c("TRIAGE_HOME", "PROJECT_ROOT", "CL_LOCAL_JSON",
                       "TRIAGE_PPI_ROOT", "TRIAGE_WORKFLOW_DIR",
                       "DATASET_CONTEXT_JSON_PATH"),
                     function(v) Sys.getenv(v, unset = ""), character(1))
  expect_identical(env_post, env_pre)
  # no stage outputs beyond preflight
  expect_false(file.exists(file.path(res$out_root, "03b_cassia")))
})
