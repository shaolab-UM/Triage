# Integration hotfix tests: canonical paths, generic context propagation,
# cross-stage DEG handoff, and manifest absolute-run-dir wiring.
library(Triage)

script_path <- function(...) test_path("..", "..", "reproducibility", "scripts", ...)
skip_if_runner_scripts_unavailable <- function() {
  testthat::skip_if_not(file.exists(script_path("run_triage.R")),
                        "Repository-level runner scripts are excluded from the R package bundle.")
}

# --- tiny valid resource root helpers (mirror test-review-fixes.R) ---------
triage_staged_resource_root_species <- function(species = "human") {
  rr <- tempfile("triage_res_")
  dir.create(file.path(rr, "ppi"), recursive = TRUE)
  dir.create(file.path(rr, "collectri"), recursive = TRUE)
  dir.create(file.path(rr, "cellmarker"), recursive = TRUE)
  code <- if (species == "human") "9606" else "10090"
  writeLines(paste0("gene1\talias1\tsource1"), file.path(rr, "ppi",
    paste0(code, ".protein.aliases.v12.0.txt")))
  writeLines("p1\tp2\t900", file.path(rr, "ppi",
    paste0(code, ".protein.physical.links.v12.0.txt")))
  saveRDS(data.frame(source = "TF1", target = "G1", weight = 1),
          file.path(rr, "collectri", paste0("collectri_", species, "_network.rds")))
  cm <- if (species == "human") "Cell_marker_Human.xlsx" else "Cell_marker_Mouse.xlsx"
  writexl::write_xlsx(data.frame(cell_name = "monocyte", marker = "LST1"),
                      file.path(rr, "cellmarker", cm))
  rr
}

triage_deg_fixture <- function(path, clusters = c("Alpha", "Beta")) {
  readr::write_csv(
    data.frame(cluster = rep(clusters, each = 2),
               gene = c("G1", "G2", "G3", "G4"),
               avg_log2FC = c(1, 0.5, 0.8, 0.3),
               p_val = c(0.01, 0.02, 0.03, 0.04),
               p_val_adj = c(0.05, 0.1, 0.15, 0.2),
               pct.1 = c(0.5, 0.4, 0.6, 0.3),
               pct.2 = c(0.1, 0.2, 0.15, 0.25)),
    path)
  path
}

# --- A: generic DEG -> 01b -> 03a uses only anonymous `cluster` ------------
test_that("01b emits a `cluster` column of anonymous ids and 03a accepts it", {
  skip_if_runner_scripts_unavailable()
  run_dir <- tempfile("triage_int_")
  dir.create(run_dir)
  deg <- triage_deg_fixture(file.path(run_dir, "markers.csv"))
  map_dir <- file.path(run_dir, "01_input_prep")
  dir.create(map_dir)
  rscript <- file.path(R.home("bin"), "Rscript")
  prep <- system2(rscript,
                  shQuote(c(script_path("pipeline", "01b_prepare_user_deg.R"),
                            "--deg", deg, "--out-dir", map_dir)),
                  stdout = "", stderr = "")
  expect_identical(as.integer(prep), 0L)
  maskdeg <- readr::read_csv(file.path(map_dir, "maskdeg.csv"),
                             show_col_types = FALSE)
  expect_identical(names(maskdeg)[1], "cluster")
  expect_false("cluster_anon" %in% names(maskdeg))
  expect_true(all(maskdeg$cluster %in% c("cluster_1", "cluster_2")))
  expect_false(any(maskdeg$cluster %in% c("Alpha", "Beta")))
  expect_false(file.exists(file.path(map_dir, "true_label.csv")))
  # 03a consumes the 01b output unchanged
  filtered <- system2(rscript,
                      shQuote(c(script_path("pipeline", "03a_filter_deg.R"),
                                "--deg", file.path(map_dir, "maskdeg.csv"),
                                "--out_dir", run_dir)),
                      stdout = "", stderr = "")
  expect_identical(as.integer(filtered), 0L)
  fdeg <- readr::read_csv(file.path(run_dir, "filtered_deg.csv"),
                          show_col_types = FALSE)
  expect_true(nrow(fdeg) > 0)
  expect_true(all(fdeg$cluster %in% c("cluster_1", "cluster_2")))
  unlink(run_dir, recursive = TRUE)
})

# --- B: spaces in cwd + default out = "results" resolve absolutely ---------
test_that("run_triage resolves default out from a working directory with spaces", {
  space_dir <- tempfile("triage space ")
  dir.create(space_dir)
  deg <- triage_deg_fixture(file.path(space_dir, "markers.csv"))
  rr <- triage_staged_resource_root_species("human")
  res <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "test-key",
      LLM_API_BASE_URL = "https://example.test/chat/completions"),
    withr::with_dir(space_dir, {
      run_triage(deg = deg, species = "human", tissue = "pancreas",
                 out = "results", dataset_name = "user_dataset",
                 preflight_only = TRUE, resource_root = rr)
    }))
  expect_true(grepl("^/", res$out_root))
  expect_identical(normalizePath(res$out_root, mustWork = FALSE),
                   normalizePath(file.path(space_dir, "results",
                                           "user_dataset", res$run_tag),
                                 mustWork = FALSE))
  # generic context file written and wired for the run
  ctx <- file.path(res$runtime_home, "dataset_context.json")
  expect_true(file.exists(ctx))
  ctx_obj <- jsonlite::fromJSON(ctx, simplifyVector = FALSE)
  expect_identical(ctx_obj$tissue[[1]], "pancreas")
  expect_identical(ctx_obj$dataset_scope[[1]], "mixed_unknown")
  expect_identical(ctx_obj$gate_mode[[1]], "flag_only")
  unlink(space_dir, recursive = TRUE)
  unlink(rr, recursive = TRUE)
})

# --- C: generic context overrides any colliding dataset name ---------------
test_that("context JSON wins over manuscript profiles for colliding names", {
  ctx <- tempfile("ctx_", fileext = ".json")
  jsonlite::write_json(list(
    species = "mouse", tissue = "hippocampus", study_context = "reference",
    dataset_scope = "mixed_unknown", scope_profile = "mixed_unknown",
    gate_mode = "flag_only", allowed_lineages = character(0)),
    ctx, auto_unbox = TRUE)
  cfg <- withr::with_envvar(c(DATASET_CONTEXT_JSON_PATH = ctx), {
    Triage:::get_dataset_config("census_immune", project_root = tempdir())
  })
  expect_identical(cfg$species, "mouse")
  expect_identical(cfg$tissue, "hippocampus")
  expect_identical(cfg$scope_profile, "mixed_unknown")
  expect_identical(cfg$gate_mode, "flag_only")
  unlink(ctx)
  # without the override, the released benchmark context is untouched
  cfg0 <- Triage:::get_dataset_config("census_immune", project_root = tempdir())
  expect_identical(cfg0$tissue, "blood")
  # mouse preflight resolves the 10090 naming contract
  rr_mouse <- triage_staged_resource_root_species("mouse")
  pf <- withr::with_envvar(
    c(DEEPSEEK_API_KEY = "test-key",
      LLM_API_BASE_URL = "https://example.test/chat/completions"),
    triage_preflight(species = "mouse", resource_root = rr_mouse))
  expect_true(pf$ok)
  unlink(rr_mouse, recursive = TRUE)
})

# --- D: stage-05 intermediate root is honored by the real resolver ---------
test_that("stage-05 writes enrichment TSVs under the requested intermediate root", {
  skip_if_runner_scripts_unavailable()
  llm_run <- script_path("llm_run.R")
  # Execute the REAL resolver and writer from llm_run.R in a subprocess:
  # parse the actual script, evaluate resolve_bioinfo_dir +
  # save_enrichment_tsv from it, write through save_enrichment_tsv, and
  # verify the file lands under the explicit TRIAGE_INTERMEDIATE_ROOT.
  driver <- tempfile("drv_", fileext = ".R")
  ok_file <- tempfile("drvok_", fileext = ".txt")
  root <- tempfile("tsv_root_")
  dir.create(root)
  tmpl <- '
    Sys.setenv(TRIAGE_INTERMEDIATE_ROOT = %1$s)
    src <- parse(%2$s)
    fns <- Filter(function(e) is.call(e) && length(e) == 3 &&
      e[[1]] == quote("<-") && identical(e[[2]], quote(save_enrichment_tsv)), src)
    resv <- Filter(function(e) is.call(e) && length(e) == 3 &&
      e[[1]] == quote("<-") && identical(e[[2]], quote(resolve_bioinfo_dir)), src)
    stopifnot(length(fns) == 1, length(resv) == 1)
    save_enrichment_tsv <- eval(fns[[1]])
    resolve_bioinfo_dir <- eval(resv[[1]])
    d <- resolve_bioinfo_dir("user_dataset_LLM_Input_Run", "cluster_1")
    stopifnot(identical(d, file.path(%1$s, "user_dataset_LLM_Input_Run",
                                   "bioinformatics_tsv", "cluster_1")))
    save_enrichment_tsv(data.frame(ID = "x", p = 0.01),
                        "go_biological_process", d)
    f <- file.path(d, "go_biological_process_full_results.tsv")
    stopifnot(file.exists(f))
    # without the env var the historical cwd-relative default is preserved
    Sys.unsetenv("TRIAGE_INTERMEDIATE_ROOT")
    stopifnot(identical(resolve_bioinfo_dir("user_dataset_LLM_Input_Run", "cluster_1"),
                        file.path("intermediate_outputs", "user_dataset_LLM_Input_Run",
                                  "bioinformatics_tsv", "cluster_1")))
    writeLines("SUBPROCESS_OK", %3$s)
  '
  writeLines(sprintf(tmpl, dQuote(root), dQuote(llm_run), dQuote(ok_file)), driver)
  rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(rscript, shQuote(c("--vanilla", driver)),
                    stdout = "", stderr = "")
  expect_identical(as.integer(status), 0L)
  expect_true(file.exists(ok_file) &&
              grepl("SUBPROCESS_OK", readLines(ok_file, warn = FALSE)[1]))
  # both orchestrators pass an explicit absolute intermediate root
  for (runner in c(script_path("run_triage.R"),
                   file.path(system.file("workflow", package = "Triage"),
                             "..", "..", "R", "triage_workflow.R"))) {
    if (!file.exists(runner)) next
    txt <- paste(readLines(runner, warn = FALSE), collapse = "\n")
    expect_true(grepl('"--intermediate_root", intermediate_root', txt,
                      fixed = TRUE), info = runner)
    expect_identical(
      length(gregexpr("paste0(dataset_name, \"_LLM_Input_Run\")", txt,
                      fixed = TRUE)[[1]]), 1L,
      info = runner)
    expect_true(grepl('intermediate_run_dir, "bioinformatics_tsv"', txt,
                      fixed = TRUE), info = runner)
    expect_true(grepl('"--intermediate_outputs_dir", intermediate_run_dir',
                      txt, fixed = TRUE), info = runner)
  }
  # stage 05 declares the option and wires the env contract
  s05 <- paste(readLines(script_path("pipeline", "05_prepare_llm_inputs.R"),
                         warn = FALSE), collapse = "\n")
  expect_true(grepl("TRIAGE_INTERMEDIATE_ROOT", s05, fixed = TRUE))
  unlink(c(root, driver), recursive = TRUE)
})

# --- E: run manifests carry an absolute run_dir -----------------------------
test_that("run manifests are written from an absolute run_dir", {
  skip_if_runner_scripts_unavailable()
  for (f in c("run_triage.R", file.path("pipeline", "05_prepare_llm_inputs.R"))) {
    txt <- paste(readLines(script_path(f), warn = FALSE), collapse = "\n")
    if (f == "run_triage.R") {
      expect_true(grepl("out_root <- normalizePath", txt, fixed = TRUE))
      expect_true(grepl("run_dir <- out_root", txt, fixed = TRUE))
    }
  }
  # absolute manifest entries resolve through the packaged helper unchanged
  absolute <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  expect_identical(
    Triage:::resolve_manifest_run_dir(tempdir(), "ds", absolute), absolute)
})

# --- F: pinned GOSemSim + clusterProfiler provenance contract ---------------
test_that("tested-revision provenance helpers and pins are wired everywhere", {
  expect_true(Triage:::.triage_clusterprofiler_sha_matches())
  expect_false(
    Triage:::.triage_clusterprofiler_sha_matches("deadbeef"))
  expect_identical(Triage:::.triage_tested_clusterprofiler_sha(),
                   "f9f0d502508cacd258ac1a1cba6d5d497b98fe6c")
  expect_identical(Triage:::.triage_tested_gosemsim_sha(),
                   "67e3da1dd3ee9d7c5067b2044fcf979e0cf6480d")
  # both preflights and stage 06b reference the tested pins
  skip_if_runner_scripts_unavailable()
  pf_txt <- paste(readLines(script_path("preflight_check.R"), warn = FALSE),
                  collapse = "\n")
  expect_true(grepl(".triage_tested_clusterprofiler_sha()", pf_txt,
                    fixed = TRUE))
  expect_true(grepl(".triage_tested_gosemsim_sha()", pf_txt,
                    fixed = TRUE))
  for (f in c("pipeline/06b_run_inter.R",
              file.path("pipeline", "06b_run_inter.R"))) {
    txt <- paste(readLines(script_path(f), warn = FALSE), collapse = "\n")
    expect_true(grepl("f9f0d502508cacd258ac1a1cba6d5d497b98fe6c", txt,
                      fixed = TRUE), info = f)
  }
  wf_txt <- paste(readLines(system.file("workflow", "README.md",
                                        package = "Triage"), warn = FALSE),
                  collapse = "\n")
  # installer carries both pins (source-level, from the repository file)
  ih_txt <- paste(readLines(test_path("..", "..", "R", "install_helpers.R"),
                            warn = FALSE), collapse = "\n")
  expect_true(grepl("67e3da1dd3ee9d7c5067b2044fcf979e0cf6480d", ih_txt,
                    fixed = TRUE))
  expect_true(grepl("remotes::install_github(\"YuLab-SMU/GOSemSim\"", ih_txt,
                    fixed = TRUE))
  expect_true(grepl("remotes::install_github(\"YuLab-SMU/clusterProfiler\"", ih_txt,
                    fixed = TRUE))
})
