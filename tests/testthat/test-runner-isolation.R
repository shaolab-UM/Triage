# Integration checks for the user-input preparation and preflight scripts.
# They only execute deterministic setup code and never run an LLM stage.

repo_root <- normalizePath(test_path("..", ".."), winslash = "/")
script_path <- function(...) file.path(repo_root, "reproducibility", "scripts", ...)
skip_if_runner_scripts_unavailable <- function() {
  testthat::skip_if_not(
    file.exists(script_path("run_triage.R")),
    "Repository-level runner scripts are excluded from the R package bundle."
  )
}
run_rscript <- function(script, args, env = character()) {
  out <- tempfile("triage-runner-test-", fileext = ".log")
  on.exit(unlink(out), add = TRUE)
  status <- system2("Rscript", c(script, args), stdout = out, stderr = out,
                    env = env)
  list(status = status, output = paste(readLines(out, warn = FALSE), collapse = "\n"))
}

test_that("user preparation anonymizes biological-looking cluster names", {
  skip_if_runner_scripts_unavailable()
  tmp <- tempfile("triage-user-prep-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  deg <- file.path(tmp, "markers.csv")
  readr::write_csv(data.frame(
    cluster = rep(c("Malignant epithelial", "Activated T cell"), each = 2),
    gene = c("EPCAM", "KRT19", "CD3D", "TRAC"),
    avg_log2FC = c(2, 1, 3, 2), p_val = 0, p_val_adj = 0,
    pct.1 = .5, pct.2 = .1
  ), deg)

  out <- file.path(tmp, "prepared")
  res <- run_rscript(script_path("pipeline", "01b_prepare_user_deg.R"),
                     c("--deg", deg, "--out-dir", out))
  expect_equal(res$status, 0L, info = res$output)

  map <- readr::read_csv(file.path(out, "cluster_map.csv"), show_col_types = FALSE)
  masked <- readr::read_csv(file.path(out, "maskdeg.csv"), show_col_types = FALSE)
  expect_equal(map$cluster_id, c("cluster_1", "cluster_2"))
  expect_true(all(c("Activated T cell", "Malignant epithelial") %in% map$cluster_label))
  expect_false("cluster" %in% names(masked))
  expect_equal(names(masked), c("cluster_anon", "gene", "avg_log2FC", "p_val",
                                "p_val_adj", "pct.1", "pct.2"))
  expect_true(all(masked$cluster_anon %in% c("cluster_1", "cluster_2")))
  expect_false(file.exists(file.path(out, "true_label.csv")))
  expect_false(any(grepl("Malignant|Activated", unlist(masked), fixed = FALSE)))
})

test_that("user preparation rejects malformed DEG input before pipeline stages", {
  skip_if_runner_scripts_unavailable()
  tmp <- tempfile("triage-user-prep-bad-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  deg <- file.path(tmp, "bad.csv")
  readr::write_csv(data.frame(cluster = "Tumor", gene = "EPCAM", avg_log2FC = 2), deg)
  res <- run_rscript(script_path("pipeline", "01b_prepare_user_deg.R"),
                     c("--deg", deg, "--out-dir", file.path(tmp, "prepared")))
  expect_false(identical(res$status, 0L))
  expect_match(res$output, "missing required column")
})

test_that("preflight rejects malformed CellMarkerDB and CollecTRI resources", {
  skip_if_runner_scripts_unavailable()
  tmp <- tempfile("triage-preflight-")
  home <- file.path(tmp, "home")
  project <- file.path(tmp, "project")
  ppi <- file.path(tmp, "ppi")
  dir.create(file.path(home, "inputs", "raw", "ontology"), recursive = TRUE)
  dir.create(file.path(home, "inputs", "raw", "cellmarker"), recursive = TRUE)
  dir.create(file.path(project, "inputs", "raw", "collectri"), recursive = TRUE)
  dir.create(ppi, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  file.copy(file.path(repo_root, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"),
            file.path(home, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"))
  writeLines("a b c", file.path(ppi, "9606.protein.aliases.v12.0.txt"))
  writeLines("a b c", file.path(ppi, "9606.protein.physical.links.v12.0.txt"))
  saveRDS(data.frame(source = "TF", target = "GENE"),
          file.path(project, "inputs", "raw", "collectri", "collectri_human_network.rds"))
  writexl::write_xlsx(data.frame(other = "value"),
                      file.path(home, "inputs", "raw", "cellmarker", "Cell_marker_Human.xlsx"))

  env <- c(
    paste0("TRIAGE_HOME=", home), paste0("PROJECT_ROOT=", project),
    paste0("TRIAGE_PPI_ROOT=", ppi), "DEEPSEEK_API_KEY=test-key",
    "LLM_API_BASE_URL=https://example.invalid/chat/completions"
  )
  res <- run_rscript(script_path("preflight_check.R"),
                     c("--species", "human", "--skip-internet"), env)
  expect_false(identical(res$status, 0L))
  expect_match(res$output, "CollecTRI RDS is invalid")
  expect_match(res$output, "CellMarkerDB spreadsheet is invalid")
})

test_that("benchmark preflight resolves context without reading evaluation files", {
  skip_if_runner_scripts_unavailable()
  res <- run_rscript(script_path("run_triage.R"),
                     c("--benchmark", "Census_immune", "--preflight-only",
                       "--out", tempfile("triage-benchmark-")))
  # Local resources are intentionally absent: preflight stops before stages.
  expect_false(identical(res$status, 0L))
  expect_match(res$output, "dataset = Census_immune")
  expect_match(res$output, "species = human")
  expect_match(res$output, "tissue = blood")
  expect_match(res$output, "study context = healthy_immune")
  expect_match(res$output, "model_inputs/maskdeg.csv")

  source_text <- paste(readLines(script_path("run_triage.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("reproducibility/primary/.*/evaluation|reference_cl", source_text))
  expect_match(source_text, "if \\(!is.null\\(opt\\$`reference-labels`\\)\\)")
  expect_match(source_text, "--true_label_csv", fixed = TRUE)
})
