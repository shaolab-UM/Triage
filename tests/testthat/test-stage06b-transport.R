# Stage-06b LLM transport routing (no network): model/endpoint reach the
# transport, the fanyi endpoint shim is self-restoring, and the
# intermediate root survives a real multisession worker.
library(Triage)

llm_run_path <- test_path("..", "..", "reproducibility", "scripts", "llm_run.R")

# --- (a) model/api_key/task routing reaches clusterProfiler::interpret -----
test_that("interpret routing passes model, api_key and task through", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    interpret = function(enrich_list, context = NULL, model = NULL,
                         api_key = NULL, task = "interpretation", ...) {
      captured$model <- model
      captured$api_key <- api_key
      captured$task <- task
      "interpret-called"
    },
    .package = "clusterProfiler")
  out <- Triage:::.triage_interpret(
    enrich_list = list(),
    context = "ctx",
    model = "deepseek-v4-flash",
    api_key = "test-key",
    base_url = Triage:::.triage_fanyi_default_endpoint,
    task = "annotation")
  expect_identical(out, "interpret-called")
  expect_identical(captured$model, "deepseek-v4-flash")
  expect_identical(captured$api_key, "test-key")
  expect_identical(captured$task, "annotation")
})

# --- (b) .triage_fanyi_transport posts to the supplied endpoint ------------
test_that("fanyi transport posts to the supplied endpoint URL", {
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    req_perform = function(req) {
      captured$url <- req$url
      httr2::response(status_code = 200,
                      headers = list("content-type" = "application/json"),
                      body = charToRaw('{"choices":[{"message":{"content":"ok"}}]}'))
    },
    .package = "httr2")
  res <- Triage:::.triage_fanyi_transport(
    list(list(role = "user", content = "hi")),
    model = "deepseek-v4-flash",
    api_key = "test-key",
    base_url = "https://custom.example/chat/completions")
  expect_identical(captured$url, "https://custom.example/chat/completions")
  expect_identical(httr2::resp_status(res), 200L)
})

# --- (c) shim routes chat_request and restores the fanyi binding -----------
test_that("endpoint shim routes fanyi and restores the original binding", {
  binding_name <- ".deepseek_query_messages"
  ns <- asNamespace("fanyi")
  orig <- get(binding_name, envir = ns)
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    req_perform = function(req) {
      captured$url <- req$url
      httr2::response(status_code = 200,
                      headers = list("content-type" = "application/json"),
                      body = charToRaw('{"choices":[{"message":{"content":"ok"}}]}'))
    },
    .package = "httr2")
  # mock clusterProfiler::interpret so it exercises fanyi::chat_request,
  # exactly the call chain interpret() uses internally
  testthat::local_mocked_bindings(
    interpret = function(enrich_list, context = NULL, model = NULL,
                         api_key = NULL, task = "interpretation", ...) {
      fanyi::chat_request("hi", model = model, api_key = api_key)
    },
    .package = "clusterProfiler")
  out <- Triage:::.triage_interpret(
    enrich_list = list(),
    model = "deepseek-v4-flash",
    api_key = "test-key",
    base_url = "https://custom.example/chat/completions",
    task = "annotation")
  expect_identical(captured$url, "https://custom.example/chat/completions")
  # fanyi::chat_request returns the parsed message content
  expect_identical(out, "ok")
  # binding restored: identical function, locked again
  expect_identical(get(binding_name, envir = ns), orig)
  expect_true(environmentIsLocked(asNamespace("fanyi")))
})

# --- (d) tested fanyi version contract --------------------------------------
test_that("fanyi version helpers expose the tested version", {
  expect_identical(Triage:::.triage_tested_fanyi_version(), "0.1.0")
  expect_identical(Triage:::.triage_fanyi_default_endpoint,
                   "https://api.deepseek.com/v1/chat/completions")
  expect_true(Triage:::.triage_fanyi_version_matches())
  expect_false(Triage:::.triage_fanyi_version_matches("9.9.9"))
})

# --- (e) multisession worker writes under the intermediate root ------------
test_that("multisession workers honor the intermediate root", {
  testthat::skip_if_not(file.exists(llm_run_path),
                        "Repository-level workflow scripts are excluded from the R package bundle.")
  testthat::skip_if_not(requireNamespace("future.apply", quietly = TRUE))
  # pull the real resolver + writer function definitions out of llm_run.R
  exprs <- parse(text = paste(readLines(llm_run_path, warn = FALSE),
                              collapse = "\n"))
  wanted <- c("resolve_bioinfo_dir", "save_enrichment_tsv")
  picked <- Filter(function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && length(e) >= 3 &&
      is.symbol(e[[2]]) && as.character(e[[2]]) %in% wanted
  }, exprs)
  expect_length(picked, 2L)
  code_txt <- paste(vapply(picked, function(e) paste(deparse(e), collapse = "\n"), character(1)), collapse = "\n")
  root <- tempfile("ms_root_")
  plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(plan), add = TRUE)
  # the root is set AFTER plan(multisession): the worker must still see it
  # because the parent exports it explicitly as a future global.
  fns_env <- new.env(parent = baseenv())
  eval(parse(text = code_txt), envir = fns_env)
  out <- withr::with_envvar(c(TRIAGE_INTERMEDIATE_ROOT = root), {
    future.apply::future_lapply(
    "cluster_1",
    function(qid) {
      d <- resolve_bioinfo_dir("user_dataset_LLM_Input_Run", qid,
                               intermediate_root = triage_intermediate_root)
      save_enrichment_tsv(data.frame(ID = "x", p = 0.01),
                          "go_biological_process", d)
      d
    },
    future.globals = list(
      resolve_bioinfo_dir = fns_env$resolve_bioinfo_dir,
      save_enrichment_tsv = fns_env$save_enrichment_tsv,
      triage_intermediate_root = Sys.getenv("TRIAGE_INTERMEDIATE_ROOT",
                                            unset = ""))
  )
  })
  expected <- file.path(root, "user_dataset_LLM_Input_Run",
                        "bioinformatics_tsv", "cluster_1",
                        "go_biological_process_full_results.tsv")
  expect_true(file.exists(expected))
  # nothing written relative to the test working directory
  expect_false(file.exists(file.path(getwd(), "intermediate_outputs",
                                     "user_dataset_LLM_Input_Run")))
})
