# API configuration UX tests (no network access).
library(Triage)

test_that("api_base_url is forwarded to the provider (NULL keeps env fallback)", {
  jin <- read_triage_input(test_path("fixtures", "cluster_1_round1.json"))
  head_fixture <- jsonlite::fromJSON(test_path("fixtures", "head_round1.json"),
                                     simplifyVector = FALSE)
  fenced <- paste0("```json\n",
                   jsonlite::toJSON(head_fixture, auto_unbox = TRUE, null = "null"),
                   "\n```")

  calls <- list()
  stub <- function(prompt_json_string, api_key, model = "deepseek-v4-flash",
                   temperature = 0, timeout_seconds = 1200, max_retries = 4,
                   retry_delay = 4, system_prompt = NULL, base_url = NULL) {
    calls[[length(calls) + 1L]] <<- list(base_url = base_url)
    list(ok = TRUE, text = fenced, usage = NULL, status = 200L,
         model = model, error = NULL, retry_count = 0L, error_class = NULL)
  }
  testthat::local_mocked_bindings(invoke_deepseek_api = stub, .package = "Triage")

  fin <- run_triage_adjudication(jin, use_api = TRUE, api_key = "test-key",
                                 api_base_url = "https://cli-endpoint.example/chat/completions",
                                 dataset_name = "census_immune",
                                 project_root = test_path("fixtures"))
  expect_true(attr(fin, "gate")$ok)
  expect_length(calls, 1L)
  expect_equal(calls[[1]]$base_url, "https://cli-endpoint.example/chat/completions")

  # NULL keeps the existing environment-based resolution inside the provider.
  fin2 <- run_triage_adjudication(jin, use_api = TRUE, api_key = "test-key",
                                  dataset_name = "census_immune",
                                  project_root = test_path("fixtures"))
  expect_true(attr(fin2, "gate")$ok)
  expect_length(calls, 2L)
  expect_null(calls[[2]]$base_url)
})

test_that("NAMESPACE keeps exactly six exports", {
  ns_lines <- readLines(system.file("NAMESPACE", package = "Triage"))
  expect_equal(sum(grepl("^export\\(", ns_lines)), 6L)
})

test_that("CLI --api-key/--api-base-url override environment values and never echo the key", {
  runner <- test_path("..", "..", "reproducibility", "scripts", "run_triage.R")
  testthat::skip_if_not(file.exists(runner),
                        "Repository-level runner scripts are excluded from the R package bundle.")
  out_file <- tempfile(fileext = ".log")
  on.exit(unlink(out_file), add = TRUE)

  env <- c(DEEPSEEK_API_KEY = "env-secret-key-value",
           LLM_API_BASE_URL = "https://env-endpoint.example/chat/completions")

  # A. CLI values take precedence over environment values.
  status <- withr::with_envvar(env, {
    system2("Rscript", c(shQuote(runner), "--benchmark", "Census_immune",
                         "--preflight-only", "--out", tempfile(),
                         "--api-key", "cli-secret-key-value",
                         "--api-base-url", "https://cli-endpoint.example/chat/completions"),
            stdout = out_file, stderr = out_file)
  })
  log <- paste(readLines(out_file, warn = FALSE), collapse = "\n")
  expect_true(grepl("https://cli-endpoint.example/chat/completions", log, fixed = TRUE))
  expect_false(grepl("https://env-endpoint.example", log, fixed = TRUE))
  # D. No API key value is printed anywhere in normal output.
  expect_false(grepl("cli-secret-key-value", log, fixed = TRUE))
  expect_false(grepl("env-secret-key-value", log, fixed = TRUE))
  expect_true(status %in% c(0L, 1L))  # preflight outcome may vary by machine; URL assert above is the contract

  # B. Environment-only configuration still works (URL observed by preflight).
  status2 <- withr::with_envvar(env, {
    system2("Rscript", c(shQuote(runner), "--benchmark", "Census_immune",
                         "--preflight-only", "--out", tempfile()),
            stdout = out_file, stderr = out_file)
  })
  log2 <- paste(readLines(out_file, warn = FALSE), collapse = "\n")
  expect_true(grepl("https://env-endpoint.example/chat/completions", log2, fixed = TRUE))
  expect_false(grepl("env-secret-key-value", log2, fixed = TRUE))
  expect_true(status2 %in% c(0L, 1L))
})
