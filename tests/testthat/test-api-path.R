# Deterministic mocked test for the Handling Editor API path.
library(Triage)
# No network access: invoke_deepseek_api is replaced by a capturing stub.

test_that("API path sends the serialized Handling Editor query (not a nested request)", {
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
    calls[[length(calls) + 1L]] <<- list(
      prompt = prompt_json_string, model = model, temperature = temperature,
      system_prompt = system_prompt)
    list(ok = TRUE, text = fenced, usage = NULL, status = 200L,
         model = model, error = NULL, retry_count = 0L, error_class = NULL)
  }
  testthat::local_mocked_bindings(invoke_deepseek_api = stub, .package = "Triage")

  expect_equal(length(calls), 0L)
  fin <- run_triage_adjudication(jin, use_api = TRUE, api_key = "test-key",
                                 model = "deepseek-v4-flash", temperature = 0,
                                 dataset_name = "census_immune",
                                 project_root = test_path("fixtures"))

  # Exactly one API call (Handling Editor only; Chief QC is not invoked).
  expect_length(calls, 1L)

  # The wrapper received the serialized Handling Editor QUERY object, not a
  # nested chat-completions request.
  q <- jsonlite::fromJSON(calls[[1]]$prompt, simplifyVector = FALSE)
  expect_false("messages" %in% names(q))
  expect_false("model" %in% names(q))
  expect_false("stream" %in% names(q))
  expect_true(is.list(q) && length(q) > 0)

  # Primary model identifier and temperature.
  expect_equal(calls[[1]]$model, "deepseek-v4-flash")
  expect_equal(calls[[1]]$temperature, 0)

  # Handling Editor system prompt passed separately.
  sp <- calls[[1]]$system_prompt
  expect_true(is.character(sp) && nzchar(sp))
  expect_identical(sp, Triage:::build_head_editor_system_prompt_compact(species_value = "human"))

  # The returned $text field was parsed into the adjudication record.
  expect_true(isTRUE(attr(fin, "gate")$ok))
  expect_equal(fin$final_decision$primary_cell_type, "classical monocyte")
})

test_that("API path fails informatively when the model returns non-JSON text", {
  jin <- read_triage_input(test_path("fixtures", "cluster_1_round1.json"))
  stub_non_json <- function(prompt_json_string, api_key, model = "deepseek-v4-flash",
                            temperature = 0, timeout_seconds = 1200, max_retries = 4,
                            retry_delay = 4, system_prompt = NULL, base_url = NULL) {
    list(ok = TRUE, text = "no json here at all", usage = NULL, status = 200L,
         model = model, error = NULL, retry_count = 0L, error_class = NULL)
  }
  testthat::local_mocked_bindings(invoke_deepseek_api = stub_non_json, .package = "Triage")
  expect_error(
    run_triage_adjudication(jin, use_api = TRUE, api_key = "test-key",
                            dataset_name = "census_immune",
                            project_root = test_path("fixtures")),
    regexp = "could not parse a JSON adjudication object")
})

test_that("API path fails informatively without a key", {
  jin <- read_triage_input(test_path("fixtures", "cluster_1_round1.json"))
  withr::with_envvar(c(DEEPSEEK_API_KEY = "", TRIAGE_LLM_API_KEY_ENV = "DEEPSEEK_API_KEY"), {
    expect_error(
      run_triage_adjudication(jin, use_api = TRUE,
                              dataset_name = "census_immune",
                              project_root = test_path("fixtures")),
      regexp = "no API key")
  })
})
