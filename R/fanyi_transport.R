# Stage-06b LLM transport contract.
#
# The tested clusterProfiler revision (4.19.4.008 @ f9f0d502...) routes
# interpret() LLM calls through fanyi::chat_request, whose DeepSeek transport
# (fanyi 0.1.0, the internal .deepseek_query_messages) hard-codes the
# DeepSeek /v1 endpoint and ignores the canonical LLM_API_BASE_URL contract
# used by every other Triage LLM stage. fanyi 0.1.1 keeps the same hard-coded
# endpoint and still has no base-url parameter, so there is no supported
# upstream mechanism. The helpers below keep stage 06b on the SAME effective
# full endpoint as every other Triage stage without changing clusterProfiler,
# fanyi, or any prompt: a temporary, locally-installed transport shim
# (the verbatim pinned fanyi transport with a configurable URL) is installed
# for the duration of one interpret() call and always restored on exit.

.triage_fanyi_default_endpoint <- "https://api.deepseek.com/v1/chat/completions"
.triage_tested_fanyi_version <- function() "0.1.0"
.triage_fanyi_version_matches <- function(version = .triage_tested_fanyi_version()) {
  if (!requireNamespace("fanyi", quietly = TRUE)) return(FALSE)
  identical(as.character(utils::packageDescription("fanyi")$Version), version)
}

# Verbatim fanyi 0.1.0 DeepSeek transport with the endpoint taken from the
# Triage contract instead of the hard-coded constant. Key/model fallbacks,
# request body, status handling, and error format are identical.
.triage_fanyi_transport <- function(messages, model = NULL, api_key = NULL,
                                    max_tokens = 4096, base_url, ...) {
  # fanyi's key/model store is only consulted when the caller did not supply
  # explicit values; an absent store must not break explicit-key calls.
  key_info <- if (is.null(model) || is.null(api_key)) {
    tryCatch(utils::getFromNamespace("get_translate_appkey", "fanyi")("dsk"),
             error = function(e) NULL)
  } else NULL
  if (is.null(model)) {
    user_model <- if (!is.null(key_info) && !is.null(key_info$user_model)) {
      key_info$user_model
    } else "deepseek-chat"
  } else {
    user_model <- model
  }
  if (is.null(api_key)) {
    api_key <- if (!is.null(key_info)) key_info$key else NULL
  }
  if (is.null(api_key) || !nzchar(api_key)) {
    stop("API key for deepseek is missing.")
  }
  url <- base_url
  body <- list(model = user_model, messages = messages, stream = FALSE,
               max_tokens = max_tokens)
  response <- httr2::req_perform(httr2::req_body_json(
    httr2::req_headers(httr2::request(url),
                       `Content-Type` = "application/json",
                       Authorization = paste("Bearer", api_key)),
    body))
  if (httr2::resp_status(response) != 200) {
    error_content <- httr2::resp_body_json(response)
    error_msg <- if (!is.null(error_content$error$message)) {
      error_content$error$message
    } else {
      httr2::resp_body_string(response)
    }
    stop(sprintf("API request failed: %s", error_msg))
  }
  response
}

.triage_interpret <- function(enrich_list, context = NULL,
                              model = "deepseek-v4-flash", api_key = NULL,
                              base_url, n_pathways = 20,
                              task = "interpretation", ...) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(".triage_interpret: the enrichment reviewer requires the pinned ",
         "clusterProfiler revision (see install_triage_dependencies()).")
  }
  if (!requireNamespace("fanyi", quietly = TRUE)) {
    stop(".triage_interpret: the enrichment reviewer requires fanyi ",
         "(install with remotes::install_version(\"fanyi\", version = \"0.1.0\")).")
  }
  if (is.null(base_url) || !nzchar(base_url)) {
    stop(".triage_interpret: LLM_API_BASE_URL is not set; stage 06b posts to ",
         "the same full chat-completions endpoint as every other Triage LLM stage.")
  }
  call_interpret <- function() {
    clusterProfiler::interpret(enrich_list, context = context,
                               n_pathways = n_pathways, model = model,
                               api_key = api_key, task = task, ...)
  }
  if (identical(base_url, .triage_fanyi_default_endpoint)) {
    return(call_interpret())
  }
  # Temporary transport shim: point the pinned fanyi DeepSeek transport at
  # the canonical endpoint for exactly this call, then always restore.
  ns <- asNamespace("fanyi")
  binding <- ".deepseek_query_messages"
  if (!exists(binding, envir = ns, inherits = FALSE)) {
    stop(".triage_interpret: the pinned fanyi DeepSeek transport ",
         "(fanyi:::.deepseek_query_messages) is unavailable; the tested ",
         "fanyi 0.1.0 transport is required for custom endpoints.")
  }
  orig <- get(binding, envir = ns, inherits = FALSE)
  shim <- function(messages, model = NULL, api_key = NULL,
                   max_tokens = 4096, ...) {
    .triage_fanyi_transport(messages, model = model, api_key = api_key,
                            max_tokens = max_tokens, base_url = base_url, ...)
  }
  was_locked <- bindingIsLocked(binding, ns)
  unlock_binding <- get("unlockBinding", envir = asNamespace("base"))
  lock_binding <- get("lockBinding", envir = asNamespace("base"))
  unlock_binding(binding, ns)
  assign(binding, shim, envir = ns)
  restore <- function() {
    assign(binding, orig, envir = ns)
    if (was_locked) lock_binding(binding, ns)
  }
  on.exit(restore(), add = TRUE)
  call_interpret()
}
