# =========================================================================
# provider_deepseek.R — DeepSeek-compatible API transport
# Extracted verbatim from scripts/pipeline/09_run_judge.R (v1.0.0 release).
# Endpoint/key/model/temperature/retry handling only; no scientific logic.
# =========================================================================

invoke_deepseek_api <- function(prompt_json_string, api_key, model="deepseek-v4-flash",
                                temperature=0, timeout_seconds=1200, max_retries=4, retry_delay=4,
                                system_prompt = NULL) {
  if (is.null(api_key) || !nzchar(api_key)) stop("Missing API key in environment.")
  
  retries <- 0
  last_err <- NULL
  last_status <- NA_integer_
  last_err_class <- "unknown"
  
  while (retries < max_retries) {
    sys_msg <- if (!is.null(system_prompt) && nzchar(as.character(system_prompt))) {
      paste(
        as.character(system_prompt),
        "Return ONLY valid JSON (no markdown fences, no extra text).",
        "Use JSON null (without quotes) for nulls.",
        "Do NOT add extra top-level keys.",
        sep = "\n"
      )
    } else {
      paste(
        "Return ONLY valid JSON (no markdown fences, no extra text).",
        "Use JSON null (without quotes) for nulls.",
        "Do NOT add extra top-level keys."
      )
    }
    req_body <- list(
      model = model,
      messages = list(
        list(role="system", content=sys_msg),
        list(role="user", content=prompt_json_string)
      ),
      temperature = temperature,
      stream = FALSE
    )
    
    resp <- tryCatch({
      httr::POST(
        url=DEEPSEEK_BASE_URL,
        httr::add_headers(
          `Content-Type`="application/json",
          `Authorization`=paste("Bearer", api_key)
        ),
        body=jsonlite::toJSON(req_body, auto_unbox=TRUE, null="null"),
        encode="raw",
        httr::timeout(timeout_seconds)
      )
    }, error=function(e) e)
    
    if (inherits(resp, "error")) {
      last_err <- as.character(resp$message %||% "http_error")
      last_err_class <- classify_http_error(NA_integer_, last_err)
      retries <- retries + 1
      Sys.sleep(calc_backoff(retry_delay, retries, last_err_class))
      next
    }
    
    status <- httr::status_code(resp)
    last_status <- status
    if (status == 200) {
      content <- httr::content(resp, as="parsed")
      out <- tryCatch(content$choices[[1]]$message$content, error=function(e) NULL)
      usage <- content$usage %||% NULL
      if (!is.null(out) && nzchar(out)) {
        return(list(ok=TRUE, text=out, usage=usage, status=status, model=model, error=NULL, retry_count=retries, error_class=NULL))
      } else {
        last_err <- "empty_content"
        last_err_class <- "empty_content"
      }
    } else {
      last_err <- tryCatch({
        txt <- httr::content(resp, as="text", encoding="UTF-8")
        paste0("status_", status, "_", substr(txt, 1, 200))
      }, error=function(e) paste0("status_", status))
      last_err_class <- classify_http_error(status, last_err)
    }
    
    retries <- retries + 1
    Sys.sleep(calc_backoff(retry_delay, retries, last_err_class))
  }
  
  list(ok=FALSE, text=NULL, usage=NULL, status=last_status, model=model, error=last_err %||% "unknown", retry_count=retries, error_class=last_err_class)
}

usage_to_fields <- function(usage) {
  if (is.null(usage) || !is.list(usage)) {
    return(list(prompt_tokens=NA_integer_, completion_tokens=NA_integer_, total_tokens=NA_integer_))
  }
  list(
    prompt_tokens = suppressWarnings(as.integer(usage$prompt_tokens %||% NA_integer_)),
    completion_tokens = suppressWarnings(as.integer(usage$completion_tokens %||% NA_integer_)),
    total_tokens = suppressWarnings(as.integer(usage$total_tokens %||% NA_integer_))
  )
}

append_jsonl <- function(path, obj) {
  ensure_dir(dirname(path))
  line <- jsonlite::toJSON(obj, auto_unbox=TRUE, null="null")
  cat(line, "\n", file=path, append=TRUE)
}

make_request_id <- function(cid, stage, round) {
  paste0(cid, "_", stage, "_", round, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", paste0(sample(c(letters, 0:9), 6, replace=TRUE), collapse=""))
}

# --------------------------
# Scheme A constraints
# --------------------------
