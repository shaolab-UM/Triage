# ts_pancreas_cluster1 — Handling Editor API call (requires API access)

Demonstrates the API path of `run_triage_adjudication()` for `TS_pancreas`
cluster 1 (released identity: B cell, CL:0000236): the Handling Editor
model is called, its raw text is
parsed, and the deterministic post-processing/normalization chain
produces the final adjudication record.

## Requirements

- An API key
- An OpenAI-compatible chat-completions endpoint (full URL including
  `/chat/completions`), e.g. `https://api.deepseek.com/chat/completions`

Supply both directly in `run_demo_api.R` (`api_key` / `api_base_url`), or
set the variables to `NULL` and provide `DEEPSEEK_API_KEY` /
`LLM_API_BASE_URL` as environment variables instead.

The demo uses the package default model unless you pass a different one
to `run_triage_adjudication(model = ...)`. This example covers only the
package API path (Handling Editor + deterministic postprocessing); the
full publication benchmark orchestration is in `reproducibility/scripts/`.

## Files

- `run_demo_api.R` — the demo script (stops with an informative message
  until you replace the placeholder key)

## Run

```r
# 1. Replace api_key (and api_base_url if needed) at the top of the script.
# 2. Then:
source("examples/ts_pancreas_cluster1/run_demo_api.R")
```

Note: this demo performs real LLM calls and is intentionally excluded
from the automated test suite.
