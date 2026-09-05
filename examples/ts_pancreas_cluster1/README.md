# ts_pancreas_cluster1 — Handling Editor API call (requires API access)

Demonstrates the API path of `run_triage_adjudication()` for `TS_pancreas`
cluster 1 (released identity: B cell, CL:0000236): the Handling Editor
model is called, its raw text is
parsed, and the deterministic post-processing/normalization chain
produces the final adjudication record.

## Requirements

- `DEEPSEEK_API_KEY` (or point `LLM_API_KEY_ENV` at another variable)
- An OpenAI-compatible chat-completions endpoint (full URL including
  `/chat/completions`) in `LLM_API_BASE_URL`, e.g.
  `https://api.deepseek.com/chat/completions`

The demo uses the package default model unless you pass a different one
to `run_triage_adjudication(model = ...)`. This example covers only the
package API path (Handling Editor + deterministic postprocessing); the
full publication benchmark orchestration is in `reproducibility/scripts/`.

## Files

- `run_demo_api.R` — the demo script (does nothing without the key)

## Run

```r
Sys.setenv(DEEPSEEK_API_KEY = "sk-...")       # or use ~/.Renviron
Sys.setenv(LLM_API_BASE_URL = "https://api.deepseek.com/chat/completions")
source("examples/ts_pancreas_cluster1/run_demo_api.R")
```

Note: this demo performs real LLM calls and is intentionally excluded
from the automated test suite. Without `DEEPSEEK_API_KEY` the script
stops with an informative message.
