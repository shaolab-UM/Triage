# pancreas_demo — API-backed adjudication (requires API access)

Demonstrates the full API path of `run_triage_adjudication()` for a
`TS_pancreas` cluster: the head editor model is called, its raw text is
parsed, and the deterministic post-processing/normalization chain
produces the final adjudication record.

## Requirements

- A DeepSeek-compatible endpoint:
  - `DEEPSEEK_API_KEY` (or point `LLM_API_KEY_ENV` at another variable)
  - `LLM_API_BASE_URL` (default: DeepSeek production endpoint)
- Model configuration (matches the final manuscript):
  - Primary workflow: all roles use `deepseek-v4-flash`, temperature 0
  - Profiling experiment (separate): Handling Editor `deepseek-reasoner`,
    Chief QC `deepseek-chat`
  - Alternatives: GPT-5.4 / Gemini-3 Pro Preview / Claude Opus 4.5,
    same model for both roles, 1 run per dataset (DeepSeek: 5 runs/dataset)

## Files

- `run_demo_api.R` — the demo script (does nothing without the key)

## Run

```r
Sys.setenv(DEEPSEEK_API_KEY = "sk-...")       # or use ~/.Renviron
Sys.setenv(LLM_API_BASE_URL = "https://api.deepseek.com")
source("examples/pancreas_demo/run_demo_api.R")
```

Note: this demo performs real LLM calls and is intentionally excluded
from the automated test suite. Without `DEEPSEEK_API_KEY` the script
stops with an informative message.
