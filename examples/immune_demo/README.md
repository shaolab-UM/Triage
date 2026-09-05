# immune_demo — deterministic, no API required

Reproduces the adjudication decision for `Census_immune` cluster 1
without any LLM calls, using a precomputed head-editor output.

## Files

- `run_demo.R` — the demo script
- `head_round1.json` — precomputed head-editor adjudication for cluster 1
  (the reviewer-input dossier is read from `reproducibility/primary/`)

## Run

From the repository root:

```r
source("examples/immune_demo/run_demo.R")
```

Expected output: `classical monocyte` / `CL:0000860`, decision category
`cassia_better`, local gate OK.

No API key is needed. The demo calls `run_triage_adjudication()` with
`use_api = FALSE`.
