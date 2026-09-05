# census_immune_cluster1_demo — deterministic, no API required

Demonstrates the deterministic adjudication core for `Census_immune`
cluster 1 (released identity: classical monocyte, CL:0000860) without any
LLM calls, by postprocessing a precomputed Handling Editor draft.

Provenance note: `head_round1.json` is supplied for demonstration and
mirrors the released decision; it is not an archived Handling Editor round
output from the original run (round outputs were not retained in the
release). `cluster_1_round1.json` is the released reviewer-input record
(byte-identical copy).

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
