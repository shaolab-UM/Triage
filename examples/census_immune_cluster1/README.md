# census_immune_cluster1 — deterministic postprocessing using a precomputed Handling Editor output

Demonstrates the deterministic adjudication core for `Census_immune`
cluster 1 (released identity: classical monocyte, CL:0000860). No API call
is made: `run_triage_adjudication(use_api = FALSE)` requires a precomputed
Handling Editor output and executes only deterministic downstream
processing (release policy, CL normalization, citation gates, local
validation gate). The Handling Editor is NOT executed.

`head_round1.json` is a precomputed Handling Editor fixture for the
deterministic package example. `cluster_1_round1.json` is the corresponding
reviewer-input record. These files support the deterministic package example;
the archived end-to-end workflow example is provided under `full_workflow/`.

## Files

- `run_demo.R` — the demo script
- `cluster_1_round1.json` — released reviewer-input record for cluster 1
- `head_round1.json` — precomputed Handling Editor draft used as
  demonstration input

## Run

From the repository root:

```r
source("examples/census_immune_cluster1/run_demo.R")
```

Expected output: `classical monocyte` / `CL:0000860`, decision category
`cassia_better`, local gate OK.
