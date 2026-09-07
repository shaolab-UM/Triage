# Reported primary benchmark example — Census_immune cluster_1

This example contains the released input and final adjudication records for
`Census_immune` cluster 1 from the primary benchmark reported in the
manuscript.

The reported result is:

- **Triage:** classical monocyte, CL:0000860, confidence 0.93

Published reference:

- CD14-positive monocyte, CL:0001054

Reference labels are used only for post-adjudication evaluation (the
`evaluation` block in the final record; they are never used by the
adjudication stages).

## Files

- `input/maskdeg.csv` — the released masked DEG input for cluster 1
  (deterministic `cluster_1` subset of
  `reproducibility/primary/Census_immune/model_inputs/maskdeg.csv`)
- `outputs/final_adjudication.json` — the released final adjudication
  record (`reproducibility/primary/Census_immune/final/cluster_1.json`),
  including the final decision, reviewer records, evidence, QC disposition
  and the evaluation block with the published reference and CL similarity
  (78.3778251389102)
- `outputs/reviewer_mapping_registry.tsv` — the released CL-Linker mapping
  rows for cluster 1 (extracted from
  `reproducibility/primary/reviewer_mapping_registry.tsv`)
- `outputs/final_summary.csv` — the cluster's row from the released primary
  benchmark summary (`reproducibility/expected_results/primary_summary.tsv`)
- `cluster_1_round1.json` — the released reviewer-input record
- `head_round1.json` — a precomputed Handling Editor fixture for the
  deterministic package demo below

## Deterministic package demo (no API key required)

```r
source("examples/census_immune_cluster1/run_demo.R")
```

This runs the deterministic package API (postprocessing of the precomputed
Handling Editor fixture; no network access) and returns the reported
identity: classical monocyte, CL:0000860.

## Reproducing the full workflow

The full workflow runner, its requirements and the released benchmark data
are described in the top-level `README.md`
(`--benchmark Census_immune --cluster-id cluster_1`).
