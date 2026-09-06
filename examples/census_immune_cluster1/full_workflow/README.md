# Full end-to-end workflow example — Census_immune cluster_1

This directory archives a **real, single-command end-to-end run** of the Triage
workflow on one cluster, executed with
`Rscript reproducibility/scripts/run_triage.R --benchmark Census_immune --cluster-id cluster_1`
and reproduced here as a compact example.

- The input is **anonymized cluster-level differential-expression data**
  (`input/maskdeg.csv`: 100 DEG rows for `cluster_1` only).
- The pipeline ran from DEG input through CASSIA annotation, candidate
  construction, evidence preparation, the in-house and enrichment reviewers,
  CL-Linker mapping, adjudication-input assembly, Handling Editor adjudication,
  Chief QC, and final summarization.
- **Reference labels are not used anywhere in adjudication.** No
  `true_label.csv` or `reference_cl.tsv` exists in this example; identity
  decisions are derived only from the DEG/evidence inputs.
- Fresh LLM executions can differ in wording and confidence from this archived
  example because model output is stochastic. The pipeline, prompts, and
  deterministic post-processing rules are fixed.

## Contents

```text
input/
  maskdeg.csv                       anonymized DEG input (cluster_1)
outputs/
  cassia_annotation.csv             stage 03b: CASSIA annotation + score
  structured_candidates.csv         stage 04a: marker-based candidate list
  in_house_summary.csv              stage 07: in-house reviewer summary
  enrichment_summary.csv            stage 07b: enrichment reviewer summary
  reviewer_mapping_registry.tsv     stage 07.5: CL-Linker reviewer mappings
  judge_input.json                  stage 08: adjudication input dossier
  final_adjudication.json           stage 09: final adjudication record
  final_summary.csv                 stage 09/10: final summary row
```

## Archived result

The final adjudication for `cluster_1` is `monocyte` (`CL:0000576`,
confidence 0.96, decision category `tie`) with release state
`release_with_audit_note` and `chief_qc_status = passed`. The more specific
candidate `classical monocyte` (`CL:0000860`) was supported by two of three
reviewers; the deterministic rules resolved the panel to the consensus-locked
parent identity recorded in `final_adjudication.json` (see its
`decision_trace` for the full reasoning trace).

Source-path provenance columns inside `reviewer_mapping_registry.tsv` were
generalized to `<run_dir>` so that no local temporary paths are embedded; all
scientific fields are unchanged.

## Reproducing

See the repository `README.md` ("Full Triage workflow from DEG input") for the
required API credentials, external resources, and the benchmark command used to
produce this run.
