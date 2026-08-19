# CellMarkerDB additional-reviewer sensitivity

## Analysis definition

The sensitivity analysis added a CellMarkerDB-derived ORA assignment as a fourth reviewer. Cluster DEGs were analyzed against CellMarkerDB. The resulting assignment entered adjudication as a non-core reviewer.

Two conditions were reported:

- CellMarkerDB added without reliability filtering.
- CellMarkerDB added with reliability filtering.

For non-core reviewer candidates, the adjudication code uses a reliability threshold of 0.55 and an evidence-alignment threshold of 0.45. If a separate alignment value is unavailable, reliability is used for the alignment check. The three core reviewers (`cassia`, `our`, `enrich`) are preserved by the gate.

The non-core reliability score in `09_run_judge.R` combines:

- valid CL mapping: weight 0.35;
- mapping quality: weight 0.30;
- ontology proximity to the lock: weight 0.25;
- label specificity: weight 0.10.

The code also marks no-CL, generic-label or large lineage-distance cases as low reliability.

## Publication-final results

The release contains two derived data files:

- `data/sensitivity/cellmarkerdb_additional_reviewer_summary.tsv`: overall mean similarity, paired delta, bootstrap CI, exact rate, manual-review rate and paired P values from Table S7.
- `data/sensitivity/cellmarkerdb_additional_reviewer_by_dataset.tsv`: the five dataset-specific paired deltas and pooled overall delta used for Fig. 3A, together with the corresponding full-Triage means.

The reported means are:

- CellMarkerDB added: 78.7712 (delta -1.2078).
- CellMarkerDB added with reliability filtering: 80.1512 (delta +0.1723).

## Reproducibility boundary

The release contains the baseline CellMarkerDB loader (`06b_run_inter.R`), the adjudication reliability/gating implementation (`09_run_judge.R`) and the publication sensitivity summaries.

The separate script that generated the additional CellMarkerDB reviewer assignments is not included in this release. Therefore the release does not claim to regenerate those fourth-reviewer assignments from raw DEG files. It does preserve the exact resource specification, gate implementation and publication-derived results needed to reconstruct the reported Fig. 3A and Table S7 values.

The released derived sensitivity results can be reproduced at the summary and figure-source-data level with the included scripts.


## Scope of the released checker

`validate_cellmarkerdb_sensitivity_summary.R` is a summary-level consistency check. It validates the retained publication source tables and does not regenerate the historical fourth-reviewer assignments from DEG files.
