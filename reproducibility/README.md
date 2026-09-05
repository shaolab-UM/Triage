# Reproducing the manuscript analyses

All scripts below are deterministic (no LLM calls) unless stated. Run
from the repository root. The `Triage` R package must be installed from
this repository (`R CMD INSTALL` of the repo root) — see the top-level
README.

## Primary benchmark

- Input directory: `reproducibility/primary/`
  (`Census_immune`, `Sikkema_lung`, `TS_kidney`, `TS_pancreas`,
  `Zheng_blood`; final adjudication records in `<dataset>/final/`,
  masked DEG inputs in `<dataset>/model_inputs/`, benchmark references in
  `<dataset>/evaluation/`)
- Script: `reproducibility/scripts/release/check_public_json_counts.R`
- Expected output: per-dataset counts 16/11/4/10/9 = 50, printed to stdout
  (exit status 0)

Full pipeline re-execution (requires external LLM API credentials):

- Script: `reproducibility/scripts/run_pipeline.sh`
- Expected output: regenerated run directory under `outputs/<dataset>/`

## External validation

- Input directory: `reproducibility/external/`
  (`Anderson_DLPFC`, `Zha_AD_mouse`, `S6K1_organoid`; cluster exclusions
  in `S6K1_organoid/exclusions.tsv`)
- Script: `reproducibility/scripts/release/validate_public_release.R`
- Expected output: per-dataset counts 18/9/15 = 42 plus structure and
  schema checks, printed to stdout (exit status 0)

## CL-Linker evaluation

- Input directory: `reproducibility/cl_linker/` (gold manual mappings and
  comparator predictions) and
  `reproducibility/primary/reviewer_mapping_registry.tsv`
- Script: `reproducibility/scripts/analysis/evaluate_cl_linker.R`
- Expected output: recomputed comparator and operational tables under
  `results/cl_linker_evaluation/` (gitignored), matching
  `reproducibility/expected_results/`

## Selector controls

- Input directory: `reproducibility/selector_controls/`
  (`selector_inputs.tsv`, `selector_expected_tableS10.tsv`)
- Script: `reproducibility/scripts/analysis/reproduce_selector_controls.R`
- Expected output: `results/selector_controls/selector_summary.tsv` and
  `selector_cluster_results.tsv` (gitignored); the retrospective oracle
  Overall value must reproduce the published Table S10 value
  (84.9083370479893), checked against
  `reproducibility/selector_controls/selector_expected_tableS10.tsv`

## Sensitivity analyses

- Input directory: `reproducibility/sensitivity/`
  (`cellmarkerdb_additional_reviewer_summary.tsv`,
  `cellmarkerdb_additional_reviewer_by_dataset.tsv`)
- Script:
  `reproducibility/scripts/sensitivity/cellmarkerdb/validate_cellmarkerdb_sensitivity_summary.R`
- Expected output: consistency checks of the summary against the
  by-dataset table, printed to stdout (exit status 0)

## Expected results

- Directory: `reproducibility/expected_results/`
- Contents: `primary_adjudication.jsonl` and
  `validation_adjudication.jsonl` (combined release records),
  `primary_summary.tsv` and `external_summary.tsv` (per-cluster released
  decisions), `selector_summary.tsv` and the CL-Linker summary tables
  (published comparator/operational outcomes)
- Use: byte-level reference for any recomputation; the release checksums
  are recorded in `MANIFEST_SHA256.tsv`
