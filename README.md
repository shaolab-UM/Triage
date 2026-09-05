# Triage

Triage is an adjudication workflow for cell-type annotation that integrates three reviewer outputs, biological evidence and Cell Ontology constraints. This repository contains the core workflow, CL-Linker utilities, publication-facing primary/validation outputs and deterministic downstream controls used in the manuscript.

## Release status

This repository contains the minimal public code and data release for Triage v1.0.0. It includes the reported primary benchmark (50 clusters), external validation (42 clusters), core Triage/CL-Linker workflow code and retained downstream evaluation/control scripts.

## Credentials

Never commit API credentials. Runtime-specific placeholder values are written as `XXXXX`. Use environment variables such as:

```text
DEEPSEEK_API_KEY
DISGENET_API_KEY
TRIAGE_HOME
PROJECT_ROOT
CL_LOCAL_JSON
TRIAGE_PPI_ROOT
```


## Repository layout

- `scripts/pipeline/` — primary workflow stages.
- `scripts/analysis/` — deterministic selector-control reproduction.
- `scripts/sensitivity/` — publication-level CellMarkerDB sensitivity check.
- `scripts/release/` — public JSON count and schema checks.
- `lib/` — CL-Linker, CL similarity and shared utilities.
- `config/` — dataset, ontology-normalization, prompt and formal adjudication settings.
- `data/primary/` — masked DEG inputs, final adjudication inputs, evaluation files and 50 public final JSONs.
- `data/validation/` — 42 public final validation JSONs.
- `data/release/` — combined primary/validation JSONL files.
- `resources/` — external-resource specifications and checksums.
- `docs/` — reproducibility, reference-boundary and audit notes.

## Primary benchmark inputs

For each of the five primary datasets:

```text
data/primary/<dataset>/model_inputs/maskdeg.csv
data/primary/<dataset>/handling_editor_round1_inputs/*_round1.json
data/primary/<dataset>/evaluation/true_label.csv
data/primary/<dataset>/evaluation/cluster_map.csv
data/primary/<dataset>/evaluation/reference_cl.tsv
data/primary/<dataset>/final/cluster_*.json
```

Model-facing inputs and evaluation references are intentionally separated. See `docs/REFERENCE_BOUNDARY.md`.

## Reported counts

Primary benchmark:

- Census immune: 16
- Sikkema lung: 11
- Tabula Sapiens kidney: 4
- Tabula Sapiens pancreas: 10
- Zheng blood: 9
- Total: 50

External validation:

- Anderson DLPFC: 18
- Zha AD mouse: 9
- S6K1 organoid: 15
- Total: 42


## Formal primary adjudication settings

The primary benchmark used `deepseek-v4-flash` as handling editor and `deepseek-chat` for Chief QC, with temperature 0. Exact adjudication settings are in `config/primary_adjudication_profile.tsv`. See `docs/REPRODUCE_PRIMARY_ADJUDICATION.md`.

## Deterministic selector controls

```bash
Rscript scripts/analysis/reproduce_selector_controls.R   --repo-root /path/to/Triage   --cl-json /path/to/CL-ontology-v2025-07-30.json
```

Expected overall values:

- Majority vote: 41.3207
- Top reviewer after percentile normalization: 61.6773
- Ontology-only control: 62.9211
- Triage: 79.9790
- Retrospective oracle: 84.9083 (reference-using maximum across the three reviewer outputs plus Triage, including the publication's label-mapped retrospective candidates; Table S10)

## CellMarkerDB sensitivity

The baseline enrichment reviewer reads CellMarkerDB through `scripts/pipeline/06b_run_inter.R`. The additional-reviewer sensitivity source values are under `data/sensitivity/`.

The CellMarkerDB spreadsheets are not redistributed. Place local copies at:

```text
inputs/raw/cellmarker/Cell_marker_Human.xlsx
inputs/raw/cellmarker/Cell_marker_Mouse.xlsx
```

Exact file checksums are in `resources/cellmarkerdb/CHECKSUMS.tsv`. See `docs/CELLMARKERDB_SENSITIVITY.md`.

## External resources

The workflow also expects Cell Ontology, STRING and CollecTRI resources under `inputs/raw/`; see `resources/README.md`.

## Release validation

After installing the required R environment:

```bash
Rscript scripts/release/check_public_json_counts.R
Rscript scripts/release/validate_public_release.R
Rscript scripts/sensitivity/cellmarkerdb/validate_cellmarkerdb_sensitivity_summary.R
```

The selector check additionally requires the Cell Ontology JSON.

## Reproducibility limitations

The retained project snapshot did not contain the original complete `renv.lock` or full `sessionInfo()` output. The software versions supported by the retained record are summarized in `sessionInfo.txt` and `docs/DEPENDENCIES.md`.


## Core workflow runner

A clean workflow wrapper is provided at:

```bash
bash scripts/run_pipeline.sh   --dataset Census_immune   --masked-deg data/primary/Census_immune/model_inputs/maskdeg.csv
```

The wrapper starts from anonymized cluster-level DEG input and runs the core reviewer, CL-Linker and adjudication stages through the final summary. Reference labels are only introduced when `--true-label-csv` is explicitly supplied for Step 11 evaluation.

See `scripts/run_pipeline.sh --help` for runtime requirements.

## CL-Linker evaluation

```bash
Rscript scripts/analysis/evaluate_cl_linker.R   --repo-root /path/to/Triage   --cl-json /path/to/CL-ontology-v2025-07-30.json
```

See `docs/CL_LINKER_EVALUATION.md`.

## Repository scope

Final manuscript figure-rendering and table-formatting scripts are not part of this public code package. The repository provides the core Triage/CL-Linker workflow, publication-facing cluster-level outputs and deterministic downstream evaluation/control scripts.


## Provenance

This repository uses `v1.0.0` as the public software release version. See `docs/PROVENANCE.md`.

## Historical model configuration

The reported Chief QC identifier `deepseek-chat` is retained as historical
configuration metadata. Its availability for a new execution depends on the API
provider or compatible gateway. See `docs/MODEL_CONFIGURATION.md`.

## Model-generated interpretive text

Reviewer rationales and final evidence summaries are model-generated
interpretive text. They may incorporate prior biological knowledge, including
canonical markers that are not present in the supplied DEG list. Such mentions
should not be interpreted as independently measured expression evidence.
Quantitative benchmark analyses do not use these prose fields.


## Model/provider provenance

Historical model identifiers and provider category are recorded in `docs/MODEL_CONFIGURATION.md` and `config/model_run_provenance.tsv`. Credentials and the exact runtime endpoint string are not embedded in the public repository.


## Publication-facing reviewer metadata

The public final JSON omits adjudication-relative reviewer fields
(`final_cl_id`, `support_class`, `is_correct` and `matches_final_label`).
Reviewer labels, reviewer CL assignments and model-generated rationale are
retained, while `final_decision` is the single authoritative publication-facing
final state.

Within reviewer objects, `reviewer_cl_id` preserves the source-facing reviewer
mapping and `cell_ontology_id` records the CL assignment in the Handling
Editor's reviewer-specific method verdict. These values may differ after
normalization or refinement. Selector analyses use the frozen reviewer CL IDs
in `data/primary/selector_inputs.tsv`.
