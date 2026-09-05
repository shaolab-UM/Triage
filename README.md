# Triage

Triage is a downstream adjudication workflow for single-cell cell-type
annotation: it begins after candidate annotations have been generated and
adjudicates them against three reviewer outputs (a CASSIA agent, an in-house
marker-based annotator, and an enrichment reviewer), biological evidence,
and Cell Ontology constraints, through a deterministic normalization/gating
chain and an LLM Handling Editor. Chief QC and final release-state
assignment run in the batch pipeline scripts.

The repository has **two components**:

1. **The `Triage` R package** — an installable package exposing a compact
   public API. It loads without any API keys; network access is only needed
   for the optional API-backed adjudication path.
2. **`reproducibility/`** — manuscript data, validation scripts and expected
   results for the Triage v1.0.0 publication release.

## Navigation

- [Installation](#installation)
- [Quick start](#quick-start)
- [What Triage does](#what-triage-does)
- [Example workflows](#example-workflows)
- [Package structure](#package-structure)
- [Reproducing the manuscript analyses](#reproducing-the-manuscript-analyses)
- [Model configuration](#model-configuration)
- [External resources](#external-resources)
- [Reproducibility notes](#reproducibility-notes)
- [Citation](#citation)
- [License](#license)

## Installation

```bash
R CMD build .
R CMD INSTALL Triage_1.0.0.tar.gz
```

Or from GitHub: `remotes::install_github("shaolab-UM/Triage")`.

## Quick start

Deterministic demonstration — postprocessing of a precomputed Handling
Editor draft, no API key required (see `examples/census_immune_cluster1_demo/`):

```r
library(Triage)

ontology <- load_triage_ontology()
jin <- read_triage_input("examples/census_immune_cluster1_demo/cluster_1_round1.json")

fin <- run_triage_adjudication(
  jin,
  ontology = ontology,
  use_api = FALSE,
  head_output = "examples/census_immune_cluster1_demo/head_round1.json",
  dataset_name = "census_immune",
  project_root = "reproducibility/primary"
)
attr(fin, "gate")$ok                      # deterministic local gate
fin$final_decision$primary_cell_type      # classical monocyte
validate_triage_result(fin, jin)
```

See the vignette: `vignette("getting-started", package = "Triage")`.

## What Triage does

For each cluster, Triage:

1. Parses and normalizes the three reviewer summaries against the Cell
   Ontology (`map_cell_ontology`, CL-Linker `run_cl_linker`).
2. Assembles the adjudication dossier (`build_adjudication_input`).
3. Runs the head-editor model over the dossier (API path) or consumes a
   precomputed head-editor output (no-API path).
4. Applies the deterministic post-processing chain: release policy,
   CL normalization, citation/label gates, and the local validation gate
   (`run_triage_adjudication`, `validate_triage_result`).

Chief QC and final release-state assignment are not part of the package
API; they run in `reproducibility/scripts/pipeline/09_run_judge.R`.

## Example workflows

- `examples/census_immune_cluster1_demo/` — deterministic no-API
  postprocessing of a precomputed Handling Editor draft for
  `Census_immune` cluster 1 (released identity: classical monocyte,
  CL:0000860). The precomputed draft is supplied for demonstration; it is
  not an archived Handling Editor round output from the original run.
- `examples/ts_pancreas_cluster1_demo/` — Handling Editor API call for
  `TS_pancreas` cluster 1 (released identity: B cell, CL:0000236; requires
  `DEEPSEEK_API_KEY`; not part of the test suite). Chief QC is not invoked
  by the package.

## Package structure

```text
R/                        public API + ontology, CL-Linker, adjudication engine
inst/extdata/ontology/    Cell Ontology JSON snapshot (v2025-07-30)
inst/prompts/, inst/schemas/, inst/config/
vignettes/getting-started.Rmd
tests/testthat/           deterministic tests (no API calls)
examples/                 census_immune_cluster1_demo (no API), ts_pancreas_cluster1_demo (API)
```

## Reproducing the manuscript analyses

Layout:

```text
reproducibility/primary/        5 primary datasets (final JSONs, model inputs, evaluation files)
reproducibility/external/       3 external-validation datasets (42 final JSONs, S6K1 exclusions)
reproducibility/cl_linker/      CL-Linker evaluation source data
reproducibility/sensitivity/    CellMarkerDB additional-reviewer sensitivity data
reproducibility/selector_controls/  selector inputs + expected Table S10 values
reproducibility/expected_results/   release JSONL records + summary tables
reproducibility/scripts/        pipeline, analysis, release and sensitivity scripts
reproducibility/config/         adjudication profile and run provenance
docs/                           methodology and boundary documentation
```

Reported counts — primary benchmark (50 clusters): Census immune 16,
Sikkema lung 11, Tabula Sapiens kidney 4, Tabula Sapiens pancreas 10,
Zheng blood 9. External validation (42 clusters): Anderson DLPFC 18,
Zha AD mouse 9, S6K1 organoid 15.

Deterministic release checks (no LLM calls):

```bash
Rscript reproducibility/scripts/release/check_public_json_counts.R
Rscript reproducibility/scripts/release/validate_public_release.R
Rscript reproducibility/scripts/analysis/evaluate_cl_linker.R
Rscript reproducibility/scripts/analysis/reproduce_selector_controls.R
```

The selector check reproduces the published Table S10 / Fig. 3C-D controls,
including the retrospective oracle value 84.908337 (Overall; see
`docs/SELECTOR_REPRODUCTION.md`).

The full primary pipeline can be rerun with:

```bash
bash reproducibility/scripts/run_pipeline.sh \
  --dataset Census_immune \
  --masked-deg reproducibility/primary/Census_immune/model_inputs/maskdeg.csv
```

Reference labels are only introduced when `--true-label-csv` is explicitly
supplied for Step 11 evaluation.

## Model configuration

The primary benchmark used `deepseek-v4-flash` for every LLM role
(CASSIA/In-house/clusterProfiler reviewer-side calls, Handling Editor and
Chief QC) at temperature 0. A separate profiling / model-comparison
experiment used `deepseek-reasoner` as Handling Editor and `deepseek-chat`
as Chief QC; alternatives (GPT-5.4, Gemini-3 Pro Preview, Claude Opus 4.5)
used the same model for both roles with 1 run per dataset (DeepSeek:
5 runs/dataset). See `docs/MODEL_CONFIGURATION.md` and
`reproducibility/config/model_run_provenance.tsv`.

Never commit API credentials; placeholder values are written as `XXXXX`
(`.env.example` documents the variables, e.g. `DEEPSEEK_API_KEY`).

## External resources

The Cell Ontology JSON ships with the package
(`system.file("extdata/ontology", package = "Triage")`; env override
`CL_LOCAL_JSON`). STRING/PPI, CollecTRI and CellMarkerDB spreadsheets are
not redistributed; place local copies under `inputs/raw/` as specified in
`resources/README.md` (checksums in `resources/CHECKSUMS.tsv`).

## Reproducibility notes

- Final adjudication JSONs, counts and expected results are frozen; the
  release checksums are recorded in `MANIFEST_SHA256.tsv`.
- Model-generated interpretive text (reviewer rationales, evidence prose)
  should not be read as independently measured expression evidence;
  quantitative benchmark analyses do not use these prose fields.
- Reviewer metadata in the public final JSONs preserves reviewer labels/CL
  assignments while `final_decision` remains the single
  publication-facing authoritative state (`docs/REFERENCE_BOUNDARY.md`,
  `docs/DATA_LAYOUT.md`).
- The retained snapshot did not include a full `renv.lock`; supported
  versions are summarized in `sessionInfo.txt` and `docs/DEPENDENCIES.md`.

## Citation

Citation information will be provided upon publication. Software release:
Triage v1.0.0 (`VERSION`, `docs/PROVENANCE.md`).

## License

MIT — see `LICENSE.md`.
