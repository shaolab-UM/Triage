# Software dependencies

Dependencies fall into two distinct groups:

1. **`Triage` package core** — declared in `DESCRIPTION` and installed with
   the package.
2. **Full workflow** — the additional R packages used by the pipeline
   stages (03a–10). These are checked by `triage_preflight()` (installed
   package) or `reproducibility/scripts/preflight_check.R` (repository
   runner) before a run starts; they are intentionally **not** declared in
   `DESCRIPTION`, so a plain package install stays lightweight. Install
   them with `install_triage_dependencies()`.

The installed-package workflow runtime (`run_triage()`) executes the same
stage scripts bundled byte-identically under `inst/workflow/` (see
`inst/workflow/README.md`); the repository `reproducibility/` tree is the
manuscript reproduction layer. External resources are resolved from the
Triage user-data directory (`tools::R_user_dir("Triage", "data")`)
populated by `setup_triage_resources()`; STRING is additionally honored
via `TRIAGE_PPI_ROOT`.

The manuscript reports R 4.5.1; the versions below were tested on
R 4.5.2. Versions are recorded as *tested* versions, not requirements.

## 1. Triage package core

Declared in `DESCRIPTION` (`Imports`):

- cellmarkeraccordion (GitHub: `TebaldiLab/cellmarkeraccordion`, declared in
  `Remotes`; tested 1.0.0)
- data.table (tested 1.18.2)
- dplyr (tested 1.2.1)
- httr (tested 1.4.8)
- jsonlite (tested 2.0.0)
- ontologyIndex
- ontologySimilarity
- purrr (tested 1.2.1)
- readr (tested 2.2.0)
- rlang (tested 1.3.0)
- stats, utils (base)
- stringdist
- stringr (tested 1.6.0)
- tibble (tested 3.3.1)

The package loads and its deterministic paths run without any API key;
network access is needed only for the optional API-backed adjudication path.

## 2. Full workflow (pipeline stages 03a–10)

Required to run `reproducibility/scripts/run_triage.R` end to end; the
preflight check reports all missing entries together before analysis:

- optparse, jsonlite, dplyr, readr, stringr, purrr, tibble, rlang
- digest, writexl, glue, data.table, tidyr, knitr, readxl
- httr, memoise, cachem, tictoc
- future, future.apply, furrr
- rio, AnnotationDbi, org.Hs.eg.db

## 3. Reviewer / enrichment packages (full workflow)

- **CASSIA** (tested 0.1.0) — CASSIA reviewer stage (03b); not on CRAN,
  install from its published source distribution.
- clusterProfiler (tested 4.19.4.8), DOSE (tested 4.4.0),
  ReactomePA (tested 1.52.0), enrichR (tested 3.4)
- decoupleR (tested 2.14.0) — CollecTRI TF-activity evidence (stage 05)
- disgenet2r (tested 1.2.4) — optional disease-evidence step
  (`DISGENET_API_KEY` optional)
- KEGG.db (tested 1.0; deprecated upstream but still loadable)

`06b_run_inter.R` contains an optional GitHub installation path for
`clusterProfiler` (`--install_clusterprofiler_github`).

## 4. Optional input preparation

- Seurat — only for `01a --mode seurat` (building a DEG table from a Seurat
  object). The CSV workflow does not require it.

## 5. External resources

The Cell Ontology JSON ships with the package/repo. STRING/PPI, CollecTRI
and CellMarkerDB files are not redistributed; expected paths and checksums
are documented in `resources/README.md` and `resources/CHECKSUMS.tsv`.

## 6. API / network

- A DeepSeek-compatible LLM endpoint. Pass your key and the **full**
  chat-completions endpoint URL (e.g.
  `https://api.deepseek.com/chat/completions`) directly on the command
  line with `--api-key` and `--api-base-url`; alternatively set the
  `DEEPSEEK_API_KEY` and `LLM_API_BASE_URL` environment variables. The
  endpoint contract is used as-is (no path appending).
- Internet access for the preflight reachability probe and, where used,
  resource downloads.

Supported R/package versions are summarized in `sessionInfo.txt`.
