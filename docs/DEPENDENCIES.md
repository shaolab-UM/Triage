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

- cellmarkeraccordion (pinned
  `TebaldiLab/cellmarkeraccordion@df19b668b26b3718bb8d760e013674ae39d811c0`
  in `Remotes`; package version label 1.0.0; the tested commit is 21
  commits after the `v1.0.0` git tag, whose target is `de85d37e…`, so the
  exact SHA is pinned rather than the tag)
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
- withr (environment restoration in `run_triage()`)

The package loads without an API key, and the deterministic Quick Start
runs offline after installation. The full workflow requires network access
for dependency and resource acquisition, literature retrieval and
LLM-backed stages.

## 2. Full workflow (pipeline stages 03a–10)

Required to run `reproducibility/scripts/run_triage.R` end to end; the
preflight check reports all missing entries together before analysis:

- optparse, jsonlite, dplyr, readr, stringr, purrr, tibble, rlang
- digest, writexl, glue, data.table, tidyr, knitr, readxl
- httr, memoise, cachem, tictoc
- future, future.apply, furrr, xml2, fs
- rio, AnnotationDbi, org.Hs.eg.db

## 3. Reviewer / enrichment packages (full workflow)

- **CASSIA** — CASSIA reviewer stage (03b); pinned to the tested revision
  `b008c0ac3dd81b2c2dff131d20f5081a58aca027`
  (`remotes::install_github("ElliotXie/CASSIA", ref = ..., subdir =
  "CASSIA_R")`), installed by `install_triage_dependencies()`. CASSIA also
  needs its Python backend (`CASSIA::setup_cassia_env()` once; the
  preflight verifies via a version-compatible lookup of CASSIA's internal
  `check_python_env()` (present in the tested revision) and compares the
  installed revision against the tested SHA).
- clusterProfiler (tested 4.19.4.8), DOSE (tested 4.4.0),
  ReactomePA (tested 1.52.0), enrichR (tested 3.4)
- decoupleR (tested 2.14.0) — CollecTRI TF-activity evidence (stage 05)
- **reactome.db** (tested 1.92.0) — required for the canonical Reactome
  enrichment dimension (stage 05); current Bioconductor release, installed
  by `install_triage_dependencies()` and checked by both preflights
  (without it stage 05 skips Reactome evidence).
- **KEGG.db** — required for the canonical KEGG evidence dimension:
  stage 05 calls `clusterProfiler::enrichKEGG(use_internal_data = TRUE)`,
  which is backed by the KEGG.db data package and stops without it (the
  stage's error handling would silently drop KEGG evidence). KEGG.db was
  removed from Bioconductor with release 3.11 ("use KEGGREST instead");
  the verified install route is the archived source tarball
  `https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz`,
  installed automatically by `install_triage_dependencies()` and checked
  (including a data-load check) by both preflights.
- disgenet2r (tested 1.2.4) — optional disease-evidence step, required
  only when `DISGENET_API_KEY` is set
  (`remotes::install_gitlab("medbio/disgenet2r")`)
- R.utils — used by `setup_triage_resources()` to decompress the STRING
  archives; installed by `install_triage_dependencies()`.
- **fanyi (tested 0.1.0)** — required for stage 06b: the enrichment
  reviewer's `clusterProfiler::interpret()` call routes its LLM transport
  through `fanyi::chat_request()`. fanyi hard-codes
  `https://api.deepseek.com/v1/chat/completions` and has no base-url
  parameter (0.1.1 is identical), so stage 06b overrides the endpoint with
  the canonical `LLM_API_BASE_URL` value via a narrow, self-restoring
  namespace shim (`Triage:::.triage_interpret()`). fanyi is pinned to the
  tested version (`remotes::install_version("fanyi", version = "0.1.0")`),
  installed by `install_triage_dependencies()` and checked by both
  preflights.

`06b_run_inter.R` contains an optional GitHub installation path for
`clusterProfiler` (`--install_clusterprofiler_github`), pinned to the
tested revisions (GOSemSim then clusterProfiler).

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
