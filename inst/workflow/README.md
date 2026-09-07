# Triage installed-package workflow runtime

This directory bundles the validated full-workflow pipeline scripts so that
the installed `Triage` package can run the complete DEG-input workflow
(`run_triage()`) without a repository checkout or shell access.

## Source relationship

Every file here is a byte-identical copy of the manuscript-reproduction
script of the same name in the repository tree:

| Installed file | Repository source |
|---|---|
| `llm_run.R` | `reproducibility/scripts/llm_run.R` |
| `pipeline/<stage>.R` | `reproducibility/scripts/pipeline/<stage>.R` |

The repository `reproducibility/` tree remains the single scientific
reference (manuscript reproduction layer); `inst/workflow/` is a thin
execution layer for installed-package users. No scientific logic is
duplicated in a divergent form — the copies are exact, and a package test
(`tests/testthat/test-workflow-runtime.R`) asserts byte identity between
the two locations whenever the repository tree is available. To update the
installed runtime, re-copy the repository scripts into `inst/workflow/`.

The only runtime difference is environmental: `run_triage()` sets
`TRIAGE_WORKFLOW_DIR` to this directory so that stage 05 sources the
bundled `llm_run.R`, and it stages external resources (Cell Ontology,
CollecTRI, CellMarkerDB) from the Triage user-data directory
(`tools::R_user_dir("Triage", "data")`) into a per-run runtime home so the
validated default path resolution applies unchanged.
