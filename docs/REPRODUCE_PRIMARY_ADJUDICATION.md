# Primary adjudication reproducibility

Two reproducibility levels are provided.

## 1. Exact API-facing payload archive

The exact 50 round-1 Handling Editor payloads from the reported primary run are stored under:

```text
data/primary/<dataset>/handling_editor_round1_inputs/
```

They contain the final system/user prompt material, reviewer inputs, biological dossier and ontology context sent into the Handling Editor stage. Their SHA-256 hashes are listed in `data/primary/handling_editor_round1_manifest.tsv`.

These files are provenance records. They are not presented as the pre-`09_run_judge.R` temporary input directory.

## 2. Regenerating adjudication inputs with the released workflow

The formal primary run created temporary judge inputs with `08_build_judge_inputs.R`, then passed them to `09_run_judge.R`. The temporary `/tmp/...` directory was not retained in the project snapshot. The reviewer mapping registry used for the primary datasets is included at:

```text
data/primary/reviewer_mapping_registry.tsv
```

The main pipeline scripts are provided so that reviewer outputs and judge inputs can be regenerated from the released workflow and required external resources.

The formal primary adjudication settings were:

```text
Handling Editor: deepseek-v4-flash
Chief QC: deepseek-v4-flash
temperature: 0
max_rounds: 3
always_run_chief: TRUE
gate_mode: fixed_k
gate_k: 2
delta_depth: 2
identity_policy: ontology_guarded
parity_mode: TRUE
enable_aggressive_lca_hard_guard: TRUE
enable_method_reliability_gate: FALSE
meta_reviewer_gate: FALSE
ols_first: FALSE
release_policy: auto
```

The machine-readable profile is `config/primary_adjudication_profile.tsv`.

LLM/API execution may not be bitwise deterministic even at temperature 0. Publication-facing final outputs are therefore included separately under `data/primary/<dataset>/final/`.
