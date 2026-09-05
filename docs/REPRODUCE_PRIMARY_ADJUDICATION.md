# Primary adjudication reproducibility

## Reproducibility boundary

The public release provides:

- **Derived benchmark inputs** — masked DEG model inputs, evaluation files
  and reviewer mapping registry under `reproducibility/primary/`.
- **Configuration** — `reproducibility/config/primary_adjudication_profile.tsv`
  and `reproducibility/config/model_run_provenance.tsv`.
- **Pipeline scripts** — the full workflow under
  `reproducibility/scripts/pipeline/`.
- **Publication-facing final outputs** — the 50 final adjudication JSONs
  under `reproducibility/primary/<dataset>/final/`.

Re-executing the model stages (Handling Editor, Chief QC) requires
external API access. LLM execution may not be bitwise deterministic even
at temperature 0, and identical generated text should not be expected;
the released final outputs remain the publication-facing record.

## Formal primary adjudication settings

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

The machine-readable profile is `reproducibility/config/primary_adjudication_profile.tsv`.

## Regenerating adjudication inputs with the released workflow

The pipeline creates judge inputs with `08_build_judge_inputs.R`, then
passes them to `09_run_judge.R`. The reviewer mapping registry used for
the primary datasets is included at:

```text
reproducibility/primary/reviewer_mapping_registry.tsv
```

The pipeline scripts regenerate reviewer outputs and judge inputs from
the released workflow and the required external resources.
