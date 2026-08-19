# Configuration notes

`dataset_config.R` preserves the configuration used by the analysis, including internal dataset aliases for the three validation datasets. Publication-facing data directories use the manuscript names:

- `screview_clean` -> `Anderson_DLPFC`
- `screview_mouse` -> `Zha_AD_mouse`
- `screview_ture` -> `S6K1_organoid`

These aliases are retained in configuration because changing model-facing context after the analysis would no longer reproduce the reported execution.
