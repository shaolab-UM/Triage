# Data layout

Each reported dataset has a publication-facing `final/` directory. Individual cleaned JSON records are stored one file per cluster. This preserves the dataset/cluster organization used by the workflow while avoiding publication of internal engineering traces.

## Primary benchmark

- Census_immune: 16 clusters
- Sikkema_lung: 11 clusters
- TS_kidney: 4 clusters
- TS_pancreas: 10 clusters
- Zheng_blood: 9 clusters

Expected total: 50 cluster JSON files.

## External validation

- Anderson_DLPFC: 18 clusters
- Zha_AD_mouse: 9 clusters
- S6K1_organoid: 15 clusters

Expected total: 42 cluster JSON files.


## Publication-facing JSON records

The released JSON files are retained publication-facing records. The historical transformation step used to prepare these cleaned public records is not included as an executable release step. The public package therefore validates the retained records directly with `reproducibility/scripts/release/check_public_json_counts.R` and `reproducibility/scripts/release/validate_public_release.R`.


## Reviewer metadata in publication-facing JSON

Reviewer objects retain source-facing information: the reviewer label, reviewer
Cell Ontology assignment and model-generated strengths/weaknesses.

The following reviewer-to-final derived fields are intentionally omitted:

- `final_cl_id`
- `support_class`
- `is_correct`
- `matches_final_label`

They are defined relative to an adjudication-state final decision and are not
required for the published benchmark endpoints. Because the public final record
is a publication-facing synchronized record, retaining adjudication-relative
derived fields can mix two decision states in one JSON object.

`final_decision` is therefore the single authoritative publication-facing final
state.

### Reviewer CL fields

The two CL fields in a reviewer object have different provenance:

- `reviewer_cl_id` preserves the CL assignment attached to the incoming reviewer
  summary before later adjudication-side normalization or refinement.
- `cell_ontology_id` is the CL assignment stored in the Handling Editor's
  reviewer-specific method verdict and may therefore differ from
  `reviewer_cl_id`.

These fields should not be used to reconstruct the selector controls. The
authoritative frozen reviewer CL IDs used by the selector analyses are provided
in `reproducibility/selector_controls/selector_inputs.tsv`.
