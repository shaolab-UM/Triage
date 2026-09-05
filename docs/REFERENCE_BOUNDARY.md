# Reference-label boundary

The public release separates model-facing records from benchmark references.

## Model-facing primary inputs

Each primary dataset contains:

- `model_inputs/maskdeg.csv`: differential-expression input with anonymized `cluster_N` identifiers.
- `final/cluster_*.json`: publication-facing final adjudication records.

Internal round-1 Handling Editor payloads are internal adjudication
inputs and are not redistributed in the public release. The release
checksum manifest (`MANIFEST_SHA256.tsv`) covers all released files.

The system prompt states that reviewer hypotheses are "not ground
truth"; this methodological phrase is not a benchmark-label input.

## Evaluation-only files

Each primary dataset stores benchmark references under `evaluation/`:

- `true_label.csv`
- `cluster_map.csv`
- `reference_cl.tsv`

These files are not part of the Handling Editor payloads.

In the pipeline source, reference labels are created during anonymization (`02a_build_true_label.R`), used to create anonymous cluster identifiers (`02b_build_maskdeg.R`) and read for scoring by `11_eval_accuracy.R`. Adjudication operates on masked clusters and reviewer/evidence inputs. Selector controls use reference CL identifiers only for retrospective scoring.
