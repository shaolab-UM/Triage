# Reference-label boundary

The public release separates model-facing records from benchmark references.

## Model-facing primary inputs

Each primary dataset contains:

- `model_inputs/maskdeg.csv`: differential-expression input with anonymized `cluster_N` identifiers.
- `handling_editor_round1_inputs/*_round1.json`: the exact API-facing round-1 Handling Editor payload retained from the formal primary run.

The 50 formal round-1 payloads were scanned during release assembly for benchmark/reference-like **field names**, local absolute paths and literal credentials. No such fields or local/secret values were detected. The system prompt itself states that reviewer hypotheses are “not ground truth”; this methodological phrase is not a benchmark-label input. Checksums and scan status are recorded in `data/primary/handling_editor_round1_manifest.tsv`.

Older intermediate judge-input files are not published as formal inputs because they are not byte-identical to the retained API-facing payloads from the reported run.

## Evaluation-only files

Each primary dataset stores benchmark references under `evaluation/`:

- `true_label.csv`
- `cluster_map.csv`
- `reference_cl.tsv`

These files are not part of the Handling Editor payloads.

In the pipeline source, reference labels are created during anonymization (`02a_build_true_label.R`), used to create anonymous cluster identifiers (`02b_build_maskdeg.R`) and read for scoring by `11_eval_accuracy.R`. Adjudication operates on masked clusters and reviewer/evidence inputs. Selector controls use reference CL identifiers only for retrospective scoring.
