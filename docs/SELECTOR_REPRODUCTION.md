# Selector-control reproduction

This release reconstructs the deterministic selector controls reported in Table S10 and Fig. 3C-D without LLM or API calls.

## Inputs

`data/primary/selector_inputs.tsv` contains the 50 primary benchmark clusters.

Native score sources are:

- CASSIA: the native score stored in the final judge input.
- In-house: the reviewer diagnostic method score from the final adjudication record.
- clusterProfiler: confidence from the final enrichment reviewer summary.

The three reviewer score streams are percentile-normalized independently using average ranks.

## Controls

`reproduce_selector_controls.R` evaluates:

1. Majority vote.
2. Reviewer selection after method-specific percentile normalization.
3. The deterministic ontology-only control.
4. Triage.
5. A retrospective best-available oracle.

The retrospective oracle is reference-using and is computed cluster by cluster as the maximum CL similarity across **CASSIA, In-house, clusterProfiler and Triage**, additionally including the publication's label-mapped retrospective candidates (mapped through Table S4; three clusters gain a stronger candidate). It is therefore a descriptive retrospective ceiling, not an independently deployable selector.

Deferred selector outputs receive CL similarity 0, matching the publication analysis.

## Expected overall values

- Majority vote: 41.3207
- Top reviewer after percentile normalization: 61.6773
- Ontology-only control: 62.9211
- Triage: 79.9790
- Retrospective oracle: 84.9083

Expected unavailable/deferred counts are 24/50 for majority vote, 9/50 for reviewer-score selection and 5/50 for the ontology-only control.

`data/primary/selector_expected_tableS10.tsv` contains the publication summary used for the final check.
