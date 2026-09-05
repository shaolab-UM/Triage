# Cleanup review items

Publication-alignment items identified during the structural cleanup
of 2026-09-05. Items 1 and 2 were resolved on 2026-09-05 against the
latest manuscript (5th revision, Tables renumbered: sensitivity =
Table S9, selector controls = Table S10). Item 3 remains open.

## 1. Selector-control table numbering — RESOLVED

The expected-values file is now
`data/primary/selector_expected_tableS10.tsv`, matching the final
manuscript numbering (selector controls = Table S10, Fig. 3C-D).
References in `README.md`, `docs/SELECTOR_REPRODUCTION.md` and
`scripts/analysis/reproduce_selector_controls.R` were updated. The
CellMarkerDB sensitivity documentation now cites Table S9 (the
sensitivity table in the final numbering).

## 2. Retrospective oracle value — RESOLVED

The manuscript (Table S10) reports a retrospective selector overall of
`84.9083` (Unresolved disagreement group `87.518`), replacing the
released `83.7732` / `84.531`. The change comes from three clusters
where the publication's retrospective selector used comparator output
labels mapped to Cell Ontology IDs (Table S4 mapping):

- Census_immune/cluster_5: clusterProfiler -> CL:0000910 (similarity 100)
- TS_pancreas/cluster_10: CASSIA -> CL:0008024 (similarity 87.6180)
- TS_pancreas/cluster_5: In-house -> CL:0000125 (similarity 4.0068)

`data/primary/selector_inputs.tsv` retrospective columns were updated
to the published values, `reproduce_selector_controls.R` now includes
these published label-mapped candidates in the retrospective ceiling,
and all consumer documents/scripts report `84.9083`. The reproduction
script passes against the updated expected values (overall
84.908337..., matching the manuscript exactly).

## 3. No R-package skeleton

The repository has no `DESCRIPTION`, `NAMESPACE` or `R/` directory. If
the manuscript positions Triage as an installable R package, a package
skeleton (or an explicit statement that the release is script-based)
is a publication-positioning issue to address separately. Not
addressed in this structural cleanup.

## Notes

- A full sweep of tracked text files found no references to the old
  root-level data directories; documentation, scripts and release
  checks already expected the `data/...` layout, so no path rewrites
  were required beyond the directory moves themselves.
