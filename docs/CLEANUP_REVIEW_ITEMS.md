# Cleanup review items

Unresolved publication-alignment items identified during the structural
cleanup of 2026-09-05. None of these were changed in the cleanup; they
require separate scientific verification before any value or numbering
is altered.

## 1. Selector-control table numbering may be stale

`data/primary/selector_expected_tableS7.tsv` is named after the old
supplementary table numbering ("Table S7"). Current manuscript table
numbering may place the selector controls at Table S10. The file
contents (including all numerical values) are unchanged in this
cleanup; renaming or renumbering requires confirmation against the
final manuscript table numbering.

## 2. Retrospective oracle value in the released selector package

The released selector material carries a retrospective oracle value of
`83.7732` (e.g. `README.md`, `docs/SELECTOR_REPRODUCTION.md`,
`scripts/analysis/reproduce_selector_controls.R` and
`data/primary/selector_expected_tableS7.tsv`). Later manuscript work
may use a different retrospective oracle value. The value was left
exactly as released; resolving the discrepancy requires scientific
verification, not a structural cleanup.

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
