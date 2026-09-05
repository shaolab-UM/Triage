# CellMarkerDB resource

The benchmark used locally downloaded CellMarkerDB spreadsheets:

- `Cell_marker_Human.xlsx`
- `Cell_marker_Mouse.xlsx`

The files were read by `reproducibility/scripts/pipeline/06b_run_inter.R`, where the `cell_name` and `marker` columns are supplied to `clusterProfiler::enricher` as `TERM2GENE`.

For release hygiene, the third-party XLSX files are not redistributed in this repository. Place local copies at:

```text
inputs/raw/cellmarker/Cell_marker_Human.xlsx
inputs/raw/cellmarker/Cell_marker_Mouse.xlsx
```

`CHECKSUMS.tsv` records the exact file sizes and SHA-256 hashes of the copies used in the benchmark snapshot.

## Baseline reviewer versus sensitivity reviewer

CellMarkerDB has two distinct roles in the analysis:

1. The baseline `clusterProfiler` reviewer uses CellMarkerDB as one gene-set source in `06b_run_inter.R`.
2. The Fig. 3A sensitivity analysis adds an additional CellMarkerDB-derived ORA assignment as a non-core fourth reviewer, with an unfiltered condition and a reliability-filtered condition.

These roles should not be conflated when describing the workflow.
