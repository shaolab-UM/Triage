# Primary benchmark public adjudication records

The primary benchmark contains 50 clusters across five datasets:

- Census immune: 16
- Sikkema lung: 11
- Tabula Sapiens kidney: 4
- Tabula Sapiens pancreas: 10
- Zheng blood: 9

Each `final/*.json` file is a publication-facing record. Reviewer and execution metadata are taken from the retained formal adjudication records. Final decision, evidence, QC disposition and reference-evaluation fields are synchronized to `P3_Primary_Adjudication` in the submitted Supplementary Data.

These files are intentionally not raw debugging traces. Internal decision-trace, fallback, shadow-pool, local-path and timestamp fields are omitted from the public release.
