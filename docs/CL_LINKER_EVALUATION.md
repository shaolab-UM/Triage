# CL-Linker evaluation

The release includes `reproducibility/scripts/analysis/evaluate_cl_linker.R`.

It checks the two evaluation layers described in the manuscript:

1. **Direct mapping on 97 unique single-identity labels.** The script recomputes the OLS top-result, OLS exact and exact/synonym comparator outcome counts from archived per-label predictions. The publication-final CL-Linker row is retained as source data in `reproducibility/cl_linker/publication_direct_mapping_summary.tsv`.
2. **Operational mapping on 140 reviewer annotations.** This is recomputed directly from `reproducibility/primary/reviewer_mapping_registry.tsv` and `reproducibility/cl_linker/manual_mapping_gold.tsv`, including coverage, exact/direct-parent mappings, overspecific mappings, incorrect mappings and abstentions. The ten mixed/ambiguous reference cases are also checked for single-ID withholding.

The separate per-label output from the 97-label CL-Linker direct evaluation was not present in the project snapshot supplied for release assembly. The repository therefore does not fabricate that missing file. The final direct-mapping CL-Linker counts and paired-comparison statistics are distributed as publication source data.

Run:

```bash
Rscript reproducibility/scripts/analysis/evaluate_cl_linker.R   --repo-root /path/to/Triage   --cl-json /path/to/CL-ontology-v2025-07-30.json
```

Expected publication values include:

- 97-label CL-Linker exact mapping: 63/97 (64.9%).
- 97-label CL-Linker exact or direct parent: 64/97 (66.0%).
- Operational set: 115/140 mapped (82.1%).
- Operational exact mapping: 102/115 released mappings.
- Operational exact or direct-parent mapping: 104/115 released mappings (90.4%).


The final per-label CL-Linker prediction file for the 97-label direct evaluation was not present in the retained project record. No replacement per-label CL-Linker prediction file is constructed for this release; the reported direct-mapping row remains publication source data.
