# External resources

The exact Cell Ontology JSON used by the reported analysis is included under `inputs/raw/ontology/`. Other large third-party resources are not redistributed; place local copies under `inputs/raw/` or override paths with the documented environment variables.

Expected layout:

```text
inputs/raw/
├── ontology/
│   └── CL-ontology-v2025-07-30.json
├── ppi/
│   ├── 9606.protein.aliases.v12.0.txt
│   ├── 9606.protein.physical.links.v12.0.txt
│   ├── 10090.protein.aliases.v12.0.txt
│   └── 10090.protein.physical.links.v12.0.txt
├── collectri/
│   ├── collectri_human_network.rds
│   └── collectri_mouse_network.rds
└── cellmarker/
    ├── Cell_marker_Human.xlsx
    └── Cell_marker_Mouse.xlsx
```

CellMarkerDB file sizes and SHA-256 hashes used in the benchmark snapshot are recorded in `resources/cellmarkerdb/CHECKSUMS.tsv`.

Exact checksums of the retained resource snapshots are provided. The retained project record did not include retrieval dates for every external resource.

Exact snapshot sizes and SHA-256 hashes for Cell Ontology, STRING and CollecTRI are recorded in `resources/CHECKSUMS.tsv`.


## Cell Ontology snapshot

The exact ontology file used for ontology-dependent scoring is included:

```text
inputs/raw/ontology/CL-ontology-v2025-07-30.json
```

File size: 1,455,986 bytes  
SHA-256: `2cb74a09d7f9bc9413e638f9263463b4e519f74808a1b8e9b622772ce0d3d3ea`

This avoids dependence on a later ontology download or an undocumented conversion step when reproducing CL similarity, selector controls and CL-Linker ontology evaluation.
