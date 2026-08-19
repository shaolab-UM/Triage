# Masked runtime values

`XXXXX` denotes a value that must be supplied by the user or deployment environment.

The public release masks runtime-specific values such as:

- API credentials;
- LLM provider base URL;
- local repository/resource paths in `.env.example`.

Scientific/reproducibility constants are intentionally **not** masked. These include model identifiers reported in the manuscript, Cell Ontology release, thresholds, dataset sizes, prompt profiles and benchmark values.

Before running the workflow, set the corresponding environment variables. In particular:

```text
DEEPSEEK_API_KEY=XXXXX
DISGENET_API_KEY=XXXXX
LLM_API_BASE_URL=XXXXX
TRIAGE_HOME=XXXXX
PROJECT_ROOT=XXXXX
CL_LOCAL_JSON=XXXXX
TRIAGE_PPI_ROOT=XXXXX
```
