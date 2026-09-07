# Model configuration

## Primary benchmark and main Triage workflow

The reported primary workflow used DeepSeek V4 Flash throughout:

- CASSIA reviewer: `deepseek-v4-flash`
- In-house reviewer: `deepseek-v4-flash`
- clusterProfiler reviewer: `deepseek-v4-flash`
- Handling Editor: `deepseek-v4-flash`
- Chief QC: `deepseek-v4-flash`
- Handling Editor and Chief QC temperature: 0

Machine-readable records:

- `reproducibility/config/primary_adjudication_profile.tsv`
- `reproducibility/config/model_run_provenance.tsv`

## Separate profiling / model-comparison experiment

A separate computational experiment (reported independently of the primary
benchmark in the manuscript) used:

- Handling Editor: `deepseek-reasoner`
- Chief QC: `deepseek-chat`

Alternative model-comparison conditions each used one model for both the
Handling Editor and Chief QC: GPT-5.4, Gemini-3 Pro Preview and Claude Opus
4.5. These identifiers describe only the profiling / model-comparison
experiment and do not define the primary workflow; no profiling result values
are affected by this repository release.

## Endpoint and provenance

The reported runs used a provider-hosted DeepSeek-compatible OpenAI Chat
Completions endpoint. The reviewer mapping registry
(`reproducibility/primary/reviewer_mapping_registry.tsv`) records
`deepseek-v4-flash` for the primary reviewer runs.

Model identifiers describe the configuration used for the reported
analysis. A new execution should use model identifiers accepted by the
configured provider and should document any substitution as a new
execution configuration.
