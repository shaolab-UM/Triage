# Historical model configuration

The reported primary benchmark used:

- Handling Editor: `deepseek-v4-flash`
- Chief QC: `deepseek-chat`
- Temperature: 0

`deepseek-chat` is retained as the historical Chief QC model identifier used by
the reported analysis. Runtime availability of a historical model name depends
on the API provider or compatible gateway used for a new execution. If a new
execution substitutes a different model identifier, that execution should be
documented as a new configuration and should not be presented as an identical
historical rerun.

## Provider and run provenance

The retained source code used a provider-hosted DeepSeek-compatible OpenAI Chat Completions API endpoint. The exact endpoint string and all credentials are intentionally masked from the public code.

For the reported primary adjudication, the retained configuration records:

- Handling Editor model identifier: `deepseek-v4-flash`
- Chief QC model identifier: `deepseek-chat`
- Temperature: 0
- Provider category: DeepSeek-compatible OpenAI Chat Completions API
- Retained formal artifact timestamps: August 2026

These model identifiers describe the historical analysis configuration. A new execution should use model identifiers accepted by the configured provider and should document any substitution as a new execution configuration.
