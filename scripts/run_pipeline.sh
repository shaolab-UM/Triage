#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Triage core workflow runner

Starts from a masked cluster-level DEG CSV and runs the public Triage workflow:
03a -> 03b -> 04a -> 05 -> 06 -> 06b -> 07 -> 07b -> 07.5 -> 08 -> 09 -> 10
Step 11 evaluation is optional when a reference-label CSV is supplied.

Usage:
  bash scripts/run_pipeline.sh \
    --dataset Census_immune \
    --masked-deg data/primary/Census_immune/model_inputs/maskdeg.csv \
    [--run-tag YYYYMMDD_HHMMSS] \
    [--workers 4] \
    [--true-label-csv path/to/true_label.csv]

Required runtime environment:
  DEEPSEEK_API_KEY
  LLM_API_BASE_URL
  TRIAGE_HOME (optional; inferred from this script if unset)
  CL_LOCAL_JSON (optional; defaults to inputs/raw/ontology/CL-ontology-v2025-07-30.json)

Notes:
  - This wrapper intentionally starts from anonymized/masked DEG input.
  - 01a/02a/02b remain available for benchmark input preparation but are not
    automatically run here, keeping evaluation labels outside the adjudication path.
  - External Cell Ontology, STRING, CollecTRI and CellMarkerDB resources must be
    placed under inputs/raw/ as documented in resources/README.md.
EOF
}

DATASET=""
MASKED_DEG=""
RUN_TAG=""
WORKERS=4
TRUE_LABEL=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset) DATASET="$2"; shift 2 ;;
    --masked-deg) MASKED_DEG="$2"; shift 2 ;;
    --run-tag) RUN_TAG="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --true-label-csv) TRUE_LABEL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$DATASET" ]] || { echo "--dataset is required" >&2; exit 2; }
[[ -n "$MASKED_DEG" ]] || { echo "--masked-deg is required" >&2; exit 2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TRIAGE_HOME="${TRIAGE_HOME:-$(cd "$SCRIPT_DIR/.." && pwd)}"
export PROJECT_ROOT="${PROJECT_ROOT:-$TRIAGE_HOME}"

MASKED_DEG="$(python3 - <<'PY' "$MASKED_DEG"
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
)"
[[ -f "$MASKED_DEG" ]] || { echo "Masked DEG file not found: $MASKED_DEG" >&2; exit 1; }

export CL_LOCAL_JSON="${CL_LOCAL_JSON:-$TRIAGE_HOME/inputs/raw/ontology/CL-ontology-v2025-07-30.json}"
[[ -f "$CL_LOCAL_JSON" ]] || {
  echo "Cell Ontology JSON not found: $CL_LOCAL_JSON" >&2
  echo "See resources/README.md." >&2
  exit 1
}

[[ -n "${DEEPSEEK_API_KEY:-}" && "${DEEPSEEK_API_KEY:-}" != "XXXXX" ]] || {
  echo "Set DEEPSEEK_API_KEY before running the API stages." >&2
  exit 1
}
[[ -n "${LLM_API_BASE_URL:-}" && "${LLM_API_BASE_URL:-}" != "XXXXX" ]] || {
  echo "Set LLM_API_BASE_URL before running the API stages." >&2
  exit 1
}

if [[ -z "$RUN_TAG" ]]; then RUN_TAG="$(date +%Y%m%d_%H%M%S)"; fi
RUN_DIR="$TRIAGE_HOME/outputs/$DATASET/$RUN_TAG"
mkdir -p "$RUN_DIR/logs"

cfg_value() {
  local key="$1"
  Rscript -e '
    args <- commandArgs(trailingOnly=TRUE)
    source(file.path(args[2], "config", "dataset_config.R"))
    x <- get_dataset_config(args[1], args[2])
    v <- x[[args[3]]]
    if (is.null(v) || length(v) == 0) quit(status=3)
    cat(as.character(v[[1]]))
  ' "$DATASET" "$TRIAGE_HOME" "$key"
}

SPECIES="$(cfg_value species)"
TISSUE="$(cfg_value tissue)"
STUDY_CONTEXT="$(cfg_value study_context)"

run_step() {
  local name="$1"; shift
  echo "[$(date '+%F %T')] START $name"
  "$@" 2>&1 | tee "$RUN_DIR/logs/${name}.log"
  echo "[$(date '+%F %T')] DONE  $name"
}

# 03a: deterministic DEG filtering
mkdir -p "$RUN_DIR/03a_filtered"
run_step 03a \
  Rscript "$TRIAGE_HOME/scripts/pipeline/03a_filter_deg.R" \
  --deg "$MASKED_DEG" \
  --out_dir "$RUN_DIR/03a_filtered"
FILTERED="$RUN_DIR/03a_filtered/filtered_deg.csv"
[[ -f "$FILTERED" ]] || { echo "03a output missing: $FILTERED" >&2; exit 1; }

# 03b: CASSIA reviewer
mkdir -p "$RUN_DIR/03b_cassia"
run_step 03b \
  Rscript "$TRIAGE_HOME/scripts/pipeline/03b_run_cassia.R" \
  --deg "$FILTERED" \
  --out_dir "$RUN_DIR/03b_cassia" \
  --tissue "$TISSUE" \
  --species "$SPECIES" \
  --study_context "$STUDY_CONTEXT" \
  --workers "$WORKERS" \
  --api_base_url "$LLM_API_BASE_URL"
mapfile -t CASSIA_FILES < <(find "$RUN_DIR/03b_cassia" -type f -path '*/01_annotation_results/annotation_cassia_FINAL_RESULTS.csv' | sort)
[[ "${#CASSIA_FILES[@]}" -eq 1 ]] || {
  echo "Expected exactly one CASSIA final CSV; found ${#CASSIA_FILES[@]}" >&2
  exit 1
}
CASSIA_CSV="${CASSIA_FILES[0]}"

# 04a: candidate generation
mkdir -p "$RUN_DIR/04a_candidates"
run_step 04a \
  Rscript "$TRIAGE_HOME/scripts/pipeline/04a_build_candidates.R" \
  --deg "$FILTERED" \
  --out_dir "$RUN_DIR/04a_candidates" \
  --ontology "$CL_LOCAL_JSON" \
  --species "$SPECIES"
CANDIDATES="$RUN_DIR/04a_candidates/cellanno/structured_candidates.csv"
[[ -f "$CANDIDATES" ]] || { echo "04a candidate file missing: $CANDIDATES" >&2; exit 1; }

# 05: biological dossier + LLM query construction.
# llm_run.R writes intermediate_outputs relative to the working directory,
# therefore run this stage from RUN_DIR.
mkdir -p "$RUN_DIR/05c_llm_queries"
(
  cd "$RUN_DIR"
  run_step 05 \
    Rscript "$TRIAGE_HOME/scripts/pipeline/05_prepare_llm_inputs.R" \
    --mode full \
    --deg_file "$FILTERED" \
    --candidates_file "$CANDIDATES" \
    --out_root "$RUN_DIR/05c_llm_queries" \
    --dataset_name "$DATASET" \
    --project_root "$TRIAGE_HOME"
)

# 06: In-house reviewer
mkdir -p "$RUN_DIR/06_llm_outputs"
run_step 06 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/06_run_llm_pipeline.R" \
  --mode full \
  --input_root "$RUN_DIR/05c_llm_queries" \
  --output_root "$RUN_DIR/06_llm_outputs" \
  --dataset_name "$DATASET" \
  --project_root "$TRIAGE_HOME" \
  --workers "$WORKERS" \
  --max_rounds 3 \
  --cl_local_json "$CL_LOCAL_JSON"

# 06b: enrichment / clusterProfiler reviewer
mkdir -p "$RUN_DIR/06b_inter"
BIOINFO="$RUN_DIR/intermediate_outputs/${DATASET}_LLM_Input_Run/bioinformatics_tsv"
run_step 06b \
  Rscript "$TRIAGE_HOME/scripts/pipeline/06b_run_inter.R" \
  --marker_csv "$FILTERED" \
  --step1_dir "$RUN_DIR/05c_llm_queries/step1_report_queries" \
  --bioinfo_dir "$BIOINFO" \
  --out_dir "$RUN_DIR/06b_inter" \
  --dataset_name "$DATASET" \
  --species "$SPECIES" \
  --workers "$WORKERS" \
  --cl_local_json "$CL_LOCAL_JSON"

# 07 / 07b: compact reviewer summaries
mkdir -p "$RUN_DIR/07_our_summary" "$RUN_DIR/07b_inter_summary"
run_step 07 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/07_our_llm_summary.R" \
  --out_root "$RUN_DIR/06_llm_outputs" \
  --dataset_name "$DATASET" \
  --out_dir "$RUN_DIR/07_our_summary"

run_step 07b \
  Rscript "$TRIAGE_HOME/scripts/pipeline/07b_inter_summary.R" \
  --in_dir "$RUN_DIR/06b_inter/final_passed" \
  --out_dir "$RUN_DIR/07b_inter_summary" \
  --dataset_name "$DATASET"

# 07.5: CL-Linker reviewer mapping registry for this run.
MANIFEST="$RUN_DIR/evidence_mapping_run.tsv"
printf 'dataset\trun_dir\n%s\t%s\n' "$DATASET" "$RUN_TAG" > "$MANIFEST"
mkdir -p "$RUN_DIR/07.5_mapping"
run_step 07_5 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/07.5_build_reviewer_mapping.R" \
  --cl_json "$CL_LOCAL_JSON" \
  --out_dir "$RUN_DIR/07.5_mapping" \
  --manifest "$MANIFEST" \
  --dataset "$DATASET"

mapfile -t REGISTRY_FILES < <(find "$RUN_DIR/07.5_mapping" -maxdepth 1 -type f -name 'reviewer_mapping_registry*.tsv' | sort)
[[ "${#REGISTRY_FILES[@]}" -eq 1 ]] || {
  echo "Expected exactly one mapping registry; found ${#REGISTRY_FILES[@]}" >&2
  exit 1
}
REGISTRY="${REGISTRY_FILES[0]}"

# 08: assemble adjudication input records
mkdir -p "$RUN_DIR/08_judge_inputs"
run_step 08 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/08_build_judge_inputs.R" \
  --cassia_csv "$CASSIA_CSV" \
  --in_house_summary_csv "$RUN_DIR/07_our_summary/summary.csv" \
  --enrichment_summary_csv "$RUN_DIR/07b_inter_summary/summary.csv" \
  --mapping_registry_csv "$REGISTRY" \
  --intermediate_outputs_dir "$RUN_DIR/intermediate_outputs/${DATASET}_LLM_Input_Run" \
  --step1_dir "$RUN_DIR/05c_llm_queries/step1_report_queries" \
  --out_dir "$RUN_DIR/08_judge_inputs" \
  --dataset_name "$DATASET"

# 09: Handling Editor + Chief QC
mkdir -p "$RUN_DIR/09_judge_outputs"
run_step 09 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/09_run_judge.R" \
  --dataset_name "$DATASET" \
  --judge_input_dir "$RUN_DIR/08_judge_inputs" \
  --out_root "$RUN_DIR/09_judge_outputs" \
  --model_head "${TRIAGE_MODEL_HEAD:-deepseek-v4-flash}" \
  --model_chief "${TRIAGE_MODEL_CHIEF:-deepseek-chat}" \
  --temperature 0 \
  --workers "$WORKERS" \
  --max_rounds 3 \
  --always_run_chief \
  --release_policy auto \
  --ols_first FALSE

# 10: publication-facing post summary
run_step 10 \
  Rscript "$TRIAGE_HOME/scripts/pipeline/10_judge_post_summary.R" \
  --out_root "$RUN_DIR/09_judge_outputs" \
  --dataset_name "$DATASET"

# 11: optional evaluation; reference data are introduced only here.
if [[ -n "$TRUE_LABEL" ]]; then
  TRUE_LABEL="$(python3 - <<'PY' "$TRUE_LABEL"
import os,sys
print(os.path.abspath(sys.argv[1]))
PY
)"
  [[ -f "$TRUE_LABEL" ]] || { echo "Reference label CSV not found: $TRUE_LABEL" >&2; exit 1; }
  mkdir -p "$RUN_DIR/11_eval"
  run_step 11 \
    Rscript "$TRIAGE_HOME/scripts/pipeline/11_eval_accuracy.R" \
    --dataset_name "$DATASET" \
    --true_label_csv "$TRUE_LABEL" \
    --cassia_csv "$CASSIA_CSV" \
    --our_csv "$RUN_DIR/07_our_summary/summary.csv" \
    --inter_csv "$RUN_DIR/07b_inter_summary/summary.csv" \
    --judge_csv "$RUN_DIR/09_judge_outputs/summary_final.csv" \
    --judge_final_dir "$RUN_DIR/09_judge_outputs/final" \
    --out_dir "$RUN_DIR/11_eval" \
    --cl_json "$CL_LOCAL_JSON"
fi

echo
echo "Triage workflow completed successfully."
echo "Run directory: $RUN_DIR"
echo "Final adjudication: $RUN_DIR/09_judge_outputs/final"
echo "Summary: $RUN_DIR/09_judge_outputs/summary_final.csv"
