# =========================================================================
# adjudication.R — deterministic adjudication core
# Adapted from scripts/pipeline/09_run_judge.R (v1.0.0 release).
# RULES, schema, prompts, normalization, release policy, deterministic adjudication rules,
# local gate. API invocation lives in provider_deepseek.R.
# =========================================================================

# Prompt profile used by the primary publication pipeline
# (config/primary_adjudication_profile.tsv: prompt_profile = compact).
PROMPT_PROFILE <- "compact"

# Pipeline verbosity flag (was a command-line option in the release pipeline)
VERBOSE <- FALSE

# Dataset configuration object (set by reproducibility pipeline scripts via
# config/dataset_config.R; NULL in standalone package use).
cfg <- NULL

vlog <- function(...) {
  if (isTRUE(VERBOSE)) message(...)
}
vcat <- function(...) {
  if (isTRUE(VERBOSE)) cat(...)
}
# Try to find script directory
pr <- "."
triageHome09 <- Sys.getenv("TRIAGE_HOME", unset = "")
if (!nzchar(triageHome09)) {
  scriptArgV09 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  triageHome09 <- if (length(scriptArgV09) > 0) {
    dirname(dirname(dirname(normalizePath(sub("^--file=", "", scriptArgV09[[1]]), winslash = "/", mustWork = FALSE))))
  } else getwd()
}
# dataset_config and cl_normalizer logic ships with the package namespace.
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

MARKER_CACHE <- new.env(parent = emptyenv())

normalize_gene_list <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  toupper(unique(x))
}

get_marker_map <- function(species_value) {
  key <- if (tolower(species_value %||% "human") == "mouse") "Mouse" else "Human"
  if (exists(key, envir = MARKER_CACHE, inherits = FALSE)) return(get(key, envir = MARKER_CACHE))
  data(accordion_marker, package = "cellmarkeraccordion")
  am_dt <- as.data.frame(accordion_marker, stringsAsFactors = FALSE)
  am_dt <- am_dt[am_dt$marker_type == "positive" & am_dt$species == key, , drop = FALSE]
  if (nrow(am_dt) == 0) {
    assign(key, list(), envir = MARKER_CACHE)
    return(list())
  }
  am_dt$CL_celltype <- tolower(trimws(am_dt$CL_celltype))
  am_dt$marker <- toupper(trimws(am_dt$marker))
  by_cell <- split(am_dt$marker, am_dt$CL_celltype)
  assign(key, by_cell, envir = MARKER_CACHE)
  by_cell
}

extract_gene_tiers <- function(jin) {
  step1 <- jin$inputs$original_evidence$step1_report_query %||% NULL
  degs <- step1$degs$data %||% list()
  in_scope <- degs$in_scope_top_genes %||% degs$geneSymbol %||% character(0)
  list(
    in_scope_top = normalize_gene_list(in_scope),
    recovery = normalize_gene_list(degs$ranks_16_50_in_scope %||% character(0)),
    out_scope = normalize_gene_list(degs$out_of_scope_top_genes %||% character(0)),
    unknown = normalize_gene_list(degs$unknown_top_genes %||% character(0))
  )
}

normalize_candidate_label <- function(x) {
  x <- normalize_label(x)
  if (is.na(x)) return(NA_character_)
  tolower(stringr::str_replace_all(x, "\\s+", " "))
}

compute_evidence_support <- function(label, marker_map, tiers,
                                    w_top = 2.0, w_recovery = 1.0, w_unknown = 0.5) {
  if (is.null(label) || is.na(label)) return(list(score = 0, top = 0, recovery = 0, unknown = 0, out = 0))
  key <- normalize_candidate_label(label)
  if (is.na(key) || !nzchar(key) || !key %in% names(marker_map)) {
    return(list(score = 0, top = 0, recovery = 0, unknown = 0, out = 0))
  }
  markers <- unique(marker_map[[key]])
  top_hits <- sum(markers %in% tiers$in_scope_top)
  rec_hits <- sum(markers %in% tiers$recovery)
  unk_hits <- sum(markers %in% tiers$unknown)
  out_hits <- sum(markers %in% tiers$out_scope)
  score <- (w_top * top_hits) + (w_recovery * rec_hits) + (w_unknown * unk_hits)
  list(score = score, top = top_hits, recovery = rec_hits, unknown = unk_hits, out = out_hits)
}

# --------------------------
# Compute evidence alignment for each method (cassia/our/enrich) vs raw dossier
# Returns alignment scores (0-1) for each method based on how well their predicted markers
# match the in_scope_top_genes and ranks_16_50_in_scope from raw dossier
# --------------------------
compute_method_evidence_alignment <- function(judge_input_obj, marker_map, tiers, cl_graph, cl_cfg) {
  inputs <- judge_input_obj$inputs %||% list()
  
  # Helper to compute alignment for a single method
  # Generic scoring: favor evidence alignment and specificity, avoid label-specific hacks
  get_alignment <- function(method_summary, method_role) {
    if (is.null(method_summary) || !is.list(method_summary)) return(list(score = 0, confidence = 0))
    
    # Get method's predicted cell type and its markers
    pred_label <- method_summary$top1_cell_type %||% method_summary$predicted_label %||% NA_character_
    if (is.na(pred_label) || !nzchar(pred_label)) return(list(score = 0, confidence = 0))
    
    # Get method's own evidence markers
    method_markers <- method_summary$evidence_markers_method %||% method_summary$evidence_markers %||% character(0)
    if (is.list(method_markers)) method_markers <- unlist(method_markers, recursive = TRUE, use.names = FALSE)
    method_markers <- normalize_gene_list(method_markers)
    
    # Compute overlap with raw dossier tiers (precision-style)
    top_overlap <- sum(method_markers %in% tiers$in_scope_top)
    recovery_overlap <- sum(method_markers %in% tiers$recovery)
    marker_count <- max(1, length(method_markers))
    top_ratio <- top_overlap / marker_count
    recovery_ratio <- recovery_overlap / marker_count

    # Specificity score via ontology depth (generic labels are shallow)
    specificity_score <- 0
    has_cl <- !is.null(cl_graph) && !is.null(cl_cfg)
    if (has_cl) {
      res <- tryCatch(normalize_cl_three_state(pred_label, "", cl_cfg), error = function(e) NULL)
      clid <- res$final_clid %||% NA_character_
      if (!is.na(clid) && nzchar(clid)) {
        depth <- get_depth_to_root(clid, cl_graph)
        specificity_score <- max(0, min(1, depth / 10))
      }
    }

    # Generic penalty for shallow labels (depth <= 2)
    generic_penalty <- if (has_cl && specificity_score <= 0.2) -0.2 else 0

    # Alignment score calculation (generic, dataset-agnostic)
    raw_score <- (0.6 * top_ratio) + (0.3 * recovery_ratio) + (0.1 * specificity_score) + generic_penalty
    normalized_score <- max(0, min(1, raw_score))

    # Confidence based on marker volume only (avoid label-specific heuristics)
    confidence <- min(1, length(method_markers) / 12)

    list(score = normalized_score, confidence = confidence, specificity_score = specificity_score)
  }
  
  cassia_align <- get_alignment(inputs$cassia_summary, "cassia")
  in_house_summary <- inputs$in_house_summary %||% inputs$our_summary
  in_house_align <- get_alignment(in_house_summary, "in_house")
  enrich_summary <- inputs$enrich_summary %||% inputs$inter_summary
  enrich_align <- get_alignment(enrich_summary, "enrichment")
  
  list(
    cassia = cassia_align,
    in_house = in_house_align,
    our = in_house_align,  # legacy alias
    enrich = enrich_align,
    max_score = max(cassia_align$score, in_house_align$score, enrich_align$score, na.rm = TRUE)
  )
}

# --------------------------
# Constants + rules (single source of truth)
# --------------------------
RULES <- list(
  BANNED_PRIMARY_TOKENS = c("mixed", "unresolved", "artifact", "ambiguous", "unknown", "other", "doublet", "multiplet", "triplet"),
  ALLOWED_DECISION_CATEGORIES = c("cassia_better", "in_house_better", "enrich_better", "tie", "third_party_override"),
  TOP_LEVEL_KEYS = c("cluster_id", "final_decision", "method_verdict", "evidence", "post_issues", "cluster_state", "audit_report", "third_party_adjudication", "manual_review_plan"),
  MANUAL_REVIEW = list(
    require_actions = TRUE,
    require_evidence_pointers = TRUE,
    allow_empty_goals = TRUE,
    allow_null_when_no_review = FALSE
  ),
  CITATION = list(
    require_any_if_allowlist_nonempty = TRUE,
    auto_fill_when_missing = FALSE
  )
)

rules_banned_tokens <- function() RULES$BANNED_PRIMARY_TOKENS
rules_allowed_decision_categories <- function() RULES$ALLOWED_DECISION_CATEGORIES
normalize_decision_category <- function(x) {
  value <- tolower(trimws(as.character(x %||% "")))
  if (identical(value, "our_better")) value <- "in_house_better"
  value
}
rules_top_level_keys <- function() RULES$TOP_LEVEL_KEYS

rules_banned_tokens_pattern <- function() {
  paste0("(", paste(rules_banned_tokens(), collapse = "|"), ")")
}

rules_decision_categories_text <- function() {
  paste(rules_allowed_decision_categories(), collapse = " | ")
}

rules_banned_tokens_text <- function() {
  paste(rules_banned_tokens(), collapse = ", ")
}

rules_manual_review_text <- function() {
  paste(
    "Manual review plan rules (authoritative):",
    "- Manual review is a LAST-RESORT release decision, not a default response to reviewer disagreement.",
    "- Set needs_manual_review=true only for unresolved release blockers: no valid final CL ID, unresolved strong mixed/doublet/contamination, unresolved cross-lineage or subtype conflict, unsupported specificity, failed deterministic hard gate, or low confidence (<0.5).",
    "- Do NOT set needs_manual_review merely because an auxiliary candidate score is unavailable, reviewers disagree but a supported dominant identity is selected, a second refinement round was not beneficial, or minor secondary signals remain.",
    "- If post_issues.needs_manual_review=true: manual_review_plan MUST exist and actions must contain >=1 item.",
    "- Each action MUST include evidence_pointers with at least one {type,value}.",
    "- goals may be empty (no FAIL if goals=[]).",
    "- If needs_manual_review=false: manual_review_plan must be priority=null, goals=[], actions=[].",
    sep = "\n"
  )
}

rules_citation_text <- function() {
  paste(
    "Citation rules (authoritative):",
    "- Validation source of truth is validate_citations() in this script.",
    "- Only evidence.cited_pmids and evidence.cited_enrichment_terms count as citations.",
    "- Citations must come from allowlists only; no inferred/new PMID/term.",
    "- When require_any=true and allowlists are non-empty, at least one allowed citation must remain after filtering.",
    sep = "\n"
  )
}

rules_banned_scope_text <- function() {
  "Banned token check scope: ONLY report_to_be_validated.final_decision.primary_cell_type."
}

rules_text_for_head <- function() {
  paste(
    rules_citation_text(),
    rules_manual_review_text(),
    paste0("Allowed decision_category values: ", rules_decision_categories_text()),
    paste0("Banned tokens (primary_cell_type only): ", rules_banned_tokens_text()),
    rules_banned_scope_text(),
    sep = "\n"
  )
}

rules_text_for_chief <- function() {
  paste(
    rules_citation_text(),
    rules_manual_review_text(),
    paste0("Allowed decision_category values: ", rules_decision_categories_text()),
    paste0("Banned tokens (primary_cell_type only): ", rules_banned_tokens_text()),
    rules_banned_scope_text(),
    "Chief QC must enforce the explicit hard constraints above plus the evidence-sufficiency checks defined in its pass_fail_policy; it must not choose a replacement label.",
    "A valid Head identity is the biological anchor being validated. Deterministic post-processing may calibrate specificity only along a direct ontology ancestor/descendant relation and only to a real supported reviewer candidate; it must not switch to a non-hierarchical identity branch.",
    sep = "\n"
  )
}

# The LLM endpoint is resolved at call time in invoke_deepseek_api()
# (LLM_API_BASE_URL env; no hard-coded provider endpoint in the package).

save_json_pretty <- function(obj, path) {
  ensure_dir(dirname(path))
  jsonlite::write_json(obj, path, auto_unbox=TRUE, pretty=TRUE, null="null")
}
read_json_safely <- function(path) jsonlite::fromJSON(path, simplifyVector=FALSE)

# --------------------------
# JSON extraction
# --------------------------
strip_fences <- function(text) {
  if (is.null(text) || text == "") return(NA_character_)
  stringr::str_replace_all(text, "^```json\\s*|\\s*```$", "")
}

extract_first_json_object_stack <- function(text, expected_cluster_id = NULL) {
  if (is.null(text) || text == "") return(NA_character_)
  text2 <- strip_fences(text)
  
  # Collect all JSON object candidates, then pick the best one.
  candidates <- character(0)
  
  ok_full <- tryCatch({ jsonlite::fromJSON(text2, simplifyVector=FALSE); TRUE }, error=function(e) FALSE)
  if (ok_full) candidates <- c(candidates, text2)
  
  chars <- strsplit(text2, "", fixed=TRUE)[[1]]
  start <- NA_integer_
  depth <- 0L
  in_str <- FALSE
  esc <- FALSE
  
  for (i in seq_along(chars)) {
    ch <- chars[[i]]
    if (in_str) {
      if (esc) esc <- FALSE
      else if (ch == "\\") esc <- TRUE
      else if (ch == "\"") in_str <- FALSE
      next
    } else {
      if (ch == "\"") { in_str <- TRUE; next }
    }
    
    if (!in_str) {
      if (ch == "{") {
        if (depth == 0L && is.na(start)) start <- i
        depth <- depth + 1L
      } else if (ch == "}") {
        if (depth > 0L) depth <- depth - 1L
        if (depth == 0L && !is.na(start)) {
          cand <- paste0(chars[start:i], collapse="")
          ok3 <- tryCatch({ jsonlite::fromJSON(cand, simplifyVector=FALSE); TRUE }, error=function(e) FALSE)
          if (ok3) candidates <- c(candidates, cand)
          start <- NA_integer_
        }
      }
    }
  }
  
  if (length(candidates) == 0) return(NA_character_)
  
  score_candidate <- function(obj) {
    if (is.null(obj) || !is.list(obj)) return(-999)
    score <- 0
    
    cid <- as.character(obj$cluster_id %||% NA_character_)
    if (!is.na(cid) && nzchar(cid)) {
      score <- score + 1
      if (stringr::str_detect(cid, "^[0-9]+$")) score <- score + 4 else score <- score - 1
    }
    if (!is.null(expected_cluster_id) && nzchar(expected_cluster_id)) {
      if (!is.na(cid) && identical(cid, expected_cluster_id) && cid != "CLUSTER_0") score <- score + 6
    }
    
    fd <- obj$final_decision %||% NULL
    if (is.list(fd)) {
      score <- score + 3
      dc <- normalize_decision_category(fd$decision_category %||% NA_character_)
      if (!is.na(dc) && dc %in% rules_allowed_decision_categories()) score <- score + 5
      if (!is.na(dc) && stringr::str_detect(dc, "\\|")) score <- score - 5
      
      fct <- as.character(fd$primary_cell_type %||% NA_character_)
      if (!is.na(fct) && nzchar(fct)) {
        score <- score + 1
        if (stringr::str_detect(stringr::str_to_lower(fct), rules_banned_tokens_pattern())) {
          score <- score - 2
        }
      }
    }
    
    if (is.list(obj$method_verdict %||% NULL)) score <- score + 1
    if (is.list(obj$post_issues %||% NULL)) score <- score + 1
    score
  }
  
  parsed <- lapply(candidates, function(s) tryCatch(jsonlite::fromJSON(s, simplifyVector=FALSE), error=function(e) NULL))
  scores <- vapply(parsed, score_candidate, numeric(1))
  best_i <- which(scores == max(scores, na.rm=TRUE))
  # tie-break: prefer later JSON object in the text
  best_i <- best_i[[length(best_i)]]
  candidates[[best_i]]
}

# --------------------------
# DeepSeek API
# --------------------------
classify_http_error <- function(status, err_msg) {
  if (!is.na(status) && status == 429) return("rate_limit")
  if (!is.na(status) && status >= 500) return("server_error")
  if (!is.na(status) && status >= 400) return("client_error")
  if (!is.null(err_msg) && nzchar(err_msg)) return("network_error")
  "unknown"
}

calc_backoff <- function(base_delay, attempt, err_class) {
  delay <- base_delay * (2^(attempt - 1))
  if (err_class == "rate_limit") delay <- delay * 2
  delay <- min(delay, 60)
  delay + stats::runif(1, 0, 1)
}

usage_to_fields <- function(usage) {
  if (is.null(usage) || !is.list(usage)) {
    return(list(prompt_tokens=NA_integer_, completion_tokens=NA_integer_, total_tokens=NA_integer_))
  }
  list(
    prompt_tokens = suppressWarnings(as.integer(usage$prompt_tokens %||% NA_integer_)),
    completion_tokens = suppressWarnings(as.integer(usage$completion_tokens %||% NA_integer_)),
    total_tokens = suppressWarnings(as.integer(usage$total_tokens %||% NA_integer_))
  )
}

append_jsonl <- function(path, obj) {
  ensure_dir(dirname(path))
  line <- jsonlite::toJSON(obj, auto_unbox=TRUE, null="null")
  cat(line, "\n", file=path, append=TRUE)
}

make_request_id <- function(cid, stage, round) {
  paste0(cid, "_", stage, "_", round, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", paste0(sample(c(letters, 0:9), 6, replace=TRUE), collapse=""))
}

# --------------------------
# Scheme A constraints
# --------------------------
judge_output_schema_text <- function() {
  paste(
    "{",
    '  "cluster_id": "CLUSTER_0",',
    '  "cluster_state": "clean | mixed | doublet_suspected | contaminated_or_ambient | transitional_state",',
    "",
    '  "final_decision": {',
    '    "primary_cell_type": "ONE primary label (free text, no mixed/doublet words)",',
    '    "secondary_signals": ["list secondary lineages/contamination signals"],',
    '    "mixture_explanation": "short text explaining mixed/doublet/contamination if present",',
    '    "final_cell_ontology_id": "CL:.... or null",',
    paste0('    "decision_category": "', rules_decision_categories_text(), '",'),
    '    "confidence_primary": 0.0',
    "  },",
    "",
    '  "method_verdict": {',
    '    "cassia": {',
    '      "predicted_cell_type": "string or null",',
    '      "cell_ontology_id": "CL:.... or null",',
    '      "strengths": ["string"],',
    '      "weaknesses": ["string"]',
    "    },",
    '    "our_method": {',
      '      "predicted_cell_type": "string or null",',
      '      "cell_ontology_id": "CL:.... or null",',
      '      "strengths": ["string"],',
      '      "weaknesses": ["string"]',
    "    },",
    '    "enrich": {',
      '      "predicted_cell_type": "string or null",',
      '      "cell_ontology_id": "CL:.... or null",',
      '      "strengths": ["string"],',
      '      "weaknesses": ["string"]',
    "    }",
    "  },",
    "",
    '  "audit_report": {',
    '    "reviewer_support": {',
    '      "cassia_supported": true,',
    '      "our_supported": false,',
    '      "enrich_supported": false',
    "    },",
    '    "flags": ["reviewer_disagreement | missing_citation | lineage_conflict | evidence_weak | other"],',
    '    "notes": "Short, audit-style notes explaining why each reviewer is or is not supported by ORIGINAL evidence."',
    "  },",
    "",
    '  "third_party_adjudication": {',
    '    "primary_cell_type": "string (required if ALL THREE reviewers are explicitly unsupported; otherwise can be null)",',
    '    "confidence": 0.0,',
    '    "evidence_pointers": [',
    '      {"type": "gene | pmid | term", "value": "GENE_X | 1234567 | TERM_X"}',
    "    ],",
    '    "why_reviewers_failed": "Short explanation of why CASSIA/in-house/enrichment are unsupported (if applicable)."',
    "  },",
    "",
    '  "evidence": {',
    '    "supporting": ["string (must mention concrete genes/pmids/terms when available)"],',
    '    "conflicting": ["string"],',
    '    "enrichment_or_literature": ["string"],',
    '    "cited_pmids": ["1234567"],',
    '    "cited_enrichment_terms": ["TERM_X"]',
    "  },",
    "",
    '  "post_issues": {',
    '    "needs_manual_review": false,',
    '    "flags": ["mixed_signature"],',
    '    "notes": "string"',
    "  },",
    "",
    '  "manual_review_plan": {',
    '    "priority": "low | medium | high",',
    '    "goals": ["what uncertainty to resolve"],',
    '    "actions": [',
    "      {",
    '        "action_type": "subcluster | qc_recheck | deg_rerun | doublet_check | ambient_check | reference_mapping | ontology_fix | visual_check | other",',
    '        "why": "Trigger rationale grounded in ORIGINAL evidence (cite genes/pmids/terms).",',
    '        "how": "Executable steps (e.g., rerun DE with covariate; re-cluster at higher resolution; run doublet detection).",',
    '        "expected_outcomes": "If outcome A then ..., if outcome B then ...",',
    '        "evidence_pointers": [',
    '          {"type": "gene | pmid | term", "value": "GENE_X | 1234567 | TERM_X"}',
    "        ]",
    "      }",
    "    ]",
    "  }",
    "}",
    sep="\n"
  )
}

build_head_editor_system_prompt_legacy <- function(species_value = "human") {
  schema_txt <- judge_output_schema_text()
  species_line <- sprintf("SPECIES: %s", species_value)
  paste(
    species_line,
    "ROLE: You are the Head Editor (Handling Editor) for a scientific single-cell annotation review.",
    "",
    "You receive THREE reviewer reports (CASSIA, in-house, clusterProfiler-based enrichment) plus the ORIGINAL BIOLOGICAL EVIDENCE dossier.",
    "Each reviewer report provides a top1_cell_type and may provide a small topk_cell_types shortlist; treat topk as plausible candidates (hints), not ground truth.",
    "You must also act as an auditor: explicitly assess whether each reviewer is supported by ORIGINAL cluster evidence (DEGs/markers, bioinformatic programs/enrichment). Literature may interpret observed evidence but must not substitute for missing cluster-level support.",
    "Adopt the mindset of a careful human expert: first judge cluster state (clean/mixed/doublet_suspected/contaminated/transitional), then choose the most plausible primary cell type while explicitly noting secondary signals/contamination.
IMPORTANT (Avoid overusing mixed/ambiguous): In real data, most clusters are not perfectly pure. Do NOT label a cluster as mixed or transitional_state just because there are a few secondary markers or mild reviewer disagreement.
Use cluster_state='clean' in the default case, and place minor contamination/activation/state signals into final_decision.secondary_signals.
Only use cluster_state='mixed' when there is STRONG, symmetric support for >=2 incompatible candidate labels AND neither dominates after integrating evidence layers. Mild disagreement or an unavailable auxiliary score is not a reason for manual review. Set needs_manual_review=true only when the conflict remains unresolved after evidence integration.",
    "CASSIA and in-house are reviewer opinions, not ground truth. Treat them equally; NEVER default to CASSIA when evidence is equivocal.",
    "Evidence provenance: identity support must originate from observations in this cluster/dossier (e.g., markers/DEGs, assay-derived programs, enrichment/pathway signals). Literature/RAG may interpret the biological meaning or discriminative value of an OBSERVED feature, but cannot supply an unobserved feature or independently create subtype/branch specificity.",
    "Objective evidence: inputs.evidence_summary.lineage may provide a soft low-dimensional lineage hint; do NOT treat it as decisive.",
    "If evidence is mixed, keep one PRIMARY cell type but list secondary_signals and a mixture_explanation; do NOT flip the primary solely because secondary markers are strong.",
    "Lineage priors are weak: do not treat any lineage hint as decisive; rely on integrated biological evidence.",
    "Dominant-identity rule: pick ONE primary label (no mixed/doublet words); secondary signals go to secondary_signals. Mark needs_manual_review only when strong mixed/doublet/contamination evidence remains unresolved and no defensible dominant identity can be released automatically.",
    "Evidence discipline: reviewer labels, candidate availability, ontology names, dataset tissue, and literature about a cell type are hypotheses/priors or interpretation aids, not cluster observations. Shared/general observations support only the identity level they establish; literature cannot create missing specificity.",
    "Hierarchy-aware granularity: distinguish hierarchical refinement from branch switching. If candidate B is a DESCENDANT of candidate A, A and B are not competing hypotheses: do NOT require evidence that excludes A. Require direct positive evidence for the biological component(s) added by B; when those added components are supported, prefer the most specific ontology-valid identity whose claimed components are all supported.",
    "Component preservation: if a proposed label contains several interpretable components and only some are unsupported, remove only the unsupported component(s). Do not collapse a directly supported state/developmental/subtype component merely because an orthogonal component is unresolved. Ontology direct_parents may be used only to navigate to an intermediate valid label, never as biological evidence.",
    "Ontology-defined refinement rule: a state, developmental, maturation, memory, activation, or other modifier MUST NOT automatically be relegated to secondary_signals. If that modifier, combined with the supported lineage, corresponds to an ONTOLOGY-VALID DESCENDANT identity AND the ORIGINAL dossier directly supports the modifier, incorporate it into the PRIMARY identity. Use secondary_signals ONLY for supported programs that are orthogonal to, transient relative to, or not represented in the selected ontology identity. Decision table: (a) state evidence present but ontology has no reasonable descendant -> secondary; (b) ontology has a descendant but no direct evidence for the added component -> do not refine; (c) ontology-valid descendant AND direct evidence for the added component -> primary must be refined.",
    "Non-hierarchical contrastive test: when alternatives are NON_HIERARCHICAL, choose a branch only if the ORIGINAL dossier contains positive evidence that meaningfully favors that branch. If no such evidence exists, do not switch branches by tissue, literature, reviewer majority, or elimination alone.",
    "Branch-stability rule: an exact reviewer majority on a CLID anchors the ontology BRANCH, not the final granularity. Hierarchical refinement away from that branch is allowed with positive evidence for the added component. But overriding that majority to a NON-HIERARCHICAL branch requires explicit, stable conflict resolution (clear branch-specific evidence in the ORIGINAL dossier, or unresolved conflict routed to manual review). Do NOT release a non-hierarchical branch flip driven by a single unstable Head proposal; preserve the majority branch or route to Chief QC / manual review.",
    "Negative-evidence rule: lack of support for a competitor, tissue-based exclusion, or absence of another lineage's markers is not positive support for the selected identity. Missing expected evidence counts against a candidate only when that feature should reasonably be observable in the supplied data; otherwise treat it as unknown.",
    "Frozen ontology context is authoritative for same/ancestor/descendant/non_hierarchical relations; never infer ontology relations from lexical similarity.",
    "Consistency constraint: reasoning text must match selected fields; do NOT state a hypothesis is rejected while selecting it as primary or secondary.",
    "Identity authority: final_decision.primary_cell_type + final_cell_ontology_id are your scientific adjudication. Deterministic downstream rules validate, audit, and route to QC; they are not intended to silently substitute a different valid biological identity.",
    "",
    "You MUST output ONLY one JSON object that EXACTLY matches this schema (same keys; no extra keys):",
    schema_txt,
    "",
    "Rules:",
    "0) Always include audit_report, third_party_adjudication, and manual_review_plan keys. If not applicable, keep fields null/empty per schema.",
    "1) Always set cluster_state first (clean/mixed/doublet_suspected/contaminated_or_ambient/transitional_state).",
    "2) Output ONE primary_cell_type (no mixed/doublet words). Secondary/contamination signals go to secondary_signals with mixture_explanation.",
    "3) decision_category must be one of the allowed categories (see Rules text below).",
    "3b) Use ECC-style review criteria when deciding decision_category (apply the SAME rubric to all three reviewers):",
    "    - Evidence coverage: did they integrate multiple OBSERVED evidence layers (markers/DEGs plus enrichment/programs when available), using literature only to interpret those observations?",
    "    - Consistency: is the conclusion internally consistent with the dossier (and do they acknowledge mixed/doublet signals when present)?",
    "    - Conservatism appropriateness: if evidence is insufficient for a very specific subtype, do they avoid over-specific labels while still choosing ONE primary type?",
    "4) Strengths/weaknesses: cite why each reviewer is stronger/weaker (observed marker/program/enrichment support, literature interpretation of observed features, handling of mixed signals).",
    "4b) audit_report.reviewer_support must explicitly say whether CASSIA, in-house, and enrichment are each supported by ORIGINAL evidence.",
    "4c) For a hierarchical refinement, evidence.supporting must state the positive observation supporting the ADDED component(s); it need not exclude the ancestor. For a non-hierarchical branch choice, evidence.supporting must state the positive observation that favors the selected branch.",
    "5) Enrichment/Literature must be referenced when present; marker lists alone are insufficient if higher-level evidence exists.",
    "6) If ALL THREE reviewers are explicitly unsupported (cassia_supported=false, our_supported=false, enrich_supported=false), you MUST provide third_party_adjudication with a primary_cell_type and evidence_pointers, and set decision_category=third_party_override.",
    "6b) If mixed/doublet/contamination is strongly supported and unresolved: keep a primary, set needs_manual_review=true, and explain the mixture. For minor or resolved secondary signals, keep needs_manual_review=false and record an audit note; do NOT output mixed/doublet words in primary_cell_type.",
    "7) No default to Cassia; treat lineage hints as weak prior only.",
    "",
    "Rules (single source of truth):",
    rules_text_for_head(),
    "",
    "Strict JSON rules:",
    "A) Use JSON null (without quotes). NEVER output \"null\" as a string.",
    "B) primary_cell_type MUST NOT contain banned tokens listed above.",
    "",
    "Hard constraint: If primary_cell_type contains a banned token, the output is INVALID.",
    "Do not output any text outside the JSON.",
    sep="\n"
  )
}


# --------------------------
# Compact Head prompt (A/B profile)
# Keeps the same schema/output contract but consolidates overlapping biology rules.
# --------------------------
build_head_editor_system_prompt_compact <- function(species_value = "human") {
  schema_txt <- judge_output_schema_text()
  paste(
    sprintf("SPECIES: %s", species_value),
    "ROLE: You are the Handling Editor for a scientific single-cell annotation review.",
    "INPUTS: three reviewer reports plus the ORIGINAL biological dossier and frozen ontology_context. Reviewer labels/top-k are hypotheses, not ground truth.",

    "DECISION ORDER (authoritative):",
    "1) Establish the dominant biological lineage/identity from ORIGINAL cluster observations. Minor contamination/activation does not by itself make a cluster mixed.",
    "2) HIERARCHICAL REFINEMENT: ancestor and descendant are not competing hypotheses. Choose the most specific ontology-valid descendant for which every ADDED biological component has direct positive support. Preserve supported components; remove only unsupported orthogonal components.",
    "3) NON-HIERARCHICAL BRANCH CHOICE: require ORIGINAL cluster-level positive evidence that meaningfully favors the selected branch. Shared/general programs do not distinguish branches.",
    "4) EVIDENCE PROVENANCE: markers/DEGs, assay-derived programs, and enrichment/pathway observations may support identity. Tissue context, reviewer votes, candidate availability, literature, or elimination of alternatives cannot create missing specificity. Literature may only interpret an actually observed feature.",
    "5) PRIMARY VS SECONDARY: an ontology-valid supported modifier may belong in the primary identity; use secondary_signals only for orthogonal/transient/contaminating programs not represented in the selected identity.",

    "ONTOLOGY / STABILITY:",
    "- ontology_context is authoritative for same/ancestor/descendant/non_hierarchical relations; never infer ontology relations from wording.",
    "- An exact reviewer CLID majority is a branch anchor, not a granularity lock. Hierarchical refinement is allowed with positive support for added components. A non-hierarchical majority override requires clear branch-specific ORIGINAL evidence; otherwise keep the majority branch or route unresolved conflict to QC/manual review.",
    "- direct_parents are navigation aids only, never biological evidence.",

    "REVIEWER AUDIT:",
    "- Apply the same evidence rubric to CASSIA, in-house, and enrichment. Explicitly mark each reviewer supported/unsupported from ORIGINAL observations.",
    "- If all three are unsupported, use third_party_adjudication; otherwise do not invent an unnecessary fourth label.",

    "CLUSTER STATE / RELEASE:",
    "- Default cluster_state='clean'. Use mixed/doublet/contaminated/transitional only for strong unresolved evidence. Keep one primary identity and record minor secondary signals separately.",
    "- needs_manual_review is last-resort for unresolved release blockers, not ordinary reviewer disagreement.",

    "OUTPUT REQUIREMENTS:",
    paste0("- decision_category must be one of: ", rules_decision_categories_text(), "."),
    paste0("- primary_cell_type must not contain: ", rules_banned_tokens_text(), "."),
    "- Cite only allowed PMIDs/enrichment terms supplied in citation_requirements; one valid citation is sufficient when citation is required.",
    "- If needs_manual_review=false, manual_review_plan must be priority=null, goals=[], actions=[].",
    "- evidence.supporting must name the positive observation supporting any hierarchical ADDED component or non-hierarchical branch choice.",
    "Return ONLY one JSON object matching this schema exactly:",
    schema_txt,
    "Use JSON null, never the string \"null\". Do not add extra top-level keys or prose outside JSON.",
    sep = "\n"
  )
}

build_head_editor_system_prompt <- function(species_value = "human") {
  profile <- tolower(trimws(as.character(PROMPT_PROFILE %||% "legacy")))
  if (identical(profile, "compact")) {
    build_head_editor_system_prompt_compact(species_value)
  } else {
    build_head_editor_system_prompt_legacy(species_value)
  }
}


# --------------------------
# Citation requirements (deterministic)
# --------------------------
sanitize_pmids <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x) && !is.character(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- as.character(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) return(character(0))
  # Extract 6-10 digit sequences; robust to embedded newlines/spaces like "2849400\n5"
  digits <- stringr::str_extract_all(x, "[0-9]{6,10}")
  digits <- unlist(digits, use.names = FALSE)
  digits <- digits[nzchar(digits)]
  unique(digits)
}

extract_allowed_pmids <- function(jin) {
  pmids <- NULL
  pmids <- jin$inputs$original_evidence$step1_report_query$evidence$retrieved_articles$articles_db$pmid %||%
    jin$inputs$original_evidence$step1_report_query$evidence$retrieved_articles$pmid %||%
    jin$inputs$original_evidence$evidence$retrieved_articles$articles_db$pmid %||%
    NULL
  sanitize_pmids(pmids)
}

extract_allowed_enrichment_terms <- function(jin, limit = 50L) {
  bio <- jin$inputs$original_evidence$step1_report_query$bioinfo$analysis_results %||% NULL
  if (is.null(bio) || !is.list(bio)) return(character(0))
  
  terms <- character(0)
  for (k in names(bio)) {
    v <- bio[[k]]
    if (is.list(v) && length(v) > 0) {
      # list of objects: pull common name fields
      if (is.list(v[[1]])) {
        for (item in v) {
          if (!is.list(item)) next
          cand <- item$term %||% item$description %||% item$name %||% item$pathway %||% item$id %||% NULL
          if (!is.null(cand)) terms <- c(terms, as.character(cand))
        }
      }
    }
  }
  terms <- unique(str_trim(terms))
  terms <- terms[nzchar(terms)]
  if (length(terms) > limit) terms <- terms[seq_len(limit)]
  terms
}

extract_citation_requirements <- function(jin) {
  pmids <- extract_allowed_pmids(jin)
  enrich <- extract_allowed_enrichment_terms(jin)
  require_any_cfg <- isTRUE(RULES$CITATION$require_any_if_allowlist_nonempty)
  require_any <- if (isTRUE(require_any_cfg)) (length(pmids) > 0 || length(enrich) > 0) else FALSE
  list(
    require_any = require_any,
    allowed_pmids = pmids,
    allowed_enrichment_terms = enrich,
    preferred_pmids = pmids,
    preferred_enrichment_terms = enrich,
    auto_fill_when_missing = isTRUE(RULES$CITATION$auto_fill_when_missing)
  )
}

validate_citations <- function(j, citation_requirements) {
  reasons <- character(0)
  if (is.null(citation_requirements) || !is.list(citation_requirements)) {
    citation_requirements <- list(require_any = FALSE, allowed_pmids = character(0), allowed_enrichment_terms = character(0))
  }
  allowed_pmids <- sanitize_pmids(citation_requirements$allowed_pmids %||% character(0))
  allowed_terms <- unique(str_trim(as.character(citation_requirements$allowed_enrichment_terms %||% character(0))))
  allowed_terms <- allowed_terms[nzchar(allowed_terms)]
  require_any <- isTRUE(citation_requirements$require_any %||% FALSE)

  cited_pmids <- character(0)
  cited_terms <- character(0)
  if (is.list(j$evidence)) {
    cited_pmids <- sanitize_pmids(j$evidence$cited_pmids %||% character(0))
    cited_terms <- j$evidence$cited_enrichment_terms %||% character(0)
    if (is.list(cited_terms) && !is.character(cited_terms)) cited_terms <- unlist(cited_terms, recursive = TRUE, use.names = FALSE)
    cited_terms <- unique(str_trim(as.character(cited_terms)))
    cited_terms <- cited_terms[nzchar(cited_terms)]
  }

  if (length(allowed_pmids) > 0 && length(cited_pmids) > 0) {
    bad <- setdiff(cited_pmids, allowed_pmids)
    if (length(bad) > 0) reasons <- c(reasons, paste0("cited_pmids_not_in_allowlist: ", paste(bad, collapse=",")))
  }
  if (length(allowed_terms) > 0 && length(cited_terms) > 0) {
    badt <- setdiff(cited_terms, allowed_terms)
    if (length(badt) > 0) reasons <- c(reasons, paste0("cited_enrichment_terms_not_in_allowlist: ", paste(badt, collapse=" | ")))
  }

  if (isTRUE(require_any) && (length(allowed_pmids) > 0 || length(allowed_terms) > 0)) {
    ok_pmid <- length(cited_pmids) > 0 && length(allowed_pmids) > 0
    ok_term <- length(cited_terms) > 0 && length(allowed_terms) > 0
    if (!(ok_pmid || ok_term)) reasons <- c(reasons, "missing_required_citation")
  }

  list(ok = length(reasons) == 0, reasons = reasons, cited_pmids = cited_pmids, cited_terms = cited_terms)
}

# --------------------------
# Core label helpers
# --------------------------
normalize_label <- function(x) {
  if (is.null(x)) return(NA_character_)
  x2 <- as.character(x)
  if (length(x2) == 0) return(NA_character_)
  x2 <- x2[[1]]
  x2 <- stringr::str_replace_all(x2, "\\s+", " ")
  x2 <- stringr::str_trim(x2)
  if (!nzchar(x2)) return(NA_character_)
  x2
}

sanitize_for_validation <- function(x) {
  if (is.null(x)) return(x)
  x2 <- normalize_label(x)
  if (is.na(x2)) return(x2)
  x2 <- stringr::str_replace(x2, "\\s*\\(.*$", "")
  x2 <- stringr::str_replace_all(x2, "[,:;]+$", "")
  x2 <- stringr::str_trim(x2)
  if (!nzchar(x2)) NA_character_ else x2
}

sanitize_for_display <- function(x) {
  normalize_label(x)
}

ensure_post_issues <- function(j) {
  if (is.null(j$post_issues) || !is.list(j$post_issues)) {
    j$post_issues <- list(needs_manual_review = FALSE, flags = list(), notes = "")
  }
  if (is.null(j$post_issues$flags)) j$post_issues$flags <- list()
  if (is.null(j$post_issues$notes)) j$post_issues$notes <- ""
  j
}


# --------------------------
# Guardrail: avoid overusing mixed/ambiguous labels
# --------------------------
downgrade_overused_cluster_state <- function(j) {
  if (is.null(j) || !is.list(j)) return(j)
  st <- as.character(j$cluster_state %||% "")
  conf <- suppressWarnings(as.numeric(j$final_decision$confidence_primary %||% NA_real_))
  decision_cat <- normalize_decision_category(j$final_decision$decision_category %||% "")
  original_state <- st
  # If the model says mixed/transitional but confidence is reasonably high OR there is a clear decision winner,
  # treat it as a mostly-clean cluster with secondary signals, and keep uncertainty in needs_manual_review/post_issues.
  if (st %in% c("mixed","transitional_state")) {
    clear_winner <- decision_cat %in% c("cassia_better","in_house_better","enrich_better","third_party_override")
    if (is.finite(conf) && conf >= 0.7 && clear_winner) {
      j$cluster_state <- "clean"
      j$post_issues <- ensure_post_issues(j)$post_issues
      note0 <- str_trim(as.character(j$post_issues$notes %||% ""))
      j$post_issues$notes <- str_trim(paste(note0, "Mixed/transitional state softened to clean; uncertainty kept in secondary_signals/manual review."))
    }
  }
  # Avoid 'ambiguous' as a default: if primary label exists and confidence is not extremely low, keep clean.
  if (st == "ambiguous") {
    if ((nzchar(j$final_decision$primary_cell_type %||% "")) && (is.finite(conf) && conf >= 0.55)) {
      j$cluster_state <- "clean"
      j$post_issues <- ensure_post_issues(j)$post_issues
      note0 <- str_trim(as.character(j$post_issues$notes %||% ""))
      j$post_issues$notes <- str_trim(paste(note0, "Ambiguous state softened to clean; usable primary label with moderate confidence."))
    }
  }
  if (!is.na(original_state) && nzchar(original_state) && original_state != (j$cluster_state %||% "")) {
    j$post_issues <- ensure_post_issues(j)$post_issues
    j$post_issues$cluster_state_raw <- original_state
    flags <- j$post_issues$flags %||% list()
    if (!is.list(flags)) flags <- as.list(flags)
    j$post_issues$flags <- unique(c(unlist(flags, recursive = TRUE, use.names = FALSE), "mixed_or_ambiguous_signals"))
  }
  j
}

ensure_enrich_verdict <- function(j, judge_input_obj) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.null(j$method_verdict) || !is.list(j$method_verdict)) j$method_verdict <- list()
  existing <- j$method_verdict$enrich %||% j$method_verdict$inter %||% NULL
  enrich_in <- judge_input_obj$inputs$enrich_summary %||% judge_input_obj$inputs$inter_summary %||% list()
  enrich_obj <- existing %||% list(
    predicted_cell_type = enrich_in$predicted_label %||% NULL,
    cell_ontology_id = enrich_in$cell_ontology_id %||% enrich_in$cl_id %||% NULL,
    strengths = list(),
    weaknesses = list()
  )
  j$method_verdict$enrich <- enrich_obj
  # `inter` mirrors the canonical `enrich` field; it is not an additional reviewer.
  j$method_verdict$inter <- enrich_obj
  j
}

ensure_in_house_alias <- function(j) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.null(j$method_verdict) || !is.list(j$method_verdict)) return(j)
  # Keep backward compatibility: preserve our_method, add in_house mirror.
  if (!is.null(j$method_verdict$our_method) && is.list(j$method_verdict$our_method)) {
    j$method_verdict$in_house <- j$method_verdict$our_method
  } else if (!is.null(j$method_verdict$in_house) && is.list(j$method_verdict$in_house)) {
    j$method_verdict$our_method <- j$method_verdict$in_house
  }
  j
}

as_json_string_array <- function(x) {
  if (is.null(x)) return(list())
  if (is.list(x)) {
    vals <- unlist(x, recursive = TRUE, use.names = FALSE)
    vals <- as.character(vals)
    vals <- vals[nzchar(vals)]
    return(as.list(vals))
  }
  vals <- as.character(x)
  vals <- vals[nzchar(vals)]
  as.list(vals)
}

force_array_fields <- function(j) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.list(j$post_issues)) {
    j$post_issues$flags <- as_json_string_array(j$post_issues$flags %||% list())
    j$post_issues$audit_flags <- as_json_string_array(j$post_issues$audit_flags %||% list())
    j$post_issues$auto_qc_flags <- as_json_string_array(j$post_issues$auto_qc_flags %||% list())
    j$post_issues$release_blockers <- as_json_string_array(j$post_issues$release_blockers %||% list())
  }
  if (is.list(j$audit_report)) {
    j$audit_report$flags <- as_json_string_array(j$audit_report$flags %||% list())
  }
  j
}

ensure_audit_report_support <- function(j) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.null(j$audit_report) || !is.list(j$audit_report)) {
    j$audit_report <- list(reviewer_support = list(), flags = list(), notes = "")
  }
  rs <- j$audit_report$reviewer_support %||% list()
  if (!is.list(rs)) rs <- list()

  # Missing reviewer-support assessments are UNKNOWN, never silently promoted to TRUE.
  # The deterministic local gate requires explicit TRUE/FALSE for all three canonical reviewers.
  if (is.null(rs$cassia_supported)) rs$cassia_supported <- NA
  if (is.null(rs$our_supported)) rs$our_supported <- NA
  if (is.null(rs$enrich_supported) && !is.null(rs$inter_supported)) rs$enrich_supported <- rs$inter_supported
  if (is.null(rs$enrich_supported)) rs$enrich_supported <- NA

  # `inter` mirrors enrichment only; it is not a fourth reviewer.
  if (is.null(rs$inter_supported)) rs$inter_supported <- rs$enrich_supported
  j$audit_report$reviewer_support <- rs
  if (is.null(j$audit_report$flags)) j$audit_report$flags <- list()
  if (is.null(j$audit_report$notes)) j$audit_report$notes <- ""
  j
}

normalize_manual_review_plan <- function(j) {
  j <- ensure_post_issues(j)
  needs <- isTRUE(j$post_issues$needs_manual_review %||% FALSE)
  if (!needs) {
    j$manual_review_plan <- list(priority = NULL, goals = list(), actions = list())
    return(j)
  }
  if (is.null(j$manual_review_plan) || !is.list(j$manual_review_plan)) {
    j$manual_review_plan <- list(priority = "medium", goals = character(0), actions = list())
  }
  j
}

extract_first_gene_from_jin <- function(jin) {
  if (is.null(jin) || !is.list(jin)) return(NA_character_)
  step1 <- jin$inputs$original_evidence$step1_report_query %||% NULL
  g <- step1$degs$data$geneSymbol %||% NULL
  if (!is.null(g)) {
    g <- as.character(g)
    g <- g[!is.na(g) & nzchar(g)]
    if (length(g) > 0) return(g[[1]])
  }
  te <- step1$candidates_for_evaluation$data$Top_Evidence_Genes %||% NULL
  if (!is.null(te)) {
    if (is.list(te)) te <- unlist(te, recursive = TRUE, use.names = FALSE)
    te <- as.character(te)
    te <- te[!is.na(te) & nzchar(te)]
    if (length(te) > 0) return(te[[1]])
  }
  NA_character_
}

filter_citations_to_allowlist <- function(j, citation_requirements) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.null(citation_requirements) || !is.list(citation_requirements)) return(j)
  allowed_pmids <- sanitize_pmids(citation_requirements$allowed_pmids %||% character(0))
  allowed_terms <- unique(str_trim(as.character(citation_requirements$allowed_enrichment_terms %||% character(0))))
  allowed_terms <- allowed_terms[nzchar(allowed_terms)]
  require_any <- isTRUE(citation_requirements$require_any)
  if (!is.list(j$evidence)) return(j)
  pm <- j$evidence$cited_pmids %||% list()
  if (is.list(pm)) pm <- unlist(pm, recursive = TRUE, use.names = FALSE)
  pm <- sanitize_pmids(pm)
  if (length(allowed_pmids) > 0) {
    pm <- pm[pm %in% allowed_pmids]
  }
  te <- j$evidence$cited_enrichment_terms %||% list()
  if (is.list(te)) te <- unlist(te, recursive = TRUE, use.names = FALSE)
  te <- unique(str_trim(as.character(te)))
  te <- te[!is.na(te) & nzchar(te)]
  if (length(allowed_terms) > 0) {
    te <- te[te %in% allowed_terms]
  }

  # Keep this function as pure filtering only. Missing-citation policy is handled
  # by validate_citations()+retry/manual-review flow, not by mutating outputs here.

  j$evidence$cited_pmids <- as.list(pm)
  j$evidence$cited_enrichment_terms <- as.list(te)
  j
}

ensure_manual_review_action <- function(j, jin, citation_requirements) {
  j <- ensure_post_issues(j)
  needs <- isTRUE(j$post_issues$needs_manual_review %||% FALSE)
  if (!needs) return(j)
  mrp <- j$manual_review_plan %||% list()
  actions <- mrp$actions %||% list()
  if (is.list(actions) && length(actions) > 0) return(j)

  infer_transition_context <- function(x) {
    fd <- x$final_decision %||% list()
    st <- as.character(x$cluster_state %||% "")
    raw_st <- as.character(x$post_issues$cluster_state_raw %||% "")
    mix_txt <- tolower(as.character(fd$mixture_explanation %||% ""))
    sec_raw <- fd$secondary_signals %||% list()
    if (!is.list(sec_raw)) sec_raw <- as.list(sec_raw)
    sec <- as.character(unlist(sec_raw, recursive = TRUE, use.names = FALSE))
    sec <- sec[!is.na(sec) & nzchar(sec)]
    flags_raw <- x$post_issues$flags %||% list()
    flags <- as.character(unlist(flags_raw, recursive = TRUE, use.names = FALSE))
    flags <- flags[nzchar(flags)]

    has_transition_state <- st %in% c("transitional_state") || raw_st %in% c("transitional_state")
    has_transition_text <- grepl("transit|trajectory|differentiat|toward|towards", mix_txt)
    is_transition_like <- isTRUE(has_transition_state) || isTRUE(has_transition_text)

    list(
      transition_like = is_transition_like,
      primary = as.character(fd$primary_cell_type %||% NA_character_),
      secondary = sec
    )
  }

  transition_ctx <- infer_transition_context(j)
  gene <- extract_first_gene_from_jin(jin)
  if (is.na(gene)) {
    allowed_terms <- citation_requirements$allowed_enrichment_terms %||% character(0)
    if (length(allowed_terms) > 0) {
      ep <- list(list(type = "term", value = as.character(allowed_terms[[1]])))
    } else {
      ep <- list()
    }
  } else {
    ep <- list(list(type = "gene", value = gene))
  }

  if (isTRUE(transition_ctx$transition_like) && length(transition_ctx$secondary) > 0) {
    to_label <- transition_ctx$secondary[[1]]
    from_label <- transition_ctx$primary %||% "current primary"
    mrp$actions <- list(list(
      action_type = "other",
      why = paste0("Potential transitional trajectory detected: ", from_label, " -> ", to_label, "."),
      how = "Run pseudotime/subcluster validation and compare lineage marker continuity across candidate states.",
      expected_outcomes = paste0("If transition is supported, keep primary anchor and document trajectory toward ", to_label, "; if not, retain primary and treat secondary as contamination/state effect."),
      evidence_pointers = ep
    ))
  } else {
    mrp$actions <- list(list(
      action_type = "subcluster",
      why = "Mixed or ambiguous signals require separating potential subpopulations.",
      how = "Re-cluster at higher resolution and re-run DEG analysis on resulting subclusters.",
      expected_outcomes = "If distinct lineages appear, annotate separately; otherwise retain primary label with mixed flag.",
      evidence_pointers = ep
    ))
  }
  j$manual_review_plan <- mrp
  j
}

standardize_transition_notes <- function(j) {
  if (is.null(j) || !is.list(j) || !is.list(j$final_decision)) return(j)
  fd <- j$final_decision
  primary <- as.character(fd$primary_cell_type %||% "")
  if (!nzchar(primary)) return(j)

  sec_raw <- fd$secondary_signals %||% list()
  if (!is.list(sec_raw)) sec_raw <- as.list(sec_raw)
  sec <- as.character(unlist(sec_raw, recursive = TRUE, use.names = FALSE))
  sec <- sec[!is.na(sec) & nzchar(sec)]

  st <- as.character(j$cluster_state %||% "")
  raw_st <- as.character(j$post_issues$cluster_state_raw %||% "")
  mix_txt <- as.character(fd$mixture_explanation %||% "")
  mix_txt_low <- tolower(mix_txt)
  is_transition_state <- st %in% c("transitional_state") || raw_st %in% c("transitional_state")
  has_transition_text <- grepl("transit|trajectory|differentiat|toward|towards", mix_txt_low)
  if (!(is_transition_state || has_transition_text)) return(j)

  if (length(sec) == 0) return(j)

  to_label <- sec[[1]]
  canonical <- paste0("Transitional signal detected: primary kept as ", primary, " for evaluation; possible trajectory toward ", to_label, ".")
  if (!nzchar(stringr::str_squish(mix_txt))) {
    j$final_decision$mixture_explanation <- canonical
  } else if (!grepl("primary kept as .* for evaluation", mix_txt_low)) {
    j$final_decision$mixture_explanation <- stringr::str_squish(paste(mix_txt, canonical))
  }

  j <- ensure_post_issues(j)
  note0 <- str_trim(as.character(j$post_issues$notes %||% ""))
  if (!grepl("possible trajectory toward", tolower(note0), fixed = TRUE)) {
    j$post_issues$notes <- str_trim(paste(note0, paste0("Transition note: primary=", primary, "; possible_to=", to_label, ".")))
  }
  j
}

flag_issue <- function(j, flag, note=NULL) {
  j <- ensure_post_issues(j)
  flags_raw <- j$post_issues$flags %||% list()
  flags_vec <- as.character(unlist(flags_raw, recursive = TRUE, use.names = FALSE))
  flags_vec <- flags_vec[nzchar(flags_vec)]
  flags_vec <- unique(c(flags_vec, as.character(flag)))
  j$post_issues$flags <- as.list(flags_vec)
  if (!is.null(note) && nzchar(note)) {
    note0 <- str_trim(as.character(j$post_issues$notes %||% ""))
    j$post_issues$notes <- str_trim(paste(note0, note))
  }
  j
}

labels_match_loose <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  if (length(a) == 0) a <- ""
  if (length(b) == 0) b <- ""
  a <- a[[1]]; b <- b[[1]]  
  if (is.na(a) || is.na(b)) return(FALSE)
  
  a <- tolower(stringr::str_trim(a))
  b <- tolower(stringr::str_trim(b))
  if (!nzchar(a) || !nzchar(b)) return(FALSE)
  
  if (identical(a, b)) return(TRUE)
  if (stringr::str_detect(a, stringr::fixed(b, ignore_case = TRUE))) return(TRUE)
  if (stringr::str_detect(b, stringr::fixed(a, ignore_case = TRUE))) return(TRUE)
  FALSE
}

# --------------------------
# 5-class reviewer support relative to the final adjudicated CL ID.
#   exact_or_equivalent      final == reviewer (same CL term)
#   compatible_parent        reviewer is an ancestor of final (compatible but broader)
#   compatible_descendant    reviewer is a descendant of final (compatible but narrower)
#   incompatible             both valid CL IDs but no hierarchy relation
#   unresolved               reviewer CL ID missing/invalid OR final CL ID missing
# --------------------------
cl_support_class <- function(final_clid, reviewer_clid, cl_graph) {
  clean <- function(x) {
    x <- as.character(x)
    if (length(x) == 0L) return(NA_character_)
    x <- x[[1]]
    x <- trimws(x)
    if (is.na(x) || !nzchar(x) || toupper(x) == "NA") return(NA_character_)
    x
  }
  f <- clean(final_clid)
  r <- clean(reviewer_clid)
  if (is.na(f) || is.na(r)) return("unresolved")
  if (identical(f, r)) return("exact_or_equivalent")
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return("unresolved")
  node_f <- cl_graph$cl[[f]]
  node_r <- cl_graph$cl[[r]]
  if (is.null(node_f) || is.null(node_r)) return("unresolved")
  anc_f <- node_f$ancestors %||% NULL
  anc_r <- node_r$ancestors %||% NULL
  anc_f_ids <- unique(c(names(anc_f), if (is.atomic(anc_f)) as.character(unname(anc_f)) else character(0)))
  anc_r_ids <- unique(c(names(anc_r), if (is.atomic(anc_r)) as.character(unname(anc_r)) else character(0)))
  if (r %in% anc_f_ids) return("compatible_parent")
  if (f %in% anc_r_ids) return("compatible_descendant")
  "incompatible"
}

# --------------------------
# Deterministic: populate method correctness relative to final primary label.
# (adds nested fields is_correct; allowed but optional)
# --------------------------
set_method_correctness_vs_final <- function(j, judge_input_obj, cl_graph = NULL) {
  if (is.null(j) || !is.list(j)) return(j)
  if (is.null(j$final_decision) || !is.list(j$final_decision)) return(j)
  primary <- sanitize_for_display(j$final_decision$primary_cell_type %||% NA_character_)
  if (is.na(primary) || !nzchar(primary)) return(j)
  final_clid <- sanitize_for_display(j$final_decision$final_cell_ontology_id %||% NA_character_)

  cassia_label <- sanitize_for_display(judge_input_obj$inputs$cassia_summary$top1_cell_type %||% NA_character_)
  in_house_summary_in <- judge_input_obj$inputs$in_house_summary %||% judge_input_obj$inputs$our_summary %||% list()
  our_label <- sanitize_for_display(in_house_summary_in$top1_cell_type %||% NA_character_)
  enrich_summary_in <- judge_input_obj$inputs$enrich_summary %||% judge_input_obj$inputs$inter_summary %||% list()
  enrich_label <- sanitize_for_display(enrich_summary_in$predicted_label %||% NA_character_)

  # CL IDs: prefer frozen registry mapping, fall back to reviewer-provided CL ID.
  cassia_clid <- judge_input_obj$inputs$cassia_summary$mapping$cell_ontology_id %||%
    judge_input_obj$inputs$cassia_summary$cell_ontology_id %||% NA_character_
  our_clid <- in_house_summary_in$mapping$cell_ontology_id %||%
    in_house_summary_in$cell_ontology_id %||% NA_character_
  enrich_clid <- enrich_summary_in$mapping$cell_ontology_id %||%
    enrich_summary_in$cell_ontology_id %||% NA_character_

  if (is.null(j$method_verdict) || !is.list(j$method_verdict)) j$method_verdict <- list()
  if (is.null(j$method_verdict$cassia) || !is.list(j$method_verdict$cassia)) j$method_verdict$cassia <- list()
  if (is.null(j$method_verdict$our_method) || !is.list(j$method_verdict$our_method)) j$method_verdict$our_method <- list()
  if (is.null(j$method_verdict$enrich) || !is.list(j$method_verdict$enrich)) j$method_verdict$enrich <- list()

  set_support <- function(mv_key, label, clid) {
    if (is.null(j$method_verdict[[mv_key]])) j$method_verdict[[mv_key]] <- list()
    j$method_verdict[[mv_key]]$matches_final_label <<- isTRUE(labels_match_loose(primary, label))
    j$method_verdict[[mv_key]]$support_class <<- cl_support_class(final_clid, clid, cl_graph)
    j$method_verdict[[mv_key]]$reviewer_cl_id <<- clean_clid(clid)
    j$method_verdict[[mv_key]]$final_cl_id <<- clean_clid(final_clid)
    if (is.null(j$method_verdict[[mv_key]]$is_correct)) {
      # is_correct: exact or compatible parent/descendant (hierarchy-aware)
      j$method_verdict[[mv_key]]$is_correct <<- j$method_verdict[[mv_key]]$support_class %in%
        c("exact_or_equivalent", "compatible_parent", "compatible_descendant")
    }
    j$method_verdict[[mv_key]]
  }
  clean_clid <- function(x) {
    x <- as.character(x)
    if (length(x) == 0L) return(NA_character_)
    x <- trimws(x[[1]])
    if (is.na(x) || !nzchar(x) || toupper(x) == "NA") return(NA_character_)
    x
  }

  set_support("cassia", cassia_label, cassia_clid)
  set_support("our_method", our_label, our_clid)
  set_support("enrich", enrich_label, enrich_clid)
  # `inter` mirrors `enrich` only.
  j$method_verdict$inter <- j$method_verdict$enrich
  j
}

# --------------------------
# Head Editor query
# --------------------------
build_head_editor_instruction <- function() {
  "Return ONLY the adjudication JSON."
}


# --------------------------
# Generic frozen-ontology context
# --------------------------
# IMPORTANT:
# - No biological label, tissue, marker, dataset, or GT is hard-coded here.
# - Relations are computed only from the frozen CL graph + reviewer candidates.
# - This context informs adjudication/QC but never supplies a ground-truth answer.

canonical_reviewer_method_ontology <- function(x) {
  x <- tolower(trimws(as.character(x %||% "")))
  if (x %in% c("our", "our_method", "in_house", "in_house_method")) return("our")
  if (x %in% c("inter", "enrich", "enrichment", "enrichment_method")) return("enrich")
  if (x %in% c("cassia")) return("cassia")
  x
}

ontology_relation_record <- function(a, b, cl_graph) {
  a <- as.character(a %||% NA_character_)
  b <- as.character(b %||% NA_character_)
  valid_a <- !is.na(a) && nzchar(a) && !is.null(cl_graph$cl[[a]])
  valid_b <- !is.na(b) && nzchar(b) && !is.null(cl_graph$cl[[b]])
  if (!isTRUE(valid_a) || !isTRUE(valid_b)) {
    return(list(
      relation = "unavailable",
      lca_clid = NA_character_,
      lca_label = NA_character_,
      a_distance_to_lca = NA_real_,
      b_distance_to_lca = NA_real_
    ))
  }

  relation <- if (identical(a, b)) {
    "same"
  } else if (isTRUE(is_ancestor_of(a, b, cl_graph))) {
    "a_is_ancestor_of_b"
  } else if (isTRUE(is_ancestor_of(b, a, cl_graph))) {
    "b_is_ancestor_of_a"
  } else {
    "non_hierarchical"
  }

  anc_a <- get_anc_map(a, cl_graph) %||% numeric(0)
  anc_b <- get_anc_map(b, cl_graph) %||% numeric(0)
  common <- intersect(names(anc_a), names(anc_b))
  lca <- NA_character_
  if (length(common) > 0L) {
    depths <- vapply(common, function(id) get_depth_to_root(id, cl_graph), numeric(1))
    lca <- common[[which.max(depths)]]
  }

  local_label <- function(id) {
    if (is.na(id) || !nzchar(id) || is.null(cl_graph$cl[[id]])) return(NA_character_)
    as.character(cl_graph$cl[[id]]$label %||% NA_character_)
  }

  list(
    relation = relation,
    lca_clid = lca,
    lca_label = local_label(lca),
    a_distance_to_lca = if (!is.na(lca) && lca %in% names(anc_a)) as.numeric(anc_a[[lca]]) else NA_real_,
    b_distance_to_lca = if (!is.na(lca) && lca %in% names(anc_b)) as.numeric(anc_b[[lca]]) else NA_real_
  )
}

reviewer_support_from_head <- function(head_out, method) {
  rs <- head_out$audit_report$reviewer_support %||% list()
  method <- canonical_reviewer_method_ontology(method)
  keys <- switch(
    method,
    cassia = c("cassia_supported"),
    our = c("our_supported", "in_house_supported", "our_method_supported"),
    enrich = c("enrich_supported", "enrichment_supported", "inter_supported"),
    character(0)
  )
  for (k in keys) {
    v <- rs[[k]] %||% NULL
    if (is.logical(v) && length(v) == 1L && !is.na(v)) return(v)
  }
  NA
}

build_reviewer_ontology_context <- function(judge_input_obj, cl_graph, cl_cfg) {
  if (is.null(cl_graph) || is.null(cl_graph$cl) || is.null(judge_input_obj) || !is.list(judge_input_obj)) {
    return(list(status = "unavailable"))
  }

  cands <- tryCatch(
    collect_candidates_from_inputs(judge_input_obj$inputs %||% list(), cl_cfg, cl_graph),
    error = function(e) list()
  )
  cands <- Filter(function(x) {
    method <- canonical_reviewer_method_ontology(x$method %||% "")
    clid <- as.character(x$clid %||% NA_character_)
    rank <- suppressWarnings(as.integer(x$rank %||% 99L))
    method %in% c("cassia", "our", "enrich") &&
      !is.na(clid) && nzchar(clid) && !is.null(cl_graph$cl[[clid]]) &&
      !is.na(rank) && rank <= 3L
  }, cands)

  # one record per reviewer/rank/CLID
  if (length(cands) > 0L) {
    keys <- vapply(cands, function(x) {
      paste(
        canonical_reviewer_method_ontology(x$method %||% ""),
        as.integer(x$rank %||% 99L),
        as.character(x$clid %||% ""),
        sep = "|"
      )
    }, character(1))
    cands <- cands[!duplicated(keys)]
  }

  local_label <- function(id) {
    if (is.na(id) || !nzchar(id) || is.null(cl_graph$cl[[id]])) return(NA_character_)
    as.character(cl_graph$cl[[id]]$label %||% NA_character_)
  }

  candidate_records <- lapply(cands, function(x) {
    cid <- as.character(x$clid %||% NA_character_)
    anc <- cl_graph$cl[[cid]]$ancestors %||% list()
    av <- suppressWarnings(as.numeric(unlist(anc, use.names = TRUE)))
    pids <- if (length(av) > 0L) names(av[!is.na(av) & av == 1]) else character(0)
    parents <- lapply(unique(pids), function(pid) {
      list(clid = pid, label = local_label(pid))
    })
    list(
      method = canonical_reviewer_method_ontology(x$method %||% ""),
      rank = as.integer(x$rank %||% 99L),
      submitted_label = as.character(x$label %||% NA_character_),
      clid = cid,
      canonical_label = local_label(cid),
      direct_parents = parents,
      mapping_status = as.character(x$coerced_status %||% ""),
      map_quality = as.integer(x$map_quality %||% 0L)
    )
  })

  # Strict majority is descriptive context only, never a forced label.
  top1 <- Filter(function(x) as.integer(x$rank %||% 99L) == 1L, cands)
  if (length(top1) > 0L) {
    methods <- vapply(top1, function(x) canonical_reviewer_method_ontology(x$method %||% ""), character(1))
    top1 <- top1[!duplicated(methods)]
  }
  n_methods <- length(top1)
  majority_clid <- NA_character_
  majority_votes <- 0L
  majority_fraction <- NA_real_
  if (n_methods > 0L) {
    ids <- vapply(top1, function(x) as.character(x$clid %||% NA_character_), character(1))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) > 0L) {
      tb <- sort(table(ids), decreasing = TRUE)
      majority_votes <- as.integer(tb[[1]])
      majority_fraction <- majority_votes / n_methods
      unique_winner <- sum(tb == tb[[1]]) == 1L
      if (isTRUE(unique_winner) && majority_votes >= 2L && majority_votes > n_methods / 2) {
        majority_clid <- names(tb)[[1]]
      }
    }
  }

  unique_ids <- unique(vapply(cands, function(x) as.character(x$clid %||% NA_character_), character(1)))
  unique_ids <- unique_ids[!is.na(unique_ids) & nzchar(unique_ids)]
  pairwise <- list()
  if (length(unique_ids) >= 2L) {
    cmb <- utils::combn(unique_ids, 2, simplify = FALSE)
    pairwise <- lapply(cmb, function(pair) {
      rel <- ontology_relation_record(pair[[1]], pair[[2]], cl_graph)
      list(
        clid_a = pair[[1]],
        label_a = local_label(pair[[1]]),
        clid_b = pair[[2]],
        label_b = local_label(pair[[2]]),
        relation = rel$relation,
        lca_clid = rel$lca_clid,
        lca_label = rel$lca_label,
        a_distance_to_lca = rel$a_distance_to_lca,
        b_distance_to_lca = rel$b_distance_to_lca
      )
    })
  }

  list(
    status = "ok",
    ontology_source = "frozen_CL_graph",
    ontology_rule = paste(
      "Ontology relations in this object are authoritative.",
      "Do not infer equivalence, ancestry, or subtype relations from wording similarity.",
      "Reviewer majority is contextual evidence only and does not determine the final identity.",
      "direct_parents are ontology navigation aids only; they are not reviewer votes or biological evidence."
    ),
    reviewer_candidates = candidate_records,
    direct_top1_majority = list(
      clid = majority_clid,
      label = local_label(majority_clid),
      votes = majority_votes,
      n_reviewers = n_methods,
      fraction = majority_fraction
    ),
    pairwise_relations = pairwise
  )
}

build_chief_ontology_context <- function(head_out, judge_input_obj, cl_graph, cl_cfg) {
  base <- build_reviewer_ontology_context(judge_input_obj, cl_graph, cl_cfg)
  if (!identical(base$status %||% "", "ok")) return(base)

  head_clid <- as.character(head_out$final_decision$final_cell_ontology_id %||% NA_character_)
  head_label <- as.character(head_out$final_decision$primary_cell_type %||% NA_character_)
  if (is.na(head_clid) || !nzchar(head_clid) || is.null(cl_graph$cl[[head_clid]])) {
    base$selected_head <- list(clid = head_clid, label = head_label, valid = FALSE)
    base$requires_discriminative_evidence <- FALSE
    return(base)
  }

  det_clid <- as.character(head_out$decision_trace$deterministic_pre_preserve_clid %||% NA_character_)
  dc_clid <- as.character(head_out$decision_trace$direct_consensus_clid %||% NA_character_)

  # Attach Head relation/support to each reviewer candidate.
  recs <- base$reviewer_candidates %||% list()
  recs <- lapply(recs, function(rec) {
    rel <- ontology_relation_record(head_clid, rec$clid, cl_graph)
    method <- canonical_reviewer_method_ontology(rec$method %||% "")
    rec$relation_to_head <- rel$relation
    rec$lca_clid <- rel$lca_clid
    rec$lca_label <- rel$lca_label
    rec$head_distance_to_lca <- rel$a_distance_to_lca
    rec$candidate_distance_to_lca <- rel$b_distance_to_lca
    rec$reviewer_supported_by_head <- reviewer_support_from_head(head_out, method)
    rec$is_direct_majority_anchor <- !is.na(dc_clid) && nzchar(dc_clid) && identical(rec$clid, dc_clid)
    rec$is_deterministic_proposal <- !is.na(det_clid) && nzchar(det_clid) && identical(rec$clid, det_clid)
    rec
  })
  base$reviewer_candidates <- recs

  local_label <- function(id) {
    if (is.na(id) || !nzchar(id) || is.null(cl_graph$cl[[id]])) return(NA_character_)
    as.character(cl_graph$cl[[id]]$label %||% NA_character_)
  }

  det_rel <- ontology_relation_record(head_clid, det_clid, cl_graph)
  dc_rel <- ontology_relation_record(head_clid, dc_clid, cl_graph)

  priority_nonhier <- Filter(function(rec) {
    isTRUE(identical(rec$relation_to_head %||% "", "non_hierarchical")) &&
      (
        isTRUE(rec$is_deterministic_proposal %||% FALSE) ||
        isTRUE(rec$is_direct_majority_anchor %||% FALSE) ||
        (as.integer(rec$rank %||% 99L) == 1L && isTRUE(rec$reviewer_supported_by_head %||% FALSE))
      )
  }, recs)

  base$selected_head <- list(
    clid = head_clid,
    label = local_label(head_clid),
    submitted_label = head_label,
    valid = TRUE
  )
  base$deterministic_proposal <- list(
    clid = det_clid,
    label = local_label(det_clid),
    relation_to_head = det_rel$relation,
    lca_clid = det_rel$lca_clid,
    lca_label = det_rel$lca_label
  )
  base$direct_majority_anchor <- list(
    clid = dc_clid,
    label = local_label(dc_clid),
    relation_to_head = dc_rel$relation,
    lca_clid = dc_rel$lca_clid,
    lca_label = dc_rel$lca_label
  )
  base$requires_discriminative_evidence <- length(priority_nonhier) > 0L
  base$priority_nonhierarchical_competitors <- priority_nonhier
  base$qc_rule <- paste(
    "When requires_discriminative_evidence=true, positive evidence that is compatible with multiple competing branches is not sufficient by itself.",
    "The selected Head identity must be justified by ORIGINAL dossier evidence that distinguishes its ontology branch or claimed specificity.",
    "If such distinguishing evidence is absent, fail with an existing biological QC tag; do not choose a replacement label."
  )
  base
}

build_head_editor_query <- function(judge_input_obj, cl_graph = NULL, cl_cfg = NULL) {
  list(
    query_id = paste0(judge_input_obj$cluster_id, "_head_editor"),
    analysis_type = "Head Editor: Final adjudication between three reviewers",
    llm_prompt = list(
      system = build_head_editor_system_prompt(
        species_value = if (exists("cfg", inherits = TRUE) && !is.null(cfg$species)) cfg$species[[1]] %||% "human" else "human"),
      species_value = if (exists("cfg", inherits = TRUE) && !is.null(cfg$species)) cfg$species[[1]] %||% "human" else "human",
      instruction = build_head_editor_instruction()
    ),
    citation_requirements = extract_citation_requirements(judge_input_obj),
    ontology_context = build_reviewer_ontology_context(judge_input_obj, cl_graph, cl_cfg),
    inputs = judge_input_obj$inputs %||% list(),
    cluster_id = judge_input_obj$cluster_id
  )
}

# --------------------------
# Chief QC (QC-only; Scheme A compatible)
# --------------------------
build_chief_editor_instructions_legacy <- function() {
  list(
    role = "You are the Chief Editor performing QC on a Head Editor adjudication JSON.",
    task = paste(
      "You will be given:",
      "1) report_to_be_validated (Head adjudication JSON),",
      "2) original_dossier (biological evidence dossier),",
      "3) ontology_context (frozen CL graph relations among the Head selection, deterministic proposal, majority anchor, and reviewer candidates, when available).",
      "",
      "Your job is QC only, NOT to re-adjudicate from scratch and NOT to choose a replacement label.",
      "Validate BOTH report integrity and whether the ORIGINAL biological evidence is sufficient for the selected primary label and its claimed specificity. Apply contrastive review without treating the mere existence of a competing candidate as evidence against the Head selection.",
      "If biological support is insufficient, FAIL with a structured reason; the Head Editor, not the Chief, performs any re-adjudication.",
      "",
      "Scheme A enforcement (single source of truth):",
      rules_text_for_chief(),
      "Scope reminders:",
      "- Only check banned tokens against report_to_be_validated.final_decision.primary_cell_type.",
      "- Only evidence.cited_pmids and evidence.cited_enrichment_terms count as citations.",
      "- Do NOT infer citations from free text.",
      "FAILURE REPORTING FORMAT (required):",
      "- Every item in failure_reasons must include the JSON field path(s) inspected and the observed value(s).",
      "- If you cannot cite a concrete field path + value, do NOT mark FAIL; emit an audit_warning instead.",
      "- Resolved reviewer disagreement, unavailable auxiliary scores, and minor secondary signals are audit items and do not require manual review. Strong unresolved mixed/doublet/artifact evidence must be represented by release_blockers + needs_manual_review=true."
      ,"- manual_review_plan.actions[].evidence_pointers may reference genes/terms and do NOT require PMID matching."
    ),
    pass_fail_policy = list(
      pass_definition = c(
        "PASS if schema is correct, internal consistency holds, and evidence section does not fabricate markers/pathways.",
        "PASS when the selected primary label and its claimed specificity have sufficient positive support in the ORIGINAL dossier. Distinguishing support may be a coherent combination of multiple observed markers/programs/enrichment features; it need not be unique, exclusive, or externally re-validated inside the dossier.",
        "PASS when resolved uncertainty is documented in audit_flags without manual review; require needs_manual_review=true only when release_blockers remain unresolved."
      ),
      strict_output_rules = c(
        "NO SELF-CORRECTION / NO NARRATIVE: Do not include 'however', 'but', 're-evaluate', 'actually', 'on second thought', or any self-correction phrasing anywhere in failure_reasons. If you are unsure, do not FAIL; put it in audit_warnings.",
        "FAIL REASON FORMAT (strict): each failure reason MUST be exactly one line in this template:",
        "TAG | FIELD_PATH | OBSERVED_VALUE | RULE",
        "TAG must be one of: MISSING_KEY, BANNED_TOKEN, BAD_CONFIDENCE_RANGE, BAD_DECISION_CATEGORY, CITATION_MISSING, CITATION_NOT_ALLOWED, MANUAL_REVIEW_PLAN_MISSING, THIRD_PARTY_REQUIRED_MISSING, EVIDENCE_INSUFFICIENT, LABEL_EVIDENCE_MISMATCH, UNSUPPORTED_SPECIFICITY, UNRESOLVED_LINEAGE_CONFLICT, OTHER_SCHEMA.",
        "Each line must correspond to exactly one fail condition. Do not combine multiple issues in one line.",
        "If validation_status is PASSED, failure_reasons MUST be an empty list []. If FAILED, failure_reasons MUST be non-empty.",
        "Do not write contradictory statements across fields; the presence of any failure_reasons implies FAILED.",
        "PASS output example:",
        "validation_status='VALIDATION PASSED'",
        "failure_reasons=[]",
        "audit_warnings=['...optional...']",
        "FAILED output example:",
        "validation_status='VALIDATION FAILED'",
        "failure_reasons=[",
        "BANNED_TOKEN | final_decision.primary_cell_type | 'mixed CELL_LABEL_C' | banned token present",
        "]"
      ),
      fail_conditions = c(
        "FAIL if required keys are missing or confidence not in [0,1].",
        "FAIL if any null is written as string \"null\".",
        "FAIL if primary_cell_type contains banned tokens (see Rules).",
        "FAIL if decision_category is not one of allowed categories (see Rules).",
        "FAIL if citation_requirements allow ANY PMIDs or enrichment terms, AND BOTH arrays are empty: evidence.cited_pmids is empty AND evidence.cited_enrichment_terms is empty.",
        "FAIL if cited_pmids or cited_enrichment_terms contain values outside the allowlist.",
        "MANUAL_REVIEW_PLAN_MISSING: FAIL only if report_to_be_validated.post_issues.needs_manual_review=true AND (manual_review_plan.actions missing OR actions length == 0). Do NOT require goals.",
        "Do NOT fail manual_review evidence_pointers for missing PMID. Gene/term pointers are valid if type/value are non-empty.",
        "THIRD_PARTY_REQUIRED_MISSING: FAIL if all three reviewers are explicitly unsupported AND third_party_adjudication is missing a primary_cell_type or evidence_pointers.",
        "EVIDENCE_INSUFFICIENT: FAIL if the selected primary label has no concrete positive support from cluster-derived observations in the ORIGINAL dossier (markers/DEGs, assay-derived programs, enrichment/pathways). Literature may establish the meaning of an observed feature but cannot substitute for a missing cluster observation.",
        "LABEL_EVIDENCE_MISMATCH: FAIL if the Head's claimed supporting evidence is materially inconsistent with the selected lineage/identity, the ORIGINAL dossier more strongly supports an incompatible lineage, OR the Head explicitly documents direct positive support for a biologically meaningful component but drops that supported component from the primary identity solely because a different orthogonal component is unresolved.",
        "UNSUPPORTED_SPECIFICITY: for a hierarchical refinement, FAIL only if the biological component(s) added beyond the ancestor lack direct positive cluster-level support; do NOT demand evidence that excludes the ancestor. For a non-hierarchical branch choice, require positive evidence that meaningfully favors the selected branch.",
        "ONTOLOGY-CONTEXT RULE: ontology_context is authoritative for ontology relations; lexical similarity is not. direct_parents are navigation aids only. requires_discriminative_evidence=true applies to non-hierarchical contrastive review, not ancestor-versus-descendant refinement.",
        "FALSIFICATION RULE: reviewer choice, candidate availability, dataset tissue, compatibility, exclusion of alternatives, or absence of another lineage's markers are not distinguishing positive evidence. Missing expected evidence is counterevidence only when it should reasonably be observable in the supplied data.",
        "Evidence-sufficiency rule: observed features need not be exclusive to one identity. A coherent set of multiple positive cluster observations that jointly favors the selected identity is sufficient; literature may explain why those OBSERVED features are discriminative, but literature/context alone cannot create missing specificity. Do not demand a single uniquely specific marker. FAIL only if no positive distinguishing pattern remains. Do not prescribe a replacement label.",
        "UNRESOLVED_LINEAGE_CONFLICT: FAIL if strong evidence supports incompatible lineages with no dominant identity, yet the Head releases a single primary label without an unresolved blocker/manual-review state.",
        "BRANCH-OVERRIDE: if the Head selects a NON-HIERARCHICAL branch that overrides an exact reviewer CLID majority, FAIL unless the ORIGINAL dossier contains clear branch-specific positive evidence justifying the override; absence of such evidence is UNSUPPORTED_SPECIFICITY/LABEL_EVIDENCE_MISMATCH. A single Head proposal without stable support must not flip a majority-anchored branch.",
        "For biological FAIL tags, cite the relevant Head field path and summarize the concrete dossier evidence in OBSERVED_VALUE/RULE. Do not invent markers or external facts.",
        "When a positive observation is actually present in the dossier, interpret its biological relevance; do not dismiss it merely because the dossier does not explicitly label it as a validated distinguishing marker.",
        "Do NOT FAIL merely because reviewers disagree, one auxiliary score is unavailable, or minor activation/contamination/secondary signals exist when a dominant identity is adequately supported.",
        "The Chief must never prescribe a replacement label. suggested_fixes may request evidence repair, broader specificity, or re-adjudication, but must not name a ground-truth answer.",
        "One citation is sufficient. Do NOT require both PMID and term, and do NOT require more than one.",
        "Citation compliance is separate from biological evidence sufficiency: an allowed PMID/term can document interpretation, but its presence alone does not provide cluster-level support for the selected identity.",
        "PASS example: dossier has PMIDs=10 terms=20; Head cited_pmids=[\"PMID_X\"], cited_terms=[] => PASS.",
        "PASS example: dossier has PMIDs=0 terms>0; Head cited_terms=[\"TERM_X\"] => PASS.",
        "FAIL example: dossier has PMIDs/terms; Head cited_pmids=[] and cited_terms=[] => FAIL."
      ),
      non_fail_audit_items = c(
        "Low confidence is NOT a fail by itself -> AUDIT WARNING.",
        "Missing ontology id -> AUDIT WARNING.",
        "Overstated wording -> AUDIT WARNING.",
        "Reviewer disagreement alone is NOT a fail if the Head resolves it with sufficient evidence.",
        "Unavailable auxiliary reviewer score/confidence is NOT a fail by itself.",
        "Minor activation, contamination, or secondary lineage signals are NOT a fail when the dominant primary identity remains sufficiently supported.",
        "Generic labels like 'GENERIC_LABEL_A', 'GENERIC_LABEL_B', 'GENERIC_LABEL_C' are NOT banned tokens. Do NOT FAIL; use audit_warnings only if you think the label is too generic.",
        "Only FAIL for explicit fail_conditions listed above. If suboptimal but not listed, record audit_warnings only."
      )
    ),
    output_format = list(
      schema = list(
        validation_status = "Either 'VALIDATION PASSED' or 'VALIDATION FAILED'",
        failure_reasons = "list[string] (Only if FAILED)",
        schema_ok = "boolean",
        consistency_ok = "boolean",
        evidence_ok = "boolean",
        audit_warnings = "list[string]",
        suggested_fixes = "list[string]",
        final_verdict_summary = "string"
      )
    )
  )
}

build_chief_system_prompt <- function() {
  ins <- build_chief_editor_instructions()
  paste(as.character(ins$role %||% ""), as.character(ins$task %||% ""), sep = "\n")
}

build_chief_editor_query <- function(cluster_id, judge_output_obj, original_dossier_obj, citation_requirements,
                                     ontology_context = NULL) {
  list(
    query_id = paste0(cluster_id, "_chief_editor_qc"),
    analysis_type = "Chief Editor: QC for Head Editor adjudication JSON",
    instructions_for_llm = build_chief_editor_instructions(),
    input_data = list(
      report_to_be_validated = judge_output_obj,
      original_dossier = original_dossier_obj,
      ontology_context = ontology_context %||% list(status = "unavailable"),
      citation_requirements = citation_requirements,
      qc_helpers = list(
        needs_manual_review = isTRUE(judge_output_obj$post_issues$needs_manual_review %||% FALSE),
        release_state = as.character(judge_output_obj$post_issues$release_state %||% ""),
        audit_flags = as.list(get_audit_flags(judge_output_obj)),
        release_blockers = as.list(get_release_blockers(judge_output_obj)),
        manual_review_actions_n = length(judge_output_obj$manual_review_plan$actions %||% list()),
        cited_pmids = sanitize_pmids(judge_output_obj$evidence$cited_pmids %||% list()),
        cited_enrichment_terms = unique(str_trim(as.character(unlist(judge_output_obj$evidence$cited_enrichment_terms %||% list(), recursive=TRUE))))
      )
    )
  )
}


# Compact Chief prompt: QC only; no re-adjudication and no replacement label.
build_chief_editor_instructions_compact <- function() {
  list(
    role = "You are the Chief Editor performing QC on a Head Editor adjudication JSON.",
    task = paste(
      "QC ONLY. Do not re-adjudicate from scratch and do not choose a replacement label.",
      "Use report_to_be_validated, ORIGINAL dossier, and ontology_context.",
      "Check four biological questions:",
      "1) Does the selected dominant lineage/identity have concrete positive ORIGINAL cluster support?",
      "2) If the selection is a hierarchical refinement, are the ADDED biological components directly supported? Do not require the descendant to exclude its ancestor.",
      "3) If the selection is non-hierarchical to a plausible competitor, is there positive ORIGINAL evidence that meaningfully favors the selected branch?",
      "4) Did the Head improperly use tissue context, reviewer votes, candidate availability, literature, or exclusion of alternatives to create missing specificity? Literature may interpret observed evidence but cannot replace it.",
      "Shared/general programs support only the level they establish. Do not demand one unique marker; a coherent set of positive observations can be sufficient.",
      "If a Head explicitly documents a supported component but drops it only because an orthogonal component is unresolved, treat that as LABEL_EVIDENCE_MISMATCH.",
      "For an exact reviewer-majority branch override, a non-hierarchical flip requires clear branch-specific ORIGINAL evidence; otherwise it is UNSUPPORTED_SPECIFICITY/LABEL_EVIDENCE_MISMATCH.",
      "Schema/citation/manual-review failures must cite a concrete field path and observed value. If you cannot point to a concrete failure, use audit_warnings rather than FAIL.",
      "Chief must never name a replacement label.",
      sep = "\n"
    ),
    pass_fail_policy = list(
      pass_definition = c(
        "PASS when schema/internal consistency are valid and the selected identity/specificity is supported under the four biological checks.",
        "Resolved reviewer disagreement, unavailable auxiliary scores, and minor secondary signals are audit items, not failures."
      ),
      strict_output_rules = c(
        "If PASSED, failure_reasons must be []. If FAILED, failure_reasons must be non-empty.",
        "Each failure reason must be: TAG | FIELD_PATH | OBSERVED_VALUE | RULE.",
        "Allowed biological tags: EVIDENCE_INSUFFICIENT, LABEL_EVIDENCE_MISMATCH, UNSUPPORTED_SPECIFICITY, UNRESOLVED_LINEAGE_CONFLICT. Existing schema/citation tags remain allowed.",
        "Do not prescribe a replacement identity."
      ),
      fail_conditions = c(
        "EVIDENCE_INSUFFICIENT: selected identity lacks concrete positive ORIGINAL cluster support.",
        "UNSUPPORTED_SPECIFICITY: added hierarchical component lacks direct support, or a non-hierarchical branch lacks positive branch-favoring evidence.",
        "LABEL_EVIDENCE_MISMATCH: Head reasoning/evidence is materially inconsistent with the selected identity, including dropping a directly supported component because another orthogonal component is unresolved.",
        "UNRESOLVED_LINEAGE_CONFLICT: incompatible lineages remain strongly supported with no dominant identity but the report releases a single label without a blocker.",
        "Tissue/context/literature/reviewer vote/exclusion alone never establishes missing specificity.",
        "Apply existing deterministic schema, citation, manual-review-plan, and third-party requirements supplied in rules_text_for_chief()."
      ),
      non_fail_audit_items = c(
        "Reviewer disagreement alone is not a failure.",
        "Low confidence, unavailable auxiliary scores, or minor activation/contamination alone are not failures.",
        "If evidence is suboptimal but no explicit fail condition is met, use audit_warnings."
      )
    ),
    output_format = list(
      schema = list(
        validation_status = "Either 'VALIDATION PASSED' or 'VALIDATION FAILED'",
        failure_reasons = "list[string] (Only if FAILED)",
        schema_ok = "boolean",
        consistency_ok = "boolean",
        evidence_ok = "boolean",
        audit_warnings = "list[string]",
        suggested_fixes = "list[string]",
        final_verdict_summary = "string"
      )
    ),
    deterministic_rules = rules_text_for_chief()
  )
}

build_chief_editor_instructions <- function() {
  profile <- tolower(trimws(as.character(PROMPT_PROFILE %||% "legacy")))
  if (identical(profile, "compact")) {
    build_chief_editor_instructions_compact()
  } else {
    build_chief_editor_instructions_legacy()
  }
}

# --------------------------
# Load the latest head editor output for a cluster (best-effort).
# --------------------------
load_latest_head_output <- function(cid, out_dir_head) {
  if (is.null(out_dir_head) || !dir.exists(out_dir_head)) return(NULL)
  pat <- paste0("^", cid, "_LLM_JUDGE_OUTPUT_round[0-9]+\\.json$")
  files <- list.files(out_dir_head, pattern = pat, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rounds <- suppressWarnings(as.integer(stringr::str_extract(basename(files), "[0-9]+$")))
  files <- files[order(rounds, decreasing = TRUE)]
  read_json_safely(files[[1]])
}

# If final outputs lost required fields, hydrate from the last head output on disk.
hydrate_final_from_head <- function(final_out, cid, head_out_last = NULL, out_dir_head = NULL) {
  if (is.null(final_out) || !is.list(final_out)) final_out <- list()
  missing_keys <- c("final_decision", "method_verdict", "evidence")
  needs_fill <- vapply(missing_keys, function(k) is.null(final_out[[k]]), logical(1))
  if (!any(needs_fill)) return(final_out)
  
  src <- NULL
  if (is.list(head_out_last) && !is.null(head_out_last$final_decision)) src <- head_out_last
  if (is.null(src)) src <- load_latest_head_output(cid, out_dir_head)
  if (!is.list(src)) return(final_out)
  
  filled <- FALSE
  for (k in missing_keys) {
    if (is.null(final_out[[k]]) && !is.null(src[[k]])) {
      final_out[[k]] <- src[[k]]
      filled <- TRUE
    }
  }
  
  if (filled) {
    final_out <- ensure_post_issues(final_out)
    flags <- final_out$post_issues$flags %||% character(0)
    final_out$post_issues$flags <- unique(c(as.character(flags), "final_repaired_from_head"))
    note0 <- str_trim(as.character(final_out$post_issues$notes %||% ""))
    final_out$post_issues$notes <- str_trim(paste(note0, "Hydrated missing fields from head output."))
  }
  
  final_out
}

# --------------------------
# Ontology helpers (no hardcoded labels)
# --------------------------

LABEL_CLID_CACHE <- new.env(parent = emptyenv())

normalize_label_for_match <- function(x) {
  x %||% "" |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("\\[.*?\\]", " ") |>
    stringr::str_replace_all("\\(.*?\\)", " ") |>
    stringr::str_replace_all("[^a-z0-9]+", " ") |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_trim()
}

label_tokens <- function(x) {
  n <- normalize_label_for_match(x)
  if (!nzchar(n)) return(character(0))
  toks <- unlist(strsplit(n, " ", fixed = TRUE), use.names = FALSE)
  toks <- toks[nzchar(toks)]
  unique(toks)
}

jaccard_sim <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  inter / uni
}

best_token_match <- function(label, cl_graph) {
  idx <- build_label_index(cl_graph)
  if (is.null(idx)) return(list(best_id = NA_character_, best_score = 0, second_score = 0, ties = 0))
  toks <- label_tokens(label)
  if (length(toks) == 0) return(list(best_id = NA_character_, best_score = 0, second_score = 0, ties = 0))

  score_map <- list()
  for (term in idx$term_index %||% list()) {
    if (is.null(term$id) || is.null(term$tokens)) next
    sc <- jaccard_sim(toks, term$tokens)
    if (sc <= 0) next
    prev <- score_map[[term$id]] %||% -Inf
    if (sc > prev) score_map[[term$id]] <- sc
  }
  if (length(score_map) == 0) return(list(best_id = NA_character_, best_score = 0, second_score = 0, ties = 0))

  ids <- names(score_map)
  scores <- vapply(ids, function(x) score_map[[x]], numeric(1))
  best_score <- max(scores)
  best_ids <- ids[scores == best_score]
  second_score <- if (length(scores[scores < best_score]) > 0) max(scores[scores < best_score]) else 0
  list(best_id = best_ids[[1]], best_score = best_score, second_score = second_score, ties = length(best_ids))
}

build_label_index <- function(cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(NULL)
  key <- "label_index"
  if (exists(key, envir = LABEL_CLID_CACHE, inherits = FALSE)) {
    return(get(key, envir = LABEL_CLID_CACHE, inherits = FALSE))
  }
  idx_label <- list()
  idx_syn <- list()
  term_index <- list()
  for (id in names(cl_graph$cl)) {
    if (!startsWith(id, "CL:")) next
    term <- cl_graph$cl[[id]]
    if (isTRUE(term$deprecated)) next
    lbl <- term$label %||% NULL
    syn <- term$synonyms %||% NULL
    add_lbl <- function(x, target, add_term_index = FALSE) {
      if (is.null(x)) return()
      x <- as.character(x)[1]
      if (is.na(x) || !nzchar(x)) return()
      n <- normalize_label_for_match(x)
      if (!nzchar(n)) return()
      if (target == "label") {
        if (is.null(idx_label[[n]])) idx_label[[n]] <<- character(0)
        idx_label[[n]] <<- unique(c(idx_label[[n]], id))
      } else {
        if (is.null(idx_syn[[n]])) idx_syn[[n]] <<- character(0)
        idx_syn[[n]] <<- unique(c(idx_syn[[n]], id))
      }
      if (isTRUE(add_term_index)) {
        term_index[[length(term_index) + 1]] <<- list(id = id, label = x, tokens = label_tokens(x))
      }
    }
    add_lbl(lbl, "label", add_term_index = TRUE)
    if (!is.null(syn)) {
      syns <- if (is.list(syn)) unlist(syn, recursive = TRUE, use.names = FALSE) else syn
      for (s in syns) add_lbl(s, "syn", add_term_index = FALSE)
    }
  }
  assign(key, list(map_label = idx_label, map_syn = idx_syn, term_index = term_index), envir = LABEL_CLID_CACHE)
  get(key, envir = LABEL_CLID_CACHE, inherits = FALSE)
}

get_label_candidates_exact <- function(label, cl_graph) {
  if (is.null(label)) return(character(0))
  lbl <- as.character(label)[1]
  if (is.na(lbl) || !nzchar(lbl)) return(character(0))
  idx <- build_label_index(cl_graph)
  if (is.null(idx)) return(character(0))
  key <- normalize_label_for_match(lbl)
  if (!nzchar(key)) return(character(0))
  cands <- idx$map_label[[key]] %||% character(0)
  unique(as.character(cands))
}

get_synonym_candidates_exact <- function(label, cl_graph) {
  if (is.null(label)) return(character(0))
  lbl <- as.character(label)[1]
  if (is.na(lbl) || !nzchar(lbl)) return(character(0))
  idx <- build_label_index(cl_graph)
  if (is.null(idx)) return(character(0))
  key <- normalize_label_for_match(lbl)
  if (!nzchar(key)) return(character(0))
  cands <- idx$map_syn[[key]] %||% character(0)
  unique(as.character(cands))
}

coerce_clid <- function(label, clid, cl_cfg, cl_graph, min_best = 0.6, min_delta = 0.1) {
  strip_parenthetical <- function(x) {
    if (is.null(x)) return(NA_character_)
    s <- as.character(x)[1]
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    s <- gsub("\\\\([^)]*\\\\)", "", s)
    s <- stringr::str_squish(s)
    if (!nzchar(s)) return(NA_character_)
    s
  }

  label_clean <- strip_parenthetical(label)
  raw <- as.character(clid %||% "")
  raw_present <- nzchar(raw)

  exact_label_cands <- get_label_candidates_exact(label_clean, cl_graph)
  exact_syn_cands <- character(0)
  if (length(exact_label_cands) == 0) {
    exact_syn_cands <- get_synonym_candidates_exact(label_clean, cl_graph)
  }
  exact_cands <- if (length(exact_label_cands) > 0) exact_label_cands else exact_syn_cands
  exact_unique <- (length(exact_cands) == 1)
  exact_ambig <- (length(exact_cands) > 1)

  token_match <- best_token_match(label_clean, cl_graph)
  token_ok <- (token_match$best_score >= min_best) && ((token_match$best_score - token_match$second_score) >= min_delta) && (token_match$ties == 1)

  candidate_id <- NA_character_
  status <- "unmapped"
  candidates <- exact_cands

  if (exact_unique) {
    candidate_id <- exact_cands[[1]]
    if (length(exact_label_cands) > 0) {
      status <- "mapped_exact"
    } else if (token_ok && token_match$best_id == candidate_id) {
      status <- "mapped_synonym"
    } else if (token_ok && token_match$best_id != candidate_id) {
      candidate_id <- token_match$best_id
      status <- "mapped_token"
    } else {
      candidate_id <- NA_character_
      status <- "ambiguous"
    }
  } else if (!exact_ambig && token_ok) {
    candidate_id <- token_match$best_id
    status <- "mapped_token"
  } else if (exact_ambig) {
    status <- "ambiguous"
  }

  # Priority: when label maps exactly to a CLID, use that consistently
  # This ensures same label always gets same CLID regardless of raw input
  if (exact_unique && !is.na(candidate_id) && nzchar(candidate_id)) {
    return(list(clid = candidate_id, status = status, candidates = candidates, coerced_clid = candidate_id))
  }
  
  if (raw_present) {
    if (!is.na(candidate_id) && nzchar(candidate_id)) {
      if (isTRUE(candidate_id == raw) || isTRUE(is_ancestor_of(raw, candidate_id, cl_graph))) {
        return(list(clid = candidate_id, status = "raw_used", candidates = candidates, coerced_clid = candidate_id))
      }
      return(list(clid = raw, status = "raw_incompatible", candidates = candidates, coerced_clid = candidate_id))
    }
    return(list(clid = raw, status = "raw_used", candidates = candidates, coerced_clid = raw))
  }

  if (!is.na(candidate_id) && nzchar(candidate_id)) {
    return(list(clid = candidate_id, status = status, candidates = candidates, coerced_clid = candidate_id))
  }
  list(clid = NA_character_, status = status, candidates = candidates, coerced_clid = NA_character_)
}

map_quality_from_status <- function(status) {
  st <- as.character(status %||% "")
  if (!nzchar(st)) return(0L)
  if (st == "frozen_input_clid") return(3L)
  if (st == "raw_incompatible") return(1L)
  if (grepl("desc", st, fixed = TRUE)) return(2L)
  if (st %in% c("mapped_exact", "mapped_synonym")) return(3L)
  if (st %in% c("mapped_token")) return(1L)
  if (grepl("raw_used", st, fixed = TRUE)) return(3L)
  0L
}

normalize_score <- function(score) {
  s <- suppressWarnings(as.numeric(score))
  if (is.na(s)) return(NA_real_)
  if (s > 1) s <- s / 100
  max(min(s, 1), 0)
}

get_support_score <- function(x) {
  if (is.null(x)) return(NA_real_)
  score <- suppressWarnings(as.numeric(x$score %||% NA_real_))
  if (!is.na(score)) return(normalize_score(score))
  cs <- suppressWarnings(as.numeric(x$raw_inputs$score %||% NA_real_))
  if (!is.na(cs)) return(normalize_score(cs))
  cs2 <- suppressWarnings(as.numeric(x$raw_inputs$candidate_score_main %||% NA_real_))
  if (!is.na(cs2)) return(normalize_score(cs2))
  NA_real_
}

get_summary_top1_label <- function(summary) {
  summary$top1_cell_type %||% summary$predicted_label %||% summary$label %||% NA_character_
}

get_summary_clid <- function(summary) {
  summary$cell_ontology_id %||% summary$cl_id %||% NA_character_
}

get_summary_topk <- function(summary, top1_label) {
  topk <- summary$topk_cell_types %||% list()
  if (is.list(topk)) topk <- unlist(topk, recursive = TRUE, use.names = FALSE)
  topk <- as.character(topk)
  topk <- topk[!is.na(topk) & nzchar(topk)]
  if (!is.na(top1_label) && nzchar(top1_label)) topk <- c(top1_label, topk)
  unique(topk)
}

collect_method_candidates <- function(method_name, summary, cl_cfg, cl_graph, max_k = 3) {
  if (is.null(summary) || !is.list(summary)) return(list())
  top1_label <- get_summary_top1_label(summary)
  topk <- get_summary_topk(summary, top1_label)
  if (length(topk) > max_k) topk <- topk[seq_len(max_k)]
  base_score <- get_support_score(summary)
  out <- list()
  for (i in seq_along(topk)) {
    lbl <- topk[[i]]
    raw_clid <- if (i == 1) get_summary_clid(summary) else NA_character_

    # Rank-1 CL IDs from the fixed stage-08 registry input are authoritative.
    # Rank>1 textual hints may be normalized, but a valid rank-1 CL ID is not reinterpreted.
    raw_valid <- !is.na(raw_clid) && nzchar(as.character(raw_clid)) &&
      !is.null(cl_graph$cl[[as.character(raw_clid)]])
    if (i == 1L && isTRUE(raw_valid)) {
      coerced <- list(
        clid = as.character(raw_clid),
        status = "frozen_input_clid",
        candidates = as.character(raw_clid),
        coerced_clid = as.character(raw_clid)
      )
    } else {
      coerced <- coerce_clid(lbl, raw_clid, cl_cfg, cl_graph)
    }

    mq <- map_quality_from_status(coerced$status)
    sc <- if (!is.na(base_score)) max(base_score - 0.05 * (i - 1), 0) else NA_real_
    out[[length(out) + 1]] <- list(
      label = lbl,
      clid = coerced$clid,
      raw_clid = raw_clid,
      coerced_status = coerced$status,
      coerced_candidates = coerced$candidates,
      map_quality = mq,
      score = sc,
      method = method_name,
      rank = i
    )
  }
  out
}

get_anc_map <- function(clid, cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(NULL)
  if (is.na(clid) || !nzchar(clid)) return(NULL)
  term <- cl_graph$cl[[clid]] %||% NULL
  anc <- term$ancestors %||% NULL
  out <- c()
  out[clid] <- 0
  if (!is.null(anc)) {
    for (k in names(anc)) {
      v <- suppressWarnings(as.numeric(anc[[k]]))
      if (!is.na(v)) out[k] <- v
    }
  }
  out
}

CHILDREN_MAP_CACHE <- new.env(parent = emptyenv())

build_children_map <- function(cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(NULL)
  key <- "children_map"
  if (exists(key, envir = CHILDREN_MAP_CACHE, inherits = FALSE)) {
    return(get(key, envir = CHILDREN_MAP_CACHE, inherits = FALSE))
  }
  children_map <- list()
  for (child_id in names(cl_graph$cl)) {
    if (!startsWith(child_id, "CL:")) next
    anc <- cl_graph$cl[[child_id]]$ancestors %||% NULL
    if (is.null(anc)) next
    d <- unlist(anc, use.names = TRUE)
    if (length(d) == 0) next
    parents <- names(d[d == 1])
    if (length(parents) == 0) next
    for (p in parents) {
      if (is.null(children_map[[p]])) children_map[[p]] <- character(0)
      children_map[[p]] <- unique(c(children_map[[p]], child_id))
    }
  }
  assign(key, children_map, envir = CHILDREN_MAP_CACHE)
  children_map
}

get_descendants <- function(root_id, max_depth_k, cl_graph) {
  if (is.null(root_id) || !nzchar(root_id)) return(data.frame())
  if (max_depth_k <= 0) return(data.frame())
  children_map <- build_children_map(cl_graph)
  if (is.null(children_map)) return(data.frame())
  seen <- setNames(0, root_id)
  queue <- list(root_id)
  while (length(queue) > 0) {
    cur <- queue[[1]]
    queue <- queue[-1]
    cur_depth <- suppressWarnings(as.numeric(seen[cur]))
    if (is.na(cur_depth)) cur_depth <- 0
    if (cur_depth >= max_depth_k) next
    kids <- children_map[[cur]] %||% character(0)
    for (k in kids) {
      if (is.na(k) || !nzchar(k)) next
      if (is.na(seen[k])) {
        seen[k] <- cur_depth + 1
        queue[[length(queue) + 1]] <- k
      }
    }
  }
  out_ids <- names(seen)
  out_ids <- out_ids[out_ids != root_id]
  if (length(out_ids) == 0) return(data.frame())
  data.frame(clid = out_ids, depth_from_root = as.numeric(seen[out_ids]), stringsAsFactors = FALSE)
}

get_depth_from_anc <- function(anc_map) {
  if (is.null(anc_map) || length(anc_map) == 0) return(0)
  max(as.numeric(anc_map), na.rm = TRUE)
}

get_depth_to_root <- function(clid, cl_graph) {
  if (is.na(clid) || !nzchar(clid)) return(0)
  anc <- get_anc_map(clid, cl_graph)
  if (is.null(anc) || length(anc) == 0) return(0)
  max(as.numeric(anc), na.rm = TRUE)
}


# --------------------------
# Release policies: v2 blocker split and v3 automated-QC escalation
# --------------------------
normalize_issue_flags <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

valid_final_clid <- function(x, cl_graph = NULL) {
  x <- as.character(x %||% "")
  if (length(x) == 0L) return(FALSE)
  x <- trimws(x[[1]])
  if (!grepl("^CL:[0-9]+$", x)) return(FALSE)
  if (!is.null(cl_graph) && !is.null(cl_graph$cl) && is.null(cl_graph$cl[[x]])) return(FALSE)
  TRUE
}

release_policy_is_v2 <- function(x) {
  identical(tolower(trimws(as.character(x %||% "legacy"))), "blocker")
}

release_policy_is_v3 <- function(x) {
  identical(tolower(trimws(as.character(x %||% "legacy"))), "auto")
}

release_policy_is_modern <- function(x) {
  release_policy_is_v2(x) || release_policy_is_v3(x)
}

get_auto_qc_flags <- function(j) {
  if (is.null(j) || !is.list(j)) return(character(0))
  pi <- j$post_issues %||% list()
  normalize_issue_flags(pi$auto_qc_flags %||% character(0))
}

get_chief_qc_status <- function(j) {
  if (is.null(j) || !is.list(j)) return("not_run")
  pi <- j$post_issues %||% list()
  x <- tolower(trimws(as.character(pi$chief_qc_status %||% "not_run")))
  if (!x %in% c("not_run", "pending", "passed", "failed")) x <- "not_run"
  x
}

# Automated-first release policy:
#   audit_flags       = information only
#   auto_qc_flags     = automatically escalate to Chief, not a human
#   release_blockers  = human review only after deterministic/Chief checks remain unresolved
apply_release_policy_v3 <- function(j, cl_graph = NULL) {
  if (is.null(j) || !is.list(j)) return(j)
  j <- ensure_post_issues(j)

  flags <- normalize_issue_flags(j$post_issues$flags %||% character(0))
  fd <- j$final_decision %||% list()
  dt <- j$decision_trace %||% list()

  final_label <- trimws(as.character(fd$primary_cell_type %||% ""))
  final_clid <- trimws(as.character(fd$final_cell_ontology_id %||% ""))
  confidence <- suppressWarnings(as.numeric(fd$confidence_primary %||% NA_real_))
  state <- tolower(trimws(as.character(j$cluster_state %||% "clean")))
  decision_cat <- normalize_decision_category(fd$decision_category %||% "")
  final_rule <- trimws(as.character(dt$final_rule %||% ""))
  support_gate_pass <- isTRUE(dt$support_gate_pass %||% FALSE)
  final_valid <- valid_final_clid(final_clid, cl_graph) && nzchar(final_label)
  final_depth <- if (final_valid && !is.null(cl_graph)) get_depth_to_root(final_clid, cl_graph) else NA_real_
  clear_decision <- decision_cat %in% c(
    "cassia_better", "in_house_better", "enrich_better",
    "third_party_override", "tie"
  )

  dominant_identity <- isTRUE(final_valid) &&
    is.finite(confidence) && confidence >= 0.65 &&
    nzchar(final_rule) && clear_decision

  supported_identity <- isTRUE(dominant_identity) &&
    (support_gate_pass || confidence >= 0.80)

  strong_identity <- isTRUE(dominant_identity) &&
    support_gate_pass && confidence >= 0.75

  chief_status <- get_chief_qc_status(j)
  chief_passed <- identical(chief_status, "passed")
  chief_failed <- identical(chief_status, "failed")

  audit_flags <- character(0)
  auto_qc_flags <- character(0)
  blockers <- character(0)

  add_audit <- function(x) audit_flags <<- unique(c(audit_flags, x))
  add_auto <- function(x) auto_qc_flags <<- unique(c(auto_qc_flags, x))
  add_blocker <- function(x) blockers <<- unique(c(blockers, x))

  always_audit <- c(
    "score_missing", "score_unavailable", "rerun_not_beneficial",
    "mixed_or_ambiguous_signals", "final_repaired_from_head",
    "dataset_tissue_missing", "llm_tissue_classifier_missing_keep_all",
    "generic_but_lock_stable", "ontology_id_mismatch_removed",
    "parity_avoid_generic_fallback", "parity_lineage_descendant_preference",
    "lock_confident_freeze", "decision_category_repaired",
    "deterministic_rerank_disagreed_with_head",
    "deterministic_hierarchical_adjustment_from_head",
    "deterministic_hierarchical_adjustment_blocked",
    "deterministic_identity_unavailable_head_preserved"
  )

  hard_structural <- c(
    "no_lock", "lock_missing", "lock_no_in_tree",
    "head_missing", "max_rounds_exceeded",
    "missing_or_invalid_clid", "missing_clid", "missing_primary_cell_type",
    "manual_review_plan_missing", "manual_review_actions_missing",
    "manual_review_action_evidence_pointers_missing",
    "third_party_adjudication_missing", "third_party_primary_cell_type_missing",
    "third_party_evidence_pointers_missing"
  )

  soft_signal_pattern <- paste(
    c(
      "^minor_",
      "activation_signal$",
      "activated_state$",
      "secondary_signal$",
      "squamous_differentiation_signal$",
      "cytotoxic_activated_state$",
      "minor_contamination_signal$",
      "minor_neutrophil_like_signal$"
    ),
    collapse = "|"
  )

  contamination_pattern <- paste(
    c(
      "possible_contamination",
      "contamination_signal",
      "contamination_suspected",
      "possible_ambient_contamination",
      "ambient_contamination",
      "possible_doublet",
      "doublet_signal"
    ),
    collapse = "|"
  )

  subtype_pattern <- paste(
    c(
      "subtype_conflict",
      "subtype_uncertainty",
      "subtype_ambiguity",
      "subtype_unresolved"
    ),
    collapse = "|"
  )

  for (flag in flags) {
    if (flag %in% always_audit) {
      add_audit(if (identical(flag, "score_missing")) "score_unavailable" else flag)

    } else if (grepl(soft_signal_pattern, flag, ignore.case = TRUE, perl = TRUE)) {
      # Minor/secondary state signals do not block a stable dominant identity.
      if (dominant_identity) add_audit(paste0(flag, "_documented"))
      else add_auto(flag)

    } else if (flag %in% hard_structural) {
      add_blocker(flag)

    } else if (flag %in% c("chief_qc_failed", "chief_qc_failed_but_local_gate_passed")) {
      if (chief_passed) add_audit("prior_chief_qc_failure_resolved")
      else add_blocker(flag)

    } else if (flag %in% c("missing_required_citation")) {
      # Citation completion is an automated QC/repair task first.
      add_auto(flag)

    } else if (flag %in% c("head_identity_outside_reviewer_candidate_pool",
                           "head_vs_majority_cross_branch_qc",
                           "head_vs_deterministic_cross_branch_qc",
                           "consensus_override_requires_qc")) {
      # Non-hierarchical identity disagreement is never
      # silently rewritten; automated Chief QC resolves the biological conflict.
      add_auto(flag)

    } else if (flag %in% c("ontology_id_mismatch", "reviewer_disagreement")) {
      # This flag describes input disagreement, not necessarily a bad final ID.
      if (dominant_identity) add_audit("input_ontology_disagreement_resolved")
      else add_auto("ontology_disagreement_requires_qc")

    } else if (flag %in% c("cross_lineage_conflict", "lineage_conflict",
                           "unresolved_cross_lineage_conflict")) {
      if (strong_identity) add_audit("cross_lineage_disagreement_resolved")
      else if (dominant_identity) add_auto("cross_lineage_resolution_check")
      else add_blocker("no_dominant_identity_cross_lineage")

    } else if (flag %in% c("sibling_conflict", "unresolved_sibling_conflict")) {
      if (strong_identity) add_audit("sibling_disagreement_resolved")
      else add_auto("sibling_resolution_check")

    } else if (flag %in% c("weak_score", "evidence_weak",
                           "weak_support_unresolved")) {
      if (supported_identity) add_audit("weak_auxiliary_support_resolved")
      else add_auto("weak_support_requires_qc")

    } else if (flag %in% c("small_margin", "small_margin_unresolved")) {
      if (strong_identity) add_audit("small_margin_resolved")
      else add_auto("small_margin_requires_qc")

    } else if (flag %in% c("high_level_lock", "unsupported_specificity")) {
      if (dominant_identity && is.finite(final_depth) && final_depth >= 3) {
        add_audit("high_level_lock_resolved_by_specific_final")
      } else {
        add_auto("specificity_requires_qc")
      }

    } else if (flag %in% c("mixed_signature", "unresolved_mixed_signature")) {
      if (state == "clean" && dominant_identity) {
        add_audit("minor_mixed_signal_resolved")
      } else if (dominant_identity) {
        add_auto("mixed_signal_requires_qc")
      } else {
        add_blocker("no_dominant_identity_mixed_signal")
      }

    } else if (grepl(contamination_pattern, flag, ignore.case = TRUE, perl = TRUE)) {
      # Contamination/ambient signals are often secondary. Route to automated QC.
      if (state == "clean" && dominant_identity) add_auto(flag)
      else if (dominant_identity) add_auto(paste0(flag, "_state_check"))
      else add_blocker("no_dominant_identity_contamination")

    } else if (grepl(subtype_pattern, flag, ignore.case = TRUE, perl = TRUE)) {
      # Prefer conservative automatic resolution; human review only if no stable identity.
      if (strong_identity) add_audit(paste0(flag, "_resolved"))
      else if (dominant_identity) add_auto(paste0(flag, "_qc"))
      else add_blocker("no_dominant_identity_subtype")

    } else {
      # Unknown flags are not silently released and are not immediately sent to a human.
      # They enter automated Chief QC first.
      add_auto(paste0("unclassified_flag:", flag))
    }
  }

  # Deterministic validity checks.
  if (!nzchar(final_label)) add_blocker("missing_primary_cell_type")
  if (!valid_final_clid(final_clid, cl_graph)) add_blocker("missing_or_invalid_final_clid")

  if (!is.finite(confidence) || confidence < 0.5) {
    add_blocker("low_confidence")
  } else if (confidence < 0.65) {
    add_auto("moderate_low_confidence")
  }

  if (state %in% c("mixed", "doublet_suspected", "contaminated_or_ambient")) {
    if (dominant_identity) {
      add_auto(paste0("cluster_state_", state, "_requires_qc"))
    } else {
      add_blocker(paste0("no_dominant_identity_cluster_state_", state))
    }
  } else if (identical(state, "transitional_state")) {
    if (dominant_identity) add_audit("transitional_state_with_dominant_identity")
    else add_auto("transitional_state_requires_qc")
  }

  if (identical(decision_cat, "third_party_override") &&
      (!is.finite(confidence) || confidence < 0.65)) {
    add_auto("third_party_override_requires_qc")
  }

  # Chief is the automated escalation layer.
  if (chief_failed) add_blocker("chief_qc_failed")

  if (chief_passed && length(auto_qc_flags) > 0L) {
    add_audit(paste0("chief_qc_passed:", auto_qc_flags))
    auto_qc_flags <- character(0)
  }

  blockers <- unique(blockers[nzchar(blockers)])
  auto_qc_flags <- unique(setdiff(auto_qc_flags[nzchar(auto_qc_flags)], blockers))
  audit_flags <- unique(setdiff(
    audit_flags[nzchar(audit_flags)],
    c(blockers, auto_qc_flags)
  ))

  # Human review is only for unresolved hard blockers.
  needs_review <- length(blockers) > 0L

  j$post_issues$flags <- as.list(flags)
  j$post_issues$audit_flags <- as.list(audit_flags)
  j$post_issues$auto_qc_flags <- as.list(auto_qc_flags)
  j$post_issues$release_blockers <- as.list(blockers)
  j$post_issues$needs_manual_review <- needs_review

  j$post_issues$release_state <- if (needs_review) {
    "manual_review"
  } else if (length(auto_qc_flags) > 0L) {
    "auto_qc_pending"
  } else if (length(audit_flags) > 0L) {
    "release_with_audit_note"
  } else {
    "release"
  }

  j$post_issues$policy_version <- "auto"
  j
}

get_release_blockers <- function(j) {
  if (is.null(j) || !is.list(j)) return(character(0))
  pi <- j$post_issues %||% list()
  normalize_issue_flags(pi$release_blockers %||% character(0))
}

get_audit_flags <- function(j) {
  if (is.null(j) || !is.list(j)) return(character(0))
  pi <- j$post_issues %||% list()
  normalize_issue_flags(pi$audit_flags %||% character(0))
}

apply_release_policy <- function(j, cl_graph = NULL, release_policy = "legacy") {
  if (is.null(j) || !is.list(j)) return(j)
  j <- ensure_post_issues(j)

  if (release_policy_is_v3(release_policy)) {
    return(apply_release_policy_v3(j, cl_graph = cl_graph))
  }

  if (!release_policy_is_v2(release_policy)) {
    j$post_issues$policy_version <- "legacy"
    j$post_issues$release_state <- if (isTRUE(j$post_issues$needs_manual_review %||% FALSE)) "manual_review" else "release"
    return(j)
  }

  flags <- normalize_issue_flags(j$post_issues$flags %||% character(0))
  fd <- j$final_decision %||% list()
  dt <- j$decision_trace %||% list()

  final_label <- trimws(as.character(fd$primary_cell_type %||% ""))
  final_clid <- trimws(as.character(fd$final_cell_ontology_id %||% ""))
  confidence <- suppressWarnings(as.numeric(fd$confidence_primary %||% NA_real_))
  state <- tolower(trimws(as.character(j$cluster_state %||% "clean")))
  decision_cat <- normalize_decision_category(fd$decision_category %||% "")
  final_rule <- trimws(as.character(dt$final_rule %||% ""))
  support_gate_pass <- isTRUE(dt$support_gate_pass %||% FALSE)
  final_valid <- valid_final_clid(final_clid, cl_graph) && nzchar(final_label)
  final_depth <- if (final_valid && !is.null(cl_graph)) get_depth_to_root(final_clid, cl_graph) else NA_real_
  clear_decision <- decision_cat %in% c("cassia_better", "in_house_better", "enrich_better", "third_party_override", "tie")
  resolved_core <- isTRUE(final_valid) && is.finite(confidence) && confidence >= 0.70 && nzchar(final_rule) && clear_decision
  resolved_supported <- isTRUE(resolved_core) && support_gate_pass
  resolved_strong <- isTRUE(resolved_supported) && confidence >= 0.75

  audit_flags <- character(0)
  blockers <- character(0)
  add_audit <- function(x) audit_flags <<- unique(c(audit_flags, x))
  add_blocker <- function(x) blockers <<- unique(c(blockers, x))

  always_audit <- c(
    "score_missing", "score_unavailable", "rerun_not_beneficial",
    "mixed_or_ambiguous_signals", "final_repaired_from_head",
    "dataset_tissue_missing", "llm_tissue_classifier_missing_keep_all",
    "generic_but_lock_stable", "ontology_id_mismatch_removed",
    "parity_avoid_generic_fallback", "parity_lineage_descendant_preference",
    "lock_confident_freeze"
  )
  always_block <- c(
    "no_lock", "lock_missing", "lock_no_in_tree", "subtype_conflict",
    "head_missing", "max_rounds_exceeded", "chief_qc_failed",
    "chief_qc_failed_but_local_gate_passed", "missing_required_citation",
    "missing_or_invalid_clid", "missing_clid", "missing_primary_cell_type",
    "manual_review_plan_missing", "manual_review_actions_missing",
    "manual_review_action_evidence_pointers_missing",
    "third_party_adjudication_missing", "third_party_primary_cell_type_missing",
    "third_party_evidence_pointers_missing"
  )

  for (flag in flags) {
    if (flag %in% always_audit) {
      add_audit(if (identical(flag, "score_missing")) "score_unavailable" else flag)
    } else if (flag %in% always_block) {
      add_blocker(flag)
    } else if (flag %in% c("ontology_id_mismatch", "reviewer_disagreement")) {
      if (resolved_core) add_audit("input_ontology_disagreement_resolved") else add_blocker("unresolved_ontology_disagreement")
    } else if (flag %in% c("cross_lineage_conflict", "lineage_conflict")) {
      if (resolved_strong) add_audit("cross_lineage_disagreement_resolved") else add_blocker("unresolved_cross_lineage_conflict")
    } else if (identical(flag, "sibling_conflict")) {
      if (resolved_strong) add_audit("sibling_disagreement_resolved") else add_blocker("unresolved_sibling_conflict")
    } else if (flag %in% c("weak_score", "evidence_weak")) {
      if (resolved_supported) add_audit("weak_auxiliary_support_resolved") else add_blocker("weak_support_unresolved")
    } else if (identical(flag, "small_margin")) {
      if (resolved_strong) add_audit("small_margin_resolved") else add_blocker("small_margin_unresolved")
    } else if (identical(flag, "high_level_lock")) {
      if (resolved_core && is.finite(final_depth) && final_depth >= 3) {
        add_audit("high_level_lock_resolved_by_specific_final")
      } else {
        add_blocker("unsupported_specificity")
      }
    } else if (identical(flag, "mixed_signature")) {
      if (state == "clean" && resolved_strong) add_audit("minor_mixed_signal_resolved") else add_blocker("unresolved_mixed_signature")
    } else {
      # Unknown flags remain conservative until explicitly classified.
      add_blocker(flag)
    }
  }

  # Deterministic blockers independent of model-generated flags.
  if (!nzchar(final_label)) add_blocker("missing_primary_cell_type")
  if (!valid_final_clid(final_clid, cl_graph)) add_blocker("missing_or_invalid_final_clid")
  if (!is.finite(confidence) || confidence < 0.5) add_blocker("low_confidence")

  if (state %in% c("mixed", "doublet_suspected", "contaminated_or_ambient")) {
    add_blocker(paste0("unresolved_cluster_state_", state))
  } else if (identical(state, "transitional_state")) {
    if (resolved_core) add_audit("transitional_state_with_dominant_identity") else add_blocker("unresolved_transitional_state")
  }

  if (identical(decision_cat, "third_party_override") && (!is.finite(confidence) || confidence < 0.65)) {
    add_blocker("low_confidence_third_party_override")
  }

  blockers <- unique(blockers[nzchar(blockers)])
  audit_flags <- unique(setdiff(audit_flags[nzchar(audit_flags)], blockers))
  needs_review <- length(blockers) > 0L

  j$post_issues$flags <- as.list(flags)
  j$post_issues$audit_flags <- as.list(audit_flags)
  j$post_issues$release_blockers <- as.list(blockers)
  j$post_issues$needs_manual_review <- needs_review
  j$post_issues$release_state <- if (needs_review) {
    "manual_review"
  } else if (length(audit_flags) > 0L) {
    "release_with_audit_note"
  } else {
    "release"
  }
  j$post_issues$policy_version <- "blocker"
  j
}

get_msca_two <- function(clid_a, clid_b, cl_graph) {
  am_a <- get_anc_map(clid_a, cl_graph)
  am_b <- get_anc_map(clid_b, cl_graph)
  if (is.null(am_a) || is.null(am_b)) return(list(id = NA_character_, depth = NA_real_))
  common <- intersect(names(am_a), names(am_b))
  if (length(common) == 0) return(list(id = NA_character_, depth = NA_real_))
  dist_sum <- vapply(common, function(id) as.numeric(am_a[[id]]) + as.numeric(am_b[[id]]), numeric(1))
  best_id <- common[which.min(dist_sum)][[1]]
  list(id = best_id, depth = get_depth_from_anc(get_anc_map(best_id, cl_graph)))
}

get_parent_clid <- function(clid, cl_graph) {
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(NA_character_)
  if (is.na(clid) || !nzchar(clid)) return(NA_character_)
  term <- cl_graph$cl[[clid]] %||% NULL
  anc <- term$ancestors %||% NULL
  if (is.null(anc)) return(NA_character_)
  parent <- names(anc)[which(as.numeric(anc) == 1)]
  if (length(parent) == 0) return(NA_character_)
  parent[[1]]
}

is_developmental_stage <- function(clid = NA_character_, label = NA_character_, cl_graph = NULL, cl_cfg = NULL, stage_root_clids = character(0), allow_token_fallback = FALSE) {
  if (!is.na(clid) && nzchar(as.character(clid)) && length(stage_root_clids) > 0 && !is.null(cl_graph)) {
    x <- as.character(clid)
    for (r in stage_root_clids) {
      if (is.na(r) || !nzchar(r)) next
      if (identical(x, r) || is_descendant_of(x, r, cl_graph)) return(TRUE)
    }
    return(FALSE)
  }
  if (!isTRUE(allow_token_fallback)) return(FALSE)
  s <- tolower(stringr::str_squish(as.character(label %||% "")))
  if (!nzchar(s)) return(FALSE)
  stringr::str_detect(s, "progenitor|precursor|stem|blast")
}

normalize_judge_final_decision_cl_preserve <- function(final_obj, cl_cfg) {
  if (is.null(final_obj) || !is.list(final_obj)) return(final_obj)
  if (!is.list(final_obj$final_decision)) return(final_obj)
  fd <- final_obj$final_decision
  res <- normalize_cl_three_state(fd$primary_cell_type %||% "", fd$final_cell_ontology_id %||% "", cl_cfg)
  if (is.null(fd$final_cell_ontology_id) || !nzchar(fd$final_cell_ontology_id)) {
    fd$final_cell_ontology_id <- res$final_clid %||% fd$final_cell_ontology_id
  }
  if (is.null(fd$primary_cell_type) || !nzchar(fd$primary_cell_type)) {
    fd$primary_cell_type <- res$final_name %||% fd$primary_cell_type
  }
  if (!is.null(fd$greedy_cell_type) || !is.null(fd$greedy_cell_ontology_id)) {
    resg <- normalize_cl_three_state(fd$greedy_cell_type %||% "", fd$greedy_cell_ontology_id %||% "", cl_cfg)
    if (is.null(fd$greedy_cell_ontology_id) || !nzchar(fd$greedy_cell_ontology_id)) {
      fd$greedy_cell_ontology_id <- resg$final_clid %||% fd$greedy_cell_ontology_id
    }
    if (is.null(fd$greedy_cell_type) || !nzchar(fd$greedy_cell_type)) {
      fd$greedy_cell_type <- resg$final_name %||% fd$greedy_cell_type
    }
  }
  final_obj$final_decision <- fd
  final_obj
}

finalize_output <- function(out, jin, cl_cfg, cl_graph, cid = NULL, out_dir_head = NULL, head_out_last = NULL,
                            release_policy = "legacy") {
  if (is.null(out) || !is.list(out)) return(out)
  if (!is.null(cid) && nzchar(as.character(cid))) out$cluster_id <- cid
  out <- set_method_correctness_vs_final(out, jin, cl_graph)
  out <- standardize_transition_notes(out)
  out <- ensure_audit_report_support(out)
  out <- ensure_enrich_verdict(out, jin)
  out <- ensure_in_house_alias(out)
  out <- normalize_judge_final_decision_cl_preserve(out, cl_cfg)
  if (!is.null(out_dir_head) && !is.null(head_out_last)) {
    out <- hydrate_final_from_head(out, out$cluster_id %||% cid, head_out_last, out_dir_head)
  }
  out <- apply_release_policy(out, cl_graph = cl_graph, release_policy = release_policy)
  out <- normalize_manual_review_plan(out)
  if (isTRUE(out$post_issues$needs_manual_review %||% FALSE)) {
    out <- ensure_manual_review_action(out, jin, extract_citation_requirements(jin))
  }
  out <- force_array_fields(out)
  out
}

postprocess_head_out <- function(head_out, jin, cite_req, cl_cfg, cl_graph, dataset_cfg = NULL, species_value = "human",
                                 min_depth_for_eligibility = 3,
                                 delta_depth = 2,
                                 gate_mode = "fixed_k",
                                 gate_k = 2,
                                 disable_aggressive_trigger = FALSE,
                                 enable_method_reliability_gate = FALSE,
                                 enable_aggressive_lca_hard_guard = TRUE,
                                 protect_lock_inputs = FALSE,
                                 meta_reviewer_gate = FALSE,
                                 parity_mode = FALSE,
                                 identity_policy = "ontology_guarded",
                                 release_policy = "legacy") {
  if (is.null(head_out) || !is.list(head_out)) return(head_out)
  normalize_stage <- function(x) {
    x <- normalize_manual_review_plan(x)
    x <- ensure_audit_report_support(x)
    x <- filter_citations_to_allowlist(x, cite_req)
    x <- ensure_manual_review_action(x, jin, cite_req)
    x <- force_array_fields(x)
    x <- ensure_enrich_verdict(x, jin)
    x <- ensure_in_house_alias(x)
    x <- normalize_judge_final_decision_cl_preserve(x, cl_cfg)
    x
  }

  choose_stage <- function(x) {
    apply_consensus_subtype_policy(
      x, jin, cl_cfg, cl_graph,
      dataset_cfg = dataset_cfg,
      species_value = species_value,
      min_depth_for_eligibility = min_depth_for_eligibility,
      delta_depth = delta_depth,
      gate_mode = gate_mode,
      gate_k = gate_k,
      disable_aggressive_trigger = disable_aggressive_trigger,
      enable_method_reliability_gate = enable_method_reliability_gate,
      enable_aggressive_lca_hard_guard = enable_aggressive_lca_hard_guard,
      protect_lock_inputs = protect_lock_inputs,
      meta_reviewer_gate = meta_reviewer_gate,
      parity_mode = parity_mode,
      identity_policy = identity_policy
    )
  }

  head_out <- normalize_stage(head_out)
  head_out <- choose_stage(head_out)
  head_out <- apply_release_policy(head_out, cl_graph = cl_graph, release_policy = release_policy)
  head_out <- normalize_manual_review_plan(head_out)
  needs_review <- isTRUE(head_out$post_issues$needs_manual_review %||% FALSE)
  if (needs_review) {
    head_out <- ensure_manual_review_action(head_out, jin, cite_req)
  }
  head_out <- force_array_fields(head_out)
  head_out
}


# --------------------------
# Reviewer-anchor + evidence adjudication policy (generalized, no hardcoded labels)
# --------------------------

normalize_summary_method_key <- function(k) {
  mk <- sub("_summary$", "", as.character(k))
  if (mk %in% c("our", "our_method", "in_house", "in_house_method")) mk <- "our"
  if (mk %in% c("inter", "enrich", "enrichment", "enrichment_method")) mk <- "enrich"
  mk
}

is_valid_method_summary <- function(x) {
  if (is.null(x) || !is.list(x)) return(FALSE)
  has_top <- !is.null(x$top1_cell_type) || !is.null(x$predicted_label)
  has_any <- has_top || !is.null(x$topk_cell_types)
  isTRUE(has_any)
}

get_method_summary_for_frozen <- function(inputs, method) {
  method <- normalize_summary_method_key(method)
  if (identical(method, "cassia")) return(inputs$cassia_summary %||% NULL)
  if (identical(method, "our")) return(inputs$in_house_summary %||% inputs$our_summary %||% NULL)
  if (identical(method, "enrich")) return(inputs$enrich_summary %||% inputs$inter_summary %||% NULL)
  NULL
}

collect_frozen_candidates_from_inputs <- function(inputs, cl_graph) {
  frozen <- inputs$frozen_reviewer_candidates %||% list()
  if (!is.list(frozen) || length(frozen) == 0L) return(list())

  out <- list()
  seen <- character(0)
  for (rec in frozen) {
    if (is.null(rec) || !is.list(rec)) next
    method <- normalize_summary_method_key(rec$method %||% "")
    if (!method %in% c("cassia", "our", "enrich")) next

    rank <- suppressWarnings(as.integer(rec$rank %||% 99L))
    label <- as.character(rec$label %||% NA_character_)
    clid <- as.character(rec$clid %||% NA_character_)
    if (is.na(rank) || rank < 1L || is.na(label) || !nzchar(label)) next
    if (is.na(clid) || !nzchar(clid) || is.null(cl_graph$cl[[clid]])) next

    key <- paste(method, rank, clid, sep = "|")
    if (key %in% seen) next
    seen <- c(seen, key)

    # Preserve the per-method score behavior. The ontology ID itself is fixed
    # in Step 08; Step 09 does not reinterpret/re-map the TOP-K label.
    sm <- get_method_summary_for_frozen(inputs, method)
    base_score <- get_support_score(sm %||% list())
    sc <- if (!is.na(base_score)) max(base_score - 0.05 * (rank - 1L), 0) else NA_real_

    mq <- suppressWarnings(as.integer(rec$map_quality %||% NA_integer_))
    if (is.na(mq)) mq <- if (rank == 1L) 3L else 2L

    out[[length(out) + 1L]] <- list(
      label = label,
      clid = clid,
      raw_clid = clid,
      coerced_status = as.character(rec$mapping_status %||% rec$mapping_source %||% "frozen_step08"),
      coerced_candidates = clid,
      map_quality = mq,
      score = sc,
      method = method,
      rank = rank,
      candidate_mapping_source = as.character(rec$mapping_source %||% "frozen_step08")
    )
  }
  out
}

collect_candidates_from_inputs <- function(inputs, cl_cfg, cl_graph) {
  frozen_payload <- inputs$frozen_reviewer_candidates %||% list()
  if (is.list(frozen_payload) && length(frozen_payload) > 0L) {
    # Formal mapping-safe path: mapped records participate in ontology logic;
    # unmapped records remain visible in reviewer summaries/reports but are never
    # silently remapped in Step 09.
    return(collect_frozen_candidates_from_inputs(inputs, cl_graph))
  }

  # Compatibility fallback for judge inputs that do not contain the fixed payload.
  summary_keys <- names(inputs)
  summary_keys <- summary_keys[grepl("_summary$", summary_keys)]
  summary_keys <- summary_keys[!summary_keys %in% c("evidence_summary")]
  canonical_priority <- c(
    cassia_summary = 1,
    in_house_summary = 2,
    our_summary = 3,
    enrich_summary = 4,
    inter_summary = 5
  )
  summary_keys <- summary_keys[order(vapply(summary_keys, function(k) {
    p <- canonical_priority[k] %||% 99L
    as.integer(p)
  }, integer(1)))]
  seen_methods <- character(0)
  candidates <- list()
  for (k in summary_keys) {
    sm <- inputs[[k]]
    if (!is_valid_method_summary(sm)) next
    mkey <- normalize_summary_method_key(k)
    if (mkey %in% seen_methods) next
    seen_methods <- c(seen_methods, mkey)
    candidates <- c(candidates, collect_method_candidates(mkey, sm, cl_cfg, cl_graph))
  }
  candidates
}

validate_frozen_reviewer_candidate_payload <- function(jin, cl_graph, path = "") {
  inputs <- jin$inputs %||% list()
  frozen <- inputs$frozen_reviewer_candidates %||% list()
  if (!is.list(frozen) || length(frozen) == 0L) {
    stop("Missing inputs$frozen_reviewer_candidates in judge input",
         if (nzchar(path)) paste0(": ", path) else "",
         ". Regenerate Step 08 with 08_build_judge_inputs_v4_3_3_frozen_topk.R.")
  }

  # Determine the reviewer set that should be present. M-series inputs explicitly
  # carry reviewer_availability; standard inputs default to all three reviewers.
  avail <- as.character(unlist(
    inputs$reviewer_availability$available_reviewers %||% c("cassia", "in_house", "enrichment"),
    recursive = TRUE, use.names = FALSE
  ))
  avail <- vapply(avail, function(x) {
    k <- normalize_summary_method_key(x)
    if (identical(k, "in_house")) "our" else k
  }, character(1))
  avail <- unique(avail[avail %in% c("cassia", "our", "enrich")])

  methods <- vapply(frozen, function(x) normalize_summary_method_key(x$method %||% ""), character(1))
  methods <- unique(methods[methods %in% c("cassia", "our", "enrich")])
  if (!setequal(methods, avail)) {
    stop("Frozen reviewer-candidate method set mismatch",
         if (nzchar(path)) paste0(" in ", path) else "",
         ": found=", paste(sort(methods), collapse = ","),
         " expected=", paste(sort(avail), collapse = ","))
  }

  for (m in avail) {
    top1 <- Filter(function(x) {
      identical(normalize_summary_method_key(x$method %||% ""), m) &&
        as.integer(x$rank %||% 99L) == 1L
    }, frozen)
    if (length(top1) != 1L) {
      stop("Frozen reviewer-candidate payload must contain exactly one TOP1 for ", m,
           if (nzchar(path)) paste0(" in ", path) else "")
    }
    cid <- as.character(top1[[1]]$clid %||% NA_character_)

    # Missing CLID is a legitimate mapping outcome (ambiguous/unmapped), not a
    # malformed judge input. Such a reviewer remains visible to Head/Chief as text
    # evidence, but cannot vote in CL-ID consensus or ontology relations.
    if (!is.na(cid) && nzchar(cid) && is.null(cl_graph$cl[[cid]])) {
      stop("Frozen TOP1 contains a non-empty CLID that is absent from the frozen ontology for ", m,
           if (nzchar(path)) paste0(" in ", path) else "",
           ": ", cid)
    }
  }
  invisible(TRUE)
}

compute_top1_lock_context <- function(candidates, enable_method_reliability_gate, protect_lock_inputs) {
  top1_candidates <- Filter(function(x) x$rank == 1, candidates)
  top1_methods <- unique(vapply(top1_candidates, function(x) as.character(x$method %||% ""), character(1)))
  top1_methods <- top1_methods[nzchar(top1_methods)]
  primary_methods <- intersect(c("cassia", "our", "enrich"), top1_methods)
  lock_source_methods <- top1_methods
  if (isTRUE(enable_method_reliability_gate) && isTRUE(protect_lock_inputs) && length(primary_methods) > 0) {
    lock_source_methods <- primary_methods
  }
  lock_source <- if (length(lock_source_methods) > 0) paste(lock_source_methods, collapse = ",") else "none"
  top1_for_lock <- Filter(function(x) (x$method %||% "") %in% lock_source_methods, top1_candidates)

  # Defensive one-vote-per-method deduplication. Aliases are canonicalized upstream
  # (in_house/our -> our; inter/enrichment -> enrich), so duplicate summaries cannot vote twice.
  if (length(top1_for_lock) > 0L) {
    m <- vapply(top1_for_lock, function(x) as.character(x$method %||% ""), character(1))
    top1_for_lock <- top1_for_lock[!duplicated(m)]
  }

  n_source_methods <- length(unique(vapply(top1_for_lock, function(x) as.character(x$method %||% ""), character(1))))
  valid_top1 <- Filter(function(x) {
    clid <- as.character(x$clid %||% NA_character_)
    method <- as.character(x$method %||% "")
    !is.na(clid) && nzchar(clid) && nzchar(method)
  }, top1_for_lock)

  top1_clids <- unique(vapply(valid_top1, function(x) as.character(x$clid), character(1)))
  top1_clids <- top1_clids[!is.na(top1_clids) & nzchar(top1_clids)]

  direct_consensus_clid <- NA_character_
  direct_consensus_votes <- 0L
  direct_consensus_fraction <- NA_real_
  direct_consensus_reason <- ""

  if (length(valid_top1) > 0L && n_source_methods > 0L) {
    clid_counts <- sort(table(vapply(valid_top1, function(x) as.character(x$clid), character(1))), decreasing = TRUE)
    if (length(clid_counts) > 0L) {
      top_n <- as.integer(clid_counts[[1]])
      second_n <- if (length(clid_counts) >= 2L) as.integer(clid_counts[[2]]) else 0L
      # TRUE majority: >=2 independent reviewers, >50% of source methods, unique winner.
      # Thus 1/3 and 2/4 are not treated as direct consensus.
      majority_ok <- top_n >= 2L && top_n > (n_source_methods / 2) && top_n > second_n
      if (isTRUE(majority_ok)) {
        direct_consensus_clid <- names(clid_counts)[[1]]
        direct_consensus_votes <- top_n
        direct_consensus_fraction <- top_n / n_source_methods
        direct_consensus_reason <- paste0("direct_top1_majority_", top_n, "of", n_source_methods)
      }
    }
  }

  list(
    top1_candidates = top1_candidates,
    top1_methods = top1_methods,
    primary_methods = primary_methods,
    lock_source_methods = lock_source_methods,
    lock_source = lock_source,
    top1_for_lock = top1_for_lock,
    top1_clids = top1_clids,
    n_source_methods = as.integer(n_source_methods),
    direct_consensus_clid = direct_consensus_clid,
    direct_consensus_votes = as.integer(direct_consensus_votes),
    direct_consensus_fraction = as.numeric(direct_consensus_fraction),
    direct_consensus_reason = direct_consensus_reason
  )
}

apply_meta_reviewer_gate <- function(candidates, candidates_raw, top1_methods, method_reliability, method_alignments, meta_reviewer_gate) {
  meta_reviewer_decisions <- list()
  if (!isTRUE(meta_reviewer_gate) || length(candidates) == 0) {
    top1_candidates <- Filter(function(x) x$rank == 1, candidates)
    return(list(
      candidates = candidates,
      top1_candidates = top1_candidates,
      top1_methods = top1_methods,
      meta_reviewer_decisions = meta_reviewer_decisions
    ))
  }

  reliability_min <- 0.55
  alignment_min <- 0.45
  core_methods <- c("cassia", "our", "enrich")
  methods_seen <- unique(vapply(candidates, function(x) as.character(x$method %||% ""), character(1)))
  methods_seen <- methods_seen[nzchar(methods_seen)]
  non_core <- setdiff(methods_seen, core_methods)
  keep_non_core <- character(0)

  for (m in non_core) {
    rel <- suppressWarnings(as.numeric(method_reliability[[m]]$score %||% 0))
    if (!is.finite(rel)) rel <- 0
    align <- suppressWarnings(as.numeric(method_alignments[[m]]$score %||% NA_real_))
    if (!is.finite(align)) align <- rel
    keep <- isTRUE(rel >= reliability_min && align >= alignment_min)
    meta_reviewer_decisions[[length(meta_reviewer_decisions) + 1]] <- list(
      method = m,
      keep = keep,
      reliability_score = rel,
      alignment_score = align,
      reliability_min = reliability_min,
      alignment_min = alignment_min
    )
    if (isTRUE(keep)) keep_non_core <- c(keep_non_core, m)
  }

  keep_methods <- unique(c(core_methods, keep_non_core))
  candidates_filtered <- Filter(function(x) (x$method %||% "") %in% keep_methods, candidates)
  candidates_out <- candidates
  if (length(candidates_filtered) > 0) {
    candidates_out <- candidates_filtered
  } else {
    candidates_out <- candidates_raw
  }
  top1_candidates_out <- Filter(function(x) x$rank == 1, candidates_out)
  top1_methods_out <- unique(vapply(top1_candidates_out, function(x) as.character(x$method %||% ""), character(1)))
  top1_methods_out <- top1_methods_out[nzchar(top1_methods_out)]

  list(
    candidates = candidates_out,
    top1_candidates = top1_candidates_out,
    top1_methods = top1_methods_out,
    meta_reviewer_decisions = meta_reviewer_decisions
  )
}

resolve_lock_from_top1_clids <- function(top1_clids, cl_graph, min_votes) {
  lock_id <- NA_character_
  lock_reason <- "no_clids"

  top1_clids <- unique(as.character(top1_clids %||% character(0)))
  top1_clids <- top1_clids[!is.na(top1_clids) & nzchar(top1_clids)]

  if (length(top1_clids) == 1L) {
    # One valid ontology term can define context. Direct majority remains a separate soft anchor.
    lock_id <- top1_clids[[1]]
    lock_reason <- "single_clid"
  } else if (length(top1_clids) > 1L) {
    anc_maps <- lapply(top1_clids, function(clid) get_anc_map(clid, cl_graph))
    votes <- list()
    dist_sum <- list()
    dist_cnt <- list()
    for (am in anc_maps) {
      if (is.null(am) || length(am) == 0) next
      for (id in names(am)) {
        if (is.null(votes[[id]])) votes[[id]] <- 0
        votes[[id]] <- votes[[id]] + 1
        if (is.null(dist_sum[[id]])) dist_sum[[id]] <- 0
        dist_sum[[id]] <- dist_sum[[id]] + as.numeric(am[[id]])
        if (is.null(dist_cnt[[id]])) dist_cnt[[id]] <- 0
        dist_cnt[[id]] <- dist_cnt[[id]] + 1
      }
    }
    if (length(votes) > 0) {
      ids <- names(votes)
      vcnt <- vapply(ids, function(id) votes[[id]], numeric(1))
      # MCA fallback requires support from at least two distinct top1 ontology terms.
      effective_min_votes <- max(2L, as.integer(ceiling(min_votes %||% 1)))
      keep <- ids[vcnt >= effective_min_votes]
      if (length(keep) > 0) {
        mean_dist <- vapply(keep, function(id) {
          ds <- dist_sum[[id]] %||% 1
          dc <- dist_cnt[[id]] %||% 1
          ds / dc
        }, numeric(1))
        depth <- vapply(keep, function(id) {
          d <- suppressWarnings(as.numeric(get_depth_to_root(id, cl_graph)))
          if (is.finite(d)) d else -Inf
        }, numeric(1))
        keep_df <- data.frame(id=keep, votes=vcnt[keep], mean_dist=mean_dist, depth=depth, stringsAsFactors=FALSE)
        # Nearest supported ontology context first; raw ancestor vote count is only a tie-breaker.
        # This avoids the structural bias that makes broad/root ancestors win by construction.
        keep_df <- keep_df[order(keep_df$mean_dist, -keep_df$votes, -keep_df$depth, keep_df$id), , drop=FALSE]
        lock_id <- keep_df$id[[1]]
        lock_reason <- "consensus_msca"
      } else {
        lock_reason <- "insufficient_shared_ontology_support"
      }
    }
  }
  list(lock_id=lock_id, lock_reason=lock_reason)
}

compute_lock_depth_stats <- function(lock_id, cl_graph) {
  lock_depth <- NA_real_
  lock_anc_count <- NA_real_
  if (!is.na(lock_id) && nzchar(lock_id)) {
    lock_anc <- get_anc_map(lock_id, cl_graph)
    if (!is.null(lock_anc)) {
      lock_depth <- max(as.numeric(lock_anc), na.rm = TRUE)
      lock_anc_count <- length(lock_anc)
    }
  }
  list(lock_depth = lock_depth, lock_anc_count = lock_anc_count)
}

is_valid_clid_for_rules <- function(x) {
  is.character(x) && length(x) > 0 && !is.na(x[[1]]) && nzchar(x[[1]]) && startsWith(x[[1]], "CL:")
}

build_reviewer_eligibility <- function(top1_for_rules, cl_graph, min_depth_for_eligibility) {
  reviewer_diagnostics <- list()
  eligible <- list()
  eligible_methods_for_guard <- list()
  for (cand in top1_for_rules) {
    raw <- as.character(cand$raw_clid %||% NA_character_)
    coerced <- as.character(cand$clid %||% NA_character_)
    clid_final <- if (is_valid_clid_for_rules(raw)) raw else if (is_valid_clid_for_rules(coerced)) coerced else NA_character_
    depth <- if (is_valid_clid_for_rules(clid_final)) get_depth_to_root(clid_final, cl_graph) else NA_real_
    eligible_for_gate <- is.finite(depth) && depth >= min_depth_for_eligibility
    exclude_reason <- ""
    if (!eligible_for_gate) {
      if (!is_valid_clid_for_rules(clid_final)) exclude_reason <- "missing_or_invalid_clid"
      else exclude_reason <- paste0("depth_below_threshold:", depth, "<", min_depth_for_eligibility)
    }
    reviewer_diagnostics[[length(reviewer_diagnostics) + 1]] <- list(
      reviewer_key = cand$method,
      reviewer_name = ifelse(cand$method == "our", "in-house", cand$method),
      label = cand$label %||% "",
      raw_clid = raw,
      coerced_clid = coerced,
      clid_final_for_rules = clid_final,
      depth = depth,
      method_score = as.numeric(cand$score %||% NA_real_),
      eligible_for_gate = isTRUE(eligible_for_gate),
      exclude_reason = exclude_reason
    )
    if (isTRUE(eligible_for_gate)) {
      eligible[[length(eligible) + 1]] <- list(method = cand$method, clid = clid_final, depth = depth)
      eligible_methods_for_guard[[length(eligible_methods_for_guard) + 1]] <- list(
        method = cand$method,
        label = cand$label %||% "",
        clid = clid_final,
        depth = depth,
        score = as.numeric(cand$score %||% NA_real_)
      )
    }
  }
  list(
    reviewer_diagnostics = reviewer_diagnostics,
    eligible = eligible,
    eligible_methods_for_guard = eligible_methods_for_guard
  )
}

pick_best_method_for_guard <- function(pool_methods) {
  valid <- Filter(function(x) is.finite(x$score), pool_methods)
  if (length(valid) == 0) return(NULL)
  ord <- order(
    vapply(valid, function(x) x$score, numeric(1)),
    vapply(valid, function(x) x$depth %||% -Inf, numeric(1)),
    vapply(valid, function(x) as.character(x$method %||% ""), character(1)),
    decreasing = TRUE
  )
  valid[[ord[[1]]]]
}

compute_supported_lca_stats <- function(eligible, gate_k, cl_graph) {
  eligible_n <- length(eligible)
  supported_lca_clid <- NA_character_
  supported_lca_depth <- NA_real_
  support_n <- 0L
  support_gate_pass <- FALSE
  if (eligible_n >= gate_k) {
    best_pair_depth <- -Inf
    for (i in seq_len(eligible_n - 1)) {
      for (j in (i + 1):eligible_n) {
        lca <- get_msca_two(eligible[[i]]$clid, eligible[[j]]$clid, cl_graph)
        d <- as.numeric(lca$depth %||% NA_real_)
        if (!is.na(lca$id) && nzchar(lca$id) && is.finite(d) && d > best_pair_depth) {
          best_pair_depth <- d
          supported_lca_clid <- lca$id
          supported_lca_depth <- d
        }
      }
    }
    if (!is.na(supported_lca_clid) && nzchar(supported_lca_clid)) {
      support_n <- sum(vapply(eligible, function(x) isTRUE(is_descendant_of(x$clid, supported_lca_clid, cl_graph)) || identical(x$clid, supported_lca_clid), logical(1)))
      support_gate_pass <- (support_n >= gate_k)
    }
  }
  list(
    eligible_n = eligible_n,
    supported_lca_clid = supported_lca_clid,
    supported_lca_depth = supported_lca_depth,
    support_n = support_n,
    support_gate_pass = support_gate_pass
  )
}

select_mode_by_depth_gap <- function(lock_reason, support_gate_pass, depth_gap, delta_depth, disable_aggressive_trigger) {
  mode_selected <- "balanced"
  mode_trigger_reason <- "Balanced mode: trigger conditions not met."
  aggressive_trigger <- !isTRUE(disable_aggressive_trigger) && identical(lock_reason, "consensus_msca") && isTRUE(support_gate_pass) && is.finite(depth_gap) && depth_gap >= delta_depth
  if (isTRUE(aggressive_trigger)) {
    mode_selected <- "aggressive"
    mode_trigger_reason <- paste0("Aggressive mode enabled: consensus_msca lock depth gap ", sprintf("%.2f", depth_gap), " >= ", delta_depth, ".")
  }
  list(mode_selected = mode_selected, mode_trigger_reason = mode_trigger_reason, aggressive_trigger = aggressive_trigger)
}

pick_best_supported_clid_from_eligible <- function(eligible) {
  if (length(eligible) == 0) return(NA_character_)
  counts <- list()
  depths <- list()
  for (e in eligible) {
    k <- e$clid
    counts[[k]] <- (counts[[k]] %||% 0) + 1
    depths[[k]] <- max(depths[[k]] %||% 0, e$depth %||% 0)
  }
  ids <- names(counts)
  if (length(ids) == 0) return(NA_character_)
  ord <- order(vapply(ids, function(k) counts[[k]] %||% 0, numeric(1)),
               vapply(ids, function(k) depths[[k]] %||% 0, numeric(1)),
               ids,
               decreasing = TRUE)
  ids[[ord[[1]]]]
}

in_lock_scope <- function(clid, lock_id, cl_graph) {
  if (is.na(lock_id) || !nzchar(lock_id)) return(TRUE)
  if (is.na(clid) || !nzchar(clid)) return(FALSE)
  isTRUE(is_ancestor_of(lock_id, clid, cl_graph)) || isTRUE(clid == lock_id)
}

compute_anchor_clid_from_candidates <- function(cands, lock_id, cl_graph) {
  topk <- Filter(function(x) {
    !is.na(x$clid) && nzchar(x$clid) && x$rank <= 3 && x$method %in% c("cassia", "our", "enrich")
  }, cands)
  if (length(topk) == 0) return(lock_id)
  anchor <- topk[[1]]$clid
  if (length(topk) > 1) {
    for (i in 2:length(topk)) {
      anchor <- get_msca_two(anchor, topk[[i]]$clid, cl_graph)$id
    }
  }
  if (is.na(anchor) || !nzchar(anchor)) lock_id else anchor
}

order_candidates_by_priority <- function(cands, cfg_override = NULL) {
  if (length(cands) == 0) return(cands)
  pool_ctx <- build_pool_ctx(cands)
  cfg_local <- cfg_override %||% list(min_margin = 0.12, score_na_policy = "effective_surrogate")
  eff_scores <- vapply(cands, function(x) {
    as.numeric(effective_score(x, pool_ctx, cfg_local)$value %||% -Inf)
  }, numeric(1))
  ord <- order(
    vapply(cands, function(x) x$method_count %||% 0, numeric(1)),
    vapply(cands, function(x) x$map_quality_rank %||% 0, numeric(1)),
    eff_scores,
    vapply(cands, function(x) x$specificity, numeric(1)),
    decreasing = TRUE
  )
  cands[ord]
}

evaluate_branch_consistency_gate <- function(candidate_clid, anchor_clid, cl_graph, anchor_support, anchor_strong_threshold, branch_gate_slack = 1) {
  if (is.na(anchor_clid) || !nzchar(anchor_clid)) return(TRUE)
  if (is.na(candidate_clid) || !nzchar(candidate_clid)) return(FALSE)
  lca <- get_msca_two(candidate_clid, anchor_clid, cl_graph)
  if (is.null(lca$id) || is.na(lca$id)) return(FALSE)
  anchor_depth <- get_depth_to_root(anchor_clid, cl_graph)
  lca_depth <- lca$depth %||% 0
  lca_ok <- lca_depth >= (anchor_depth - branch_gate_slack)
  if (!lca_ok) return(FALSE)

  anchor_strength <- anchor_support[[anchor_clid]] %||% 0
  evidence_weak <- anchor_strength < (anchor_strong_threshold %||% Inf)
  if (evidence_weak) {
    candidate_depth <- get_depth_to_root(candidate_clid, cl_graph)
    depth_gap <- candidate_depth - anchor_depth
    depth_gap_cap <- 1
    if (depth_gap > depth_gap_cap) return(FALSE)
  }
  TRUE
}

evaluate_candidate_support <- function(candidate, pool_all, all_candidates, lock_id, min_margin, strong_score, cl_graph) {
  reasons <- character(0)
  review_flags <- character(0)
  scores <- vapply(pool_all, function(x) x$score_adj %||% NA_real_, numeric(1))
  scores <- scores[!is.na(scores)]
  top1 <- if (length(scores) > 0) max(scores) else NA_real_
  top2 <- if (length(scores) > 1) sort(scores, decreasing = TRUE)[[2]] else NA_real_
  margin <- if (length(scores) >= 2) (top1 - top2) else NA_real_
  if (!is.na(candidate$score_adj) && !is.na(candidate$score) && candidate$score < strong_score) review_flags <- c(review_flags, "weak_score")
  if (!is.na(margin) && margin < min_margin) review_flags <- c(review_flags, "small_margin")
  cross_conflict <- FALSE
  if (!is.na(lock_id) && nzchar(lock_id) && length(scores) > 0) {
    outside <- Filter(function(x) !in_lock_scope(x$clid, lock_id, cl_graph), all_candidates)
    out_scores <- vapply(outside, function(x) x$score %||% NA_real_, numeric(1))
    out_scores <- out_scores[!is.na(out_scores)]
    if (length(out_scores) > 0 && is.finite(top1) && max(out_scores) >= (top1 - min_margin)) {
      cross_conflict <- TRUE
    }
  }
  if (cross_conflict) review_flags <- c(review_flags, "cross_lineage_conflict")

  pass <- FALSE
  if (!is.na(candidate$method_count) && candidate$method_count >= 0.5) pass <- TRUE
  if (!pass && !is.na(candidate$map_quality) && candidate$map_quality >= 2 && !is.na(candidate$score)) pass <- TRUE
  if (!pass && isTRUE(candidate$top_ranked %||% FALSE)) pass <- TRUE

  list(pass = pass, reasons = reasons, review_flags = review_flags, margin = margin, cross_conflict = cross_conflict)
}

validate_policy_cfg <- function(policy_cfg) {
  allowed_top <- c("trigger", "override", "pool")
  if (!is.list(policy_cfg)) stop("policy_cfg must be a list")
  top_keys <- names(policy_cfg)
  if (is.null(top_keys) || length(top_keys) != 3L || !setequal(top_keys, allowed_top)) {
    stop("policy_cfg must contain exactly trigger, override, pool")
  }

  allowed_trigger <- c("score_margin_unsteady_max", "method_gap_unsteady_max", "evidence_gap_unsteady_max")
  allowed_override <- c("min_margin", "score_na_policy")
  allowed_pool <- c("soft_map_quality_min", "shadow_per_lineage")

  check_keys <- function(x, allowed, prefix) {
    if (!is.list(x)) stop(prefix, " must be a list")
    keys <- names(x) %||% character(0)
    unknown <- setdiff(keys, allowed)
    if (length(unknown) > 0) stop(prefix, " unknown keys: ", paste(unknown, collapse = ","))
  }

  check_nonneg_num <- function(x, key) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
      stop("policy_cfg invalid numeric non-negative value for ", key)
    }
  }

  check_keys(policy_cfg$trigger, allowed_trigger, "policy_cfg$trigger")
  check_keys(policy_cfg$override, allowed_override, "policy_cfg$override")
  check_keys(policy_cfg$pool, allowed_pool, "policy_cfg$pool")

  check_nonneg_num(policy_cfg$trigger$score_margin_unsteady_max, "trigger$score_margin_unsteady_max")
  check_nonneg_num(policy_cfg$trigger$method_gap_unsteady_max, "trigger$method_gap_unsteady_max")
  check_nonneg_num(policy_cfg$trigger$evidence_gap_unsteady_max, "trigger$evidence_gap_unsteady_max")
  check_nonneg_num(policy_cfg$override$min_margin, "override$min_margin")
  check_nonneg_num(policy_cfg$pool$soft_map_quality_min, "pool$soft_map_quality_min")
  check_nonneg_num(policy_cfg$pool$shadow_per_lineage, "pool$shadow_per_lineage")
  if (!is.character(policy_cfg$override$score_na_policy) || length(policy_cfg$override$score_na_policy) != 1L ||
      !(policy_cfg$override$score_na_policy %in% c("require_finite_new", "effective_surrogate"))) {
    stop("policy_cfg$override$score_na_policy must be one of: require_finite_new,effective_surrogate")
  }
  invisible(TRUE)
}

build_pool_ctx <- function(effective_pool) {
  pool <- effective_pool %||% list()
  raw <- vapply(pool, function(x) {
    xm <- suppressWarnings(as.numeric(x$method_count %||% 0))
    xe <- suppressWarnings(as.numeric(x$evidence_support %||% 0))
    xm <- if (is.finite(xm)) xm else 0
    xe <- if (is.finite(xe)) xe else 0
    max(0, xm) + max(0, xe)
  }, numeric(1))
  raw <- raw[is.finite(raw)]
  percentile_fn <- function(v) {
    if (!is.finite(v)) return(NA_real_)
    n <- length(raw)
    if (n == 0) return(0.5)
    lt <- sum(raw < v)
    eq <- sum(raw == v)
    (lt + 0.5 * eq) / n
  }
  list(raw = raw, percentile_fn = percentile_fn, pool = pool)
}

effective_score <- function(x, pool_ctx, cfg_override) {
  s <- suppressWarnings(as.numeric(x$score_adj %||% NA_real_))
  if (is.finite(s)) {
    return(list(value = s, source = "adj", raw_support = NA_real_))
  }
  xm <- suppressWarnings(as.numeric(x$method_count %||% 0))
  xe <- suppressWarnings(as.numeric(x$evidence_support %||% 0))
  xm <- if (is.finite(xm)) xm else 0
  xe <- if (is.finite(xe)) xe else 0
  raw_support <- max(0, xm) + max(0, xe)
  list(value = as.numeric(pool_ctx$percentile_fn(raw_support)), source = "surrogate", raw_support = raw_support)
}

compute_lock_instability <- function(cand_table, cfg_trigger) {

  if (length(cand_table) == 0) {
    return(list(
      unsteady = TRUE,
      reasons = c("no_candidate"),
      margin = NA_real_,
      method_gap = NA_real_,
      ev_gap = NA_real_
    ))
  }

  pool_ctx <- build_pool_ctx(cand_table)
  cfg_override_local <- list(min_margin = cfg_trigger$score_margin_unsteady_max, score_na_policy = "effective_surrogate")
  scores <- vapply(cand_table, function(x) {
    effective_score(x, pool_ctx, cfg_override_local)$value
  }, numeric(1))
  ord <- order(scores, decreasing = TRUE, na.last = TRUE)

  top1 <- cand_table[[ord[[1]]]]
  top2 <- if (length(ord) >= 2) cand_table[[ord[[2]]]] else NULL

  s1 <- as.numeric(effective_score(top1, pool_ctx, cfg_override_local)$value %||% NA_real_)
  s2 <- if (is.null(top2)) NA_real_ else as.numeric(effective_score(top2, pool_ctx, cfg_override_local)$value %||% NA_real_)
  m1 <- as.numeric(top1$method_count %||% 0)
  m2 <- if (is.null(top2)) 0 else as.numeric(top2$method_count %||% 0)
  e1 <- as.numeric(top1$evidence_support %||% 0)
  e2 <- if (is.null(top2)) 0 else as.numeric(top2$evidence_support %||% 0)

  margin <- if (is.finite(s1) && is.finite(s2)) (s1 - s2) else NA_real_
  method_gap <- m1 - m2
  ev_gap <- e1 - e2

  reasons <- character(0)
  if (!is.finite(margin) || margin <= cfg_trigger$score_margin_unsteady_max) reasons <- c(reasons, "small_margin")
  if (!is.finite(method_gap) || method_gap <= cfg_trigger$method_gap_unsteady_max) reasons <- c(reasons, "small_method_gap")
  if (!is.finite(ev_gap) || ev_gap <= cfg_trigger$evidence_gap_unsteady_max) reasons <- c(reasons, "small_evidence_gap")

  list(
    unsteady = length(reasons) > 0,
    reasons = unique(reasons),
    margin = margin,
    method_gap = method_gap,
    ev_gap = ev_gap
  )
}

should_enter_rescue <- function(chosen_generic, instability) {
  isTRUE(chosen_generic) && isTRUE(instability$unsteady)
}

should_enter_non_generic_rescue <- function(effective_pool, base, cfg, cl_graph) {
  if (is.null(effective_pool) || length(effective_pool) == 0) {
    return(list(ok = FALSE, reason = "empty_pool", checks = list()))
  }
  if (is.null(base) || is.null(cl_graph) || !nzchar(as.character(base$clid %||% ""))) {
    return(list(ok = FALSE, reason = "missing_base_or_graph", checks = list()))
  }

  pool_ctx <- build_pool_ctx(effective_pool)
  cfg_override <- cfg$override %||% list(min_margin = 0.12, score_na_policy = "effective_surrogate")
  ord <- order(
    vapply(effective_pool, function(x) as.numeric(effective_score(x, pool_ctx, cfg_override)$value %||% -Inf), numeric(1)),
    vapply(effective_pool, function(x) as.numeric(x$method_count %||% 0), numeric(1)),
    vapply(effective_pool, function(x) as.numeric(x$evidence_support %||% 0), numeric(1)),
    decreasing = TRUE
  )
  best <- effective_pool[[ord[[1]]]]
  ov <- should_override(base, best, cfg_override, pool_ctx = pool_ctx, cl_graph = cl_graph)

  is_desc <- isTRUE(ov$checks$candidate_is_descendant %||% FALSE)
  consistency_ok <- isTRUE(ov$checks$consistency_ok %||% FALSE)
  stricter_desc_margin_ok <- isTRUE(ov$checks$stricter_desc_margin_ok %||% FALSE)
  ok <- isTRUE(ov$ok) && is_desc && consistency_ok && stricter_desc_margin_ok

  list(
    ok = ok,
    reason = if (ok) "pass" else "best_candidate_failed_strict_descendant_gate",
    checks = ov$checks,
    best_candidate = list(
      clid = as.character(best$clid %||% ""),
      score_adj = as.numeric(best$score_adj %||% NA_real_),
      method_count = as.numeric(best$method_count %||% 0),
      evidence_support = as.numeric(best$evidence_support %||% 0)
    )
  )
}

build_candidate_pools <- function(cand_ordered,
                                  gate_decisions,
                                  chosen,
                                  cl_graph,
                                  cfg_pool,
                                  require_deeper = FALSE,
                                  is_generic_fn,
                                  is_stage_fn) {
  chosen_clid <- chosen$clid %||% NA_character_
  chosen_depth <- get_depth_to_root(chosen_clid, cl_graph)

  check_desc_eligibility <- function(x, stage_mode = "strict") {
    clid <- x$clid %||% NA_character_
    label <- x$label %||% ""
    if (is.na(clid) || !nzchar(clid)) return(list(ok = FALSE, reason = "missing_clid", stage_flag = FALSE))
    if (!is_descendant_of(clid, chosen_clid, cl_graph)) return(list(ok = FALSE, reason = "not_descendant", stage_flag = FALSE))
    if (isTRUE(require_deeper)) {
      xd <- get_depth_to_root(clid, cl_graph)
      if (!is.finite(xd) || !is.finite(chosen_depth) || xd <= chosen_depth) return(list(ok = FALSE, reason = "require_deeper_fail", stage_flag = FALSE))
    }
    if (isTRUE(is_generic_fn(clid, label))) return(list(ok = FALSE, reason = "generic_filtered", stage_flag = FALSE))
    stage_flag <- isTRUE(is_stage_fn(clid, label))
    # In shadow mode: keep stage-flagged candidates but mark them; prefer non-stage if available
    # In strict mode: filter out stage-flagged candidates
    if (isTRUE(stage_flag) && identical(stage_mode, "strict")) {
      return(list(ok = FALSE, reason = "stage_filtered", stage_flag = TRUE))
    }
    list(ok = TRUE, reason = "pass", stage_flag = stage_flag)
  }

  is_valid_desc <- function(x) {
    isTRUE(check_desc_eligibility(x, stage_mode = "strict")$ok)
  }

  rank_pool <- function(pool) {
    if (length(pool) == 0) return(pool)
    pool_ctx <- build_pool_ctx(pool)
    cfg_local <- list(min_margin = 0.12, score_na_policy = "effective_surrogate")
    ord <- order(
      vapply(pool, function(x) as.numeric(effective_score(x, pool_ctx, cfg_local)$value %||% -Inf), numeric(1)),
      vapply(pool, function(x) x$method_count %||% 0, numeric(1)),
      vapply(pool, function(x) x$evidence_support %||% 0, numeric(1)),
      vapply(pool, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
      decreasing = TRUE
    )
    pool[ord]
  }

  candidate_from_decision <- function(x) x$candidate %||% list()

  gate_pool <- Filter(function(x) is_valid_desc(x), lapply(Filter(function(d) isTRUE(d$kept), gate_decisions), candidate_from_decision))
  soft_pool <- Filter(function(x) {
    if (!is_valid_desc(x)) return(FALSE)
    as.numeric(x$map_quality %||% 0) >= cfg_pool$soft_map_quality_min
  }, cand_ordered)

  lineage_key <- function(clid) {
    anc <- get_anc_map(clid, cl_graph)
    if (is.null(anc) || length(anc) == 0) return(clid)
    if (!is.finite(chosen_depth)) return(clid)
    anc_ids <- names(anc)
    anc_depth <- as.numeric(anc)
    idx <- which(anc_depth == (chosen_depth + 1))
    if (length(idx) > 0) return(as.character(anc_ids[idx[[1]]]))
    clid
  }

  rejected_all <- Filter(function(d) !isTRUE(d$kept), gate_decisions)
  rejected_desc <- Filter(function(d) {
    cand <- d$candidate %||% list()
    clid <- cand$clid %||% NA_character_
    if (is.na(clid) || !nzchar(clid)) return(FALSE)
    isTRUE(is_descendant_of(clid, chosen_clid, cl_graph))
  }, rejected_all)

  shadow_desc_filter_counts <- list(
    missing_clid = 0L,
    not_descendant = 0L,
    require_deeper_fail = 0L,
    generic_filtered = 0L,
    pass = 0L
  )
  shadow_attempted <- lapply(rejected_desc, candidate_from_decision)
  shadow_raw <- list()
  shadow_stage_flags <- logical(0)
  shadow_stage_flagged_n <- 0L
  for (cand in shadow_attempted) {
    chk <- check_desc_eligibility(cand, stage_mode = "shadow")
    key <- as.character(chk$reason %||% "not_descendant")
    # Only count in filter_counts if ok=FALSE (these are actual filter reasons)
    if (!isTRUE(chk$ok)) {
      if (is.null(shadow_desc_filter_counts[[key]])) shadow_desc_filter_counts[[key]] <- 0L
      shadow_desc_filter_counts[[key]] <- as.integer(shadow_desc_filter_counts[[key]] %||% 0L) + 1L
    } else {
      # Candidate passed eligibility check
      shadow_desc_filter_counts$pass <- as.integer(shadow_desc_filter_counts$pass %||% 0L) + 1L
      if (isTRUE(chk$stage_flag)) {
        shadow_stage_flagged_n <- shadow_stage_flagged_n + 1L
      }
      shadow_raw[[length(shadow_raw) + 1L]] <- cand
      shadow_stage_flags <- c(shadow_stage_flags, isTRUE(chk$stage_flag))
    }
  }
  if (length(shadow_raw) > 0) {
    shadow_ctx <- build_pool_ctx(shadow_raw)
    cfg_local <- list(min_margin = 0.12, score_na_policy = "effective_surrogate")
    ord_shadow <- order(
      vapply(shadow_raw, function(x) as.numeric(effective_score(x, shadow_ctx, cfg_local)$value %||% -Inf), numeric(1)),
      vapply(shadow_raw, function(x) x$method_count %||% 0, numeric(1)),
      vapply(shadow_raw, function(x) x$evidence_support %||% 0, numeric(1)),
      vapply(shadow_raw, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
      decreasing = TRUE
    )
    shadow_raw <- shadow_raw[ord_shadow]
    shadow_stage_flags <- shadow_stage_flags[ord_shadow]
  }
  seen_key <- character(0)
  per_lineage_cap <- max(1L, as.integer(cfg_pool$shadow_per_lineage))
  key_count <- list()
  shadow_pool <- list()
  shadow_non_stage <- list()
  shadow_stage_only <- list()
  if (length(shadow_raw) > 0) {
    for (i in seq_along(shadow_raw)) {
      cand <- shadow_raw[[i]]
      stage_flag <- if (length(shadow_stage_flags) >= i) isTRUE(shadow_stage_flags[[i]]) else FALSE
      if (isTRUE(stage_flag)) {
        shadow_stage_only[[length(shadow_stage_only) + 1L]] <- cand
      } else {
        shadow_non_stage[[length(shadow_non_stage) + 1L]] <- cand
      }
    }
  }

  shadow_ranked_for_pick <- if (length(shadow_non_stage) > 0) shadow_non_stage else shadow_stage_only
  if (length(shadow_non_stage) == 0 && length(shadow_stage_only) > 0) {
    shadow_ranked_for_pick <- list(shadow_stage_only[[1]])
  }

  for (cand in shadow_ranked_for_pick) {
    k <- lineage_key(cand$clid %||% "")
    n <- key_count[[k]] %||% 0L
    if (n >= per_lineage_cap) next
    shadow_pool[[length(shadow_pool) + 1L]] <- cand
    key_count[[k]] <- n + 1L
    seen_key <- c(seen_key, k)
  }

  reject_stage_counts <- table(vapply(rejected_desc, function(d) as.character(d$reject_stage %||% "unknown"), character(1)))
  reject_stage_counts <- as.list(as.integer(reject_stage_counts))
  names(reject_stage_counts) <- names(table(vapply(rejected_desc, function(d) as.character(d$reject_stage %||% "unknown"), character(1))))

  shadow_empty_reason <- "not_empty"
  if (length(shadow_pool) == 0) {
    if (length(rejected_all) == 0) {
      shadow_empty_reason <- "no_rejects"
    } else if (length(rejected_desc) == 0) {
      shadow_empty_reason <- "no_descendant_rejects"
    } else {
      shadow_empty_reason <- "descendant_filtered_out"
    }
  }

  list(
    gate_pool = rank_pool(gate_pool),
    soft_pool = rank_pool(soft_pool),
    shadow_pool = shadow_pool,
    meta = list(
      gate_pool_n = length(gate_pool),
      soft_pool_n = length(soft_pool),
      shadow_pool_n = length(shadow_pool),
      gate_reject_total_n = length(rejected_all),
      gate_reject_descendant_n = length(rejected_desc),
      gate_reject_descendant_by_stage = reject_stage_counts,
      shadow_kept_n = length(shadow_pool),
      shadow_empty_reason = shadow_empty_reason,
      shadow_desc_filter_counts = shadow_desc_filter_counts,
      shadow_stage_flagged_n = shadow_stage_flagged_n,
      shadow_desc_ref = list(attempted_n = length(shadow_attempted), kept_n = length(shadow_pool))
    )
  )
}

compute_base_structure_guard <- function(base, candidate, cl_graph, cfg_override) {
  base_clid <- as.character(base$clid %||% "")
  cand_clid <- as.character(candidate$clid %||% "")
  base_depth <- NA_real_
  candidate_depth <- NA_real_
  base_subtree_size <- 0L
  base_subtree_ratio <- NA_real_
  candidate_is_descendant <- FALSE
  base_is_ultra_generic <- FALSE

  if (!is.null(cl_graph) && nzchar(base_clid)) {
    base_depth <- as.numeric(get_depth_to_root(base_clid, cl_graph))
    max_desc_depth <- as.integer(cfg_override$ultra_generic_desc_depth_k %||% 64L)
    desc_df <- get_descendants(base_clid, max_desc_depth, cl_graph)
    base_subtree_size <- as.integer(nrow(desc_df))
    total_nodes <- as.integer(length(cl_graph$cl %||% list()))
    if (is.finite(total_nodes) && total_nodes > 0) {
      base_subtree_ratio <- as.numeric(base_subtree_size / total_nodes)
    }
    if (nzchar(cand_clid)) {
      candidate_is_descendant <- isTRUE(is_descendant_of(cand_clid, base_clid, cl_graph))
      candidate_depth <- as.numeric(get_depth_to_root(cand_clid, cl_graph))
    }
  }

  ultra_depth_max <- as.numeric(cfg_override$ultra_generic_depth_max %||% 2)
  ultra_subtree_ratio_min <- as.numeric(cfg_override$ultra_generic_subtree_ratio_min %||% 0.02)
  ultra_subtree_size_min <- as.integer(cfg_override$ultra_generic_subtree_size_min %||% 50L)

  base_is_ultra_generic <- is.finite(base_depth) && (base_depth <= ultra_depth_max) &&
    ((is.finite(base_subtree_ratio) && base_subtree_ratio >= ultra_subtree_ratio_min) ||
      (is.finite(base_subtree_size) && base_subtree_size >= ultra_subtree_size_min))

  list(
    base_depth = as.numeric(base_depth),
    candidate_depth = as.numeric(candidate_depth),
    base_subtree_size = as.integer(base_subtree_size),
    base_subtree_ratio = as.numeric(base_subtree_ratio),
    base_is_ultra_generic = isTRUE(base_is_ultra_generic),
    candidate_is_descendant = isTRUE(candidate_is_descendant)
  )
}

should_override <- function(base, candidate, cfg_override, pool_ctx = NULL, cl_graph = NULL) {
  b_method <- as.numeric(base$method_count %||% 0)
  b_ev <- as.numeric(base$evidence_support %||% 0)
  b_score <- as.numeric(base$score_adj %||% NA_real_)

  c_method <- as.numeric(candidate$method_count %||% 0)
  c_ev <- as.numeric(candidate$evidence_support %||% 0)
  c_score <- as.numeric(candidate$score_adj %||% NA_real_)

  support_up <- (c_method > b_method) || (c_ev > b_ev)
  support_both_up <- (c_method > b_method) && (c_ev > b_ev)

  score_na_case <- "none"
  if (!is.finite(b_score) && !is.finite(c_score)) {
    score_na_case <- "both_na"
  } else if (!is.finite(b_score) && is.finite(c_score)) {
    score_na_case <- "base_na"
  } else if (is.finite(b_score) && !is.finite(c_score)) {
    score_na_case <- "new_na"
  }

  score_ok <- FALSE
  score_ok_reason <- ""
  if (is.null(pool_ctx)) {
    pool_ctx <- build_pool_ctx(list(base, candidate))
  }
  es_base <- effective_score(base, pool_ctx, cfg_override)
  es_new <- effective_score(candidate, pool_ctx, cfg_override)
  score_source_base <- as.character(es_base$source %||% "surrogate")
  score_source_new <- as.character(es_new$source %||% "surrogate")
  effective_base_score <- as.numeric(es_base$value %||% NA_real_)
  effective_new_score <- as.numeric(es_new$value %||% NA_real_)

  struct_guard <- compute_base_structure_guard(base, candidate, cl_graph, cfg_override)
  candidate_is_descendant <- isTRUE(struct_guard$candidate_is_descendant)
  base_is_ultra_generic <- isTRUE(struct_guard$base_is_ultra_generic)
  support_not_down_both <- !((c_method < b_method) && (c_ev < b_ev))

  consistency_ok <- TRUE
  consistency_ok_reason <- "non_descendant"
  if (isTRUE(candidate_is_descendant) && !isTRUE(base_is_ultra_generic)) {
    consistency_ok <- isTRUE(support_both_up)
    consistency_ok_reason <- if (isTRUE(consistency_ok)) {
      "descendant_requires_both_support_pass"
    } else {
      "descendant_requires_both_support_fail"
    }
  } else if (isTRUE(candidate_is_descendant) && isTRUE(base_is_ultra_generic)) {
    consistency_ok <- TRUE
    consistency_ok_reason <- "ultra_generic_descendant_permissive"
  }

  # Pool-size aware gating with stricter conditions for small pools
  pool_size <- length(pool_ctx$raw %||% pool_ctx$pool %||% list())
  pool_size_1_strong_margin_used <- FALSE
  tie_depth_gain_used <- FALSE
  support_gate_ok <- isTRUE(support_up)

  # For pool_size == 1: require descendant relation + stronger score margin
  # and allow support_up to be either method OR evidence under this strong margin.
  # For pool_size == 2: require stronger margin.
  if (pool_size == 1) {
    pool_size_1_strong_margin_used <- TRUE
    strong_margin <- as.numeric(cfg_override$pool_size_1_strong_margin %||% 0.35)
    if (!isTRUE(candidate_is_descendant)) {
      return(list(
        ok = FALSE,
        reasons = c("small_pool_requires_descendant"),
        checks = list(
          support_up = isTRUE(support_up),
          support_not_down_both = isTRUE(support_not_down_both),
          score_ok = FALSE,
          consistency_ok = FALSE,
          consistency_ok_reason = "pool_size_1_requires_descendant",
          score_na_case = score_na_case,
          score_ok_reason = "pool_size_1_requires_descendant",
          score_source_base = score_source_base,
          score_source_new = score_source_new,
          effective_base_score = as.numeric(effective_base_score),
          effective_new_score = as.numeric(effective_new_score),
          candidate_is_descendant = isTRUE(candidate_is_descendant),
          pool_size = pool_size,
          pool_size_1_strong_margin_used = TRUE,
          tie_depth_gain_used = FALSE,
          base_depth = struct_guard$base_depth,
          candidate_depth = struct_guard$candidate_depth,
          base_subtree_size = struct_guard$base_subtree_size,
          base_is_ultra_generic = struct_guard$base_is_ultra_generic,
          stricter_desc_margin_ok = TRUE
        )
      ))
    }

    strong_margin_ok <- is.finite(effective_new_score) && is.finite(effective_base_score) &&
      (effective_new_score >= (effective_base_score + strong_margin))
    if (!isTRUE(strong_margin_ok)) {
      depth_gain_min <- as.integer(cfg_override$depth_gain_min %||% 1L)
      tie_depth_gain_ok <- isTRUE(base_is_ultra_generic) && isTRUE(candidate_is_descendant) &&
        is.finite(effective_new_score) && is.finite(effective_base_score) &&
        (effective_new_score >= effective_base_score) &&
        is.finite(struct_guard$candidate_depth) && is.finite(struct_guard$base_depth) &&
        (struct_guard$candidate_depth >= (struct_guard$base_depth + depth_gain_min)) &&
        isTRUE(support_up)

      if (isTRUE(tie_depth_gain_ok)) {
        tie_depth_gain_used <- TRUE
        consistency_ok <- TRUE
        consistency_ok_reason <- "pool_size_1_tie_allow_depth_gain"
        support_gate_ok <- isTRUE(support_up)
      } else {
        return(list(
          ok = FALSE,
          reasons = c("small_pool_strong_margin_fail"),
          checks = list(
            support_up = isTRUE(support_up),
            support_not_down_both = isTRUE(support_not_down_both),
            score_ok = FALSE,
            consistency_ok = FALSE,
            consistency_ok_reason = "pool_size_1_strong_margin_fail",
            score_na_case = score_na_case,
            score_ok_reason = paste0("pool_size_1_requires_strong_margin_", strong_margin),
            score_source_base = score_source_base,
            score_source_new = score_source_new,
            effective_base_score = as.numeric(effective_base_score),
            effective_new_score = as.numeric(effective_new_score),
            candidate_is_descendant = isTRUE(candidate_is_descendant),
            pool_size = pool_size,
            pool_size_1_strong_margin_used = TRUE,
            tie_depth_gain_used = FALSE,
            base_depth = struct_guard$base_depth,
            candidate_depth = struct_guard$candidate_depth,
            base_subtree_size = struct_guard$base_subtree_size,
            base_is_ultra_generic = struct_guard$base_is_ultra_generic,
            stricter_desc_margin_ok = TRUE
          )
        ))
      }
    }

    # For pool_size==1 and strong margin pass, permit descendant rescue with support_up (OR) only.
    if (!isTRUE(tie_depth_gain_used)) {
      consistency_ok <- TRUE
      consistency_ok_reason <- "pool_size_1_strong_margin_override"
      support_gate_ok <- isTRUE(support_up)
    }
  } else if (pool_size == 2) {
    # Small pool: require stronger margin
    min_margin_small_pool <- as.numeric(cfg_override$min_margin_small_pool %||% 0.2)
    if (is.finite(effective_new_score) && is.finite(effective_base_score) &&
        (effective_new_score < (effective_base_score + min_margin_small_pool))) {
      return(list(
        ok = FALSE,
        reasons = c("small_pool_margin_fail"),
        checks = list(
          support_up = isTRUE(support_up),
          support_not_down_both = isTRUE(support_not_down_both),
          score_ok = FALSE,
          consistency_ok = isTRUE(consistency_ok),
          consistency_ok_reason = consistency_ok_reason,
          score_na_case = score_na_case,
          score_ok_reason = paste0("pool_size_2_requires_margin_", min_margin_small_pool),
          score_source_base = score_source_base,
          score_source_new = score_source_new,
          effective_base_score = as.numeric(effective_base_score),
          effective_new_score = as.numeric(effective_new_score),
          candidate_is_descendant = isTRUE(candidate_is_descendant),
          pool_size = pool_size,
          pool_size_1_strong_margin_used = FALSE,
          tie_depth_gain_used = FALSE,
          base_depth = struct_guard$base_depth,
          candidate_depth = struct_guard$candidate_depth,
          base_subtree_size = struct_guard$base_subtree_size,
          base_is_ultra_generic = struct_guard$base_is_ultra_generic,
          stricter_desc_margin_ok = TRUE
        )
      ))
    }
  }
  # pool_size >= 3: normal flow, no hard block

  if (identical(cfg_override$score_na_policy, "require_finite_new")) {
    if (!is.finite(c_score)) {
      score_ok <- FALSE
      score_ok_reason <- "new_non_finite_blocked"
    } else if (!is.finite(b_score)) {
      score_ok <- TRUE
      score_ok_reason <- "base_non_finite_new_finite_pass"
    } else {
      score_ok <- (c_score >= (b_score - cfg_override$min_margin))
      score_ok_reason <- if (isTRUE(score_ok)) "finite_margin_pass" else "finite_margin_fail"
    }
  } else if (identical(cfg_override$score_na_policy, "effective_surrogate")) {
    score_ok <- is.finite(effective_new_score) && is.finite(effective_base_score) &&
      (effective_new_score >= (effective_base_score - cfg_override$min_margin))
    score_ok_reason <- if (isTRUE(score_ok)) "effective_margin_pass" else "effective_margin_fail"
  }
  if (isTRUE(tie_depth_gain_used)) {
    tie_score_ok <- is.finite(effective_new_score) && is.finite(effective_base_score) &&
      (effective_new_score >= effective_base_score)
    score_ok <- isTRUE(tie_score_ok)
    score_ok_reason <- if (isTRUE(tie_score_ok)) {
      "pool_size_1_tie_allow_depth_gain"
    } else {
      "pool_size_1_tie_score_not_worse_fail"
    }
  }

  stricter_desc_margin_ok <- TRUE
  stricter_desc_margin <- as.numeric(cfg_override$descendant_strict_margin %||% 0.1)
  if (isTRUE(candidate_is_descendant) && !isTRUE(base_is_ultra_generic)) {
    stricter_desc_margin_ok <- is.finite(effective_new_score) && is.finite(effective_base_score) &&
      (effective_new_score >= (effective_base_score + stricter_desc_margin))
  }

  reasons <- character(0)
  if (!isTRUE(support_gate_ok)) reasons <- c(reasons, "support_failed")
  if (!isTRUE(score_ok)) reasons <- c(reasons, "score_failed")
  if (!isTRUE(consistency_ok)) reasons <- c(reasons, "consistency_failed")
  if (!isTRUE(stricter_desc_margin_ok)) reasons <- c(reasons, "descendant_margin_failed")
  if (isTRUE(candidate_is_descendant) && !isTRUE(base_is_ultra_generic) &&
      (!isTRUE(consistency_ok) || !isTRUE(stricter_desc_margin_ok))) {
    reasons <- c(reasons, "lock_confident_freeze")
  }

  final_ok <- isTRUE(support_gate_ok) && isTRUE(score_ok) && isTRUE(consistency_ok) && isTRUE(stricter_desc_margin_ok)
  list(
    ok = final_ok,
    reasons = reasons,
    checks = list(
      support_up = isTRUE(support_up),
      support_not_down_both = isTRUE(support_not_down_both),
      support_gate_ok = isTRUE(support_gate_ok),
      consistency_ok = isTRUE(consistency_ok),
      consistency_ok_reason = consistency_ok_reason,
      score_ok = isTRUE(score_ok),
      score_na_case = score_na_case,
      score_ok_reason = score_ok_reason,
      score_source_base = score_source_base,
      score_source_new = score_source_new,
      effective_base_score = as.numeric(effective_base_score),
      effective_new_score = as.numeric(effective_new_score),
      candidate_is_descendant = isTRUE(candidate_is_descendant),
      base_depth = as.numeric(struct_guard$base_depth),
      candidate_depth = as.numeric(struct_guard$candidate_depth),
      base_subtree_size = as.integer(struct_guard$base_subtree_size),
      base_is_ultra_generic = isTRUE(struct_guard$base_is_ultra_generic),
      stricter_desc_margin_ok = isTRUE(stricter_desc_margin_ok),
      pool_size = as.integer(pool_size),
      pool_size_1_strong_margin_used = isTRUE(pool_size_1_strong_margin_used),
      tie_depth_gain_used = isTRUE(tie_depth_gain_used)
    )
  )
}

apply_post_fallback <- function(state, pools, cfg, rule_tag, cl_graph = NULL) {
  out <- list(
    updated_state = state,
    logs = list(
      gate_pool_n = pools$meta$gate_pool_n %||% 0L,
      soft_pool_n = pools$meta$soft_pool_n %||% 0L,
      shadow_pool_n = pools$meta$shadow_pool_n %||% 0L,
      gate_reject_total_n = pools$meta$gate_reject_total_n %||% 0L,
      gate_reject_descendant_n = pools$meta$gate_reject_descendant_n %||% 0L,
      gate_reject_descendant_by_stage = pools$meta$gate_reject_descendant_by_stage %||% list(),
      shadow_kept_n = pools$meta$shadow_kept_n %||% 0L,
      shadow_empty_reason = pools$meta$shadow_empty_reason %||% "not_empty",
      shadow_desc_filter_counts = pools$meta$shadow_desc_filter_counts %||% list(),
      shadow_desc_ref = pools$meta$shadow_desc_ref %||% list(attempted_n = 0L, kept_n = 0L),
      shadow_pool_used = FALSE,
      why_not_override = "no_candidate",
      override_checks = list(support_up = FALSE, support_not_down_both = TRUE, support_gate_ok = FALSE, consistency_ok = TRUE, consistency_ok_reason = "", score_ok = FALSE, score_na_case = "none", score_ok_reason = "", score_source_base = "", score_source_new = "", effective_base_score = NA_real_, effective_new_score = NA_real_, candidate_is_descendant = FALSE, base_depth = NA_real_, candidate_depth = NA_real_, base_subtree_size = 0L, base_is_ultra_generic = FALSE, stricter_desc_margin_ok = TRUE, pool_size = 0L, pool_size_1_strong_margin_used = FALSE, tie_depth_gain_used = FALSE),
      override_applied = FALSE,
      override_reason = ""
    )
  )

  candidate <- NULL
  if (length(pools$gate_pool) > 0) {
    candidate <- pools$gate_pool[[1]]
  } else if (length(pools$soft_pool) > 0) {
    candidate <- pools$soft_pool[[1]]
  } else if (length(pools$shadow_pool) > 0) {
    out$logs$shadow_pool_used <- TRUE
    candidate <- pools$shadow_pool[[1]]
  } else {
    out$logs$why_not_override <- "shadow_pool_empty"
    return(out)
  }

  pool_ctx <- build_pool_ctx(c(pools$soft_pool %||% list(), pools$shadow_pool %||% list(), pools$gate_pool %||% list()))
  ov <- should_override(state$base, candidate, cfg$override, pool_ctx = pool_ctx, cl_graph = cl_graph)
  out$logs$override_checks <- ov$checks
  if (isTRUE(ov$ok)) {
    out$updated_state$chosen <- candidate
    out$updated_state$final_rule <- paste0("iter9_post_fallback_", rule_tag)
    out$updated_state$anchor_fallback_happened <- FALSE
    out$logs$override_applied <- TRUE
    out$logs$why_not_override <- "overridden"
    out$logs$override_reason <- rule_tag
    return(out)
  }

  if (isTRUE("lock_confident_freeze" %in% (ov$reasons %||% character(0)))) {
    out$logs$why_not_override <- "lock_confident_freeze"
  } else if (isTRUE(out$logs$shadow_pool_used)) {
    out$logs$why_not_override <- "shadow_pool_used_but_failed"
  } else if (!isTRUE(ov$checks$support_up)) {
    out$logs$why_not_override <- "support_failed"
  } else if (!isTRUE(ov$checks$consistency_ok)) {
    out$logs$why_not_override <- "consistency_failed"
  } else {
    out$logs$why_not_override <- "score_failed"
  }
  out
}

apply_consensus_subtype_policy <- function(head_out, judge_input_obj, cl_cfg, cl_graph,
                                           dataset_cfg = NULL,
                                           min_votes = 1.0,
                                           min_margin = 0.12,
                                           strong_score = 0.6,
                                           expand_depth_k = 3,
                                           score_lambda = 2.5,
                                           lock_beta = 1.5,
                                           conflict_delta = 3,
                                           specificity_eps = 0.05,
                                           branch_gate_slack = 1,
                                           min_depth_for_eligibility = 3,
                                           delta_depth = 2,
                                           gate_mode = "fixed_k",
                                            gate_k = 2,
                                            disable_aggressive_trigger = FALSE,
                                           enable_method_reliability_gate = FALSE,
                                            enable_aggressive_lca_hard_guard = TRUE,
                                             protect_lock_inputs = FALSE,
                                             meta_reviewer_gate = FALSE,
                                             parity_mode = FALSE,
                                             identity_policy = "ontology_guarded",
                                             species_value = "human",
                                            policy_cfg = NULL) {
  if (is.null(head_out) || !is.list(head_out)) return(head_out)
  if (is.null(judge_input_obj) || !is.list(judge_input_obj)) return(head_out)
  if (is.null(cl_graph) || is.null(cl_graph$cl)) return(head_out)
  identity_policy <- tolower(trimws(as.character(identity_policy %||% "ontology_guarded")))
  if (!identity_policy %in% c("ontology_guarded", "head_preserve", "legacy_rerank")) {
    stop("identity_policy must be one of: ontology_guarded, head_preserve, legacy_rerank")
  }
  inputs <- judge_input_obj$inputs %||% list()
  candidates <- collect_candidates_from_inputs(inputs, cl_cfg, cl_graph)
  candidates_raw <- candidates

  count_by_method <- function(cands) {
    if (length(cands) == 0) return(list())
    methods <- vapply(cands, function(x) as.character(x$method %||% "other"), character(1))
    methods[!nzchar(methods)] <- "other"
    tb <- table(methods)
    out <- as.list(as.integer(tb))
    names(out) <- names(tb)
    out
  }

  if (length(candidates) == 0) return(head_out)

  n_raw_candidates_by_method <- count_by_method(candidates_raw)

  tissue_filter_applied <- FALSE
  tissue_filter_reason <- ""
  tissue_filter_labels <- character(0)
  tissue_filter_n_before <- as.integer(length(candidates))
  tissue_filter_n_after <- as.integer(length(candidates))
  tissue_filter_results <- list()
  tissue_filter_summary <- list(n_in = 0L, n_out = 0L, n_uncertain = 0L, suspect_leak = FALSE, suspect_leak_reason = "")

  normalize_token <- function(x) {
    x <- tolower(stringr::str_squish(as.character(x %||% "")))
    x <- gsub("[^a-z0-9 ]", " ", x, perl = TRUE)
    stringr::str_squish(x)
  }
  normalize_verdict <- function(v) {
    x <- normalize_token(v)
    if (identical(x, "in tissue")) x <- "in_tissue"
    if (identical(x, "out of tissue")) x <- "out_of_tissue"
    if (!x %in% c("in_tissue", "out_of_tissue", "uncertain")) x <- "uncertain"
    x
  }
  extract_llm_tissue_payload <- function(h, tissue_prior) {
    if (is.null(h) || !is.list(h)) return(NULL)
    payload <- h$tissue_candidate_filter %||% h$tissue_filter %||% NULL
    if (is.null(payload) || !is.list(payload)) return(NULL)
    p_tissue <- as.character(payload$tissue_prior %||% "")
    if (nzchar(tissue_prior) && nzchar(p_tissue) && !identical(normalize_token(tissue_prior), normalize_token(p_tissue))) {
      return(NULL)
    }
    payload
  }

  cfg_tissue <- as.character((dataset_cfg %||% list())$tissue %||% character(0))
  cfg_tissue <- unique(cfg_tissue[!is.na(cfg_tissue) & nzchar(cfg_tissue)])
  tissue_prior <- as.character(cfg_tissue[[1]] %||% "")
  tissue_filter_labels <- cfg_tissue

  if (!nzchar(tissue_prior)) {
    tissue_filter_reason <- "dataset_tissue_missing"
  } else {
    llm_payload <- extract_llm_tissue_payload(head_out, tissue_prior)
    if (is.null(llm_payload) || !is.list(llm_payload$results)) {
      tissue_filter_reason <- "llm_tissue_classifier_missing_keep_all"
    } else {
      raw_results <- llm_payload$results
      cand_results <- lapply(seq_along(candidates), function(i) {
        cand <- candidates[[i]]
        cid <- paste0("cand_", i)
        clabel <- as.character(cand$label %||% "")
        hit <- NULL
        for (r in raw_results) {
          rid <- as.character(r$id %||% "")
          rlabel <- as.character(r$label %||% "")
          if (identical(rid, cid) || (nzchar(clabel) && nzchar(rlabel) && identical(normalize_token(clabel), normalize_token(rlabel)))) {
            hit <- r
            break
          }
        }
        verdict <- normalize_verdict(hit$verdict %||% "uncertain")
        conf <- suppressWarnings(as.numeric(hit$confidence %||% 0.5))
        if (!is.finite(conf)) conf <- 0.5
        reason <- as.character(hit$reason %||% "LLM classifier did not provide a candidate-specific reason.")
        list(id = cid, label = clabel, verdict = verdict, confidence = conf, reason = reason)
      })

      verdicts <- vapply(cand_results, function(x) as.character(x$verdict %||% "uncertain"), character(1))
      tissue_filter_results <- cand_results
      tissue_filter_summary <- llm_payload$summary %||% list(
        n_in = as.integer(sum(verdicts == "in_tissue")),
        n_out = as.integer(sum(verdicts == "out_of_tissue")),
        n_uncertain = as.integer(sum(verdicts == "uncertain")),
        suspect_leak = FALSE,
        suspect_leak_reason = ""
      )

      keep_idx <- which(verdicts != "out_of_tissue")
      if (length(keep_idx) > 0) {
        candidates <- candidates[keep_idx]
        tissue_filter_applied <- TRUE
        tissue_filter_reason <- "candidate_filter_by_llm_tissue_classifier"
        tissue_filter_n_after <- as.integer(length(candidates))
      } else {
        tissue_filter_reason <- "all_out_of_tissue_keep_all_for_safety"
      }
    }
  }

  if (is.null(policy_cfg)) {
    policy_cfg <- list(
      trigger = list(
        score_margin_unsteady_max = min_margin,
        method_gap_unsteady_max = 0,
        evidence_gap_unsteady_max = 0
      ),
      override = list(
        min_margin = min_margin,
        score_na_policy = "effective_surrogate"
      ),
      pool = list(
        soft_map_quality_min = 0,
        shadow_per_lineage = 1
      )
    )
  }
  validate_policy_cfg(policy_cfg)
  policy_cfg$trigger$score_margin_unsteady_max <- policy_cfg$override$min_margin

  is_generic_label <- function(clid, lbl) {
    x <- tolower(stringr::str_squish(as.character(lbl %||% "")))
    if (!nzchar(x)) return(TRUE)

    cid <- as.character(clid %||% "")
    if (!is.null(cl_graph) && nzchar(cid)) {
      d <- suppressWarnings(as.numeric(get_depth_to_root(cid, cl_graph)))
      max_desc_depth <- as.integer(policy_cfg$override$ultra_generic_desc_depth_k %||% 64L)
      desc_df <- get_descendants(cid, max_desc_depth, cl_graph)
      subtree_size <- as.integer(nrow(desc_df))
      total_nodes <- as.integer(length(cl_graph$cl %||% list()))
      subtree_ratio <- if (is.finite(total_nodes) && total_nodes > 0) as.numeric(subtree_size / total_nodes) else NA_real_
      ultra_depth_max <- as.numeric(policy_cfg$override$ultra_generic_depth_max %||% 2)
      ultra_subtree_ratio_min <- as.numeric(policy_cfg$override$ultra_generic_subtree_ratio_min %||% 0.02)
      ultra_subtree_size_min <- as.integer(policy_cfg$override$ultra_generic_subtree_size_min %||% 50L)
      if (is.finite(d)) {
        return(d <= ultra_depth_max &&
          ((is.finite(subtree_ratio) && subtree_ratio >= ultra_subtree_ratio_min) ||
             (is.finite(subtree_size) && subtree_size >= ultra_subtree_size_min)))
      }
    }

    x <- gsub("\\([^)]*\\)", "", x)
    x <- stringr::str_squish(x)
    # Text fallback only for placeholder generic labels.
    identical(x, "cell")
  }

  normalize_candidate_clid <- function(cand) {
    # Never remap a valid fixed rank-1 reviewer CL ID at this stage.
    raw <- as.character(cand$raw_clid %||% NA_character_)
    is_real_top1 <- as.integer(cand$rank %||% 99L) == 1L &&
      as.character(cand$method %||% "") != "expanded"
    raw_valid <- !is.na(raw) && nzchar(raw) && !is.null(cl_graph$cl[[raw]])
    if (isTRUE(is_real_top1) && isTRUE(raw_valid)) {
      cand$clid <- raw
      cand$coerced_status <- "frozen_input_clid"
      cand$map_quality <- max(as.numeric(cand$map_quality %||% 0), 3)
      return(cand)
    }

    label <- cand$label %||% NA_character_
    res <- tryCatch(normalize_cl_three_state(label, "", cl_cfg), error = function(e) NULL)
    norm_clid <- res$final_clid %||% NA_character_
    if (!is.na(norm_clid) && nzchar(norm_clid)) {
      if (is.na(cand$clid) || !nzchar(cand$clid)) {
        cand$clid <- norm_clid
      } else if (!identical(cand$clid, norm_clid)) {
        cand$raw_clid <- cand$clid
        cand$clid <- norm_clid
        cand$map_quality <- max(0, (cand$map_quality %||% 0) - 1)
        cand$coerced_status <- paste0(cand$coerced_status %||% "", "|label_clid_norm")
      }
    }
    cand
  }

  stage_root_clids <- get0("STAGE_ROOT_CLIDS", ifnotfound = character(0), inherits = TRUE)
  allow_stage_token_fallback <- isTRUE(get0("ALLOW_STAGE_TOKEN_FALLBACK", ifnotfound = FALSE, inherits = TRUE))
  if (is.null(stage_root_clids)) stage_root_clids <- character(0)
  stage_root_clids <- unique(as.character(stage_root_clids))
  stage_root_clids <- stage_root_clids[!is.na(stage_root_clids) & nzchar(stage_root_clids)]
  if (!is.null(cl_graph) && !is.null(cl_graph$cl)) {
    stage_root_clids <- stage_root_clids[stage_root_clids %in% names(cl_graph$cl)]
  }
  if (length(stage_root_clids) == 0 && !isTRUE(allow_stage_token_fallback)) {
    stop("STAGE_ROOT_CLIDS is empty or unavailable; strict CL-based stage policy cannot run. Define STAGE_ROOT_CLIDS in dataset_config.R or enable ALLOW_STAGE_TOKEN_FALLBACK explicitly.")
  }
  if (length(stage_root_clids) == 0 && isTRUE(allow_stage_token_fallback)) {
    warning("STAGE_ROOT_CLIDS is empty; falling back to token-based stage detection because ALLOW_STAGE_TOKEN_FALLBACK=TRUE.")
  }

  head_out <- ensure_post_issues(head_out)
  head_out <- normalize_judge_final_decision_cl_preserve(head_out, cl_cfg)
  if (is.null(head_out$audit_report) || !is.list(head_out$audit_report)) {
    head_out$audit_report <- list(reviewer_support = list(), flags = list(), notes = "")
  }

  coerce_notes <- character(0)
  for (cand in candidates) {
    raw <- cand$raw_clid %||% ""
    coerced <- cand$clid %||% ""
    status <- cand$coerced_status %||% ""
    if (nzchar(raw) || nzchar(coerced)) {
      coerce_notes <- c(coerce_notes, paste0(cand$method, ":", status, ":", raw, "->", coerced))
    }
  }

  # Determine lock from top1 clids
  top1_ctx <- compute_top1_lock_context(candidates, enable_method_reliability_gate, protect_lock_inputs)
  top1_candidates <- top1_ctx$top1_candidates
  top1_methods <- top1_ctx$top1_methods
  primary_methods <- top1_ctx$primary_methods
  lock_source_methods <- top1_ctx$lock_source_methods
  lock_source <- top1_ctx$lock_source
  top1_for_lock <- top1_ctx$top1_for_lock
  top1_clids <- top1_ctx$top1_clids
  top1_for_rules <- Filter(function(x) x$rank == 1, candidates)

  lock_ctx <- resolve_lock_from_top1_clids(top1_clids, cl_graph, min_votes)
  lock_id <- lock_ctx$lock_id
  lock_id_original <- lock_id
  lock_reason <- lock_ctx$lock_reason

  # Lock generation is unchanged; corrections are applied only after fallback resolution.

  lock_depth_ctx <- compute_lock_depth_stats(lock_id, cl_graph)
  lock_depth <- lock_depth_ctx$lock_depth
  lock_anc_count <- lock_depth_ctx$lock_anc_count

  eligible_ctx <- build_reviewer_eligibility(top1_for_rules, cl_graph, min_depth_for_eligibility)
  reviewer_diagnostics <- eligible_ctx$reviewer_diagnostics
  eligible <- eligible_ctx$eligible
  eligible_methods_for_guard <- eligible_ctx$eligible_methods_for_guard
  eligibility_drop_reasons_top <- list()
  if (length(reviewer_diagnostics) > 0) {
    rs <- vapply(reviewer_diagnostics, function(x) as.character(x$exclude_reason %||% ""), character(1))
    rs <- rs[nzchar(rs)]
    if (length(rs) > 0) {
      tb <- sort(table(rs), decreasing = TRUE)
      eligibility_drop_reasons_top <- as.list(as.integer(tb))
      names(eligibility_drop_reasons_top) <- names(tb)
    }
  }

  best_eligible_method <- pick_best_method_for_guard(eligible_methods_for_guard)
  best_other_method <- pick_best_method_for_guard(Filter(function(x) !identical(x$method, "cassia"), eligible_methods_for_guard))

  support_ctx <- compute_supported_lca_stats(eligible, gate_k, cl_graph)
  eligible_n <- support_ctx$eligible_n
  supported_lca_clid <- support_ctx$supported_lca_clid
  supported_lca_depth <- support_ctx$supported_lca_depth
  support_n <- support_ctx$support_n
  support_gate_pass <- support_ctx$support_gate_pass

  depth_gap <- if (is.finite(supported_lca_depth) && is.finite(lock_depth)) supported_lca_depth - lock_depth else NA_real_
  mode_ctx <- select_mode_by_depth_gap(lock_reason, support_gate_pass, depth_gap, delta_depth, disable_aggressive_trigger)
  mode_selected <- mode_ctx$mode_selected
  mode_trigger_reason <- mode_ctx$mode_trigger_reason
  aggressive_trigger <- mode_ctx$aggressive_trigger

  final_rule <- "balanced_default"
  anchor_fallback_happened <- FALSE
  aggressive_lca_hard_guard_applied <- FALSE
  post_fallback_override_applied <- FALSE
  post_fallback_override_reason <- ""
  post_fallback_why_not_override <- ""
  post_fallback_rule1_gate_n <- 0L
  post_fallback_rule1_soft_n <- 0L
  post_fallback_rule2_gate_n <- 0L
  post_fallback_rule2_soft_n <- 0L
  lock_unsteady <- FALSE
  lock_unsteady_reasons <- ""
  gate_pool_n <- 0L
  soft_pool_n <- 0L
  shadow_pool_n <- 0L
  gate_reject_total_n <- 0L
  gate_reject_descendant_n <- 0L
  gate_reject_descendant_by_stage <- list()
  shadow_kept_n <- 0L
  shadow_empty_reason <- "no_rejects"
  shadow_desc_filter_counts <- list(missing_clid = 0L, not_descendant = 0L, require_deeper_fail = 0L, generic_filtered = 0L, pass = 0L)
  shadow_stage_flagged_n <- 0L
  shadow_desc_ref <- list(attempted_n = 0L, kept_n = 0L)
  shadow_pool_used <- FALSE
  override_checks <- list(support_up = FALSE, support_not_down_both = TRUE, support_gate_ok = FALSE, consistency_ok = TRUE, consistency_ok_reason = "", score_ok = FALSE, score_na_case = "none", score_ok_reason = "", score_source_base = "", score_source_new = "", effective_base_score = NA_real_, effective_new_score = NA_real_, candidate_is_descendant = FALSE, base_depth = NA_real_, candidate_depth = NA_real_, base_subtree_size = 0L, base_is_ultra_generic = FALSE, stricter_desc_margin_ok = TRUE, pool_size = 0L, pool_size_1_strong_margin_used = FALSE, tie_depth_gain_used = FALSE)
  score_na_case <- "none"
  score_ok_reason <- ""
  consistency_ok <- TRUE
  consistency_ok_reason <- ""
  score_source_base <- ""
  score_source_new <- ""
  effective_base_score <- NA_real_
  effective_new_score <- NA_real_
  base_depth <- NA_real_
  candidate_depth <- NA_real_
  base_subtree_size <- 0L
  base_is_ultra_generic <- FALSE
  pool_size_1_strong_margin_used <- FALSE
  instability_pool_n <- 0L
  branch_gate_relaxed_due_to_missing_support <- NULL
  n_with_valid_clid <- 0L
  n_after_meta_reviewer_gate <- 0L
  n_after_anchor_constraints <- 0L
  n_in_cand_table <- 0L
  n_in_cand_ordered <- 0L
  guardrail_trace <- list(triggered = FALSE, gap = NA_real_, gap_ratio = NA_real_, best_method = NA_character_, best_score = NA_real_, judge_score = NA_real_)
  cassia_override_trace <- list(triggered = FALSE, reason = "disabled_in_top2_variant", cassia_clid = NA_character_, cassia_depth = NA_real_, max_other_depth = NA_real_, chosen_override_clid = NA_character_)
  label_tiebreak_trace <- list(triggered = FALSE, reason = "", candidate_labels = list(), chosen_label = NA_character_)

  if (!is.na(lock_id) && nzchar(lock_id)) {
    candidates <- lapply(candidates, function(cand) {
      cands <- cand$coerced_candidates %||% character(0)
      if (length(cands) <= 1) return(cand)
      cands <- cands[!is.na(cands) & nzchar(cands)]
      if (length(cands) == 0) return(cand)

      desc <- Filter(function(x) is_descendant_of(x, lock_id, cl_graph), cands)
      if (length(desc) > 0) {
        depths <- vapply(desc, function(x) get_depth_to_root(x, cl_graph), numeric(1))
        best <- desc[which.max(depths)][[1]]
      } else if (lock_id %in% cands) {
        best <- lock_id
      } else {
        dists <- vapply(cands, function(x) {
          get_ontology_distance(lock_id, x, cl_graph)
        }, numeric(1))
        best <- cands[which.min(dists)][[1]]
      }
      if (!is.na(best) && nzchar(best)) cand$clid <- best
      cand
    })
  }

  pool <- candidates
  if (!is.na(lock_id) && nzchar(lock_id)) {
    in_lock_pool <- Filter(function(x) in_lock_scope(x$clid, lock_id, cl_graph), candidates)
    if (length(in_lock_pool) == 0) {
      head_out$post_issues$needs_manual_review <- TRUE
      head_out$post_issues$flags <- unique(c(as.character(head_out$post_issues$flags %||% character(0)), "lock_no_in_tree"))
    }
  }

  if (length(pool) == 0) return(head_out)

  num_candidates_before <- length(pool)

  if (!is.na(lock_id) && nzchar(lock_id)) {
    expanded <- get_descendants(lock_id, expand_depth_k, cl_graph)
    if (nrow(expanded) > 0) {
      expanded_candidates <- lapply(seq_len(nrow(expanded)), function(i) {
        cid <- expanded$clid[[i]]
        list(
          label = local_lookup_by_clid(cid, cl_cfg) %||% NA_character_,
          clid = cid,
          raw_clid = NA_character_,
          coerced_status = "expanded",
          coerced_candidates = character(0),
          map_quality = 0,
          score = NA_real_,
          method = "expanded",
          rank = 99,
          depth_from_lock = expanded$depth_from_root[[i]]
        )
      })
      pool <- c(pool, expanded_candidates)
    }
  }

  num_candidates_after <- length(pool)

  # Normalize CLID by label to avoid CLID/label mismatches
  pool <- lapply(pool, normalize_candidate_clid)
  n_with_valid_clid <- as.integer(sum(vapply(pool, function(x) !is.na(x$clid) && nzchar(x$clid), logical(1))))

  # Compute evidence alignment for each method vs raw dossier
  # This allows us to weight methods by how well they match the actual evidence
  marker_map <- get_marker_map(species_value)
  tiers <- extract_gene_tiers(judge_input_obj)
  method_alignments <- compute_method_evidence_alignment(judge_input_obj, marker_map, tiers, cl_graph, cl_cfg)

  # Optional method reliability gate (default OFF to preserve historical results)
  method_reliability <- list()
  for (m in top1_methods) {
    method_reliability[[m]] <- list(score = 1, reasons = character(0), is_primary = m %in% c("cassia", "our", "enrich"))
  }
  if (isTRUE(enable_method_reliability_gate) && length(top1_methods) > 0) {
    top1_by_method <- list()
    for (m in top1_methods) {
      cm <- Filter(function(x) (x$method %||% "") == m && (x$rank %||% 99L) == 1L, candidates)
      if (length(cm) > 0) top1_by_method[[m]] <- cm[[1]]
    }
    for (m in top1_methods) {
      if (m %in% c("cassia", "our", "enrich")) {
        method_reliability[[m]] <- list(score = 1, reasons = character(0), is_primary = TRUE)
        next
      }
      cand <- top1_by_method[[m]] %||% list()
      m_clid <- as.character(cand$clid %||% NA_character_)
      m_mq <- suppressWarnings(as.numeric(cand$map_quality %||% NA_real_))
      m_lbl <- as.character(cand$label %||% "")
      has_clid <- !is.na(m_clid) && nzchar(m_clid)
      map_score <- if (!is.na(m_mq)) max(0, min(1, m_mq / 3)) else 0
      generic_flag <- is_generic_label(m_clid, m_lbl)
      dist_to_lock <- NA_real_
      dist_score <- 0
      if (has_clid && !is.na(lock_id) && nzchar(lock_id)) {
        dist_to_lock <- suppressWarnings(as.numeric(get_ontology_distance(lock_id, m_clid, cl_graph)))
        if (is.finite(dist_to_lock)) {
          dist_score <- max(0, 1 - min(dist_to_lock, 8) / 8)
        }
      }
      rel <- (0.35 * ifelse(has_clid, 1, 0)) + (0.30 * map_score) + (0.25 * dist_score) + (0.10 * ifelse(generic_flag, 0, 1))
      reasons <- character(0)
      if (!has_clid) reasons <- c(reasons, "no_clid")
      if (generic_flag) reasons <- c(reasons, "generic_label")
      if (is.finite(dist_to_lock) && dist_to_lock >= 6) reasons <- c(reasons, "lineage_far_from_lock")
      if (length(reasons) > 0) rel <- 0.01
      method_reliability[[m]] <- list(score = max(0, min(1, rel)), reasons = unique(reasons), is_primary = FALSE)
    }
  }

  # Optional hard gate for non-core methods (default OFF to preserve historical results).
  # Purpose: automatically drop newly injected methods that look noisy for this cluster.
  meta_gate_res <- apply_meta_reviewer_gate(candidates, candidates_raw, top1_methods, method_reliability, method_alignments, meta_reviewer_gate)
  candidates <- meta_gate_res$candidates
  n_after_meta_reviewer_gate <- as.integer(length(candidates))
  top1_candidates <- meta_gate_res$top1_candidates
  top1_methods <- meta_gate_res$top1_methods
  meta_reviewer_decisions <- meta_gate_res$meta_reviewer_decisions
  
  # method_count per clid - weighted by evidence alignment
  # Keep this dynamic so future reviewer methods are included automatically.
  real_input_methods <- unique(vapply(candidates, function(x) as.character(x$method %||% ""), character(1)))
  real_input_methods <- real_input_methods[nzchar(real_input_methods)]
  method_weights <- list(
    cassia = method_alignments$cassia$score %||% 0,
    our = method_alignments$our$score %||% 0,
    enrich = method_alignments$enrich$score %||% 0
  )
  if (isTRUE(enable_method_reliability_gate)) {
    for (m in top1_methods) {
      if (m %in% c("cassia", "our", "enrich")) next
      rel <- suppressWarnings(as.numeric(method_reliability[[m]]$score %||% 0.01))
      base_align <- suppressWarnings(as.numeric(method_alignments[[m]]$score %||% NA_real_))
      if (!is.finite(base_align)) base_align <- rel
      method_weights[[m]] <- max(0, min(1, base_align)) * max(0, min(1, rel))
    }
  }
  
  method_count <- list()
  for (cand in pool) {
    clid <- cand$clid %||% NA_character_
    if (is.na(clid) || !nzchar(clid)) next
    method_name <- cand$method %||% ""
    # Only count real input methods, not "expanded" or other synthetic sources
    if (!method_name %in% real_input_methods) next
    key <- clid
    if (is.null(method_count[[key]])) method_count[[key]] <- 0
    weight <- method_weights[[method_name]] %||% 0
    method_count[[key]] <- method_count[[key]] + weight
  }

  # method_count_all per clid - includes expanded/synthetic sources
  method_count_all <- list()
  for (cand in pool) {
    clid <- cand$clid %||% NA_character_
    if (is.na(clid) || !nzchar(clid)) next
    method_name <- cand$method %||% ""
    key <- clid
    if (is.null(method_count_all[[key]])) method_count_all[[key]] <- list()
    method_count_all[[key]] <- unique(c(method_count_all[[key]], method_name))
  }

  # Use already-computed marker_map and tiers from method alignment calculation
  tier_anchor_min_hits <- 2
  tier_anchor_max_top_hits <- 1
  anchor_boost <- 0.5
  contamination_penalty_weight <- 0.5
  guardrail_gap <- 2
  guardrail_penalty_weight <- 0.5

  cand_table <- lapply(pool, function(cand) {
    clid <- cand$clid %||% NA_character_
    # method_count stores weighted evidence alignment (numeric)
    mcnt <- if (!is.na(clid) && !is.null(method_count[[clid]])) method_count[[clid]] else 0
    # method_count_all stores list of all methods including expanded
    mcnt_all <- if (!is.na(clid) && !is.null(method_count_all[[clid]])) length(method_count_all[[clid]]) else max(1, mcnt)
    mq_base <- as.numeric(cand$map_quality %||% 0)
    mq_rank <- max(0, mq_base - (1.5 * (cand$rank %||% 1L - 1L)))
    depth_to_root <- if (!is.na(clid) && nzchar(clid)) get_depth_to_root(clid, cl_graph) else 0
    dist_to_lock <- if (!is.na(lock_id) && nzchar(lock_id) && !is.na(clid) && nzchar(clid)) {
      get_ontology_distance(lock_id, clid, cl_graph)
    } else NA_real_
    in_lock_candidate <- if (!is.na(lock_id) && nzchar(lock_id) && !is.na(clid) && nzchar(clid)) {
      isTRUE(is_ancestor_of(lock_id, clid, cl_graph)) || isTRUE(clid == lock_id)
    } else FALSE
    depth_from_lock <- if (!is.null(cand$depth_from_lock)) {
      cand$depth_from_lock
    } else if (isTRUE(in_lock_candidate) && is.finite(dist_to_lock)) {
      dist_to_lock
    } else 0
    apply_depth_penalty <- is.na(cand$score) || (!is.na(cand$score) && cand$score < strong_score) || (mq_base < 2)
    penalty_alpha <- score_lambda
    penalty_beta <- score_lambda * lock_beta
    dist_pen <- if (is.finite(dist_to_lock)) dist_to_lock else (lock_depth %||% 5) + 5
    lock_penalty <- 0
    if (!is.na(lock_id) && nzchar(lock_id)) {
      if (isTRUE(in_lock_candidate)) {
        lock_penalty <- if (apply_depth_penalty) penalty_alpha * depth_from_lock else 0
      } else {
        lock_penalty <- (if (apply_depth_penalty) penalty_alpha * depth_from_lock else 0) + (penalty_beta * dist_pen)
      }
    }
    evidence_stats <- compute_evidence_support(cand$label, marker_map, tiers)
    anchor_bonus <- if (evidence_stats$recovery >= tier_anchor_min_hits && evidence_stats$top <= tier_anchor_max_top_hits) anchor_boost else 0
    score_adj <- if (!is.na(cand$score)) cand$score - lock_penalty + anchor_bonus else NA_real_
    rel_score <- suppressWarnings(as.numeric(method_reliability[[cand$method %||% ""]]$score %||% 1))
    if (!is.finite(rel_score)) rel_score <- 1
    rel_score <- max(0, min(1, rel_score))
    rel_penalty <- 0
    if (isTRUE(enable_method_reliability_gate) && !((cand$method %||% "") %in% c("cassia", "our", "enrich"))) {
      rel_penalty <- (1 - rel_score) * 1.5
      if (is.na(score_adj)) {
        score_adj <- -rel_penalty
      } else {
        score_adj <- score_adj - rel_penalty
      }
    }
    list(
      label = cand$label,
      clid = clid,
      score = cand$score,
      score_adj = score_adj,
      method = cand$method,
      rank = cand$rank,
      method_count = mcnt,
      method_count_all = mcnt_all,
      map_quality = mq_base,
      map_quality_rank = mq_rank,
      specificity = depth_to_root,
      depth_from_lock = depth_from_lock,
      dist_to_lock = dist_to_lock,
      depth_to_root = depth_to_root,
      lock_penalty = lock_penalty,
      evidence_support = evidence_stats$score,
      evidence_top_hits = evidence_stats$top,
      evidence_recovery_hits = evidence_stats$recovery,
      evidence_unknown_hits = evidence_stats$unknown,
      evidence_out_hits = evidence_stats$out,
      anchor_bonus = anchor_bonus,
      reliability_score = rel_score,
      reliability_penalty = rel_penalty,
      reliability_reasons = method_reliability[[cand$method %||% ""]]$reasons %||% character(0)
    )
  })

  max_support <- max(vapply(cand_table, function(x) x$evidence_support, numeric(1)), na.rm = TRUE)
  if (!is.finite(max_support)) max_support <- 0

  cand_table <- lapply(cand_table, function(cand) {
    guardrail_penalty <- 0
    if (max_support > 0 && cand$evidence_support <= (max_support - guardrail_gap)) {
      guardrail_penalty <- guardrail_penalty_weight * (max_support - cand$evidence_support)
    }
    contamination_penalty <- contamination_penalty_weight * (cand$evidence_out_hits %||% 0)
    cand$score_adj <- if (!is.na(cand$score_adj)) cand$score_adj - guardrail_penalty - contamination_penalty else cand$score_adj
    cand$guardrail_penalty <- guardrail_penalty
    cand$contamination_penalty <- contamination_penalty
    cand
  })

  # Lineage consistency + specificity balance (generic, no hardcoded biology)
  anchor_support_threshold <- if (max_support > 0) max_support * 0.6 else NA_real_
  strong_support_threshold <- if (max_support > 0) max_support * 0.8 else NA_real_
  anchor_candidates <- Filter(function(x) {
    !is.na(x$clid) && nzchar(x$clid) && is.finite(x$evidence_support) &&
      max_support > 0 && x$evidence_support >= anchor_support_threshold &&
      ((x$method_count %||% 0) >= min_votes || (x$map_quality %||% 0) >= 2)
  }, cand_table)
  anchor_clids <- unique(vapply(anchor_candidates, function(x) x$clid, character(1)))
  strong_anchor_candidates <- Filter(function(x) {
    !is.na(x$clid) && nzchar(x$clid) && is.finite(x$evidence_support) &&
      max_support > 0 && x$evidence_support >= strong_support_threshold
  }, cand_table)
  strong_anchor_clids <- unique(vapply(strong_anchor_candidates, function(x) x$clid, character(1)))
  dc_anchor_clid <- top1_ctx$direct_consensus_clid %||% NA_character_

  # Direct majority is a SOFT reviewer anchor. It is not automatically promoted to
  # strong_anchor_clids and therefore cannot by itself trigger the hard branch override.
  if (!is.na(dc_anchor_clid) && nzchar(as.character(dc_anchor_clid %||% ""))) {
    anchor_clids <- unique(c(anchor_clids, as.character(dc_anchor_clid)))
  }

  anchor_support <- list()
  if (length(anchor_candidates) > 0) {
    for (a in anchor_candidates) anchor_support[[a$clid]] <- max(anchor_support[[a$clid]] %||% 0, a$evidence_support %||% 0)
  }
  if (!is.na(dc_anchor_clid) && nzchar(as.character(dc_anchor_clid %||% ""))) {
    dc_ev <- vapply(Filter(function(x) {
      !is.na(x$clid) && nzchar(x$clid) && identical(as.character(x$clid), as.character(dc_anchor_clid))
    }, cand_table), function(x) as.numeric(x$evidence_support %||% 0), numeric(1))
    dc_ev <- dc_ev[is.finite(dc_ev)]
    dc_ev_max <- if (length(dc_ev) > 0L) max(dc_ev) else 0
    anchor_support[[as.character(dc_anchor_clid)]] <- max(anchor_support[[as.character(dc_anchor_clid)]] %||% 0, dc_ev_max)
  }

  anchor_depths <- list()
  if (length(anchor_clids) > 0) {
    for (a in anchor_clids) anchor_depths[[a]] <- get_depth_to_root(a, cl_graph) %||% 0
  }
  lineage_penalty_weight <- 0.4
  over_specific_weight <- 0.2
  over_conservative_weight <- 0.2
  anchor_strong_threshold <- if (max_support > 0) max_support * 0.8 else Inf

  cand_table <- lapply(cand_table, function(cand) {
    lineage_penalty <- 0
    over_specific_penalty <- 0
    over_conservative_penalty <- 0
    anchor_dist_min <- NA_real_
    if (length(anchor_clids) > 0 && !is.na(cand$clid) && nzchar(cand$clid)) {
      dists <- vapply(anchor_clids, function(a) get_ontology_distance(cand$clid, a, cl_graph), numeric(1))
      dists <- dists[is.finite(dists)]
      if (length(dists) > 0) {
        anchor_dist_min <- min(dists)
        lineage_penalty <- max(0, anchor_dist_min - 1) * lineage_penalty_weight

         # Over-specific: candidate much deeper than nearest strong anchor with weak evidence
         # Depth-cap penalty now tied to anchor strength: weak anchor allows more depth, strong anchor is strict
         nearest_anchor <- anchor_clids[[which.min(dists)]]
         anchor_depth <- anchor_depths[[nearest_anchor]] %||% 0
         anchor_strength <- anchor_support[[nearest_anchor]] %||% 0
         depth_gap <- (cand$depth_to_root %||% 0) - anchor_depth
         if (depth_gap >= 2 && (cand$evidence_support %||% 0) < anchor_strong_threshold) {
           strength_factor <- max(0.5, anchor_strength / anchor_strong_threshold)
           over_specific_penalty <- (depth_gap * over_specific_weight) / strength_factor
         }

        # Over-conservative: candidate is ancestor of a strong anchor
        if (!is.na(cand$clid) && nzchar(cand$clid)) {
          strong_anchors <- anchor_clids[vapply(anchor_clids, function(a) (anchor_support[[a]] %||% 0) >= anchor_strong_threshold, logical(1))]
          if (length(strong_anchors) > 0) {
            anchor_depths_desc <- vapply(strong_anchors, function(a) {
              if (isTRUE(is_ancestor_of(cand$clid, a, cl_graph))) anchor_depths[[a]] %||% 0 else NA_real_
            }, numeric(1))
            anchor_depths_desc <- anchor_depths_desc[is.finite(anchor_depths_desc)]
            if (length(anchor_depths_desc) > 0) {
              depth_gap_desc <- max(anchor_depths_desc) - (cand$depth_to_root %||% 0)
              if (depth_gap_desc >= 2) {
                over_conservative_penalty <- depth_gap_desc * over_conservative_weight
          }
        }
      }
    }
  }  # Descendant override remains disabled.
    }

    total_lineage_penalty <- lineage_penalty + over_specific_penalty + over_conservative_penalty
    cand$score_adj <- if (!is.na(cand$score_adj)) cand$score_adj - total_lineage_penalty else cand$score_adj
    cand$lineage_penalty <- lineage_penalty
    cand$over_specific_penalty <- over_specific_penalty
    cand$over_conservative_penalty <- over_conservative_penalty
    cand$anchor_dist_min <- anchor_dist_min
    cand
  })

  # Evidence anchor hard constraint: if strong anchors exist, only keep candidates within distance 1
  cand_table_all <- cand_table
  if (length(strong_anchor_clids) > 0) {
    cand_table_strong <- Filter(function(cand) {
      if (is.na(cand$clid) || !nzchar(cand$clid)) return(FALSE)
      dists <- vapply(strong_anchor_clids, function(a) get_ontology_distance(cand$clid, a, cl_graph), numeric(1))
      dists <- dists[is.finite(dists)]
      length(dists) > 0 && min(dists) <= 1
    }, cand_table)
    if (length(cand_table_strong) > 0) {
      cand_table <- cand_table_strong
    } else {
      cand_table <- cand_table_all
    }
  }
  n_after_anchor_constraints <- as.integer(length(cand_table))
  n_in_cand_table <- as.integer(length(cand_table))

  # If strong anchors span multiple branches, cap specificity to avoid over-specific selection
  mixed_lineage <- FALSE
  if (length(strong_anchor_clids) > 1) {
    dist_mat <- c()
    for (i in seq_along(strong_anchor_clids)) {
      for (j in seq_along(strong_anchor_clids)) {
        if (i >= j) next
        d <- get_ontology_distance(strong_anchor_clids[[i]], strong_anchor_clids[[j]], cl_graph)
        if (is.finite(d)) dist_mat <- c(dist_mat, d)
      }
    }
    if (length(dist_mat) > 0 && max(dist_mat) > 2) mixed_lineage <- TRUE
  }
  if (mixed_lineage && length(anchor_depths) > 0) {
    min_anchor_depth <- min(vapply(anchor_depths, function(x) x %||% 0, numeric(1)))
    depth_cap <- min_anchor_depth + 1
    cand_table <- lapply(cand_table, function(cand) {
      if (!is.na(cand$depth_to_root) && cand$depth_to_root > depth_cap) {
        cand$score_adj <- if (!is.na(cand$score_adj)) cand$score_adj - (cand$depth_to_root - depth_cap) else cand$score_adj
      }
      cand
    })
  }

  # If a strong anchor exists, do not allow ancestors far above it to win (prevents over-conservative fallback)
  if (length(strong_anchor_clids) > 0) {
    cand_table <- Filter(function(cand) {
      if (is.na(cand$clid) || !nzchar(cand$clid)) return(TRUE)
      for (a in strong_anchor_clids) {
        if (isTRUE(is_ancestor_of(cand$clid, a, cl_graph))) {
          depth_gap <- (anchor_depths[[a]] %||% 0) - (cand$depth_to_root %||% 0)
          if (depth_gap >= 2) return(FALSE)
        }
      }
      TRUE
    }, cand_table)
    if (length(cand_table) == 0) cand_table <- cand_table_all
  }

  greedy_target <- NA_character_
  greedy_label <- NA_character_
  
  # Extract LLM's anchor: the original final_decision before post-processing
  # This is the LLM's informed decision that we should respect
  llm_anchor_clid <- as.character(head_out$final_decision$final_cell_ontology_id %||% NA_character_)
  llm_anchor_label <- as.character(head_out$final_decision$primary_cell_type %||% NA_character_)
  head_identity_valid <- !is.na(llm_anchor_clid) && nzchar(llm_anchor_clid) &&
    llm_anchor_clid %in% names(cl_graph$cl) && !is.na(llm_anchor_label) && nzchar(llm_anchor_label)
  force_preserve_head_identity <- identical(identity_policy, "head_preserve") && isTRUE(head_identity_valid)
  ontology_guarded_identity <- identical(identity_policy, "ontology_guarded") && isTRUE(head_identity_valid)
  head_policy_active <- isTRUE(force_preserve_head_identity) || isTRUE(ontology_guarded_identity)
  # Internal alias for the full-preserve identity mode.
  preserve_head_identity <- isTRUE(force_preserve_head_identity)
  
  # Extract reviewer_support from LLM's audit_report to identify unsupported methods
  # Methods marked as unsupported should NOT have their CLIDs selected as greedy_target
  reviewer_support <- head_out$audit_report$reviewer_support %||% list()
  normalize_support_key <- function(k) {
    kk <- sub("_supported$", "", as.character(k))
    if (identical(kk, "our_method")) kk <- "our"
    kk
  }
  unsupported_methods <- character(0)
  if (is.list(reviewer_support) && length(reviewer_support) > 0) {
    for (k in names(reviewer_support)) {
      if (!grepl("_supported$", k)) next
      if (isFALSE(reviewer_support[[k]])) unsupported_methods <- c(unsupported_methods, normalize_support_key(k))
    }
  }
  unsupported_methods <- unique(unsupported_methods)
  
  # Build map of CLIDs to originating methods (ALL candidates from real input methods)
  # This helps identify if a CLID only comes from unsupported methods
  clid_to_methods <- list()
  for (cand in pool) {
    method_name <- cand$method %||% ""
    # Only track real input methods
    if (!method_name %in% real_input_methods) next
    if (!is.na(cand$clid) && nzchar(cand$clid)) {
      key <- cand$clid
      if (is.null(clid_to_methods[[key]])) clid_to_methods[[key]] <- character(0)
      clid_to_methods[[key]] <- unique(c(clid_to_methods[[key]], method_name))
    }
  }
  
  # Helper: check if a CLID should be excluded because it ONLY comes from unsupported methods
  is_unsupported_only <- function(clid) {
    if (is.na(clid) || !nzchar(clid)) return(FALSE)
    methods_for_clid <- clid_to_methods[[clid]] %||% character(0)
    if (length(methods_for_clid) == 0) return(FALSE)  # expanded/synthetic candidates are OK
    # If ALL methods that produced this CLID are unsupported AND alignment is weak, exclude it
    if (!all(methods_for_clid %in% unsupported_methods)) return(FALSE)
    align_scores <- vapply(methods_for_clid, function(m) {
      if (!is.null(method_alignments[[m]])) method_alignments[[m]]$score %||% 0 else 0
    }, numeric(1))
    all(align_scores < 0.3)
  }
  
  if (length(cand_table) > 0) {
    # method_count is now weighted sum (0-3 max), use lower threshold
    c2 <- Filter(function(x) !is.na(x$clid) && x$method_count >= 0.5, cand_table)
    if (length(c2) > 0) {
      # GUARD 1: Filter out CLIDs that ONLY come from unsupported methods
      # This uses the LLM's own evidence assessment (reviewer_support)
      c2_supported <- Filter(function(x) !is_unsupported_only(x$clid), c2)
      if (length(c2_supported) == 0) {
        # All candidates are from unsupported methods; fall back to original pool
        c2_supported <- c2
      }
      
      # Sort by score_adj first, with specificity as a weak tiebreaker.
      # This prevents deep but biologically wrong types from winning
      ord <- order(vapply(c2_supported, function(x) x$score_adj %||% -Inf, numeric(1)),
                   vapply(c2_supported, function(x) x$specificity, numeric(1)),
                   decreasing = TRUE)
      best_candidate <- c2_supported[[ord[[1]]]]
      
      # GUARD 2: Don't override LLM's decision with a progenitor/precursor type
      # The LLM's decision (e.g., "MATURE_CELL_TYPE") is mature;
      # don't replace it with developmental types (e.g., "PROGENITOR_LIKE_TYPE")
      if (!is.na(llm_anchor_clid) && nzchar(llm_anchor_clid)) {
      llm_is_progenitor <- is_developmental_stage(clid = llm_anchor_clid, label = llm_anchor_label, cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)
      best_is_progenitor <- is_developmental_stage(clid = best_candidate$clid, label = best_candidate$label, cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)
        
        # If LLM chose mature type but greedy wants progenitor, keep LLM's choice
        if (!llm_is_progenitor && best_is_progenitor) {
          # Check if LLM's choice is in the candidate pool
          llm_in_pool <- any(vapply(c2_supported, function(x) identical(x$clid, llm_anchor_clid), logical(1)))
          if (llm_in_pool) {
            greedy_target <- llm_anchor_clid
            greedy_label <- llm_anchor_label
          } else {
            # LLM's choice not in multi-method pool; use best non-progenitor
        non_prog <- Filter(function(x) !is_developmental_stage(clid = x$clid, label = x$label, cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback), c2_supported)
            if (length(non_prog) > 0) {
              ord2 <- order(vapply(non_prog, function(x) x$score_adj %||% -Inf, numeric(1)),
                           vapply(non_prog, function(x) x$specificity, numeric(1)),
                           decreasing = TRUE)
              greedy_target <- non_prog[[ord2[[1]]]]$clid
              greedy_label <- non_prog[[ord2[[1]]]]$label
            } else {
              greedy_target <- best_candidate$clid
              greedy_label <- best_candidate$label
            }
          }
        } else {
          greedy_target <- best_candidate$clid
          greedy_label <- best_candidate$label
        }
      } else {
        greedy_target <- best_candidate$clid
        greedy_label <- best_candidate$label
      }
    }
  }

  cand_filtered <- cand_table
  consensus_mode <- "all"
  multi_topk <- Filter(function(x) (x$method_count %||% 0) >= min_votes, cand_table)
  if (length(multi_topk) > 0) {
    cand_filtered <- multi_topk
    consensus_mode <- "topk"
  } else {
    multi_extended <- Filter(function(x) (x$method_count_all %||% 0) >= min_votes, cand_table)
    if (length(multi_extended) > 0) {
      cand_filtered <- multi_extended
      consensus_mode <- "extend_topk"
    }
  }
  high_quality <- any(vapply(cand_table, function(x) (x$map_quality %||% 0) >= 2, logical(1)))
  if (high_quality) {
    best_high <- max(vapply(cand_table, function(x) if ((x$map_quality %||% 0) >= 2 && !is.na(x$score_adj)) x$score_adj else -Inf, numeric(1)))
    cand_filtered <- Filter(function(x) {
      if ((x$map_quality %||% 0) >= 2) return(TRUE)
      if ((x$map_quality %||% 0) == 1 && !is.na(x$score_adj) && is.finite(best_high) && x$score_adj >= (best_high - min_margin)) return(TRUE)
      FALSE
    }, cand_table)
  }

  if (length(cand_filtered) == 0) cand_filtered <- cand_table

  # Keep the soft direct-majority anchor visible even if an intermediate filter drops it.
  # Reintroduction only restores candidacy; it does not force the final label.
  direct_consensus_anchor_reintroduced <- FALSE
  if (!is.na(dc_anchor_clid) && nzchar(as.character(dc_anchor_clid %||% ""))) {
    dc_id <- as.character(dc_anchor_clid)
    has_dc <- any(vapply(cand_filtered, function(x) {
      !is.na(x$clid) && identical(as.character(x$clid), dc_id)
    }, logical(1)))
    if (!isTRUE(has_dc)) {
      dc_top1 <- Filter(function(x) {
        !is.na(x$clid) && identical(as.character(x$clid), dc_id) &&
          as.integer(x$rank %||% 99L) == 1L &&
          as.character(x$method %||% "") %in% real_input_methods
      }, cand_table_all)
      if (length(dc_top1) > 0L) {
        cand_filtered <- c(cand_filtered, dc_top1)
        direct_consensus_anchor_reintroduced <- TRUE
      }
    }
  }

  if (length(cand_filtered) > 1L) {
    dedup_key <- vapply(cand_filtered, function(x) paste(
      as.character(x$method %||% ""), as.integer(x$rank %||% 99L),
      as.character(x$clid %||% ""), sep = "|"
    ), character(1))
    cand_filtered <- cand_filtered[!duplicated(dedup_key)]
  }

  cand_ordered <- order_candidates_by_priority(cand_filtered, policy_cfg$override)
  n_in_cand_ordered <- as.integer(length(cand_ordered))
  if (length(cand_ordered) > 0) {
    cand_ordered[[1]]$top_ranked <- TRUE
  }

  anchor_clid <- compute_anchor_clid_from_candidates(candidates, lock_id, cl_graph)
  drop_reasons <- names(eligibility_drop_reasons_top %||% list())
  has_input_drop <- any(grepl("missing_or_invalid_clid", drop_reasons)) ||
    any(grepl("^depth_below_threshold:", drop_reasons))
  fragile_anchor_case <- FALSE
  branch_gate_relaxed_due_to_missing_support <- isTRUE(fragile_anchor_case)
  rejected_by_gate <- character(0)
  gate_decisions <- lapply(cand_ordered, function(cand) {
    is_desc <- FALSE
    if (!is.na(lock_id) && nzchar(lock_id) && !is.na(cand$clid) && nzchar(cand$clid)) {
      is_desc <- isTRUE(is_descendant_of(cand$clid, lock_id, cl_graph)) || identical(cand$clid, lock_id)
    }
    ok <- TRUE
    if (!isTRUE(fragile_anchor_case && isTRUE(is_desc))) {
      ok <- evaluate_branch_consistency_gate(cand$clid, anchor_clid, cl_graph, anchor_support, anchor_strong_threshold, branch_gate_slack)
    } else {
      if (is_generic_label(cand$clid, cand$label %||% "")) ok <- FALSE
    }
    lineage_key <- if (!is.na(cand$clid) && nzchar(cand$clid)) {
      p <- get_parent_clid(cand$clid, cl_graph)
      if (!is.na(p) && nzchar(p)) p else cand$clid
    } else {
      ""
    }
    stage <- if (isTRUE(ok)) "none" else "branch_gate"
    reason <- if (isTRUE(ok)) {
      if (isTRUE(fragile_anchor_case && isTRUE(is_desc))) "branch_gate_relaxed_narrow" else "none"
    } else {
      "branch_gate_reject"
    }
    if (!ok) {
      rejected_by_gate <<- c(rejected_by_gate, paste0(cand$label %||% "", "|", cand$clid %||% ""))
    }
    list(
      kept = isTRUE(ok),
      reject_stage = stage,
      reject_reason = reason,
      is_descendant = isTRUE(is_desc),
      lineage_key = lineage_key,
      candidate = cand
    )
  })
  cand_ordered_gate <- lapply(Filter(function(d) isTRUE(d$kept), gate_decisions), function(d) d$candidate)

  # ---- Evidence-qualified consensus override helpers ----
  # IMPORTANT: do NOT enforce this gate on greedy_target. greedy_target can be NA when
  # reviewer scores are unavailable, while a later fallback / method-alignment path can
  # still select a real final candidate. The gate is evaluated once, late, against the
  # ACTUAL chosen_clid after all candidate-selection and decision-alignment steps.
  consensus_override_attempted <- FALSE
  consensus_override_allowed <- NA
  consensus_override_reason <- ""
  override_candidate_clid <- NA_character_
  override_relation_to_anchor <- ""
  contrastive_shared_score <- NA_real_
  contrastive_candidate_unique_score <- NA_real_
  contrastive_anchor_unique_score <- NA_real_
  contrastive_margin <- NA_real_
  contrastive_candidate_total_score <- NA_real_
  contrastive_anchor_total_score <- NA_real_
  contrastive_candidate_real_methods <- character(0)
  consensus_override_forced_anchor <- FALSE
  consensus_override_gate_stage <- "final_pre_canonicalization"
  contrastive_override_min_margin <- 0.5

  score_observed_markers <- function(markers) {
    markers <- normalize_gene_list(markers)
    if (length(markers) == 0L) {
      return(list(score = 0, top = 0L, recovery = 0L, unknown = 0L, genes = character(0)))
    }
    top_genes <- intersect(markers, tiers$in_scope_top %||% character(0))
    rec_genes <- intersect(markers, tiers$recovery %||% character(0))
    unk_genes <- intersect(markers, tiers$unknown %||% character(0))
    score <- (2.0 * length(top_genes)) + (1.0 * length(rec_genes)) + (0.5 * length(unk_genes))
    list(
      score = as.numeric(score),
      top = as.integer(length(top_genes)),
      recovery = as.integer(length(rec_genes)),
      unknown = as.integer(length(unk_genes)),
      genes = unique(c(top_genes, rec_genes, unk_genes))
    )
  }

  markers_for_clid <- function(clid) {
    if (is.null(clid) || is.na(clid) || !nzchar(as.character(clid))) return(character(0))
    lbl <- local_lookup_by_clid(as.character(clid), cl_cfg) %||% ""
    key <- normalize_candidate_label(lbl)
    if (is.na(key) || !nzchar(key) || !key %in% names(marker_map)) return(character(0))
    normalize_gene_list(marker_map[[key]] %||% character(0))
  }

  evaluate_final_consensus_override <- function(candidate_clid) {
    res <- list(
      attempted = FALSE,
      allowed = NA,
      reason = "",
      candidate_clid = as.character(candidate_clid %||% NA_character_),
      relation = "",
      shared_score = NA_real_,
      candidate_unique_score = NA_real_,
      anchor_unique_score = NA_real_,
      margin = NA_real_,
      candidate_total_score = NA_real_,
      anchor_total_score = NA_real_,
      candidate_real_methods = character(0)
    )

    anchor <- as.character(dc_anchor_clid %||% NA_character_)
    cand <- as.character(candidate_clid %||% NA_character_)
    if (is.na(anchor) || !nzchar(anchor) || is.na(cand) || !nzchar(cand)) return(res)
    if (identical(cand, anchor)) {
      res$relation <- "same_clid"
      res$reason <- "final_matches_direct_consensus_anchor"
      return(res)
    }

    related <- isTRUE(is_descendant_of(cand, anchor, cl_graph)) ||
      isTRUE(is_ancestor_of(cand, anchor, cl_graph))
    if (related) {
      res$relation <- "ancestor_or_descendant"
      res$reason <- "hierarchical_refinement_or_broadening"
      return(res)
    }

    # Non-hierarchical / sibling-branch change relative to direct majority.
    res$attempted <- TRUE
    res$relation <- "incompatible_branch"
    res$candidate_real_methods <- as.character(clid_to_methods[[cand]] %||% character(0))

    cand_markers <- markers_for_clid(cand)
    anchor_markers <- markers_for_clid(anchor)
    if (length(cand_markers) == 0L || length(anchor_markers) == 0L) {
      # Marker reference is insufficient to make a deterministic override decision.
      # Leave the chosen label intact and escalate to automated Chief QC instead of
      # silently forcing the majority anchor.
      res$allowed <- NA
      res$reason <- "contrastive_marker_reference_unavailable_requires_qc"
      return(res)
    }

    shared_markers <- intersect(cand_markers, anchor_markers)
    cand_unique_markers <- setdiff(cand_markers, anchor_markers)
    anchor_unique_markers <- setdiff(anchor_markers, cand_markers)

    shared_ev <- score_observed_markers(shared_markers)
    cand_unique_ev <- score_observed_markers(cand_unique_markers)
    anchor_unique_ev <- score_observed_markers(anchor_unique_markers)

    res$shared_score <- shared_ev$score
    res$candidate_unique_score <- cand_unique_ev$score
    res$anchor_unique_score <- anchor_unique_ev$score
    res$candidate_total_score <- shared_ev$score + cand_unique_ev$score
    res$anchor_total_score <- shared_ev$score + anchor_unique_ev$score
    res$margin <- cand_unique_ev$score - anchor_unique_ev$score

    has_real_candidate_source <- length(res$candidate_real_methods) > 0L
    override_ok <- isTRUE(has_real_candidate_source) &&
      is.finite(res$candidate_unique_score) &&
      is.finite(res$anchor_unique_score) &&
      (res$candidate_unique_score > (res$anchor_unique_score + contrastive_override_min_margin))

    res$allowed <- isTRUE(override_ok)
    if (!has_real_candidate_source) {
      res$reason <- "override_not_justified_synthetic_or_unreviewed_candidate"
    } else if (isTRUE(override_ok)) {
      res$reason <- paste0(
        "contrastive_override_allowed: candidate_unique=", sprintf("%.2f", res$candidate_unique_score),
        " > anchor_unique=", sprintf("%.2f", res$anchor_unique_score),
        " + margin=", sprintf("%.2f", contrastive_override_min_margin),
        "; shared_observed=", sprintf("%.2f", res$shared_score)
      )
    } else {
      res$reason <- paste0(
        "override_not_justified: candidate_unique=", sprintf("%.2f", res$candidate_unique_score),
        " <= anchor_unique=", sprintf("%.2f", res$anchor_unique_score),
        " + margin=", sprintf("%.2f", contrastive_override_min_margin),
        "; shared_observed=", sprintf("%.2f", res$shared_score)
      )
    }
    res
  }

  chosen <- NULL
  chosen_support <- list(pass = FALSE, reasons = c("no_candidate"), margin = NA_real_)
  specificity_tiebreak <- FALSE

  if (!is.na(greedy_target) && nzchar(greedy_target)) {
    greedy_pool <- Filter(function(x) x$clid == greedy_target, cand_ordered_gate)
    if (length(greedy_pool) > 0) {
      greedy_best <- greedy_pool[[1]]
      support <- evaluate_candidate_support(greedy_best, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
      if (support$pass) {
        chosen <- greedy_best
        chosen_support <- support
      }
    }
  }

  if (is.null(chosen)) {
    for (cand in cand_ordered_gate) {
      if (is.na(cand$clid) && !is.na(lock_id) && nzchar(lock_id)) next
      support <- evaluate_candidate_support(cand, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
      if (support$pass) {
        chosen <- cand
        chosen_support <- support
        break
      }
    }
  }

  if (is.null(chosen) && !is.na(anchor_clid) && nzchar(anchor_clid)) {
    fallback_clid <- if (!is.na(lock_id) && nzchar(lock_id)) lock_id else anchor_clid
    if (identical(mode_selected, "aggressive")) {
      if (!is.na(supported_lca_clid) && nzchar(supported_lca_clid)) {
        fallback_clid <- supported_lca_clid
        final_rule <- "aggressive_fallback_lca"
      } else {
      fallback_clid <- pick_best_supported_clid_from_eligible(eligible)
        final_rule <- "aggressive_pick_best_supported"
      }
    }

    # Parity guardrail: avoid collapsing to very generic fallback labels when there are
    # specific, supported descendants available in the gated candidate pool.
    if (isTRUE(parity_mode) && !is.na(fallback_clid) && nzchar(fallback_clid)) {
      fallback_depth <- get_depth_to_root(fallback_clid, cl_graph)
      generic_fallback <- is_generic_label(fallback_clid, local_lookup_by_clid(fallback_clid, cl_cfg) %||% "") ||
        (is.finite(fallback_depth) && fallback_depth <= 2)
      if (generic_fallback) {
        specific_pool <- Filter(function(x) {
          if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
          if (!is_descendant_of(x$clid, fallback_clid, cl_graph)) return(FALSE)
          if ((x$map_quality %||% 0) < 2) return(FALSE)
      if (is_generic_label(x$clid, x$label %||% "")) return(FALSE)
          if (is_developmental_stage(clid = x$clid, label = x$label %||% "", cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)) return(FALSE)
          TRUE
        }, cand_ordered_gate)
        if (length(specific_pool) > 0) {
          ord_sp <- order(
            vapply(specific_pool, function(x) x$method_count %||% 0, numeric(1)),
            vapply(specific_pool, function(x) x$score_adj %||% -Inf, numeric(1)),
            vapply(specific_pool, function(x) x$evidence_support %||% 0, numeric(1)),
            decreasing = TRUE
          )
          chosen <- specific_pool[[ord_sp[[1]]]]
          chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
          anchor_fallback_happened <- FALSE
          final_rule <- "parity_avoid_generic_fallback"
        } else {
          # Relaxed rescue path: if strict branch gate filtered everything out,
          # still avoid collapsing to generic fallback by taking the best
          # non-generic descendant from pre-gate ordered candidates.
          specific_pool_relaxed <- Filter(function(x) {
            if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
            if (!is_descendant_of(x$clid, fallback_clid, cl_graph)) return(FALSE)
            if ((x$map_quality %||% 0) < 1) return(FALSE)
            if (is_generic_label(x$clid, x$label %||% "")) return(FALSE)
            if (is_developmental_stage(clid = x$clid, label = x$label %||% "", cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)) return(FALSE)
            TRUE
          }, cand_ordered)
          if (length(specific_pool_relaxed) > 0) {
            ord_relaxed <- order(
              vapply(specific_pool_relaxed, function(x) x$method_count %||% 0, numeric(1)),
              vapply(specific_pool_relaxed, function(x) x$score_adj %||% -Inf, numeric(1)),
              vapply(specific_pool_relaxed, function(x) x$evidence_support %||% 0, numeric(1)),
              vapply(specific_pool_relaxed, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
              decreasing = TRUE
            )
            chosen <- specific_pool_relaxed[[ord_relaxed[[1]]]]
            chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
            anchor_fallback_happened <- FALSE
            final_rule <- "parity_relaxed_descendant_rescue"
          }
        }
      }

      # Additional parity guardrail for lineage-level fallbacks:
      # even if fallback is not generic, prefer a significantly deeper
      # same-lineage descendant when support is comparable.
      if (is.null(chosen) && is.finite(fallback_depth)) {
        lineage_pool <- Filter(function(x) {
          if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
          if (!is_descendant_of(x$clid, fallback_clid, cl_graph)) return(FALSE)
          if ((x$map_quality %||% 0) < 2) return(FALSE)
          if (is_generic_label(x$clid, x$label %||% "")) return(FALSE)
          if (is_developmental_stage(clid = x$clid, label = x$label %||% "", cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)) return(FALSE)
          x_depth <- get_depth_to_root(x$clid, cl_graph)
          if (!is.finite(x_depth) || x_depth <= fallback_depth) return(FALSE)
          TRUE
        }, cand_ordered)

        if (length(lineage_pool) > 0) {
          # Keep only descendants with comparable support / score, then pick most specific.
          fallback_score <- suppressWarnings(max(vapply(Filter(function(x) {
            !is.na(x$clid) && nzchar(x$clid) && identical(x$clid, fallback_clid)
          }, cand_ordered), function(x) as.numeric(x$score_adj %||% -Inf), numeric(1)), na.rm = TRUE))

          lineage_pool <- Filter(function(x) {
            x_score <- as.numeric(x$score_adj %||% -Inf)
            x_support <- as.numeric(x$support_n %||% 0)
            score_ok <- !is.finite(fallback_score) || (is.finite(x_score) && x_score >= (fallback_score - 2 * min_margin))
            support_ok <- x_support >= 1
            isTRUE(score_ok) && isTRUE(support_ok)
          }, lineage_pool)

          if (length(lineage_pool) > 0) {
            ord_lineage <- order(
              vapply(lineage_pool, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
              vapply(lineage_pool, function(x) x$method_count %||% 0, numeric(1)),
              vapply(lineage_pool, function(x) x$score_adj %||% -Inf, numeric(1)),
              decreasing = TRUE
            )
            chosen <- lineage_pool[[ord_lineage[[1]]]]
            chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
            anchor_fallback_happened <- FALSE
            final_rule <- "parity_lineage_descendant_preference"
          }
        }
      }
    }
    # Before fallback, prefer the highest-scoring supported input candidate when available.
    # This prevents a generic fallback from replacing a supported reviewer candidate.
    if (is.null(chosen)) {
      # Filter to only real input methods (not expanded from ontology)
      real_method_candidates <- Filter(function(x) {
        method_name <- x$method %||% ""
        method_name %in% real_input_methods && !is.na(x$clid) && nzchar(x$clid)
      }, cand_ordered)

      # Filter to candidates that are NOT only from unsupported methods
      # (i.e., at least one supporting method)
      supported_real_candidates <- Filter(function(x) {
        !is_unsupported_only(x$clid)
      }, real_method_candidates)

      if (length(supported_real_candidates) > 0) {
        # Pick the highest scoring supported input candidate
        scores <- vapply(supported_real_candidates, function(x) {
          s <- as.numeric(x$score_adj %||% NA_real_)
          if (is.finite(s)) s else -Inf
        }, numeric(1))
        best_idx <- which.max(scores)
        if (length(best_idx) == 0L || !is.finite(scores[[best_idx]])) {
          # deterministic fallback: keep candidate order priority when scores are unavailable
          best_idx <- 1L
        }
        chosen <- supported_real_candidates[[best_idx]]
        chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
        anchor_fallback_happened <- FALSE
        final_rule <- "supported_input_candidate"
      } else {
        # No supported input candidates - use existing fallback logic
        anchor_pool <- Filter(function(x) x$clid == fallback_clid, cand_ordered)
        if (length(anchor_pool) > 0) {
          chosen <- anchor_pool[[1]]
          chosen_support <- list(pass = TRUE, reasons = c("anchor_fallback"), margin = NA_real_)
        } else if (isTRUE(enable_aggressive_lca_hard_guard) &&
                   identical(mode_selected, "aggressive") &&
                   isTRUE(support_gate_pass) &&
                   !is.na(supported_lca_clid) && nzchar(supported_lca_clid) &&
                   !is.na(fallback_clid) && nzchar(fallback_clid) &&
                   identical(supported_lca_clid, fallback_clid)) {
          chosen <- list(
            label = local_lookup_by_clid(supported_lca_clid, cl_cfg) %||% NA_character_,
            clid = supported_lca_clid
          )
          chosen_support <- list(pass = TRUE, reasons = c("aggressive_supported_lca_hard_guard"), margin = NA_real_)
          anchor_fallback_happened <- FALSE
          aggressive_lca_hard_guard_applied <- TRUE
          final_rule <- "aggressive_supported_lca_hard_guard"
        } else {
          chosen <- list(label = local_lookup_by_clid(fallback_clid, cl_cfg) %||% NA_character_, clid = fallback_clid)
          chosen_support <- list(pass = TRUE, reasons = c("anchor_fallback"), margin = NA_real_)
          anchor_fallback_happened <- TRUE
          if (identical(mode_selected, "balanced")) final_rule <- "balanced_anchor_fallback"
        }
      }
    }
  }

  if (!is.null(chosen) && identical(final_rule, "balanced_anchor_fallback") && !is.na(chosen$clid) && nzchar(chosen$clid)) {
    lock_present <- !is.na(lock_id) && nzchar(lock_id)
    chosen_depth <- get_depth_to_root(chosen$clid, cl_graph)
    chosen_generic <- is_generic_label(chosen$clid, local_lookup_by_clid(chosen$clid, cl_cfg) %||% "") ||
      (is.finite(chosen_depth) && chosen_depth <= 2)

    stage_fn <- function(clid, label) {
      is_developmental_stage(
        clid = clid,
        label = label,
        cl_graph = cl_graph,
        cl_cfg = cl_cfg,
        stage_root_clids = stage_root_clids,
        allow_token_fallback = allow_stage_token_fallback
      )
    }

    pools_rule1 <- build_candidate_pools(
      cand_ordered = cand_ordered,
      gate_decisions = gate_decisions,
      chosen = chosen,
      cl_graph = cl_graph,
      cfg_pool = policy_cfg$pool,
      require_deeper = FALSE,
      is_generic_fn = is_generic_label,
      is_stage_fn = stage_fn
    )
    effective_pool <- c(pools_rule1$soft_pool %||% list(), pools_rule1$shadow_pool %||% list())
    instability_pool_n <- as.integer(length(effective_pool))

    instability <- if (isTRUE(lock_present)) compute_lock_instability(effective_pool, policy_cfg$trigger) else {
      list(unsteady = FALSE, reasons = c("lock_missing"), margin = NA_real_, method_gap = NA_real_, ev_gap = NA_real_)
    }
    lock_unsteady <- isTRUE(instability$unsteady)
    lock_unsteady_reasons <- paste(instability$reasons %||% character(0), collapse = "|")

    non_generic_probe <- list(ok = FALSE, reason = "not_checked", checks = list())
    if (!isTRUE(chosen_generic) && isTRUE(lock_present) && isTRUE(lock_unsteady)) {
      non_generic_probe <- should_enter_non_generic_rescue(effective_pool, chosen, policy_cfg, cl_graph)
    }
    non_generic_rescue_ok <- isTRUE(non_generic_probe$ok)

    # Special case: depth_gap == 1 with supported_lca being deeper - prefer specificity
    depth_gap_1_override <- FALSE
    if (isTRUE(lock_present) && isTRUE(lock_unsteady) && 
        is.finite(depth_gap) && depth_gap == 1 &&
        !is.na(supported_lca_clid) && nzchar(supported_lca_clid) &&
        is.finite(supported_lca_depth) && is.finite(chosen_depth) &&
        supported_lca_depth > chosen_depth) {
      # Check if supported_lca is a descendant of chosen (more specific)
      if (is_descendant_of(supported_lca_clid, chosen$clid, cl_graph)) {
        # Prefer the more specific supported_lca over the generic lock
        chosen <- list(
          label = local_lookup_by_clid(supported_lca_clid, cl_cfg) %||% NA_character_, 
          clid = supported_lca_clid
        )
        chosen_support <- list(pass = TRUE, reasons = c("depth_gap_1_refinement"), margin = NA_real_)
        anchor_fallback_happened <- FALSE
        final_rule <- "balanced_depth_gap_1_refinement"
        post_fallback_override_applied <- TRUE
        post_fallback_override_reason <- "depth_gap_1_refinement"
        depth_gap_1_override <- TRUE
      }
    }
    
    if (!isTRUE(depth_gap_1_override)) {
      if (!isTRUE(chosen_generic) && !isTRUE(non_generic_rescue_ok)) {
        post_fallback_why_not_override <- "not_generic"
      } else if (!isTRUE(lock_present)) {
        post_fallback_why_not_override <- "lock_missing"
      } else if (!isTRUE(lock_unsteady)) {
        post_fallback_why_not_override <- "generic_but_lock_stable"
      }
    }

    rescue_enter <- isTRUE(should_enter_rescue(chosen_generic, instability)) || isTRUE(non_generic_rescue_ok)
    if (isTRUE(rescue_enter) && nzchar(lock_id %||% "")) {
      state <- list(
        chosen = chosen,
        chosen_support = chosen_support,
        final_rule = final_rule,
        anchor_fallback_happened = anchor_fallback_happened,
        base = chosen
      )

      rule1_result <- apply_post_fallback(state, pools_rule1, policy_cfg, "non_generic", cl_graph = cl_graph)

      post_fallback_rule1_gate_n <- as.integer(pools_rule1$meta$gate_pool_n %||% 0L)
      post_fallback_rule1_soft_n <- as.integer(pools_rule1$meta$soft_pool_n %||% 0L)
      gate_pool_n <- post_fallback_rule1_gate_n
      soft_pool_n <- post_fallback_rule1_soft_n
      shadow_pool_n <- as.integer(pools_rule1$meta$shadow_pool_n %||% 0L)
      gate_reject_total_n <- as.integer(pools_rule1$meta$gate_reject_total_n %||% 0L)
      gate_reject_descendant_n <- as.integer(pools_rule1$meta$gate_reject_descendant_n %||% 0L)
      gate_reject_descendant_by_stage <- pools_rule1$meta$gate_reject_descendant_by_stage %||% list()
      shadow_kept_n <- as.integer(pools_rule1$meta$shadow_kept_n %||% 0L)
      shadow_stage_flagged_n <- as.integer(pools_rule1$meta$shadow_stage_flagged_n %||% 0L)
      shadow_empty_reason <- as.character(pools_rule1$meta$shadow_empty_reason %||% "no_rejects")
      shadow_desc_filter_counts <- pools_rule1$meta$shadow_desc_filter_counts %||% shadow_desc_filter_counts
      shadow_desc_ref <- pools_rule1$meta$shadow_desc_ref %||% shadow_desc_ref
      shadow_pool_used <- isTRUE(rule1_result$logs$shadow_pool_used)
      override_checks <- rule1_result$logs$override_checks
      score_na_case <- as.character(override_checks$score_na_case %||% "none")
      score_ok_reason <- as.character(override_checks$score_ok_reason %||% "")
      consistency_ok <- isTRUE(override_checks$consistency_ok %||% TRUE)
      consistency_ok_reason <- as.character(override_checks$consistency_ok_reason %||% "")
      score_source_base <- as.character(override_checks$score_source_base %||% "")
      score_source_new <- as.character(override_checks$score_source_new %||% "")
      effective_base_score <- as.numeric(override_checks$effective_base_score %||% NA_real_)
      effective_new_score <- as.numeric(override_checks$effective_new_score %||% NA_real_)
      base_depth <- as.numeric(override_checks$base_depth %||% NA_real_)
      candidate_depth <- as.numeric(override_checks$candidate_depth %||% NA_real_)
      base_subtree_size <- as.integer(override_checks$base_subtree_size %||% 0L)
      base_is_ultra_generic <- isTRUE(override_checks$base_is_ultra_generic %||% FALSE)
      pool_size_1_strong_margin_used <- isTRUE(override_checks$pool_size_1_strong_margin_used %||% FALSE)

      if (isTRUE(rule1_result$logs$override_applied)) {
        chosen <- rule1_result$updated_state$chosen
        chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
        final_rule <- rule1_result$updated_state$final_rule
        anchor_fallback_happened <- FALSE
        post_fallback_override_applied <- TRUE
        post_fallback_override_reason <- "rule1"
        post_fallback_why_not_override <- "overridden"
      } else {
        post_fallback_why_not_override <- rule1_result$logs$why_not_override
        if (!isTRUE(chosen_generic) && isTRUE(non_generic_rescue_ok) && identical(post_fallback_why_not_override, "not_generic")) {
          post_fallback_why_not_override <- "non_generic_rescue_attempt_failed"
        }

      pools_rule2 <- build_candidate_pools(
        cand_ordered = cand_ordered,
        gate_decisions = gate_decisions,
        chosen = chosen,
        cl_graph = cl_graph,
        cfg_pool = policy_cfg$pool,
          require_deeper = TRUE,
          is_generic_fn = is_generic_label,
          is_stage_fn = stage_fn
        )
        rule2_result <- apply_post_fallback(state, pools_rule2, policy_cfg, "descendant", cl_graph = cl_graph)

        post_fallback_rule2_gate_n <- as.integer(pools_rule2$meta$gate_pool_n %||% 0L)
        post_fallback_rule2_soft_n <- as.integer(pools_rule2$meta$soft_pool_n %||% 0L)
        gate_pool_n <- gate_pool_n + post_fallback_rule2_gate_n
        soft_pool_n <- soft_pool_n + post_fallback_rule2_soft_n
        shadow_pool_n <- shadow_pool_n + as.integer(pools_rule2$meta$shadow_pool_n %||% 0L)
        gate_reject_total_n <- gate_reject_total_n + as.integer(pools_rule2$meta$gate_reject_total_n %||% 0L)
        gate_reject_descendant_n <- gate_reject_descendant_n + as.integer(pools_rule2$meta$gate_reject_descendant_n %||% 0L)
        stage_r2 <- pools_rule2$meta$gate_reject_descendant_by_stage %||% list()
        if (length(stage_r2) > 0) {
          for (k in names(stage_r2)) {
            gate_reject_descendant_by_stage[[k]] <- (gate_reject_descendant_by_stage[[k]] %||% 0L) + as.integer(stage_r2[[k]] %||% 0L)
          }
        }
        shadow_kept_n <- shadow_kept_n + as.integer(pools_rule2$meta$shadow_kept_n %||% 0L)
        shadow_stage_flagged_n <- shadow_stage_flagged_n + as.integer(pools_rule2$meta$shadow_stage_flagged_n %||% 0L)
        f2 <- pools_rule2$meta$shadow_desc_filter_counts %||% list()
        if (length(f2) > 0) {
          for (k in names(f2)) {
            shadow_desc_filter_counts[[k]] <- as.integer(shadow_desc_filter_counts[[k]] %||% 0L) + as.integer(f2[[k]] %||% 0L)
          }
        }
        r2_ref <- pools_rule2$meta$shadow_desc_ref %||% list(attempted_n = 0L, kept_n = 0L)
        shadow_desc_ref <- list(
          attempted_n = as.integer(shadow_desc_ref$attempted_n %||% 0L) + as.integer(r2_ref$attempted_n %||% 0L),
          kept_n = as.integer(shadow_desc_ref$kept_n %||% 0L) + as.integer(r2_ref$kept_n %||% 0L)
        )
        if (shadow_kept_n > 0) {
          shadow_empty_reason <- "not_empty"
        } else {
          r2_empty <- as.character(pools_rule2$meta$shadow_empty_reason %||% "")
          if (nzchar(r2_empty) && !identical(r2_empty, "not_empty")) shadow_empty_reason <- r2_empty
        }
        shadow_pool_used <- isTRUE(shadow_pool_used) || isTRUE(rule2_result$logs$shadow_pool_used)
        override_checks <- rule2_result$logs$override_checks
        score_na_case <- as.character(override_checks$score_na_case %||% "none")
        score_ok_reason <- as.character(override_checks$score_ok_reason %||% "")
        consistency_ok <- isTRUE(override_checks$consistency_ok %||% TRUE)
        consistency_ok_reason <- as.character(override_checks$consistency_ok_reason %||% "")
        score_source_base <- as.character(override_checks$score_source_base %||% "")
        score_source_new <- as.character(override_checks$score_source_new %||% "")
        effective_base_score <- as.numeric(override_checks$effective_base_score %||% NA_real_)
        effective_new_score <- as.numeric(override_checks$effective_new_score %||% NA_real_)
        base_depth <- as.numeric(override_checks$base_depth %||% NA_real_)
        candidate_depth <- as.numeric(override_checks$candidate_depth %||% NA_real_)
        base_subtree_size <- as.integer(override_checks$base_subtree_size %||% 0L)
        base_is_ultra_generic <- isTRUE(override_checks$base_is_ultra_generic %||% FALSE)
        pool_size_1_strong_margin_used <- isTRUE(override_checks$pool_size_1_strong_margin_used %||% FALSE)

        if (isTRUE(rule2_result$logs$override_applied)) {
          chosen <- rule2_result$updated_state$chosen
          chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
          final_rule <- rule2_result$updated_state$final_rule
          anchor_fallback_happened <- FALSE
          post_fallback_override_applied <- TRUE
          post_fallback_override_reason <- "rule2"
          post_fallback_why_not_override <- "overridden"
        } else {
          post_fallback_why_not_override <- rule2_result$logs$why_not_override
        }
      }
    } else if (isTRUE(chosen_generic) && isTRUE(lock_present) && !isTRUE(lock_unsteady) && !nzchar(post_fallback_why_not_override)) {
      post_fallback_why_not_override <- "lock_unsteady_false"
    } else if (!isTRUE(chosen_generic) && isTRUE(lock_present) && isTRUE(lock_unsteady) && isTRUE(non_generic_rescue_ok) && !nzchar(post_fallback_why_not_override)) {
      post_fallback_why_not_override <- "non_generic_rescue_candidate_found"
    }
  }

  # Aggressive refinement: when fallback lands on a broad supported LCA, prefer
  # a deeper descendant in the same supported lineage if evidence remains comparable.
  if (!is.null(chosen) && identical(mode_selected, "aggressive") && !is.na(chosen$clid) && nzchar(chosen$clid)) {
    chosen_depth <- get_depth_to_root(chosen$clid, cl_graph)
    chosen_score_adj <- as.numeric(chosen$score_adj %||% NA_real_)
    chosen_evidence <- as.numeric(chosen$evidence_support %||% 0)
    refine_pool <- Filter(function(x) {
      if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
      if (!is_descendant_of(x$clid, chosen$clid, cl_graph)) return(FALSE)
      if ((x$map_quality %||% 0) < 2) return(FALSE)
      x_depth <- get_depth_to_root(x$clid, cl_graph)
      if (!is.finite(x_depth) || x_depth <= chosen_depth) return(FALSE)
      x_score <- as.numeric(x$score_adj %||% NA_real_)
      if (is.finite(chosen_score_adj) && is.finite(x_score) && x_score < (chosen_score_adj - 2 * min_margin)) return(FALSE)
      x_ev <- as.numeric(x$evidence_support %||% 0)
      if (x_ev < (chosen_evidence - 1)) return(FALSE)
      TRUE
    }, cand_ordered_gate)

    if (length(refine_pool) > 0) {
      ord_ref <- order(
        vapply(refine_pool, function(x) !is_developmental_stage(clid = x$clid, label = x$label %||% "", cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback), logical(1)),
        vapply(refine_pool, function(x) x$method_count %||% 0, numeric(1)),
        vapply(refine_pool, function(x) x$evidence_support %||% 0, numeric(1)),
        vapply(refine_pool, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
        vapply(refine_pool, function(x) x$score_adj %||% -Inf, numeric(1)),
        decreasing = TRUE
      )
      refined <- refine_pool[[ord_ref[[1]]]]
      chosen <- refined
      chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
      final_rule <- "aggressive_descendant_refine"
    }
  }

  # Final hard override: if strong anchors exist, force choice within anchor branch
  if (length(strong_anchor_clids) > 0) {
    in_strong_branch <- function(clid) {
      if (is.na(clid) || !nzchar(clid)) return(FALSE)
      dists <- vapply(strong_anchor_clids, function(a) get_ontology_distance(clid, a, cl_graph), numeric(1))
      dists <- dists[is.finite(dists)]
      length(dists) > 0 && min(dists) <= 1
    }
    if (is.null(chosen) || !in_strong_branch(chosen$clid)) {
      strong_pool <- Filter(function(x) {
        if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
        if (!in_strong_branch(x$clid)) return(FALSE)
        if (!evaluate_branch_consistency_gate(x$clid, anchor_clid, cl_graph, anchor_support, anchor_strong_threshold, branch_gate_slack)) return(FALSE)
        TRUE
      }, cand_table_all)
      if (length(strong_pool) > 0) {
        ords <- order(vapply(strong_pool, function(x) x$evidence_support %||% -Inf, numeric(1)),
                      vapply(strong_pool, function(x) x$score_adj %||% -Inf, numeric(1)),
                      vapply(strong_pool, function(x) x$method_count %||% 0, numeric(1)),
                      decreasing = TRUE)
        chosen <- strong_pool[[ords[[1]]]]
        chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
      }
    }
  }

  if (!is.null(chosen) && length(cand_ordered) >= 2) {
    c1 <- cand_ordered[[1]]
    c2 <- cand_ordered[[2]]
    if (!is.na(c1$map_quality) && !is.na(c2$map_quality) && c1$map_quality >= 2 && c2$map_quality >= 2) {
      if (isTRUE(c1$method_count == c2$method_count) && isTRUE(c1$map_quality_rank == c2$map_quality_rank)) {
        if (!is.na(c1$score_adj) && !is.na(c2$score_adj) && abs(c1$score_adj - c2$score_adj) <= specificity_eps) {
          if (is_ancestor_of(c1$clid, c2$clid, cl_graph)) {
            chosen <- c2
            chosen_support <- evaluate_candidate_support(c2, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
            specificity_tiebreak <- TRUE
          } else if (is_ancestor_of(c2$clid, c1$clid, cl_graph)) {
            chosen <- c1
      chosen_support <- evaluate_candidate_support(c1, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
            specificity_tiebreak <- TRUE
          }
        }
      }
    }
  }

  if (is.null(chosen) && length(cand_ordered_gate) > 0) {
    chosen <- cand_ordered_gate[[1]]
      chosen_support <- evaluate_candidate_support(chosen, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
    low_quality <- all(vapply(cand_ordered_gate, function(x) (x$map_quality %||% 0) <= 1, logical(1)))
    weak_scores <- all(vapply(cand_ordered_gate, function(x) is.na(x$score) || (x$score < strong_score), logical(1)))
    cross_conflict_any <- isTRUE(chosen_support$cross_conflict %||% FALSE)
    if (low_quality && weak_scores && cross_conflict_any) {
      chosen_support$reasons <- unique(c(chosen_support$reasons, "no_candidate_passed_gate"))
    }
  }

  # Descendant override is disabled.
  # if (!is.null(chosen) && !is.na(lock_id) && nzchar(lock_id)) {
  if (FALSE) {
    if (!is.na(chosen$clid) && nzchar(chosen$clid) && identical(chosen$clid, lock_id)) {
      subtype_candidates <- Filter(function(x) {
        if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
        if (!is_descendant_of(x$clid, lock_id, cl_graph)) return(FALSE)
        if (is_developmental_stage(clid = x$clid, label = x$label, cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)) return(FALSE)
        if (is_unsupported_only(x$clid)) return(FALSE)
        TRUE
      }, cand_ordered_gate)

      # Relaxed pool: if gate is too strict, still allow descendant override from
      # pre-gate candidates with basic quality.
      if (length(subtype_candidates) == 0) {
        subtype_candidates <- Filter(function(x) {
          if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
          if (!is_descendant_of(x$clid, lock_id, cl_graph)) return(FALSE)
          if ((x$map_quality %||% 0) < 1) return(FALSE)
          if (is_generic_label(x$clid, x$label %||% "")) return(FALSE)
          if (is_developmental_stage(clid = x$clid, label = x$label, cl_graph = cl_graph, cl_cfg = cl_cfg, stage_root_clids = stage_root_clids, allow_token_fallback = allow_stage_token_fallback)) return(FALSE)
          TRUE
        }, cand_ordered)
      }

      if (length(subtype_candidates) > 0) {
        subtype_real <- Filter(function(x) (x$method_count %||% 0) >= min_votes, subtype_candidates)
        subtype_pool <- if (length(subtype_real) > 0) subtype_real else {
          Filter(function(x) (x$method_count_all %||% 0) >= min_votes, subtype_candidates)
        }
        if (length(subtype_pool) == 0) subtype_pool <- subtype_candidates

        ord_sub <- order(
          vapply(subtype_pool, function(x) x$score_adj %||% -Inf, numeric(1)),
          vapply(subtype_pool, function(x) x$map_quality_rank %||% 0, numeric(1)),
          vapply(subtype_pool, function(x) x$specificity, numeric(1)),
          decreasing = TRUE
        )
        best_sub <- subtype_pool[[ord_sub[[1]]]]
        score_ok <- is.na(chosen$score_adj) || is.na(best_sub$score_adj) || best_sub$score_adj >= (chosen$score_adj - 2 * min_margin)
        if (isTRUE(score_ok) && (best_sub$map_quality %||% 0) >= 1) {
          support_sub <- evaluate_candidate_support(best_sub, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
          # Accept when gate passes OR when descendant has tangible support and is more specific.
          chosen_depth <- get_depth_to_root(chosen$clid, cl_graph)
          sub_depth <- get_depth_to_root(best_sub$clid, cl_graph)
          best_score_adj <- as.numeric(best_sub$score_adj %||% NA_real_)
          best_ev <- as.numeric(best_sub$evidence_support %||% 0)
          soft_support <- ((best_sub$method_count %||% 0) >= 1) || (best_ev >= 1) || is.finite(best_score_adj)
          # Allow soft-rejected descendants to enter the override check.
          # only if they are non-generic, more specific, and near ancestor score.
          soft_near <- !is.finite(chosen$score_adj) || !is.finite(best_score_adj) || (best_score_adj >= (chosen$score_adj - 2.5 * min_margin))
    relaxed_accept <- soft_support && isTRUE(soft_near) && is.finite(sub_depth) && (!is.finite(chosen_depth) || sub_depth > chosen_depth) && !is_generic_label(best_sub$clid, best_sub$label %||% "")
          if (support_sub$pass || is.na(chosen$score_adj) || relaxed_accept) {
            chosen <- best_sub
            chosen_support <- support_sub
            chosen_support$reasons <- unique(c(chosen_support$reasons, "iter5_descendant_override"))
            final_rule <- "iter5_descendant_override"
          }
        }
      }
    }
  }

  score_missing_flag <- FALSE
  if (length(cand_ordered) >= 2) {
    c1 <- cand_ordered[[1]]
    c2 <- cand_ordered[[2]]
    if (isTRUE(c1$method_count == c2$method_count) && isTRUE(c1$map_quality_rank == c2$map_quality_rank)) {
      if ((is.na(c1$score) && !is.na(c2$score)) || (!is.na(c1$score) && is.na(c2$score))) {
        score_missing_flag <- TRUE
      }
    }
  }

  if (is.null(chosen)) return(head_out)

  chosen_clid <- chosen$clid
  chosen_label <- chosen$label
  if (is.na(chosen_label) || !nzchar(chosen_label)) {
    chosen_label <- local_lookup_by_clid(chosen_clid, cl_cfg) %||% head_out$final_decision$primary_cell_type
  }

  set_choice <- function(clid, label = NULL, rule_name = NULL) {
    chosen_clid <<- as.character(clid %||% NA_character_)
    if (is.null(label) || is.na(label) || !nzchar(as.character(label))) {
      chosen_label <<- as.character(local_lookup_by_clid(chosen_clid, cl_cfg) %||% chosen_label)
    } else {
      chosen_label <<- as.character(label)
    }
    if (!is.null(rule_name) && nzchar(as.character(rule_name))) final_rule <<- as.character(rule_name)
    invisible(NULL)
  }

  # Decision-category metadata parsing (does not change the biological label).
  decision_cat <- normalize_decision_category(head_out$final_decision$decision_category %||% "")
  preferred_method <- NA_character_
  if (decision_cat == "cassia_better") preferred_method <- "cassia"
  if (decision_cat == "in_house_better") preferred_method <- "in_house"
  if (decision_cat == "enrich_better") preferred_method <- "enrich"
  if (decision_cat == "third_party_override") preferred_method <- "third_party"
  decision_category_v2 <- if (!is.na(preferred_method) && nzchar(preferred_method)) "method_better" else "tie_or_override"

  method_key_aliases <- function(method_key) {
    mk <- as.character(method_key %||% "")
    if (mk == "our") return(c("our", "in_house", "our_method"))
    if (mk == "inter") mk <- "enrich"
    if (mk == "enrich") return(c("enrich"))
    c(mk)
  }

  pick_top1_by_method <- function(method_key, require_gate_pass = TRUE) {
    src <- if (isTRUE(require_gate_pass)) cand_ordered_gate else cand_ordered
    aliases <- method_key_aliases(method_key)
    cands <- Filter(function(x) {
      as.character(x$method %||% "") %in% aliases &&
        as.integer(x$rank %||% 99L) == 1L &&
        !is.na(x$clid) && nzchar(x$clid)
    }, src)
    if (length(cands) == 0) return(NULL)

    ord <- order(
      vapply(cands, function(x) as.numeric(x$score_adj %||% -Inf), numeric(1)),
      vapply(cands, function(x) as.numeric(x$method_count %||% 0), numeric(1)),
      vapply(cands, function(x) as.numeric(x$evidence_support %||% 0), numeric(1)),
      vapply(cands, function(x) as.numeric(x$score %||% -Inf), numeric(1)),
      vapply(cands, function(x) as.numeric(x$specificity %||% 0), numeric(1)),
      -vapply(cands, function(x) as.integer(x$rank %||% 99L), integer(1)),
      decreasing = TRUE
    )
    cand <- cands[[ord[[1]]]]
    if (!isTRUE(require_gate_pass)) return(cand)
    support <- evaluate_candidate_support(cand, cand_table, candidates, lock_id, min_margin, strong_score, cl_graph)
    if (!isTRUE(support$pass)) return(NULL)
    cand
  }

  # Frozen top1 lookup from the pre-hard-filter candidate table, for consistency/audit only.
  pick_top1_by_method_all <- function(method_key) {
    aliases <- method_key_aliases(method_key)
    cands <- Filter(function(x) {
      as.character(x$method %||% "") %in% aliases &&
        as.integer(x$rank %||% 99L) == 1L &&
        !is.na(x$clid) && nzchar(x$clid)
    }, cand_table_all)
    if (length(cands) == 0L) return(NULL)
    cands[[1]]
  }

  # decision_category is descriptive metadata and does not change the label.
  # Do not force the biological label to match a possibly inconsistent LLM category.
  # Final category consistency is checked after all biological selection/contrastive gates.


  # Top2 same-lineage promotion (simple, CLID-first):
  # when locked on a broad consensus ancestor, allow deeper rank>1 candidates
  # in the same lineage if score is close enough.
  if (!is.na(chosen_clid) && nzchar(chosen_clid) && !is.na(lock_id) && nzchar(lock_id) &&
      identical(chosen_clid, lock_id) && identical(lock_reason, "consensus_msca") && isTRUE(support_gate_pass)) {
    chosen_depth <- get_depth_to_root(chosen_clid, cl_graph)
    chosen_score_vals <- vapply(Filter(function(x) !is.na(x$clid) && identical(x$clid, chosen_clid), cand_ordered), function(x) as.numeric(x$score %||% NA_real_), numeric(1))
    chosen_score_vals <- chosen_score_vals[is.finite(chosen_score_vals)]
    chosen_score_raw <- if (length(chosen_score_vals) > 0) max(chosen_score_vals) else NA_real_
    chosen_top_hits <- suppressWarnings(max(vapply(Filter(function(x) !is.na(x$clid) && identical(x$clid, chosen_clid), cand_ordered), function(x) as.numeric(x$evidence_top_hits %||% 0), numeric(1)), na.rm = TRUE))
    if (!is.finite(chosen_top_hits)) chosen_top_hits <- 0

    promo_pool <- Filter(function(x) {
      if (is.na(x$clid) || !nzchar(x$clid)) return(FALSE)
      if (!isTRUE(is_descendant_of(x$clid, lock_id, cl_graph))) return(FALSE)
      if ((x$rank %||% 1L) <= 1L) return(FALSE)
      # Do not promote synthetic expanded descendants.
      if (identical(as.character(x$method %||% ""), "expanded")) return(FALSE)

      # Require a real CLID signal for promotion; avoid label-only coerced promotion.
      raw_clid <- as.character(x$raw_clid %||% "")
      coerce_status <- as.character(x$coerce_status %||% "")
      if (!nzchar(raw_clid) && coerce_status %in% c("mapped_exact", "mapped_synonym", "mapped_token")) return(FALSE)

      x_depth <- get_depth_to_root(x$clid, cl_graph)
      if (!is.finite(x_depth) || x_depth <= chosen_depth) return(FALSE)
      x_score_raw <- as.numeric(x$score %||% NA_real_)
      x_top_hits <- as.numeric(x$evidence_top_hits %||% 0)
      if (is.finite(chosen_score_raw) && is.finite(x_score_raw) && x_score_raw < (chosen_score_raw - min_margin)) return(FALSE)
      if (is.finite(chosen_top_hits) && chosen_top_hits > 0 && is.finite(x_top_hits) && x_top_hits < (0.6 * chosen_top_hits)) return(FALSE)
      TRUE
    }, cand_ordered_gate)

    if (length(promo_pool) > 0) {
      ord <- order(
        vapply(promo_pool, function(x) as.numeric(x$score %||% -Inf), numeric(1)),
        vapply(promo_pool, function(x) as.numeric(x$evidence_top_hits %||% 0), numeric(1)),
        vapply(promo_pool, function(x) get_depth_to_root(x$clid, cl_graph), numeric(1)),
        vapply(promo_pool, function(x) as.integer(x$rank %||% 99L), integer(1)),
        decreasing = TRUE
      )
      pick <- promo_pool[[ord[[1]]]]
      set_choice(pick$clid, pick$label, "same_lineage_top2_promotion")
      guardrail_trace <- list(
        triggered = TRUE,
        gap = NA_real_,
        gap_ratio = NA_real_,
        best_method = pick$method %||% NA_character_,
        best_score = as.numeric(pick$score %||% NA_real_),
        judge_score = chosen_score_raw
      )
    }
  }

  # ---- Final consensus-override gate ----
  # Evaluate the ACTUAL final candidate, not greedy_target. This catches fallback,
  # post-fallback override, strong-anchor rescue, method-alignment, and top2-promotion paths.
  pre_consensus_final_clid <- if (isTRUE(head_policy_active)) llm_anchor_clid else as.character(chosen_clid %||% NA_character_)
  pre_consensus_final_label <- if (isTRUE(head_policy_active)) llm_anchor_label else as.character(chosen_label %||% NA_character_)
  co_eval <- evaluate_final_consensus_override(pre_consensus_final_clid)

  consensus_override_attempted <- isTRUE(co_eval$attempted)
  consensus_override_allowed <- co_eval$allowed
  consensus_override_reason <- as.character(co_eval$reason %||% "")
  override_candidate_clid <- as.character(co_eval$candidate_clid %||% NA_character_)
  override_relation_to_anchor <- as.character(co_eval$relation %||% "")
  contrastive_shared_score <- suppressWarnings(as.numeric(co_eval$shared_score %||% NA_real_))
  contrastive_candidate_unique_score <- suppressWarnings(as.numeric(co_eval$candidate_unique_score %||% NA_real_))
  contrastive_anchor_unique_score <- suppressWarnings(as.numeric(co_eval$anchor_unique_score %||% NA_real_))
  contrastive_margin <- suppressWarnings(as.numeric(co_eval$margin %||% NA_real_))
  contrastive_candidate_total_score <- suppressWarnings(as.numeric(co_eval$candidate_total_score %||% NA_real_))
  contrastive_anchor_total_score <- suppressWarnings(as.numeric(co_eval$anchor_total_score %||% NA_real_))
  contrastive_candidate_real_methods <- as.character(co_eval$candidate_real_methods %||% character(0))

  if (nzchar(Sys.getenv("DEBUG_OVERRIDE"))) {
    cat(
      "[OVERRIDE-DEBUG] final_candidate=", pre_consensus_final_clid,
      " anchor=", as.character(dc_anchor_clid %||% NA_character_),
      " relation=", override_relation_to_anchor,
      " attempted=", consensus_override_attempted,
      " allowed=", if (is.na(consensus_override_allowed)) "NA" else as.character(consensus_override_allowed),
      " cand_unique=", if (is.na(contrastive_candidate_unique_score)) "NA" else sprintf("%.2f", contrastive_candidate_unique_score),
      " anchor_unique=", if (is.na(contrastive_anchor_unique_score)) "NA" else sprintf("%.2f", contrastive_anchor_unique_score),
      " shared=", if (is.na(contrastive_shared_score)) "NA" else sprintf("%.2f", contrastive_shared_score),
      " reason=", consensus_override_reason,
      "\n", sep = ""
    )
  }

  if (isTRUE(consensus_override_attempted)) {
    if (isTRUE(head_policy_active)) {
      # Protected policies are discrete and threshold-free at the identity level:
      # any non-hierarchical disagreement between a valid Head identity and a direct
      # reviewer majority is sent to automated Chief QC. Contrastive marker scores
      # remain diagnostics only and never rewrite identity.
      head_out <- flag_issue(
        head_out,
        "head_vs_majority_cross_branch_qc",
        paste0("Valid Head identity ", pre_consensus_final_clid,
               " differs non-hierarchically from direct-majority anchor ",
               as.character(dc_anchor_clid %||% NA_character_),
               "; identity preserved and routed to automated Chief QC. ",
               "Contrastive diagnostics: ", consensus_override_reason)
      )
    } else {
      # Ablation path: retain the evidence-qualified rerank behavior.
      if (isTRUE(consensus_override_allowed)) {
        final_rule <- if (identical(final_rule, "balanced_default")) "contrastive_consensus_override" else final_rule
      } else if (identical(consensus_override_allowed, FALSE)) {
        anchor_id <- as.character(dc_anchor_clid %||% NA_character_)
        anchor_pool <- Filter(function(x) {
          !is.na(x$clid) && identical(as.character(x$clid), anchor_id) &&
            as.integer(x$rank %||% 99L) == 1L &&
            as.character(x$method %||% "") %in% real_input_methods &&
            !is_unsupported_only(anchor_id)
        }, cand_table_all)

        if (length(anchor_pool) > 0L) {
          ord_anchor <- order(
            -vapply(anchor_pool, function(x) as.integer(x$rank %||% 99L), integer(1)),
            vapply(anchor_pool, function(x) as.numeric(x$evidence_support %||% 0), numeric(1)),
            vapply(anchor_pool, function(x) as.numeric(x$map_quality %||% 0), numeric(1)),
            decreasing = TRUE
          )
          anchor_choice <- anchor_pool[[ord_anchor[[1]]]]
          set_choice(anchor_id, anchor_choice$label %||% local_lookup_by_clid(anchor_id, cl_cfg),
                     "consensus_anchor_retained_contrastive_gate")
          chosen_support <- evaluate_candidate_support(anchor_choice, cand_table_all, candidates, lock_id, min_margin, strong_score, cl_graph)
          consensus_override_forced_anchor <- TRUE
        } else {
          head_out <- flag_issue(
            head_out,
            "consensus_override_requires_qc",
            paste0("Direct-majority anchor ", anchor_id,
                   " could not safely replace non-hierarchical candidate ", pre_consensus_final_clid,
                   "; automated Chief QC required. ", consensus_override_reason)
          )
        }
      } else {
        head_out <- flag_issue(
          head_out,
          "consensus_override_requires_qc",
          paste0("Non-hierarchical candidate ", pre_consensus_final_clid,
                 " differs from direct-majority anchor ", as.character(dc_anchor_clid %||% NA_character_),
                 "; deterministic contrastive evidence unavailable. ", consensus_override_reason)
        )
      }
    }
  }

  # ---- Ontology-guarded identity resolution ----
  # The Head is the biological anchor. Deterministic post-processing may adjust
  # specificity only when ALL of the following are true:
  #   (1) Head and deterministic proposal are in a direct ancestor/descendant relation;
  #   (2) the deterministic target is a real reviewer candidate (not ontology expansion);
  #   (3) the target is not supported only by reviewers that the Head explicitly rejected.
  # Non-hierarchical changes are prohibited and routed to Chief QC.
  deterministic_pre_preserve_clid <- as.character(chosen_clid %||% NA_character_)
  deterministic_pre_preserve_label <- as.character(chosen_label %||% NA_character_)
  head_identity_restore_applied <- FALSE
  head_identity_candidate_present <- FALSE
  head_identity_preserved_final <- FALSE
  hierarchical_adjustment_applied <- FALSE
  hierarchical_adjustment_candidate_present <- FALSE
  hierarchical_adjustment_candidate_supported <- FALSE
  head_deterministic_relation <- if (isTRUE(head_identity_valid)) "deterministic_unavailable" else "head_invalid"
  identity_resolution_action <- if (isTRUE(head_identity_valid)) "pending" else "deterministic_only_head_invalid"

  if (isTRUE(head_policy_active)) {
    head_pool <- Filter(function(x) {
      !is.na(x$clid) && identical(as.character(x$clid), llm_anchor_clid)
    }, cand_table_all)
    head_identity_candidate_present <- length(head_pool) > 0L

    det_clid <- as.character(deterministic_pre_preserve_clid %||% NA_character_)
    det_valid <- !is.na(det_clid) && nzchar(det_clid) && det_clid %in% names(cl_graph$cl)

    det_real_pool <- if (isTRUE(det_valid)) {
      Filter(function(x) {
        !is.na(x$clid) && identical(as.character(x$clid), det_clid) &&
          as.character(x$method %||% "") %in% real_input_methods &&
          !identical(as.character(x$method %||% ""), "expanded")
      }, cand_table_all)
    } else {
      list()
    }
    hierarchical_adjustment_candidate_present <- length(det_real_pool) > 0L
    hierarchical_adjustment_candidate_supported <- isTRUE(det_valid) &&
      isTRUE(hierarchical_adjustment_candidate_present) &&
      !is_unsupported_only(det_clid)

    restore_head <- function(rule_name, issue_flag = NULL, issue_note = NULL) {
      head_identity_restore_applied <<- !identical(
        as.character(deterministic_pre_preserve_clid %||% NA_character_),
        llm_anchor_clid
      )
      chosen_clid <<- llm_anchor_clid
      chosen_label <<- llm_anchor_label
      head_identity_preserved_final <<- TRUE
      final_rule <<- rule_name
      consensus_override_forced_anchor <<- FALSE
      if (!is.null(issue_flag) && nzchar(issue_flag)) {
        head_out <<- flag_issue(head_out, issue_flag, issue_note %||% issue_flag)
      }
    }

    if (!isTRUE(det_valid)) {
      head_deterministic_relation <- "deterministic_unavailable"
      identity_resolution_action <- "head_preserved_deterministic_unavailable"
      restore_head(
        "head_editor_preserved_deterministic_unavailable",
        "deterministic_identity_unavailable_head_preserved",
        paste0("Deterministic post-processing did not produce a valid frozen CLID; valid Head identity ",
               llm_anchor_clid, " was preserved.")
      )

    } else if (identical(det_clid, llm_anchor_clid)) {
      head_deterministic_relation <- "exact"
      identity_resolution_action <- "head_and_deterministic_agree"
      chosen_clid <- llm_anchor_clid
      chosen_label <- llm_anchor_label
      head_identity_preserved_final <- TRUE
      final_rule <- "head_editor_deterministic_agree"

    } else {
      det_is_descendant <- isTRUE(is_ancestor_of(llm_anchor_clid, det_clid, cl_graph))
      det_is_ancestor <- isTRUE(is_ancestor_of(det_clid, llm_anchor_clid, cl_graph))

      if (isTRUE(det_is_descendant) || isTRUE(det_is_ancestor)) {
        head_deterministic_relation <- if (isTRUE(det_is_descendant)) {
          "deterministic_descendant_of_head"
        } else {
          "deterministic_ancestor_of_head"
        }

        if (isTRUE(force_preserve_head_identity)) {
          identity_resolution_action <- "head_preserve_forced"
          restore_head(
            "head_editor_preserved_over_hierarchical_rerank",
            "deterministic_rerank_disagreed_with_head",
            paste0("Deterministic proposal ", det_clid, " is ontology-hierarchical with Head ",
                   llm_anchor_clid, " but identity_policy=head_preserve forces the Head identity.")
          )

        } else if (isTRUE(ontology_guarded_identity) &&
                   isTRUE(head_identity_candidate_present) &&
                   isTRUE(hierarchical_adjustment_candidate_supported)) {
          # Keep the deterministic proposal only as a within-lineage specificity calibration.
          hierarchical_adjustment_applied <- TRUE
          identity_resolution_action <- if (isTRUE(det_is_descendant)) {
            "hierarchical_refinement_allowed"
          } else {
            "hierarchical_broadening_allowed"
          }
          head_identity_preserved_final <- FALSE
          final_rule <- if (isTRUE(det_is_descendant)) {
            "ontology_guarded_hierarchical_refinement"
          } else {
            "ontology_guarded_hierarchical_broadening"
          }
          head_out <- flag_issue(
            head_out,
            "deterministic_hierarchical_adjustment_from_head",
            paste0("Allowed ontology-hierarchical specificity calibration from Head ",
                   llm_anchor_clid, " to real supported reviewer candidate ", det_clid,
                   "; relation=", head_deterministic_relation, ".")
          )

        } else {
          # The relation is hierarchical, but provenance/support is insufficient to let
          # deterministic heuristics alter the valid Head identity.
          identity_resolution_action <- "hierarchical_adjustment_blocked_preserve_head"
          restore_head(
            "head_editor_preserved_hierarchical_adjustment_blocked",
            "deterministic_hierarchical_adjustment_blocked",
            paste0("Blocked ontology-hierarchical deterministic proposal ", det_clid,
                   " because it was not a supported real reviewer candidate or the Head identity ",
                   "was outside the reviewer candidate pool; preserved Head ", llm_anchor_clid, ".")
          )
        }

      } else {
        head_deterministic_relation <- "non_hierarchical"
        identity_resolution_action <- "head_preserved_cross_branch_qc"
        restore_head(
          "head_editor_preserved_cross_branch_guard",
          "head_vs_deterministic_cross_branch_qc",
          paste0("Deterministic proposal ", det_clid,
                 " is non-hierarchical with valid Head identity ", llm_anchor_clid,
                 "; cross-branch replacement was blocked and routed to automated Chief QC.")
        )
        # Diagnostic marker.
        head_out <- flag_issue(
          head_out,
          "deterministic_rerank_disagreed_with_head",
          paste0("Deterministic candidate selection proposed ", det_clid,
                 " but ontology guard preserved Head ", llm_anchor_clid, ".")
        )
      }
    }

    # If the final identity is the Head identity, retain Head support diagnostics.
    if (identical(as.character(chosen_clid %||% NA_character_), llm_anchor_clid)) {
      if (length(head_pool) > 0L) {
        ord_head <- order(
          -vapply(head_pool, function(x) as.integer(x$rank %||% 99L), integer(1)),
          vapply(head_pool, function(x) as.numeric(x$map_quality %||% 0), numeric(1)),
          decreasing = TRUE
        )
        head_choice <- head_pool[[ord_head[[1]]]]
        chosen_support <- evaluate_candidate_support(
          head_choice, cand_table_all, candidates, lock_id, min_margin, strong_score, cl_graph
        )
      } else {
        chosen_support <- list(
          pass = TRUE,
          reasons = c("head_identity_outside_reviewer_candidate_pool"),
          review_flags = c("head_identity_outside_reviewer_candidate_pool"),
          margin = NA_real_
        )
        head_out <- flag_issue(
          head_out,
          "head_identity_outside_reviewer_candidate_pool",
          paste0("Valid Head identity ", llm_anchor_clid,
                 " was not present in the reviewer candidate pool; identity preserved and sent to automated QC.")
        )
      }
    }
  } else {
    # Under reranking or an invalid Handling Editor output, deterministic selection remains authoritative.
    head_identity_preserved_final <- isTRUE(head_identity_valid) &&
      identical(as.character(chosen_clid %||% NA_character_), llm_anchor_clid)
    if (isTRUE(head_identity_valid) && !is.na(deterministic_pre_preserve_clid) &&
        nzchar(deterministic_pre_preserve_clid)) {
      if (identical(deterministic_pre_preserve_clid, llm_anchor_clid)) {
        head_deterministic_relation <- "exact"
      } else if (isTRUE(is_ancestor_of(llm_anchor_clid, deterministic_pre_preserve_clid, cl_graph))) {
        head_deterministic_relation <- "deterministic_descendant_of_head"
      } else if (isTRUE(is_ancestor_of(deterministic_pre_preserve_clid, llm_anchor_clid, cl_graph))) {
        head_deterministic_relation <- "deterministic_ancestor_of_head"
      } else {
        head_deterministic_relation <- "non_hierarchical"
      }
    }
  }

  # Simple label refinement rule: if multiple labels map to the same chosen CLID,
  # keep the finer textual label (more specific wording) without changing CLID.
  normalize_specific_label <- function(lbl) {
    x <- as.character(lbl %||% "")
    if (!nzchar(x)) return("")
    x <- gsub("\\([^)]*\\)", "", x)
    stringr::str_squish(x)
  }
  if (!is.na(chosen_clid) && nzchar(chosen_clid)) {
    same_clid_labels <- vapply(Filter(function(x) {
      !is.na(x$clid) && nzchar(x$clid) && identical(x$clid, chosen_clid) && !is.na(x$label) && nzchar(x$label)
    }, pool), function(x) as.character(x$label), character(1))
    same_clid_labels <- unique(c(as.character(chosen_label %||% ""), same_clid_labels))
    same_clid_labels <- same_clid_labels[nzchar(same_clid_labels)]
    if (length(same_clid_labels) > 1) {
      norm <- vapply(same_clid_labels, normalize_specific_label, character(1))
      tok_n <- vapply(norm, function(z) if (!nzchar(z)) 0L else length(strsplit(tolower(z), "\\s+")[[1]]), integer(1))
      char_n <- nchar(norm, type = "chars")
      ord <- order(tok_n, char_n, norm, decreasing = TRUE)
      chosen_label <- same_clid_labels[[ord[[1]]]]
      label_tiebreak_trace <- list(
        triggered = TRUE,
        reason = "Selected the most specific label among candidates mapped to the same CLID.",
        candidate_labels = as.list(same_clid_labels),
        chosen_label = chosen_label
      )
      if (identical(final_rule, "balanced_default")) final_rule <- "same_clid_finer_label"
    }
  }

  # Canonicalize final primary label to CL label for stable wording.
  canonical_primary <- local_lookup_by_clid(chosen_clid, cl_cfg)
  if (!is.na(canonical_primary) && nzchar(canonical_primary)) {
    chosen_label <- canonical_primary
    if (identical(final_rule, "balanced_default") || identical(final_rule, "same_clid_finer_label") || identical(final_rule, "atomic_primary_label")) {
      final_rule <- "canonical_primary_label"
    }
  }

  # Enforce atomic primary label (single cell type phrase) within chosen CLID.
  is_multi_label_text <- function(lbl) {
    s <- tolower(stringr::str_squish(as.character(lbl %||% "")))
    if (!nzchar(s)) return(FALSE)
    grepl("\\s*,\\s*|\\s*;\\s*|\\s*\\|\\s*|\\s*/\\s*|\\s+and\\s+", s)
  }
  if (!is.na(chosen_clid) && nzchar(chosen_clid) && is_multi_label_text(chosen_label)) {
    atomic_pool <- Filter(function(x) {
      if (is.na(x$clid) || !nzchar(x$clid) || !identical(x$clid, chosen_clid)) return(FALSE)
      lbl <- as.character(x$label %||% "")
      nzchar(lbl) && !is_multi_label_text(lbl)
    }, pool)
    if (length(atomic_pool) > 0) {
      ord_atomic <- order(
        vapply(atomic_pool, function(x) as.numeric(x$score %||% -Inf), numeric(1)),
        vapply(atomic_pool, function(x) as.numeric(x$evidence_top_hits %||% 0), numeric(1)),
        vapply(atomic_pool, function(x) nchar(as.character(x$label %||% ""), type = "chars"), numeric(1)),
        decreasing = TRUE
      )
      chosen_label <- as.character(atomic_pool[[ord_atomic[[1]]]]$label %||% chosen_label)
      label_tiebreak_trace <- list(
        triggered = TRUE,
        reason = "Replaced multi-label primary with the best atomic label under the same CLID.",
        candidate_labels = as.list(vapply(atomic_pool, function(x) as.character(x$label %||% ""), character(1))),
        chosen_label = chosen_label
      )
      if (identical(final_rule, "balanced_default") || identical(final_rule, "same_clid_finer_label")) {
        final_rule <- "atomic_primary_label"
      }
    } else {
      canonical <- local_lookup_by_clid(chosen_clid, cl_cfg)
      if (!is.na(canonical) && nzchar(canonical)) {
        chosen_label <- canonical
        if (identical(final_rule, "balanced_default") || identical(final_rule, "same_clid_finer_label")) {
          final_rule <- "atomic_primary_label"
        }
      }
    }
  }

  # Normalize and deduplicate secondary signals; drop entries equivalent to primary.
  normalize_signal <- function(x) {
    s <- tolower(stringr::str_squish(as.character(x %||% "")))
    s <- gsub("\\([^)]*\\)", "", s)
    s <- gsub("[^a-z0-9+ -]", "", s, perl = TRUE)
    s <- gsub("\\bcells\\b", "cell", s)
    stringr::str_squish(s)
  }
  sec_raw <- head_out$final_decision$secondary_signals %||% list()
  if (!is.list(sec_raw)) sec_raw <- as.list(sec_raw)
  sec_vals <- as.character(unlist(sec_raw, recursive = TRUE, use.names = FALSE))
  sec_vals <- sec_vals[nzchar(sec_vals)]
  if (length(sec_vals) > 0) {
    primary_norm <- normalize_signal(chosen_label)
    kept <- character(0)
    seen <- character(0)
    for (sv in sec_vals) {
      nsv <- normalize_signal(sv)
      if (!nzchar(nsv)) next
      if (identical(nsv, primary_norm)) next
      if (nsv %in% seen) next
      kept <- c(kept, sv)
      seen <- c(seen, nsv)
    }
    head_out$final_decision$secondary_signals <- as.list(kept)
  } else {
    head_out$final_decision$secondary_signals <- list()
  }

  sibling_conflict <- FALSE
  if (length(cand_ordered) >= 2) {
    c1 <- cand_ordered[[1]]
    c2 <- cand_ordered[[2]]
    p1 <- get_parent_clid(c1$clid, cl_graph)
    p2 <- get_parent_clid(c2$clid, cl_graph)
    if (!is.na(p1) && !is.na(p2) && p1 == p2) {
      if (is.finite(chosen_support$margin) && chosen_support$margin < min_margin) sibling_conflict <- TRUE
    }
  }

  conflict_fallback <- FALSE
  if (length(cand_ordered) >= 2 && !is.na(lock_id) && nzchar(lock_id)) {
    c1 <- cand_ordered[[1]]
    c2 <- cand_ordered[[2]]
    if (!is.na(c1$clid) && !is.na(c2$clid)) {
      related <- is_ancestor_of(c1$clid, c2$clid, cl_graph) || is_ancestor_of(c2$clid, c1$clid, cl_graph)
      if (!related) {
        msca <- get_msca_two(c1$clid, c2$clid, cl_graph)
        if (!is.na(msca$id) && is_descendant_of(msca$id, lock_id, cl_graph)) {
          dist <- get_ontology_distance(c1$clid, c2$clid, cl_graph)
          if (!is.infinite(dist) && dist >= 2) {
            diff_adj <- abs((c1$score_adj %||% 0) - (c2$score_adj %||% 0))
            if (diff_adj < conflict_delta) conflict_fallback <- TRUE
          }
        }
      }
    }
  }

  risk_flags <- character(0)
  review_reasons <- chosen_support$review_flags %||% character(0)
  risk_flags <- unique(c(risk_flags, review_reasons))
  if (!chosen_support$pass) risk_flags <- unique(c(risk_flags, chosen_support$reasons))
  if (is.na(lock_id) || !nzchar(lock_id)) risk_flags <- unique(c(risk_flags, "no_lock"))
  if (is.finite(lock_anc_count) && lock_anc_count <= 1) risk_flags <- unique(c(risk_flags, "high_level_lock"))

  # ontology_id_mismatch based on pairwise ontology distance of top1 CLIDs
  top1_cands <- Filter(function(x) x$rank == 1 && !is.na(x$clid) && nzchar(x$clid), candidates)
  if (length(top1_cands) >= 2) {
    flag_mismatch <- FALSE
    for (i in seq_len(length(top1_cands) - 1)) {
      for (j in (i + 1):length(top1_cands)) {
        a <- top1_cands[[i]]$clid
        b <- top1_cands[[j]]$clid
        if (is.na(a) || is.na(b) || !nzchar(a) || !nzchar(b) || a == b) next
        d_a <- get_depth_from_anc(get_anc_map(a, cl_graph))
        d_b <- get_depth_from_anc(get_anc_map(b, cl_graph))
        msca <- get_msca_two(a, b, cl_graph)
        d_m <- msca$depth
        if (is.na(d_m)) next
        dist <- d_a + d_b - 2 * d_m
        qa <- as.numeric(top1_cands[[i]]$map_quality %||% 0)
        qb <- as.numeric(top1_cands[[j]]$map_quality %||% 0)
        if (qa >= 2 && qb >= 2 && dist >= 2) flag_mismatch <- TRUE
        if ((qa == 1 || qb == 1) && dist >= 4) flag_mismatch <- TRUE
      }
    }
    if (isTRUE(flag_mismatch)) risk_flags <- unique(c(risk_flags, "ontology_id_mismatch"))
  }
  if (sibling_conflict) risk_flags <- unique(c(risk_flags, "sibling_conflict"))

  if (conflict_fallback) {
    conflict_labels <- character(0)
    if (length(cand_ordered) >= 2) {
      conflict_labels <- c(cand_ordered[[1]]$label %||% NA_character_, cand_ordered[[2]]$label %||% NA_character_)
      conflict_labels <- conflict_labels[!is.na(conflict_labels) & nzchar(conflict_labels)]
    }
    if (length(conflict_labels) > 0) {
      sec <- head_out$final_decision$secondary_signals %||% list()
      if (is.character(sec)) sec <- as.list(sec)
      sec <- unique(c(as.character(sec), conflict_labels))
      head_out$final_decision$secondary_signals <- sec
    }
    head_out$post_issues$needs_manual_review <- TRUE
    head_out$post_issues$flags <- unique(c(as.character(head_out$post_issues$flags %||% character(0)), "subtype_conflict"))
  } else if (length(risk_flags) > 0) {
    head_out$post_issues$needs_manual_review <- TRUE
    head_out$post_issues$flags <- unique(c(as.character(head_out$post_issues$flags %||% character(0)), risk_flags))
  }

  if (score_missing_flag) {
    head_out$post_issues$needs_manual_review <- TRUE
    head_out$post_issues$flags <- unique(c(as.character(head_out$post_issues$flags %||% character(0)), "score_missing"))
  }

  alternatives <- character(0)
  for (cand in cand_ordered) {
    if (!is.na(cand$label) && nzchar(cand$label) && cand$label != chosen_label) {
      alternatives <- c(alternatives, cand$label)
    }
    if (length(alternatives) >= 3) break
  }
  if (length(alternatives) > 0) {
    sec <- head_out$final_decision$secondary_signals %||% list()
    if (is.character(sec)) sec <- as.list(sec)
    sec <- unique(c(as.character(sec), alternatives))
    head_out$final_decision$secondary_signals <- sec
  }

  final_method_lock_applied <- FALSE
  final_method_lock_method <- ""
  final_method_lock_from_clid <- NA_character_
  final_method_lock_to_clid <- NA_character_
  final_method_lock_recommendation <- NULL
  lock_method <- if (decision_cat == "cassia_better") "cassia" else if (decision_cat == "in_house_better") "in_house" else if (decision_cat == "enrich_better") "enrich" else ""
  if (nzchar(lock_method)) {
    lock_cand <- pick_top1_by_method(lock_method, require_gate_pass = TRUE)
    if (is.null(lock_cand)) lock_cand <- pick_top1_by_method(lock_method, require_gate_pass = FALSE)
    if (!is.null(lock_cand) && !is.na(lock_cand$clid) && nzchar(lock_cand$clid)) {
      final_method_lock_method <- lock_method
      final_method_lock_from_clid <- as.character(chosen_clid %||% NA_character_)
      final_method_lock_to_clid <- as.character(lock_cand$clid %||% NA_character_)
      final_method_lock_recommendation <- list(
        method = lock_method,
        from_clid = final_method_lock_from_clid,
        to_clid = final_method_lock_to_clid,
        to_label = as.character(lock_cand$label %||% local_lookup_by_clid(lock_cand$clid, cl_cfg) %||% NA_character_),
        note = "diagnostic_only"
      )
    }
  }

  tissue_hard_constraint_applied <- FALSE
  tissue_hard_constraint_tokens <- character(0)
  tissue_hard_constraint_roots <- as.character(tissue_filter_labels %||% character(0))
  tissue_hard_constraint_from_clid <- NA_character_
  tissue_hard_constraint_to_clid <- NA_character_
  tissue_hard_constraint_reason <- if (isTRUE(tissue_filter_applied)) {
    "deprecated_replaced_by_tissue_candidate_filter"
  } else {
    paste0("deprecated_filter_not_applied:", tissue_filter_reason)
  }



  # decision_category is metadata and does not change the label.
  # A *_better category is only valid when that reviewer's frozen top1 is ontology-compatible
  # with the ACTUAL final CLID. Otherwise repair category->tie rather than changing biology.
  decision_category_before_consistency <- as.character(decision_cat %||% "")
  decision_category_consistency_repaired <- FALSE
  decision_category_consistency_reason <- ""
  decision_category_preferred_top1_clid <- NA_character_
  decision_category_relation_to_final <- ""

  clear_method_key <- if (decision_cat == "cassia_better") {
    "cassia"
  } else if (decision_cat == "in_house_better") {
    "our"
  } else if (decision_cat == "enrich_better") {
    "enrich"
  } else {
    ""
  }

  if (nzchar(clear_method_key)) {
    pref_top1 <- pick_top1_by_method_all(clear_method_key)
    if (is.null(pref_top1) || is.na(pref_top1$clid) || !nzchar(as.character(pref_top1$clid))) {
      decision_category_consistency_repaired <- TRUE
      decision_category_consistency_reason <- "preferred_method_top1_unavailable"
      decision_category_relation_to_final <- "unavailable"
    } else {
      pref_clid <- as.character(pref_top1$clid)
      decision_category_preferred_top1_clid <- pref_clid
      if (identical(pref_clid, as.character(chosen_clid))) {
        decision_category_relation_to_final <- "exact"
      } else if (isTRUE(is_ancestor_of(pref_clid, as.character(chosen_clid), cl_graph)) ||
                 isTRUE(is_ancestor_of(as.character(chosen_clid), pref_clid, cl_graph))) {
        decision_category_relation_to_final <- "hierarchical"
      } else {
        decision_category_relation_to_final <- "incompatible"
        decision_category_consistency_repaired <- TRUE
        decision_category_consistency_reason <- paste0(
          "preferred_method_top1_incompatible_with_final:", pref_clid, "->", as.character(chosen_clid)
        )
      }
    }

    if (isTRUE(decision_category_consistency_repaired)) {
      decision_cat <- "tie"
      decision_category_v2 <- "tie_or_override"
      preferred_method <- NA_character_
      head_out$final_decision$decision_category <- "tie"
      head_out$final_decision$preferred_method <- NULL
      head_out <- flag_issue(head_out, "decision_category_repaired", decision_category_consistency_reason)
    }
  } else if (identical(decision_cat, "tie")) {
    preferred_method <- NA_character_
    head_out$final_decision$preferred_method <- NULL
  }

  # Final output normalization pass (single-point guarantee):
  # - primary_cell_type uses canonical CL label for chosen CLID
  # - secondary_signals forced to list and deduplicated against primary
  normalize_signal_final <- function(x) {
    s <- tolower(stringr::str_squish(as.character(x %||% "")))
    s <- gsub("\\([^)]*\\)", "", s)
    s <- gsub("[^a-z0-9+ -]", "", s, perl = TRUE)
    s <- gsub("\\bcells\\b", "cell", s)
    stringr::str_squish(s)
  }
  canonical_primary <- local_lookup_by_clid(chosen_clid, cl_cfg)
  if (!is.na(canonical_primary) && nzchar(canonical_primary)) {
    chosen_label <- canonical_primary
  }
  sec_vals <- head_out$final_decision$secondary_signals %||% list()
  if (!is.list(sec_vals)) sec_vals <- as.list(sec_vals)
  sec_vals <- as.character(unlist(sec_vals, recursive = TRUE, use.names = FALSE))
  sec_vals <- sec_vals[nzchar(sec_vals)]
  if (length(sec_vals) > 0) {
    pnorm <- normalize_signal_final(chosen_label)
    kept <- character(0)
    seen <- character(0)
    for (sv in sec_vals) {
      nsv <- normalize_signal_final(sv)
      if (!nzchar(nsv)) next
      if (identical(nsv, pnorm)) next
      if (nsv %in% seen) next
      kept <- c(kept, sv)
      seen <- c(seen, nsv)
    }
    head_out$final_decision$secondary_signals <- as.list(kept)
  } else {
    head_out$final_decision$secondary_signals <- list()
  }

  audit_note <- paste0(
    "lock_id=", lock_id %||% "NA",
    "; lock_reason=", lock_reason,
    "; greedy_target=", greedy_target %||% "NA",
    "; chosen=", chosen_label %||% "NA",
    "; chosen_clid=", chosen_clid %||% "NA",
    "; margin=", if (is.na(chosen_support$margin)) "NA" else sprintf("%.4f", chosen_support$margin),
    "; reasons=", paste(chosen_support$reasons, collapse = "|"),
    "; review=", paste(chosen_support$review_flags %||% character(0), collapse = "|"),
    "; alternatives=", paste(alternatives, collapse = "|")
  )
  if (length(coerce_notes) > 0) {
    audit_note <- paste(audit_note, "coerce:", paste(coerce_notes, collapse = ";"))
  }
  head_out$audit_report$notes <- str_trim(paste(as.character(head_out$audit_report$notes %||% ""), audit_note))
  if (length(rejected_by_gate) > 0) {
    rej_txt <- paste(head(rejected_by_gate, 5), collapse = "; ")
    head_out$audit_report$notes <- str_trim(paste(head_out$audit_report$notes, paste0("branch_gate_reject: ", rej_txt)))
  }

  head_out$final_decision$primary_cell_type <- chosen_label
  head_out$final_decision$final_cell_ontology_id <- chosen_clid
  head_out$final_decision$decision_category <- decision_cat
  if (!is.na(preferred_method) && nzchar(preferred_method)) {
    head_out$final_decision$preferred_method <- preferred_method
  } else {
    head_out$final_decision$preferred_method <- NULL
  }

  enabled_methods <- unique(vapply(Filter(function(x) x$rank == 1, candidates), function(x) x$method, character(1)))
  head_out$method_aliases <- list(our_method = "in_house")
  head_out$judge_config <- list(
    enabled_methods = as.list(enabled_methods),
    mode_selected = mode_selected,
    script_version = "v1.0.0",
    prompt_profile = as.character(PROMPT_PROFILE %||% "legacy"),
    run_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  head_out$judge_params <- list(
    min_depth_for_eligibility = min_depth_for_eligibility,
    delta_depth = delta_depth,
    disable_aggressive_trigger = isTRUE(disable_aggressive_trigger),
    gate_mode = gate_mode,
    gate_k = gate_k,
    enable_method_reliability_gate = isTRUE(enable_method_reliability_gate),
    enable_aggressive_lca_hard_guard = isTRUE(enable_aggressive_lca_hard_guard),
    protect_lock_inputs = isTRUE(protect_lock_inputs),
    meta_reviewer_gate = isTRUE(meta_reviewer_gate),
    identity_policy = identity_policy,
    min_votes = min_votes,
    min_margin = min_margin,
    strong_score = strong_score,
    expand_depth_k = expand_depth_k,
    score_lambda = score_lambda,
    lock_beta = lock_beta,
    conflict_delta = conflict_delta,
    specificity_eps = specificity_eps,
    branch_gate_slack = branch_gate_slack
  )
  head_out$decision_trace <- list(
    mode_selected = mode_selected,
    mode_trigger_reason = mode_trigger_reason,
    lock_id_original = lock_id_original,
    lock_id_final_used = lock_id,
    lock_id = lock_id,
    lock_source = lock_source,
    lock_source_methods = as.list(lock_source_methods),
    direct_consensus_clid = as.character(top1_ctx$direct_consensus_clid %||% NA_character_),
    direct_consensus_votes = as.integer(top1_ctx$direct_consensus_votes %||% 0L),
    direct_consensus_fraction = as.numeric(top1_ctx$direct_consensus_fraction %||% NA_real_),
    direct_consensus_reason = as.character(top1_ctx$direct_consensus_reason %||% ""),
    reviewer_anchor_id = as.character(if (!is.na(dc_anchor_clid) && nzchar(as.character(dc_anchor_clid %||% ""))) dc_anchor_clid else lock_id),
    reviewer_anchor_source = if (!is.na(dc_anchor_clid) && nzchar(as.character(dc_anchor_clid %||% ""))) "direct_top1_majority" else "ontology_mca_context",
    identity_policy = identity_policy,
    reviewer_candidate_source = if (length(judge_input_obj$inputs$frozen_reviewer_candidates %||% list()) > 0L) "step08_frozen_reviewer_candidates" else "legacy_summary_mapping",
    head_identity_valid = isTRUE(head_identity_valid),
    head_identity_clid = as.character(llm_anchor_clid %||% NA_character_),
    head_identity_label = as.character(llm_anchor_label %||% NA_character_),
    head_identity_preserved = isTRUE(head_identity_preserved_final),
    head_identity_restore_applied = isTRUE(head_identity_restore_applied),
    head_identity_candidate_present = isTRUE(head_identity_candidate_present),
    head_deterministic_relation = as.character(head_deterministic_relation %||% ""),
    identity_resolution_action = as.character(identity_resolution_action %||% ""),
    hierarchical_adjustment_applied = isTRUE(hierarchical_adjustment_applied),
    hierarchical_adjustment_candidate_present = isTRUE(hierarchical_adjustment_candidate_present),
    hierarchical_adjustment_candidate_supported = isTRUE(hierarchical_adjustment_candidate_supported),
    deterministic_pre_preserve_clid = deterministic_pre_preserve_clid,
    deterministic_pre_preserve_label = deterministic_pre_preserve_label,
    consensus_override_attempted = isTRUE(consensus_override_attempted),
    consensus_override_allowed = if (is.na(consensus_override_allowed)) NA else isTRUE(consensus_override_allowed),
    consensus_override_reason = as.character(consensus_override_reason %||% ""),
    override_candidate_clid = as.character(override_candidate_clid %||% NA_character_),
    override_relation_to_anchor = as.character(override_relation_to_anchor %||% ""),
    contrastive_shared_score = if (is.na(contrastive_shared_score)) NA_real_ else as.numeric(contrastive_shared_score),
    contrastive_candidate_unique_score = if (is.na(contrastive_candidate_unique_score)) NA_real_ else as.numeric(contrastive_candidate_unique_score),
    contrastive_anchor_unique_score = if (is.na(contrastive_anchor_unique_score)) NA_real_ else as.numeric(contrastive_anchor_unique_score),
    contrastive_margin = if (is.na(contrastive_margin)) NA_real_ else as.numeric(contrastive_margin),
    contrastive_candidate_total_score = if (is.na(contrastive_candidate_total_score)) NA_real_ else as.numeric(contrastive_candidate_total_score),
    contrastive_anchor_total_score = if (is.na(contrastive_anchor_total_score)) NA_real_ else as.numeric(contrastive_anchor_total_score),
    contrastive_candidate_real_methods = as.list(as.character(contrastive_candidate_real_methods %||% character(0))),
    contrastive_override_min_margin = as.numeric(contrastive_override_min_margin),
    consensus_override_forced_anchor = isTRUE(consensus_override_forced_anchor),
    consensus_override_gate_stage = consensus_override_gate_stage,
    pre_consensus_final_clid = pre_consensus_final_clid,
    pre_consensus_final_label = pre_consensus_final_label,
    direct_consensus_anchor_reintroduced = isTRUE(direct_consensus_anchor_reintroduced),
    decision_category_before_consistency = decision_category_before_consistency,
    decision_category_after_consistency = decision_cat,
    decision_category_consistency_repaired = isTRUE(decision_category_consistency_repaired),
    decision_category_consistency_reason = decision_category_consistency_reason,
    decision_category_preferred_top1_clid = decision_category_preferred_top1_clid,
    decision_category_relation_to_final = decision_category_relation_to_final,
    lock_reason = lock_reason,
    lock_depth = lock_depth,
    supported_lca_clid = supported_lca_clid,
    supported_lca_depth = supported_lca_depth,
    depth_gap = depth_gap,
    eligible_n = eligible_n,
    support_n = support_n,
    support_gate_mode = gate_mode,
    support_gate_threshold = gate_k,
    support_gate_threshold_effective = as.integer(gate_k),
    support_gate_pass = isTRUE(support_gate_pass),
    branch_gate_relaxed_due_to_missing_support = if (isTRUE(branch_gate_relaxed_due_to_missing_support)) TRUE else NULL,
    decision_category_original = decision_cat,
    decision_category_v2 = decision_category_v2,
    preferred_method = preferred_method,
    final_rule = final_rule,
    final_method_lock_applied = isTRUE(final_method_lock_applied),
    final_method_lock_method = final_method_lock_method,
    final_method_lock_from_clid = final_method_lock_from_clid,
    final_method_lock_to_clid = final_method_lock_to_clid,
    final_method_lock_recommendation = final_method_lock_recommendation,
    tissue_filter_applied = isTRUE(tissue_filter_applied),
    tissue_filter_reason = tissue_filter_reason,
    tissue_filter_labels = as.list(as.character(tissue_filter_labels %||% character(0))),
    tissue_filter_n_before = as.integer(tissue_filter_n_before),
    tissue_filter_n_after = as.integer(tissue_filter_n_after),
    tissue_filter_results = lapply(tissue_filter_results, function(x) {
      list(
        id = as.character(x$id %||% ""),
        label = as.character(x$label %||% ""),
        verdict = as.character(x$verdict %||% "uncertain"),
        confidence = as.numeric(x$confidence %||% 0.5),
        reason = as.character(x$reason %||% "")
      )
    }),
    tissue_filter_summary = tissue_filter_summary,
    tissue_hard_constraint_applied = isTRUE(tissue_hard_constraint_applied),
    tissue_hard_constraint_tokens = as.list(as.character(tissue_hard_constraint_tokens %||% character(0))),
    tissue_hard_constraint_roots = as.list(as.character(tissue_hard_constraint_roots %||% character(0))),
    tissue_hard_constraint_from_clid = tissue_hard_constraint_from_clid,
    tissue_hard_constraint_to_clid = tissue_hard_constraint_to_clid,
    tissue_hard_constraint_reason = tissue_hard_constraint_reason,

    aggressive_lca_hard_guard_applied = isTRUE(aggressive_lca_hard_guard_applied),
    method_reliability = lapply(names(method_reliability), function(m) {
      r <- method_reliability[[m]]
      list(method = m, score = as.numeric(r$score %||% NA_real_), is_primary = isTRUE(r$is_primary), reasons = as.list(as.character(r$reasons %||% character(0))))
    }),
    meta_reviewer_decisions = meta_reviewer_decisions,
    anchor_fallback_happened = anchor_fallback_happened,
    post_fallback_override_applied = post_fallback_override_applied,
    post_fallback_override_reason = post_fallback_override_reason,
    post_fallback_why_not_override = post_fallback_why_not_override,
    post_fallback_rule1_gate_n = as.integer(post_fallback_rule1_gate_n),
    post_fallback_rule1_soft_n = as.integer(post_fallback_rule1_soft_n),
    post_fallback_rule2_gate_n = as.integer(post_fallback_rule2_gate_n),
    post_fallback_rule2_soft_n = as.integer(post_fallback_rule2_soft_n),
    lock_unsteady = isTRUE(lock_unsteady),
    lock_unsteady_reasons = lock_unsteady_reasons,
    instability_pool_n = as.integer(instability_pool_n),
    gate_pool_n = as.integer(gate_pool_n),
    soft_pool_n = as.integer(soft_pool_n),
    shadow_pool_n = as.integer(shadow_pool_n),
    gate_reject_total_n = as.integer(gate_reject_total_n),
    gate_reject_descendant_n = as.integer(gate_reject_descendant_n),
    gate_reject_descendant_by_stage = gate_reject_descendant_by_stage,
    shadow_kept_n = as.integer(shadow_kept_n),
    shadow_stage_flagged_n = as.integer(shadow_stage_flagged_n),
    shadow_empty_reason = shadow_empty_reason,
    shadow_desc_filter_counts = as.character(jsonlite::toJSON(shadow_desc_filter_counts %||% list(), auto_unbox = TRUE)),
    shadow_desc_ref = shadow_desc_ref,
    shadow_pool_used = isTRUE(shadow_pool_used),
    support_up = isTRUE(override_checks$support_up %||% FALSE),
    score_ok = isTRUE(override_checks$score_ok %||% FALSE),
    score_na_case = score_na_case,
    consistency_ok = isTRUE(consistency_ok),
    consistency_ok_reason = consistency_ok_reason,
    score_ok_reason = score_ok_reason,
    score_source_base = score_source_base,
    score_source_new = score_source_new,
    effective_base_score = effective_base_score,
    effective_new_score = effective_new_score,
    base_depth = as.numeric(base_depth),
    candidate_depth = as.numeric(candidate_depth),
    base_subtree_size = as.integer(base_subtree_size),
    base_is_ultra_generic = isTRUE(base_is_ultra_generic),
    pool_size = as.integer(override_checks$pool_size %||% 0L),
    pool_size_1_strong_margin_used = isTRUE(pool_size_1_strong_margin_used),
    override_checks = override_checks,
    guardrail = guardrail_trace,
    cassia_too_generic_override = cassia_override_trace,
    label_tiebreak = label_tiebreak_trace,
    reviewer_diagnostics = reviewer_diagnostics
    ,n_raw_candidates_by_method = n_raw_candidates_by_method
    ,n_with_valid_clid = as.integer(n_with_valid_clid)
    ,n_after_meta_reviewer_gate = as.integer(n_after_meta_reviewer_gate)
    ,n_after_anchor_constraints = as.integer(n_after_anchor_constraints)
    ,n_in_cand_table = as.integer(n_in_cand_table)
    ,n_in_cand_ordered = as.integer(n_in_cand_ordered)
    ,eligibility_drop_reasons_top = eligibility_drop_reasons_top
  )
  cand_dump <- lapply(cand_ordered, function(cand) {
    methods_real <- character(0)
    if (!is.na(cand$clid) && nzchar(cand$clid)) {
      methods_real <- clid_to_methods[[cand$clid]] %||% character(0)
    }
    list(
      label = cand$label,
      clid = cand$clid,
      score = cand$score,
      score_adj = cand$score_adj,
      method = cand$method,
      rank = cand$rank,
      method_count = cand$method_count,
      method_count_all = cand$method_count_all,
      methods_supporting = methods_real,
      map_quality = cand$map_quality,
      map_quality_rank = cand$map_quality_rank,
      specificity = cand$specificity,
      depth_from_lock = cand$depth_from_lock,
      dist_to_lock = cand$dist_to_lock,
      depth_to_root = cand$depth_to_root,
      lock_penalty = cand$lock_penalty,
      evidence_support = cand$evidence_support,
      evidence_top_hits = cand$evidence_top_hits,
      evidence_recovery_hits = cand$evidence_recovery_hits,
      evidence_unknown_hits = cand$evidence_unknown_hits,
      evidence_out_hits = cand$evidence_out_hits,
      anchor_bonus = cand$anchor_bonus,
      guardrail_penalty = cand$guardrail_penalty,
      contamination_penalty = cand$contamination_penalty
      ,reliability_score = cand$reliability_score
      ,reliability_penalty = cand$reliability_penalty
      ,reliability_reasons = as.list(as.character(cand$reliability_reasons %||% character(0)))
    )
  })
  all_real_top1_dump <- lapply(Filter(function(cand) {
    as.integer(cand$rank %||% 99L) == 1L &&
      as.character(cand$method %||% "") %in% real_input_methods
  }, cand_table_all), function(cand) {
    list(
      method = as.character(cand$method %||% ""),
      label = as.character(cand$label %||% ""),
      clid = as.character(cand$clid %||% NA_character_),
      raw_clid = as.character(cand$raw_clid %||% NA_character_),
      coerced_status = as.character(cand$coerced_status %||% ""),
      map_quality = as.numeric(cand$map_quality %||% NA_real_),
      evidence_support = as.numeric(cand$evidence_support %||% 0),
      score = as.numeric(cand$score %||% NA_real_),
      score_adj = as.numeric(cand$score_adj %||% NA_real_)
    )
  })

  attr(head_out, "diag") <- list(
    lock_clid = lock_id,
    num_candidates_before = num_candidates_before,
    num_candidates_after = num_candidates_after,
    chosen_pred_main_clid = chosen_clid,
    chosen_pred_subtypes = alternatives,
    expanded_added = num_candidates_after - num_candidates_before,
    specificity_tiebreak = specificity_tiebreak,
    score_missing_flag = score_missing_flag,
      candidate_dump = list(
        lock_id = lock_id,
        lock_reason = lock_reason,
        direct_consensus_clid = as.character(top1_ctx$direct_consensus_clid %||% NA_character_),
        direct_consensus_votes = as.integer(top1_ctx$direct_consensus_votes %||% 0L),
        direct_consensus_fraction = as.numeric(top1_ctx$direct_consensus_fraction %||% NA_real_),
        direct_consensus_reason = as.character(top1_ctx$direct_consensus_reason %||% ""),
        min_votes = min_votes,
        min_margin = min_margin,
        strong_score = strong_score,
      score_lambda = score_lambda,
      consensus_mode = consensus_mode,
      chosen_clid = chosen_clid,
      chosen_label = chosen_label,
      chosen_score = chosen$score,
      chosen_score_adj = chosen$score_adj,
      direct_consensus_anchor_reintroduced = isTRUE(direct_consensus_anchor_reintroduced),
      all_real_top1 = all_real_top1_dump,
      candidates = cand_dump
    )
  )
  head_out
}


run_single_cluster_test <- function(input_json_path,
                                    out_root = NULL,
                                    expected_substring = NULL,
                                    cl_cfg = NULL,
                                    cl_graph = NULL,
                                    expand_depth_k = 3L) {
  jin <- read_json_safely(input_json_path)
  if (is.null(jin) || !is.list(jin)) stop("Invalid input JSON: ", input_json_path)
  cid <- as.character(jin$cluster_id %||% NA_character_)

  if (is.null(cl_cfg) || is.null(cl_graph)) {
    cl_cfg <- make_cl_cfg(Sys.getenv("CL_LOCAL_JSON", unset = file.path(Sys.getenv("TRIAGE_HOME", unset = getwd()), "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json")), prefer_ols = FALSE, cache_dir = ".ols_cache")
    cl_graph <- load_cl_graph(cl_cfg)
  }

  head_out <- NULL
  if (!is.null(out_root) && nzchar(out_root)) {
    out_head <- file.path(out_root, "head_editor_outputs")
    head_out <- load_latest_head_output(cid, out_head)
  }
  if (is.null(head_out)) stop("Missing head output for cluster: ", cid)

  cite_req <- extract_citation_requirements(jin)
  head_out <- postprocess_head_out(
    head_out, jin, cite_req, cl_cfg, cl_graph,
    dataset_cfg = if (exists("cfg", inherits = TRUE)) cfg else NULL,
    species_value = if (exists("cfg", inherits = TRUE) && !is.null(cfg$species)) cfg$species[[1]] %||% "human" else "human",
    release_policy = "auto"
  )
  chosen_label <- as.character(head_out$final_decision$primary_cell_type %||% NA_character_)
  chosen_clid <- as.character(head_out$final_decision$final_cell_ontology_id %||% NA_character_)
  needs_review <- isTRUE(head_out$post_issues$needs_manual_review %||% FALSE)

  final_sim <- head_out
  final_sim <- normalize_manual_review_plan(final_sim)
  final_sim <- force_array_fields(final_sim)
  final_sim <- ensure_enrich_verdict(final_sim, jin)
  final_sim <- normalize_judge_final_decision_cl_preserve(final_sim, cl_cfg)
  final_label <- as.character(final_sim$final_decision$primary_cell_type %||% NA_character_)
  final_clid <- as.character(final_sim$final_decision$final_cell_ontology_id %||% NA_character_)
  if (!identical(chosen_label, final_label) || !identical(chosen_clid, final_clid)) {
    cat("[WARN] final post-processing altered choice\n")
    cat("  chosen=", chosen_label, " / ", chosen_clid, "\n", sep="")
    cat("  final =", final_label, " / ", final_clid, "\n", sep="")
  }

  inputs <- jin$inputs %||% list()
  cand_top1 <- c(
    collect_method_candidates("cassia", inputs$cassia_summary %||% list(), cl_cfg, cl_graph, max_k = 1),
    collect_method_candidates("our", inputs$our_summary %||% list(), cl_cfg, cl_graph, max_k = 1),
    collect_method_candidates("enrich", inputs$enrich_summary %||% inputs$inter_summary %||% list(), cl_cfg, cl_graph, max_k = 1)
  )

  cat("method_clids:\n")
  for (cand in cand_top1) {
    cat(" - ", cand$method, ": raw=", cand$raw_clid %||% "NA",
        ", coerced=", cand$clid %||% "NA",
        ", status=", cand$coerced_status %||% "NA",
        ", map_q=", cand$map_quality %||% NA, "\n", sep = "")
  }

  candidates <- c(
    collect_method_candidates("cassia", inputs$cassia_summary %||% list(), cl_cfg, cl_graph),
    collect_method_candidates("our", inputs$our_summary %||% list(), cl_cfg, cl_graph),
    collect_method_candidates("enrich", inputs$enrich_summary %||% inputs$inter_summary %||% list(), cl_cfg, cl_graph)
  )

  top1_candidates <- Filter(function(x) x$rank == 1, candidates)
  top1_clids <- unique(vapply(top1_candidates, function(x) x$clid %||% NA_character_, character(1)))
  top1_clids <- top1_clids[!is.na(top1_clids) & nzchar(top1_clids)]

  lock_id <- NA_character_
  if (length(top1_clids) == 1) {
    lock_id <- top1_clids[[1]]
  } else if (length(top1_clids) > 1) {
    anc_maps <- lapply(top1_clids, function(clid) get_anc_map(clid, cl_graph))
    votes <- list()
    dist_sum <- list()
    dist_cnt <- list()
    for (am in anc_maps) {
      if (is.null(am) || length(am) == 0) next
      for (id in names(am)) {
        if (is.null(votes[[id]])) votes[[id]] <- 0
        votes[[id]] <- votes[[id]] + 1
        if (is.null(dist_sum[[id]])) dist_sum[[id]] <- 0
        dist_sum[[id]] <- dist_sum[[id]] + as.numeric(am[[id]])
        if (is.null(dist_cnt[[id]])) dist_cnt[[id]] <- 0
        dist_cnt[[id]] <- dist_cnt[[id]] + 1
      }
    }
    ids <- names(votes)
    vcnt <- vapply(ids, function(id) votes[[id]], numeric(1))
    keep <- ids[vcnt >= 2]
    if (length(keep) > 0) {
      mean_dist <- vapply(keep, function(id) {
        ds <- dist_sum[[id]] %||% 1
        dc <- dist_cnt[[id]] %||% 1
        ds / dc
      }, numeric(1))
      keep_df <- data.frame(id = keep, votes = vcnt[keep], mean_dist = mean_dist, stringsAsFactors = FALSE)
      keep_df <- keep_df[order(-keep_df$votes, keep_df$mean_dist), , drop = FALSE]
      lock_id <- keep_df$id[[1]]
    }
  }

  vcat("lock_id=", lock_id %||% "NA", "\n", sep = "")
  vcat("expand_depth_k=", as.integer(expand_depth_k), "\n", sep = "")

  in_lock <- function(clid) {
    if (is.na(lock_id) || !nzchar(lock_id)) return(TRUE)
    if (is.na(clid) || !nzchar(clid)) return(FALSE)
    isTRUE(is_ancestor_of(lock_id, clid, cl_graph)) || isTRUE(clid == lock_id)
  }

  pool <- candidates
  scores <- vapply(pool, function(x) x$score %||% NA_real_, numeric(1))
  scores <- scores[!is.na(scores)]
  if (length(scores) > 0) {
    top1 <- max(scores)
    top2 <- if (length(scores) > 1) sort(scores, decreasing = TRUE)[[2]] else NA_real_
    margin <- if (is.finite(top2)) (top1 - top2) else NA_real_
    vcat("top2_margin=", if (is.na(margin)) "NA" else sprintf("%.4f", margin), "\n", sep = "")
  } else {
    vcat("top2_margin=NA\n")
  }

  greedy_target <- NA_character_
  method_count <- list()
  for (cand in pool) {
    clid <- cand$clid %||% NA_character_
    if (is.na(clid) || !nzchar(clid)) next
    if (is.null(method_count[[clid]])) method_count[[clid]] <- list()
    method_count[[clid]] <- unique(c(method_count[[clid]], cand$method))
  }
  cand_table <- lapply(pool, function(cand) {
    clid <- cand$clid %||% NA_character_
    mcnt <- if (!is.na(clid) && !is.null(method_count[[clid]])) length(method_count[[clid]]) else 1L
    spec <- if (!is.na(clid) && nzchar(clid)) get_depth_to_root(clid, cl_graph) else 0
    list(clid = clid, score = cand$score, method_count = mcnt, specificity = spec)
  })
  c2 <- Filter(function(x) !is.na(x$clid) && x$method_count >= 1, cand_table)
  if (length(c2) > 0) {
    ord <- order(vapply(c2, function(x) x$score %||% -Inf, numeric(1)),
                 vapply(c2, function(x) x$specificity, numeric(1)),
                 decreasing = TRUE)
    greedy_target <- c2[[ord[[1]]]]$clid
  }
  vcat("greedy_target=", greedy_target %||% "NA", "\n", sep = "")

  vcat("cluster_id=", cid, "\n")
  vcat("chosen_label=", chosen_label, "\n")
  vcat("chosen_clid=", chosen_clid, "\n")
  vcat("needs_manual_review=", needs_review, "\n")
  diag <- attr(head_out, "diag")
  if (!is.null(diag) && is.list(diag)) {
    vcat("expanded_added=", as.character(diag$expanded_added %||% "NA"), "\n", sep = "")
    vcat("specificity_tiebreak=", as.character(diag$specificity_tiebreak %||% "NA"), "\n", sep = "")
    vcat("score_missing_flag=", as.character(diag$score_missing_flag %||% "NA"), "\n", sep = "")
  }
  flags <- head_out$post_issues$flags %||% list()
  if (is.list(flags)) flags <- unlist(flags, recursive = TRUE, use.names = FALSE)
  flags <- as.character(flags)
  flags <- flags[nzchar(flags)]
  vcat("review_flags=", if (length(flags) > 0) paste(flags, collapse = ";") else "", "\n", sep = "")
  audit_note <- as.character(head_out$audit_report$notes %||% "")
  vcat("audit_notes=", audit_note, "\n", sep = "")

  if (!is.null(expected_substring) && nzchar(expected_substring)) {
    if (is.na(chosen_label) || !grepl(expected_substring, chosen_label, fixed = TRUE)) {
      vlog("[WARN] expected_substring not found in chosen_label")
    }
  }

  invisible(head_out)
}

# --------------------------
# Local hard gate (deterministic)
# UPDATED: removed any marker_panels-based hard constraints
# --------------------------
has_valid_evidence_pointers <- function(ep) {
  if (is.null(ep)) return(FALSE)
  if (!is.list(ep) || length(ep) == 0) return(FALSE)
  if (!is.list(ep[[1]])) return(FALSE)
  ok <- vapply(ep, function(x) {
    t <- as.character(x$type %||% "")
    v <- as.character(x$value %||% "")
    nzchar(t) && nzchar(v)
  }, logical(1))
  any(ok)
}

judge_local_gate <- function(j, citation_requirements = NULL, judge_input_obj = NULL) {
  reasons <- character(0)
  j <- ensure_audit_report_support(j)
  
  req_top <- rules_top_level_keys()
  missing_top <- req_top[!req_top %in% names(j)]
  if (length(missing_top) > 0) reasons <- c(reasons, paste0("missing_top_keys: ", paste(missing_top, collapse=",")))
  
  if (!is.list(j$final_decision)) reasons <- c(reasons, "final_decision_not_object")
  else {
    fct <- as.character(j$final_decision$primary_cell_type %||% "")
    if (!nzchar(fct)) reasons <- c(reasons, "missing_primary_cell_type")
    if (nzchar(fct) && grepl(banned_tokens_pattern(), tolower(fct))) reasons <- c(reasons, "primary_cell_type_contains_banned_token")
    
    conf <- suppressWarnings(as.numeric(j$final_decision$confidence_primary %||% NA_real_))
    if (is.na(conf) || conf < 0 || conf > 1) reasons <- c(reasons, "confidence_primary_out_of_range")
    
    dc <- as.character(j$final_decision$decision_category %||% "")
    if (!(dc %in% rules_allowed_decision_categories())) reasons <- c(reasons, "decision_category_invalid")
  }
  
  if (!is.list(j$method_verdict)) reasons <- c(reasons, "method_verdict_not_object")

  if (!is.list(j$audit_report)) {
    reasons <- c(reasons, "audit_report_not_object")
  } else {
    rs <- j$audit_report$reviewer_support %||% NULL
    cs <- rs$cassia_supported %||% NULL
    os <- rs$our_supported %||% NULL
    es <- rs$enrich_supported %||% NULL
    if (!is.logical(cs) || length(cs) != 1 || is.na(cs) ||
        !is.logical(os) || length(os) != 1 || is.na(os) ||
        !is.logical(es) || length(es) != 1 || is.na(es)) {
      reasons <- c(reasons, "audit_report_invalid_reviewer_support")
    }
  }
  
  # Citation requirement: must cite at least one allowed PMID and/or enrichment term when provided.
  if (!is.list(j$evidence)) {
    reasons <- c(reasons, "evidence_not_object")
  } else {
    citation_check <- validate_citations(j, citation_requirements)
    if (length(citation_check$reasons) > 0) reasons <- c(reasons, citation_check$reasons)
  }

  rs <- j$audit_report$reviewer_support %||% NULL
  cs <- rs$cassia_supported %||% NULL
  os <- rs$our_supported %||% NULL
  es <- rs$enrich_supported %||% NULL
  all_reviewers_unsupported <- isFALSE(cs) && isFALSE(os) && isFALSE(es)

  if (all_reviewers_unsupported) {
    tpa <- j$third_party_adjudication %||% NULL
    if (!is.list(tpa)) {
      reasons <- c(reasons, "third_party_adjudication_missing")
    } else {
      tpa_label <- as.character(tpa$primary_cell_type %||% "")
      if (!nzchar(tpa_label)) reasons <- c(reasons, "third_party_primary_cell_type_missing")
      if (!has_valid_evidence_pointers(tpa$evidence_pointers %||% NULL)) reasons <- c(reasons, "third_party_evidence_pointers_missing")
    }
    dc <- as.character(j$final_decision$decision_category %||% "")
    if (dc != "third_party_override") reasons <- c(reasons, "decision_category_not_third_party_override")
  }

  needs_review <- isTRUE(j$post_issues$needs_manual_review %||% FALSE)
  if (needs_review) {
    mrp <- j$manual_review_plan %||% NULL
    if (!is.list(mrp)) {
      reasons <- c(reasons, "manual_review_plan_missing")
    } else {
      actions <- mrp$actions %||% NULL
      if (!is.list(actions) || length(actions) == 0) {
        reasons <- c(reasons, "manual_review_actions_missing")
      } else if (isTRUE(RULES$MANUAL_REVIEW$require_evidence_pointers)) {
        for (a in actions) {
          if (!has_valid_evidence_pointers(a$evidence_pointers %||% NULL)) reasons <- c(reasons, "manual_review_action_evidence_pointers_missing")
        }
      }
    }
  }
  
  list(ok = length(reasons)==0, reasons = reasons)
}

validate_rules_consistency <- function() {
  reasons <- character(0)
  head_prompt <- build_head_editor_system_prompt(
    species_value = if (exists("cfg", inherits = TRUE) && !is.null(cfg$species)) cfg$species[[1]] %||% "human" else "human")
  chief_task <- build_chief_editor_instructions()$task %||% ""
  chief_rules_text <- rules_text_for_chief() %||% ""
  chief_combined <- paste(chief_task, chief_rules_text, sep = "\n")
  schema_txt <- judge_output_schema_text()
  dec_text <- rules_decision_categories_text()
  ban_text <- rules_banned_tokens_text()
  if (!stringr::str_detect(head_prompt, stringr::fixed(dec_text))) reasons <- c(reasons, "head_prompt_missing_decision_categories")
  if (!stringr::str_detect(chief_combined, stringr::fixed(dec_text))) reasons <- c(reasons, "chief_prompt_missing_decision_categories")
  if (!stringr::str_detect(head_prompt, stringr::fixed(ban_text))) reasons <- c(reasons, "head_prompt_missing_banned_tokens")
  if (!stringr::str_detect(chief_combined, stringr::fixed(ban_text))) reasons <- c(reasons, "chief_prompt_missing_banned_tokens")
  if (!stringr::str_detect(schema_txt, stringr::fixed(dec_text))) reasons <- c(reasons, "schema_missing_decision_categories")
  list(ok = length(reasons) == 0, reasons = reasons)
}

should_run_chief <- function(head_out, cite_req, always_run=FALSE, gate_ok=TRUE,
                             release_policy="legacy") {
  if (isTRUE(always_run)) return(list(run=TRUE, reason="always_run"))
  if (!isTRUE(gate_ok)) return(list(run=TRUE, reason="local_gate_failed"))
  if (!chief_precheck_ok(head_out, cite_req)) return(list(run=TRUE, reason="precheck_failed"))
  conf <- suppressWarnings(as.numeric(head_out$final_decision$confidence_primary %||% NA_real_))
  needs <- isTRUE(head_out$post_issues$needs_manual_review %||% FALSE)
  if (!is.na(conf) && conf < 0.5) return(list(run=TRUE, reason="low_confidence"))

  if (release_policy_is_v3(release_policy)) {
    blockers <- get_release_blockers(head_out)
    auto_qc <- get_auto_qc_flags(head_out)
    if (needs || length(blockers) > 0L) {
      return(list(run=TRUE, reason="release_blockers_present"))
    }
    if (length(auto_qc) > 0L) {
      return(list(run=TRUE, reason="automated_qc_required"))
    }
    return(list(run=FALSE, reason="fast_path_audit_only_or_clean"))
  }

  if (release_policy_is_v2(release_policy)) {
    blockers <- get_release_blockers(head_out)
    if (needs || length(blockers) > 0L) return(list(run=TRUE, reason="release_blockers_present"))
    return(list(run=FALSE, reason="fast_path_audit_only_or_clean"))
  }

  flags_raw <- normalize_issue_flags(head_out$post_issues$flags %||% list())
  if (needs) return(list(run=TRUE, reason="needs_review"))
  if (length(flags_raw) > 0L) return(list(run=TRUE, reason="flags_present"))
  list(run=FALSE, reason="fast_path")
}

chief_precheck_ok <- function(head_out, cite_req) {
  if (is.null(head_out) || !is.list(head_out)) return(FALSE)
  citation_check <- validate_citations(head_out, cite_req)
  if (!isTRUE(citation_check$ok)) return(FALSE)
  if (isTRUE(head_out$post_issues$needs_manual_review %||% FALSE) && isTRUE(RULES$MANUAL_REVIEW$require_actions)) {
    actions <- head_out$manual_review_plan$actions %||% list()
    if (!is.list(actions) || length(actions) == 0) return(FALSE)
  }
  TRUE
}

attach_chief_qc_failure <- function(head_out, chief_out) {
  if (is.null(head_out) || !is.list(head_out)) return(head_out)
  if (is.null(head_out$post_issues) || !is.list(head_out$post_issues)) {
    head_out$post_issues <- list(needs_manual_review = TRUE, flags = character(0), notes = "")
  }
  head_out$post_issues$needs_manual_review <- TRUE
  flags <- head_out$post_issues$flags %||% character(0)
  flags <- unique(c(as.character(flags), "chief_qc_failed"))
  head_out$post_issues$flags <- flags
  
  fr <- NULL
  if (!is.null(chief_out) && is.list(chief_out)) fr <- chief_out$failure_reasons %||% NULL
  fr_txt <- ""
  if (!is.null(fr)) {
    if (is.character(fr)) fr_txt <- paste(fr, collapse = " | ")
    else fr_txt <- tryCatch(jsonlite::toJSON(fr, auto_unbox = TRUE, null = "null"), error = function(e) "")
  }
  
  note0 <- str_trim(as.character(head_out$post_issues$notes %||% ""))
  note1 <- if (nzchar(fr_txt)) paste0("Chief QC failed: ", fr_txt) else "Chief QC failed."
  head_out$post_issues$notes <- str_trim(paste(note0, note1))
  head_out
}

# --------------------------
# Run one cluster
# --------------------------

# Banned-token pattern helper (09_run_judge.R lines 507-508)
banned_tokens_pattern <- function() rules_banned_tokens_pattern()
