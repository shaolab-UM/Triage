# =============================================================
# 00_cl_linker.R — CL-Linker: shared ontology normalization layer
#
# Design goals:
#   Candidate retrieval: high recall — maximize the chance that the correct CL term enters the candidate set
#   Final mapping:      high precision — release a unique CL ID only when evidence is sufficient
#
# Intended use (shared across workflows and not restricted to Triage):
#   - Triage (09 judge)
#   - ontology-aware deterministic selector (15)
#   - benchmark scoring (11)
#   - future user-supplied annotation results
#
# Output states:
#   mapped | ambiguous | unmapped   (ambiguous/unmapped are safety states rather than execution failures)
#
# Layers:
#   Layer 0: reviewer-provided CL ID (validate and prioritize when present)
#   Layer 1: deterministic dictionary (canonical label + synonyms + normalization)
#   Layer 2: versioned alias dictionary (cl_aliases.tsv, optional)
# =============================================================


# ---- Versioned alias table (v1, aligned with cl_normalizer canonicalization) ----
CL_ALIASES <- c(
  # neural
  "excitatory neuron" = "glutamatergic neuron",
  "cortical excitatory neuron" = "glutamatergic neuron",
  "inhibitory neuron" = "GABAergic neuron",
  "gaba neuron" = "GABAergic neuron",
  "neuronal" = "neuron",
  "neuronal progenitor" = "neural progenitor cell",
  "neuronal progenitor cell" = "neural progenitor cell",
  "intermediate progenitor" = "neural progenitor cell",
  "intermediate progenitor cell" = "neural progenitor cell",
  "radial glia" = "radial glial cell",
  "radial glia-1" = "radial glial cell",
  "radial glia-2" = "radial glial cell",
  "cycling progenitor" = "neural progenitor cell",
  "choroid" = "choroid plexus epithelial cell",
  "mesenchyme" = "mesenchymal cell",
  "retinal pigment epithelium" = "retinal pigment epithelial cell",
  "retinal progenitor" = "retinal progenitor cell",
  # immune
  "t lymphocyte" = "T cell",
  "t lymphocytes" = "T cell",
  "b lymphocyte" = "B cell",
  "b lymphocytes" = "B cell",
  "natural killer cells" = "natural killer cell",
  "nk cell" = "natural killer cell",
  "nk cells" = "natural killer cell",
  "plasma cells" = "plasma cell",
  "endothelial cells" = "endothelial cell",
  "mesothelial cells" = "mesothelial cell",
  "macrophage (myeloid lineage" = "macrophage",
  "monocyte (myeloid lineage)" = "monocyte",
  "cytotoxic ab t lymphocyte" = "cytotoxic T cell",
  # general
  "stem cell" = "stem cell",
  "progenitor cell" = "progenitor cell"
)

alias_lookup <- function(label, idx) {
  if (is.na(label) || !nzchar(label)) return(NA_character_)
  key <- idx$norm_form(label)
  for (i in seq_along(CL_ALIASES)) {
    alias_key <- idx$norm_form(names(CL_ALIASES)[i])
    if (identical(key, alias_key)) {
      target <- unname(CL_ALIASES[i])
      # normalize the target before index lookup
      nk <- idx$norm_form(target)
      if (exists(nk, envir = idx$norm_idx, inherits = FALSE)) {
        return(get(nk, envir = idx$norm_idx))
      }
    }
  }
  NA_character_
}

# ---- Build Cell Ontology indexes ----
build_cl_index <- function(cl_json_path) {
  cl <- jsonlite::fromJSON(cl_json_path, simplifyVector = FALSE)

  # label -> CL ID index (canonical label)
  label_idx <- new.env(parent = emptyenv())
  # normalized form -> CL ID (conservative normalization: case, pluralization, punctuation, Greek letters and CD markers)
  norm_idx <- new.env(parent = emptyenv())
  # synonym -> CL ID
  syn_idx <- new.env(parent = emptyenv())
  # CL ID -> info
  clid_info <- new.env(parent = emptyenv())

  # ---- conservative normalization (without changing biological meaning) ----
  norm_form <- function(x) {
    x <- tolower(x)
    # Greek letters
    x <- gsub("\u03b1", "alpha", x)   # α
    x <- gsub("\u03b2", "beta", x)    # β
    x <- gsub("\u03b3", "gamma", x)   # γ
    x <- gsub("\u03b4", "delta", x)   # δ
    # CD markers: CD4+ / CD4 positive / CD4pos
    x <- gsub("cd([0-9]+)[+]?", "cd\\1positive", x)
    x <- gsub("cd([0-9]+)pos", "cd\\1positive", x)
    # convert hyphens, underscores and slashes to spaces
    x <- gsub("[-_/]", " ", x)
    # normalize simple trailing plural s, e.g. cells -> cell
    x <- gsub("cells\\b", "cell", x)
    x <- gsub("lymphocytes\\b", "lymphocyte", x)
    x <- gsub("monocytes\\b", "monocyte", x)
    x <- gsub("neurons\\b", "neuron", x)
    x <- gsub("cells$", "cell", x)
    # normalize repeated whitespace and non-alphanumeric characters
    x <- gsub("[^a-z0-9]", "", x)
    x
  }

  for (id in names(cl)) {
    term <- cl[[id]]
    lbl <- term$label %||% ""
    syns <- term$synonyms %||% list()
    if (is.list(syns)) syns <- unlist(syns, recursive = TRUE, use.names = FALSE)
    syns <- as.character(syns)

    info <- list(
      cl_id = id,
      canonical_label = lbl,
      synonyms = syns,
      obsolete = isTRUE(term$obsolete %||% FALSE),
      replacement = term$replacement %||% NA_character_
    )
    assign(id, info, envir = clid_info)

    if (nzchar(lbl)) {
      assign(lbl, id, envir = label_idx)                   # canonical exact
      nk <- norm_form(lbl)
      if (!exists(nk, envir = norm_idx, inherits = FALSE)) assign(nk, id, envir = norm_idx)
    }
    for (s in syns) {
      if (!nzchar(s)) next
      if (!exists(s, envir = syn_idx, inherits = FALSE)) assign(s, id, envir = syn_idx)
      nk <- norm_form(s)
      if (!exists(nk, envir = norm_idx, inherits = FALSE)) assign(nk, id, envir = norm_idx)
    }
  }

  list(
    cl = cl,
    label_idx = label_idx,
    norm_idx = norm_idx,
    syn_idx = syn_idx,
    clid_info = clid_info,
    norm_form = norm_form
  )
}

# ---- Layer 0: reviewer-provided CL ID ----
# validate the provided CL ID (must exist in the ontology and not be obsolete)
validate_provided_clid <- function(clid, idx, raw_label = NA_character_) {
  if (is.na(clid) || !nzchar(clid)) return(list(ok = FALSE, reason = "missing"))
  if (!startsWith(clid, "CL:")) return(list(ok = FALSE, reason = "invalid_format"))
  info <- idx$clid_info[[clid]]
  if (is.null(info)) return(list(ok = FALSE, reason = "not_in_ontology"))
  if (isTRUE(info$obsolete)) return(list(ok = FALSE, reason = "obsolete"))
  # : label-consistency check (the canonical label for the provided CL ID must be compatible with the raw label)
  # token-level check: label canon , (stricter than character-substring matching)
  if (!is.na(raw_label) && nzchar(raw_label)) {
    canon <- info$canonical_label %||% ""
    words_of <- function(x) sort(unique(strsplit(tolower(x), "[^a-z0-9]+")[[1]]))
    wl <- words_of(raw_label)
    wc <- words_of(canon)
    compatible <- length(wl) > 0 && length(wc) > 0 &&
                  (all(wl %in% wc) || all(wc %in% wl))
    if (!compatible) {
      return(list(ok = FALSE, reason = paste0("label_mismatch: '", raw_label, "' vs canonical '", canon, "'")))
    }
  }
  list(ok = TRUE, reason = "provided_valid", canonical_label = info$canonical_label)
}

# ---- Layer 1: deterministic dictionary ----
# Return: list(clid, method, status, canonical_label)
dict_lookup <- function(label, idx) {
  if (is.na(label) || !nzchar(label)) return(list(clid = NA_character_, method = "unmapped", status = "unmapped"))

  # 1. canonical exact
  if (exists(label, envir = idx$label_idx, inherits = FALSE)) {
    return(list(clid = get(label, envir = idx$label_idx), method = "canonical_exact", status = "mapped"))
  }
  # 2. synonym exact
  if (exists(label, envir = idx$syn_idx, inherits = FALSE)) {
    return(list(clid = get(label, envir = idx$syn_idx), method = "synonym_exact", status = "mapped"))
  }
  # 3. conservative normalization
  nk <- idx$norm_form(label)
  if (exists(nk, envir = idx$norm_idx, inherits = FALSE)) {
    return(list(clid = get(nk, envir = idx$norm_idx), method = "normalized", status = "mapped"))
  }
  # 4. trailing explanatory parenthesis (try the full label first, then strip one level only)
  no_paren <- sub("\\s*\\(.*\\)\\s*$", "", label)
  if (nzchar(no_paren) && no_paren != label) {
    nk_np <- idx$norm_form(no_paren)
    if (exists(nk_np, envir = idx$norm_idx, inherits = FALSE)) {
      return(list(clid = get(nk_np, envir = idx$norm_idx), method = "parenthetical_core", status = "mapped"))
    }
  }
  # 4b. strip explanatory parentheses (which may occur inside the label, "excitatory (glutamatergic) neurons").
  #     general lexical rule: remove parenthetical content and retry normalized matching; without label-specific hard-coding。
  all_noparen <- gsub("\\([^)]*\\)", " ", label, perl = TRUE)
  all_noparen <- stringr::str_squish(all_noparen)
  if (nzchar(all_noparen) && all_noparen != label) {
    nk_anp <- idx$norm_form(all_noparen)
    if (exists(nk_anp, envir = idx$norm_idx, inherits = FALSE)) {
      return(list(clid = get(nk_anp, envir = idx$norm_idx), method = "parenthetical_all", status = "mapped"))
    }
  }
  # 5. first-token fallback (first-token fallback label — try once without recursion)
  first <- strsplit(trimws(label), "\\s+")[[1]][1]
  if (!is.na(first) && nzchar(first) && first != label) {
    nk_first <- idx$norm_form(first)
    if (exists(nk_first, envir = idx$norm_idx, inherits = FALSE)) {
      return(list(clid = get(nk_first, envir = idx$norm_idx), method = "first_token", status = "mapped_low_conf"))
    }
  }
  list(clid = NA_character_, method = "unmapped", status = "unmapped")
}


# ---- Step 2: pure-R lexical top-k candidates (after deterministic mapping fails) ----
# use Jaro-Winkler and Levenshtein similarity to score the normalized label against all canonical CL labels
lexical_candidates <- function(label, idx, top_k = 5, min_score = 0.55) {
  if (!requireNamespace("stringdist", quietly = TRUE)) return(data.frame())
  nk <- idx$norm_form(label)
  if (!nzchar(nk)) return(data.frame())

  # collect normalized forms of all canonical labels
  all_labels <- ls(idx$norm_idx, all.names = TRUE)
  # compute Jaro-Winkler similarity (useful for short labels)
  jw <- stringdist::stringsim(nk, all_labels, method = "jw", p = 0.1)
  # normalize Levenshtein similarity
  lv <- stringdist::stringsim(nk, all_labels, method = "lv")
  score <- pmax(jw, lv)
  ord <- order(score, decreasing = TRUE)
  top <- head(ord, top_k * 3)  # retrieve extra candidates before filtering

  res <- lapply(top, function(i) {
    if (score[i] < min_score) return(NULL)
    clid <- get(all_labels[i], envir = idx$norm_idx)
    data.frame(
      cl_id = clid,
      canonical_label = idx$clid_info[[clid]]$canonical_label %||% all_labels[i],
      lexical_score = round(score[i], 3),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  # deduplicate (keep the highest score for each CL ID)
  df <- df[order(-df$lexical_score), ]
  df <- df[!duplicated(df$cl_id), ]
  head(df, top_k)
}

# ---- Main entry point: complete CL-Linker workflow ----
# Input: label (free text) + provided_clid ( NA)
# Output: complete mapping audit
cl_link <- function(label, provided_clid = NA_character_, idx) {
  raw_label <- as.character(label %||% NA_character_)
  if (is.na(raw_label) || !nzchar(raw_label)) {
    return(list(raw_label = NA_character_, parsed_core_label = NA_character_,
                canonical_label = NA_character_, cl_id = NA_character_,
                mapping_method = "empty", mapping_status = "unmapped",
                candidate_cl_ids = character(0), provided_cl_id_valid = FALSE))
  }

  # Layer 0: provided CL ID
  if (!is.na(provided_clid) && nzchar(provided_clid)) {
    v <- validate_provided_clid(provided_clid, idx, raw_label)
    if (v$ok) {
      return(list(raw_label = raw_label,
                  parsed_core_label = raw_label,
                  canonical_label = v$canonical_label %||% raw_label,
                  cl_id = provided_clid,
                  mapping_method = "provided_valid_clid",
                  mapping_status = "mapped",
                  candidate_cl_ids = provided_clid,
                  provided_cl_id_valid = TRUE))
    }
  }

  # Layer 1: dictionary
  d <- dict_lookup(raw_label, idx)
  if (is.na(d$clid)) {
    # general lexical fallback: strip explanatory parentheses and retry ( "excitatory (glutamatergic) neurons").
    stripped <- stringr::str_squish(gsub("\\([^)]*\\)", " ", raw_label, perl = TRUE))
    if (nzchar(stripped) && stripped != raw_label) {
      d2 <- dict_lookup(stripped, idx)
      if (!is.na(d2$clid)) {
        d <- d2
        d$method <- "parenthetical_all"
        d$status <- "mapped"
      } else {
        a2 <- alias_lookup(stripped, idx)
        if (!is.na(a2)) {
          d <- list(clid = a2, method = "alias_parenthetical_all", status = "mapped")
        }
      }
    }
  }
  if (is.na(d$clid)) {
    a <- alias_lookup(raw_label, idx)
    if (!is.na(a)) {
      d <- list(clid = a, method = "alias", status = "mapped")
    }
  }
  list(raw_label = raw_label,
       parsed_core_label = if (d$status == "unmapped") NA_character_ else raw_label,
       canonical_label = if (!is.na(d$clid)) idx$clid_info[[d$clid]]$canonical_label %||% raw_label else NA_character_,
       cl_id = d$clid,
       mapping_method = d$method,
       mapping_status = d$status,
       candidate_cl_ids = if (!is.na(d$clid)) d$clid else character(0),
       provided_cl_id_valid = FALSE)
}

# ---- Batch linking ----
cl_link_batch <- function(labels, provided_clids = NULL, idx) {
  if (is.null(provided_clids)) provided_clids <- rep(NA_character_, length(labels))
  lapply(seq_along(labels), function(i) {
    cl_link(labels[i], provided_clids[i], idx)
  })
}
