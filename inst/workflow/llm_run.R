###################################################################
# LLM INPUT PREPARATION AND VALIDATION UTILITIES
# Shared functions for structured LLM input preparation, validation and citation handling
###################################################################

# --- 0. Environment setup and package loading ---
suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(glue); library(stringr); library(readr); library(dplyr); library(purrr); library(memoise); library(cachem)
  library(clusterProfiler); library(org.Hs.eg.db); library(DOSE); library(ReactomePA); library(enrichR); library(data.table); library(tidyr)
  library(knitr)
  library(future); library(future.apply); library(tictoc)
  library(decoupleR); library(tibble); library(org.Hs.eg.db); library(AnnotationDbi)
})
# org.Mm.eg.db is a mouse-only annotation package; it is loaded
# conditionally (mouse analysis loads and validates it explicitly).

# DisGeNET disease evidence is OPTIONAL: it requires DISGENET_API_KEY and the
# disgenet2r package. Without the key the workflow degrades gracefully —
# the disease_disgenet evidence dimension is skipped with one informational
# message and every other evidence dimension is unchanged.
disgenet_key <- Sys.getenv("DISGENET_API_KEY", unset = "")
if (nzchar(disgenet_key)) {
  if (!requireNamespace("disgenet2r", quietly = TRUE)) {
    stop("llm_run: DISGENET_API_KEY is set but the 'disgenet2r' package is ",
         "missing. Install it with: remotes::install_gitlab(\"medbio/disgenet2r\")")
  }
  suppressPackageStartupMessages(library(disgenet2r))
} else {
  message("[INFO] DISGENET_API_KEY not set; DisGeNET disease evidence is skipped (all other evidence dimensions unchanged).")
}

# --- Project root for relative resource files (important for future workers) ---
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = Sys.getenv("TRIAGE_HOME", unset = getwd()))

# --- Optional CL normalization (for candidate CLID enrichment) ---
# CL normalization logic ships with the Triage package namespace
# (loaded by the sourcing scripts); CL_NORMALIZER_PATH is an optional
# override for advanced users only.
cl_normalizer_path <- Sys.getenv("CL_NORMALIZER_PATH", unset = "")
if (nzchar(cl_normalizer_path) && file.exists(cl_normalizer_path)) {
  source(cl_normalizer_path)
}

CL_LOCAL_JSON <- Sys.getenv("CL_LOCAL_JSON", unset = file.path(PROJECT_ROOT, "inputs", "raw", "ontology", "CL-ontology-v2025-07-30.json"))
CL_CFG_CACHE <- NULL
CL_GRAPH_CACHE <- NULL


get_cl_cfg <- function() {
  if (!is.null(CL_CFG_CACHE)) return(CL_CFG_CACHE)
  if (!nzchar(CL_LOCAL_JSON) || !file.exists(CL_LOCAL_JSON)) return(NULL)
  CL_CFG_CACHE <<- Triage:::make_cl_cfg(CL_LOCAL_JSON, prefer_ols = FALSE, cache_dir = "")
  CL_CFG_CACHE
}

get_cl_graph <- function() {
  if (!is.null(CL_GRAPH_CACHE)) return(CL_GRAPH_CACHE)
  if (!nzchar(CL_LOCAL_JSON) || !file.exists(CL_LOCAL_JSON)) return(NULL)
  CL_GRAPH_CACHE <<- jsonlite::fromJSON(CL_LOCAL_JSON, simplifyVector = FALSE)
  CL_GRAPH_CACHE
}

available_cores <- future::availableCores()
default_workers <- if (nzchar(Sys.getenv("TRIAGE_WORKERS", unset = ""))) as.integer(Sys.getenv("TRIAGE_WORKERS")) else max(1L, floor(available_cores / 2L))
plan(multisession, workers = default_workers)
cat(glue::glue("--- Parallel processing enabled (workers = {default_workers}) ---\n"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_primary_context_value <- function(x, fallback = NULL) {
  if (is.null(x) || length(x) == 0) return(fallback)
  x[[1]]
}

extract_terms_recursive <- function(x) {
  terms <- character(0)
  if (is.null(x)) return(terms)
  if (is.character(x)) {
    return(x)
  }
  if (is.list(x)) {
    if (!is.null(x$term)) terms <- c(terms, as.character(x$term))
    if (!is.null(x$name)) terms <- c(terms, as.character(x$name))
    if (!is.null(x$pathway)) terms <- c(terms, as.character(x$pathway))
    if (!is.null(x$description)) terms <- c(terms, as.character(x$description))
    for (v in x) {
      terms <- c(terms, extract_terms_recursive(v))
    }
  }
  unique(terms[!is.na(terms) & nzchar(terms)])
}

build_program_evidence <- function(bioinfo_results) {
  terms <- extract_terms_recursive(bioinfo_results)
  if (length(terms) == 0) return(list())
  terms_lc <- tolower(terms)
  pick <- function(pattern) {
    idx <- stringr::str_detect(terms_lc, pattern)
    unique(terms[idx])
  }
   programs <- list(
     epithelial_polarity_transport = pick("polarity|apical|basal|basolateral|tight junction|adherens junction|cell junction|transport|transporter|channel|junction|epithelial"),
     secretory_program = pick("secretory|secretion|vesicle|exocytosis|endocytosis|lysosome|phagocytosis"),
     basement_ecm_program = pick("basement membrane|extracellular matrix|ecm|cell adhesion|collagen|laminin|integrin|focal adhesion"),
     immune_activation = pick("immune|inflammatory|cytokine|interferon|antigen|leukocyte|innate|adaptive"),
     neural_differentiation = pick("neurogenesis|neuron differentiation|synapse|axon|dendrite|neuronal"),
     cell_cycle = pick("cell cycle|mitotic|g2/m|s phase|proliferation|checkpoint"),
     metabolic_program = pick("metabolic|mitochondrial|oxidative|respiratory|glycolysis|melanogenesis|pigment"),
     reproductive_germline_program = pick("spermatogenesis|sperm|acrosome|flagellum|gamete|germ cell|meiotic|meiosis|fertilization")
   )
  programs <- programs[vapply(programs, length, integer(1)) > 0]
  programs
}

build_program_evidence_from_markers <- function(degs_state_df, degs_df) {
  collect_genes <- function(df, n_max = 50L) {
    if (is.null(df)) return(character(0))
    df <- tibble::as_tibble(df)
    if (!"geneSymbol" %in% colnames(df)) return(character(0))
    genes <- as.character(df$geneSymbol)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    unique(utils::head(genes, n_max))
  }

  genes_state <- collect_genes(degs_state_df, n_max = 50L)
  genes_16_50 <- character(0)
  if (!is.null(degs_df) && "geneSymbol" %in% colnames(degs_df)) {
    genes <- as.character(degs_df$geneSymbol)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) >= 16) {
      genes_16_50 <- genes[16:min(50, length(genes))]
    }
  }

  genes <- unique(c(genes_state, genes_16_50))
  if (length(genes) == 0) return(list(programs = list(), support_level = list()))

  hits <- function(pattern) sum(grepl(pattern, genes, ignore.case = TRUE))
  support <- function(n) if (n >= 3) "strong" else if (n >= 2) "moderate" else if (n >= 1) "weak" else "none"

  epithelial_hits <- hits("^KRT|CLDN|OCLN|TJP|EPCAM|MUC|ITGA|ITGB|CDH|LAMA|COL")
  transport_hits <- hits("^SLC|^ABCA|^ABCB|^ABCC|^ABCD|ATP6V|VAMP|RAB|SNX")
  junction_hits <- hits("^GJA|^GJB|TJ|CLDN|OCLN|TJP|JAM")
  secretory_hits <- hits("SEC|COPA|COPB|COPG|EXOC|RAB|VAMP|STX|SNAP")

  programs <- list(
    epithelial_program_markers = genes[grepl("^KRT|CLDN|OCLN|TJP|EPCAM|MUC|ITGA|ITGB|CDH|LAMA|COL", genes, ignore.case = TRUE)],
    transport_program_markers = genes[grepl("^SLC|^ABCA|^ABCB|^ABCC|^ABCD|ATP6V|VAMP|RAB|SNX", genes, ignore.case = TRUE)],
    junction_program_markers = genes[grepl("^GJA|^GJB|TJ|CLDN|OCLN|TJP|JAM", genes, ignore.case = TRUE)],
    secretory_program_markers = genes[grepl("SEC|COPA|COPB|COPG|EXOC|RAB|VAMP|STX|SNAP", genes, ignore.case = TRUE)]
  )
  programs <- programs[vapply(programs, length, integer(1)) > 0]

  support_level <- list(
    epithelial = support(epithelial_hits + junction_hits),
    transport = support(transport_hits),
    secretory = support(secretory_hits)
  )

  list(programs = programs, support_level = support_level)
}

###################################################################
# SECTION 1: Five-dimension hybrid analysis module (MODIFIED)
###################################################################

# Explicit intermediate-output root (backward compatible). When
# TRIAGE_INTERMEDIATE_ROOT is set to an absolute path, enrichment TSVs are
# written under that root regardless of worker working directories; when it
# is unset, the historical working-directory-relative location is preserved.
resolve_bioinfo_dir <- function(run_name_prefix, masked_qid, intermediate_root = NULL) {
  # The explicit root is preferred: multisession workers are spawned when the
  # future plan is created, which happens BEFORE stage 05 exports
  # TRIAGE_INTERMEDIATE_ROOT, so the parent-resolved value must be forwarded
  # as a future global rather than read from the worker environment.
  explicit_root <- if (!is.null(intermediate_root) && nzchar(intermediate_root)) {
    intermediate_root
  } else {
    Sys.getenv("TRIAGE_INTERMEDIATE_ROOT", unset = "")
  }
  if (nzchar(explicit_root)) {
    file.path(explicit_root, run_name_prefix, "bioinformatics_tsv", masked_qid)
  } else {
    file.path("intermediate_outputs", run_name_prefix, "bioinformatics_tsv", masked_qid)
  }
}

save_enrichment_tsv <- function(result_obj, analysis_name, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  file_path <- file.path(output_dir, glue::glue("{analysis_name}_full_results.tsv"))
  
  if (is.null(result_obj)) {
    cat(glue::glue("   -> Skipping TSV for {analysis_name}: Result is NULL.\n"))
    return(invisible(NULL))
  }
  
  result_df <- NULL
  
  # Support enrichResult objects
  if (inherits(result_obj, "enrichResult")) {
    result_df <- as.data.frame(result_obj@result)
    
  } else if (isS4(result_obj) && "results" %in% slotNames(result_obj) && is.data.frame(result_obj@results)) {
    result_df <- as.data.frame(result_obj@results)
    
  } else if (is.data.frame(result_obj) || tibble::is_tibble(result_obj)) {
    result_df <- as.data.frame(result_obj)
    
  } else if (is.character(result_obj)) {
    result_df <- tibble::tibble(term = as.character(result_obj))
  }
  
  if (is.null(result_df) || nrow(result_df) == 0) {
    cat(glue::glue("   -> Skipping TSV for {analysis_name}: No data rows to write.\n"))
  } else {
    tryCatch({
      readr::write_tsv(result_df, file_path)
      cat(glue::glue("   -> Successfully saved full results for {analysis_name} to {file_path}\n"))
    }, error = function(e) {
      warning(glue::glue("Failed to write TSV for {analysis_name}. Error: {e$message}"))
    })
  }
}

get_species_resources <- function(species_name) {
  species_name <- tolower(species_name)
  if (species_name == "human") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    return(list(
      org_db = org.Hs.eg.db,
      kegg_organism = "hsa",
      reactome_organism = "human",
      gene_symbol_vocab = "HGNC"
    ))
  } else if (species_name == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Mouse analysis requires the 'org.Mm.eg.db' package. Please install it.")
    }
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    return(list(
      org_db = org.Mm.eg.db,
      kegg_organism = "mmu",
      reactome_organism = "mouse",
      gene_symbol_vocab = "MGI"
    ))
  } else {
    stop(glue::glue("Unsupported species: '{species_name}'. Supported species are 'human' and 'mouse'."))
  }
}

convert_genes_to_entrez <- function(gene_symbols, org_db) {
  if (is.null(gene_symbols) || length(gene_symbols) == 0) return(character(0))
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
  if (length(gene_symbols) == 0) return(character(0))
  mapIds(org_db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>%
    na.omit() %>% unique()
}

summarize_enrich <- function(res, n) {
  # The local Reactome implementation returns a character vector filtered and ordered by p value
  if (is.character(res)) {
    if (length(res) == 0) return(character(0))
    return(head(res, n))
  }
  if (is.null(res) || !inherits(res, "enrichResult") || nrow(res@result) == 0) return(character(0))
  significant_terms <- res@result %>% dplyr::filter(p.adjust < 0.05)
  if (nrow(significant_terms) == 0) return(character(0))
  significant_terms %>% head(n) %>% pull(Description)
}

analyze_ppi_from_local_file <- function(gene_list, ppi_files, n_hubs) {
  tryCatch({
    aliases_dt <- fread(ppi_files$aliases, col.names = c("id", "alias", "src"))
    symbol_to_string_id <- aliases_dt[src == "BioMart_HUGO", .(id, alias)]
    input_string_ids <- symbol_to_string_id[alias %in% gene_list, id]
    if (length(input_string_ids) < 2) return("Fewer than 2 genes mapped.")
    
    interactions_dt <- fread(ppi_files$interactions, col.names = c("p1", "p2", "score"))
    subnet_dt <- interactions_dt[p1 %in% input_string_ids & p2 %in% input_string_ids]
    if (nrow(subnet_dt) == 0) return("No interactions found.")
    
    degree_counts <- as.data.table(table(c(subnet_dt$p1, subnet_dt$p2)))
    setnames(degree_counts, c("id", "degree"))
    hub_gene_names <- degree_counts %>%
      arrange(desc(degree)) %>%
      head(n_hubs) %>%
      left_join(symbol_to_string_id, by = "id") %>%
      pull(alias) %>%
      na.omit() %>%
      paste(collapse = ", ")
    
    if (nchar(hub_gene_names) > 0) paste("Top hub genes:", hub_gene_names) else "No hubs identified."
  }, error = function(e) {
    warning("Local PPI failed: ", e$message)
    "Local PPI analysis failed."
  })
}

# ✅ DisGeNET (FREE/ACADEMIC SAFE): avoid disease_enrichment endpoint (often 403 on ACADEMIC)
# Use gene2disease(database="CURATED") per gene + local aggregation.

# ---- Local Reactome enrichment(replacing the network-dependent ReactomePA::enrichPathway call, <6s)----
enrich_reactome_local <- function(gene_entrez, pvalueCutoff = 0.05) {
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    warning("reactome.db is unavailable; skipping Reactome enrichment")
    return(character(0))
  }
  suppressPackageStartupMessages(library(reactome.db))
  path2ext <- as.list(reactomePATHID2EXTID)
  path2name <- as.list(reactomePATHID2NAME)
  bg <- unique(unlist(path2ext))
  n_bg <- length(bg); n_sel <- length(gene_entrez)
  ids <- unique(as.character(gene_entrez))
  res <- lapply(names(path2ext), function(pid) {
    pgenes <- path2ext[[pid]]
    overlap <- sum(ids %in% pgenes)
    if (overlap == 0) return(NULL)
    p <- phyper(overlap - 1, length(pgenes), n_bg - length(pgenes), n_sel, lower.tail = FALSE)
    list(pathway = unname(path2name[[pid]]), overlap = overlap, size = length(pgenes), pvalue = p)
  })
  res <- do.call(rbind, lapply(res, function(x) if (is.null(x)) NULL else as.data.frame(x, stringsAsFactors = FALSE)))
  if (is.null(res) || nrow(res) == 0) return(character(0))
  res <- res[order(res$pvalue), ]
  sig <- res[res$pvalue < pvalueCutoff, ]
  if (nrow(sig) == 0) return(character(0))
  gsub("^Homo sapiens: ", "", sig$pathway)
}

analyze_disgenet_enrichment <- function(gene_list, n, vocabulary = "HGNC") {
  api_key <- Sys.getenv("DISGENET_API_KEY", unset = NA_character_)
  
  if (is.na(api_key) || nchar(api_key) < 30) {
    return(character(0))
  }
  if (length(gene_list) < 2) return(character(0))
  
  db_to_use <- "CURATED"
  
  g2d_list <- lapply(unique(gene_list), function(g) {
    suppressWarnings(
      tryCatch(
        withCallingHandlers(
          disgenet2r::gene2disease(
            gene       = g,
            vocabulary = vocabulary,
            database   = db_to_use
          ),
          warning = function(w) {
            if (grepl("timeout|Timeout|HTTP|network|connection", conditionMessage(w), ignore.case = TRUE)) {
              invokeRestart("muffleWarning")
            }
          }
        ),
        error = function(e) NULL
      )
    )
  })
  
  g2d_df <- suppressWarnings(dplyr::bind_rows(g2d_list))
  if (is.null(g2d_df) || nrow(g2d_df) == 0) return(character(0))
  
  disease_name_col <- intersect(
    c("disease_name", "diseaseName", "disease_name.x", "disease"),
    names(g2d_df)
  )
  disease_id_col <- intersect(
    c("diseaseid", "diseaseId", "disease_id"),
    names(g2d_df)
  )
  
  col_to_use <- if (length(disease_name_col) > 0) disease_name_col[1] else
    if (length(disease_id_col) > 0) disease_id_col[1] else NA_character_
  
  if (is.na(col_to_use)) return(character(0))
  
  out <- g2d_df %>%
    dplyr::mutate(.disease = .data[[col_to_use]]) %>%
    dplyr::filter(!is.na(.disease), .disease != "") %>%
    dplyr::count(.disease, name = "n_genes_supporting") %>%
    dplyr::arrange(dplyr::desc(n_genes_supporting)) %>%
    dplyr::slice_head(n = as.integer(n)) %>%
    dplyr::pull(.disease)
  
  if (length(out) == 0) character(0) else as.character(out)
}

run_decoupleR_gsea <- function(marker_df, organism = "human", n_tfs = 15) {
  tryCatch({
    if (!all(c("gene", "avg_log2FC") %in% names(marker_df))) {
      stop("Marker data frame must contain 'gene' and 'avg_log2FC' columns.")
    }
    
    net_filename <- switch(
      organism,
      human = file.path(PROJECT_ROOT, "inputs", "raw", "collectri", "collectri_human_network.rds"),
      mouse = file.path(PROJECT_ROOT, "inputs", "raw", "collectri", "collectri_mouse_network.rds"),
      stop(paste0("Unsupported organism '", organism, "'."))
    )
    
    net_path <- if (grepl("^(/|[A-Za-z]:)", net_filename)) {
      normalizePath(net_filename, mustWork = FALSE)
    } else {
      normalizePath(file.path(PROJECT_ROOT, net_filename), mustWork = FALSE)
    }
    
    if (!file.exists(net_path)) stop(paste("Local network file not found:", net_path))
    net <- readRDS(net_path)
    
    ranked_genes <- marker_df %>%
      dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::select(gene, avg_log2FC) %>%
      tibble::deframe()
    
    tf_activities <- decoupleR::run_fgsea(mat = ranked_genes, net = net, minsize = 5)
    
    activated_tfs <- tf_activities %>%
      dplyr::filter(score > 0) %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      head(n_tfs) %>%
      pull(source)
    
    inhibited_tfs <- tf_activities %>%
      dplyr::filter(score < 0) %>%
      dplyr::arrange(score) %>%
      head(n_tfs) %>%
      pull(source)
    
    list(
      activated = activated_tfs,
      inhibited = inhibited_tfs,
      raw_activities = tf_activities
    )
  }, error = function(e) {
    warning("decoupleR GSEA failed: ", e$message)
    list(activated = character(0), inhibited = character(0), raw_activities = tibble::tibble())
  })
}

###################################################################
# <<< MODIFIED: SETUP_DYNAMIC_CACHES >>>
###################################################################
setup_dynamic_caches <- function(dataset_name, config, dimension_definitions, ppi_files, species_resources) {
  main_cache_dir <- file.path(Sys.getenv("TRIAGE_CACHE_ROOT", unset = "api_caches"), dataset_name)
  if (!dir.exists(main_cache_dir)) dir.create(main_cache_dir, recursive = TRUE)
  
  perform_bioinformatics_analysis_impl <- function(marker_df_key, config, ppi_files, prompt_context_config, resources) {
    if (isTRUE(config$cache_only)) return(list())
    marker_df <- unserialize(marker_df_key)
    
    gene_list_for_other_analyses <- marker_df %>%
      head(config$n_deg_for_other_analyses) %>%
      pull(gene)
    
    entrez_ids <- convert_genes_to_entrez(gene_list_for_other_analyses, resources$org_db)
    
    if (length(entrez_ids) < 5) {
      warning("Not enough genes with Entrez IDs for analysis.")
      return(list())
    }
    
    run_enrich <- function(analysis_func, ...) {
      tryCatch({
        res <- suppressMessages(analysis_func(...))
        if (inherits(res, "enrichResult") && "ONTOLOGY" %in% slotNames(res) && res@ontology %in% c("BP", "CC", "MF")) {
          res <- suppressMessages(clusterProfiler::simplify(res, cutoff = 0.7))
        }
        res
      }, error = function(e) {
        warning(paste("Enrichment analysis failed:", e$message))
        NULL
      })
    }
    
  current_species <- get_primary_context_value(prompt_context_config$species, "human")
    disgenet_entities_to_use <- gene_list_for_other_analyses
    disgenet_vocabulary_to_use <- "HGNC"
    
    if (current_species == "mouse") {
      cat("  -> [DisGeNET] Species is mouse, converting gene symbols to Entrez IDs for query.\n")
      disgenet_entities_to_use <- entrez_ids
      disgenet_vocabulary_to_use <- "ENTREZ"
    }
    
    list(
      go_bp = run_enrich(clusterProfiler::enrichGO, gene = entrez_ids, ont = "BP", pvalueCutoff = 0.05, OrgDb = resources$org_db),
      go_cc = run_enrich(clusterProfiler::enrichGO, gene = entrez_ids, ont = "CC", pvalueCutoff = 0.05, OrgDb = resources$org_db),
      go_mf = run_enrich(clusterProfiler::enrichGO, gene = entrez_ids, ont = "MF", pvalueCutoff = 0.05, OrgDb = resources$org_db),
      kegg = run_enrich(clusterProfiler::enrichKEGG, gene = entrez_ids, organism = resources$kegg_organism, pvalueCutoff = 0.05, use_internal_data = TRUE),
      reactome = tryCatch(enrich_reactome_local(entrez_ids, 0.05), error = function(e) character(0)),
      disease_do = run_enrich(DOSE::enrichDO, gene = entrez_ids, pvalueCutoff = 0.05),
      
      decoupleR_gsea = tryCatch({
        run_decoupleR_gsea(marker_df, organism = current_species, n_tfs = config$n_tf)
      }, error = function(e) NULL),
      
      ppi_hubs = analyze_ppi_from_local_file(gene_list_for_other_analyses, ppi_files, config$n_hubs),
      
      disease_disgenet = tryCatch({
        analyze_disgenet_enrichment(
          gene_list = disgenet_entities_to_use,
          n = config$n_disease_each,
          vocabulary = disgenet_vocabulary_to_use
        )
      }, error = function(e) character(0))
    )
  }
  
  query_litsense_api_impl <- function(gene_set_key, species, tissue, study_context, keywords_key, config) {
    if (isTRUE(config$cache_only)) return(tibble::tibble())
    suppressPackageStartupMessages({
      library(httr); library(jsonlite); library(glue); library(stringr); library(dplyr); library(tibble)
    })
    
    gene_set <- unlist(stringr::str_split(gene_set_key, ","))
    keywords <- unlist(stringr::str_split(keywords_key, ","))
    
    if (length(gene_set) == 0) return(tibble::tibble())
    
    all_terms <- c(gene_set, species, tissue, keywords)
    final_query_logic <- paste(all_terms, collapse = " ")
    
    encoded_query_string <- URLencode(final_query_logic)
    base_url <- "https://www.ncbi.nlm.nih.gov/research/litsense-api/api/"
    full_url <- paste0(base_url, "?query=", encoded_query_string, "&rerank=true")
    
  on.exit({ cat("   - Pausing for 1.1s to respect API rate limit.\n"); Sys.sleep(1.1) })
    
    tryCatch({
      cat(glue::glue("Calling LitSense API for dimension '{keywords_key}'\n"))
      cat(glue::glue("  -> [USER LOGIC: AND] Final URL: {full_url}\n"))
      
      request_config <- httr::config(connecttimeout = 20, timeout = 60)
      response <- httr::RETRY(
        "GET", url = full_url, config = request_config, times = 7,
        pause_base = 1, pause_cap = 60, pause_min = 1,
        terminate_on = c(400, 401, 403, 404), quiet = FALSE
      )
      
      if (httr::status_code(response) != 200) {
    warning(glue::glue("LitSense failed with status {httr::status_code(response)}"))
        return(tibble::tibble())
      }
      
      content_text <- httr::content(response, "text", encoding = "UTF-8")
      if (!jsonlite::validate(content_text)) {
        warning("LitSense returned non-JSON content.")
        return(tibble::tibble())
      }
      
      parsed <- jsonlite::fromJSON(content_text)
      
      if (is.null(parsed) || !is.data.frame(parsed) || nrow(parsed) == 0 || !"score" %in% names(parsed)) {
        return(tibble::tibble())
      }
      
      primary_threshold <- config$litsense_primary_threshold %||% 0.6
      fallback_threshold <- config$litsense_fallback_threshold %||% 0.5
      
      high_score_results <- parsed %>% tibble::as_tibble() %>% dplyr::filter(score >= primary_threshold)
      
      if (nrow(high_score_results) > 0) {
        cat(glue::glue("  -> Success! Found {nrow(high_score_results)} articles with score >= {primary_threshold}\n"))
        return(high_score_results)
      } else {
        cat(glue::glue("  -> No articles found with score >= {primary_threshold}. Attempting fallback with threshold >= {fallback_threshold}...\n"))
        fallback_results <- parsed %>% tibble::as_tibble() %>% dplyr::filter(score >= fallback_threshold)
        
        if (nrow(fallback_results) > 0) {
          cat(glue::glue("  -> Fallback Success! Found {nrow(fallback_results)} articles with score >= {fallback_threshold}\n"))
          return(fallback_results)
        } else {
          cat("  -> Fallback failed. No articles found above the lower threshold either.\n")
          return(tibble::tibble())
        }
      }
    }, error = function(e) {
      warning(paste("LitSense API call failed:", e$message))
      tibble::tibble()
    })
  }

  # =========================================================
  # Europe PMC PRIMARY: search -> cache fullTextXML -> snippets
  # Output schema aligned to LitSense: pmid/pmcid/title/text/score
  # =========================================================
  eupmc_safe_get_text <- function(url, timeout_sec = 60) {
    res <- httr::GET(url, httr::timeout(timeout_sec))
    httr::stop_for_status(res)
    httr::content(res, as = "text", encoding = "UTF-8")
  }

  eupmc_urlenc <- function(x) URLencode(x, reserved = TRUE)

  eupmc_normalize_pmcid <- function(x) {
    x <- as.character(x %||% "")
    m <- stringr::str_extract(x, "PMC\\d+")
    ifelse(is.na(m), "", m)
  }

  eupmc_build_query <- function(gene_set, species, tissue, data_type, keywords) {
    gene_set <- gene_set[!is.na(gene_set) & nzchar(gene_set)]
    keywords <- keywords[!is.na(keywords) & nzchar(keywords)]
    if (length(gene_set) == 0) return("")
    gene_clause <- glue::glue("(({paste0('TITLE_ABS:', gene_set) %>% paste(collapse=' OR ')}) OR ({paste0('BODY:', gene_set) %>% paste(collapse=' OR ')}))")
    sc_terms <- c(
      '"single cell"', '"single-cell"', 'scRNAseq', '"scRNA-seq"',
      '"single nucleus"', '"single-nucleus"', 'snRNAseq', '"snRNA-seq"',
      'transcriptome', 'organoid'
    )
    sc_clause <- glue::glue("({paste(sc_terms, collapse=' OR ')})")
    ctx_terms <- c(species, tissue, data_type, keywords)
    ctx_terms <- ctx_terms[!is.na(ctx_terms) & nzchar(ctx_terms)]
    ctx_clause <- if (length(ctx_terms) > 0) glue::glue("({paste(shQuote(ctx_terms), collapse=' OR ')})") else ""
    pieces <- c(gene_clause, sc_clause, ctx_clause, "sort_date:y")
    paste(pieces[pieces != ""], collapse = " AND ")
  }

  eupmc_search <- function(query, page_size = 30, result_type = "core", timeout_sec = 60) {
    base <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    url <- glue::glue("{base}?query={eupmc_urlenc(query)}&format=json&resultType={result_type}&pageSize={page_size}")
    txt <- eupmc_safe_get_text(url, timeout_sec = timeout_sec)
    dat <- jsonlite::fromJSON(txt, flatten = TRUE)
    res_raw <- dat$resultList$result
    if (is.null(res_raw)) return(tibble::tibble())
    res <- tibble::as_tibble(res_raw)
    if (nrow(res) == 0) return(tibble::tibble())
    getcol <- function(df, col, default="") if (col %in% names(df)) df[[col]] else rep(default, nrow(df))
    tibble::tibble(
      rank  = seq_len(nrow(res)),
      pmid  = as.character(getcol(res, "pmid", "")),
      pmcid = as.character(getcol(res, "pmcid", "")),
      doi   = as.character(getcol(res, "doi", "")),
      title = as.character(getcol(res, "title", "")),
      year  = as.character(getcol(res, "pubYear", "")),
      abstract = as.character(getcol(res, "abstractText", ""))
    ) %>%
      dplyr::mutate(pmcid_norm = eupmc_normalize_pmcid(pmcid)) %>%
      dplyr::distinct(pmid, pmcid_norm, .keep_all = TRUE) %>%
      dplyr::arrange(rank)
  }

  eupmc_fulltext_xml_url <- function(pmcid_norm) {
    glue::glue("https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid_norm}/fullTextXML")
  }

  eupmc_cache_dir <- function(main_cache_dir) {
    d <- file.path(main_cache_dir, "eupmc_fulltext_cache")
    fs::dir_create(d, recurse = TRUE)
    d
  }

  eupmc_cache_one <- function(main_cache_dir, pmcid_norm, pmid="", doi="", query_used="", force=FALSE, timeout_sec=60) {
    cache_dir <- eupmc_cache_dir(main_cache_dir)
    xml_path  <- file.path(cache_dir, glue::glue("{pmcid_norm}.fullTextXML.xml"))
    meta_path <- file.path(cache_dir, glue::glue("{pmcid_norm}.meta.json"))
    if (!force && file.exists(xml_path)) {
      return(list(ok=TRUE, cached=TRUE, xml_path=xml_path, meta_path=meta_path))
    }
    url <- eupmc_fulltext_xml_url(pmcid_norm)
    xml_text <- tryCatch(eupmc_safe_get_text(url, timeout_sec = timeout_sec), error=function(e) NULL)
    if (is.null(xml_text) || nchar(xml_text) < 200) {
      return(list(ok=FALSE, cached=FALSE, reason="fullTextXML fetch failed/empty", url=url))
    }
    writeLines(xml_text, xml_path, useBytes = TRUE)
    jsonlite::write_json(list(
      pmcid = pmcid_norm, pmid = pmid, doi = doi,
      fetched_from = url,
      fetched_at_utc = format(Sys.time(), tz="UTC", usetz=TRUE),
      query_used = query_used
    ), meta_path, auto_unbox = TRUE, pretty = TRUE)
    list(ok=TRUE, cached=FALSE, xml_path=xml_path, meta_path=meta_path)
  }

  eupmc_extract_blocks <- function(xml_text) {
    doc <- xml2::read_xml(xml_text)
    xml2::xml_ns_strip(doc)
    xpaths <- list(
      abstract = ".//abstract//p",
      body = ".//body//p",
      fig = ".//fig//caption//p",
      table_caption = ".//table-wrap//caption//p",
      table_foot = ".//table-wrap-foot//p"
    )
    purrr::imap_dfr(xpaths, function(xp, sec) {
      nodes <- as.list(xml2::xml_find_all(doc, xp))
      if (length(nodes) == 0) return(tibble::tibble(section=character(), text=character()))
      tibble::tibble(
        section = sec,
        text = purrr::map_chr(nodes, ~ xml2::xml_text(.x, trim = TRUE))
      )
    }) %>%
      dplyr::mutate(text = stringr::str_replace_all(text, "\\s+", " ")) %>%
      dplyr::filter(nchar(text) > 40) %>%
      dplyr::distinct(section, text)
  }

  eupmc_split_sentences <- function(text) {
    s <- stringr::str_split(text, "(?<=[\\.\\?\\!])\\s+(?=[A-Z0-9\\(])", simplify = FALSE)[[1]]
    s <- s[!is.na(s) & nchar(stringr::str_trim(s)) > 0]
    s
  }

  eupmc_marker_hits <- function(text, markers) {
    markers <- markers[!is.na(markers) & markers != ""]
    if (length(markers) == 0) return(character(0))
    up <- stringr::str_to_upper(text)
    unique(purrr::keep(markers, function(m) {
      stringr::str_detect(up, glue::glue("\\b{stringr::str_to_upper(m)}\\b"))
    }))
  }

  eupmc_snippets_from_xml <- function(xml_path, pmid, pmcid, title, gene_markers, dim_markers, topk=10, require_gene_hit=TRUE) {
    xml_text <- tryCatch(paste(readLines(xml_path, warn = FALSE), collapse = "\n"),
                         error = function(e) NA_character_)
    if (is.na(xml_text) || nchar(xml_text) < 200) return(tibble::tibble())
    blocks <- eupmc_extract_blocks(xml_text)
    if (nrow(blocks) == 0) return(tibble::tibble())
    sents <- blocks %>%
      dplyr::mutate(sent = purrr::map(text, eupmc_split_sentences)) %>%
      tidyr::unnest(sent) %>%
      dplyr::transmute(section = section, text = stringr::str_replace_all(sent, "\\s+", " ")) %>%
      dplyr::filter(nchar(text) >= 60) %>%
      dplyr::distinct(section, text)
    if (nrow(sents) == 0) return(tibble::tibble())
    scored <- sents %>%
      dplyr::mutate(
        gene_hits = purrr::map(text, ~ eupmc_marker_hits(.x, gene_markers)),
        dim_hits  = purrr::map(text, ~ eupmc_marker_hits(.x, dim_markers)),
        n_gene = purrr::map_int(gene_hits, length),
        n_dim  = purrr::map_int(dim_hits, length),
        sec_bonus = dplyr::case_when(
          section %in% c("fig","table_caption","table_foot") ~ 2L,
          section == "abstract" ~ 1L,
          TRUE ~ 0L
        ),
        score = n_gene + 0.5 * n_dim + sec_bonus
      ) %>%
      { if (isTRUE(require_gene_hit)) dplyr::filter(., n_gene >= 1) else dplyr::filter(., n_gene >= 1 | n_dim >= 1) } %>%
      dplyr::arrange(dplyr::desc(score), dplyr::desc(n_gene), dplyr::desc(n_dim)) %>%
      dplyr::slice_head(n = topk) %>%
      dplyr::transmute(
        pmid = as.character(pmid),
        pmcid = as.character(pmcid),
        title = as.character(title %||% ""),
        text = as.character(text),
        score = as.numeric(score),
        annotations = gene_hits
      )
    scored
  }

  eupmc_snippets_from_abstract <- function(pmid, pmcid, title, abstract, gene_markers, dim_markers, topk=5, require_gene_hit=TRUE) {
    text <- paste(title %||% "", abstract %||% "")
    text <- stringr::str_replace_all(text, "\\s+", " ")
    if (!nzchar(text) || nchar(text) < 40) return(tibble::tibble())
    sents <- tibble::tibble(section = "abstract", text = eupmc_split_sentences(text)) %>%
      dplyr::filter(nchar(text) >= 60) %>%
      dplyr::distinct(section, text)
    if (nrow(sents) == 0) return(tibble::tibble())
    scored <- sents %>%
      dplyr::mutate(
        gene_hits = purrr::map(text, ~ eupmc_marker_hits(.x, gene_markers)),
        dim_hits  = purrr::map(text, ~ eupmc_marker_hits(.x, dim_markers)),
        n_gene = purrr::map_int(gene_hits, length),
        n_dim  = purrr::map_int(dim_hits, length),
        sec_bonus = 1L,
        score = n_gene + 0.5 * n_dim + sec_bonus
      ) %>%
      { if (isTRUE(require_gene_hit)) dplyr::filter(., n_gene >= 1) else dplyr::filter(., n_gene >= 1 | n_dim >= 1) } %>%
      dplyr::arrange(dplyr::desc(score), dplyr::desc(n_gene), dplyr::desc(n_dim)) %>%
      dplyr::slice_head(n = topk) %>%
      dplyr::transmute(
        pmid = as.character(pmid),
        pmcid = as.character(pmcid),
        title = as.character(title %||% ""),
        text = as.character(text),
        score = as.numeric(score),
        annotations = gene_hits
      )
    scored
  }

  query_eupmc_primary_impl <- function(gene_set_key, species, tissue, data_type, study_context, keywords_key, config, main_cache_dir) {
    if (isTRUE(config$cache_only)) return(tibble::tibble())
    if (!requireNamespace("xml2", quietly = TRUE)) stop("Need xml2")
    if (!requireNamespace("furrr", quietly = TRUE)) stop("Need furrr")
    suppressPackageStartupMessages({
      library(dplyr); library(tibble); library(purrr); library(stringr); library(tidyr); library(glue); library(fs); library(xml2); library(future); library(furrr)
    })
    gene_set <- unlist(stringr::str_split(gene_set_key, ","))
    keywords <- unlist(stringr::str_split(keywords_key, ","))
    gene_set <- gene_set[!is.na(gene_set) & nzchar(gene_set)]
    keywords <- keywords[!is.na(keywords) & nzchar(keywords)]
    if (length(gene_set) == 0) return(tibble::tibble())
    q <- eupmc_build_query(
      gene_set = gene_set,
      species = species %||% "",
      tissue  = tissue %||% "",
      data_type = data_type %||% "",
      keywords = keywords
    )
    if (is.na(q) || !nzchar(q)) return(tibble::tibble())
    max_hits   <- as.integer(config$eupmc_max_hits   %||% 30)
    max_papers <- as.integer(config$eupmc_max_papers %||% 8)
    topk_snips <- as.integer(config$eupmc_topk_snips %||% 10)
    timeout_sec <- as.integer(config$eupmc_timeout_sec %||% 60)
    n_dl <- as.integer(config$eupmc_download_workers %||% 8)
    n_parse <- as.integer(config$eupmc_parse_workers %||% 4)

    cat(glue::glue("Calling Europe PMC (PRIMARY) for dimension '{keywords_key}'\n"))
    cat(glue::glue("  -> Query: {q}\n"))
    hits <- eupmc_search(q, page_size = max_hits, result_type = "core", timeout_sec = timeout_sec)
    if (is.null(hits) || !is.data.frame(hits) || nrow(hits) == 0) return(tibble::tibble())
    hits <- hits %>% dplyr::slice_head(n = max_papers)
    hits_with_pmcid <- hits %>% dplyr::filter(nzchar(pmcid_norm))
    hits_no_pmcid <- hits %>% dplyr::filter(!nzchar(pmcid_norm))

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_dl)
    dl <- furrr::future_pmap(
      hits_with_pmcid,
      function(rank, pmid, pmcid, doi, title, year, pmcid_norm) {
        c(
          list(pmid = pmid, pmcid_norm = pmcid_norm, title = title),
          eupmc_cache_one(main_cache_dir, pmcid_norm, pmid = pmid, doi = doi, query_used = q,
                          force = FALSE, timeout_sec = timeout_sec)
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )

    ok <- purrr::keep(dl, ~ isTRUE(.x$ok) && nzchar(.x$xml_path))

    future::plan(future::multisession, workers = n_parse)
    sn_fulltext <- tibble::tibble()
    if (length(ok) > 0) {
      sn_fulltext <- furrr::future_map_dfr(
        ok,
        function(x) {
          eupmc_snippets_from_xml(
            xml_path = x$xml_path,
            pmid = x$pmid %||% "",
            pmcid = x$pmcid_norm %||% "",
            title = x$title %||% "",
            gene_markers = gene_set,
            dim_markers = keywords,
            topk = topk_snips,
            require_gene_hit = isTRUE(config$eupmc_require_gene_hit)
          )
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    }
    sn_abstract <- tibble::tibble()
    if (nrow(hits_no_pmcid) > 0) {
      sn_abstract <- hits_no_pmcid %>%
        dplyr::mutate(snippets = purrr::pmap(
          list(pmid, pmcid, title, abstract),
          ~ eupmc_snippets_from_abstract(..1, ..2, ..3, ..4, gene_set, keywords, topk = 3,
                                         require_gene_hit = isTRUE(config$eupmc_require_gene_hit))
        )) %>%
        dplyr::select(snippets) %>%
        tidyr::unnest(snippets)
    }
    dplyr::bind_rows(sn_fulltext, sn_abstract)
  }
  
  list(
    perform_bioinformatics_analysis = memoise::memoise(
      purrr::partial(perform_bioinformatics_analysis_impl, resources = species_resources),
      cache = cachem::cache_disk(file.path(main_cache_dir, "bioinfo_cache"))
    ),
    query_eupmc_primary_cached = memoise::memoise(
      purrr::partial(query_eupmc_primary_impl, main_cache_dir = main_cache_dir),
      cache = cachem::cache_disk(file.path(main_cache_dir, "eupmc_primary_cache"))
    ),
    query_litsense_api_cached = memoise::memoise(
      query_litsense_api_impl,
      cache = cachem::cache_disk(file.path(main_cache_dir, "litsense_cache"))
    ),
    generate_smart_context = function(...) ""
  )
}

get_evidence_partitioned <- function(full_gene_set_key, context, config, dimension_definitions, cached_functions, prompt_context_config) {
  full_gene_set <- unlist(stringr::str_split(full_gene_set_key, ",")) %>% unique() %>% sort()
  cat(glue::glue("--- Retrieving literature evidence for a STABILIZED gene set of {length(full_gene_set)} genes...\n"))
  
  enhanced_context <- c(context, unlist(prompt_context_config))
  min_articles <- 3
  dimension_names <- names(dimension_definitions)
  if (length(dimension_names) == 0) {
    cat("  -> WARNING: No dimension definitions. Literature evidence will be empty.\n")
    return(list(
      articles_db = list(pmid = character(0), text = character(0), score = numeric(0)),
      relevance_map = list()
    ))
  }
  
  sanitize_litsense_output <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
      return(tibble::tibble(pmid = character(), title = character(), text = character(), score = numeric()))
    expected_cols <- c("pmid", "pmcid", "title", "text", "score", "section", "annotations")
    for (col in expected_cols) if (!col %in% names(df)) df[[col]] <- NA
    df %>%
      dplyr::mutate(
        pmid = as.character(pmid),
        title = as.character(title),
        text = as.character(text),
        score = as.numeric(score)
      ) %>%
      dplyr::filter(!is.na(pmid), nchar(pmid) > 0) %>%
      dplyr::select(any_of(c("pmid", "pmcid", "title", "text", "score", "section", "annotations"))) %>%
      tibble::as_tibble()
  }
  
  n_candidates_to_fetch <- 10
  
  retrieve_candidates <- function(genes, ctx, source = c("eupmc", "litsense"), pctx = prompt_context_config) {
    source <- match.arg(source)
    sorted_genes <- sort(unique(genes))
    stable_chunk_key <- paste(sorted_genes, collapse = ",")
    if (length(dimension_names) == 0) return(list())
    lapply(dimension_names, function(d) {
      keywords_key <- paste(dimension_definitions[[d]], collapse = ",")
      raw_df <- if (source == "eupmc") {
        cached_functions$query_eupmc_primary_cached(
          gene_set_key = stable_chunk_key,
          species = get_primary_context_value(pctx$species, "human"),
          tissue  = get_primary_context_value(pctx$tissue, ""),
          data_type = get_primary_context_value(pctx$data_type, "single-cell RNA"),
          study_context = get_primary_context_value(pctx$study_context, ""),
          keywords_key = keywords_key,
          config = config
        )
      } else {
        cached_functions$query_litsense_api_cached(
          gene_set_key = stable_chunk_key,
          species = get_primary_context_value(pctx$species, "human"),
          tissue = get_primary_context_value(pctx$tissue, ""),
          study_context = get_primary_context_value(pctx$study_context, ""),
          keywords_key = keywords_key,
          config = config
        )
      }
      safe_df <- sanitize_litsense_output(raw_df)
      safe_df %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::slice_head(n = n_candidates_to_fetch) %>%
        dplyr::select(any_of(c("pmid", "pmcid", "title", "text", "score")))
    }) %>% stats::setNames(dimension_names)
  }

  has_any_rows <- function(x) {
    any(vapply(x, function(df) is.data.frame(df) && nrow(df) > 0, logical(1)))
  }
  
  retrieval_tier <- "strict"
  safe_pull <- function(x, d) {
    if (is.null(x) || !is.list(x)) return(tibble::tibble())
    if (!d %in% names(x)) return(tibble::tibble())
    x[[d]]
  }
  if (length(full_gene_set) > config$genes_per_sample) {
    chunks <- split(full_gene_set, ceiling(seq_along(full_gene_set) / config$genes_per_sample))
    all_evidence_list <- lapply(chunks, function(chunk) retrieve_candidates(chunk, enhanced_context, source = "eupmc"))
    initial_retrieval <- purrr::map(dimension_names, function(d) {
      all_evidence_list %>%
        purrr::map_dfr(~ safe_pull(.x, d)) %>%
        dplyr::distinct(pmid, .keep_all = TRUE) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::slice_head(n = n_candidates_to_fetch)
    }) %>% stats::setNames(dimension_names)
  } else {
    initial_retrieval <- retrieve_candidates(full_gene_set, enhanced_context, source = "eupmc")
  }

  total_rows <- sum(vapply(initial_retrieval, function(df) if (is.data.frame(df)) nrow(df) else 0L, integer(1)))
  if (total_rows < min_articles) {
    retrieval_tier <- "relaxed"
    relaxed_pctx <- prompt_context_config
    relaxed_pctx$study_context <- ""
    relaxed_pctx$user_notes <- ""
    if (length(full_gene_set) > config$genes_per_sample) {
      all_evidence_list <- lapply(chunks, function(chunk) retrieve_candidates(chunk, enhanced_context, source = "eupmc", pctx = relaxed_pctx))
      initial_retrieval <- purrr::map(dimension_names, function(d) {
        all_evidence_list %>%
          purrr::map_dfr(~ safe_pull(.x, d)) %>%
          dplyr::distinct(pmid, .keep_all = TRUE) %>%
          dplyr::arrange(dplyr::desc(score)) %>%
          dplyr::slice_head(n = n_candidates_to_fetch)
      }) %>% stats::setNames(dimension_names)
    } else {
      initial_retrieval <- retrieve_candidates(full_gene_set, enhanced_context, source = "eupmc", pctx = relaxed_pctx)
    }
    total_rows <- sum(vapply(initial_retrieval, function(df) if (is.data.frame(df)) nrow(df) else 0L, integer(1)))
    if (total_rows < min_articles) {
      retrieval_tier <- "litsense"
      cat("  -> Europe PMC insufficient. Falling back to LitSense...\n")
      if (length(full_gene_set) > config$genes_per_sample) {
        all_evidence_list <- lapply(chunks, function(chunk) retrieve_candidates(chunk, enhanced_context, source = "litsense", pctx = relaxed_pctx))
        initial_retrieval <- purrr::map(dimension_names, function(d) {
          all_evidence_list %>%
            purrr::map_dfr(~ safe_pull(.x, d)) %>%
            dplyr::distinct(pmid, .keep_all = TRUE) %>%
            dplyr::arrange(dplyr::desc(score)) %>%
            dplyr::slice_head(n = n_candidates_to_fetch)
        }) %>% stats::setNames(dimension_names)
      } else {
        initial_retrieval <- retrieve_candidates(full_gene_set, enhanced_context, source = "litsense", pctx = relaxed_pctx)
      }
      total_rows <- sum(vapply(initial_retrieval, function(df) if (is.data.frame(df)) nrow(df) else 0L, integer(1)))
      if (total_rows < 1) retrieval_tier <- "empty"
    }
  }
  
  final_evidence <- list()
  allocated_pmids <- c()
  
  for (dim_name in names(initial_retrieval)) {
    unique_top_articles <- initial_retrieval[[dim_name]] %>%
      dplyr::filter(!pmid %in% allocated_pmids) %>%
      dplyr::slice_head(n = 5)
    final_evidence[[dim_name]] <- unique_top_articles
    allocated_pmids <- c(allocated_pmids, unique_top_articles$pmid)
  }
  
  all_articles_pool <- dplyr::bind_rows(initial_retrieval, .id = "dimension") %>%
    dplyr::distinct(pmid, .keep_all = TRUE) %>%
    dplyr::filter(!pmid %in% allocated_pmids) %>%
    dplyr::arrange(dplyr::desc(score))
  
  for (dim_name in names(final_evidence)) {
    needed <- 5 - nrow(final_evidence[[dim_name]])
    if (needed > 0 && nrow(all_articles_pool) > 0) {
      articles_to_add <- all_articles_pool %>% dplyr::slice_head(n = needed)
      if (nrow(articles_to_add) > 0) {
        final_evidence[[dim_name]] <- dplyr::bind_rows(final_evidence[[dim_name]], articles_to_add)
        all_articles_pool <- all_articles_pool %>% dplyr::filter(!pmid %in% articles_to_add$pmid)
      }
    }
    final_evidence[[dim_name]] <- final_evidence[[dim_name]] %>% dplyr::slice_head(n = 5)
  }
  
  all_unique_articles_df <- dplyr::bind_rows(final_evidence) %>%
    dplyr::distinct(pmid, .keep_all = TRUE)
  
  articles_db_subset <- all_unique_articles_df %>%
    dplyr::select(any_of(c("pmid", "text", "score")))
  
  articles_db_columnar <- as.list(articles_db_subset)
  relevance_map <- purrr::map(final_evidence, ~ .x$pmid)
  
  cat("  - Reformatting literature evidence into columnar database structure.\n")
  
  list(
    articles_db = articles_db_columnar,
    relevance_map = relevance_map,
    retrieval_tier = retrieval_tier
  )
}

###################################################################
# SECTION 2: Main generation function (MODIFIED PROMPT)
###################################################################

run_simple_llm_call <- function(prompt_text, model = "deepseek-chat", max_tokens = as.integer(Sys.getenv("DEEPSEEK_MAX_TOKENS", unset = 2000))) {
  api_key <- Sys.getenv("DEEPSEEK_API_KEY")
  base_url <- Sys.getenv("LLM_API_BASE_URL", unset = Sys.getenv("CASSIA_API_BASE_URL", unset = "XXXXX"))
  
  req_body <- list(
    model = model,
    messages = list(list(role = "user", content = prompt_text)),
    temperature = 0,
    top_p = 1,
    presence_penalty = 0,
    frequency_penalty = 0,
    max_tokens = max_tokens
  )
  
  response <- httr::POST(
    url = base_url,
    httr::add_headers(
      `Content-Type` = "application/json",
      `Authorization` = paste("Bearer", api_key)
    ),
    body = jsonlite::toJSON(req_body, auto_unbox = TRUE),
    encode = "raw"
  )
  
  if (httr::status_code(response) == 200) {
    content <- httr::content(response, "parsed")
    return(content$choices[[1]]$message$content)
  } else {
    warning(paste("Simple LLM call for Step 0 failed with status:", httr::status_code(response)))
    return("Pre-analysis failed.")
  }
}

build_expert_report_instructions_core <- function(prompt_context_list, use_context = TRUE, use_neuro_rules = FALSE) {
  context_items <- prompt_context_list
  context_items$user_notes <- NULL
  context_items$dataset_scope <- NULL
  context_items <- unlist(context_items)
  context_items <- context_items[!is.na(context_items)]
  context_description <- if (isTRUE(use_context)) {
    paste(names(context_items), context_items, sep = ": ", collapse = ", ")
  } else {
    "not provided"
  }
  tissue_raw <- if (isTRUE(use_context)) {
    tolower(get_primary_context_value(prompt_context_list$tissue, fallback = ""))
  } else {
    ""
  }
  dataset_scope <- tolower(get_primary_context_value(prompt_context_list$dataset_scope, fallback = "mixed_unknown"))
  scope_profile <- get_primary_context_value(prompt_context_list$scope_profile, fallback = dataset_scope)
  neuro_tissues <- c("brain", "cortex", "dlpfc", "hippocampus", "retina", "eye", "spinal cord", "cerebellum")
  is_neuro_context <- isTRUE(use_neuro_rules) &&
    any(vapply(neuro_tissues, function(x) grepl(x, tissue_raw, fixed = TRUE), logical(1)))
  
  # ==========================================================
  # ECC breakdown schema; the LLM returns integer component scores
  # ==========================================================
  ecc_breakdown_schema <- list(
    marker_support      = "integer (0-2: 0=absent, 1=weak, 2=strong/canonical)",
    functional_support  = "integer (0-2: 0=none, 1=weak/generic, 2=specific match)",
    literature_support  = "integer (0-2: 0=no_snippet, 1=indirect, 2=direct_snippet)",
    candidate_agreement = "integer (0-2: 0=conflict/absent, 1=weak_agree, 2=strong_agree)"
  )
  
  # ==========================================================
  # Structured reasoning fields for explicit contradiction tracking
  # ==========================================================
  structured_reasoning_schema <- list(
    evidence_contradictions = "string (Explicitly list any Top DEGs or Bioinfo results that CONTRADICT this lineage assignment. If none, write 'No hard contradictions'.)",
    overall_assessment = "string",
    lineage_identity = "list[string]",
    functional_capability = "list[string]",
    cellular_state = "list[string]"
  )
  
  context_rule <- if (isTRUE(use_context)) {
    glue::glue("CONTEXT SUMMARY: {context_description}. Scope profile: {scope_profile}. Context is a soft prior used to break ties, except when scope-adaptive gating applies; in that case, core identity must stay within context-consistent broad lineages.")
  } else {
    "TISSUE-BLIND MODE: Do not assume tissue/organ; rely only on DEGs, bioinfo, and evidence."
  }

  base_rules <- c(
    "Return a single valid JSON object in the requested schema; do NOT output validation_status / pass/fail / major/minor. Subtype Level 2 must include core_identity and phenotypic_label.",
    "For Main Type and Subtype Level 1, populate candidate_cell_type with your inferred lineage label. For Subtype Level 2, use core_identity.candidate_cell_type.",
    "EVIDENCE SOURCE RULE: Use in_scope_top_genes as primary evidence for lineage.",
    "LAYERED WORKFLOW RULE: (1) Use Top50 + recovery tier to score main lineage candidates (do NOT hard-lock from a single tier). (2) Then refine subtype only if evidence is sufficient. (3) If subtype evidence is insufficient, stop at the parent lineage (do NOT fall back to generic 'cell').",
    "EVIDENCE COVERAGE RULE: Use all available ranked in-scope genes (Top50 + recovery tier) to score candidates; do NOT rely on Top15 alone.",
    "WEIGHTED CANDIDATE RANKING: Compute a weighted evidence score for each plausible lineage/subtype (Top50=2 points per canonical marker, recovery tier=1 point per canonical marker). Unknown markers may contribute weak support (weight <=0.5) but cannot define lineage alone. List 3-5 candidates ranked by total score, then choose the highest-scoring in-scope candidate.",
    "SUBTYPE THRESHOLD RULE: Do NOT finalize a subtype unless there are >=2 subtype-specific markers supporting it. A single marker may be noted as a hint but cannot define the subtype.",
    "EVIDENCE ROUTING RULE: Top50 = strong evidence, recovery tier = weak evidence, unknown = weak support only (cannot define lineage alone), out-of-scope = contradiction only.",
    "CONFIDENCE CALIBRATION RULE: Set confidence_score_breakdown strictly from evidence. Do NOT force main_type confidence to exceed subtype confidence; subtype may be higher if evidence is stronger.",
    "SCOPE-AWARE TOP15 RULE: In restricted scopes (e.g., immune_enriched), Top15 may only lock a lineage if it contains >=2 canonical scope-consistent markers AND no more than 1 canonical out-of-scope marker. Otherwise Top15 is hint-only and you MUST consult ranks_16_50_in_scope.",
    "WEIGHTED EVIDENCE RULE: Use a soft weighting (Top50 canonical = 2 points, recovery tier canonical = 1 point). Prefer the highest-scoring in-scope lineage. Out-of-scope canonical markers count as contradictions only, not as evidence to switch lineages. Only leave scope if Top50+recovery tier have no coherent in-scope anchors.",
    "CANDIDATE CONSISTENCY RULE: Candidate list is for consistency checking only. If candidates conflict with evidence, ignore them and stay with evidence-derived lineage.",
    "OPTIMALITY PREFERENCE RULE: Among lineage-consistent annotations that would pass Step 2 validation, prefer the most specific subtype supported by evidence. Remaining at a higher-level category is suboptimal when finer resolution is justified.",
    "STAGE REFINEMENT COMMIT RULE: If subtype_level_1 structured_reasoning explicitly identifies a more specific differentiation stage within the same lineage (e.g., 'stage-2 within lineage-A'), you MUST use that refined stage as subtype_level_1.candidate_cell_type. It is NOT allowed to keep a broader parent label at subtype_level_1 when a specific stage is supported by the evidence.",
    "STAGE–LINEAGE DISTINCTION RULE: Labels that primarily describe developmental or cell-cycle stage do not imply a lineage change and should inherit the parent lineage unless contradicted by markers.",
    "MANDATORY IDENTITY REFINEMENT (Program-level): If in_scope_top_genes establish a clear lineage, and secondary evidence (ranks_16_50_in_scope) shows a coherent, identity-defining program, you MUST refine to the most specific compatible identity within that lineage. Do NOT remain at a generic sibling identity when the program-level evidence supports refinement.",
    "PROGRAM EVIDENCE FIELD: Use program_evidence (from enrichment terms) to assess secondary programs for subtype refinement. Treat it as structured evidence, not as primary lineage evidence.",
    "USER-NOTES USAGE RULE: Use user_notes only as a soft prior; do not treat them as lineage evidence.",
    "STEP-2 VALIDATION GATE (MANDATORY): All greedy, organ-specific terminal hypotheses MUST be explicitly validated in Step 2. Validation requires POSITIVE evidence for the defining functional programs of that terminal identity. If validation PASSES, the organ-specific terminal identity may be assigned; if validation FAILS, the model MUST downgrade.",
    "FAIL-SAFE DOWNGRADE RULE: When validation of an organ-specific terminal identity fails, the model MUST (1) keep core_identity at the most specific conservative lineage supported by evidence, (2) preserve the greedy hypothesis ONLY as an association in subtype_level_2.phenotypic_label using non-committal language (e.g., '<organ>-associated', '<terminal>-like'), and (3) explicitly state that the organ-specific hypothesis was considered but not validated.",
    "ASSOCIATIVE ORGAN-AFFINITY RULE (MANDATORY): If evidence suggests an organ/tissue affinity but does not meet the positive program-level validation for a terminal identity, you MUST express the affinity ONLY in subtype_level_2.phenotypic_label using '-like' / 'associated' wording. Do NOT put organ names into core_identity.candidate_cell_type.",
    "HIERARCHY CONSISTENCY HARD RULE: subtype_level_1.candidate_cell_type MUST be a biologically compatible child (or same-level refinement) of main_type.candidate_cell_type. Do NOT output mutually incompatible broad lineages across hierarchy levels. If ambiguity spans multiple broad lineages, choose a higher-level parent that can contain them, then refine conservatively.",
    "IMPORTANT CONSTRAINTS: Organ-specific terminal identities MUST NOT be assigned without functional-program validation. Association terms MUST NOT appear in core_identity.candidate_cell_type. Missing markers alone MUST NOT be used as validation failure; validation depends on positive program-level evidence.",
    "WEAK PRIOR RULE (APPLY LAST): candidates_for_evaluation and Top_Evidence_Genes are weak/noisy hints only. Use them only after forming an independent hypothesis from degs/program_evidence/bioinfo.",
    "MAIN TYPE HARD BOUNDARY: Main Type must prioritize in_scope_top_genes when they provide a coherent, scope-consistent lineage. If Top15 is weak, contaminated, or scope-inconsistent, you MUST consult ranks_16_50_in_scope to choose the best-supported in-scope lineage rather than outputting mixed/ambiguous.",
    "COARSE-LINEAGE FIRST RULE: First decide the broad lineage class using ONLY in_scope_top_genes. Use candidates_for_evaluation only as a weak prior for tie-breaking.",
    "TIE-BREAK RULE (PROGRAM COLLISION): If multiple candidates share similar programs but differ in coarse lineage, choose the coarse lineage best supported by in_scope_top_genes. If still ambiguous, output a higher-level core identity rather than forcing a specific subtype.",
    "LINEAGE-SPECIFIC SECONDARY PROGRAM RULE: When a strong core lineage program is established from top-ranked DEGs, you MUST explicitly evaluate whether secondary evidence (from ranks 16–50 or degs_state_support) supports a more specialized identity within that lineage. If supported, refine the subtype accordingly; if not, remain at the conservative lineage level.",
    "SPECIALIZED PROGRAM RULE: If a specialized lineage program is strong and secondary evidence supports a more specific identity within that lineage, prefer that identity; if ambiguous, stop at a higher-level parent rather than forcing a terminal subtype.",
    "LOW-RANK MARKER LIMIT: markers ranked in the recovery tier may be auxiliary only; they may not override strong Top50 evidence unless Top50 lacks coherent lineage anchors.",
    "OUT-OF-SCOPE MARKER RULE: Use out_of_scope_top_genes ONLY for contamination notes in phenotypic_label/evidence_contradictions. They MUST NOT influence core_identity.",
    "UNKNOWN MARKER RULE: unknown_top_genes have no lineage mapping; treat them as contamination-only and do NOT use them to change core identity.",
    "We are not asking the model to ignore evidence; we are reclassifying evidence into in-scope vs contamination channels. The LLM must remain evidence-based but only in-scope evidence can affect typing.",
    "GLOBAL SUMMARY RULE: Use global_summary only for consistency checking; it must NOT override in_scope_top_genes for core identity.",
    "TOP-15 SCOPE OVERRIDE RULE: When out_of_scope_top_genes conflict with dataset_scope but ranks_16_50_in_scope contain coherent in-scope lineage markers, you MUST prefer the in-scope lineage and treat out_of_scope_top_genes as contamination.",
    "COHERENT ANCHOR SCOPE RULE: In restricted scopes, coherent anchors may come from in_scope_top_genes or ranks_16_50_in_scope; out_of_scope_top_genes are contamination-only and cannot override in-scope anchors.",
    "DOMINANT LINEAGE PRESERVATION RULE: If multiple canonical markers support a coherent lineage, keep that lineage even when a small number of conflicting markers exist; treat conflicts as contamination and record them in evidence_contradictions.",
    "CONFLICT HANDLING RULE: Out-of-scope markers should be recorded as contradictions and reflected in phenotypic_label; they should not overturn core identity unless they form a coherent alternative lineage.",
    "CONTAMINATION THRESHOLD RULE: Override the core lineage only when in-scope lineage markers are weak or absent across in_scope_top_genes and ranks_16_50_in_scope.",
    "IN-SCOPE TIE-BREAK RULE: When multiple in-scope lineages are plausible, prefer the lineage supported by multiple evidence layers (in_scope_top_genes + ranks_16_50_in_scope + program_evidence + literature) over a lineage supported by a single layer.",
    "LAST-RESORT AMBIGUITY RULE: Do NOT output mixed/ambiguous/unknown/other labels or rejection schemas unless instructions explicitly mark FINAL ROUND (round guardrail or validator_feedback). Before that, pick the best-supported lineage and express uncertainty only in phenotypic_label/evidence_contradictions.",
    "SCHEMA-IDENTITY CONSISTENCY RULE: open_world_summary.best_cell_type must align with the broad lineage of core identity.",
       "OPEN-WORLD OUTPUT: fill open_world_summary.best_cell_type, evidence_genes (from in_scope_top_genes), and missing_candidate_suspects (<=3).",
    "2-LAYER DECISION: candidate_cell_type is the core lineage noun; all state/quality terms go to phenotypic_label only.",
    context_rule,
    "CLID RULE: cell_ontology_id must be 'CL:########' or '' (empty). Do not fabricate.",
    "ECC OUTPUT RULE: fill confidence_score_breakdown with integers (marker/function/literature/candidate_agreement, 0-2). Do not compute totals. In structured_reasoning.overall_assessment, the first sentence must be 'ECC breakdown: marker=X/2; function=Y/2; literature=Z/2; candidate_agreement=W/2.'",
    "ITERATION FEEDBACK RULE: If input_data.validator_feedback is provided (list of strings), you MUST revise your output to address these items. Treat them as highest-priority edits. Do NOT change the coarse lineage unless the feedback explicitly says the lineage is wrong. Keep all other mandatory rules unchanged."
  )

  non_neuro_fallback_rules <- c(
    "NON-NEURAL FALLBACK: If lineage anchors are weak, choose the most conservative plausible broad lineage and put dominant state into phenotypic_label; do not force a neural identity. Use 'cell (ambiguous lineage)' only as a FINAL ROUND last resort."
  )

  strategy_rules <- switch(
    dataset_scope,
    immune_enriched = c(
      "DATASET SCOPE: immune_enriched. Prefer conservative in-scope core identities; avoid organ-specific terminal identities.",
      "If evidence is weak or conflicting, output a higher-level in-scope immune lineage rather than forcing a subtype.",
      "Do NOT fall back to 'cell' or non-immune lineages unless Top15+16-50 lack any immune canonical anchors.",
      "Reject only if two incompatible lineages both have strong Top15 anchors; otherwise keep a conservative in-scope core."
    ),
    whole_tissue = c(
      "DATASET SCOPE: whole_tissue. Organ-specific terminal identities are allowed ONLY when program_evidence provides positive support; otherwise keep the parent lineage and express affinity in phenotypic_label.",
      "If evidence is weak, stay at a conservative lineage level rather than forcing an organ-specific subtype."
    ),
    mixed_unknown = c(
    "DATASET SCOPE: mixed_unknown. Prefer the best-supported broad lineage; only output an ambiguous lineage as a FINAL ROUND last resort after scanning in_scope_top_genes and ranks_16_50_in_scope."
    ),
    c("DATASET SCOPE: mixed_unknown. Prefer the best-supported broad lineage; only output an ambiguous lineage as a FINAL ROUND last resort after scanning in_scope_top_genes and ranks_16_50_in_scope.")
  )

  neuro_only_rules <- c(
    "LINEAGE OVERRIDE RULE: If Top DEGs provide strong, consistent lineage markers, prioritize that lineage when dataset_scope is NOT restricted (e.g., mixed_unknown/whole_tissue). For restricted scopes, treat out-of-scope Top DEGs as contamination unless in-scope markers are absent across Top 50.",
    "PROGENITOR/STATE CAUTION: If progenitor or transient-state markers dominate, choose a higher-level lineage label rather than forcing a mature subtype.",
    "ORGANOID DIVERSITY RULE: Organoid datasets may contain off-target lineages beyond the primary tissue. If DEGs strongly support a non-primary lineage, do NOT force a primary-tissue label; follow the DEG evidence."
  )

  extra_rules <- c(
    strategy_rules,
    if (is_neuro_context) {
      neuro_only_rules
    } else if (isTRUE(use_context)) {
      non_neuro_fallback_rules
    } else {
      character(0)
    }
  )

  instructions <- list(
    role = "You are a world-class computational biologist and Knowledge Engineer. Your primary task is to produce a definitive, expert-level hierarchical annotation in a structured JSON format for a knowledge graph.",
    task = "Your task is to analyze the provided dossier and produce a hierarchical JSON. A key innovation in your output is the separation of a cell's stable 'core identity' from its transient 'phenotypic label'.",
    expert_workflow = list(
      description = "To construct the JSON data payload, you MUST follow this mandatory, step-by-step process designed to decouple identity from phenotype.",
      step_0_5 = list(
        title = "Step 0.5: Evidence Consistency Sanity Check",
        instruction = paste(
          "Before any analysis, you MUST perform a critical sanity check on the top-ranked DEGs.",
          "Scope-aware rule: In restricted scopes, out_of_scope_top_genes are contamination-only and must NOT drive lineage decisions.",
          "If Top15 is weak/contaminated, consult ranks_16_50_in_scope to confirm the most plausible in-scope lineage.",
          "Only flag a severe inconsistency if BOTH Top15 and ranks_16_50_in_scope lack coherent in-scope lineage anchors."
        )
      ),
      step_1 = list(
        title = "Step 1: Evidence Deconstruction",
        instruction = paste(
          "List and group the key markers into TWO explicit buckets:",
          "(A) Functional/Pathway markers (what processes they indicate),",
          "(B) Cell-type markers (lineage-specific markers).",
          "Use DEG rank ordering: earlier genes are more important."
        )
      ),
      step_2 = list(
        title = "Step 2: Determine and Document Core Identities",
        instruction = paste(
          "Cross-reference lineage markers with known scRNA-seq knowledge.",
          "Determine the most probable GENERAL cell type.",
          "Then list the TOP 3 plausible subtypes and select the most likely one for `Subtype Level 1`.",
          "Explain why alternatives are less likely."
        )
      ),
      step_3 = list(
        title = "Step 3: Synthesize the Phenotypic Label for the Deepest Subtype",
        instruction = "For `Subtype Level 2`, synthesize a phenotypic label integrating identity + function + state, without altering the core lineage."
      ),
      step_4 = list(
        title = "Step 4: Populate the Decoupled JSON Structure",
        instruction = "Populate the JSON fields and provide a concise summary of the reasoning."
      )
    ),
    final_output_rules = c(
      base_rules,
      extra_rules
    )
  )
  
  instructions$output_format <- list(
    description = "Return answer as a single, valid JSON object.",
    hierarchical_evaluation_schema = list(
      format_description = "Standard hierarchical annotation output.",
      schema = list(
        open_world_summary = list(
          best_cell_type = "string (may be outside candidates)",
          evidence_genes = "list[string] (from Top DEGs)",
          missing_candidate_suspects = "list[string] (<=3; if candidates missed likely lineage)"
        ),
        main_type_schema = list(
          candidate_cell_type = "string",
          cell_ontology_id = "string",
          confidence_score_breakdown = ecc_breakdown_schema,
          structured_reasoning = structured_reasoning_schema
        ),
        subtype_level_1_schema = list(
          candidate_cell_type = "string",
          cell_ontology_id = "string",
          confidence_score_breakdown = ecc_breakdown_schema,
          structured_reasoning = structured_reasoning_schema
        ),
        subtype_level_2_schema = list(
          core_identity = list(
            candidate_cell_type = "string",
            cell_ontology_id = "string"
          ),
          phenotypic_label = "string",
          confidence_score_breakdown = ecc_breakdown_schema,
          structured_reasoning = structured_reasoning_schema
        )
      )
    ),
    quality_control_rejection_schema = list(
      format_description = "ONLY use this format if Step 0.5 fails.",
      schema = list(
        evaluation_status = "string (Must be 'Rejected - Technical Artifact')",
        rejection_reason = "string",
        probable_artifact_type = "string",
        rejection_summary = "string",
        evidence_for_constituent_identities = "list[object]"
      )
    ),
    candidate_rejection_schema = list(
      format_description = "ONLY use this format if data is consistent, but ALL provided candidates are incorrect.",
      schema = list(
        evaluation_status = "string (Must be 'All Candidates Rejected')",
        rejection_reasoning = "string",
        new_hypothesis = "object"
      )
    )
  )
  
  instructions
}

build_expert_report_instructions_tissue_specific <- function(prompt_context_list) {
  tissue_raw <- tolower(get_primary_context_value(prompt_context_list$tissue, fallback = ""))
  neuro_tissues <- c("brain", "cortex", "dlpfc", "hippocampus", "retina", "eye", "spinal cord", "cerebellum")
  is_neuro_context <- any(vapply(neuro_tissues, function(x) grepl(x, tissue_raw, fixed = TRUE), logical(1)))
  build_expert_report_instructions_core(prompt_context_list, use_context = TRUE, use_neuro_rules = is_neuro_context)
}

build_expert_report_instructions_tissue_blind <- function(prompt_context_list) {
  build_expert_report_instructions_core(prompt_context_list, use_context = FALSE, use_neuro_rules = FALSE)
}

build_expert_report_instructions <- function(prompt_context_list) {
  tissue_raw <- tolower(get_primary_context_value(prompt_context_list$tissue, fallback = ""))
  tissue_specific <- nzchar(tissue_raw) && !grepl("unknown|unspecified|mixed|broad|multiple|multi|whole|general|na", tissue_raw)
  if (isTRUE(tissue_specific)) {
    build_expert_report_instructions_tissue_specific(prompt_context_list)
  } else {
    build_expert_report_instructions_tissue_blind(prompt_context_list)
  }
}

###################################################################
# Citation-adder instructions and citation-repair query
###################################################################
build_citation_adder_instructions <- function() {
  list(
    role = "You are an expert citation specialist. Your sole purpose is to add precise and relevant literature citations to a pre-written biological analysis JSON report.",
    task = paste(
      "You will be given a JSON report draft and a library of literature snippets.",
      "Your job is to review EACH sentence within every `structured_reasoning` block.",
      "For EACH sentence, search the entire `literature_evidence` library.",
      "If a claim is directly and specifically supported, append the correct `[pmid:PMID_NUMBER]` tag.",
      "If no direct support is found, leave the sentence unchanged."
    ),
    rules = c(
      "**Rule 0 (Key-Subject Match REQUIRED):** You may ONLY add a citation if the snippet explicitly contains the sentence's key subject. Key subject = (a) the cell type/subtype label, OR (b) at least one KEY gene named in the sentence, OR (c) the exact key process/state phrase claimed. If none of these appear, DO NOT cite.",
    "**Rule 0.5 (No Over-Claim Citations):** Do NOT cite a paper for a sentence that asserts subtype/lineage specificity unless the snippet explicitly contains (i) the same lineage/subtype term used in the sentence AND (ii) the same gene/process phrase claimed.",
      "**Rule 1 (Prioritize Citing Interpretations):** Focus citations on biological interpretations/conclusions. For sentences that merely restate dossier results (GO/KEGG/etc.), citations are optional and should only be used if the snippet adds direct interpretive support.",
      "**Rule 2 (Ensure Specific Support):** A citation is only valid if the snippet specifically supports the core subject of the sentence. Do not over-generalize.",
      "**Rule 3 (Process Sentence-by-Sentence):** For each sentence, scan all provided articles to find the best match. You may add multiple citation tags if multiple articles provide distinct evidence.",
      "**Rule 4 (DO NOT CHANGE CONTENT):** No modification of arguments, scores, or conclusions. Only append `[pmid:PMID_NUMBER]` tags.",
      "**Rule 5 (If No Evidence, Do Nothing):** If no direct and specific support exists, leave the sentence as-is without a citation.",
      "**Rule 6 (OUTPUT FORMAT):** Output must be the complete, valid JSON object from `report_to_be_corrected`, with citations added.",
      "ABSOLUTE FINAL RULE: Only add precise, well-justified `[pmid:...]` tags. Return the entire JSON object."
    )
  )
}

generate_citation_fix_query <- function(step1_query_object) {
  citation_adder_instructions <- build_citation_adder_instructions()
  list(
    query_id = glue::glue("{step1_query_object$query_id}_citation_fix"),
    analysis_type = "Step 1.5: Automated Citation Correction",
    instructions_for_llm = citation_adder_instructions,
    input_data = list(
      literature_evidence = step1_query_object$input_data$cluster_dossier$evidence,
      report_to_be_corrected_source = list(
        description = "The JSON object to be corrected is the direct output from the 'Step 1' query.",
        source_query_id = step1_query_object$query_id
      )
    )
  )
}

###################################################################
# Validator instructions (CASSIA-LIKE, CORE-IDENTITY PASS/FAIL)
###################################################################
  build_validator_instructions <- function() {
    list(
      role = "You are a strict but practical QC validator for LLM-generated cell type annotations.",
      task = paste(
        "You will be given an input_data object with:",
        "1) report_to_be_validated (the generated annotation JSON),",
        "2) original_dossier (DEGs + bioinfo + literature evidence).",
        "Your job is to decide PASS/FAIL using a CASSIA-like philosophy:",
        "PASS if the CORE LINEAGE / CORE IDENTITY is supported and there are no hard contradictions.",
      "DO NOT fail solely because of weak/incorrect citations or over-strong wording in phenotypic details—log those as AUDIT WARNINGS instead.",
      "GREEDY HYPOTHESIS VALIDATION: If the report proposes an organ-specific terminal identity, you MUST check for positive evidence of its defining functional programs. If missing, validation fails and the identity must be downgraded in the report (to a conservative lineage with association-only phenotypic_label).",
        "IMPORTANT: validation_status MUST reflect ONLY core identity (main/sub1/sub2.core_identity)."
      ),
    
    pass_fail_policy = list(
      pass_definition = c(
        "PASS if main_type and subtype_level_1 (and subtype_level_2.core_identity) are consistent with DEG markers AND not contradicted by bioinfo.",
        "PASS even if phenotype/state claims are weakly supported, as long as they are not directly contradictory.",
        "PASS if a dominant lineage is plausible despite minor conflicting anchors; record conflicts as audit warnings rather than FAIL.",
        "validation_status = VALIDATION FAILED ONLY when core identity fails."
      ),
      fail_conditions = c(
        "FAIL if the proposed CORE IDENTITY (main/sub1/sub2.core_identity) contradicts the DEG evidence.",
        "FAIL if the report asserts key markers that are NOT present in the provided DEG list as if they were observed.",
        "FAIL if the report contains direct contradictions against the dossier (e.g., claims pathway enrichment that does not appear in bioinfo lists).",
        "FAIL if there is SCHEMA CONTRADICTION: main_type/subtype_level_1/subtype_level_2.core_identity and open_world_summary.best_cell_type imply different broad lineages.",
        "FAIL if the report assigns an organ-specific terminal identity without positive program-level evidence supporting its defining function (greedy hypothesis validation failed).",
        "FAIL if main_type/subtype_level_1 relies on low-rank markers (>20 in degs) as PRIMARY evidence while in_scope_top_genes support a different lineage/program.",
        "FAIL if main_type/subtype_level_1 is NOT supported by any lineage-consistent markers within in_scope_top_genes or ranks_16_50_in_scope. If in_scope_top_genes are empty or lack coherent anchors, allow ranks_16_50_in_scope to pick the best-supported in-scope lineage and do NOT FAIL unless no plausible lineage exists across the evidence summary.",
        "FAIL if core_identity implies a coarse lineage incompatible with dataset_scope/tissue constraints when the context indicates a restricted lineage scope."
      ),
      
      non_fail_audit_items = c(
        "Over-strong marker language for non-canonical genes -> AUDIT WARNING.",
        "Weak/indirect citations for phenotype/state claims -> AUDIT WARNING.",
        "Wrong or weak citation-to-claim match (PMID exists but does not support the specific claim) -> AUDIT WARNING.",
        "If a PMID appears that is not in the evidence library, add an AUDIT WARNING; do NOT FAIL.",
        "Uncited statements that are clearly traceable to bioinfo fields -> AUDIT WARNING + suggest adding explicit source pointer.",
        "Confidence calibration/hierarchy issues that do not change identity -> AUDIT WARNING (recommend fix).",
        "Greedy organ-specific hypotheses may be preserved only as association language in phenotypic_label when validation fails.",
        "Conflicting anchors from mutually exclusive lineages when a dominant lineage is still plausible -> AUDIT WARNING; require documentation in evidence_contradictions.",
        "Mixed lineage is a FINAL ROUND last resort only: if mixed signals persist after in_scope_top_genes + ranks_16_50_in_scope review, allow a higher-level main_type and record mixed interpretation in subtype_level_2/phenotypic_label. Before final round, require a best-supported lineage choice and record uncertainty only in phenotypic_label.",
        "If out_of_scope_top_genes are present, treat them as contamination and do NOT use them to fail core identity.",
        "In restricted scopes, do not treat out_of_scope_top_genes as coherent lineage evidence if in_scope markers exist in ranks_16_50_in_scope; log as contamination instead.",
        "unknown_top_genes are contamination-only; do NOT use them to override core identity.",
        "If in_scope_top_genes are weak but not empty, prefer a conservative in-scope parent over unknown/other and record uncertainty in phenotypic_label.",
        "Global summary is for consistency checks only; do NOT override in_scope_top_genes with global summary signals."
      )
    ),
    
      rules = c(
        "**0) Rank-Aware Marker Check:** Confirm that the core identity is supported by in_scope_top_genes when possible; if in_scope_top_genes are insufficient, allow ranks_16_50_in_scope to select the best-supported in-scope lineage and avoid mixed/ambiguous unless FINAL ROUND.",
        "**0a) Scope-Aware Top15 Rule:** In restricted scopes, Top15 can only lock a lineage if it has >=2 canonical scope-consistent markers and <=1 canonical out-of-scope marker; otherwise treat Top15 as hint-only and defer to ranks_16_50_in_scope.",
        "**0b) Weighted Evidence Rule:** Use soft weighting (Top15 canonical = 2 points, ranks_16_50 canonical = 1 point). Prefer the highest-scoring in-scope lineage. Out-of-scope canonical markers count as contradictions only, not as evidence to switch lineages.",
        "**0c) Contamination Override:** If out_of_scope_top_genes are present but ranks_16_50_in_scope contain coherent in-scope lineage markers with supporting bioinfo/literature, prefer the in-scope lineage and treat out_of_scope_top_genes as contamination.",
        "**0d) Dominant Lineage Preservation:** If multiple canonical markers support a coherent lineage, keep that lineage unless an alternative lineage has stronger, coherent support across in_scope_top_genes + ranks_16_50_in_scope.",
        "**0.5) Schema Consistency Check:** Ensure all core identity fields and open_world_summary agree on the same broad lineage.",
        "**1) Core Identity Check (Primary):** Verify main_type, subtype_level_1, and subtype_level_2.core_identity against DEGs. Decide whether lineage is correct.",
        "**2) Contradiction Check (Hard):** Identify any direct contradictions against DEGs/bioinfo. Contradictions can trigger FAIL.",
        "**3) Citation Integrity Check:** Do NOT require every claim to be cited. Do NOT FAIL for citation issues. Log mismatches as audit warnings.",
        "**4) Bioinfo Traceability:** If the report mentions KEGG/Reactome/TF results, check they exist in original_dossier.bioinfo. Missing traceability -> audit warning, not fail.",
        "**5) Confidence Logic:** If confidence scores violate hierarchy constraints or look uncalibrated, flag as audit warning unless it fundamentally misrepresents identity.",
        "**6) Output:** Always return JSON in the schema below."
      ),
    
    output_format = list(
      description = "Return a single valid JSON object. PASS/FAIL is for core identity only; audit is separate.",
      schema = list(
        validation_status = "string (Either 'VALIDATION PASSED' or 'VALIDATION FAILED')",
        failure_reasons = "list[string] (Only include if FAILED; concise, evidence-based)",
        core_identity_verdict = list(
          main_type_ok = "boolean",
          subtype_level_1_ok = "boolean",
          subtype_level_2_core_ok = "boolean",
          core_identity_summary = "string (1-3 sentences explaining why identity is or isn't supported)"
        ),
        verdicts = list(
          core_identity = "PASS|FAIL",
          marker_rank_support = "PASS|FAIL",
          schema_consistency = "PASS|FAIL",
          bioinfo_traceability = "PASS|WARN",
          confidence_calibration = "PASS|WARN"
        ),
        audit_warnings = "list[string] (Non-fatal issues: weak citations, overstated marker claims, missing bioinfo pointers, etc.)",
        suggested_fixes = "list[string] (Concrete edits to improve the report without changing core identity)",
        final_verdict_summary = "string (One sentence conclusion; PASS/FAIL reflects core identity only)",
        revision_directive = list(
          description = "REQUIRED when VALIDATION FAILED. Classify the failure type to enable automated triage.",
          schema = list(
            validation_outcome = "string (pass | revise | hard_reject)",
            failure_type = "string (A_conservative_ok | B_overspecific | C_wrong_lineage | null if passed)",
            action = "string (accept_as_is | downgrade_to_parent | stop)",
            suggested_parent = "string | null (CL label of suggested parent if downgrade needed)",
            suggested_parent_clid = "string | null (CL:XXXXXXX format if known)",
            notes = "string (brief explanation of the classification)"
          ),
          classification_rules = c(
            "A_conservative_ok: Prediction is correct direction but not specific enough (i.e., a valid ancestor of the evidence-supported identity). Action: accept_as_is.",
            "B_overspecific: Prediction is correct direction but too specific (i.e., a descendant more specific than evidence supports). Action: downgrade_to_parent.",
            "C_wrong_lineage: Prediction is wrong direction entirely (incompatible broad lineage). Action: stop."
          )
        )
      )
    )
  )
}

infer_coarse_lineage <- function(label, clid = "") {
  lineage_from_name <- function(nm) {
    x <- tolower(as.character(nm %||% ""))
    if (!nzchar(x)) return(NA_character_)
    if (stringr::str_detect(x, "epithel")) return("epithelial")
    if (stringr::str_detect(x, "endothel")) return("endothelial")
    if (stringr::str_detect(x, "fibroblast|stromal|mesenchym")) return("stromal_mesenchymal")
    if (stringr::str_detect(x, "smooth muscle|myocyte|muscle")) return("muscle")
  if (stringr::str_detect(x, "glia|astro|olig|microglia|ependym|neur|neuron|glutamatergic|gaba|interneuron")) return("neural_glial")
    if (stringr::str_detect(x, "immune|lymph|t cell|b cell|nk|myeloid|mono|macro|dendritic")) return("immune")
    if (stringr::str_detect(x, "eryth|platelet|megakaryo")) return("erythroid_megakaryocytic")
    NA_character_
  }

  if (nzchar(clid)) {
    cl <- get_cl_graph()
    if (!is.null(cl) && !is.null(cl[[clid]])) {
      ids <- c(clid)
      anc <- cl[[clid]]$ancestors %||% NULL
      if (!is.null(anc)) {
        ids <- unique(c(ids, names(unlist(anc, use.names = TRUE))))
      }
      for (id in ids) {
        nm <- cl[[id]]$name %||% ""
        lin <- lineage_from_name(nm)
        if (!is.na(lin)) return(lin)
      }
    }
  }

  lin <- lineage_from_name(label)
  if (!is.na(lin)) return(lin)
  "unknown"
}

###################################################################
# Build the LLM input data object with configurable DEG-head size
###################################################################
build_llm_input_data_object <- function(qid,
                                        degs_df,
                                        n_deg_head,
                                        degs_state_df = NULL,
                                        bioinfo_results,
                                        evidence_results,
                                        candidates_df,
                                        prompt_context_list,
                                        smart_context = NULL,
                                        top_k_candidates = 10,
                                        include_preanalysis = FALSE) {
  
  candidates_df <- tibble::as_tibble(candidates_df)
  degs_df <- tibble::as_tibble(degs_df)
  
  degs_df_subset <- degs_df %>%
    utils::head(as.integer(n_deg_head)) %>%
    dplyr::select(any_of(c("geneSymbol", "avg_log2FC", "geneType")))
  
  degs_columnar_data <- as.list(degs_df_subset)


  # Optional: extra DEG evidence for state/quality (not for Main/Sub1 lineage)
  degs_state_columnar_data <- NULL

  
  top_modules <- "NA"
  
  candidates_df_subset <- candidates_df %>%
    dplyr::mutate(
      final_score = round(as.numeric(final_score), 4),
      rank = as.integer(rank)
    )

  if (!"cell_ontology_id" %in% names(candidates_df_subset)) {
    candidates_df_subset$cell_ontology_id <- NA_character_
  }
  if (!"cl_id" %in% names(candidates_df_subset)) {
    candidates_df_subset$cl_id <- NA_character_
  }

  candidates_df_subset <- candidates_df_subset %>%
    dplyr::mutate(
      cell_ontology_id = dplyr::coalesce(
        as.character(.data$cell_ontology_id),
        as.character(.data$cl_id),
        ""
      )
    ) %>%
    dplyr::arrange(rank, dplyr::desc(final_score)) %>%
    utils::head(top_k_candidates) %>%
    dplyr::select(
      Candidate_Cell_Type = cell_type,
      Cell_Ontology_ID = cell_ontology_id,
      Original_Rank = rank,
      Final_Score = final_score,
      Top_Evidence_Genes = top_evidence_genes,
      Rescue_Status = rescue_status
    )

  cl_cfg <- get_cl_cfg()
  if (!is.null(cl_cfg)) {
    norm_clids <- purrr::map_chr(candidates_df_subset$Candidate_Cell_Type, function(lbl) {
      m <- Triage:::normalize_cl_three_state(lbl, "", cl_cfg)
      as.character(m$final_clid %||% "")
    })
    candidates_df_subset$Cell_Ontology_ID <- ifelse(
      nzchar(candidates_df_subset$Cell_Ontology_ID),
      candidates_df_subset$Cell_Ontology_ID,
      norm_clids
    )
  }

  candidates_df_subset <- candidates_df_subset %>%
    dplyr::mutate(
      Candidate_Cell_Type = as.character(Candidate_Cell_Type),
      Cell_Ontology_ID = as.character(Cell_Ontology_ID),
      Coarse_Lineage = purrr::map2_chr(Candidate_Cell_Type, Cell_Ontology_ID, function(lbl, clid) {
        out <- tryCatch(infer_coarse_lineage(lbl, clid), error = function(e) "unknown")
        if (is.null(out) || !nzchar(out)) "unknown" else as.character(out)
      })
    )

  allowed_lineages <- prompt_context_list$allowed_lineages %||% NULL

  normalize_lineage_for_scope <- function(lineage, allowed) {
    allowed <- allowed %||% character(0)
    if ("neural_glial" %in% allowed && lineage %in% c("neural", "glial")) return("neural_glial")
    lineage
  }

  candidates_df_all <- candidates_df_subset

  if (!is.null(allowed_lineages) && "Coarse_Lineage" %in% names(candidates_df_subset)) {
    candidates_df_subset <- candidates_df_subset %>%
      dplyr::mutate(Coarse_Lineage = vapply(.data$Coarse_Lineage, normalize_lineage_for_scope, character(1), allowed = allowed_lineages)) %>%
      dplyr::filter(Coarse_Lineage %in% allowed_lineages)
  }
  
  map_gene_lineages <- function(candidates_df, allowed) {
    gene_map <- list()
    if (is.null(candidates_df) || nrow(candidates_df) == 0) return(gene_map)
    for (i in seq_len(nrow(candidates_df))) {
      genes_raw <- as.character(candidates_df$Top_Evidence_Genes[i] %||% "")
      lineage <- as.character(candidates_df$Coarse_Lineage[i] %||% "unknown")
      lineage <- normalize_lineage_for_scope(lineage, allowed)
      if (!nzchar(genes_raw)) next
      genes <- trimws(unlist(strsplit(genes_raw, "[,;]")))
      genes <- genes[nzchar(genes)]
      for (g in genes) {
        if (is.null(gene_map[[g]])) {
          gene_map[[g]] <- unique(lineage)
        } else {
          gene_map[[g]] <- unique(c(gene_map[[g]], lineage))
        }
      }
    }
    gene_map
  }

  apply_scope_gate <- function(top_genes, gene_map, allowed) {
    top_genes <- as.character(top_genes %||% character(0))
    allowed <- allowed %||% character(0)
    allowed_effective <- setdiff(allowed, "unknown")
    if (length(allowed) == 0) {
      return(list(
        in_scope_top15 = unique(top_genes),
        out_of_scope_top15 = character(0),
        unknown_top15 = character(0),
        out_of_scope_summary = list()
      ))
    }
    in_scope <- character(0)
    out_scope <- character(0)
    unknown_scope <- character(0)
    oos_summary <- list()
    for (g in top_genes) {
      lineages <- gene_map[[g]] %||% character(0)
      lineages <- vapply(lineages, normalize_lineage_for_scope, character(1), allowed = allowed)
      in_allowed <- length(lineages) > 0 && any(lineages %in% allowed_effective)
      if (length(lineages) == 0 || all(lineages == "unknown")) {
        unknown_scope <- c(unknown_scope, g)
      } else if (in_allowed) {
        in_scope <- c(in_scope, g)
      } else {
        out_scope <- c(out_scope, g)
        oos_summary[[g]] <- lineages
      }
    }
    list(
      in_scope_top15 = unique(in_scope),
      out_of_scope_top15 = unique(out_scope),
      unknown_top15 = unique(unknown_scope),
      out_of_scope_summary = oos_summary
    )
  }

  candidates_columnar_data <- as.list(candidates_df_subset)
  
  retrieval_tier <- evidence_results$retrieval_tier %||% "strict"
  evidence_payload <- list(
    articles_db = evidence_results$articles_db %||% list(),
    relevance_map = evidence_results$relevance_map %||% list()
  )

  new_context_description <- glue::glue(
    "Analysis of Cluster {qid}, population: {prompt_context_list$tissue}.

 Notes on inputs (in priority order):
  - PRIMARY: 'degs' contains top positive DEGs (avg_log2FC > 0) and is the primary evidence for lineage assignment (Main/Sub1).
  - SCOPE-AWARE RULE: in_scope_top_genes are the primary tier; if they are weak/contaminated, you MUST consult ranks_16_50_in_scope (recovery tier) to choose the best-supported in-scope lineage. out_of_scope_top_genes are contamination-only and must NOT influence core identity.
- SECONDARY: 'program_evidence' summarizes structured secondary programs for refinement.
 - SUPPORTING: 'bioinfo' (enrichment) and 'evidence' (RAG) should be used for phenotype/state description and for consistency checking.
 - retrieval_tier: {retrieval_tier} (strict|relaxed|empty)
 - WEAK PRIOR (LAST): 'candidates_for_evaluation' is a weak/noisy prior. Use it only AFTER forming an independent hypothesis.
 - WEAK PRIOR (LAST): 'Top_Evidence_Genes' inside candidates are hints ONLY; do NOT treat them as primary evidence.
  - Evidence routing: in_scope_top_genes and ranks_16_50_in_scope are for core identity. out_of_scope_top_genes are contamination-only. unknown_top_genes may be used as weak evidence when in-scope evidence is sparse, but must be weighted lower.
 - We are not asking the model to ignore evidence; we are reclassifying evidence into in-scope vs contamination channels. The LLM must remain evidence-based but only in-scope evidence can affect typing.
"
  )
  
  if (isTRUE(include_preanalysis) && !is.null(smart_context) && nzchar(smart_context)) {
    new_context_description <- paste0(new_context_description, "\n\n[Debug Pre-analysis]\n", smart_context)
  }
  
  program_terms <- build_program_evidence(bioinfo_results)
  program_markers <- build_program_evidence_from_markers(degs_state_df, degs_df)

  gene_map <- map_gene_lineages(candidates_df_all, allowed_lineages)
  top15_genes <- degs_df$geneSymbol %||% character(0)
  if (length(top15_genes) > 15) top15_genes <- top15_genes[1:15]
  select_in_scope_markers <- function(degs_df, gene_map, allowed, limit = 15L, allow_unknown = FALSE) {
    if (is.null(allowed) || length(allowed) == 0) {
      genes <- degs_df$geneSymbol %||% character(0)
      return(utils::head(genes, limit))
    }
    allowed_effective <- setdiff(allowed, "unknown")
    out <- character(0)
    genes <- degs_df$geneSymbol %||% character(0)
    for (g in genes) {
      lineages <- gene_map[[g]] %||% character(0)
      lineages <- vapply(lineages, normalize_lineage_for_scope, character(1), allowed = allowed)
      if (length(lineages) == 0 || all(lineages == "unknown")) {
        if (isTRUE(allow_unknown)) out <- c(out, g)
        next
      }
      if (any(lineages %in% allowed_effective)) out <- c(out, g)
      if (length(out) >= limit) break
    }
    out
  }

  scope_split <- apply_scope_gate(top15_genes, gene_map, allowed_lineages)
  top50_df <- degs_df %>% utils::head(50L) %>% dplyr::select(any_of(c("geneSymbol", "avg_log2FC", "geneType")))
  in_scope_top_limit <- suppressWarnings(as.integer(Sys.getenv("IN_SCOPE_TOP_LIMIT", "50")))
  if (is.na(in_scope_top_limit) || in_scope_top_limit <= 0) in_scope_top_limit <- 50L
  in_scope_recovery_limit <- suppressWarnings(as.integer(Sys.getenv("IN_SCOPE_RECOVERY_LIMIT", "100")))
  if (is.na(in_scope_recovery_limit) || in_scope_recovery_limit <= 0) in_scope_recovery_limit <- 100L
  allow_unknown_in_scope <- Sys.getenv("IN_SCOPE_ALLOW_UNKNOWN", "1") == "1"

  in_scope_markers <- select_in_scope_markers(top50_df, gene_map, allowed_lineages, limit = in_scope_top_limit, allow_unknown = allow_unknown_in_scope)
  in_scope_16_50 <- select_in_scope_markers(top50_df, gene_map, allowed_lineages, limit = in_scope_recovery_limit, allow_unknown = allow_unknown_in_scope)
  in_scope_16_50 <- setdiff(in_scope_16_50, in_scope_markers)

  build_global_summary <- function(degs_df, gene_map, allowed) {
    genes <- degs_df$geneSymbol %||% character(0)
    allowed <- allowed %||% character(0)
    if (length(genes) == 0) return(list())
    lineage_counts <- list()
    oos_count <- 0L
    for (g in genes) {
      lineages <- gene_map[[g]] %||% character(0)
      lineages <- vapply(lineages, normalize_lineage_for_scope, character(1), allowed = allowed)
      if (length(lineages) == 0 || length(allowed) == 0) {
        lineage_counts$unknown <- (lineage_counts$unknown %||% 0L) + 1L
        next
      }
      in_allowed <- lineages[lineages %in% allowed]
      if (length(in_allowed) > 0) {
        for (l in unique(in_allowed)) {
          lineage_counts[[l]] <- (lineage_counts[[l]] %||% 0L) + 1L
        }
      } else {
        oos_count <- oos_count + 1L
      }
    }
    total <- length(genes)
    list(
      total_genes = total,
      in_scope_count = total - oos_count,
      out_of_scope_count = oos_count,
      in_scope_ratio = if (total > 0) round((total - oos_count) / total, 3) else NA_real_,
      lineage_counts = lineage_counts
    )
  }

  top100_df <- degs_df %>% utils::head(100L) %>% dplyr::select(any_of(c("geneSymbol", "avg_log2FC", "geneType")))
  global_summary <- build_global_summary(top100_df, gene_map, allowed_lineages)

  list(
    cluster_dossier = list(
      context = new_context_description,
      degs = list(
        description = "Representative in-scope markers for lineage inference. Use ONLY in_scope_top_genes/ranks_16_50_in_scope for Main/Sub1 lineage; out_of_scope_top_genes and unknown_top_genes are contamination-only.",
        data = list(
          geneSymbol = in_scope_markers
        ),
        in_scope_top_genes = in_scope_markers,
        ranks_16_50_in_scope = in_scope_16_50,
        out_of_scope_top_genes = scope_split$out_of_scope_top15,
        unknown_top_genes = scope_split$unknown_top15,
        out_of_scope_summary = scope_split$out_of_scope_summary
      ),
      scope_profile = prompt_context_list$scope_profile %||% NULL,
      degs_state_support = NULL,
      program_evidence = list(
        description = "Structured secondary program evidence derived from enrichment terms and marker-family patterns.",
        programs_from_terms = program_terms,
        programs_from_markers = program_markers$programs %||% list(),
        support_level = program_markers$support_level %||% list()
      ),
      global_summary = list(
        description = "Global summary derived from top-ranked genes/enrichments for consistency checking only.",
        evidence_lineage_summary = global_summary
      ),
      bioinfo = list(
        description = "Bioinformatics analysis.",
        analysis_results = bioinfo_results
      ),
      evidence = list(
        description = "Literature evidence (RAG).",
        retrieval_tier = retrieval_tier,
        retrieved_articles = evidence_payload
      ),
      candidates_for_evaluation = list(
        description = glue::glue("Weak prior candidates (marker-library-based). Use only for post-hoc consistency checks. Prominent modules: {top_modules}."),
        data = candidates_columnar_data
      )
    )
  )
}

###################################################################
# DEG subset selection helper (unchanged)
###################################################################
select_markers_focus <- function(deg_df, n_pos = 20, n_neg = 20, n_abs = 20) {
  deg_df <- tibble::as_tibble(deg_df)
  
  required_cols <- c("geneSymbol", "avg_log2FC")
  if (!all(required_cols %in% colnames(deg_df))) {
    stop("DEG table must contain columns: geneSymbol, avg_log2FC")
  }
  
  pos_df <- deg_df %>% dplyr::filter(is.finite(avg_log2FC)) %>% dplyr::arrange(dplyr::desc(avg_log2FC)) %>% utils::head(n_pos)
  neg_df <- deg_df %>% dplyr::filter(is.finite(avg_log2FC), avg_log2FC < 0) %>% dplyr::arrange(avg_log2FC) %>% utils::head(n_neg)
  abs_df <- deg_df %>%
    dplyr::filter(is.finite(avg_log2FC)) %>%
    dplyr::mutate(abs_fc = abs(avg_log2FC)) %>%
    dplyr::arrange(dplyr::desc(abs_fc)) %>%
    utils::head(n_abs) %>%
    dplyr::select(-abs_fc)
  
  dplyr::bind_rows(pos_df, neg_df, abs_df) %>%
    dplyr::distinct(geneSymbol, .keep_all = TRUE)
}

###################################################################
# MAIN: generate_expert_report_query
###################################################################
generate_expert_report_query <- function(processed_deg_data,
                                         all_candidates_df,
                                         clusters_to_process,
                                         n_deg_per_cluster,
                                         ppi_files,
                                         run_name,
                                         analysis_config,
                                         dimension_definitions,
                                         prompt_context_config,
                                         top_k_candidates = 10,
                                         include_preanalysis = FALSE) {
  
  cat(glue::glue("\n\n>>> Starting STEP 1 Query Generation for {run_name} (Full Expert Review Mode with decoupleR)\n"))
  `%+%` <- function(a, b) paste0(a, b)
  
  cluster_map <- setNames(clusters_to_process, clusters_to_process)
  cluster_anon_values <- clusters_to_process
  primary_species <- tolower(get_primary_context_value(prompt_context_config$species, "human"))
  future_packages <- c(
    "dplyr", "tibble", "purrr", "stringr", "glue", "readr", "rlang",
    "magrittr", "tictoc", "AnnotationDbi", "clusterProfiler", "ReactomePA",
    "DOSE", "enrichR", "data.table", "tidyr", "decoupleR",
    "org.Hs.eg.db", "httr", "jsonlite", "memoise",
    "cachem", "knitr", "future", "future.apply"
  )
  # disgenet2r is optional: only attach it in workers when it is installed
  # and its API key is configured (the same condition as the main session).
  if (nzchar(Sys.getenv("DISGENET_API_KEY", unset = "")) &&
      requireNamespace("disgenet2r", quietly = TRUE)) {
    future_packages <- c(future_packages, "disgenet2r")
  }
  if (primary_species == "mouse") {
    future_packages <- c(future_packages, "org.Mm.eg.db")
  }
  
  globals_to_export <- list(
    "all_candidates_df" = all_candidates_df,
    "ppi_files" = ppi_files,
    "%+%" = `%+%`,
    "analysis_config" = analysis_config,
    "dimension_definitions" = dimension_definitions,
    "prompt_context_config" = prompt_context_config,
    "disgenet_api_key" = Sys.getenv("DISGENET_API_KEY"),
    "cluster_map" = cluster_map,
    "CL_CFG_CACHE" = CL_CFG_CACHE,
    "CL_GRAPH_CACHE" = CL_GRAPH_CACHE,
    "CL_LOCAL_JSON" = CL_LOCAL_JSON,
    "PROJECT_ROOT" = PROJECT_ROOT,
    "get_cl_graph" = get_cl_graph,
    "cl_normalizer_path" = cl_normalizer_path,
    
    "convert_genes_to_entrez" = convert_genes_to_entrez,
    "summarize_enrich" = summarize_enrich,
    "get_species_resources" = get_species_resources,
    "save_enrichment_tsv" = save_enrichment_tsv,
    "resolve_bioinfo_dir" = resolve_bioinfo_dir,
    # Read at globals-build time (inside the query-generation function,
    # AFTER stage 05 exported TRIAGE_INTERMEDIATE_ROOT); workers do not see
    # environment variables exported after the future plan was created.
    "triage_intermediate_root" = Sys.getenv("TRIAGE_INTERMEDIATE_ROOT", unset = ""),
    "analyze_ppi_from_local_file" = analyze_ppi_from_local_file,
    "analyze_disgenet_enrichment" = analyze_disgenet_enrichment,
    "run_decoupleR_gsea" = run_decoupleR_gsea,
    "enrich_reactome_local" = enrich_reactome_local,
    "setup_dynamic_caches" = setup_dynamic_caches,
    "get_evidence_partitioned" = get_evidence_partitioned,
    "get_primary_context_value" = get_primary_context_value,
    "get_cl_cfg" = get_cl_cfg,
    "extract_terms_recursive" = extract_terms_recursive,
    "build_program_evidence" = build_program_evidence,
    "build_program_evidence_from_markers" = build_program_evidence_from_markers,
    "infer_coarse_lineage" = infer_coarse_lineage,
    "build_llm_input_data_object" = build_llm_input_data_object,
    "build_expert_report_instructions_core" = build_expert_report_instructions_core,
    "build_expert_report_instructions_tissue_specific" = build_expert_report_instructions_tissue_specific,
    "build_expert_report_instructions_tissue_blind" = build_expert_report_instructions_tissue_blind,
    "build_expert_report_instructions" = build_expert_report_instructions,
    
    "build_citation_adder_instructions" = build_citation_adder_instructions,
    "generate_citation_fix_query" = generate_citation_fix_query,
    "build_validator_instructions" = build_validator_instructions,
    "select_markers_focus" = select_markers_focus
  )
  
  llm_input_list <- future.apply::future_lapply(
    X = clusters_to_process,
    FUN = function(qid, deg_data, n_deg_head, run_name_prefix) {
      
      n_deg_head_local <- as.integer(n_deg_head)
      
      tictoc::tic(glue::glue("Total processing for Cluster {qid}"))
      
      step_tag <- "init"
      tryCatch({
        # CL normalization comes from the Triage package namespace; the
        # optional CL_NORMALIZER_PATH override was sourced at load time.
        NULL
        
        if (!is.null(disgenet_api_key) && nzchar(disgenet_api_key)) {
          Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
        }
        
        cat(glue::glue("\n--- Processing Cluster: {qid} ---\n"))
        step_tag <- "mask_cluster"
        masked_qid <- qid
        
        step_tag <- "species_resources"
        current_species <- get_primary_context_value(prompt_context_config$species, "human")
        species_resources <- get_species_resources(current_species)
        
        step_tag <- "setup_caches"
        run_name_for_cache <- file.path(run_name_prefix, masked_qid)
        cached_functions <- setup_dynamic_caches(
          dataset_name = run_name_for_cache,
          config = analysis_config,
          dimension_definitions = dimension_definitions,
          ppi_files = ppi_files,
          species_resources = species_resources
        )
        
        # --- Step 1/4: Data Preparation ---
        tictoc::tic("  -> Step 1/4: Data Preparation, Gene Annotation & Filtering")
        
        annotate_gene_type_from_orgdb <- function(gene_symbols_vector, org_db) {
          if (is.null(gene_symbols_vector) || length(gene_symbols_vector) == 0) {
            return(character(0))
          }
          gene_symbols_vector <- gene_symbols_vector[!is.na(gene_symbols_vector) & nzchar(gene_symbols_vector)]
          if (length(gene_symbols_vector) == 0) {
            return(character(0))
          }
          gene_info <- tryCatch({
            AnnotationDbi::select(
              org_db,
              keys = gene_symbols_vector,
              columns = c("SYMBOL", "GENETYPE"),
              keytype = "SYMBOL"
            )
          }, error = function(e) {
            warning(glue::glue("AnnotationDbi::select query failed for species {current_species}. Returning empty annotation."))
            return(data.frame(SYMBOL = character(0), GENETYPE = character(0)))
          })
          
          gene_map <- gene_info %>%
            dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
            dplyr::rename(geneSymbol = SYMBOL, geneType = GENETYPE)
          
          tibble(geneSymbol = gene_symbols_vector) %>%
            left_join(gene_map, by = "geneSymbol") %>%
            mutate(geneType = ifelse(is.na(geneType), "unknown", geneType)) %>%
            pull(geneType)
        }
        
        step_tag <- "degs_prepare"
        degs_for_cluster_df_full <- deg_data %>%
          dplyr::filter(cluster_name == qid)
        
        # Lineage DEGs (for Main/Sub1): top positive markers only
        degs_for_dossier_lineage_raw <- degs_for_cluster_df_full %>%
          dplyr::filter(is.finite(avg_log2FC), avg_log2FC > 0) %>%
          dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
          utils::head(3L * as.integer(n_deg_head_local))

        # Additional DEG evidence for phenotype/state (NOT for Main/Sub1 lineage)
        degs_for_dossier_state_raw <- select_markers_focus(
          degs_for_cluster_df_full,
          n_pos = n_deg_head_local,
          n_neg = n_deg_head_local,
          n_abs = n_deg_head_local
        )

        degs_for_dossier_lineage <- degs_for_dossier_lineage_raw %>%
          dplyr::mutate(geneType = annotate_gene_type_from_orgdb(geneSymbol, species_resources$org_db))

        degs_for_dossier_state <- degs_for_dossier_state_raw %>%
          dplyr::mutate(geneType = annotate_gene_type_from_orgdb(geneSymbol, species_resources$org_db))
        
        step_tag <- "candidates"
        candidates_for_cluster <- all_candidates_df %>% dplyr::filter(Query_ID == qid)
        
        valid_genes_for_litsense <- degs_for_dossier_lineage %>%
          dplyr::filter(geneType == "protein-coding") %>%
          pull(geneSymbol) %>%
          unique()
        
        cat(glue::glue(
          "    -> LitSense Gene Filtering: Started with {nrow(degs_for_dossier_lineage)} genes, ",
          "{length(valid_genes_for_litsense)} are identified as protein-coding and kept for RAG query.\n"
        ))
        
        tictoc::toc()
        
        # --- Step 2/4: Bioinformatics ---
        tictoc::tic("  -> Step 2/4: Bioinformatics Analysis (with decoupleR)")
        
        step_tag <- "bioinfo"
        marker_df_for_decoupleR <- degs_for_cluster_df_full %>%
          dplyr::select(gene = geneSymbol, avg_log2FC)
        
        stable_marker_df_key <- serialize(marker_df_for_decoupleR, NULL)
        
        bioinfo <- cached_functions$perform_bioinformatics_analysis(
          marker_df_key = stable_marker_df_key,
          config = analysis_config,
          ppi_files = ppi_files,
          prompt_context_config = prompt_context_config
        )
        
        tictoc::toc()
        
        # --- Step 2.2: Save TSV + Summarize ---
        tictoc::tic("  -> Step 2.2/4: Saving TSV files & Summarizing for Dossier")
        
        tsv_output_dir <- resolve_bioinfo_dir(run_name_prefix, masked_qid,
                                             intermediate_root = triage_intermediate_root)
        
        save_enrichment_tsv(bioinfo$go_bp, "go_biological_process", tsv_output_dir)
        save_enrichment_tsv(bioinfo$go_cc, "go_cellular_component", tsv_output_dir)
        save_enrichment_tsv(bioinfo$go_mf, "go_molecular_function", tsv_output_dir)
        save_enrichment_tsv(bioinfo$kegg, "kegg_pathways", tsv_output_dir)
        save_enrichment_tsv(bioinfo$reactome, "reactome_pathways", tsv_output_dir)
        save_enrichment_tsv(bioinfo$disease_do, "disease_ontology", tsv_output_dir)
        save_enrichment_tsv(bioinfo$disease_disgenet, "disease_disgenet", tsv_output_dir)
        
        if (!is.null(bioinfo$decoupleR_gsea)) {
          save_enrichment_tsv(bioinfo$decoupleR_gsea$raw_activities, "decoupleR_gsea_all_tfs", tsv_output_dir)
        }
        
        bioinfo_for_llm <- list(
          biological_processes = summarize_enrich(bioinfo$go_bp, analysis_config$n_bp),
          cellular_components  = summarize_enrich(bioinfo$go_cc, analysis_config$n_cc),
          molecular_functions  = summarize_enrich(bioinfo$go_mf, analysis_config$n_mf),
          kegg_pathways        = summarize_enrich(bioinfo$kegg, analysis_config$n_pathway_each),
          reactome_pathways    = summarize_enrich(bioinfo$reactome, analysis_config$n_pathway_each),
          regulatory_mechanisms = if (!is.null(bioinfo$decoupleR_gsea)) {
            list(activated = bioinfo$decoupleR_gsea$activated, inhibited = bioinfo$decoupleR_gsea$inhibited)
          } else {
            list(activated = character(0), inhibited = character(0))
          },
          protein_interaction_hubs = bioinfo$ppi_hubs,
          disease_associations_do = summarize_enrich(bioinfo$disease_do, analysis_config$n_disease_each),
          disease_associations_disgenet = if (!is.null(bioinfo$disease_disgenet)) {
            head(as.character(bioinfo$disease_disgenet), analysis_config$n_disease_each)
          } else {
            character(0)
          }
        )
        
        tictoc::toc()
        
        # --- Step 3/4: Literature Retrieval (RAG) ---
        tictoc::tic("  -> Step 3/4: Literature Retrieval (RAG)")
        
        step_tag <- "evidence"
        if (length(valid_genes_for_litsense) == 0) {
          cat("    -> WARNING: No valid gene symbols. LitSense evidence will be empty.\n")
          evidence <- list(
            articles_db = list(pmid = character(0), text = character(0), score = numeric(0)),
            relevance_map = list()
          )
        } else {
          stable_litsense_key <- paste(sort(valid_genes_for_litsense), collapse = ",")
          evidence <- get_evidence_partitioned(
            full_gene_set_key = stable_litsense_key,
          context = glue::glue("Cluster {qid} DEGs"),
            config = analysis_config,
            dimension_definitions = dimension_definitions,
            cached_functions = cached_functions,
            prompt_context_config = prompt_context_config
          )
        }
        
        tictoc::toc()
        
        smart_context <- NULL
        if (isTRUE(include_preanalysis)) {
          tictoc::tic("  -> Step 0/5: Generating Smart Context")
          smart_context <- cached_functions$generate_smart_context(
            deg_data_for_prompt = degs_for_dossier_lineage,
            candidates_for_prompt = candidates_for_cluster,
            prompt_context_config = prompt_context_config
          )
          tictoc::toc()
        }
        
        # --- Step 4/4: Final Object Construction ---
        tictoc::tic("  -> Step 4/4: Final Object Construction")
        
        candidates_for_cluster_raw <- all_candidates_df %>% dplyr::filter(Query_ID == qid)
        candidates_for_cluster <- tibble::as_tibble(candidates_for_cluster_raw)
        
        step_tag <- "build_input"
        llm_input_data <- build_llm_input_data_object(
          qid = masked_qid,
          degs_df = degs_for_dossier_lineage,
          n_deg_head = n_deg_head_local,
          degs_state_df = degs_for_dossier_state,
          bioinfo_results = bioinfo_for_llm,
          evidence_results = evidence,
          candidates_df = candidates_for_cluster,
          prompt_context_list = prompt_context_config,
          smart_context = smart_context,
          top_k_candidates = top_k_candidates,
          include_preanalysis = include_preanalysis
        )
        
        expert_instructions <- build_expert_report_instructions(prompt_context_config)
        
        step_tag <- "final_object"
        final_json_object <- list(
          query_id = masked_qid,
          analysis_type = "Step 1: Structured Cell Type Evaluation Report Generation",
          instructions_for_llm = expert_instructions,
          input_data = llm_input_data
        )
        
        tictoc::toc()
        tictoc::toc()
        
        final_json_object
        
      }, error = function(e) {
        cat(glue::glue("\n\n!!! ERROR processing Cluster: {qid} - Skipping. !!!\n"))
        cat(glue::glue("  -> Step: {step_tag}\n"))
        cat(glue::glue("  -> Error: {e$message}\n\n"))
        print(rlang::trace_back())
        NULL
      })
    },
    deg_data = processed_deg_data,
    n_deg_head = n_deg_per_cluster,
    run_name_prefix = run_name,
    future.globals = globals_to_export,
    future.packages = future_packages,
    future.seed = TRUE
  )
  
  final_flattened_list <- llm_input_list[!sapply(llm_input_list, is.null)]
  if (length(final_flattened_list) > 0) {
    names(final_flattened_list) <- sapply(final_flattened_list, `[[`, "query_id")
  }
  final_flattened_list
}

###################################################################
# Post-process confidence in R: aggregate component scores, enforce hierarchy and apply optional caps
###################################################################

compute_ecc_confidence <- function(breakdown) {
  if (is.null(breakdown) || !is.list(breakdown)) return(NA_real_)
  keys <- c("marker_support","functional_support","literature_support","candidate_agreement")
  vals <- vapply(keys, function(k) {
    v <- breakdown[[k]]
    if (is.null(v) || is.na(v)) return(NA_real_)
    suppressWarnings(as.numeric(v))
  }, numeric(1))
  if (any(is.na(vals))) return(NA_real_)
  vals <- pmin(pmax(vals, 0), 2)
  sum(vals) / 8
}

postprocess_report_confidence <- function(report_json) {
  if (is.null(report_json) || !is.list(report_json)) return(report_json)
  
  levels <- c("main_type", "subtype_level_1", "subtype_level_2")
  
  for (lv in levels) {
    if (!is.list(report_json[[lv]])) next
    bd <- report_json[[lv]]$confidence_score_breakdown
    conf <- compute_ecc_confidence(bd)
    report_json[[lv]]$confidence <- conf
  }
  
  # hierarchy constraint removed: allow subtype confidence to exceed main if evidence is stronger
  
  report_json
}

enforce_scope_coarse_gate <- function(report_json, prompt_context_list) {
  if (is.null(report_json) || !is.list(report_json)) return(report_json)

  allowed <- prompt_context_list$allowed_lineages %||% NULL
  if (is.null(allowed) || length(allowed) == 0) return(report_json)

  gate_mode <- tolower(prompt_context_list$gate_mode %||% "hard")

  normalize_conservative_parent <- function(lineage) {
    if (!nzchar(lineage) || lineage == "unknown") return("cell (ambiguous lineage)")
    if (grepl("cell", lineage, ignore.case = TRUE)) return(lineage)
    paste(lineage, "cell")
  }

  normalize_lineage_for_scope <- function(lineage, allowed) {
    allowed <- allowed %||% character(0)
    if ("neural_glial" %in% allowed && lineage %in% c("neural", "glial")) return("neural_glial")
    lineage
  }

  get_lin <- function(lbl, clid = "") {
    out <- tryCatch(infer_coarse_lineage(lbl %||% "", clid %||% ""), error = function(e) "unknown")
    out <- normalize_lineage_for_scope(out, allowed)
    if (!nzchar(out)) "unknown" else out
  }

  main_lin <- get_lin(report_json$main_type$candidate_cell_type, report_json$main_type$cell_ontology_id)
  sub1_lin <- get_lin(report_json$subtype_level_1$candidate_cell_type, report_json$subtype_level_1$cell_ontology_id)
  sub2_lin <- get_lin(report_json$subtype_level_2$core_identity$candidate_cell_type, report_json$subtype_level_2$core_identity$cell_ontology_id)

  bad <- !(main_lin %in% allowed && sub1_lin %in% allowed && sub2_lin %in% allowed)

  if (isTRUE(bad)) {
    if (identical(gate_mode, "flag_only")) {
      report_json$subtype_level_2$phenotypic_label <- paste(
        report_json$subtype_level_2$phenotypic_label %||% "",
        "| NOTE: scope gate flag (out-of-scope core identity detected).",
        sep = " "
      )
      return(report_json)
    }

    pick_in_scope <- function(...) {
      vals <- c(...)
      vals <- vals[!is.na(vals) & nzchar(vals)]
      for (v in vals) {
        lin <- get_lin(v, "")
        if (lin %in% allowed) return(v)
      }
      ""
    }

    candidate_parent <- pick_in_scope(
      report_json$subtype_level_2$core_identity$candidate_cell_type,
      report_json$subtype_level_1$candidate_cell_type,
      report_json$main_type$candidate_cell_type,
      report_json$open_world_summary$best_cell_type
    )
    conservative <- if (nzchar(candidate_parent)) candidate_parent else "cell (ambiguous lineage)"

    main_oos <- !(main_lin %in% allowed)
    sub1_oos <- !(sub1_lin %in% allowed)
    sub2_oos <- !(sub2_lin %in% allowed)

    if (main_oos) {
      report_json$main_type$candidate_cell_type <- conservative
      report_json$main_type$cell_ontology_id <- ""
    }
    if (sub1_oos) {
      report_json$subtype_level_1$candidate_cell_type <- conservative
      report_json$subtype_level_1$cell_ontology_id <- ""
    }
    if (sub2_oos) {
      report_json$subtype_level_2$core_identity$candidate_cell_type <- conservative
      report_json$subtype_level_2$core_identity$cell_ontology_id <- ""
    }

    report_json$subtype_level_2$phenotypic_label <- paste(
      report_json$subtype_level_2$phenotypic_label %||% "",
      "| NOTE: scope-coarse gate triggered; core identity downgraded to a conservative parent.",
      sep = " "
    )
  }

  report_json
}
