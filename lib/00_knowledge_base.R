# =============================================================
# 00_knowledge_base.R — knowledge-base loading and scoring functions（reusable with caching）
# 04a ；06/08/09/11 Cell Ontology loading
# =============================================================

# ---- Cell Ontology loading（with caching：RDS cache for faster repeated access）----
load_cl_ontology <- function(json_path, cache_path = NULL) {
  if (!file.exists(json_path)) stop("Ontology file does not exist: ", json_path)
  if (!is.null(cache_path) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  ontology_data <- jsonlite::fromJSON(json_path)

  parents_list_by_name <- list()
  all_names_map <- list()
  for (term_id in names(ontology_data)) {
    child_name <- ontology_data[[term_id]]$label
    if (is.null(child_name)) next
    all_names_map[[term_id]] <- child_name
    ancestors <- ontology_data[[term_id]]$ancestors
    if (!is.null(ancestors) && length(ancestors) > 0) {
      for (ancestor_id in names(ancestors)) {
        if (ancestors[[ancestor_id]] == 1) {
          parent_name <- ontology_data[[ancestor_id]]$label
          if (!is.null(parent_name)) {
            parents_list_by_name[[child_name]] <- parent_name
            break
          }
        }
      }
    }
  }

  res <- list(
    parents = parents_list_by_name,
    name = unlist(all_names_map),
    raw = ontology_data
  )
  if (!is.null(cache_path)) saveRDS(res, cache_path)
  res
}

# ---- Load the CellMarkerAccordion knowledge base（with caching）----
load_accordion_kb <- function(species_filter = "Human", cache_path = NULL) {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  if (!requireNamespace("cellmarkeraccordion", quietly = TRUE)) {
    stop("The cellmarkeraccordion package is required")
  }
  data(accordion_marker, package = "cellmarkeraccordion", envir = environment())
  am_dt <- data.table::as.data.table(accordion_marker)

  cell_to_genes_db <- am_dt %>%
    dplyr::filter(marker_type == "positive", species == species_filter) %>%
    dplyr::group_by(CL_celltype) %>%
    dplyr::summarise(markers = list(unique(marker)), .groups = "drop") %>%
    { stats::setNames(.$markers, .$CL_celltype) }

  cell_to_negative_genes_db <- am_dt %>%
    dplyr::filter(marker_type == "negative", species == species_filter) %>%
    dplyr::group_by(CL_celltype) %>%
    dplyr::summarise(negative_markers = list(unique(marker)), .groups = "drop") %>%
    { stats::setNames(.$negative_markers, .$CL_celltype) }

  ref_sizes <- sapply(cell_to_genes_db, length)
  all_markers_flat <- unlist(cell_to_genes_db, use.names = FALSE)
  gene_counts <- table(all_markers_flat)
  sps_scores <- 1 / gene_counts

  kb <- list(
    am_dt = am_dt,
    cell_to_genes_db = cell_to_genes_db,
    cell_to_negative_genes_db = cell_to_negative_genes_db,
    ref_sizes = ref_sizes,
    sps_scores = sps_scores,
    species = species_filter
  )
  if (!is.null(cache_path)) saveRDS(kb, cache_path)
  kb
}

# ---- V4 scoring: specificity-weighted（primary method）----
annotate_specificity_weighted_V4 <- function(unknown_degs, kb, top_n = 50, penalty_factor = 2, w1 = 0.85, w2 = 0.15) {
  cell_to_genes_db <- kb$cell_to_genes_db
  cell_to_negative_genes_db <- kb$cell_to_negative_genes_db
  ref_sizes <- kb$ref_sizes
  sps_scores <- kb$sps_scores
  unknown_gene_set <- unique(unknown_degs)

  relevant_cells <- kb$am_dt %>%
    dplyr::filter(marker %in% unknown_gene_set, species == kb$species, marker_type == "positive") %>%
    dplyr::pull(CL_celltype) %>%
    unique()

  results_list <- lapply(relevant_cells, function(node_name) {
    node_genes <- cell_to_genes_db[[node_name]]
    if (is.null(node_genes)) return(NULL)
    intersection_genes <- intersect(unknown_gene_set, node_genes)
    if (length(intersection_genes) == 0) return(NULL)
    neg_genes <- cell_to_negative_genes_db[[node_name]]
    if (is.null(neg_genes)) neg_genes <- character(0)
    w_score <- sum(sps_scores[names(sps_scores) %in% intersection_genes], na.rm = TRUE)
    comb_score <- ((w_score / length(unknown_gene_set) + 1e-10)^w1) *
      ((w_score / ref_sizes[[node_name]] + 1e-10)^w2)
    final_ratio <- comb_score -
      (penalty_factor * (length(intersect(unknown_gene_set, neg_genes)) / length(unknown_gene_set)))
    if (final_ratio > 0) {
      list(cell_type = node_name, final_score = final_ratio,
           num_matching_genes = length(intersection_genes),
           evidence_genes = list(intersection_genes))
    } else NULL
  })

  res <- as.data.frame(do.call(rbind, Filter(Negate(is.null), results_list)))
  if (!is.null(res) && nrow(res) > 0) {
    for (i in c(1, 2, 3)) res[[i]] <- unlist(res[[i]])
    res <- res[order(-as.numeric(res$final_score)), ]
    res$rank <- 1:nrow(res)
    return(head(res, top_n))
  }
  data.frame()
}

# ---- V1 scoring: simple-ratio（fallback when V4 fails）----
annotate_simple_ratio_V1 <- function(unknown_degs, kb, top_n = 50, penalty_factor = 2) {
  cell_to_genes_db <- kb$cell_to_genes_db
  cell_to_negative_genes_db <- kb$cell_to_negative_genes_db
  unknown_gene_set <- unique(unknown_degs)

  relevant_cells <- kb$am_dt %>%
    dplyr::filter(marker %in% unknown_gene_set, species == kb$species, marker_type == "positive") %>%
    dplyr::pull(CL_celltype) %>%
    unique()

  results_list <- lapply(relevant_cells, function(node_name) {
    node_genes <- cell_to_genes_db[[node_name]]
    if (is.null(node_genes)) return(NULL)
    intersection_genes <- intersect(unknown_gene_set, node_genes)
    neg_genes <- cell_to_negative_genes_db[[node_name]]
    if (is.null(neg_genes)) neg_genes <- character(0)
    final_ratio <- (length(intersection_genes) / length(unknown_gene_set)) -
      (penalty_factor * (length(intersect(unknown_gene_set, neg_genes)) / length(unknown_gene_set)))
    if (final_ratio > 0) {
      list(cell_type = node_name, final_score = final_ratio,
           num_matching_genes = length(intersection_genes),
           evidence_genes = list(intersection_genes))
    } else NULL
  })

  res <- as.data.frame(do.call(rbind, Filter(Negate(is.null), results_list)))
  if (!is.null(res) && nrow(res) > 0) {
    for (i in c(1, 2, 3)) res[[i]] <- unlist(res[[i]])
    res <- res[order(-as.numeric(res$final_score)), ]
    res$rank <- 1:nrow(res)
    return(head(res, top_n))
  }
  data.frame()
}
