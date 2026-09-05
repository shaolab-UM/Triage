
# =============================================================
# 00_cl_similarity.R — CL ontology similarity utilities (shared library)
# 11_eval_accuracy.R 12_simple_selectors.R
# Use the same scoring scale for all methods: CL similarity (0-100)
# =============================================================

build_cl_similarity <- function(cl_json_path) {
  cl <- jsonlite::fromJSON(cl_json_path, simplifyVector = FALSE)

  # Infer direct parents from ancestors (distance==1)
  cl_ont <- NULL
  if (requireNamespace("ontologyIndex", quietly = TRUE)) {
    ids_all <- names(cl)
    parents_all <- lapply(ids_all, function(id) {
      anc <- cl[[id]]$ancestors %||% NULL
      if (is.null(anc)) return(character(0))
      d <- unlist(anc, use.names = TRUE)
      if (length(d) == 0) return(character(0))
      names(d[d == 1])
    })
    names(parents_all) <- ids_all
    nm_all <- ids_all
    cl_ont <- ontologyIndex::ontology_index(id = ids_all, name = nm_all, parents = parents_all)
  }

  sim_cache <- new.env(parent = emptyenv())

  get_sim_pair_cached <- function(a, b) {
    key <- paste(a, b, sep = "||")
    if (exists(key, envir = sim_cache, inherits = FALSE)) return(get(key, envir = sim_cache))
    if (is.null(cl_ont) || !requireNamespace("ontologySimilarity", quietly = TRUE)) {
      assign(key, NA_real_, envir = sim_cache)
      assign(paste(b, a, sep = "||"), NA_real_, envir = sim_cache)
      return(NA_real_)
    }
    g <- tryCatch(
      ontologySimilarity::get_sim_grid(ontology = cl_ont, term_sets = list(c(a), c(b))),
      error = function(e) NULL
    )
    if (is.null(g)) {
      assign(key, NA_real_, envir = sim_cache)
      assign(paste(b, a, sep = "||"), NA_real_, envir = sim_cache)
      return(NA_real_)
    }
    s <- suppressWarnings(as.numeric(g[1, 2]))
    if (length(s) == 0 || is.na(s)) s <- NA_real_
    if (!is.na(s) && s <= 1) s <- s * 100
    assign(key, s, envir = sim_cache)
    assign(paste(b, a, sep = "||"), s, envir = sim_cache)
    s
  }

  is_valid_clid <- function(x) {
    !is.na(x) && nzchar(x) && startsWith(x, "CL:") && !is.null(cl[[x]])
  }

  score_similarity <- function(pred, truth) {
    if (!is_valid_clid(pred) || !is_valid_clid(truth)) return(0)
    if (pred == truth) return(100)
    s <- get_sim_pair_cached(pred, truth)
    if (!is.na(s)) return(s)
    0
  }

  list(
    cl = cl,
    cl_ont = cl_ont,
    score_similarity = score_similarity,
    is_valid_clid = is_valid_clid
  )
}
