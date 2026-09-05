# =============================================================
# 00_utils.R — reusable general utility functions
# 01a/02a/02b/03a/03b/04a source
# =============================================================

# NULL/empty/NA-aware coalescing operator (matches runtime semantics of the
# cl_normalizer definition that was active when the release pipeline ran)
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (length(a) == 0) return(b)
  if (length(a) == 1) {
    if (is.na(a)) return(b)
    if (is.character(a) && !nzchar(a)) return(b)
  }
  a
}

# ---- Locate the Triage project root（ source lib ）----
triage_home <- function() {
  home <- Sys.getenv("TRIAGE_HOME", unset = "")
  if (nzchar(home)) return(home)
  # Infer from the script location: scripts/xx.R ->
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) > 0) {
    sp <- sub("^--file=", "", script_arg[[1]])
    return(dirname(dirname(normalizePath(sp, winslash = "/", mustWork = FALSE))))
  }
  getwd()
}

# ---- Noise-gene filtering shared by stages 03 and 04 ----
is_noise_gene <- function(g) {
  g <- toupper(g)
  stringr::str_detect(g, "^MT-") | stringr::str_detect(g, "^RPL") | stringr::str_detect(g, "^RPS") |
    stringr::str_detect(g, "^AC[0-9]") | stringr::str_detect(g, "^LINC") | stringr::str_detect(g, "AS1$")
}

# ---- Fixed gene-exclusion list: ribosomal, mitochondrial and housekeeping genes ----
genes_to_exclude <- c(
  # Ribosomal proteins
  "RPS2","RPS3","RPS3A","RPS4X","RPS4Y1","RPS5","RPS6","RPS7","RPS8","RPS9","RPS10",
  "RPS11","RPS12","RPS13","RPS14","RPS15","RPS15A","RPS16","RPS17","RPS18","RPS19",
  "RPS20","RPS21","RPS23","RPS24","RPS25","RPS26","RPS27","RPS27A","RPS28","RPS29",
  "RPSA","RPS19BP1","RPS27L","RPL3","RPL4","RPL5","RPL6","RPL7","RPL7A","RPL8",
  "RPL9","RPL10","RPL10A","RPL11","RPL12","RPL13","RPL13A","RPL14","RPL15","RPL17",
  "RPL18","RPL18A","RPL19","RPL21","RPL22","RPL23","RPL23A","RPL24","RPL26","RPL27",
  "RPL27A","RPL28","RPL29","RPL30","RPL31","RPL32","RPL34","RPL35","RPL35A","RPL36",
  "RPL36A","RPL37","RPL37A","RPL38","RPL39","RPLP0","RPLP1","RPLP2","RPL22L1",
  "RPL26L1","RPL39L","RPL7L1",
  # Mitochondrial genes
  "MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6","MT-CYB",
  "MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8",
  # Housekeeping genes
  "ACTB","B2M","GAPDH","GUSB","HMBS","HPRT1","PGK1","PPIA","TBP","TFRC","UBC","YWHAZ"
)

# ---- Read a DEG file（detect tab/CSV delimiters automatically）----
read_deg <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  first_line <- readLines(path, n = 1)
  delim <- if (length(first_line) > 0 && grepl("\t", first_line)) "\t" else ","
  df <- readr::read_delim(path, delim = delim, show_col_types = FALSE)
  req_cols <- c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("DEG file is missing columns: ", paste(missing_cols, collapse = ", "),
         " | observed columns: ", paste(names(df), collapse = ", "))
  }
  df
}

# ---- Reference-label normalization（ + ）----
standardize_true_label <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- stringr::str_replace(x, "\\s*[-_ ]\\d+$", "")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

# ---- Filter DEGs： + （03a ）----
filter_deg <- function(deg_df) {
  deg_df %>%
    dplyr::mutate(
      cluster = as.character(cluster),
      gene = as.character(gene),
      geneSymbol_u = toupper(gene)
    ) %>%
    dplyr::filter(
      !geneSymbol_u %in% toupper(genes_to_exclude),
      !is_noise_gene(geneSymbol_u)
    ) %>%
    dplyr::select(-geneSymbol_u)
}

# ---- Subset selection：Optionally retain the first N clusters（for validation）----
subset_clusters <- function(deg_df, n_clusters = NULL) {
  if (is.null(n_clusters) || n_clusters <= 0) return(deg_df)
  clus <- sort(unique(as.character(deg_df$cluster)))
  keep <- head(clus, n_clusters)
  deg_df %>% dplyr::filter(cluster %in% keep)
}
