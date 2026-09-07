#' Install Triage workflow dependencies
#'
#' Installs the R packages needed by the full manuscript workflow scripts
#' (stages 03a-10 under `reproducibility/scripts/pipeline/`). The Triage
#' package itself (including `cellmarkeraccordion`) is installed
#' automatically by `remotes::install_github("shaolab-UM/Triage")` via the
#' package `Remotes` field; this helper covers the remaining
#' full-workflow dependencies and classifies them into CRAN,
#' Bioconductor, and GitHub sources so the user does not have to.
#'
#' Behavior and limitations (stated explicitly):
#' \itemize{
#'   \item CRAN packages are installed with `install.packages()`.
#'   \item Bioconductor packages are installed with
#'     `BiocManager::install()` (`BiocManager` is installed automatically
#'     if missing).
#'   \item CASSIA is installed from its official source
#'     `remotes::install_github("ElliotXie/CASSIA", subdir = "CASSIA_R")`.
#'     CASSIA additionally requires a working Python environment for its
#'     bundled Python pipeline. This function never launches a Python
#'     environment setup automatically; it prints instructions instead.
#'   \item `KEGG.db` is deprecated upstream. An install attempt is made,
#'     and if it fails the user is told how to proceed honestly rather
#'     than silently skipping.
#'   \item `disgenet2r` is optional (only used for the optional
#'     DISGENET disease-evidence step and requires a
#'     `DISGENET_API_KEY`). It is reported, not installed, because no
#'     authoritative automated install source is recorded in this
#'     repository.
#'   \item Seurat is only needed for the optional Seurat-object input
#'     mode of the DEG preparation script; it is reported, not installed.
#' }
#'
#' @param species "human" or "mouse"; selects the organism annotation
#'   package (`org.Hs.eg.db` vs `org.Mm.eg.db`).
#' @param install_cassia logical; install CASSIA from its official
#'   GitHub source.
#'
#' @return Invisibly, a named list: `$installed` (missing packages that
#'   were installed), `$failed` (install attempts that failed), and
#'   `$notes` (reported optional items).
#' @export
install_triage_dependencies <- function(species = "human",
                                        install_cassia = TRUE) {
  species <- match.arg(tolower(species), c("human", "mouse"))

  cran_pkgs <- c("optparse", "digest", "writexl", "glue", "memoise",
                 "cachem", "tictoc", "future", "future.apply", "furrr",
                 "rio", "tidyr", "knitr", "readxl", "enrichR")
  bioc_pkgs <- c("AnnotationDbi",
                 if (species == "mouse") c("org.Mm.eg.db", "org.Hs.eg.db")
                 else "org.Hs.eg.db",
                 "clusterProfiler", "DOSE", "ReactomePA", "decoupleR")

  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]

  failed <- character(0)
  installed <- character(0)

  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    utils::install.packages(missing_cran, quiet = TRUE)
    installed <- c(installed, missing_cran)
    still <- missing_cran[!vapply(missing_cran, requireNamespace,
                                  logical(1), quietly = TRUE)]
    failed <- c(failed, still)
  }
  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager", quiet = TRUE)
    }
    message("Installing Bioconductor packages: ",
            paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, quiet = TRUE, update = FALSE,
                         ask = FALSE)
    installed <- c(installed, missing_bioc)
    still <- missing_bioc[!vapply(missing_bioc, requireNamespace,
                                  logical(1), quietly = TRUE)]
    failed <- c(failed, still)
  }

  # KEGG.db is deprecated upstream; attempt, then report honestly.
  if (!requireNamespace("KEGG.db", quietly = TRUE)) {
    message("Attempting KEGG.db (deprecated upstream; may fail) ...")
    ok <- tryCatch({
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        utils::install.packages("BiocManager", quiet = TRUE)
      }
      BiocManager::install("KEGG.db", quiet = TRUE, update = FALSE,
                           ask = FALSE)
      requireNamespace("KEGG.db", quietly = TRUE)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) {
      installed <- c(installed, "KEGG.db")
    } else {
      failed <- c(failed, "KEGG.db")
      message("KEGG.db could not be installed (deprecated upstream). ",
              "Install a local copy manually if the preflight check ",
              "reports it missing.")
    }
  }

  notes <- character(0)
  if (!requireNamespace("disgenet2r", quietly = TRUE)) {
    notes <- c(notes,
               paste0("disgenet2r is OPTIONAL (only for the optional ",
                      "DISGENET disease-evidence step, which requires a ",
                      "DISGENET_API_KEY). No automated install source is ",
                      "recorded in this repository; install it manually ",
                      "only if you plan to use that step."))
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    notes <- c(notes,
               paste0("Seurat is OPTIONAL (only for the Seurat-object ",
                      "input mode of the DEG preparation script)."))
  }

  if (install_cassia && !requireNamespace("CASSIA", quietly = TRUE)) {
    message("Installing CASSIA from its official source ",
            "(ElliotXie/CASSIA, subdir CASSIA_R) ...")
    cassia_ok <- tryCatch({
      if (!requireNamespace("remotes", quietly = TRUE)) {
        utils::install.packages("remotes", quiet = TRUE)
      }
      remotes::install_github("ElliotXie/CASSIA", subdir = "CASSIA_R",
                              quiet = TRUE)
      requireNamespace("CASSIA", quietly = TRUE)
    }, error = function(e) FALSE)
    if (isTRUE(cassia_ok)) {
      installed <- c(installed, "CASSIA")
    } else {
      failed <- c(failed, "CASSIA")
    }
  }
  message(paste0(
    "CASSIA NOTE: CASSIA runs its annotation pipeline through a bundled ",
    "Python implementation (via reticulate). A working Python ",
    "environment must be available. This function does NOT set up ",
    "Python automatically; follow the Python setup documented in the ",
    "CASSIA repository before running the CASSIA stage."))
  if (!requireNamespace("CASSIA", quietly = TRUE) && !install_cassia) {
    notes <- c(notes, "CASSIA not installed (install_cassia = FALSE); ",
               "stage 03b requires it.")
  }

  invisible(list(installed = installed, failed = failed, notes = notes))
}

#' Download external workflow resources
#'
#' Downloads the external resources needed by the full workflow into a
#' standard per-user data directory (`tools::R_user_dir("Triage",
#' "data")`), and prints exact instructions for wiring them into the
#' workflow repository.
#'
#' Automation is intentionally limited to sources that are verifiably
#' official and machine-readable:
#' \itemize{
#'   \item STRING v12.0: downloaded from the official STRING download
#'     server (`stringdb-downloads.org`, as linked from
#'     https://string-db.org/cgi/download); the `.txt.gz` files are
#'     downloaded and decompressed automatically.
#'   \item CollecTRI: retrieved programmatically with
#'     `decoupleR::get_collectri()` (the official programmatic source
#'     shipped with decoupleR) and saved as an RDS.
#'   \item CellMarkerDB: NOT downloaded automatically. This repository
#'     documents a manual download; see `resources/cellmarkerdb/README.md`.
#' }
#'
#' Retained benchmark checksums (`resources/CHECKSUMS.tsv` and
#' `resources/cellmarkerdb/CHECKSUMS.tsv`) describe the exact snapshot
#' used for the manuscript and are separate from these current
#' functional downloads.
#'
#' Wiring notes printed by this function:
#' \itemize{
#'   \item Installed-package workflow (`run_triage()`): no wiring needed —
#'     resources under the Triage user-data directory are resolved
#'     automatically at run time.
#'   \item Repository workflow (`Rscript
#'     reproducibility/scripts/run_triage.R`): set `TRIAGE_PPI_ROOT` to the
#'     STRING directory and copy the CollecTRI RDS and CellMarkerDB xlsx
#'     into `<repository>/inputs/raw/collectri/` and
#'     `<repository>/inputs/raw/cellmarker/`.
#' }
#'
#' @param species "human" or "mouse".
#' @param download_string logical; download the two STRING v12.0 files.
#' @param download_collectri logical; retrieve CollecTRI via decoupleR.
#' @param cellmarker_file optional path to a locally downloaded CellMarkerDB
#'   spreadsheet (`Cell_marker_Human.xlsx` / `Cell_marker_Mouse.xlsx`).
#'   CellMarkerDB has no machine-readable official download endpoint, so it
#'   cannot be fetched automatically; supply the file you downloaded
#'   manually (e.g. `cellmarker_file = file.choose()`) and it is registered
#'   under the Triage user-data directory with the canonical file name.
#'
#' @return Invisibly, the target directory (a character path).
#' @export
setup_triage_resources <- function(species = "human",
                                   download_string = TRUE,
                                   download_collectri = TRUE,
                                   cellmarker_file = NULL) {
  species <- match.arg(tolower(species), c("human", "mouse"))
  ppi_code <- if (species == "human") "9606" else "10090"
  organism <- if (species == "human") "human" else "mouse"

  base_dir <- tools::R_user_dir("Triage", "data")
  ppi_dir <- file.path(base_dir, "ppi")
  collectri_dir <- file.path(base_dir, "collectri")
  cellmarker_dir <- file.path(base_dir, "cellmarker")
  dir.create(ppi_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(collectri_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cellmarker_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(download_string)) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      stop("setup_triage_resources: the 'R.utils' package is required ",
           "to decompress the STRING archives. Install it with ",
           "install.packages(\"R.utils\").")
    }
    for (fname in c(paste0(ppi_code, ".protein.aliases.v12.0.txt"),
                    paste0(ppi_code, ".protein.physical.links.v12.0.txt"))) {
      subdir <- if (grepl("aliases", fname)) "protein.aliases.v12.0" else
        "protein.physical.links.v12.0"
      url <- sprintf(
        "https://stringdb-downloads.org/download/%s/%s.gz", subdir, fname)
      dest_gz <- file.path(ppi_dir, paste0(fname, ".gz"))
      message("Downloading ", url)
      utils::download.file(url, dest_gz, mode = "wb", quiet = TRUE)
      R.utils::gunzip(dest_gz, remove = TRUE, overwrite = TRUE)
    }
  }

  if (isTRUE(download_collectri)) {
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
      stop("setup_triage_resources: the 'decoupleR' package is required ",
           "to retrieve CollecTRI. Install it with ",
           "install_triage_dependencies() or ",
           "BiocManager::install(\"decoupleR\").")
    }
    message("Retrieving CollecTRI via decoupleR::get_collectri() ...")
    net <- decoupleR::get_collectri(organism = organism,
                                    split_complexes = FALSE)
    saveRDS(net, file.path(collectri_dir,
                           paste0("collectri_", species, "_network.rds")))
  }

  if (!is.null(cellmarker_file) && nzchar(cellmarker_file)) {
    if (!file.exists(cellmarker_file)) {
      stop("setup_triage_resources: cellmarker_file not found: ",
           cellmarker_file)
    }
    canonical <- if (species == "human") "Cell_marker_Human.xlsx" else
      "Cell_marker_Mouse.xlsx"
    dest <- file.path(cellmarker_dir, canonical)
    file.copy(cellmarker_file, dest, overwrite = TRUE)
    message("CellMarkerDB registered as ", dest)
  }

  message("Resources directory: ", base_dir)
  message("WIRING:")
  message("  - Installed-package workflow (run_triage()): no wiring needed; ",
          "resources under ", base_dir, " are resolved automatically.")
  message("  - Repository workflow (Rscript reproducibility/scripts/run_triage.R):")
  message("    1. STRING: set TRIAGE_PPI_ROOT=\"", ppi_dir, "\".")
  message("    2. CollecTRI: copy ",
          file.path(collectri_dir,
                    paste0("collectri_", species, "_network.rds")),
          " into <repository>/inputs/raw/collectri/.")
  message("    3. CellMarkerDB: copy the xlsx into ",
          "<repository>/inputs/raw/cellmarker/.")

  invisible(base_dir)
}
