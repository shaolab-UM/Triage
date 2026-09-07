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
#'     if missing). `reactome.db` (Reactome enrichment, stage 05) is part
#'     of the current Bioconductor release.
#'   \item `KEGG.db` is required for canonical KEGG evidence (stage 05
#'     `enrichKEGG(use_internal_data = TRUE)` is backed by KEGG.db) but
#'     was removed from Bioconductor with release 3.11. It is installed
#'     from the verified archived source tarball
#'     `https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz`.
#'   \item CASSIA is installed at the tested revision
#'     `b008c0ac3dd81b2c2dff131d20f5081a58aca027`
#'     (`remotes::install_github("ElliotXie/CASSIA", ref = ..., subdir =
#'     "CASSIA_R")`), not at a moving HEAD. CASSIA additionally requires
#'     a working Python environment for its bundled Python pipeline.
#'     Loading CASSIA at this revision may itself invoke CASSIA's own
#'     `setup_cassia_env()` when its environment is absent (its
#'     `.onLoad()` does this); this helper does not call
#'     `setup_cassia_env()` itself. The workflow preflight verifies
#'     readiness via a version-compatible lookup of CASSIA's internal
#'     `check_python_env()` (present in the tested revision) and checks the
#'     installed CASSIA provenance against the tested revision.
#'   \item `disgenet2r` is optional (only used for the optional
#'     DISGENET disease-evidence step and requires a
#'     `DISGENET_API_KEY`). It is reported, not installed; if a key is
#'     set, install it with
#'     `remotes::install_gitlab("medbio/disgenet2r")`.
#'   \item Seurat is only needed for the optional Seurat-object input
#'     mode of the DEG preparation script; it is reported, not installed.
#' }
#'
#' @param species "human" or "mouse"; selects the organism annotation
#'   package (`org.Hs.eg.db` vs `org.Mm.eg.db`).
#' @param install_cassia logical; install CASSIA at the tested revision.
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
                 "rio", "tidyr", "knitr", "readxl", "enrichR", "xml2", "fs",
                 "R.utils")
  bioc_pkgs <- c("AnnotationDbi",
                 if (species == "mouse") c("org.Mm.eg.db", "org.Hs.eg.db")
                 else "org.Hs.eg.db",
                 "clusterProfiler", "DOSE", "ReactomePA", "decoupleR",
                 "reactome.db")

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

  notes <- character(0)
  # KEGG.db is required for canonical KEGG evidence (stage 05 uses
  # clusterProfiler::enrichKEGG(use_internal_data = TRUE), which is backed
  # by the KEGG.db data package). KEGG.db was removed from Bioconductor
  # with release 3.11 ("use KEGGREST instead"), so the verified install
  # route is the archived Bioconductor source tarball.
  if (!requireNamespace("KEGG.db", quietly = TRUE)) {
    message("Installing KEGG.db from the archived Bioconductor source ",
            "(removed from the current Bioconductor release; URL verified) ...")
    kegg_ok <- tryCatch({
      if (!requireNamespace("remotes", quietly = TRUE)) {
        utils::install.packages("remotes", quiet = TRUE)
      }
      remotes::install_url(
        "https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz",
        quiet = TRUE)
      requireNamespace("KEGG.db", quietly = TRUE)
    }, error = function(e) FALSE)
    if (isTRUE(kegg_ok)) {
      installed <- c(installed, "KEGG.db")
    } else {
      failed <- c(failed, "KEGG.db")
    }
  }
  if (!requireNamespace("disgenet2r", quietly = TRUE)) {
    notes <- c(notes,
               paste0("disgenet2r is OPTIONAL (only for the optional ",
                      "DISGENET disease-evidence step, which requires a ",
                      "DISGENET_API_KEY). Verified install command: ",
                      "remotes::install_gitlab(\"medbio/disgenet2r\"). ",
                      "Without it the workflow skips DisGeNET evidence ",
                      "and every other evidence dimension is unchanged."))
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    notes <- c(notes,
               paste0("Seurat is OPTIONAL (only for the Seurat-object ",
                      "input mode of the DEG preparation script)."))
  }

  if (install_cassia && !requireNamespace("CASSIA", quietly = TRUE)) {
    message("Installing CASSIA at the tested revision ",
            "(ElliotXie/CASSIA @ b008c0ac3dd81b2c2dff131d20f5081a58aca027, ",
            "subdir CASSIA_R) ...")
    cassia_ok <- tryCatch({
      if (!requireNamespace("remotes", quietly = TRUE)) {
        utils::install.packages("remotes", quiet = TRUE)
      }
      remotes::install_github("ElliotXie/CASSIA",
                              ref = "b008c0ac3dd81b2c2dff131d20f5081a58aca027",
                              subdir = "CASSIA_R", quiet = TRUE)
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
    "environment must be available. Loading CASSIA (or calling its ",
    "checkers) may itself trigger CASSIA's own environment setup at its ",
    "tested revision: CASSIA's .onLoad() can invoke setup_cassia_env() ",
    "when its environment is absent. This helper does not call ",
    "setup_cassia_env() itself, but it cannot prevent CASSIA from doing ",
    "so at load time. The workflow preflight verifies readiness with ",
    "a version-compatible lookup of CASSIA's internal check_python_env ",
    "(present in the tested revision)."))
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
#'   \item CellMarkerDB: NOT downloaded automatically. The official source
#'     is the CellMarker 2.0 site
#'     (<http://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/>, Download page
#'     `CellMarker_download.html`); download
#'     `Cell_marker_Human.xlsx` / `Cell_marker_Mouse.xlsx` manually — see
#'     `resources/cellmarkerdb/README.md`.
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
#'   CellMarkerDB has no machine-readable download contract, so it is not
#'   fetched automatically; download the species file from the official
#'   CellMarker 2.0 download page
#'   (`http://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/CellMarker_download.html`)
#'   and supply the path (e.g. `cellmarker_file = file.choose()`); it is
#'   registered under the Triage user-data directory with the canonical
#'   file name. See `resources/cellmarkerdb/README.md`.
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
