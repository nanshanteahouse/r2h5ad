# =============================================================================
# install_deps.R - One-time dependency installer for r2h5ad
# =============================================================================
# Usage:
#   bash -c "source ~/miniforge3/etc/profile.d/conda.sh && conda activate atacseq-archr && Rscript install_deps.R"
#
# Installs all missing R and Python packages needed for r2h5ad conversion.
# Idempotent: skips packages that are already installed.
#
# Note: R packages are installed via conda (conda-forge) where available,
# with SeuratDisk from GitHub as the exception.
# =============================================================================

# --- Configure CRAN mirror ---------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Helper: install R package if missing (tries CRAN, then conda) -----------
install_if_missing <- function(pkg, conda_pkg = NULL, source = "cran", repo = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[SKIP] %s already installed (v%s)\n",
                pkg, as.character(packageVersion(pkg))))
    return(invisible(FALSE))
  }
  cat(sprintf("[INSTALL] %s ...\n", pkg))

  if (source == "cran") {
    # Try CRAN first
    result <- tryCatch({
      install.packages(pkg, Ncpus = max(1, parallel::detectCores() %/% 2))
      TRUE
    }, error = function(e) {
      cat(sprintf("[WARN]  CRAN install failed: %s\n", conditionMessage(e)))
      FALSE
    })

    # If CRAN fails and conda package name provided, try conda
    if (!result && !is.null(conda_pkg)) {
      cat(sprintf("[INFO]  Trying conda install: %s\n", conda_pkg))
      exit_code <- system2("conda", c("install", "-y", "-c", "conda-forge", conda_pkg))
      if (exit_code != 0) {
        stop(sprintf("Both CRAN and conda install failed for %s", pkg))
      }
    }
  } else if (source == "github") {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", Ncpus = max(1, parallel::detectCores() %/% 2))
    }
    remotes::install_github(repo, upgrade = "never")
  }
  cat(sprintf("[DONE] %s installed\n", pkg))
  return(invisible(TRUE))
}

# --- Helper: install Python package if missing -------------------------------
install_python_if_missing <- function(pkg) {
  installed <- system2("python3", c("-c",
    sprintf("import importlib.util; exit(0 if importlib.util.find_spec('%s') else 1)", pkg)),
    stdout = FALSE, stderr = FALSE)
  exit_code <- attr(installed, "status") %||% installed

  if (identical(exit_code, 0L)) {
    system2("python3", c("-c", sprintf(
      "import %s; print('[SKIP] Python %s already installed (v' + getattr(%s, '__version__', '?') + ')')", pkg, pkg, pkg)))
    return(invisible(FALSE))
  }
  cat(sprintf("[INSTALL] Python package %s ...\n", pkg))
  system2("pip", c("install", pkg))
  cat(sprintf("[DONE] Python %s installed\n", pkg))
  return(invisible(TRUE))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Main --------------------------------------------------------------------
cat("\n========================================\n")
cat("r2h5ad Dependency Installer\n")
cat("========================================\n\n")

cat("Using conda environment. Ensure it is activated before running this script.\n\n")

cat("--- R Packages ---\n")
# Seurat + qs via conda (more reliable for compiled dependencies)
install_if_missing("Seurat", conda_pkg = "r-seurat")
install_if_missing("qs", conda_pkg = "r-qs")
install_if_missing("hdf5r", conda_pkg = "r-hdf5r")  # needed by SeuratDisk
# SeuratDisk from GitHub
install_if_missing("SeuratDisk", source = "github", repo = "mojaveazure/seurat-disk")

cat("\n--- Python Packages ---\n")
install_python_if_missing("anndata")
install_python_if_missing("scanpy")

cat("\n========================================\n")
cat("Dependency installation complete.\n")
cat("========================================\n")

# --- Verification ------------------------------------------------------------
cat("\n--- Verification ---\n")
cat("R packages:\n")
for (pkg in c("Seurat", "SeuratDisk", "qs", "hdf5r")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-12s v%s  [OK]\n", pkg, as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("  %-12s [MISSING]\n", pkg))
  }
}

cat("Python packages:\n")
for (pkg in c("anndata", "scanpy")) {
  code <- sprintf(
    "import importlib.util; spec=importlib.util.find_spec('%s'); print(f'  %-12s [OK]' if spec else '  %-12s [MISSING]')", pkg, pkg, pkg)
  system2("python3", c("-c", shQuote(code)))
}

cat("\nDone.\n")
