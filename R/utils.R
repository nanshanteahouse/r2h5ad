# =============================================================================
# utils.R - Shared utility functions for r2h5ad
# =============================================================================
# Source this file from other R scripts:
#   source(file.path(dirname(sys.frame(1)$ofile), "utils.R"))
#
# Functions:
#   load_rds(path)       - Load an .rds file via readRDS
#   load_qs(path)        - Load a .qs file via qs::qread
#   load_object(path)    - Auto-dispatch based on file extension
#   is_seurat(obj)       - Check if object inherits Seurat class
#   is_sce(obj)          - Check if object inherits SingleCellExperiment
#   win_to_wsl(path)     - Convert Windows path to WSL path
#   wsl_to_win(path)     - Convert WSL path to Windows path
#   log_info(msg)        - Print INFO-level message to stderr
#   log_warn(msg)        - Print WARN-level message to stderr
#   log_error(msg)       - Print ERROR-level message to stderr
#   log_debug(msg)       - Print DEBUG-level message to stderr (when verbose)
# =============================================================================

suppressPackageStartupMessages({
  library(methods)
})

# --- Verbose flag ------------------------------------------------------------
VERBOSE <- Sys.getenv("R2H5AD_VERBOSE", unset = "false") == "true"

# --- Logging ----------------------------------------------------------------
log_info  <- function(msg) cat(sprintf("[INFO]  %s\n", msg), file = stderr())
log_warn  <- function(msg) cat(sprintf("[WARN]  %s\n", msg), file = stderr())
log_error <- function(msg) cat(sprintf("[ERROR] %s\n", msg), file = stderr())
log_debug <- function(msg) {
  if (VERBOSE) cat(sprintf("[DEBUG] %s\n", msg), file = stderr())
}

# --- Object loading ----------------------------------------------------------

#' Load an RDS file
#' @param path Path to .rds file
#' @return Deserialized R object
load_rds <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  log_info(sprintf("Loading RDS: %s", path))
  obj <- tryCatch(
    readRDS(path),
    error = function(e) {
      stop(sprintf("Failed to read RDS file: %s\n  Reason: %s", path, e$message))
    }
  )
  log_info(sprintf("Object size in memory: %.1f MB", object.size(obj) / 1e6))
  return(obj)
}

#' Load a QS file
#' @param path Path to .qs file
#' @return Deserialized R object
load_qs <- function(path) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("Package 'qs' is required. Run: install.packages('qs')")
  }
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  log_info(sprintf("Loading QS: %s", path))
  obj <- tryCatch(
    qs::qread(path),
    error = function(e) {
      stop(sprintf("Failed to read QS file: %s\n  Reason: %s", path, e$message))
    }
  )
  log_info(sprintf("Object size in memory: %.1f MB", object.size(obj) / 1e6))
  return(obj)
}

#' Auto-detect format and load
#' @param path Path to .rds or .qs file
#' @return Deserialized R object
load_object <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(load_rds(path))
  } else if (ext == "qs") {
    return(load_qs(path))
  } else {
    stop(sprintf("Unsupported file extension: .%s (expected .rds or .qs)", ext))
  }
}

# --- Object type detection ---------------------------------------------------

#' Check if object is a Seurat object
is_seurat <- function(obj) {
  inherits(obj, "Seurat")
}

#' Check if object is a SingleCellExperiment
is_sce <- function(obj) {
  inherits(obj, "SingleCellExperiment")
}

#' Check if object is a SummarizedExperiment
is_summarized_experiment <- function(obj) {
  inherits(obj, "SummarizedExperiment")
}

#' Get object dimensions (cells x features)
#' @param obj R object
#' @return Named list with cells and features counts
get_dims <- function(obj) {
  if (is_seurat(obj)) {
    counts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
                       error = function(e) NULL)
    if (is.null(counts) && requireNamespace("Seurat", quietly = TRUE)) {
      # Fallback: try with Seurat v4 API
      counts <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"),
                         error = function(e) NULL)
    }
    if (!is.null(counts)) {
      return(list(cells = ncol(counts), features = nrow(counts)))
    }
    return(list(cells = ncol(obj), features = nrow(obj)))
  } else if (is_sce(obj) || is_summarized_experiment(obj)) {
    return(list(cells = ncol(obj), features = nrow(obj)))
  } else if (is.matrix(obj) || inherits(obj, "Matrix") || inherits(obj, "dgCMatrix")) {
    return(list(cells = ncol(obj), features = nrow(obj)))
  } else if (is.list(obj) && !is.null(obj[["counts"]])) {
    mat <- obj[["counts"]]
    if (is.matrix(mat) || inherits(mat, "Matrix")) {
      return(list(cells = ncol(mat), features = nrow(mat)))
    }
    return(list(cells = 0, features = 0))
  }
  return(list(cells = 0, features = 0))
}

#' Extract count matrix from an object
#' @param obj R object
#' @param assay_name Name of assay to extract (default: "RNA")
#' @return Matrix (dgCMatrix or matrix)
extract_counts <- function(obj, assay_name = "RNA") {
  if (is_seurat(obj)) {
    # Try Seurat v5 layer API first
    counts <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts"),
      error = function(e) NULL
    )
    if (is.null(counts)) {
      # Fallback: Seurat v4 slot API
      counts <- tryCatch(
        Seurat::GetAssayData(obj, assay = assay_name, slot = "counts"),
        error = function(e) NULL
      )
    }
    if (is.null(counts)) {
      # Try other assay names
      assays_avail <- if (inherits(obj, "Seurat")) names(obj@assays) else character(0)
      stop(sprintf(
        "Could not extract counts from assay '%s'. Available assays: %s",
        assay_name, paste(assays_avail, collapse = ", ")
      ))
    }
    return(counts)
  } else if (is_sce(obj) || is_summarized_experiment(obj)) {
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      return(SingleCellExperiment::counts(obj))
    }
    if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      return(SummarizedExperiment::assay(obj, "counts"))
    }
    # Generic: try obj@assays$counts or obj[["counts"]]
    if (is.list(obj@assays) && !is.null(obj@assays$counts)) {
      return(obj@assays$counts)
    }
    stop("Cannot extract counts from SingleCellExperiment/SummarizedExperiment")
  } else if (is.matrix(obj) || inherits(obj, "Matrix") || inherits(obj, "dgCMatrix")) {
    return(obj)
  } else if (is.list(obj) && !is.null(obj[["counts"]])) {
    return(obj[["counts"]])
  } else if (is.list(obj) && !is.null(obj[["mat"]])) {
    return(obj[["mat"]])
  }
  stop("Cannot determine count matrix from object of class: ", paste(class(obj), collapse = ", "))
}

# --- Path conversion ---------------------------------------------------------

#' Convert Windows path to WSL path
#' D:\Projects\file.rds -> /mnt/d/Projects/file.rds
win_to_wsl <- function(path) {
  if (grepl("^[A-Za-z]:", path)) {
    drive <- tolower(substr(path, 1, 1))
    rest  <- substr(path, 3, nchar(path))
    rest  <- gsub("\\\\", "/", rest)
    return(sprintf("/mnt/%s%s", drive, rest))
  }
  return(path)
}

#' Convert WSL path to Windows path
#' /mnt/d/Projects/file.rds -> D:/Projects/file.rds
wsl_to_win <- function(path) {
  if (grepl("^/mnt/[a-z]/", path)) {
    drive <- toupper(substr(path, 6, 6))
    rest  <- substr(path, 7, nchar(path))
    return(sprintf("%s:%s", drive, rest))
  }
  return(path)
}
