# =============================================================================
# detect_format.R - Auto-detect RDS/QS file format and object class
# =============================================================================
# Usage:
#   Rscript detect_format.R <input_file>
#
# Output: JSON object on stdout with keys:
#   file_type    - "rds" or "qs"
#   object_class - Seurat, SingleCellExperiment, dgCMatrix, matrix, list, unknown
#   assays       - [string] list of assay names (for Seurat/SCE objects)
#   slots        - [string] list of available slots/layers
#   dims         - { cells: int, features: int }
#   file_size_mb - float
# =============================================================================

# Source utils from same directory (robust with Rscript)
script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript detect_format.R <input_file>")
}
input_file <- args[1]

# --- Detect file type --------------------------------------------------------
ext <- tolower(tools::file_ext(input_file))
if (!ext %in% c("rds", "qs")) {
  stop(sprintf("Unsupported file extension: .%s (expected .rds or .qs)", ext))
}

# --- Get file size -----------------------------------------------------------
file_size_bytes <- file.info(input_file)$size
file_size_mb <- round(file_size_bytes / 1e6, 2)

# --- Load object (minimal) ---------------------------------------------------
obj <- load_object(input_file)

# --- Inspect object ----------------------------------------------------------
obj_class <- class(obj)[1]
assays    <- character(0)
slots     <- character(0)

if (is_seurat(obj)) {
  obj_type <- "Seurat"
  # Seurat v5 assay names
  assays <- tryCatch(names(obj@assays), error = function(e) character(0))
  # Try to enumerate layers/slots
  if (length(assays) > 0) {
    slots <- tryCatch({
      a1 <- obj[[assays[1]]]
      if (inherits(a1, "Assay5")) {
        names(a1@layers)
      } else {
        c("counts", "data", "scale.data")
      }
    }, error = function(e) c("counts", "data"))
  }
} else if (is_sce(obj)) {
  obj_type <- "SingleCellExperiment"
  # Try to get assay names
  assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
} else if (is_summarized_experiment(obj)) {
  obj_type <- "SummarizedExperiment"
  assays <- tryCatch(SummarizedExperiment::assayNames(obj), error = function(e) character(0))
} else if (inherits(obj, "dgCMatrix") || inherits(obj, "Matrix")) {
  obj_type <- "dgCMatrix"
} else if (is.matrix(obj)) {
  obj_type <- if (is(obj, "matrix")) "matrix" else "matrix"
} else if (is.list(obj)) {
  # Check if it looks like a known structure
  sub_names <- names(obj)
  if (!is.null(sub_names)) {
    obj_type <- sprintf("list[%s]", paste(head(sub_names, 5), collapse = ","))
  } else {
    obj_type <- "list"
  }
} else {
  obj_type <- "unknown"
}

# --- Get dimensions ----------------------------------------------------------
dims <- get_dims(obj)

# --- Warn about large files --------------------------------------------------
if (file_size_mb > 2000) {
  log_warn(sprintf("Large file (%.0f MB). RDS files > 2GB may be unstable.", file_size_mb))
}

# --- Output JSON -------------------------------------------------------------
result <- list(
  file_type    = ext,
  object_class = obj_type,
  assays       = as.list(assays),
  slots        = as.list(slots),
  dims         = dims,
  file_size_mb = file_size_mb
)

cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = FALSE))
