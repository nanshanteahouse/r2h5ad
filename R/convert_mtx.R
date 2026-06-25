# =============================================================================
# convert_mtx.R - Scheme C: Universal MTX export to h5ad
# =============================================================================
# Usage:
#   Rscript convert_mtx.R <input_file> [output_file.h5ad] [assay_name]
#
# Works with any object that has a count matrix (Seurat, SCE, bare matrix, list).
# Steps:
#   1. Load object
#   2. Extract count matrix
#   3. Write MTX + features.tsv + barcodes.tsv to temp dir
#   4. Call Python to read MTX and save as h5ad (via scanpy)
#   5. Optionally attach cell metadata to the AnnData object
#   6. Clean up temp files
# =============================================================================

# Source utils from same directory
script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript convert_mtx.R <input_file> [output_file.h5ad] [assay_name]")
}

input_file  <- args[1]
output_file <- if (length(args) >= 2) args[2] else {
  file.path(dirname(input_file),
            paste0(tools::file_path_sans_ext(basename(input_file)), ".h5ad"))
}
assay_name  <- if (length(args) >= 3) args[3] else "RNA"

# --- Ensure Matrix package ---------------------------------------------------
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required for MTX export.")
}

# --- Create temp directory ---------------------------------------------------
tmpdir <- file.path(dirname(output_file), paste0(".tmp_mtx_", basename(tools::file_path_sans_ext(output_file))))
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

# Clean up on exit
on.exit({
  unlink(tmpdir, recursive = TRUE)
  log_debug("Cleaned up temp directory")
})

# --- Load object -------------------------------------------------------------
log_info(sprintf("Loading: %s", input_file))
obj <- load_object(input_file)
dims <- get_dims(obj)
log_info(sprintf("Object: %d cells x %d features", dims$cells, dims$features))

# --- Extract count matrix ----------------------------------------------------
log_info(sprintf("Extracting count matrix (assay: %s)...", assay_name))
counts <- extract_counts(obj, assay_name)

# Ensure it's a sparse Matrix
if (!inherits(counts, "Matrix") && !inherits(counts, "dgCMatrix")) {
  log_info("Converting to dgCMatrix...")
  counts <- as(counts, "dgCMatrix")
}

# --- Get feature names and barcodes ------------------------------------------
features <- rownames(counts)
barcodes <- colnames(counts)

if (is.null(features) || length(features) == 0) {
  features <- paste0("gene_", seq_len(nrow(counts)))
  log_warn("No feature names found; using gene_1..gene_N")
}
if (is.null(barcodes) || length(barcodes) == 0) {
  barcodes <- paste0("cell_", seq_len(ncol(counts)))
  log_warn("No barcode names found; using cell_1..cell_N")
}

# --- Write MTX files ---------------------------------------------------------
mtx_path      <- file.path(tmpdir, "matrix.mtx")
features_path <- file.path(tmpdir, "features.tsv")
barcodes_path <- file.path(tmpdir, "barcodes.tsv")

log_debug("Writing matrix.mtx...")
Matrix::writeMM(counts, file = mtx_path)

log_debug("Writing features.tsv...")
write.table(features, file = features_path, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

log_debug("Writing barcodes.tsv...")
write.table(barcodes, file = barcodes_path, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

mtx_size_mb <- round(file.info(mtx_path)$size / 1e6, 2)
log_info(sprintf("MTX written: %.2f MB", mtx_size_mb))

# --- Write cell metadata (if available) --------------------------------------
meta_path <- file.path(tmpdir, "metadata.csv")
has_meta  <- FALSE

if (is_seurat(obj)) {
  meta <- tryCatch(obj[[]], error = function(e) NULL)
  if (!is.null(meta) && ncol(meta) > 0) {
    log_debug(sprintf("Exporting %d metadata columns for %d cells", ncol(meta), nrow(meta)))
    write.csv(meta, file = meta_path, row.names = TRUE)
    has_meta <- TRUE
  }
} else if (is_sce(obj) || is_summarized_experiment(obj)) {
  meta <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(e) NULL)
  if (!is.null(meta) && ncol(meta) > 0) {
    write.csv(meta, file = meta_path, row.names = TRUE)
    has_meta <- TRUE
  }
} else if (is.list(obj) && !is.null(obj[["meta"]])) {
  meta <- obj[["meta"]]
  if (is.data.frame(meta) && ncol(meta) > 0) {
    write.csv(meta, file = meta_path, row.names = TRUE)
    has_meta <- TRUE
  }
}

# --- Call Python to assemble h5ad --------------------------------------------
log_info("Assembling h5ad from MTX...")

# Build Python script directly with path values (avoid sprintf+shQuote double quoting)
py_lines <- c(
  'import sys, os',
  'import scanpy as sc',
  'import pandas as pd',
  '',
  sprintf('mtx_dir = "%s"', tmpdir),
  sprintf('output = "%s"', output_file),
  sprintf('has_meta = %s', if (has_meta) "True" else "False"),
  '',
  '# Read MTX (the input is a directory with matrix.mtx)',
  'adata = sc.read(mtx_dir + "/matrix.mtx", cache=False).T',
  'adata.var_names = pd.read_csv(mtx_dir + "/features.tsv", header=None)[0].values',
  'adata.obs_names = pd.read_csv(mtx_dir + "/barcodes.tsv", header=None)[0].values',
  '',
  '# Attach metadata if available',
  'if has_meta:',
  '    meta_path = os.path.join(mtx_dir, "metadata.csv")',
  '    if os.path.exists(meta_path):',
  '        meta = pd.read_csv(meta_path, index_col=0)',
  '        for col in meta.columns:',
  '            adata.obs[col] = meta[col].reindex(adata.obs_names)',
  '',
  '# Write h5ad',
  'adata.write(output)',
  'print(f"Shape: {adata.shape}")',
  'obs_cols = list(adata.obs.columns) if adata.obs.shape[1] > 0 else "none"',
  'var_cols = list(adata.var.columns) if adata.var.shape[1] > 0 else "none"',
  'print(f"Obs columns: {obs_cols}")',
  'print(f"Var columns: {var_cols}")'
)

py_script_path <- file.path(tmpdir, "assemble_h5ad.py")
writeLines(py_lines, py_script_path)
log_debug(sprintf("Python script written to %s", py_script_path))

# Run the Python script
exit_code <- system2("python3", py_script_path, stdout = TRUE, stderr = TRUE)
status    <- attr(exit_code, "status")

if (is.null(status)) status <- 0L

if (is.character(exit_code)) {
  cat(exit_code, sep = "\n")
}

if (!identical(status, 0L)) {
  stop("Python assembly step failed with exit code ", status)
}

# --- Final validation --------------------------------------------------------
if (!file.exists(output_file)) {
  stop("Output file was not created: ", output_file)
}

output_size_mb <- round(file.info(output_file)$size / 1e6, 2)
log_info(sprintf("Conversion complete: %s (%.1f MB)", output_file, output_size_mb))
log_info(sprintf("Dimensions: %d cells x %d features", dims$cells, dims$features))

cat(sprintf("output=%s\ncells=%d\nfeatures=%d\nsize_mb=%.2f\n",
            output_file, dims$cells, dims$features, output_size_mb))
