#!/usr/bin/env Rscript
# =============================================================================
# assemble_seurat.R — Build Seurat object from MTX + CSV intermediates
# =============================================================================
#
# Reads the intermediates produced by extract_h5ad.py and constructs a
# Seurat object, then saves as .rds or .qs.
#
# Usage:
#   Rscript assemble_seurat.R <work_dir> <output.rds|output.qs>
#
# Expects these files in <work_dir>:
#   matrix.mtx       — count matrix (features × cells, Market Matrix)
#   barcodes.tsv     — cell barcodes
#   features.tsv     — feature/gene names
#   cell_meta.csv    — obs metadata (optional)
#   var_meta.csv     — var metadata (optional)
#   obsm_pca.csv     — PCA embeddings (optional)
#   obsm_umap.csv    — UMAP embeddings (optional)
#   obsm_tsne.csv    — t-SNE embeddings (optional)
# =============================================================================

script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript assemble_seurat.R <work_dir> <output.rds|output.qs>")
}
work_dir    <- args[1]
output_file <- args[2]
out_ext     <- tolower(tools::file_ext(output_file))

log_info(sprintf("Assembling Seurat from %s", work_dir))

# --- Read MTX ----------------------------------------------------------------
mtx_path <- file.path(work_dir, "matrix.mtx")
if (!file.exists(mtx_path)) {
  stop(sprintf("matrix.mtx not found in %s", work_dir))
}
counts <- Matrix::readMM(mtx_path)  # features × cells
log_info(sprintf("Matrix loaded: %d x %d", nrow(counts), ncol(counts)))

# --- Read barcodes and features ----------------------------------------------
barcodes_path <- file.path(work_dir, "barcodes.tsv")
features_path <- file.path(work_dir, "features.tsv")

if (!file.exists(barcodes_path)) stop("barcodes.tsv not found")
if (!file.exists(features_path)) stop("features.tsv not found")

barcodes <- readLines(barcodes_path)
features <- readLines(features_path)

log_info(sprintf("Cells: %d, Genes: %d", length(barcodes), length(features)))

# --- Set dimnames ------------------------------------------------------------
if (length(features) != nrow(counts)) {
  log_warn(sprintf("features count (%d) != matrix rows (%d)", length(features), nrow(counts)))
  if (length(features) > nrow(counts)) features <- features[1:nrow(counts)]
}
if (length(barcodes) != ncol(counts)) {
  log_warn(sprintf("barcodes count (%d) != matrix cols (%d)", length(barcodes), ncol(counts)))
  if (length(barcodes) > ncol(counts)) barcodes <- barcodes[1:ncol(counts)]
}

rownames(counts) <- features
colnames(counts) <- barcodes

# --- Create Seurat object ----------------------------------------------------
counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")
obj <- CreateSeuratObject(counts = counts, project = "h5ad2r", assay = "RNA")
log_info(sprintf("Seurat object created: %d cells x %d features", ncol(obj), nrow(obj)))

# --- Attach cell metadata (obs) ----------------------------------------------
cell_meta_path <- file.path(work_dir, "cell_meta.csv")
if (file.exists(cell_meta_path)) {
  meta <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)
  common_cells <- intersect(rownames(meta), colnames(obj))
  if (length(common_cells) > 0) {
    meta <- meta[common_cells, , drop = FALSE]
    cols_added <- 0
    for (col in colnames(meta)) {
      if (col %in% colnames(obj@meta.data)) next
      obj@meta.data[[col]] <- NA
      obj@meta.data[common_cells, col] <- meta[[col]]
      cols_added <- cols_added + 1
    }
    log_info(sprintf("Added %d cell metadata columns", cols_added))
  }
}

# --- Attach feature metadata (var) -------------------------------------------
var_meta_path <- file.path(work_dir, "var_meta.csv")
if (file.exists(var_meta_path)) {
  var_meta <- read.csv(var_meta_path, row.names = 1, check.names = FALSE)
  common_genes <- intersect(rownames(var_meta), rownames(obj))
  if (length(common_genes) > 0 && ncol(var_meta) > 0) {
    var_meta <- var_meta[common_genes, , drop = FALSE]
    for (col in colnames(var_meta)) {
      obj[["RNA"]][[col]] <- var_meta[[col]]
    }
    log_info(sprintf("Added %d feature metadata columns", ncol(var_meta)))
  }
}
# --- Attach reductions (obsm) ------------------------------------------------
obsm_specs <- list(
  pca  = list(file = "obsm_pca.csv",  key = "PC_",  name = "pca"),
  umap = list(file = "obsm_umap.csv", key = "UMAP_", name = "umap"),
  tsne = list(file = "obsm_tsne.csv", key = "tSNE_", name = "tsne")
)

for (spec in obsm_specs) {
  csv_path <- file.path(work_dir, spec$file)
  if (!file.exists(csv_path)) next
  df <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  # Align to cell order
  aligned <- mat[colnames(obj), , drop = FALSE]
  if (any(is.na(aligned))) {
    log_warn(sprintf("%d cells missing in %s reduction, filling 0",
                      sum(is.na(aligned[, 1])), spec$name))
    aligned[is.na(aligned)] <- 0
  }
  dreduc <- CreateDimReducObject(
    embeddings = aligned,
    key        = spec$key,
    assay      = "RNA"
  )
  obj[[spec$name]] <- dreduc
  log_info(sprintf("Added reduction '%s': %d dims", spec$name, ncol(aligned)))
}

# --- Save --------------------------------------------------------------------
log_info(sprintf("Saving %s...", basename(output_file)))
if (out_ext == "qs") {
  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("qs package is required for .qs output.")
  }
  qs::qsave(obj, output_file)
} else {
  saveRDS(obj, output_file)
}

output_mb <- round(file.info(output_file)$size / 1e6, 2)
log_info(sprintf("Saved: %s (%.1f MB)", basename(output_file), output_mb))

# --- JSON summary ------------------------------------------------------------
summary <- list(
  status      = "success",
  method      = "mtx",
  output      = output_file,
  output_mb   = output_mb,
  cells       = ncol(obj),
  features    = nrow(obj),
  assays      = as.list(names(obj@assays)),
  reductions  = as.list(names(obj@reductions))
)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = FALSE))
