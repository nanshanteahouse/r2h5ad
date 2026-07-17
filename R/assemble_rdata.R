#!/usr/bin/env Rscript
# =============================================================================
# assemble_rdata.R — Build .Rdata workspace from MTX + CSV intermediates
# =============================================================================
#
# Reads intermediates from extract_h5ad.py and saves a workspace dump (.Rdata)
# containing named R objects:
#   counts       — dgCMatrix (genes × cells)
#   cell_meta    — data.frame (cell-level metadata)
#   feature_meta — data.frame (feature-level metadata)
#   pca          — matrix (cells × PCs)
#   umap         — matrix (cells × 2)
#   tsne         — matrix (cells × 2)
#
# Usage:
#   Rscript assemble_rdata.R <work_dir> <output.Rdata|output.RData|output.rda>
# =============================================================================

script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript assemble_rdata.R <work_dir> <output.Rdata>")
}
work_dir    <- args[1]
output_file <- args[2]

log_info(sprintf("Assembling .Rdata from %s", work_dir))

# --- Read MTX ----------------------------------------------------------------
mtx_path <- file.path(work_dir, "matrix.mtx")
if (!file.exists(mtx_path)) stop("matrix.mtx not found")
counts <- as(Matrix::readMM(mtx_path), "dgCMatrix")

barcodes_path <- file.path(work_dir, "barcodes.tsv")
features_path <- file.path(work_dir, "features.tsv")

barcodes <- readLines(barcodes_path)
features <- readLines(features_path)

if (length(features) != nrow(counts)) {
  if (length(features) > nrow(counts)) features <- features[1:nrow(counts)]
}
if (length(barcodes) != ncol(counts)) {
  if (length(barcodes) > ncol(counts)) barcodes <- barcodes[1:ncol(counts)]
}
rownames(counts) <- features
colnames(counts) <- barcodes

log_info(sprintf("Counts matrix: %d genes × %d cells", nrow(counts), ncol(counts)))
gc()

# --- Read cell metadata ------------------------------------------------------
cell_meta <- NULL
cell_meta_path <- file.path(work_dir, "cell_meta.csv")
if (file.exists(cell_meta_path)) {
  cell_meta <- read.csv(cell_meta_path, row.names = 1, check.names = FALSE)
  log_info(sprintf("Cell metadata: %d cols × %d cells", ncol(cell_meta), nrow(cell_meta)))
}

# --- Read feature metadata ---------------------------------------------------
feature_meta <- NULL
var_meta_path <- file.path(work_dir, "var_meta.csv")
if (file.exists(var_meta_path)) {
  feature_meta <- read.csv(var_meta_path, row.names = 1, check.names = FALSE)
  log_info(sprintf("Feature metadata: %d cols × %d genes", ncol(feature_meta), nrow(feature_meta)))
}

# --- Read reductions ---------------------------------------------------------
pca <- NULL
umap <- NULL
tsne <- NULL

pca_path <- file.path(work_dir, "obsm_pca.csv")
if (file.exists(pca_path)) {
  pca <- as.matrix(read.csv(pca_path, row.names = 1, check.names = FALSE))
  log_info(sprintf("PCA: %d cells × %d PCs", nrow(pca), ncol(pca)))
}

umap_path <- file.path(work_dir, "obsm_umap.csv")
if (file.exists(umap_path)) {
  umap <- as.matrix(read.csv(umap_path, row.names = 1, check.names = FALSE))
  log_info(sprintf("UMAP: %d cells × %d dims", nrow(umap), ncol(umap)))
}

tsne_path <- file.path(work_dir, "obsm_tsne.csv")
if (file.exists(tsne_path)) {
  tsne <- as.matrix(read.csv(tsne_path, row.names = 1, check.names = FALSE))
  log_info(sprintf("t-SNE: %d cells × %d dims", nrow(tsne), ncol(tsne)))
}

# --- Save .Rdata -------------------------------------------------------------
save_vars <- c("counts")
if (!is.null(cell_meta))    save_vars <- c(save_vars, "cell_meta")
if (!is.null(feature_meta)) save_vars <- c(save_vars, "feature_meta")
if (!is.null(pca))           save_vars <- c(save_vars, "pca")
if (!is.null(umap))          save_vars <- c(save_vars, "umap")
if (!is.null(tsne))          save_vars <- c(save_vars, "tsne")

save(list = save_vars, file = output_file)
output_mb <- round(file.info(output_file)$size / 1e6, 2)

log_info(sprintf("Saved .Rdata: %s (%.1f MB) with objects: %s",
                  basename(output_file), output_mb,
                  paste(save_vars, collapse = ", ")))

# --- JSON summary ------------------------------------------------------------
summary <- list(
  status       = "success",
  method       = "rdata",
  output       = output_file,
  output_mb    = output_mb,
  cells        = ncol(counts),
  features     = nrow(counts),
  objects      = save_vars
)
rm(counts)
gc()
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = FALSE))
