#!/usr/bin/env Rscript
#
# convert_rdata.R — Extract objects from .Rdata / .RData files into MTX + CSV
#
# This script is called by rdata2h5ad.sh. It loads an .Rdata file, identifies
# the count matrix and optional metadata objects, and exports them as:
#   - matrix.mtx      (sparse count matrix, Market Matrix format)
#   - barcodes.tsv    (cell barcodes)
#   - features.tsv    (gene names)
#   - cell_meta.csv   (cell-level metadata, e.g. cluster, sample)
#   - obsm_*.csv      (dimensionality reductions: PCA, t-SNE, UMAP)
#
# Usage (via bash wrapper):
#   Rscript convert_rdata.R <input.Rdata> <outdir> [options]
#
# Options (passed as --key=value):
#   count-object=NAME       Name of the count matrix object (auto-detected if omitted)
#   filter-object=NAME      Object whose row names define which cells to keep
#   pca-object=NAME         PCA coordinates object (cells × PCs)
#   tsne-object=NAME        t-SNE coordinates object (cells × 2)
#   umap-object=NAME        UMAP coordinates object (cells × 2)
#   verbose=true            Enable debug logging
# =============================================================================

# ── Parse named arguments ──────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat('{"status":"error","message":"Usage: convert_rdata.R <input.Rdata> <outdir> [--key=value ...]"}')
  quit(status = 1)
}

input_file <- args[1]
outdir     <- args[2]

# Parse --key=value arguments
opts <- list()
for (arg in args[-(1:2)]) {
  if (grepl("^--", arg)) {
    kv <- sub("^--", "", arg)
    parts <- strsplit(kv, "=")[[1]]
    if (length(parts) == 2) {
      opts[[parts[1]]] <- parts[2]
    }
  }
}

verbose <- identical(opts$verbose, "true")
log_info <- function(msg) cat(sprintf("[INFO] %s\n", msg))
log_debug <- function(msg) if (verbose) cat(sprintf("[DEBUG] %s\n", msg))

suppressPackageStartupMessages(library(Matrix))

# ── Load .Rdata ────────────────────────────────────────────────────────────
if (!file.exists(input_file)) {
  cat(sprintf('{"status":"error","message":"File not found: %s"}', input_file))
  quit(status = 1)
}

# .Rdata uses load(), which dumps objects into the current environment
# We capture the environment before and after to identify loaded objects
.before <- ls()
load(input_file)
.after <- ls()
loaded_objs <- setdiff(.after, .before)
# Remove internal names
loaded_objs <- setdiff(loaded_objs, c(".before", ".after", "args", "opts",
                                       "input_file", "outdir", "verbose",
                                       "arg", "kv", "parts", "log_info",
                                       "log_debug", "suppressPackageStartupMessages"))

log_info(sprintf("Loaded %s → %d objects: %s",
                  basename(input_file), length(loaded_objs),
                  paste(loaded_objs, collapse = ", ")))

# ── Collect object info ────────────────────────────────────────────────────
obj_info <- list()
for (nm in loaded_objs) {
  obj <- get(nm)
  cls <- class(obj)[1]
  info <- list(name = nm, class = cls)

  if (is.data.frame(obj) || is.matrix(obj) || inherits(obj, "Matrix")) {
    info$nrow <- nrow(obj)
    info$ncol <- ncol(obj)
    info$size_mb <- as.numeric(round(object.size(obj) / 1e6, 1))
    # Heuristic: count matrix if colnames look like cell barcodes
    info$has_rownames <- !is.null(rownames(obj)) && length(rownames(obj)) > 0
    info$has_colnames <- !is.null(colnames(obj)) && length(colnames(obj)) > 0
    if (info$has_colnames) {
      first_col <- colnames(obj)[1]
      info$is_cell_like <- grepl("^[A-Za-z0-9]+[_-]", first_col) &&
                           nchar(first_col) >= 8
    } else {
      info$is_cell_like <- FALSE
    }
  } else {
    info$nrow <- 0
    info$ncol <- 0
    info$size_mb <- as.numeric(round(object.size(obj) / 1e6, 1))
    info$has_rownames <- !is.null(rownames(obj))
    info$has_colnames <- !is.null(colnames(obj))
    info$is_cell_like <- FALSE
  }
  obj_info[[nm]] <- info
}

# ── Identify count matrix ──────────────────────────────────────────────────
count_obj_name <- opts[["count-object"]]
if (is.null(count_obj_name) || !count_obj_name %in% loaded_objs) {
  # Auto-detect: prefer large numeric matrices with cell-like column names
  candidates <- list()
  for (nm in loaded_objs) {
    info <- obj_info[[nm]]
    obj <- get(nm)
    # Must be numeric data.frame or matrix
    if (!(is.data.frame(obj) || is.matrix(obj) || inherits(obj, "Matrix"))) next
    # Must have substantial size (at least 100 rows AND 100 cols)
    if (info$nrow < 100 || info$ncol < 100) next
    # Heuristic score
    score <- 0
    # Named cell barcodes → strong signal
    if (info$is_cell_like) score <- score + 100
    # Many columns → count matrix typically has the most columns
    score <- score + log10(info$ncol) * 10
    # Many rows → genes
    score <- score + log10(info$nrow) * 5
    # Name contains suggestive keywords
    nm_lower <- tolower(nm)
    if (grepl("count|dge|matrix|expr|umi|data", nm_lower)) score <- score + 50
    # Prefer data.frame (common for old R data) over matrix
    if (is.data.frame(obj)) score <- score + 10
    candidates[[nm]] <- score
  }

  if (length(candidates) == 0) {
    cat('{"status":"error","message":"No suitable count matrix found in .Rdata. Use --count-object to specify."}')
    quit(status = 1)
  }

  count_obj_name <- names(sort(unlist(candidates), decreasing = TRUE))[1]
  log_info(sprintf("Auto-detected count matrix: '%s'", count_obj_name))
} else {
  log_info(sprintf("Using specified count matrix: '%s'", count_obj_name))
}

counts_df <- get(count_obj_name)
log_info(sprintf("Count matrix '%s': %d x %d, %.1f MB",
                  count_obj_name, nrow(counts_df), ncol(counts_df),
                  obj_info[[count_obj_name]]$size_mb))

# ── Identify cell filter object ────────────────────────────────────────────
filter_obj_name <- opts[["filter-object"]]
filter_cells <- NULL
if (!is.null(filter_obj_name) && filter_obj_name %in% loaded_objs) {
  filter_obj <- get(filter_obj_name)
  filter_cells <- rownames(filter_obj)
  n_common <- sum(filter_cells %in% colnames(counts_df))
  log_info(sprintf("Filtering to %d cells from '%s' (%d / %d match colnames)",
                    length(filter_cells), filter_obj_name, n_common, length(filter_cells)))
} else {
  # Auto-detect filter: find objects whose rownames are a subset of count colnames
  for (nm in loaded_objs) {
    if (nm == count_obj_name) next
    obj <- get(nm)
    rn <- rownames(obj)
    if (is.null(rn) || length(rn) == 0) next
    n_match <- sum(rn %in% colnames(counts_df))
    ratio <- n_match / length(rn)
    if (ratio > 0.5 && n_match >= 100) {
      log_info(sprintf("Auto-detected cell filter: '%s' (%d / %d cells match)",
                        nm, n_match, length(rn)))
      filter_cells <- rn
      break
    }
  }
}

# ── Build sparse matrix ────────────────────────────────────────────────────
log_info("Converting to sparse matrix...")

if (!is.null(filter_cells)) {
  # Only keep cells that are in both the filter and the count matrix
  cells_keep <- intersect(filter_cells, colnames(counts_df))
} else {
  cells_keep <- colnames(counts_df)
}
n_cells <- length(cells_keep)
n_genes <- nrow(counts_df)
log_info(sprintf("Keeping %d cells out of %d total", n_cells, ncol(counts_df)))

# Convert to sparse dgCMatrix
# Strategy: build triplet (i, j, x) from non-zero entries per column
pts <- proc.time()
i_list <- vector("list", n_cells)
j_list <- vector("list", n_cells)
x_list <- vector("list", n_cells)
total_nz <- 0

for (col_idx in seq_len(n_cells)) {
  cell_name <- cells_keep[col_idx]
  col_vals <- counts_df[[cell_name]]
  nz_idx <- which(col_vals != 0)
  n_nz <- length(nz_idx)
  total_nz <- total_nz + n_nz
  if (n_nz > 0) {
    i_list[[col_idx]] <- nz_idx
    j_list[[col_idx]] <- rep(col_idx, n_nz)
    x_list[[col_idx]] <- col_vals[nz_idx]
  }
  if (col_idx %% 10000 == 0) {
    elapsed <- (proc.time() - pts)[3]
    log_info(sprintf("  col %d / %d (%.0fs elapsed, %d non-zero)", 
                      col_idx, n_cells, elapsed, total_nz))
  }
}

log_info(sprintf("Building sparse matrix with %d triplets...", total_nz))
i <- unlist(i_list, use.names = FALSE)
j <- unlist(j_list, use.names = FALSE)
x <- unlist(x_list, use.names = FALSE)
rm(i_list, j_list, x_list)
gc()

sparse_counts <- sparseMatrix(
  i = i, j = j, x = x,
  dims = c(n_genes, n_cells),
  dimnames = list(rownames(counts_df), cells_keep),
  giveCsparse = TRUE,
  index1 = TRUE
)
rm(i, j, x)
n_nonzero <- length(sparse_counts@x)
log_info(sprintf("Sparse matrix: %d x %d, %.1f%%%% non-zero (%d entries)",
                  nrow(sparse_counts), ncol(sparse_counts),
                  100 * n_nonzero / prod(dim(sparse_counts)),
                  n_nonzero))

# ── Write MTX + barcodes + features ────────────────────────────────────────
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_info(sprintf("Writing output to %s", outdir))

writeMM(sparse_counts, file.path(outdir, "matrix.mtx"))

# Cell names
writeLines(cells_keep, file.path(outdir, "barcodes.tsv"))

# Gene names
writeLines(rownames(sparse_counts), file.path(outdir, "features.tsv"))

# sparse_counts kept for summary below
gc()

# ── Export metadata objects (PCA, t-SNE, UMAP, cell metadata) ──────────────
meta_dfs <- list()

# PCA object
pca_obj_name <- opts[["pca-object"]]
if (!is.null(pca_obj_name) && pca_obj_name %in% loaded_objs) {
  pca_obj <- get(pca_obj_name)
  if (is.data.frame(pca_obj) || is.matrix(pca_obj)) {
    df <- as.matrix(pca_obj)  # keep as matrix, simpler
    if (!is.null(rownames(pca_obj))) {
      rownames(df) <- rownames(pca_obj)
    }
    colnames(df) <- paste0("PC_", seq_len(ncol(df)))
    write.csv(df, file.path(outdir, "obsm_pca.csv"), row.names = TRUE, quote = FALSE)
    log_info(sprintf("Exported PCA: %d cells x %d PCs", nrow(df), ncol(df)))
  }
}

# Auto-detect PCA object
if (is.null(opts[["pca-object"]])) {
  for (nm in loaded_objs) {
    if (nm == count_obj_name) next
    nm_lower <- tolower(nm)
    if (grepl("pca", nm_lower)) {
      obj <- get(nm)
      if ((is.data.frame(obj) || is.matrix(obj)) && ncol(obj) >= 2 && ncol(obj) <= 200) {
        df <- as.matrix(obj)
        if (!is.null(rownames(obj))) rownames(df) <- rownames(obj)
        colnames(df) <- paste0("PC_", seq_len(ncol(df)))
        write.csv(df, file.path(outdir, "obsm_pca.csv"), row.names = TRUE, quote = FALSE)
        log_info(sprintf("Auto-detected PCA: '%s' (%d cells x %d PCs)", nm, nrow(df), ncol(df)))
        break
      }
      }
    }
  }

# t-SNE object
tsne_obj_name <- opts[["tsne-object"]]
if (!is.null(tsne_obj_name) && tsne_obj_name %in% loaded_objs) {
  tsne_obj <- get(tsne_obj_name)
  if (is.data.frame(tsne_obj) || is.matrix(tsne_obj)) {
    df <- as.matrix(tsne_obj)
    if (!is.null(rownames(tsne_obj))) rownames(df) <- rownames(tsne_obj)
    colnames(df) <- paste0("tSNE_", seq_len(ncol(df)))
    write.csv(df, file.path(outdir, "obsm_tsne.csv"), row.names = TRUE, quote = FALSE)
    log_info(sprintf("Exported t-SNE: %d cells", nrow(df)))
  }
}

# Auto-detect t-SNE object
if (is.null(opts[["tsne-object"]])) {
  for (nm in loaded_objs) {
    if (nm == count_obj_name) next
    nm_lower <- tolower(nm)
    if (grepl("tsne|t.sne", nm_lower)) {
      obj <- get(nm)
      if ((is.data.frame(obj) || is.matrix(obj)) && ncol(obj) %in% c(2, 3)) {
        df <- as.matrix(obj)
        if (!is.null(rownames(obj))) rownames(df) <- rownames(obj)
        colnames(df) <- paste0("tSNE_", seq_len(ncol(df)))
        write.csv(df, file.path(outdir, "obsm_tsne.csv"), row.names = TRUE, quote = FALSE)
        log_info(sprintf("Auto-detected t-SNE: '%s' (%d cells)", nm, nrow(df)))
        break
      }
    }
  }
}

# UMAP object
umap_obj_name <- opts[["umap-object"]]
if (!is.null(umap_obj_name) && umap_obj_name %in% loaded_objs) {
  umap_obj <- get(umap_obj_name)
  if (is.data.frame(umap_obj) || is.matrix(umap_obj)) {
    df <- as.matrix(umap_obj)
    if (!is.null(rownames(umap_obj))) rownames(df) <- rownames(umap_obj)
    colnames(df) <- paste0("UMAP_", seq_len(ncol(df)))
    write.csv(df, file.path(outdir, "obsm_umap.csv"), row.names = TRUE, quote = FALSE)
    log_info(sprintf("Exported UMAP: %d cells", nrow(df)))
  }
}

# Auto-detect UMAP
if (is.null(opts[["umap-object"]])) {
  for (nm in loaded_objs) {
    if (nm == count_obj_name) next
    nm_lower <- tolower(nm)
    if (grepl("umap", nm_lower)) {
      obj <- get(nm)
      if ((is.data.frame(obj) || is.matrix(obj)) && ncol(obj) %in% c(2, 3)) {
        df <- as.matrix(obj)
        if (!is.null(rownames(obj))) rownames(df) <- rownames(obj)
        colnames(df) <- paste0("UMAP_", seq_len(ncol(df)))
        write.csv(df, file.path(outdir, "obsm_umap.csv"), row.names = TRUE, quote = FALSE)
        log_info(sprintf("Auto-detected UMAP: '%s' (%d cells)", nm, nrow(df)))
        break
      }
    }
  }
}

# ── Write object summary JSON (for bash wrapper to parse) ──────────────────
summary_list <- list(
  status = "success",
  count_object = count_obj_name,
  filter_object = if (!is.null(filter_obj_name)) filter_obj_name else "auto",
  n_genes = n_genes,
  n_cells = n_cells,
  total_cells = ncol(counts_df),
  n_nonzero = n_nonzero,
  objects = obj_info
)

rm(sparse_counts, counts_df)
gc()

cat(jsonlite::toJSON(summary_list, auto_unbox = TRUE, pretty = FALSE))
