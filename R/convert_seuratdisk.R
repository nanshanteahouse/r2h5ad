# =============================================================================
# convert_seuratdisk.R - Scheme A: Seurat object to h5ad via SeuratDisk
# =============================================================================
# Usage:
#   Rscript convert_seuratdisk.R <input_file> [output_file.h5ad]
#
# Converts a Seurat object (RDS/QS) to h5ad using the SeuratDisk pipeline:
#   Seurat object -> .h5Seurat -> .h5ad
#
# This is the recommended path for Seurat objects — preserves all assays,
# metadata, and dimensional reductions.
# =============================================================================

# Source utils from same directory
script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript convert_seuratdisk.R <input_file> [output_file.h5ad]")
}

input_file  <- args[1]
output_file <- if (length(args) >= 2) args[2] else {
  file.path(dirname(input_file),
            paste0(tools::file_path_sans_ext(basename(input_file)), ".h5ad"))
}

# --- Check prerequisites -----------------------------------------------------
for (pkg in c("Seurat", "SeuratDisk")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Run: Rscript install_deps.R", pkg))
  }
}

# --- SeuratDisk compatibility patch: translate slot= -> layer= ---------------
# SeuratDisk is unmaintained and uses the slot= argument which was made defunct
# in SeuratObject 5.3.0+. We monkey-patch GetAssayData/SetAssayData to translate
# slot= -> layer= so SaveH5Seurat works with modern SeuratObject.
.r2h5ad_patch_sd <- local({
  if (packageVersion("SeuratObject") >= "5.3.0") {
    log_info("Patching GetAssayData/SetAssayData slot->layer for SeuratDisk compat...")
    .orig_get <- SeuratObject::GetAssayData
    .orig_set <- SeuratObject::SetAssayData

    patched_get <- function(object, layer = "data", slot, ...) {
      if (!missing(slot)) layer <- slot
      .orig_get(object = object, layer = layer, ...)
    }
    patched_set <- function(object, layer, slot, ...) {
      if (!missing(slot)) layer <- slot
      .orig_set(object = object, layer = layer, ...)
    }

    # Patch in SeuratObject namespace (for direct callers)
    assignInNamespace("GetAssayData", patched_get, ns = "SeuratObject")
    assignInNamespace("SetAssayData", patched_set, ns = "SeuratObject")

    # Patch in SeuratDisk import environment (for SeuratDisk internal calls)
    # SeuratDisk imports these from SeuratObject but byte-compilation can cache
    # the original references. Unlocking the import bindings and replacing them
    # ensures SeuratDisks SaveH5Seurat path uses our patched versions.
    sd_imports <- parent.env(asNamespace("SeuratDisk"))
    for (fname in c("GetAssayData", "SetAssayData")) {
      if (exists(fname, envir = sd_imports, inherits = FALSE)) {
        patched <- if (fname == "GetAssayData") patched_get else patched_set
        if (bindingIsLocked(fname, sd_imports)) unlockBinding(fname, sd_imports)
        assign(fname, patched, envir = sd_imports)
        lockBinding(fname, sd_imports)
      }
    }
    TRUE
  } else {
    FALSE
  }
})
# --- Load object -------------------------------------------------------------
log_info(sprintf("Loading: %s", input_file))
obj <- load_object(input_file)

if (!is_seurat(obj)) {
  stop(sprintf(
    "convert_seuratdisk.R only handles Seurat objects. Got: %s\nUse convert_mtx.R for non-Seurat objects.",
    paste(class(obj), collapse = ", ")
  ))
}

dims <- get_dims(obj)
log_info(sprintf("Object: %d cells x %d features", dims$cells, dims$features))

# --- Compatibility: downgrade Assay5 -> Assay for SeuratDisk -----------------
# SeuratDisk has not been updated for Seurat v5's Assay5 format.
# We downgrade to v3 Assay to work around this.
suppressMessages(library(Seurat))
assay_names <- names(obj@assays)
for (aname in assay_names) {
  assay_obj <- obj[[aname]]
  if (inherits(assay_obj, "Assay5")) {
    log_info(sprintf("Downgrading assay '%s' from Assay5 to Assay for SeuratDisk compatibility...", aname))
    # Extract count matrix and data matrix using v5 API
    counts_mat <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = aname, layer = "counts"),
      error = function(e) {
        SeuratObject::GetAssayData(obj, assay = aname, layer = "data")
      }
    )
    data_mat <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = aname, layer = "data"),
      error = function(e) NULL
    )
    # Replace assay with v3 style
    suppressWarnings({
      obj[[aname]] <- Seurat::CreateAssayObject(counts = counts_mat)
    })
    # Copy data layer back if it existed
    if (!is.null(data_mat)) {
      tryCatch({
        obj[[aname]] <- SeuratObject::SetAssayData(obj[[aname]], layer = "data", new.data = data_mat)
      }, error = function(e) log_warn(sprintf("Could not copy data layer: %s", e$message)))
    }
  }
}

# --- Convert to h5Seurat -----------------------------------------------------
# Create temp file in the output directory (same filesystem for rename)
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tmp_h5seurat <- file.path(output_dir, paste0(".tmp_", basename(tools::file_path_sans_ext(output_file)), ".h5Seurat"))
tmp_h5ad     <- paste0(tools::file_path_sans_ext(tmp_h5seurat), ".h5ad")

# Clean up any leftover temp files
unlink(tmp_h5seurat)
unlink(tmp_h5ad)

log_info("Step 1/2: Saving as h5Seurat intermediate...")
tryCatch({
  SeuratDisk::SaveH5Seurat(obj, filename = tmp_h5seurat, overwrite = TRUE)
}, error = function(e) {
  stop(sprintf("SaveH5Seurat failed: %s\nThis may happen with Seurat v5 objects if SeuratDisk is outdated.\nTry: remotes::install_github('mojaveazure/seurat-disk')", e$message))
})
log_debug(sprintf("h5Seurat saved: %s (%.1f MB)",
                  tmp_h5seurat, file.info(tmp_h5seurat)$size / 1e6))

log_info("Step 2/2: Converting h5Seurat to h5ad...")
tryCatch({
  SeuratDisk::Convert(tmp_h5seurat, dest = "h5ad", overwrite = TRUE)
}, error = function(e) {
  stop(sprintf("Convert to h5ad failed: %s", e$message))
})

# The Convert function creates the h5ad next to the h5Seurat file
# Move it to the desired output path
if (tmp_h5ad != output_file) {
  log_debug(sprintf("Moving %s -> %s", tmp_h5ad, output_file))
  file.rename(tmp_h5ad, output_file)
}

# --- Cleanup ----------------------------------------------------------------
unlink(tmp_h5seurat)
unlink(tmp_h5ad)

# --- Report ------------------------------------------------------------------
output_size_mb <- round(file.info(output_file)$size / 1e6, 2)
log_info(sprintf("Conversion complete: %s (%.1f MB)", output_file, output_size_mb))
log_info(sprintf("Dimensions: %d cells x %d features", dims$cells, dims$features))

cat(sprintf("output=%s\ncells=%d\nfeatures=%d\nsize_mb=%.2f\n",
            output_file, dims$cells, dims$features, output_size_mb))
