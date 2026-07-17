#!/usr/bin/env Rscript
# =============================================================================
# convert_h5ad_seuratdisk.R — h5ad → Seurat via SeuratDisk (reverse path)
# =============================================================================
#
# Converts a .h5ad AnnData file to an .rds or .qs Seurat object using the
# SeuratDisk intermediate format:
#   h5ad  ──Convert()──►  .h5Seurat  ──LoadH5Seurat()──►  Seurat  ──►  .rds/.qs
#
# Usage:
#   Rscript convert_h5ad_seuratdisk.R <input.h5ad> <output.rds|output.qs>
#
# NOTE: LoadH5Seurat calls GetAssayData(slot=) which is defunct in
# SeuratObject >= 5.3.0. The monkey-patch below handles SaveH5Seurat
# (forward) but LoadH5Seurat's internal AssembleAssay path uses S3 dispatch
# to GetAssayData.Assay which bypasses the patched generic. If it fails,
# the caller (h5ad2r.sh) falls back to the MTX intermediate path.
# =============================================================================

script_dir <- dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
utils_path <- file.path(script_dir, "utils.R")
if (file.exists(utils_path)) source(utils_path)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert_h5ad_seuratdisk.R <input.h5ad> <output.rds|output.qs>")
}
input_file  <- args[1]
output_file <- args[2]

out_ext <- tolower(tools::file_ext(output_file))

# --- SeuratDisk monkey-patch (same as convert_seuratdisk.R) --------------------
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

    assignInNamespace("GetAssayData", patched_get, ns = "SeuratObject")
    assignInNamespace("SetAssayData", patched_set, ns = "SeuratObject")

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

# --- Convert h5ad → h5Seurat --------------------------------------------------
output_dir <- dirname(output_file)
tmp_h5seurat <- file.path(output_dir, paste0(".tmp_", basename(tools::file_path_sans_ext(output_file)), ".h5Seurat"))

log_info(sprintf("Converting %s -> h5Seurat...", basename(input_file)))

tryCatch({
  SeuratDisk::Convert(
    source    = input_file,
    dest      = "h5seurat",
    assay     = "RNA",
    overwrite = TRUE,
    verbose   = FALSE
  )

  expected_h5seurat <- paste0(tools::file_path_sans_ext(input_file), ".h5seurat")
  if (file.exists(expected_h5seurat)) {
    file.rename(expected_h5seurat, tmp_h5seurat)
    actual_h5seurat <- tmp_h5seurat
  } else if (file.exists(tmp_h5seurat)) {
    actual_h5seurat <- tmp_h5seurat
  } else {
    stop("h5Seurat file was not created")
  }
}, error = function(e) {
  stop(sprintf("Convert(h5ad→h5Seurat) failed: %s", e$message))
})

log_info(sprintf("h5Seurat created: %s (%.1f MB)",
                  basename(actual_h5seurat),
                  file.info(actual_h5seurat)$size / 1e6))

# --- Load Seurat object -------------------------------------------------------
log_info("Loading h5Seurat into Seurat object...")
obj <- tryCatch(
  SeuratDisk::LoadH5Seurat(actual_h5seurat),
  error = function(e) {
    unlink(actual_h5seurat)
    stop(sprintf("LoadH5Seurat failed: %s", e$message))
  }
)

log_info(sprintf("Seurat object: %d cells x %d features",
                  ncol(obj), nrow(obj)))
log_info(sprintf("Assays: %s", paste(names(obj@assays), collapse = ", ")))

# --- Save output --------------------------------------------------------------
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

# --- Cleanup ------------------------------------------------------------------
unlink(actual_h5seurat)

# --- JSON summary -------------------------------------------------------------
summary <- list(
  status       = "success",
  method       = "seuratdisk",
  output       = output_file,
  output_mb    = output_mb,
  cells        = ncol(obj),
  features     = nrow(obj),
  assays       = as.list(names(obj@assays)),
  reductions   = as.list(names(obj@reductions))
)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = FALSE))
