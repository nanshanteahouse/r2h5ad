# AGENTS.md — r2h5ad

Shell+R+Python CLI for converting single-cell RDS/QS/Rdata files to h5ad (AnnData).

## Entrypoint

`bash r2h5ad.sh <input> [output] [options]` — single entry point for `.rds`, `.qs`, `.Rdata` files. Auto-detects format by extension. Do not run `R/*.R` scripts directly.

## Architecture

```
r2h5ad.sh  ── detect_format.R  ── Seurat? ──yes── convert_seuratdisk.R  →  .h5Seurat  →  .h5ad
             │                   │
             .rds/.qs ──── load_object()    no
             │                                │
             │                                └── convert_mtx.R  →  matrix.mtx  →  Python assembly  →  .h5ad
             │
             .Rdata ─── convert_rdata.R  →  MTX + metadata  →  assemble_h5ad.py  →  .h5ad
```

**Three conversion paths:**
- **SeuratDisk** (primary for RDS/QS): intermediate `.h5Seurat` file → `SeuratDisk::Convert()` → rename to target. Preserves assays, metadata, reductions.
- **MTX export** (fallback for RDS/QS): writes MTX + features.tsv + barcodes.tsv to temp dir → Python script calls `scanpy.read().T` → attaches metadata CSV → writes h5ad.
- **Rdata export** (for .Rdata/.RData/.rda): calls `convert_rdata.R` to load workspace dump, extract MTX + optional embeddings → `assemble_h5ad.py` to build h5ad.

## Critical quirks

### Assay5 → Assay downgrade + SeuratObject compat patch
SeuratDisk does **not** support Seurat v5's `Assay5` format. `convert_seuratdisk.R` automatically downgrades assay layers in-memory (`Assay5@layers` → `v3 CreateAssayObject`) before `SaveH5Seurat`. This is lossy for multi-layer assays — keep an eye on layer preservation if extending.
Also, SeuratDisk is **unmaintained** and uses `slot=` which was **defunct in SeuratObject 5.3.0+**. `convert_seuratdisk.R` monkey-patches `GetAssayData`/`SetAssayData` in SeuratDisk's import env to translate `slot=` → `layer=`.

### SeuratDisk is GitHub-only
Not on CRAN or conda. Must install after environment setup:
```bash
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
```

### MTX transpose
`convert_mtx.R` line 137 calls `sc.read(...).T` because scanpy reads MTX as features×cells but AnnData expects cells×features. If the h5ad shape is transposed, check here first.

### 2GB RDS limit
`detect_format.R` warns for RDS files >2GB. Use QS format for large objects — it's faster and more stable.

.qs.gz / .rds.gz auto-decompression
`detect_format.R` and `load_object()` in `utils.R` auto-decompress `.qs.gz` / `.rds.gz`. QS uses `gunzip -c` to a temp file; RDS reads directly from a `gzfile` connection.

### Conda auto-detection
The shell script probes these paths in order: `~/miniforge3`, `~/miniconda3`, `~/anaconda3`, `/opt/miniforge3`, `/opt/conda`. Set `R2H5AD_CONDA_ENV` to override the env name. If no conda is found, it falls back to system R/Python on `$PATH`.

### Windows/WSL paths
Both `r2h5ad.sh` (bash) and `R/utils.R` (R) have `win_to_wsl()` and `wsl_to_win()` converters. Input paths like `D:\data\obj.rds` are automatically translated to `/mnt/d/data/obj.rds`.

## R scripts

All in `R/`. Each sources `utils.R` from the same directory. `utils.R` provides:
- `load_object(path)` — dispatches by extension (`.rds`, `.qs`, `.rds.gz`, `.qs.gz`)
- `is_seurat(obj)`, `is_sce(obj)`, `is_summarized_experiment(obj)`
- `get_dims(obj)` — returns `list(cells=, features=)`
- `extract_counts(obj, assay)` — handles Seurat v5 (`layer="counts"`) with v4 fallback (`slot="counts"`)

## Dependencies

| R | Python |
|---|---|
| Seurat (≥5.0), SeuratDisk (GitHub, unmaintained), qs, jsonlite, Matrix, hdf5r | anndata (≥0.12), scanpy (≥1.12) |

## Dev commands

```bash
# Setup
conda env create -f environment.yml && conda activate r2h5ad
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

# Run
bash r2h5ad.sh input.rds output.h5ad

# Verify output
python3 -c "import scanpy; print(scanpy.read_h5ad('output.h5ad'))"

# Debug
bash r2h5ad.sh input.rds --verbose --no-cleanup

# Force MTX path (non-Seurat or if SeuratDisk fails)
bash r2h5ad.sh input.rds --method mtx
```

## No tests, no CI

This is a single-commit tool repo. No test suite, no CI workflows. Verify manually after changes.
