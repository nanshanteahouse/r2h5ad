# r2h5ad — RDS / QS / Rdata to h5ad Converter

Convert single-cell data from R-native formats (`.rds`, `.qs`, `.rds.gz`, `.qs.gz`, `.Rdata`, `.RData`) to AnnData `.h5ad`, preserving metadata wherever possible.

## Quick start

```bash
# 1. One-time setup
conda env create -f environment.yml
conda activate r2h5ad
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

# 2. Convert
bash r2h5ad.sh your_data.rds output.h5ad

# 3. Verify
python3 -c "import scanpy; print(scanpy.read_h5ad('output.h5ad'))"
```

## How it works

input.rds/.rds.gz  ──► detect_format.R ──► Seurat? ──yes──► convert_seuratdisk.R ──► .h5Seurat ──► .h5ad
                     │                   │
input.qs/.qs.gz  ────► load_object()       no
                     │                       │
                     │                       └──► convert_mtx.R ──► matrix.mtx ──► .h5ad
                     │                                (universal fallback)
                     │
input.Rdata/.rda  ──► convert_rdata.R ──► MTX + metadata ──► assemble_h5ad.py ──► .h5ad
                              (legacy workspace dump)

| Method | Triggered for | Preserves |
|--------|---------------|-----------|
| **SeuratDisk** (primary for RDS/QS) | Seurat objects | All assays, metadata, dimensions ⭐ |
| **MTX export** (fallback for RDS/QS) | Everything else (SCE, matrix, list) | Count matrix + cell metadata |
| **Rdata export** (for .Rdata) | Workspace dump objects | Count matrix + PCA/t-SNE/UMAP metadata |

## Usage

```bash
bash r2h5ad.sh <input_file> [output_file] [options]
```

### Options

| Flag | Description |
|------|-------------|
| `--method seuratdisk\|mtx` | Force a specific conversion method (serialized objects only) |
| `--force` | Overwrite existing output |
| `--verbose` | Print detailed debug logs |
| `--assay NAME` | Assay to extract (default: `RNA`, for MTX fallback) |
| `--no-cleanup` | Keep temp files on error (for debugging) |
| `--skip-deps-check` | Skip pre-flight dependency verification |
| `--list` | List objects in .Rdata and exit (Rdata only) |
| `--count-object NAME` | Count matrix object name (Rdata only, auto-detected) |
| `--pca-object / --tsne-object / --umap-object NAME` | Embedding object names (Rdata only) |
### Examples

```bash
# Basic Seurat RDS conversion
bash r2h5ad.sh data/obj.rds

# QS file with custom output path
bash r2h5ad.sh data/obj.qs results/processed.h5ad --force

# Force MTX fallback (e.g. for SingleCellExperiment)
bash r2h5ad.sh data/sce.rds --method mtx

# Debug a failing conversion
bash r2h5ad.sh data/obj.rds --verbose --no-cleanup

# Convert .Rdata workspace dump
bash r2h5ad.sh data/obj.Rdata

# List objects inside an .Rdata file
bash r2h5ad.sh data/obj.Rdata --list

# Specify count matrix in Rdata with embeddings
bash r2h5ad.sh data/obj.Rdata --count-object expr --pca-object pca_result

# From Windows host (via WSL)
wsl bash D:/Projects/r2h5ad/r2h5ad.sh D:/data/obj.rds D:/output.h5ad

## Environment & dependencies

### Option A — Conda (recommended)

```bash
conda env create -f environment.yml
conda activate r2h5ad
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
```

If you already have a conda env with the needed packages, point the tool at it:

```bash
export R2H5AD_CONDA_ENV=my_existing_env
bash r2h5ad.sh ...
```

The tool automatically finds conda installations at common paths (`~/miniforge3`, `~/miniconda3`, `~/anaconda3`, `/opt/conda`).

### Option B — Manual (no conda)

R packages (any R ≥4.0):
```r
install.packages(c("Seurat", "qs", "jsonlite", "Matrix"))
remotes::install_github("mojaveazure/seurat-disk")
```

Python packages:
```bash
pip install anndata scanpy
```

The tool only requires `Rscript` and `python3` on `$PATH` — no conda dependency.

### Required packages

| Ecosystem | Package | Purpose |
|-----------|---------|---------|
| R | `Seurat` (≥5.0) | Object loading / assay handling |
| R | `SeuratDisk` | Primary conversion path (Seurat → h5ad) |
| R | `qs` | QS file format support |
| R | `jsonlite` | JSON output for format detection |
| R | `Matrix` | Sparse matrix handling |
| R | `hdf5r` | HDF5 backend (SeuratDisk dependency) |
| Python | `anndata` (≥0.12) | h5ad read/write |
| Python | `scanpy` (≥1.12) | MTX loading (fallback path) |

| Object class / format | Primary method | Fallback |
|-------------|---------------|----------|
| `Seurat` (v3/v4/v5) | SeuratDisk | MTX export |
| `SingleCellExperiment` | — | MTX export |
| `SummarizedExperiment` | — | MTX export |
| `dgCMatrix` / `Matrix` (raw) | — | MTX export |
| Named list (with `$counts`) | — | MTX export |
| `.Rdata` / `.RData` / `.rda` workspace dump | Rdata export (MTX) | — |
Seurat v5 `Assay5` objects are automatically downgraded to v3 `Assay` for SeuratDisk compatibility.

## Troubleshooting

**"Missing R packages: Seurat SeuratDisk qs"**
→ Run the dependency setup (Option A or B above), or use `--skip-deps-check` if you're sure they're installed.

**"SaveH5Seurat failed: slot deprecated in SeuratObject 5.3.0+"**
→ SeuratObject defuncted the `slot=` argument. The tool auto-patches SeuratDisk's internal calls, but if errors persist reinstall SeuratDisk:
`Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'`

**"Rscript not found on PATH"**
→ R is not installed or not on `$PATH`. Install R (≥4.0) or activate a conda environment.

**Large file warnings (>2GB)**
→ RDS has an internal 2GB stability limit. Use QS format for large objects; the tool handles both identically.

**File not found on Windows**
→ Use forward slashes: `wsl bash D:/Projects/r2h5ad/r2h5ad.sh ...`

**".qs.gz / .rds.gz support"**
→ The tool auto-decompresses `.qs.gz` and `.rds.gz` files. Use them directly as input (no manual decompression needed).

## Related

- [SeuratDisk](https://mojaveazure.github.io/seurat-disk/)
- [Scanpy file formats](https://scanpy.readthedocs.io/en/stable/api/scanpy.read.html)
- Pipeline reference: `rds_format_support.md` (in this directory)
