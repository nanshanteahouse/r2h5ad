# r2h5ad — Bidirectional single-cell format converter

Convert between R-native single-cell formats (`.rds`, `.qs`, `.Rdata`) and AnnData `.h5ad`, preserving metadata wherever possible.

## Quick start

```bash
# 1. One-time setup
conda env create -f environment.yml
conda activate r2h5ad
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

# 2. Forward: R → h5ad
bash r2h5ad.sh your_data.rds output.h5ad

# 3. Reverse: h5ad → R
bash h5ad2r.sh output.h5ad roundtrip.rds

# 4. Verify
python3 -c "import scanpy; print(scanpy.read_h5ad('output.h5ad'))"
Rscript -e "obj <- readRDS('roundtrip.rds'); cat(class(obj)[1], ncol(obj), 'cells x', nrow(obj), 'features\n')"
```

## How it works

### Forward: R → h5ad (`r2h5ad.sh`)

```
input.rds/.rds.gz  ──► detect_format.R ──► Seurat? ──yes──► convert_seuratdisk.R ──► .h5Seurat ──► .h5ad
                     │                   │
input.qs/.qs.gz  ────► load_object()       no
                     │                       │
                     │                       └──► convert_mtx.R ──► matrix.mtx ──► .h5ad
                     │                                (universal fallback)
                     │
input.Rdata/.rda  ──► convert_rdata.R ──► MTX + metadata ──► assemble_h5ad.py ──► .h5ad
                              (legacy workspace dump)
```

### Reverse: h5ad → R (`h5ad2r.sh`)

```
input.h5ad  ──►  h5ad2r.sh
                 │
                 ├── .rds/.qs output:  convert_h5ad_seuratdisk.R  ──►  .h5Seurat  ──►  Seurat  ──►  .rds/.qs
                 │                     (SeuratDisk primary path, falls back to MTX on failure)
                 │
                 └── .Rdata output:  extract_h5ad.py  ──►  MTX+CSV  ──►  assemble_rdata.R  ──►  .Rdata
```

| Direction | Method | Triggered for | Preserves |
|-----------|--------|---------------|-----------|
| R → h5ad | **SeuratDisk** | Seurat objects | All assays, metadata, reductions ⭐ |
| R → h5ad | **MTX export** | Everything else (SCE, matrix, list) | Count matrix + cell metadata |
| R → h5ad | **Rdata export** | .Rdata workspace dump | Count matrix + PCA/t-SNE/UMAP metadata |
| h5ad → R | **SeuratDisk reverse** | .rds/.qs output (primary) | Counts, data, metadata, reductions ⭐ |
| h5ad → R | **MTX intermediate** | .rds/.qs output (fallback) | Count matrix + metadata + reductions |
| h5ad → R | **Rdata export** | .Rdata output | Count matrix + metadata + reductions |

### Information preservation (round-trip)

| Component | Forward (R→h5ad) | Reverse (h5ad→R) |
|-----------|-----------------|-------------------|
| Raw counts (`raw.X`) | ✅ Stored in `raw.X` | ✅ Maps to Seurat `counts` |
| Normalized data (`X`) | ✅ Stored in `X` | ✅ Maps to Seurat `data` |
| Scaled data | ✅ (if dense) | ✅ (if dense `X`) |
| Cell metadata (`obs`) | ✅ | ✅ → Seurat `meta.data` |
| Feature metadata (`var`) | ✅ | ✅ → Seurat `meta.features` |
| PCA/UMAP/t-SNE (`obsm`) | ✅ | ✅ → Seurat reductions |
| NN graphs (`obsp`) | ✅ | ✅ → Seurat graphs |
| Seurat commands/tools log | ❌ Lost | ❌ Not recovered |
| `uns` non-standard data | ⚠️ Partial | ⚠️ Partial |
| Assay5 multi-layers | ❌ Downgraded to v3 | ❌ Not recovered |

## Usage

### R → h5ad (`r2h5ad.sh`)

```bash
bash r2h5ad.sh <input_file> [output_file] [options]
```

Options: `--method seuratdisk|mtx`, `--force`, `--verbose`, `--assay NAME`, `--no-cleanup`, `--skip-deps-check`, `--list`, `--count-object`, `--pca-object`, `--tsne-object`, `--umap-object`

### h5ad → R (`h5ad2r.sh`)

```bash
bash h5ad2r.sh <input.h5ad> [output.rds|output.qs|output.Rdata] [options]
```

| Flag | Description |
|------|-------------|
| `--method seuratdisk\|mtx` | Force conversion method (default: auto, SeuratDisk first) |
| `--force` | Overwrite existing output |
| `--verbose` | Print detailed debug logs |
| `--processed` | Use `adata.X` instead of `adata.raw.X` (default: prefer raw.X) |
| `--assay NAME` | Assay name (default: `RNA`) |
| `--no-cleanup` | Keep temp files on error (for debugging) |
| `--skip-deps-check` | Skip pre-flight dependency verification |

Output format is auto-detected by extension: `.rds` → Seurat RDS, `.qs` → Seurat QS, `.Rdata`/.RData/.rda → workspace dump.

### Examples

```bash
# Forward conversions
bash r2h5ad.sh data/obj.rds
bash r2h5ad.sh data/obj.qs results/processed.h5ad --force
bash r2h5ad.sh data/sce.rds --method mtx                      # Force MTX fallback
bash r2h5ad.sh data/obj.Rdata --list                           # List objects in .Rdata
bash r2h5ad.sh data/obj.Rdata --count-object expr --pca-object pca_result

# Reverse conversions
bash h5ad2r.sh data/processed.h5ad                             # .rds output (default)
bash h5ad2r.sh data/obj.h5ad output.qs --force                 # .qs output
bash h5ad2r.sh data/obj.h5ad output.Rdata --processed         # .Rdata with processed X
bash h5ad2r.sh data/obj.h5ad --method mtx --verbose            # Force MTX path

# Round-trip
bash r2h5ad.sh original.rds tmp.h5ad
bash h5ad2r.sh tmp.h5ad roundtrip.rds

# Debug
bash r2h5ad.sh data/obj.rds --verbose --no-cleanup
bash h5ad2r.sh data/obj.h5ad --verbose --no-cleanup

# From Windows host (via WSL)
wsl bash D:/Projects/r2h5ad/r2h5ad.sh D:/data/obj.rds D:/output.h5ad
wsl bash D:/Projects/r2h5ad/h5ad2r.sh D:/data/obj.h5ad D:/output.rds
```

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
bash h5ad2r.sh ...
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

### Required packages

| Ecosystem | Package | Purpose |
|-----------|---------|---------|
| R | `Seurat` (≥5.0) | Object loading / assay handling |
| R | `SeuratDisk` | Primary conversion path (bidirectional h5Seurat) |
| R | `qs` | QS file format support |
| R | `jsonlite` | JSON output |
| R | `Matrix` | Sparse matrix handling |
| R | `hdf5r` | HDF5 backend (SeuratDisk dependency) |
| Python | `anndata` (≥0.12) | h5ad read/write |
| Python | `scanpy` (≥1.12) | MTX loading (fallback path) |

### Supported formats

| Object class / format | Forward (R→h5ad) | Reverse (h5ad→R) |
|-------------|:---:|:---:|
| `Seurat` (v3/v4/v5) | SeuratDisk / MTX | — |
| `SingleCellExperiment` | MTX export | — |
| `SummarizedExperiment` | MTX export | — |
| `dgCMatrix` / `Matrix` (raw) | MTX export | — |
| Named list (with `$counts`) | MTX export | — |
| `.Rdata` workspace dump | Rdata export (MTX) | — |
| `.h5ad` (AnnData) | — | SeuratDisk / MTX / Rdata |

Seurat v5 `Assay5` objects are automatically downgraded to v3 `Assay` for SeuratDisk compatibility.

## Troubleshooting

**"Missing R packages: Seurat SeuratDisk qs"**
→ Run the dependency setup (Option A or B above), or use `--skip-deps-check`.

**"slot deprecated in SeuratObject 5.3.0+"**
→ The tool auto-patches SeuratDisk's internal calls. If errors persist, reinstall:
`Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'`

**"Rscript not found on PATH"**
→ Install R (≥4.0) or activate a conda environment.

**Large file warnings (>2GB)**
→ RDS has an internal 2GB stability limit. Use QS format for large objects.

**File not found on Windows**
→ Use forward slashes: `wsl bash D:/Projects/r2h5ad/r2h5ad.sh ...`

**Reverse conversion fails with SeuratDisk**
→ Try `--method mtx` to use the MTX intermediate path instead.

## Related

- [SeuratDisk](https://mojaveazure.github.io/seurat-disk/)
- [Scanpy file formats](https://scanpy.readthedocs.io/en/stable/api/scanpy.read.html)
- Detailed architecture: `AGENTS.md`
