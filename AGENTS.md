# AGENTS.md — r2h5ad

Shell+R+Python CLI for bidirectional conversion between R-native single-cell formats (RDS/QS/Rdata) and h5ad (AnnData).

## Entrypoints

| Direction | Command | Input | Output |
|---|---|---|---|
| R → h5ad | `bash r2h5ad.sh <input> [output] [options]` | `.rds`, `.qs`, `.rds.gz`, `.qs.gz`, `.Rdata`, `.RData`, `.rda` | `.h5ad` |
| h5ad → R | `bash h5ad2r.sh <input> [output] [options]` | `.h5ad` | `.rds`, `.qs`, `.Rdata`, `.RData`, `.rda` |

Auto-detects format by extension. Do not run `R/*.R` scripts directly.

## Forward architecture (R → h5ad)

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

## Reverse architecture (h5ad → R)

```
h5ad2r.sh ── output ext? ─── .rds/.qs ─── convert_h5ad_seuratdisk.R  →  .h5Seurat  →  Seurat  →  .rds/.qs
           │               │
           │               │  (fallback on SeuratDisk failure)
           │               └── extract_h5ad.py  →  MTX+CSV  →  assemble_seurat.R  →  .rds/.qs
           │
           └── .Rdata ─── extract_h5ad.py  →  MTX+CSV  →  assemble_rdata.R  →  .Rdata
```

**Three reverse paths:**
- **SeuratDisk reverse** (primary for .rds/.qs): `SeuratDisk::Convert(h5ad → h5Seurat)` → `LoadH5Seurat()` → `saveRDS()` / `qs::qsave()`. Same SeuratDisk monkey-patch applies.
- **MTX intermediate** (fallback): `extract_h5ad.py` reads h5ad, exports MTX + metadata CSVs → `assemble_seurat.R` builds Seurat object → saves as .rds/.qs.
- **Rdata export** (for .Rdata/.RData/.rda): `extract_h5ad.py` → MTX + metadata → `assemble_rdata.R` dumps named objects (counts, cell_meta, pca, umap, tsne, etc.) into .Rdata workspace.

**Information preservation (h5ad → R):**
- `X` → `data`; `raw.X` → `counts`; dense `X` → `scale.data`
- `obs` → `meta.data`; `var` → `meta.features` (feature metadata)
- `obsm['X_pca', 'X_umap', 'X_tsne']` → Seurat reductions
- `obsp` → Seurat graphs
- Lost: `uns` non-standard data, original Seurat commands/tools log, Assay5 multi-layers (already downgraded in forward)

## Critical quirks

### Assay5 → Assay downgrade + SeuratObject compat patch
SeuratDisk does **not** support Seurat v5's `Assay5` format. `convert_seuratdisk.R` automatically downgrades assay layers in-memory (`Assay5@layers` → `v3 CreateAssayObject`) before `SaveH5Seurat`. This is lossy for multi-layer assays — keep an eye on layer preservation if extending.
Also, SeuratDisk is **unmaintained** and uses `slot=` which was **defunct in SeuratObject 5.3.0+**. `convert_seuratdisk.R` monkey-patches `GetAssayData`/`SetAssayData` in SeuratDisk's import env to translate `slot=` → `layer=`.

### Same patch needed for LoadH5Seurat (reverse path)
`convert_h5ad_seuratdisk.R` applies the **identical** `slot=` → `layer=` monkey-patch before calling `SeuratDisk::Convert()` and `LoadH5Seurat()`. Both operations trigger internal SeuratDisk function calls that use the defunct `slot=` argument.

### SeuratDisk is GitHub-only
Not on CRAN or conda. Must install after environment setup:
```bash
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'
```

### MTX transpose
**Forward:** `convert_mtx.R` line 137 calls `sc.read(...).T` because scanpy reads MTX as features×cells but AnnData expects cells×features.
**Reverse:** `extract_h5ad.py` transposes `X.T` before writing MTX, because R expects features×cells. If the h5ad shape is transposed, check here first.

### 2GB RDS limit
`detect_format.R` warns for RDS files >2GB. Use QS format for large objects — it's faster and more stable.

### .qs.gz / .rds.gz auto-decompression (forward only)
`detect_format.R` and `load_object()` in `utils.R` auto-decompress `.qs.gz` / `.rds.gz`. QS uses `gunzip -c` to a temp file; RDS reads directly from a `gzfile` connection.

### Conda auto-detection
The shell script probes these paths in order: `~/miniforge3`, `~/miniconda3`, `~/anaconda3`, `/opt/miniforge3`, `/opt/conda`. Set `R2H5AD_CONDA_ENV` to override the env name. If no conda is found, it falls back to system R/Python on `$PATH`.

### Windows/WSL paths
Both shell scripts (bash) and `R/utils.R` (R) have `win_to_wsl()` and `wsl_to_win()` converters. Input paths like `D:\data\obj.rds` are automatically translated to `/mnt/d/data/obj.rds`.

## R scripts

All in `R/`. Each sources `utils.R` from the same directory. `utils.R` provides:
- `load_object(path)` — dispatches by extension (`.rds`, `.qs`, `.rds.gz`, `.qs.gz`)
- `is_seurat(obj)`, `is_sce(obj)`, `is_summarized_experiment(obj)`
- `get_dims(obj)` — returns `list(cells=, features=)`
- `extract_counts(obj, assay)` — handles Seurat v5 (`layer="counts"`) with v4 fallback (`slot="counts"`)

### Forward scripts
| Script | Purpose |
|---|---|
| `detect_format.R` | Inspect RDS/QS object class, assays, dims |
| `convert_seuratdisk.R` | Seurat → .h5Seurat → .h5ad |
| `convert_mtx.R` | Any R object → MTX → h5ad (universal fallback) |
| `convert_rdata.R` | .Rdata workspace → MTX + metadata |
| `assemble_h5ad.py` | Python: MTX + CSV → h5ad |

### Reverse scripts
| Script | Purpose |
|---|---|
| `convert_h5ad_seuratdisk.R` | h5ad → .h5Seurat → Seurat → .rds/.qs |
| `extract_h5ad.py` | Python: h5ad → MTX + metadata CSVs |
| `assemble_seurat.R` | MTX + CSV → Seurat → .rds/.qs |
| `assemble_rdata.R` | MTX + CSV → named objects → .Rdata |

## Dependencies

| R | Python |
|---|---|
| Seurat (≥5.0), SeuratDisk (GitHub, unmaintained), qs, jsonlite, Matrix, hdf5r | anndata (≥0.12), scanpy (≥1.12) |

## Dev commands

```bash
# Setup
conda env create -f environment.yml && conda activate r2h5ad
Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

# Forward: R → h5ad
bash r2h5ad.sh input.rds output.h5ad

# Reverse: h5ad → R
bash h5ad2r.sh input.h5ad output.rds

# Verify forward output
python3 -c "import scanpy; print(scanpy.read_h5ad('output.h5ad'))"

# Verify reverse output
Rscript -e "obj <- readRDS('output.rds'); cat(class(obj)[1], ncol(obj), 'cells x', nrow(obj), 'features\n')"

# Debug
bash r2h5ad.sh input.rds --verbose --no-cleanup
bash h5ad2r.sh input.h5ad --method mtx --verbose --no-cleanup

# Force MTX path (non-Seurat or if SeuratDisk fails)
bash r2h5ad.sh input.rds --method mtx
bash h5ad2r.sh input.h5ad --method mtx

# Round-trip test
bash r2h5ad.sh original.rds tmp.h5ad
bash h5ad2r.sh tmp.h5ad roundtrip.rds
```

## No tests, no CI

This is a single-commit tool repo. No test suite, no CI workflows. Verify manually after changes.
