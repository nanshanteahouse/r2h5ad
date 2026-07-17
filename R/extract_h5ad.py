#!/usr/bin/env python3
"""
extract_h5ad.py — Export AnnData h5ad to MTX + CSV intermediates

Reads h5ad with scanpy/anndata and writes intermediate files for
R-side assembly (assemble_seurat.R / assemble_rdata.R).

Usage:
    python3 R/extract_h5ad.py <input.h5ad> <outdir> [--raw] [--assay RNA]

Exports:
    matrix.mtx       — count matrix (features × cells, Market Matrix)
    barcodes.tsv     — cell barcodes
    features.tsv     — feature/gene names
    cell_meta.csv    — obs metadata
    var_meta.csv     — var metadata
    obsm_pca.csv     — PCA embeddings (if present)
    obsm_umap.csv    — UMAP embeddings (if present)
    obsm_tsne.csv    — t-SNE embeddings (if present)
    summary.json     — shape, keys, column info
"""
import sys
import os
import json
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import issparse, coo_matrix


def main():
    parser = argparse.ArgumentParser(
        description="Export AnnData h5ad to MTX+CSV intermediates for R assembly"
    )
    parser.add_argument("input", help="Input .h5ad file")
    parser.add_argument("outdir", help="Output directory for intermediates")
    parser.add_argument("--processed", action="store_true",
                        help="Use adata.X instead of adata.raw.X (default: use raw.X if available)")
    parser.add_argument("--assay", default="RNA",
                        help="Assay name for summary (default: RNA)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    print(f"[INFO] Loaded: {args.input} — {adata.shape[0]} cells × {adata.shape[1]} genes")

    # ── Determine X matrix ─────────────────────────────────────────────────
    if args.processed or adata.raw is None:
        X = adata.X
        print(f"[INFO] Using adata.X{' (no raw.X available)' if adata.raw is None else ' (--processed)'}")
    else:
        X = adata.raw.X
        print("[INFO] Using adata.raw.X")
    n_cells, n_genes = adata.shape
    nnz = X.nnz if issparse(X) else np.count_nonzero(X)
    print(f"[INFO] Matrix: {n_cells} cells × {n_genes} genes, {nnz} non-zero")

    # ── Write MTX (transpose to features×cells for R) ──────────────────────
    mtx_path = os.path.join(args.outdir, "matrix.mtx")
    X_t = X.T  # genes × cells
    if not issparse(X_t):
        X_t = coo_matrix(X_t)
    else:
        X_t = X_t.tocoo()
    mmwrite(mtx_path, X_t)
    mtx_mb = round(os.path.getsize(mtx_path) / 1e6, 2)
    print(f"[INFO] MTX written: {mtx_mb} MB")

    # ── Write barcodes.tsv ─────────────────────────────────────────────────
    barcodes = list(adata.obs_names)
    barcodes_path = os.path.join(args.outdir, "barcodes.tsv")
    with open(barcodes_path, "w") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    # ── Write features.tsv ─────────────────────────────────────────────────
    features = list(adata.var_names)
    features_path = os.path.join(args.outdir, "features.tsv")
    with open(features_path, "w") as f:
        for ft in features:
            f.write(f"{ft}\n")

    print(f"[INFO] Cells: {len(barcodes)}, Genes: {len(features)}")

    # ── Write cell_meta.csv (obs) ──────────────────────────────────────────
    if adata.obs.shape[1] > 0:
        meta_path = os.path.join(args.outdir, "cell_meta.csv")
        adata.obs.to_csv(meta_path, index=True)
        print(f"[INFO] cell_meta.csv: {len(adata.obs.columns)} columns × {len(adata.obs)} cells")

    # ── Write var_meta.csv (var) ───────────────────────────────────────────
    if adata.var.shape[1] > 0:
        var_path = os.path.join(args.outdir, "var_meta.csv")
        adata.var.to_csv(var_path, index=True)
        print(f"[INFO] var_meta.csv: {len(adata.var.columns)} columns × {len(adata.var)} genes")

    # ── Write obsm (reductions) ────────────────────────────────────────────
    obsm_exported = []
    obsm_map = {
        "X_pca": ("obsm_pca.csv", "PC_"),
        "X_umap": ("obsm_umap.csv", "UMAP_"),
        "X_tsne": ("obsm_tsne.csv", "tSNE_"),
    }

    for obsm_key, (csv_name, col_prefix) in obsm_map.items():
        if obsm_key not in adata.obsm:
            # also check alternative keys (sometimes stored lowercase or differently)
            alt_keys = [k for k in adata.obsm.keys() if k.lower().replace("_", "") == obsm_key.lower().replace("_", "")]
            if alt_keys:
                obsm_key = alt_keys[0]
            else:
                continue

        mat = adata.obsm[obsm_key]
        n_cols = mat.shape[1]
        colnames = [f"{col_prefix}{i+1}" for i in range(n_cols)]
        df = pd.DataFrame(mat, index=barcodes, columns=colnames)
        csv_path = os.path.join(args.outdir, csv_name)
        df.to_csv(csv_path, index=True)
        obsm_exported.append(obsm_key)
        print(f"[INFO] {csv_name}: {df.shape}")

    # ── Write summary.json ─────────────────────────────────────────────────
    summary = {
        "source": os.path.abspath(args.input),
        "shape": [n_cells, n_genes],
        "nnz": int(nnz),
        "has_raw": adata.raw is not None,
        "assay": args.assay,
        "obsm_keys": list(adata.obsm.keys()),
        "obsm_exported": obsm_exported,
        "obs_columns": list(adata.obs.columns) if adata.obs.shape[1] > 0 else [],
        "var_columns": list(adata.var.columns) if adata.var.shape[1] > 0 else [],
        "mtx_mb": mtx_mb,
    }
    summary_path = os.path.join(args.outdir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] summary.json written")

    return 0


if __name__ == "__main__":
    sys.exit(main())
