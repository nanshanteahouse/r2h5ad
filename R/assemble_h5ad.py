#!/usr/bin/env python3
"""
assemble_h5ad.py — Build h5ad from Rdata conversion intermediates

Reads MTX + CSV files produced by convert_rdata.R and assembles an AnnData h5ad.
Called by rdata2h5ad.sh.

Usage:
    python3 assemble_h5ad.py <work_dir> <output.h5ad>
"""
import sys
import os
import json
import pandas as pd
import numpy as np
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix


def assemble(work_dir: str, output_h5ad: str):
    """Assemble h5ad from intermediates in work_dir."""

    # ── Read MTX ───────────────────────────────────────────────────────
    mtx_path = os.path.join(work_dir, "matrix.mtx")
    if not os.path.exists(mtx_path):
        raise FileNotFoundError(f"Missing matrix.mtx in {work_dir}")
    counts = mmread(mtx_path).T.tocsr()
    print(f"[INFO] Matrix: {counts.shape} ({counts.nnz} non-zero)")

    # ── Read barcodes and features ─────────────────────────────────────
    barcodes_path = os.path.join(work_dir, "barcodes.tsv")
    features_path = os.path.join(work_dir, "features.tsv")

    if os.path.exists(barcodes_path):
        barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)[0].tolist()
    elif os.path.exists(os.path.join(work_dir, "barcodes.csv")):
        barcodes = pd.read_csv(os.path.join(work_dir, "barcodes.csv"), header=None)[0].tolist()
    else:
        barcodes = [f"cell_{i}" for i in range(counts.shape[0])]

    if os.path.exists(features_path):
        genes = pd.read_csv(features_path, sep="\t", header=None)[0].tolist()
    elif os.path.exists(os.path.join(work_dir, "features.csv")):
        genes = pd.read_csv(os.path.join(work_dir, "features.csv"), header=None)[0].tolist()
    else:
        genes = [f"gene_{i}" for i in range(counts.shape[1])]

    # Verify dimensions match
    n_cells, n_genes = counts.shape
    if len(barcodes) != n_cells:
        print(f"[WARN] barcodes count ({len(barcodes)}) != matrix cells ({n_cells})")
        if len(barcodes) > n_cells:
            barcodes = barcodes[:n_cells]
        else:
            barcodes = barcodes + [f"cell_{i}" for i in range(len(barcodes), n_cells)]
    if len(genes) != n_genes:
        print(f"[WARN] features count ({len(genes)}) != matrix genes ({n_genes})")
        if len(genes) > n_genes:
            genes = genes[:n_genes]
        else:
            genes = genes + [f"gene_{i}" for i in range(len(genes), n_genes)]

    print(f"[INFO] Cells: {len(barcodes)}, Genes: {len(genes)}")

    # ── Create AnnData ─────────────────────────────────────────────────
    adata = ad.AnnData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    # Ensure var_names are unique
    adata.var_names_make_unique()
    # Raw counts
    adata.raw = adata

    # ── Attach obsm (PCA, t-SNE, UMAP) ────────────────────────────────
    obsm_files = {
        "X_pca": "obsm_pca.csv",
        "X_tsne": "obsm_tsne.csv",
        "X_umap": "obsm_umap.csv",
    }

    for obsm_key, csv_name in obsm_files.items():
        csv_path = os.path.join(work_dir, csv_name)
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path, index_col=0)
            # Align rows to our barcode order
            idx = pd.Index(barcodes)
            matched = df.index.intersection(idx)
            if len(matched) == 0:
                print(f"[WARN] No matching cells for {csv_name}, skipping")
                continue
            reindexed = df.reindex(barcodes).values.astype(np.float32)
            # Handle NaN (cells without metadata)
            if np.any(np.isnan(reindexed)):
                print(f"[WARN] {np.isnan(reindexed).sum()} NaN values in {obsm_key}, filling with 0")
                reindexed = np.nan_to_num(reindexed, nan=0.0)
            adata.obsm[obsm_key] = reindexed
            print(f"[INFO] Added obsm '{obsm_key}': {reindexed.shape}")

    # ── Attach obs metadata ────────────────────────────────────────────
    cell_meta_path = os.path.join(work_dir, "cell_meta.csv")
    if os.path.exists(cell_meta_path):
        meta_df = pd.read_csv(cell_meta_path, index_col=0)
        for col in meta_df.columns:
            if col in adata.obs:
                continue
            aligned = meta_df[col].reindex(barcodes)
            if aligned.isna().sum() < len(barcodes) * 0.5:  # at least 50% filled
                adata.obs[col] = aligned.values
                print(f"[INFO] Added obs column '{col}' ({aligned.notna().sum()} / {len(barcodes)} cells)")

    # ── Write summary for bash wrapper ─────────────────────────────────
    summary = {
        "shape": list(adata.shape),
        "obsm_keys": list(adata.obsm.keys()),
        "obs_columns": list(adata.obs.columns),
    }

    # ── Save ───────────────────────────────────────────────────────────
    print(f"[INFO] Writing {output_h5ad}...")
    adata.write(output_h5ad, compression="gzip")
    print(f"[INFO] Done: {output_h5ad}")

    # Write summary JSON
    summary_path = os.path.join(work_dir, "h5ad_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f)

    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 assemble_h5ad.py <work_dir> <output.h5ad>")
        sys.exit(1)
    sys.exit(assemble(sys.argv[1], sys.argv[2]))
