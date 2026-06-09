#!/usr/bin/env python3
"""
Export raw counts from h5ad as CSV for GRNBoost2.
Output: cells x genes CSV, gene names as columns (hyphens → underscores),
cell IDs as index. No unnamed index column.
"""
import argparse
import scanpy as sc
import pandas as pd
import scipy.sparse

def main():
    parser = argparse.ArgumentParser(description='Export raw counts from h5ad to CSV')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--use_raw', action='store_true', default=False,
                        help='Use adata.raw.X instead of adata.X')
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    if args.use_raw and adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    # Convert sparse to dense if needed
    if scipy.sparse.issparse(X):
        X = X.toarray()

    # Hyphen → underscore (Seurat convention → standard)
    gene_names = var_names.str.replace('-', '_')

    df = pd.DataFrame(X, index=adata.obs_names, columns=gene_names)
    df.to_csv(args.output)
    print(f"Exported: {df.shape[0]} cells x {df.shape[1]} genes → {args.output}")

if __name__ == '__main__':
    main()
