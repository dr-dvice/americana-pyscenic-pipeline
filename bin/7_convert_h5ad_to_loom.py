#!/usr/bin/env python3
"""
Convert h5ad to loom format for pySCENIC.
Gene names: hyphens → underscores (Seurat → standard format).
"""
import argparse
import scanpy as sc
import scipy.sparse
import numpy as np
import loompy

def main():
    parser = argparse.ArgumentParser(description='Convert h5ad to loom for pySCENIC')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output', required=True, help='Output loom file')
    args = parser.parse_args()

    print(f"Loading h5ad: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"Shape: {adata.shape} (cells x genes)")

    adata.var_names = adata.var_names.str.replace('-', '_')
    print("Converted gene names: hyphens -> underscores")

    # Convert sparse to dense for loompy (requires ndarray)
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    row_attrs = {"Gene": adata.var_names.values}
    col_attrs = {"CellID": adata.obs_names.values}

    loompy.create(args.output, X.T, row_attrs, col_attrs)
    print(f"Saved: {args.output}")

if __name__ == '__main__':
    main()
