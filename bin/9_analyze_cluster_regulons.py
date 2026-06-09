#!/usr/bin/env python3
"""
Analyze regulon activity at cluster level from pySCENIC AUCell output.
Generates cluster activity CSV, top regulons per cluster, and heatmaps.
"""
import argparse
import os
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

ROMAN = ['I','II','III','IV','V','VI','VII','VIII','IX','X',
         'XI','XII','XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX']

def make_display_names(ids, fly_map, mean_auc=None):
    gene_ids  = [r.replace('(+)', '').replace('(-)', '') for r in ids]
    raw_names = [fly_map.get(g, g) for g in gene_ids]
    groups = defaultdict(list)
    for i, name in enumerate(raw_names):
        groups[name].append(i)
    result = list(raw_names)
    for name, idxs in groups.items():
        if len(idxs) > 1:
            if mean_auc:
                idxs = sorted(idxs, key=lambda i: mean_auc.get(ids[i], 0), reverse=True)
            for rank, i in enumerate(idxs):
                result[i] = f"{name}-{ROMAN[rank]}"
    return {ids[i]: result[i] for i in range(len(ids))}

def save_heatmap(cluster_regulons, regulon_ids, display, title, filepath):
    plot_df = cluster_regulons[regulon_ids].rename(columns=display).T
    plot_df = plot_df.loc[plot_df.mean(axis=1).sort_values(ascending=True).index]
    plt.figure(figsize=(20, 14))
    sns.heatmap(plot_df, cmap='RdYlBu_r', center=0,
                cbar_kws={'label': 'AUC Score'}, yticklabels=True, xticklabels=True)
    plt.title(title, fontsize=14)
    plt.xlabel('Cluster', fontsize=12)
    plt.ylabel('Regulon', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze regulon activity per cluster')
    parser.add_argument('--auc_loom', required=True, help='AUCell output loom file')
    parser.add_argument('--h5ad', required=True, help='Atlas h5ad file')
    parser.add_argument('--cluster_names', required=True, help='Cluster names CSV (cluster, cell_type)')
    parser.add_argument('--orthologs', required=True, help='Species-to-fly ortholog TSV')
    parser.add_argument('--species_name', default='Species', help='Species name for heatmap titles')
    parser.add_argument('--cluster_col', default='seurat_clusters', help='Column in adata.obs for cluster IDs')
    parser.add_argument('--cluster_name_col', default='cell_type', help='Column in cluster_names CSV for cell type labels')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load AUCell scores
    with h5py.File(args.auc_loom, 'r') as f:
        auc_struct = f['col_attrs/RegulonsAUC'][:]
        cell_ids   = f['col_attrs/CellID'][:].astype(str)

    regulon_names = list(auc_struct.dtype.names)
    auc_mtx = np.array([auc_struct[n] for n in regulon_names]).T
    auc_df  = pd.DataFrame(auc_mtx, index=cell_ids, columns=regulon_names)
    print(f"AUC matrix: {auc_df.shape}")

    # Fly ortholog map (best bitscore per species gene)
    orthologs = pd.read_csv(args.orthologs, sep='\t')
    orthologs = orthologs[orthologs['pident'] >= 50.0]
    orthologs = orthologs.sort_values('bitscore', ascending=False).drop_duplicates('sp_transcript')
    fly_map   = dict(zip(orthologs['sp_transcript'], orthologs['fly_gene']))

    # Load cluster annotations
    adata = sc.read_h5ad(args.h5ad)
    # Ensure cell ID format consistency (hyphens → underscores if needed)
    adata.obs_names = adata.obs_names.str.replace('-', '_')
    clusters_raw = adata.obs[args.cluster_col]

    # Detect whether clusters are integer IDs (need mapping) or string labels (use directly)
    try:
        clusters_int = clusters_raw.astype(int)
        # Integer cluster IDs — map via cluster_names CSV
        cluster_map_df = pd.read_csv(args.cluster_names)
        cluster_dict = dict(zip(cluster_map_df['cluster'], cluster_map_df[args.cluster_name_col]))
        cluster_labels = clusters_int.map(cluster_dict).fillna(clusters_int.astype(str))
        print(f"Cluster IDs: integer (mapped {len(cluster_dict)} names from CSV)")
    except (ValueError, TypeError):
        # String cluster labels — use directly, CSV is optional
        cluster_labels = clusters_raw.astype(str)
        print(f"Cluster IDs: string labels (using directly, {cluster_labels.nunique()} unique)")

    common_cells = auc_df.index.intersection(cluster_labels.index)
    print(f"Matched cells: {len(common_cells)} / {len(auc_df)}")

    auc_matched = auc_df.loc[common_cells].copy()
    auc_matched['cluster'] = cluster_labels.loc[common_cells]

    cluster_regulons = auc_matched.groupby('cluster').mean()
    print(f"Clusters: {len(cluster_regulons)}, Regulons: {cluster_regulons.shape[1]}")

    cluster_regulons.to_csv(f"{args.output_dir}/cluster_regulon_activity.csv")

    # Display names with duplicate fly name dedup (roman numerals)
    mean_auc_per_reg = cluster_regulons.mean(axis=0).to_dict()
    regulon_display  = make_display_names(regulon_names, fly_map, mean_auc=mean_auc_per_reg)

    with open(f"{args.output_dir}/top_regulons_per_cluster.txt", 'w') as f:
        for cluster in cluster_regulons.index:
            top = cluster_regulons.loc[cluster].nlargest(10)
            f.write(f"\n=== {cluster} ===\n")
            for i, (reg, score) in enumerate(top.items(), 1):
                f.write(f"{i}. {regulon_display.get(reg, reg)}: {score:.4f}\n")

    # Heatmaps
    top_var  = cluster_regulons.var(axis=0).nlargest(50).index
    top_mean = cluster_regulons.mean(axis=0).nlargest(50).index

    save_heatmap(cluster_regulons, top_var, regulon_display,
                 f'Top 50 Variable Regulons Across Clusters — {args.species_name}',
                 f'{args.output_dir}/cluster_regulon_heatmap_variance.png')
    save_heatmap(cluster_regulons, top_mean, regulon_display,
                 f'Top 50 Active Regulons (Mean AUC) Across Clusters — {args.species_name}',
                 f'{args.output_dir}/cluster_regulon_heatmap_mean_auc.png')

    print(f"Outputs in {args.output_dir}/")

if __name__ == '__main__':
    main()
