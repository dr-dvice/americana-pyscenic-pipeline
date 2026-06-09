#!/usr/bin/env python3
"""
Annotate regulon IDs with gene names from GFF and ortholog mapping.
"""
import argparse
import os
import pandas as pd
import re

def parse_gff_annotations(gff_file):
    gene_annotations = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            attributes = fields[8]
            gene_id_match = re.search(r'ID=gene-([^;]+)', attributes)
            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)
            name_match = re.search(r'Name=([^;]+)', attributes)
            gene_name = name_match.group(1) if name_match else None
            desc_match = re.search(r'description=([^;]+)', attributes)
            description = desc_match.group(1) if desc_match else None
            gene_annotations[gene_id] = {
                'name': gene_name if gene_name else gene_id,
                'description': description if description else 'uncharacterized'
            }
    return gene_annotations

def annotate_regulon(regulon_id, gene_annotations, ortholog_map):
    gene_id = regulon_id.replace('(+)', '').replace('(-)', '')
    gff_info = gene_annotations.get(gene_id, {'name': gene_id, 'description': 'not found'})
    gene_name = gff_info['name']
    description = gff_info['description']
    fly_gene = ortholog_map.get(gene_id, 'no ortholog')

    if gene_name != gene_id:
        annotation = f"{gene_name} ({fly_gene})"
    else:
        annotation = f"{gene_id} ({fly_gene})"

    return {
        'regulon_id': regulon_id,
        'gene_id': gene_id,
        'gene_name': gene_name,
        'fly_ortholog': fly_gene,
        'description': description,
        'annotation': annotation
    }

def main():
    parser = argparse.ArgumentParser(description='Annotate regulons with gene names and orthologs')
    parser.add_argument('--gff', required=True, help='GFF file for gene annotations')
    parser.add_argument('--orthologs', required=True, help='Species-to-fly ortholog TSV')
    parser.add_argument('--cluster_activity', required=True, help='Cluster regulon activity CSV')
    parser.add_argument('--top_regulons', required=True, help='Top regulons per cluster TXT')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Parsing GFF for gene names...")
    gene_annotations = parse_gff_annotations(args.gff)
    print(f"  Found {len(gene_annotations)} gene annotations")

    print("Loading ortholog mapping...")
    orthologs = pd.read_csv(args.orthologs, sep='\t')
    orthologs = orthologs.sort_values('bitscore', ascending=False).drop_duplicates('sp_transcript')
    ortholog_map = dict(zip(orthologs['sp_transcript'], orthologs['fly_gene']))
    print(f"  Found {len(ortholog_map)} ortholog mappings")

    print("Annotating cluster regulon activity matrix...")
    cluster_activity = pd.read_csv(args.cluster_activity, index_col=0)

    regulon_annotations = pd.DataFrame([
        annotate_regulon(reg, gene_annotations, ortholog_map) for reg in cluster_activity.columns
    ])
    regulon_annotations.set_index('regulon_id', inplace=True)

    regulon_annotations.to_csv(f"{args.output_dir}/regulon_annotations.csv")
    print(f"Saved: {args.output_dir}/regulon_annotations.csv")

    cluster_activity_annotated = cluster_activity.copy()
    cluster_activity_annotated.columns = [regulon_annotations.loc[col, 'annotation'] for col in cluster_activity.columns]
    cluster_activity_annotated.to_csv(f"{args.output_dir}/cluster_regulon_activity_annotated.csv")
    print(f"Saved: {args.output_dir}/cluster_regulon_activity_annotated.csv")

    # Annotate top regulons — generalized regex for any gene ID format
    print("Annotating top regulons per cluster...")
    with open(args.top_regulons, 'r') as f_in, \
         open(f"{args.output_dir}/top_regulons_per_cluster_annotated.txt", 'w') as f_out:
        for line in f_in:
            line = line.rstrip()
            # Match any non-whitespace string followed by (+) or (-)
            match = re.search(r'(\S+\([+-]\))', line)
            if match:
                regulon_id = match.group(1)
                if regulon_id in regulon_annotations.index:
                    annot = regulon_annotations.loc[regulon_id]
                    annotated_line = line.replace(
                        regulon_id,
                        f"{regulon_id} [{annot['gene_name']} / {annot['fly_ortholog']}]"
                    )
                    f_out.write(annotated_line + '\n')
                else:
                    f_out.write(line + '\n')
            else:
                f_out.write(line + '\n')

    print(f"Saved: {args.output_dir}/top_regulons_per_cluster_annotated.txt")

    # Top 20 overall
    print("Creating summary table...")
    mean_activity = cluster_activity.mean(axis=0).sort_values(ascending=False)
    top_overall = mean_activity.head(20)

    summary_df = pd.DataFrame({
        'regulon_id': top_overall.index,
        'mean_AUC': top_overall.values
    })
    summary_df = summary_df.merge(regulon_annotations, left_on='regulon_id', right_index=True)
    summary_df = summary_df[['regulon_id', 'gene_name', 'fly_ortholog', 'description', 'mean_AUC']]
    summary_df.to_csv(f"{args.output_dir}/top20_regulons_overall.csv", index=False)
    print(f"Saved: {args.output_dir}/top20_regulons_overall.csv")

    print("\nTop 20 regulons overall:")
    print(summary_df.to_string(index=False))

if __name__ == '__main__':
    main()
