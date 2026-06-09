#!/usr/bin/env python3
"""
Convert fly motif annotations to species-specific using ortholog mapping.
Replaces fly gene names with species orthologs in the motif annotation table.
"""
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create species-specific motif annotations from fly')
    parser.add_argument('--fly_tbl', required=True, help='Fly v10nr motif annotation TBL file')
    parser.add_argument('--orthologs', required=True, help='Species-to-fly ortholog TSV')
    parser.add_argument('--output', required=True, help='Output motif annotation TBL')
    parser.add_argument('--min_pident', type=float, default=50.0, help='Minimum percent identity (default: 50)')
    args = parser.parse_args()

    print("Loading ortholog mapping...")
    orthologs = pd.read_csv(args.orthologs, sep='\t')
    print(f"  Total orthologs: {len(orthologs)}")

    orthologs_filtered = orthologs[orthologs['pident'] >= args.min_pident].copy()
    print(f"  After pident >= {args.min_pident}%: {len(orthologs_filtered)}")

    # Build fly_gene -> list of species genes mapping
    fly_to_species = {}
    for _, row in orthologs_filtered.iterrows():
        fly_gene = row['fly_gene']
        species_gene = row['sp_transcript']
        if fly_gene not in fly_to_species:
            fly_to_species[fly_gene] = []
        fly_to_species[fly_gene].append(species_gene)

    print(f"\nFly genes with species orthologs: {len(fly_to_species)}")

    print("\nLoading fly motif annotations...")
    fly_annot = pd.read_csv(args.fly_tbl, sep='\t')
    has_hash = fly_annot.columns[0].startswith('#')
    fly_annot.columns = fly_annot.columns.str.replace('^#', '', regex=True)
    print(f"  Total fly annotation rows: {len(fly_annot)}")

    print("\nConverting annotations...")
    species_rows = []
    matched_genes = set()
    unmatched_genes = set()

    for _, row in fly_annot.iterrows():
        fly_gene = row['gene_name']
        if fly_gene in fly_to_species:
            matched_genes.add(fly_gene)
            for species_gene in fly_to_species[fly_gene]:
                new_row = row.copy()
                new_row['gene_name'] = species_gene
                species_rows.append(new_row)
        else:
            unmatched_genes.add(fly_gene)

    species_annot = pd.DataFrame(species_rows)

    print(f"\nResults:")
    print(f"  Fly genes matched: {len(matched_genes)}")
    print(f"  Fly genes unmatched: {len(unmatched_genes)}")
    print(f"  Original fly annotation rows: {len(fly_annot)}")
    print(f"  New species annotation rows: {len(species_annot)}")
    print(f"  Unique species genes in annotations: {species_annot['gene_name'].nunique()}")

    if has_hash:
        species_annot.columns = ['#' + species_annot.columns[0]] + list(species_annot.columns[1:])
    species_annot.to_csv(args.output, sep='\t', index=False)
    print(f"\nSaved to {args.output}")

if __name__ == '__main__':
    main()
