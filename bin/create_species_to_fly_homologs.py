#!/usr/bin/env python3
"""
Build species_to_fly_homologs.tsv from BLAST output + header/annotation files.
Output format matches what 3_create_motif_annotations.py, 9_analyze_cluster_regulons.py,
10_annotate_regulons.py, and 11_create_regulon_spreadsheet.py expect:
    sp_transcript, sp_gene_name, sp_product, sp_chr, fly_gene, fly_product, fly_chr,
    pident, length, evalue, bitscore
"""
import argparse
import sys
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description='Create species-to-fly ortholog table from BLAST maps')
    parser.add_argument('--blast_map', required=True,
                        help='BLAST outfmt6 file: species protein → fly protein')
    parser.add_argument('--sp_headers', required=True,
                        help='TSV mapping species protein IDs to gene IDs (col1=protein, col2=gene)')
    parser.add_argument('--fly_headers', required=True,
                        help='TSV mapping fly NP_ accessions to gene names (col1=NP_, col2=gene)')
    parser.add_argument('--gff_key', required=True,
                        help='Species GFF key TSV with columns: Gene, rna_type, product, chromosome')
    parser.add_argument('--output', required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    # Load BLAST map (outfmt 6: 12 standard columns)
    blast = pd.read_csv(args.blast_map, sep='\t', header=None,
                        names=['sp_protein', 'fly_np', 'pident', 'length', 'mismatch',
                               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    print(f"BLAST hits: {len(blast)}, unique species proteins: {blast['sp_protein'].nunique()}")

    # Map species protein ID → gene ID
    sp_headers = pd.read_csv(args.sp_headers, sep='\t', header=None, names=['sp_protein', 'sp_transcript'])
    prot_to_gene = dict(zip(sp_headers['sp_protein'], sp_headers['sp_transcript']))
    blast['sp_transcript'] = blast['sp_protein'].map(prot_to_gene)
    unmapped_sp = blast['sp_transcript'].isna().sum()
    if unmapped_sp > 0:
        print(f"WARNING: {unmapped_sp} BLAST hits have no gene mapping ({100*unmapped_sp/len(blast):.1f}%)")
    blast = blast.dropna(subset=['sp_transcript'])
    # Convert hyphens to underscores in gene IDs (Seurat compatibility)
    blast['sp_transcript'] = blast['sp_transcript'].str.replace('-', '_', regex=False)

    # Map fly NP_ → gene name
    fly_headers = pd.read_csv(args.fly_headers, sep='\t', header=None, names=['fly_np', 'fly_gene'])
    np_to_gene = dict(zip(fly_headers['fly_np'], fly_headers['fly_gene']))
    blast['fly_gene'] = blast['fly_np'].map(np_to_gene)
    unmapped_fly = blast['fly_gene'].isna().sum()
    if unmapped_fly > 0:
        print(f"WARNING: {unmapped_fly} fly NP_ accessions unmapped ({100*unmapped_fly/len(blast):.1f}%)")
    blast = blast.dropna(subset=['fly_gene'])

    # Load GFF key for species gene annotation
    gff_key = pd.read_csv(args.gff_key, sep='\t')
    # Handle R-style quoted column names
    gff_key.columns = gff_key.columns.str.strip('"')
    # Detect gene ID column (supports both 'Gene' and 'transcript_id')
    gene_col = 'Gene' if 'Gene' in gff_key.columns else 'transcript_id'
    gff_key[gene_col] = gff_key[gene_col].str.replace('-', '_', regex=False)
    gff_lookup = gff_key.set_index(gene_col)

    blast['sp_gene_name'] = blast['sp_transcript'].map(
        lambda g: gff_lookup.loc[g, 'product'] if g in gff_lookup.index else '')
    blast['sp_product'] = blast['sp_gene_name']  # product = description from GFF
    blast['sp_chr'] = blast['sp_transcript'].map(
        lambda g: str(gff_lookup.loc[g, 'chromosome']) if g in gff_lookup.index else '')
    # Simplify sp_gene_name: use product but truncate if too long
    blast['sp_gene_name'] = blast['sp_gene_name'].fillna('uncharacterized')
    blast['sp_product'] = blast['sp_product'].fillna('')
    blast['sp_chr'] = blast['sp_chr'].fillna('')

    # Fly product/chr left empty (consistent with pharaonis script)
    blast['fly_product'] = ''
    blast['fly_chr'] = ''

    out = blast[['sp_transcript', 'sp_gene_name', 'sp_product', 'sp_chr',
                 'fly_gene', 'fly_product', 'fly_chr',
                 'pident', 'length', 'evalue', 'bitscore']]

    out.to_csv(args.output, sep='\t', index=False)
    print(f"\nWritten: {args.output}")
    print(f"Rows: {len(out)}")
    print(f"Unique species genes: {out['sp_transcript'].nunique()}")
    print(f"Unique fly genes: {out['fly_gene'].nunique()}")
    sys.stdout.flush()


if __name__ == '__main__':
    main()
