#!/usr/bin/env python3
"""
Extract TSS regions from GFF file.
Default: 5kb upstream + 2kb downstream of TSS.
For + strand: [TSS-upstream, TSS+downstream] where TSS = gene start
For - strand: [TSS-downstream, TSS+upstream] where TSS = gene end

One region per gene (deduplicated). Gene names that pandas would silently
convert to NaN are renamed with a '_gene' suffix to prevent cisTarget
AssertionError in update_scores_for_motif_or_track.
"""
# This script was designed for Schistocerca so it might need adjustment to work with other species
import argparse
import re
import sys

# pandas silently converts these to NaN when reading cbust output
PANDAS_NAN_KEYWORDS = {
    'nan', 'na', 'null', 'none', 'n/a', 'n.a.', 'n.a', '#n/a', '#na',
}

def extract_gene_name(attributes):
    """Extract gene name from GFF attributes. Tries Name= first, then gene_id."""
    match = re.search(r'Name=([^;]+)', attributes)
    if match:
        return match.group(1)
    match = re.search(r'gene_id "([^"]+)"', attributes)
    if match:
        return match.group(1)
    match = re.search(r'ID=gene-([^;]+)', attributes)
    if match:
        return match.group(1)
    return None

def process_gff(gff_file, output_bed, upstream, downstream, scaffold_filter=None):
    seen = set()
    written = 0
    nan_renamed = 0

    with open(gff_file) as f, open(output_bed, 'w') as out:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            if scaffold_filter and chrom != scaffold_filter:
                continue

            gene_name = extract_gene_name(attributes)
            if not gene_name:
                continue

            # Guard against pandas NaN keywords
            if gene_name.lower() in PANDAS_NAN_KEYWORDS:
                gene_name = gene_name + '_gene'
                nan_renamed += 1

            # Deduplicate on gene name (one region per gene)
            if gene_name in seen:
                continue
            seen.add(gene_name)

            if strand == '+':
                region_start = max(1, start - upstream)
                region_end = start + downstream
            else:
                region_start = max(1, end - downstream)
                region_end = end + upstream

            out.write(f"{chrom}\t{region_start}\t{region_end}\t{gene_name}\t0\t{strand}\n")
            written += 1

    print(f"Wrote {written} regions to {output_bed} (upstream={upstream}, downstream={downstream})")
    if nan_renamed:
        print(f"  WARNING: Renamed {nan_renamed} gene(s) that match pandas NaN keywords")
    sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract TSS regions from GFF')
    parser.add_argument('gff', help='Input GFF file')
    parser.add_argument('output', help='Output BED file')
    parser.add_argument('--upstream', type=int, default=5000, help='Upstream bp (default: 5000)')
    parser.add_argument('--downstream', type=int, default=2000, help='Downstream bp (default: 2000)')
    parser.add_argument('--scaffold', default=None, help='Filter for specific scaffold')
    args = parser.parse_args()

    process_gff(args.gff, args.output, args.upstream, args.downstream, args.scaffold)
