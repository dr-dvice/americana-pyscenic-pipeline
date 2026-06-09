#!/usr/bin/env python3
"""
Strip ::chrom:start-end(strand) suffixes that bedtools getfasta -name appends.
Validates that no :: remain after cleaning.
"""
import sys

if len(sys.argv) != 3:
    sys.exit("Usage: 2_clean_fasta_headers.py <input.fa> <output.fa>")

input_fa, output_fa = sys.argv[1], sys.argv[2]
cleaned = 0

with open(input_fa) as fin, open(output_fa, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            line = line.split('::')[0] + '\n'
            cleaned += 1
        fout.write(line)

# Validate
with open(output_fa) as f:
    for line in f:
        if line.startswith('>') and '::' in line:
            sys.exit(f"ERROR: :: suffix still present: {line.strip()}")

print(f"Cleaned {cleaned} headers → {output_fa}")
sys.stdout.flush()
