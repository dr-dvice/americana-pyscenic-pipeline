#!/usr/bin/env python3
"""
Combine partial cisTarget scores databases and generate rankings.
Replicates the ranking algorithm from cistarget_db.py (seeded RNG tiebreaker).

Usage: python combine_partial_dbs.py <output_prefix> <total_parts> [seed]
Example: python combine_partial_dbs.py outputs/full_db/sa_v1 103 0
"""
import sys
import glob
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as pf


def load_and_combine_partials(output_prefix, total_parts):
    parts = []
    for i in range(1, total_parts + 1):
        path = f"{output_prefix}.part_{i:04d}_of_{total_parts:04d}.motifs_vs_regions.scores.feather"
        df = pf.read_table(path).to_pandas()
        df = df.set_index("regions")
        parts.append(df)
        print(f"  Loaded part {i}/{total_parts}: {df.shape[1]} motifs")

    combined = pd.concat(parts, axis=1)
    print(f"Combined scores: {combined.shape[0]} regions x {combined.shape[1]} motifs")
    return combined


def create_rankings(scores_mvr, seed):
    """Replicate cistarget_db.py ranking: descending score, seeded random tiebreak."""
    rng = np.random.default_rng(seed=seed)
    n_motifs, n_regions = scores_mvr.shape
    rankings = np.empty(scores_mvr.shape, dtype=np.int32)

    for i in range(n_motifs):
        scores = scores_mvr[i, :]
        perm = rng.permutation(n_regions)
        rank_arr = np.empty(n_regions, dtype=np.int32)
        rank_arr[perm[(-scores)[perm].argsort()]] = np.arange(n_regions, dtype=np.int32)
        rankings[i, :] = rank_arr

    return rankings


def write_feather(df, path):
    # ctxcore needs the identifier column ("motifs"/"regions") as a real column
    # (not a pandas index) AND as the LAST column: it scans from the last column
    # and its guard `not index_column_idx` wrongly rejects idx 0. So reset the
    # index to a column, move it to the end, and drop pandas index metadata.
    idx_name = df.index.name
    out = df.reset_index()
    out = out[[c for c in out.columns if c != idx_name] + [idx_name]]
    table = pa.Table.from_pandas(out, preserve_index=False)
    pf.write_feather(table, path, version=2)
    print(f"  Wrote: {path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python combine_partial_dbs.py <output_prefix> <total_parts> [seed]")
        sys.exit(1)

    output_prefix = sys.argv[1]
    total_parts = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    # 1. Load and combine partial scores
    print("Loading partial databases...")
    combined = load_and_combine_partials(output_prefix, total_parts)

    # 2. Write combined scores: regions_vs_motifs (regions as rows, motifs as columns)
    print("Writing combined scores (regions vs motifs)...")
    combined.index.name = "regions"
    write_feather(combined, f"{output_prefix}.regions_vs_motifs.scores.feather")

    # 3. Write transposed scores: motifs_vs_regions
    print("Writing combined scores (motifs vs regions)...")
    transposed = combined.T
    transposed.index.name = "motifs"
    write_feather(transposed, f"{output_prefix}.motifs_vs_regions.scores.feather")

    # 4. Create rankings (motifs as rows, rank each motif's scores across regions)
    print(f"Creating rankings (seed={seed})...")
    rankings_np = create_rankings(transposed.to_numpy(), seed)
    rankings_df = pd.DataFrame(rankings_np, index=transposed.index, columns=transposed.columns)
    rankings_df.index.name = "motifs"

    # 5. Write rankings: regions_vs_motifs (index=motifs, columns=regions).
    # ctxcore parses "regions_vs_motifs" as column_kind=regions, row_kind=motifs,
    # so the index MUST be named "motifs"
    print("Writing rankings (regions vs motifs)...")
    write_feather(rankings_df, f"{output_prefix}.regions_vs_motifs.rankings.feather")

    print(f"\nDone. Rankings database: {output_prefix}.regions_vs_motifs.rankings.feather")


if __name__ == "__main__":
    main()
