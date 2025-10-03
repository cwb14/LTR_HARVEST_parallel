#!/usr/bin/env python3
"""
cluster_ltr.py

Cluster LTR-RT entries by overlap (transitive). Designed for LTR_parallel output:
  Columns (0-based):
    3: s(lLTR)
    7: e(rLTR)
   11: chr

- Ignores lines starting with '#'
- Outputs original lines, with clusters separated by a line '###'
- Lines remain in their original input order inside each cluster.

Usage:
    python cluster_ltr.py input.tsv > output.tsv
    or
    cat input.tsv | python cluster_ltr.py -   # use '-' to read stdin
"""

import sys
import argparse

# Column indices for this specific input format
CHR_COL = 11      # 'chr'
LLTR_START_COL = 3  # s(lLTR)
RLTR_END_COL = 7    # e(rLTR)

def parse_args():
    p = argparse.ArgumentParser(description="Cluster LTR-RTs by transitive overlap (LTR_parallel format).")
    p.add_argument("infile", help="Input TSV file (use '-' for stdin)")
    return p.parse_args()

def read_records(f):
    """
    Returns list of tuples: (orig_index, chrom, start, end, line)
    - Skips empty lines
    - Skips comment lines starting with '#'
    - Skips lines with non-integer coords or not enough columns
    """
    recs = []
    for i, raw in enumerate(f):
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        # Need at least up to the columns we access
        max_idx = max(CHR_COL, LLTR_START_COL, RLTR_END_COL)
        if len(parts) <= max_idx:
            continue
        chrom = parts[CHR_COL]
        try:
            start = int(parts[LLTR_START_COL])
            end   = int(parts[RLTR_END_COL])
        except ValueError:
            # skip lines where coordinates aren't integers
            continue
        recs.append((i, chrom, start, end, line))
    return recs

def cluster_by_chrom(records):
    """
    records: list of (orig_index, chrom, start, end, line)
    Returns list of clusters:
      [ [ (orig_index, chrom, start, end, line), ... ], ... ]
    Clusters ordered by the first appearance (lowest orig_index) among their members.
    """
    from collections import defaultdict

    chrom_map = defaultdict(list)
    for rec in records:
        chrom_map[rec[1]].append(rec)

    all_clusters = []

    for chrom, recs in chrom_map.items():
        # sort by start for merging to find overlapping/transitive groups
        recs_sorted = sorted(recs, key=lambda r: (r[2], r[3]))
        i = 0
        n = len(recs_sorted)
        while i < n:
            # start a new cluster with recs_sorted[i]
            cluster_members = [recs_sorted[i]]
            cur_end = recs_sorted[i][3]
            i += 1
            # extend cluster by any recs whose start <= cur_end (transitive merging)
            while i < n and recs_sorted[i][2] <= cur_end:
                cluster_members.append(recs_sorted[i])
                if recs_sorted[i][3] > cur_end:
                    cur_end = recs_sorted[i][3]
                i += 1

            # restore original input order inside the cluster
            cluster_members_sorted_by_input = sorted(cluster_members, key=lambda r: r[0])
            all_clusters.append(cluster_members_sorted_by_input)

    # Sort clusters across chromosomes by first original index for deterministic output
    all_clusters.sort(key=lambda cl: cl[0][0] if cl else float('inf'))
    return all_clusters

def main():
    args = parse_args()

    if args.infile == "-":
        f = sys.stdin
    else:
        try:
            f = open(args.infile, "r")
        except Exception as e:
            sys.stderr.write(f"Error opening {args.infile}: {e}\n")
            sys.exit(1)

    records = read_records(f)
    if f is not sys.stdin:
        f.close()

    if not records:
        return

    clusters = cluster_by_chrom(records)

    out_lines = []
    for ci, cluster in enumerate(clusters):
        for rec in cluster:
            out_lines.append(rec[4])
        if ci != len(clusters) - 1:
            out_lines.append("###")

    sys.stdout.write("\n".join(out_lines) + ("\n" if out_lines else ""))

if __name__ == "__main__":
    main()
