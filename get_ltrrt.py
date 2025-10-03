#!/usr/bin/env python3
import sys, gzip, argparse

def smart_open(path, mode="rt"):
    if path == "-" or path is None:
        return sys.stdin if "r" in mode else sys.stdout
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def load_fasta_into_memory(fasta_path):
    """Return dict: chrom -> sequence."""
    seqs = {}
    name, out = None, []
    with smart_open(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(out).upper()
                name = line[1:].strip().split()[0]
                out = []
            else:
                out.append(line.strip())
        if name:
            seqs[name] = "".join(out).upper()
    return seqs

def slice_seq(seq, start0, end0):
    start0 = max(0, start0)
    end0 = min(len(seq), end0)
    return seq[start0:end0].upper()

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract LTR-RT sequences to FASTA from a table and a genome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("-t","--table", required=True, help="Input table (use '-' for stdin)")
    p.add_argument("-g","--genome", required=True, help="Genome FASTA")
    p.add_argument("-o","--out", default="-", help="Output FASTA (use '-' for stdout)")
    p.add_argument("--chr-col", type=int, default=12, help="1-based column index for chromosome")
    p.add_argument("--start-col", type=int, default=1, help="1-based column index for s(ret)")
    p.add_argument("--end-col", type=int, default=2, help="1-based column index for e(ret)")
    p.add_argument("--delimiter", default=None, help="Field delimiter (default: whitespace, use '\\t' for tabs)")
    p.add_argument("--zero-based", action="store_true", help="Coordinates are 0-based half-open (default: 1-based inclusive)")
    return p.parse_args()

def main():
    args = parse_args()
    delim = None if args.delimiter in (None, "", "\\s+") else ("\t" if args.delimiter == "\\t" else args.delimiter)

    # Load genome into memory
    seqs = load_fasta_into_memory(args.genome)
    if not seqs:
        sys.exit("[error] No sequences loaded from FASTA.")

    out = smart_open(args.out, "wt")

    with smart_open(args.table, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or line == "###":
                continue
            fields = line.split(delim) if delim else line.split()
            try:
                chrom = fields[args.chr_col-1]
                s = int(float(fields[args.start-col-1])) if False else int(float(fields[args.start_col-1]))
                e = int(float(fields[args.end_col-1]))
            except Exception as ex:
                print(f"[warn] skipping line (parse error): {line}\n  -> {ex}", file=sys.stderr)
                continue

            if args.zero_based:
                start0, end0 = min(s,e), max(s,e)
                start_print, end_print = start0, end0
            else:
                start0, end0 = min(s,e)-1, max(s,e)  # 1-based inclusive -> 0-based half-open
                start_print, end_print = min(s,e), max(s,e)

            if chrom not in seqs:
                print(f"[warn] chromosome not in FASTA: {chrom}", file=sys.stderr)
                continue

            seq = slice_seq(seqs[chrom], start0, end0)
            if not seq:
                print(f"[warn] empty seq {chrom}:{start_print}-{end_print}", file=sys.stderr)
                continue

            # Header format: >Chr_start_end
            header = f">{chrom}_{start_print}_{end_print}"
            out.write(header + "\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()
