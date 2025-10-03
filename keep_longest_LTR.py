#!/usr/bin/env python3
import sys

def avg_ltr_len(fields):
    # fields are 0-based; l(lLTR)=col 6 -> idx 5, l(rLTR)=col 9 -> idx 8
    return (float(fields[5]) + float(fields[8])) / 2.0

def flush_cluster(best_line, best_avg):
    if best_line is not None:
        sys.stdout.write(f"{best_line}\t{best_avg:.6f}\n")

def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else "-"
    fh = sys.stdin if infile == "-" else open(infile, "r")

    best_line = None
    best_avg = None

    for raw in fh:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            # Skip blank/comment lines (except the cluster marker handled below)
            if line.strip() == "###":
                # Treat '###' as a cluster boundary too (in case comment-char present)
                pass
            else:
                continue

        if line.strip() == "###":
            # End of a cluster → output winner and reset
            flush_cluster(best_line, best_avg)
            best_line, best_avg = None, None
            continue

        # Data line: whitespace/tab-delimited
        fields = line.split()
        # Compute average of l(lLTR) and l(rLTR)
        avg = avg_ltr_len(fields)

        if best_avg is None or avg > best_avg:
            best_line, best_avg = line, avg

    # Flush the last cluster if file doesn't end with ###
    flush_cluster(best_line, best_avg)

    if fh is not sys.stdin:
        fh.close()

if __name__ == "__main__":
    main()
