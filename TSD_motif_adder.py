#!/usr/bin/env python3
import argparse
from collections import namedtuple
import sys
import gzip

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def load_fasta(path):
    """Load FASTA into {seq_id: UPPERCASE sequence}."""
    seqs = {}
    with open_maybe_gzip(path, "rt") as fh:
        name, chunks = None, []
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    if not seqs:
        sys.exit(f"[ERROR] No sequences loaded from {path}")
    return seqs

def safe_slice(seq, start0, end0_excl):
    """Return seq[start0:end0_excl] if within bounds, else ''. """
    if start0 < 0 or end0_excl > len(seq) or start0 >= end0_excl:
        return ""
    return seq[start0:end0_excl]

def parse_scn_line(line):
    """
    Expect 12 columns like:
    s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
    """
    parts = line.strip().split()
    if len(parts) < 12:
        raise ValueError("SCN line with <12 columns")
    vals = {
        "s_ret": int(parts[0]),
        "e_ret": int(parts[1]),
        "s_lLTR": int(parts[3]),
        "e_lLTR": int(parts[4]),
        "s_rLTR": int(parts[6]),
        "e_rLTR": int(parts[7]),
        "chr": parts[11],
        "rest": parts,  # keep original tokens for passthrough
    }
    return vals

def find_tsd(chrom_seq, s_ret, e_ret, vic, tsd_len):
    """
    Search ±vic around both element bounds for a perfect TSD of length tsd_len.
    Positions are 1-based inclusive in the SCN; convert to 0-based here.
    Returns tsd_string or 'NA'.
    """
    best = None
    for left_off in range(-vic, vic + 1):
        left_start0 = (s_ret + left_off - tsd_len) - 1
        left_end0 = (s_ret + left_off) - 1
        left = safe_slice(chrom_seq, left_start0, left_end0)
        if len(left) != tsd_len:
            continue
        for right_off in range(-vic, vic + 1):
            right_start0 = e_ret + right_off  # (e_ret + 1 + ro) - 1
            right_end0 = right_start0 + tsd_len
            right = safe_slice(chrom_seq, right_start0, right_end0)
            if len(right) != tsd_len:
                continue
            if left == right:
                score = abs(left_off) + abs(right_off)
                cand = (score, left)
                if best is None or cand < best:
                    best = cand
    return best[1] if best else "NA"

def get_motif(chrom_seq, s_lLTR, e_rLTR):
    """
    Report dinucleotides at LTR boundaries:
      start of left LTR (expected 'TG'), end of right LTR (expected 'CA').
    Returns 'XX-YY' or 'NA'.
    """
    left_start0 = s_lLTR - 1
    left = safe_slice(chrom_seq, left_start0, left_start0 + 2)
    right_start0 = e_rLTR - 2
    right = safe_slice(chrom_seq, right_start0, right_start0 + 2)
    if len(left) == 2 and len(right) == 2:
        return f"{left}-{right}"
    return "NA"

RunInfo = namedtuple("RunInfo", "kind start0 end0 length seq")

def longest_poly_run(seq, allowed_set):
    """
    Return (start0, end0_excl, length, seq) of the longest contiguous run of bases
    fully in allowed_set within 'seq'. Ties broken by leftmost (smallest start).
    If none, return None.
    """
    best = None
    i, n = 0, len(seq)
    while i < n:
        if seq[i] in allowed_set:
            j = i + 1
            while j < n and seq[j] in allowed_set:
                j += 1
            length = j - i
            cand = (length, i, j, seq[i:j])
            if best is None or cand > (best[0], -best[1], best[2], best[3]):  # prefer longer, then leftmost
                best = (length, i, j, seq[i:j])
            i = j
        else:
            i += 1
    if best is None:
        return None
    length, i, j, s = best
    return (i, j, length, s)

def find_ppt(chrom, chrom_seq, s_lLTR, e_lLTR, s_rLTR, e_rLTR, ppt_min, ppt_window):
    """
    Strand-agnostic PPT finder:
      - Window 1: just right of left LTR end  -> [e_lLTR .. e_lLTR + ppt_window - 1]
      - Window 2: just left  of right LTR start -> [s_rLTR - ppt_window .. s_rLTR - 1]
    In each window, find the longest run of purines (A/G) and pyrimidines (C/T).
    Choose the single best (longest). Require length >= ppt_min.
    Return pretty string or 'NA'.
    """
    pur = set("AG")
    pyr = set("CT")

    best_overall = None  # (length, chr_start1, chr_end1, type, seq)

    # Window 1: right of left LTR
    win1_start0 = e_lLTR  # e_lLTR is 1-based inclusive; +1 -> index e_lLTR, but slice is [e_lLTR, e_lLTR+ppt_window)
    win1_end0 = e_lLTR + ppt_window
    w1 = safe_slice(chrom_seq, win1_start0, win1_end0)
    if w1:
        for kind, allowed in (("PUR", pur), ("PYR", pyr)):
            res = longest_poly_run(w1, allowed)
            if res:
                i0, j0, L, s = res
                chr_start1 = win1_start0 + i0 + 1  # back to 1-based inclusive
                chr_end1 = win1_start0 + j0        # inclusive
                cand = (L, chr_start1, chr_end1, kind, s)
                if best_overall is None or cand > best_overall:
                    best_overall = cand

    # Window 2: left of right LTR
    win2_start0 = max(0, s_rLTR - 1 - ppt_window)  # start index inclusive
    win2_end0 = s_rLTR - 1                          # end index exclusive
    w2 = safe_slice(chrom_seq, win2_start0, win2_end0)
    if w2:
        for kind, allowed in (("PUR", pur), ("PYR", pyr)):
            res = longest_poly_run(w2, allowed)
            if res:
                i0, j0, L, s = res
                chr_start1 = win2_start0 + i0 + 1  # 1-based
                chr_end1 = win2_start0 + j0        # inclusive
                cand = (L, chr_start1, chr_end1, kind, s)
                if best_overall is None or cand > best_overall:
                    best_overall = cand

    if (best_overall is None) or (best_overall[0] < ppt_min):
        return "NA"

    L, chr_start1, chr_end1, kind, s = best_overall
    return f"{kind};len={L};loc={chrom}:{chr_start1}-{chr_end1};seq={s}"

def main():
    ap = argparse.ArgumentParser(
        description="Add TSD, Motif, and PPT columns to LTRharvest/FINDER .scn output (strand-agnostic PPT)."
    )
    ap.add_argument("--genome", required=True, help="Genome FASTA (optionally .gz)")
    ap.add_argument("--scn", required=True, help="Input .scn file (optionally .gz)")
    ap.add_argument("-vic", type=int, default=2, help="Vicinity (±bp) to search for TSD (default: 2)")
    ap.add_argument("--tsd-len", type=int, default=5, help="TSD length to search (default: 5)")
    ap.add_argument("--ppt-min", type=int, default=8, help="Minimum poly-tract length to report PPT (default: 8)")
    ap.add_argument("--ppt-window", type=int, default=50, help="Window (bp) flanking each LTR junction to search PPT (default: 50)")
    args = ap.parse_args()

    seqs = load_fasta(args.genome)

    header_augmented = False
    with open_maybe_gzip(args.scn, "rt") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                print(line)
                continue
            if line.startswith("#"):
                print(line)
                # Insert a one-time note about appended columns right after the "predictions are..." header block
                if (not header_augmented) and "predictions are reported" in line:
                    print("# Added columns: TSD Motif PPT")
                    header_augmented = True
                continue

            # Data line
            try:
                rec = parse_scn_line(line)
            except Exception:
                # Pass-through unparseable lines
                print(line)
                continue

            chrom = rec["chr"]
            if chrom not in seqs:
                # Try simple Chr/chr toggle
                alt = ("chr" + chrom[3:]) if chrom.startswith("Chr") else ("Chr" + chrom[3:] if chrom.startswith("chr") else chrom)
                if alt in seqs:
                    chrom = alt
                else:
                    tsd = "NA"
                    motif = "NA"
                    ppt = "NA"
                    print("\t".join(map(str, rec["rest"] + [tsd, motif, ppt])))
                    continue

            chrom_seq = seqs[chrom]

            tsd = find_tsd(
                chrom_seq,
                rec["s_ret"],
                rec["e_ret"],
                args.vic,
                args.tsd_len
            )
            motif = get_motif(chrom_seq, rec["s_lLTR"], rec["e_rLTR"])
            ppt = find_ppt(
                chrom, chrom_seq,
                rec["s_lLTR"], rec["e_lLTR"],
                rec["s_rLTR"], rec["e_rLTR"],
                args.ppt_min, args.ppt_window
            )

            print("\t".join(map(str, rec["rest"] + [tsd, motif, ppt])))

if __name__ == "__main__":
    main()
