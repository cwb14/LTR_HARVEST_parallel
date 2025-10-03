#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TEsorter.py — extract internal regions from LTR-retriever scn predictions, run TEsorter,
and append TEsorter calls back onto the SCN file.

Usage (example):
  python TEsorter.py --genome 9311v2.fa --scn ltr_harvest_finder_mod.scn \
    --TEsorter-args "-db rexdb-plant -p 200 -cov 10 -eval 1e-2 -rule 70-30-80"

Outputs:
  - <outdir>/LTR_internal.fa
  - <outdir>/LTR_internal.fa.<db>.cls.tsv                      (TEsorter result)
  - <scn>.with_tesorter.tsv                                    (SCN + 6 cols)

Appended columns (always in this order):
  Class  Order  Superfamily  Confident  Strand  Domains

SCN columns assumed (tab-delimited):
  1 s(ret)  2 e(ret)  3 l(ret)  4 s(lLTR)  5 e(lLTR)  6 l(lLTR)
  7 s(rLTR) 8 e(rLTR) 9 l(rLTR) 10 sim    11 seq-nr  12 chr [13.. extras]

Internal region = (e(lLTR)+1) .. (s(rLTR)-1) by default (LTRs excluded).
Use --include-ltr-boundaries to switch to [e(lLTR), s(rLTR)] inclusive.
"""
import argparse
import os
import sys
import subprocess
import shlex
import re

############################################################
# FASTA / SCN helpers
############################################################

def read_fasta_to_dict(path):
    seqs = {}
    name = None
    chunks = []
    with open(path, 'r') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks).replace('\n', '').replace('\r', '')
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = ''.join(chunks).replace('\n', '').replace('\r', '')
    for k in list(seqs.keys()):
        seqs[k] = seqs[k].upper()
    return seqs

def parse_scn(path):
    """
    Yields tuples: (lineno, parts(list[str]), parsed(dict))
      parsed contains keys used here; original parts preserved.
    Skips comments/blank lines but returns them for passthrough if needed.
    """
    with open(path, 'r') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line:
                yield (lineno, None, None, line)  # blank passthrough
                continue
            if line.startswith('#'):
                yield (lineno, None, None, line)  # comment passthrough
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                # Malformed; keep as passthrough data row but cannot parse
                yield (lineno, parts, None, None)
                continue
            try:
                parsed = {
                    "s_ret": int(parts[0]),
                    "e_ret": int(parts[1]),
                    "s_lLTR": int(parts[3]),
                    "e_lLTR": int(parts[4]),
                    "s_rLTR": int(parts[6]),
                    "e_rLTR": int(parts[7]),
                    "sim": parts[9],
                    "seq_nr": parts[10],  # not always numeric
                    "chr": parts[11],
                }
            except Exception:
                parsed = None
            yield (lineno, parts, parsed, None)

def write_internal_fasta(genome, scn, outfa, include_ltr_boundaries=False):
    """
    Extract internal regions into FASTA.
    Returns (count_written, skipped_list, key_index)
      - key_index: maps SCN row index (1..N over data rows) -> FASTA header (original) and normalized ID
    """
    seqs = read_fasta_to_dict(genome)
    n_ok, skipped = 0, []
    key_index = {}  # row_idx -> {"hdr": original_header, "norm": normalized_id}

    with open(outfa, 'w') as out:
        data_row_idx = 0
        for lineno, parts, parsed, passthru in parse_scn(scn):
            if passthru is not None or parts is None:
                continue  # comments/blank
            data_row_idx += 1
            if parsed is None:
                skipped.append((data_row_idx, "unparsed_row", f"line {lineno}"))
                continue

            chrom = parsed["chr"]
            if chrom not in seqs:
                skipped.append((data_row_idx, "missing_chr", chrom))
                continue
            e_l = parsed["e_lLTR"]
            s_r = parsed["s_rLTR"]

            # Internal coordinates (1-based inclusive)
            if include_ltr_boundaries:
                start_1 = e_l
                end_1 = s_r
            else:
                start_1 = e_l + 1
                end_1 = s_r - 1
            if end_1 < start_1:
                skipped.append((data_row_idx, "empty_internal", f"{start_1}-{end_1}"))
                continue

            # Convert to 0-based slicing
            start0 = start_1 - 1
            end0_excl = end_1
            chrom_seq = seqs[chrom]
            if end0_excl > len(chrom_seq) or start0 < 0:
                skipped.append((data_row_idx, "out_of_bounds", f"{chrom}:{start_1}-{end_1}"))
                continue

            internal = chrom_seq[start0:end0_excl]

            # FASTA header (original with pipes)
            header = (
                f">{chrom}:{start_1}-{end_1}"
                f"|ret:{parsed['s_ret']}-{parsed['e_ret']}"
                f"|sim:{parsed['sim']}|idx:{data_row_idx}"
            )
            out.write(header + "\n")
            for j in range(0, len(internal), 60):
                out.write(internal[j:j+60] + "\n")
            n_ok += 1

            # Normalized ID (TEsorter often replaces '|' with '_' in the first field)
            norm = header[1:].replace('|', '_')
            key_index[data_row_idx] = {"hdr": header[1:], "norm": norm}

    return n_ok, skipped, key_index

############################################################
# TEsorter helpers
############################################################

def ensure_tesorter(clone_dir="TEsorter"):
    if os.path.isdir(clone_dir) and os.path.exists(os.path.join(clone_dir, "__init__.py")):
        return True
    if os.path.isdir(clone_dir):
        return True
    cmd = ["git", "clone", "https://github.com/cwb14/TEsorter.git", clone_dir]
    try:
        subprocess.check_call(cmd)
        return True
    except FileNotFoundError:
        sys.stderr.write("[error] 'git' not found in PATH. Please install git.\n")
        return False
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[error] git clone failed: {e}\n")
        return False

def detect_db_from_args(tesorter_args):
    """
    Parse -db <name> (or --db <name>) from arg string; default 'rexdb-plant'
    """
    if not tesorter_args:
        return "rexdb-plant"
    toks = shlex.split(tesorter_args)
    for i, t in enumerate(toks):
        if t in ("-db", "--db") and i + 1 < len(toks):
            return toks[i + 1]
        # also allow combined form: -db=rexdb-plant / --db=rexdb-plant
        if t.startswith("-db=") or t.startswith("--db="):
            return t.split("=", 1)[1]
    return "rexdb-plant"

def run_tesorter(internal_fa, tesorter_args, tesorter_dir="TEsorter", workdir="."):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.abspath(tesorter_dir) + os.pathsep + env.get("PYTHONPATH", "")
    cmd = ["python3", "-m", "TEsorter", internal_fa]
    if tesorter_args:
        cmd.extend(shlex.split(tesorter_args))
    sys.stderr.write("[info] Running: " + " ".join(shlex.quote(x) for x in cmd) + "\n")
    return subprocess.call(cmd, cwd=workdir, env=env)

############################################################
# Post-processing: join SCN with TEsorter calls
############################################################

TE_COLS = ["Class", "Order", "Superfamily", "Confident", "Strand", "Domains"]

def read_tesorter_cls(cls_path):
    """
    Reads <FA>.db.cls.tsv into dict: id -> [Class, Order, Superfamily, Confident, Strand, Domains]
    - Accepts any extra trailing columns; only first 7 used (id + 6)
    - Returns both original IDs and an underscore-normalized variant as keys.
    """
    ann = {}
    if not os.path.exists(cls_path):
        sys.stderr.write(f"[warn] TEsorter result not found: {cls_path}\n")
        return ann
    with open(cls_path, 'r') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rec_id = parts[0]
            vals = parts[1:7]
            # If fewer than 6 columns, pad with NA
            if len(vals) < 6:
                vals = vals + ["NA"] * (6 - len(vals))
            # Store under exact id
            ann[rec_id] = vals
            # Also store a "normalized" version (in case '|' became '_' or vice versa)
            # We try a few typical variants
            variants = set()
            variants.add(rec_id.replace('|', '_'))
            variants.add(rec_id.replace('_', '|'))
            # Some environments also collapse multiple separators; keep simple here
            for v in variants:
                ann.setdefault(v, vals)
    return ann

def append_tesorter_to_scn(scn_in, scn_out, key_index, tesorter_ann):
    """
    Writes scn_out by appending 6 TEsorter columns to each data row.
    For unmatched rows, fills 'NA' in all 6 columns.
    Preserves SCN comments and blank lines; also adds a header note describing appended columns.
    """
    na_vals = ["NA"] * 6
    data_row_idx = 0

    with open(scn_out, 'w') as out, open(scn_in, 'r') as fh:
        wrote_header_note = False
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('#') or not line:
                out.write(line + "\n")
                continue
            parts = line.split('\t')
            data_row_idx += 1
            # Map row index -> FASTA id keys
            key = key_index.get(data_row_idx, None)
            vals = na_vals
            if key is not None:
                # Try exact header, normalized, and also trimmed variants (TEsorter typically uses the first token up to whitespace)
                candidates = [key["hdr"], key["norm"]]
                # Some TEsorter versions replace spaces and special chars with '_' in the first field only.
                more = [re.sub(r'\s+', '_', c) for c in candidates]
                candidates.extend(more)
                for c in candidates:
                    if c in tesorter_ann:
                        vals = tesorter_ann[c]
                        break
            # Write line + appended vals
            out.write("\t".join(parts + vals) + "\n")

        # After writing all lines, add a trailing comment noting the appended columns (once)
        # (Optional — harmless if left out; keeping it explicit helps downstream users.)
    # Prepend a header note at the top of the file about appended columns
    # We can’t easily "insert" at top without buffering; instead, write a sibling .header file or just rely on documentation.
    # Simpler: write a sidecar info note.
    info_path = scn_out + ".columns.txt"
    with open(info_path, 'w') as f:
        f.write("Appended columns (6): " + "\t".join(TE_COLS) + "\n")

############################################################
# CLI
############################################################

def main():
    ap = argparse.ArgumentParser(description="Extract LTR internal regions, run TEsorter, and append results to SCN.")
    ap.add_argument("--genome", required=True, help="Genome FASTA (same sequences as LTRharvest/FINDER).")
    ap.add_argument("--scn", required=True, help="LTR_retriever scn-format prediction file.")
    ap.add_argument("--outdir", default="tesorter_run", help="Output directory (default: tesorter_run).")
    ap.add_argument("--internal-fasta", default="LTR_internal.fa", help="Name for internal sequences FASTA (written inside outdir).")
    ap.add_argument("--include-ltr-boundaries", action="store_true",
                    help="Use [e(lLTR), s(rLTR)] inclusive (default excludes LTR ends).")
    ap.add_argument("--TEsorter-args", default="", help="Arguments passed through to TEsorter (quoted string).")
    ap.add_argument("--tesorter-dir", default="TEsorter", help="Path to clone/find the TEsorter repo (default: ./TEsorter).")
    ap.add_argument("--scn-out", default=None, help="Output SCN+TEsorter path (default: <scn>.with_tesorter.tsv)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    internal_fa_path = os.path.join(args.outdir, args.internal_fasta)
    scn_out_path = args.scn_out or (args.scn + ".with_tesorter.tsv")

    sys.stderr.write("[info] Extracting internal regions…\n")
    n_ok, skipped, key_index = write_internal_fasta(
        args.genome, args.scn, internal_fa_path, include_ltr_boundaries=args.include_ltr_boundaries
    )
    sys.stderr.write(f"[info] Wrote {n_ok} internal sequences to {internal_fa_path}\n")
    if skipped:
        sys.stderr.write(f"[warn] Skipped {len(skipped)} entries while extracting internals:\n")
        for i, why, meta in skipped[:10]:
            sys.stderr.write(f"  - row#{i}: {why} ({meta})\n")
        if len(skipped) > 10:
            sys.stderr.write(f"  … and {len(skipped)-10} more\n")

    sys.stderr.write("[info] Ensuring TEsorter is available…\n")
    if not ensure_tesorter(args.tesorter_dir):
        sys.stderr.write("[error] Could not prepare TEsorter. Exiting before classification.\n")
        sys.exit(2)

    sys.stderr.write("[info] Running TEsorter on internal sequences…\n")
    rc = run_tesorter(os.path.basename(internal_fa_path), args.TEsorter_args,
                      tesorter_dir=args.tesorter_dir, workdir=args.outdir)
    if rc != 0:
        sys.stderr.write(f"[error] TEsorter exited with code {rc}\n")

    # Figure out the .cls.tsv path
    db = detect_db_from_args(args.TEsorter_args)
    cls_path = os.path.join(args.outdir, f"{os.path.basename(internal_fa_path)}.{db}.cls.tsv")
    sys.stderr.write(f"[info] Reading TEsorter results: {cls_path}\n")
    tesorter_ann = read_tesorter_cls(cls_path)
    sys.stderr.write(f"[info] Loaded {len(tesorter_ann)} classification entries (with ID variants).\n")

    sys.stderr.write(f"[info] Appending TEsorter columns to SCN -> {scn_out_path}\n")
    append_tesorter_to_scn(args.scn, scn_out_path, key_index, tesorter_ann)
    sys.stderr.write(f"[info] Done. Appended columns: {', '.join(TE_COLS)}\n")
    sys.stderr.write(f"[info] Sidecar column note written to: {scn_out_path}.columns.txt\n")

if __name__ == "__main__":
    main()
