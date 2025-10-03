#!/usr/bin/env python3
import argparse, os, sys, shutil, tempfile, subprocess, textwrap, shlex
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

VERSION = "v1.2-py-min-mask+finder"

def which(cmd):
    return shutil.which(cmd)

def run(cmd, timeout=None, cwd=None):
    # shell=False; pass as list
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=timeout, text=True, cwd=cwd
    )
    return proc

def read_fasta_lengths(fa_path):
    lengths = {}
    order = {}
    seq_order = 0
    with open(fa_path) as f:
        name = None
        seq_chunks = []
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = sum(len(s) for s in seq_chunks)
                    order[name] = seq_order; seq_order += 1
                name = line.strip()[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            lengths[name] = sum(len(s) for s in seq_chunks)
            order[name] = seq_order
    return lengths, order

def chunk_windows(seq_len, size, overlap):
    assert size > 0
    step = max(1, size - overlap)
    starts = list(range(0, seq_len, step))
    windows = []
    for i, s in enumerate(starts, 1):
        e = min(seq_len, s + size)
        if windows and e == windows[-1][2]:  # avoid 0-length duplicate at end
            continue
        windows.append((i, s, e))
        if e == seq_len:
            break
    return windows

def write_chunk_fa(out_path, header, seq):
    with open(out_path, "w") as w:
        w.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            w.write(seq[i:i+80] + "\n")

def read_fasta_sequences(fa_path):
    # return dict name->sequence (no whitespace)
    seqs = {}
    with open(fa_path) as f:
        name = None
        buf = []
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line.strip()[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

def map_coords_back(reverse_flag, len_chunk, coord_adj, coords):
    # coords is an even-length list of positions (start,end, start,end, ...)
    out = []
    for i, val in enumerate(coords):
        if reverse_flag == 0:
            out.append(val + coord_adj)
        else:
            if i % 2 == 0:
                out.append(coord_adj + len_chunk - coords[i+1])
            else:
                out.append(coord_adj + len_chunk - coords[i-1])
    return out

def ensure_symlink(src, dst):
    src = str(Path(src).resolve())
    dst_p = Path(dst)
    if dst_p.exists():
        return
    os.symlink(src, dst)

# -------------- Masking helpers --------------

# ---------- FASTA utilities ----------

def fasta_iter(path):
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line.strip()[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            yield name, "".join(buf)

def write_fasta_record(fo, name, seq, width=80):
    fo.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fo.write(seq[i:i+width] + "\n")

# ---------- miniprot helpers ----------

def ensure_miniprot(tools_dir) -> str:
    """Return path to a runnable miniprot; build into tools_dir if not on PATH."""
    mp = which("miniprot")
    if mp:
        return mp
    tools_dir = Path(tools_dir).resolve()
    repo_dir = tools_dir / "miniprot"
    exe = repo_dir / "miniprot"
    tools_dir.mkdir(parents=True, exist_ok=True)
    if not exe.exists():
        try:
            if not repo_dir.exists():
                r = run(["git", "clone", "https://github.com/lh3/miniprot", str(repo_dir)])
                if r.returncode != 0:
                    raise RuntimeError(f"git clone failed: {r.stderr.strip()}")
            r = run(["make"], cwd=str(repo_dir))
            if r.returncode != 0:
                raise RuntimeError(f"make failed: {r.stderr.strip()}")
        except Exception as e:
            raise RuntimeError(f"Failed to install miniprot: {e}")
    if not (exe.exists() and os.access(exe, os.X_OK)):
        raise RuntimeError("miniprot binary not found or not executable after build.")
    return str(exe)

def parse_mrna_intervals_from_gff_text(gff_text):
    """Return dict: chrom -> list[(start,end)] for type == mRNA (1-based, inclusive)."""
    intervals = {}
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        typ = parts[2]
        if typ != "mRNA":
            continue
        chrom = parts[0]
        try:
            s = int(parts[3]); e = int(parts[4])
        except Exception:
            continue
        if s > e:
            s, e = e, s
        intervals.setdefault(chrom, []).append((s, e))
    return intervals

def merge_intervals_per_chrom(ivd):
    """Merge overlapping/adjacent intervals for each chrom independently."""
    out = {}
    for chrom, ivals in ivd.items():
        if not ivals:
            continue
        ivals.sort(key=lambda x: x[0])
        m = [list(ivals[0])]
        for s, e in ivals[1:]:
            if s <= m[-1][1] + 0:  # adjacent ok
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[chrom] = [(s, e) for s, e in m]
    return out

def hardmask_fasta_by_intervals(in_fasta, intervals_by_chrom, out_fasta, bed_out):
    """Mask [start,end] (1-based, inclusive) with 'N' and write BED (0-based, half-open)."""
    with open(out_fasta, "w") as fo, open(bed_out, "w") as bo:
        for contig, seq in fasta_iter(in_fasta):
            if contig not in intervals_by_chrom or not intervals_by_chrom[contig]:
                write_fasta_record(fo, contig, seq)
                continue
            L = len(seq)
            sl = list(seq)
            for s, e in intervals_by_chrom[contig]:
                s0 = max(1, min(s, L))
                e0 = max(1, min(e, L))
                if s0 > e0:
                    s0, e0 = e0, s0
                for i in range(s0 - 1, e0):
                    sl[i] = "N"
                bo.write(f"{contig}\t{s0-1}\t{e0}\n")
            write_fasta_record(fo, contig, "".join(sl))

def run_miniprot_gene_mask(in_fasta, protein_faa, out_prefix, threads, miniprot_path,
                           outn=1, outs=0.8, outc=0.5):
    """
    Run miniprot, capture GFF text, mask mRNA spans, return paths:
      (gff_path, bed_path, masked_fasta)
    """
    gff_path = f"{out_prefix}.genic.gff"
    bed_path = f"{out_prefix}.genic.mask.bed"
    masked_fasta = f"{out_prefix}.genic_masked.fa"

    cmd = [
        miniprot_path,
        "--gff-only",
        "-t", str(threads),
        in_fasta,
        protein_faa,
        "-P", out_prefix,
        "--outn", str(outn),
        "--outs", str(outs),
        "--outc", str(outc),
    ]
    r = run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"miniprot failed:\n{r.stderr.strip()}")
    gff_text = r.stdout
    with open(gff_path, "w") as gfh:
        gfh.write(gff_text)

    intervals = parse_mrna_intervals_from_gff_text(gff_text)
    intervals = merge_intervals_per_chrom(intervals)
    if not intervals:
        # Nothing to mask; copy input
        shutil.copyfile(in_fasta, masked_fasta)
        Path(bed_path).write_text("")  # empty bed
        return gff_path, bed_path, masked_fasta

    hardmask_fasta_by_intervals(in_fasta, intervals, masked_fasta, bed_path)
    return gff_path, bed_path, masked_fasta

def ensure_tool(tool_name, tools_dir):
    """
    Ensure a low-complexity tool binary exists and return its path.
    If not on PATH, try to auto-install into tools_dir by cloning & compiling.
    """
    # Already in PATH?
    p = which(tool_name)
    if p:
        return p

    tools_dir = Path(tools_dir).resolve()
    tools_dir.mkdir(parents=True, exist_ok=True)

    if tool_name == "sdust":
        repo = "https://github.com/lh3/sdust.git"
        subdir = tools_dir / "sdust"
        exe = subdir / "sdust"
        if not exe.exists():
            clone_build(repo, subdir, ["make"])
        return str(exe) if exe.exists() else None

    if tool_name == "trf-mod":
        # Executable name is "trf-mod" in repo TRF-mod
        repo = "https://github.com/lh3/TRF-mod.git"
        subdir = tools_dir / "TRF-mod"
        exe = subdir / "trf-mod"
        if not exe.exists():
            clone_build(repo, subdir, ["make", "-f", "compile.mak"])
        return str(exe) if exe.exists() else None

    if tool_name == "longdust":
        repo = "https://github.com/lh3/longdust.git"
        subdir = tools_dir / "longdust"
        exe = subdir / "longdust"
        if not exe.exists():
            clone_build(repo, subdir, ["make"])
        return str(exe) if exe.exists() else None

    return None

def clone_build(repo_url, dest_dir, build_cmd):
    try:
        if not dest_dir.exists():
            r = run(["git", "clone", repo_url, str(dest_dir)])
            if r.returncode != 0:
                print(f"[WARN] git clone failed for {repo_url}: {r.stderr.strip()}", file=sys.stderr)
                return
        r = run(build_cmd, cwd=str(dest_dir))
        if r.returncode != 0:
            print(f"[WARN] build failed in {dest_dir}: {r.stderr.strip()}", file=sys.stderr)
    except Exception as e:
        print(f"[WARN] exception during install of {repo_url}: {e}", file=sys.stderr)

def read_bed_intervals(bed_path, chrom_expected):
    """
    Parse BED-like file, return list of [start,end) for matching chrom.
    Ignores malformed lines. Accepts TRF-mod extended columns.
    """
    ivals = []
    if not Path(bed_path).exists():
        return ivals
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom != chrom_expected:
                # Some tools echo original header name; we generated chunks named {chrom}_subN.
                # For safety, accept either exact chunk header or original chrom if tool transforms names.
                if chrom != chrom_expected.split("_sub")[0]:
                    continue
            try:
                s = int(parts[1]); e = int(parts[2])
                if e > s:
                    ivals.append([s, e])
            except:
                continue
    return ivals

def merge_intervals(ivals):
    if not ivals:
        return []
    ivals.sort(key=lambda x: x[0])
    out = [ivals[0]]
    for s, e in ivals[1:]:
        if s <= out[-1][1]:  # overlap/adjacent
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

def gene_mask_from_existing_gff(in_fasta, gff_in_path, out_prefix):
    """
    Use an existing miniprot GFF to hard-mask mRNA spans on the full genome.
    Returns (gff_path_out, bed_path, masked_fasta).
    """
    gff_path_out = f"{out_prefix}.genic.gff"
    bed_path = f"{out_prefix}.genic.mask.bed"
    masked_fasta = f"{out_prefix}.genic_masked.fa"

    # ensure output GFF exists adjacent to other outputs (copy for consistency)
    if str(gff_in_path) != gff_path_out:
        shutil.copyfile(gff_in_path, gff_path_out)
    else:
        # still ensure it's present; nothing to do
        pass

    gff_text = Path(gff_in_path).read_text()
    intervals = parse_mrna_intervals_from_gff_text(gff_text)
    intervals = merge_intervals_per_chrom(intervals)

    if not intervals:
        # Nothing to mask; copy input
        shutil.copyfile(in_fasta, masked_fasta)
        Path(bed_path).write_text("")  # empty bed
        return gff_path_out, bed_path, masked_fasta

    hardmask_fasta_by_intervals(in_fasta, intervals, masked_fasta, bed_path)
    return gff_path_out, bed_path, masked_fasta


def hardmask_fasta_in_place(fa_path, intervals):
    """
    Replace bases in intervals with 'N'. Assumes single-sequence FASTA.
    """
    if not intervals:
        return
    # read
    header = None
    seq_chunks = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:].split()[0]
            else:
                seq_chunks.append(line.strip())
    seq = list("".join(seq_chunks))
    L = len(seq)
    for s, e in intervals:
        s = max(0, min(L, s))
        e = max(0, min(L, e))
        for i in range(s, e):
            seq[i] = "N"
    # write back (80 cols)
    with open(fa_path, "w") as w:
        w.write(f">{header}\n")
        s = "".join(seq)
        for i in range(0, len(s), 80):
            w.write(s[i:i+80] + "\n")

# ---------- LTR_FINDER helpers ----------

def ensure_ltr_finder(tools_dir):
    """
    Clone LTR_FINDER_parallel if needed and return (ltr_finder_bin, convert_script) paths.
    Prefers the repository layout the user specified.
    """
    tools_dir = Path(tools_dir).resolve()
    repo_dir = tools_dir / "LTR_FINDER_parallel"
    if not repo_dir.exists():
        r = run(["git", "clone", "https://github.com/cwb14/LTR_FINDER_parallel.git", str(repo_dir)])
        if r.returncode != 0:
            raise RuntimeError(f"git clone LTR_FINDER_parallel failed: {r.stderr.strip()}")

    # Try the given canonical path first
    lf_bin = repo_dir / "bin" / "LTR_FINDER.x86_64-1.0.7" / "ltr_finder"
    conv = repo_dir / "bin" / "convert_ltr_finder2.pl"

    # Fallback: search for ltr_finder under bin/
    if not lf_bin.exists():
        for p in (repo_dir / "bin").rglob("ltr_finder"):
            lf_bin = p
            break

    if not lf_bin.exists():
        raise RuntimeError("Could not locate ltr_finder binary under LTR_FINDER_parallel/bin")

    if not conv.exists():
        raise RuntimeError("Could not locate convert_ltr_finder2.pl under LTR_FINDER_parallel/bin")

    if not os.access(lf_bin, os.X_OK):
        raise RuntimeError(f"ltr_finder is not executable: {lf_bin}")

    return str(lf_bin), str(conv)

# -------------- Main --------------

def main():
    ap = argparse.ArgumentParser(
        description="Run GenomeTools LTRharvest over a genome in parallel chunks (with optional low-complexity and genic hardmasking), and optionally run LTR_FINDER per chunk with conversion to LTRharvest format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("-seq", required=True, help="Input FASTA")
    ap.add_argument("-size", type=int, default=5_000_000, help="Chunk size (bp)")
    ap.add_argument("--overlap", type=int, default=0, help="Overlap between consecutive chunks (bp)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="Threads for parallel processing")
    ap.add_argument("--timeout", type=int, default=1500, help="Per-chunk timeout (seconds)")
    ap.add_argument("--gt-path", default="", help="Path to `gt` (GenomeTools). If empty, will use $PATH")
    ap.add_argument("--annotator", default="LTR_HARVEST_parallel", help="Annotator name in GFF for LTRharvest results")
    ap.add_argument(
        "--ltrharvest-args",
        default="",
        help=("Arguments string forwarded to `gt ltrharvest`. "
              "If omitted, runs `gt ltrharvest` with no additional options. "
              "Tip: see `gt ltrharvest -help`.")
    )
    # LTR_FINDER options
    ap.add_argument("--run-ltrfinder", action="store_true",
                    help="Also run LTR_FINDER on each chunk and merge its calls.")
    ap.add_argument("--ltrfinder-args", default="",
                    help="Arguments string forwarded to `ltr_finder` (placed before the input FASTA).")
    ap.add_argument("--ltrfinder-path", default="",
                    help="Path to an ltr_finder binary; if empty, will clone into ./tools/LTR_FINDER_parallel and use bundled binary.")
    ap.add_argument("--annotator-finder", default="LTR_FINDER_parallel",
                    help="Annotator name in GFF for LTR_FINDER results")
    ap.add_argument("--keep", action="store_true", help="Keep chunk folder (do not delete)")
    ap.add_argument("--workdir", default="", help="Optional working directory (default: alongside input)")
    # New: masking toggles
    ap.add_argument("--sdust", action="store_true", help="Hard-mask sdust regions per chunk before LTRharvest/LTR_FINDER")
    ap.add_argument("--trf", action="store_true", help="Hard-mask TRF-mod regions per chunk before LTRharvest/LTR_FINDER")
    ap.add_argument("--longdust", action="store_true", help="Hard-mask longdust regions per chunk before LTRharvest/LTR_FINDER")
    # Gene masking (whole-genome, pre-partition) via miniprot
    ap.add_argument("--gene-mask", action="store_true",
                    help="Hard-mask genic (mRNA) spans on the full genome using miniprot before partitioning.")
    ap.add_argument("--protein", default="",
                    help="Protein FASTA to guide genic masking (required with --gene-mask).")
    ap.add_argument("--miniprot-path", default="",
                    help="Path to a miniprot binary; if empty, will attempt to build into ./tools/miniprot.")
    ap.add_argument("--mp-outn", type=int, default=1, help="miniprot --outn")
    ap.add_argument("--mp-outs", type=float, default=0.8, help="miniprot --outs")
    ap.add_argument("--mp-outc", type=float, default=0.5, help="miniprot --outc")
    ap.add_argument(
        "--gene-gff",
        default="",
        help="Use an existing miniprot GFF to hard-mask genic (mRNA) spans on the full genome and skip running miniprot."
    )

    args = ap.parse_args()

    # Locate gt
    gt_bin = args.gt_path.strip()
    if gt_bin:
        gt_exe = Path(gt_bin)
        if gt_exe.is_dir() or gt_exe.name == "gt":
            if gt_exe.is_dir():
                gt = str(gt_exe / "gt")
            else:
                gt = str(gt_exe)
        else:
            gt = gt_bin
    else:
        gt = which("gt") or ""
    if not gt:
        sys.exit("ERROR: Could not find GenomeTools `gt` in PATH. Use --gt-path.")
    if not os.access(gt, os.X_OK):
        sys.exit(f"ERROR: `gt` is not executable at: {gt}")
    print(f"Using LTRharvest: {gt}", file=sys.stderr)

    fasta = Path(args.seq).resolve()
    if not fasta.exists():
        sys.exit(f"ERROR: FASTA not found: {fasta}")

    # Prepare workspace (do this early so we have a place to build tools & write masked FASTA)
    if args.workdir:
        workbase = Path(args.workdir).resolve()
        workbase.mkdir(parents=True, exist_ok=True)
    else:
        workbase = Path.cwd()
    chunk_dir = workbase / f"{fasta.name}.harvest"
    chunk_dir.mkdir(exist_ok=True)
    tools_dir = workbase / "tools"
    tools_dir.mkdir(exist_ok=True)

    # LTR_FINDER setup if requested
    ltr_finder_bin = ""
    converter_pl = ""
    if args.run_ltrfinder:
        if args.ltrfinder_path.strip():
            ltr_finder_bin = args.ltrfinder_path.strip()
            if not Path(ltr_finder_bin).exists():
                sys.exit(f"ERROR: --ltrfinder-path provided but not found: {ltr_finder_bin}")
            # Try to find converter alongside the repo root if user provided custom path
            # Assume default layout unless user overrides by PATHing the script themselves.
            default_conv = tools_dir / "LTR_FINDER_parallel" / "bin" / "convert_ltr_finder2.pl"
            converter_pl = str(default_conv) if default_conv.exists() else ""
        if not ltr_finder_bin:
            try:
                ltr_finder_bin, converter_pl = ensure_ltr_finder(tools_dir)
            except Exception as e:
                sys.exit(f"ERROR setting up LTR_FINDER_parallel: {e}")
        if not converter_pl or not Path(converter_pl).exists():
            sys.exit("ERROR: convert_ltr_finder2.pl could not be located. If using --ltrfinder-path, also ensure the converter script is present under tools/LTR_FINDER_parallel/bin.")

        print(f"Using LTR_FINDER: {ltr_finder_bin}", file=sys.stderr)
        print(f"Using converter : {converter_pl}", file=sys.stderr)

    # (NEW) Whole-genome genic hardmasking (either via existing GFF or by running miniprot)
    fasta_for_processing = fasta  # may be replaced if gene masking requested
    out_prefix = str(workbase / fasta.name)  # keep outputs adjacent to work

    if args.gene_gff:
        gff_in = Path(args.gene_gff).resolve()
        if not gff_in.exists():
            sys.exit(f"ERROR: --gene-gff provided but file not found: {gff_in}")
        if args.gene_mask:
            print("[WARN] Both --gene-gff and --gene-mask provided; using --gene-gff and skipping miniprot.", file=sys.stderr)
        try:
            gff_p, bed_p, masked_fa = gene_mask_from_existing_gff(
                str(fasta), str(gff_in), out_prefix
            )
            print(f"[GENE-MASK] existing GFF: {gff_p}", file=sys.stderr)
            print(f"[GENE-MASK] masked BED : {bed_p}", file=sys.stderr)
            print(f"[GENE-MASK] masked FASTA: {masked_fa}", file=sys.stderr)
            fasta_for_processing = Path(masked_fa)
        except Exception as e:
            sys.exit(f"ERROR during gene masking from existing GFF: {e}")

    elif args.gene_mask:
        if not args.protein:
            sys.exit("ERROR: --gene-mask requested but --protein FASTA not provided.")
        try:
            mp_bin = args.miniprot_path.strip() or ensure_miniprot(tools_dir)
        except Exception as e:
            sys.exit(f"ERROR: {e}")
        print(f"Using miniprot: {mp_bin}", file=sys.stderr)

        try:
            gff_p, bed_p, masked_fa = run_miniprot_gene_mask(
                str(fasta), args.protein, out_prefix,
                threads=args.threads, miniprot_path=mp_bin,
                outn=args.mp_outn, outs=args.mp_outs, outc=args.mp_outc
            )
            print(f"[GENE-MASK] miniprot GFF: {gff_p}", file=sys.stderr)
            print(f"[GENE-MASK] masked BED : {bed_p}", file=sys.stderr)
            print(f"[GENE-MASK] masked FASTA: {masked_fa}", file=sys.stderr)
            fasta_for_processing = Path(masked_fa)
        except Exception as e:
            sys.exit(f"ERROR during gene masking: {e}")
    # else: no gene masking requested

    # Outputs (keep original base names)
    out_scn = Path(f"{fasta.name}.harvest.combine.scn")
    out_gff = Path(f"{fasta.name}.harvest.combine.gff3")

    # Read lengths + order & sequences from the FASTA actually used for downstream processing
    lengths, order_map = read_fasta_lengths(fasta_for_processing)
    seqs = read_fasta_sequences(fasta_for_processing)

    # Prepare workspace (redundant but safe)
    if args.workdir:
        workbase = Path(args.workdir).resolve()
        workbase.mkdir(parents=True, exist_ok=True)
    else:
        workbase = Path.cwd()
    chunk_dir = workbase / f"{fasta.name}.harvest"
    chunk_dir.mkdir(exist_ok=True)

    # Where to install tools if needed
    tools_dir = workbase / "tools"
    tools_dir.mkdir(exist_ok=True)

    # Resolve selected maskers + ensure binaries
    selected_maskers = []
    if args.sdust:
        sd = ensure_tool("sdust", tools_dir)
        if sd: selected_maskers.append(("sdust", sd))
        else: print("[WARN] sdust not available; skipping sdust masking", file=sys.stderr)
    if args.trf:
        trf = ensure_tool("trf-mod", tools_dir)
        if trf: selected_maskers.append(("trf", trf))
        else: print("[WARN] TRF-mod not available; skipping TRF-mod masking", file=sys.stderr)
    if args.longdust:
        ld = ensure_tool("longdust", tools_dir)
        if ld: selected_maskers.append(("longdust", ld))
        else: print("[WARN] longdust not available; skipping longdust masking", file=sys.stderr)

    # Build chunks (with overlap)
    list_entries = []  # [(chunk_basename, chrom, piece_index, start0, end0, len_chunk)]
    for chrom, seq in seqs.items():
        L = len(seq)
        wins = chunk_windows(L, args.size, args.overlap)
        for piece_idx, s0, e0 in wins:
            chunk_base = f"{chrom}_sub{piece_idx}"
            chunk_fa = chunk_dir / f"{chunk_base}.fa"
            write_chunk_fa(chunk_fa, chunk_base, seq[s0:e0])
            list_entries.append((chunk_base, chrom, piece_idx, s0, e0, e0 - s0))

    # Worker for each chunk
    def worker(entry):
        (chunk_base, chrom, piece_idx, s0, e0, len_chunk) = entry
        chunk_fa = chunk_dir / f"{chunk_base}.fa"
        masked_fa = chunk_fa  # overwrite in place once masked
        idxname = chunk_dir / chunk_base
        scn_path_harvest = chunk_dir / f"{chunk_base}.harvest.scn"
        scn_path_finder  = chunk_dir / f"{chunk_base}.ltrfinder.scn"

        # Optional masking
        if selected_maskers:
            all_ivals = []
            for name, exe in selected_maskers:
                try:
                    if name == "sdust":
                        out_bed = chunk_dir / f"{chunk_base}_sdust.bed"
                        r = run([exe, str(chunk_fa), "-t", "12"])
                        if r.returncode != 0:
                            return (chunk_base, False, f"sdust failed: {r.stderr.strip()}")
                        out_bed.write_text(r.stdout)
                        ivals = read_bed_intervals(out_bed, chunk_base)
                        all_ivals.extend(ivals)
                    elif name == "trf":
                        out_bed = chunk_dir / f"{chunk_base}_trf.bed"
                        r = run([exe, str(chunk_fa), "-b5", "-g5", "-s30", "-p200"])
                        if r.returncode != 0:
                            return (chunk_base, False, f"TRF-mod failed: {r.stderr.strip()}")
                        out_bed.write_text(r.stdout)
                        ivals = read_bed_intervals(out_bed, chunk_base)
                        all_ivals.extend(ivals)
                    elif name == "longdust":
                        out_bed = chunk_dir / f"{chunk_base}_longdust.bed"
                        r = run([exe, str(chunk_fa)])
                        if r.returncode != 0:
                            return (chunk_base, False, f"longdust failed: {r.stderr.strip()}")
                        out_bed.write_text(r.stdout)
                        ivals = read_bed_intervals(out_bed, chunk_base)
                        all_ivals.extend(ivals)
                except Exception as e:
                    return (chunk_base, False, f"{name} exception: {e}")

            merged = merge_intervals(all_ivals)
            try:
                hardmask_fasta_in_place(masked_fa, merged)
            except Exception as e:
                return (chunk_base, False, f"masking failed: {e}")

        # suffixerator
        cmd1 = [gt, "suffixerator",
                "-db", str(masked_fa),
                "-indexname", str(idxname),
                "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna"]
        r1 = run(cmd1, timeout=args.timeout)
        if r1.returncode != 0:
            return (chunk_base, False, f"suffixerator failed: {r1.stderr.strip()}")

        # ltrharvest
        lh_tokens = shlex.split(args.ltrharvest_args) if args.ltrharvest_args.strip() else []
        cmd2 = [gt, "ltrharvest", "-index", str(idxname)] + lh_tokens
        r2 = run(cmd2, timeout=args.timeout)
        if r2.returncode != 0:
            return (chunk_base, False, f"ltrharvest failed: {r2.stderr.strip()}")

        with open(scn_path_harvest, "w") as w:
            w.write(r2.stdout)

        # Best-effort cleanup of large index files
        for ext in (".des",".esq",".lcp",".llv",".md5",".prj",".sds",".ssp",".suf"):
            p = chunk_dir / f"{chunk_base}{ext}"
            try:
                if p.exists(): p.unlink()
            except Exception:
                pass

        # LTR_FINDER (optional)
        if args.run_ltrfinder:
            try:
                raw_out = chunk_dir / f"{chunk_base}.ltrfinder.raw"
                lf_tokens = shlex.split(args.ltrfinder_args) if args.ltrfinder_args.strip() else []
                cmd_lf = [ltr_finder_bin] + lf_tokens + [str(masked_fa)]
                rlf = run(cmd_lf, timeout=args.timeout)
                if rlf.returncode != 0:
                    return (chunk_base, False, f"ltr_finder failed: {rlf.stderr.strip()}")
                raw_out.write_text(rlf.stdout)

                # Convert to LTRharvest format (SCN)
                # The converter prints to STDOUT
                cmd_conv = ["perl", converter_pl, str(raw_out)]
                rcv = run(cmd_conv, timeout=args.timeout)
                if rcv.returncode != 0:
                    return (chunk_base, False, f"convert_ltr_finder2.pl failed: {rcv.stderr.strip()}")
                scn_path_finder.write_text(rcv.stdout)

            except Exception as e:
                return (chunk_base, False, f"LTR_FINDER/convert exception: {e}")

        return (chunk_base, True, "")

    # Run in parallel
    print(f"Chunk count: {len(list_entries)}", file=sys.stderr)
    if selected_maskers:
        used = ",".join(n for n,_ in selected_maskers)
        print(f"Masking enabled: {used}", file=sys.stderr)
    if args.run_ltrfinder:
        print("LTR_FINDER will also be run per chunk.", file=sys.stderr)

    status = {}
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        fut2entry = {ex.submit(worker, e): e for e in list_entries}
        for fut in as_completed(fut2entry):
            base = fut2entry[fut][0]
            ok = False; msg = ""
            try:
                _base, ok, msg = fut.result()
            except Exception as e:
                msg = str(e)
            status[base] = (ok, msg)
            if ok:
                print(f"[OK] {base}", file=sys.stderr)
            else:
                print(f"[FAILED] {base} :: {msg}", file=sys.stderr)

    # Write headers & sequence-region lines in GFF
    with open(out_gff, "w") as gff:
        gff.write("##gff-version   3\n")
        for chrom, L in lengths.items():
            gff.write(f"##sequence-region   {chrom} 1 {L}\n")
        chrom_order_section = "".join(f"#{c}\n" for c in lengths.keys())
        gff.write(chrom_order_section)

    # Combine .scn (harvest first, then finder) and produce GFF features
    lh_args_display = args.ltrharvest_args if args.ltrharvest_args.strip() else "(none)"
    lf_args_display = args.ltrfinder_args if args.ltrfinder_args.strip() else "(none)"
    header = textwrap.dedent(f"""\
        #LTR_parallel -seq {fasta.name} -size {args.size} -time {args.timeout} -threads {args.threads}
        # LTRharvest args= {lh_args_display}
        # LTR_FINDER args= {lf_args_display if args.run_ltrfinder else '(not run)'}
        # Python version={VERSION}
        # predictions are reported in the following way
        # s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
        # where:
        # s = starting position
        # e = ending position
        # l = length
        # ret = LTR-retrotransposon
        # lLTR = left LTR
        # rLTR = right LTR
        # sim = similarity
        # seq-nr = sequence order
    """)

    count = 1
    with open(out_scn, "w") as outscn, open(out_gff, "a") as gff:
        outscn.write(header)

        def merge_one_scn(scn_file, annotator_name, reverse_hint=False):
            nonlocal count
            if not scn_file.exists():
                return
            with open(scn_file) as scn:
                for line in scn:
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    # Attempt to read seq_id (0 forward, 1 reverse) as in LTRharvest SCN
                    try:
                        seq_id_field = parts[10]
                        seq_id = int(seq_id_field)
                    except Exception:
                        seq_id = 0  # assume forward if missing

                    try:
                        start = int(parts[0]); end = int(parts[1]); ele_len = int(parts[2])
                        ltrt_s = int(parts[3]); lltr_e = int(parts[4]); lltr = int(parts[5])
                        rltr_s = int(parts[6]); rltr_e = int(parts[7]); rltr = int(parts[8])
                        sim = parts[9]
                        chrom = parts[11] if len(parts) > 11 else None
                    except Exception:
                        continue

                    # Figure out which chunk this SCN belongs to to recover coord_adj/len_chunk.
                    # We encoded scn_file name as <chrom>_sub<i>.{harvest|ltrfinder}.scn; use that.
                    # Since we need the original chrom and offsets, rely on list_entries.
                    # (This loop is called per chunk below, so we pass coord context when calling.)
                    return (start, end, ele_len, ltrt_s, lltr_e, lltr, rltr_s, rltr_e, rltr, sim, seq_id, chrom)

            return None

        # Iterate chunk-by-chunk to preserve coordinate adjustments
        for (chunk_base, chrom, piece_idx, s0, e0, len_chunk) in list_entries:
            coord_adj = s0

            # Helper to process a single SCN file (harvest or finder)
            def process_scn_file(scn_path, annotator_name):
                nonlocal count
                if not scn_path.exists():
                    return
                with open(scn_path) as scn:
                    for line in scn:
                        if not line.strip() or line.startswith("#"):
                            continue
                        parts = line.strip().split()
                        # Parse SCN columns; tolerate slightly different tails
                        try:
                            sim = parts[9]
                        except Exception:
                            continue
                        # seq_id for strand (0 forward, 1 reverse) if present
                        seq_id = 0
                        for v in reversed(parts[-2:]):
                            try:
                                seq_id = int(v); break
                            except: continue

                        try:
                            start = int(parts[0]); end = int(parts[1]); ele_len = int(parts[2])
                            ltrt_s = int(parts[3]); lltr_e = int(parts[4]); lltr = int(parts[5])
                            rltr_s = int(parts[6]); rltr_e = int(parts[7]); rltr = int(parts[8])
                        except Exception:
                            continue

                        coords = [start, end, ltrt_s, lltr_e, rltr_s, rltr_e]
                        mapped = map_coords_back(0 if seq_id == 0 else 1, len_chunk, coord_adj, coords)
                        start, end, ltrt_s, lltr_e, rltr_s, rltr_e = mapped

                        outscn.write(f"{start} {end} {ele_len} {ltrt_s} {lltr_e} {lltr} {rltr_s} {rltr_e} {rltr} {sim} {order_map[chrom]} {chrom}\n")

                        strand = "+" if seq_id == 0 else "-"
                        gff.write(f"{chrom}\t{annotator_name}\trepeat_region\t{start}\t{end}\t.\t{strand}\t.\tID=repeat_region{count}\n")
                        gff.write(f"{chrom}\t{annotator_name}\tLTR_retrotransposon\t{start}\t{end}\t.\t{strand}\t.\tID=LTR_retrotransposon{count};Parent=repeat_region{count};ltr_identity={sim};seq_number={order_map[chrom]}\n")
                        gff.write(f"{chrom}\t{annotator_name}\tlong_terminal_repeat\t{start}\t{lltr_e}\t.\t{strand}\t.\tParent=LTR_retrotransposon{count}\n")
                        gff.write(f"{chrom}\t{annotator_name}\tlong_terminal_repeat\t{rltr_s}\t{end}\t.\t{strand}\t.\tParent=LTR_retrotransposon{count}\n")
                        gff.write("###\n")

                        count += 1

            # Harvest SCN for this chunk
            scn_file_h = chunk_dir / f"{chunk_base}.harvest.scn"
            process_scn_file(scn_file_h, args.annotator)

            # Finder SCN for this chunk (if present)
            if args.run_ltrfinder:
                scn_file_f = chunk_dir / f"{chunk_base}.ltrfinder.scn"
                process_scn_file(scn_file_f, args.annotator_finder)

    if not args.keep:
        try:
            shutil.rmtree(chunk_dir)
        except Exception:
            pass

    print(f"Done. Outputs:\n  {out_scn}\n  {out_gff}", file=sys.stderr)

if __name__ == "__main__":
    main()
