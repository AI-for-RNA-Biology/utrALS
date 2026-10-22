#!/usr/bin/env python3
"""
Count genomic sites per sliding window along a target gene.

Supports two input formats for sites:
  • BED    – tab-separated, no header; first 6 columns used
             (chrom  start  end  name  score  strand)
  • GFF/GTF – loaded with genomemanager.gtf_utils.read_gtf
  Format is auto-detected from the file extension (.gff / .gtf → GFF/GTF;
  everything else → BED).  Extra columns are silently ignored, so the
  29-column RMBase m6A BED file works out of the box.

Sub-commands
------------
single   – process a single sites file and write one output CSV
count    – print the number of CLIP files found in --iclip-dir / --eclip-dir
prepare  – load the annotation GTF, cache target gene coordinates to JSON
run      – process one CLIP file (--task-id or $SLURM_ARRAY_TASK_ID) and
           write <clip_name>_w{window_size}.csv to --output-dir

Output CSV columns
------------------
  window_start  – position relative to the gene start (nt)
  window_end    – position relative to the gene start (nt)
  n_sites       – number of sites whose midpoint falls in the window
  n_sites_norm  – n_sites divided by the total number of sites in the whole
                  CLIP experiment (genome-wide, before subsetting to the
                  gene), scaled to sites per million total sites (SPM). This
                  makes density profiles comparable across CLIP experiments
                  that differ in total number of called sites.
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from genomemanager.gtf_utils import read_gtf


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_gene_from_annotation(annotation_path: str, gene_name: str) -> dict:
    """Return genomic coordinates of *gene_name* from a GTF annotation."""
    print(f"Loading annotation: {annotation_path}", flush=True)
    gtf = read_gtf(annotation_path)
    hits = gtf[(gtf["gene_name"] == gene_name) & (gtf["feature"] == "gene")]
    if hits.empty:
        sys.exit(f"ERROR: gene '{gene_name}' not found in {annotation_path}")
    row = hits.iloc[0]
    return {
        "seqname": str(row["seqname"]),
        "start":   int(row["start"]),
        "end":     int(row["end"]),
        "strand":  str(row["strand"]),
    }


def gene_cache_path(output_dir: str, gene_name: str) -> str:
    return os.path.join(output_dir, f"_gene_cache_{gene_name}.json")


def load_gene_cached(output_dir: str, gene_name: str,
                     annotation_path: str = None) -> dict:
    """Load gene coords from JSON cache (fast) or full annotation (slow)."""
    cache_file = gene_cache_path(output_dir, gene_name)
    if os.path.exists(cache_file):
        with open(cache_file) as fh:
            gene = json.load(fh)
        print(f"Gene info loaded from cache: {cache_file}", flush=True)
        return gene
    if annotation_path:
        return load_gene_from_annotation(annotation_path, gene_name)
    sys.exit(f"ERROR: gene cache not found at {cache_file} "
             "and --annotation was not provided")


def load_sites(sites_path: str) -> pd.DataFrame:
    """Load sites from a BED or GFF/GTF file; returns seqname/start/end/strand."""
    ext = os.path.splitext(sites_path)[1].lower()
    if ext in (".gff", ".gtf"):
        df = read_gtf(sites_path)
        return df[["seqname", "start", "end", "strand"]].copy()
    else:
        bed_cols = ["seqname", "start", "end", "name", "score", "strand"]
        return pd.read_csv(sites_path, sep="\t", header=None,
                           usecols=range(6), names=bed_cols)


def compute_density(gene: dict, sites: pd.DataFrame,
                    window_size: int, total_sites: int = None) -> pd.DataFrame:
    """Return a DataFrame with site counts per window.

    total_sites, if given, is the total number of sites in the whole CLIP
    experiment (i.e. len(sites) before subsetting to the gene); it is used
    to add a sites-per-million (SPM) normalized column so that density
    profiles can be compared across CLIP experiments of different depth.
    """
    in_gene = sites[
        (sites["seqname"] == gene["seqname"]) &
        (sites["strand"]  == gene["strand"])  &
        (sites["start"]   >= gene["start"])   &
        (sites["end"]     <= gene["end"])
    ]
    print(f"  → {len(in_gene):,} sites overlap the gene", flush=True)

    gene_len = gene["end"] - gene["start"]
    midpoints = ((in_gene["start"] + in_gene["end"]) / 2 - gene["start"]).values

    bins = np.arange(0, gene_len + window_size, window_size)
    counts, edges = np.histogram(midpoints, bins=bins)

    density = pd.DataFrame({
        "window_start": edges[:-1].astype(int),
        "window_end":   edges[1:].astype(int),
        "n_sites":      counts,
    })

    if total_sites is None:
        total_sites = len(sites)
    density["n_sites_norm"] = density["n_sites"] / total_sites * 1e6 if total_sites else 0.0

    return density


def get_all_clip_files(iclip_dir=None, eclip_dir=None):
    """Return sorted list of (filepath, label) for all .gff CLIP files."""
    entries = []
    if iclip_dir:
        for fp in sorted(glob.glob(os.path.join(iclip_dir, "*.gff"))):
            entries.append((fp, "iCLIP"))
    if eclip_dir:
        for fp in sorted(glob.glob(os.path.join(eclip_dir, "HepG2", "*.gff"))):
            entries.append((fp, "HepG2"))
        for fp in sorted(glob.glob(os.path.join(eclip_dir, "K652", "*.gff"))):
            entries.append((fp, "K652"))
    return entries


def clip_display_name(filepath: str, source: str) -> str:
    """Derive a short human-readable label from the filename and source."""
    basename = os.path.basename(filepath).replace(".gff", "")
    if source == "iCLIP":
        if "all" in basename:
            return basename.split("_all-")[1].split("-")[0] + "_iCLIP"
        return basename.split("_")[2].split("-")[0] + "_iCLIP"
    if source in ("HepG2", "K652"):
        return basename.split("_ENC")[0] + f"_{source}"
    return basename


# ── Sub-commands ──────────────────────────────────────────────────────────────

def cmd_single(args):
    """Process a single BED/GFF sites file and write one density CSV."""
    gene = load_gene_from_annotation(args.annotation, args.gene_name)
    print(f"Loading sites: {args.sites_file}", flush=True)
    sites = load_sites(args.sites_file)
    print(f"  → {len(sites):,} sites loaded", flush=True)
    density = compute_density(gene, sites, args.window_size)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    density.to_csv(args.output, index=False)
    gene_len = gene["end"] - gene["start"]
    print(f"Gene length : {gene_len:,} nt")
    print(f"Window size : {args.window_size} nt  →  {len(density)} windows")
    print(f"Output saved: {args.output}")


def cmd_count(args):
    """Print the total number of CLIP files (stdout, for shell capture)."""
    files = get_all_clip_files(args.iclip_dir, args.eclip_dir)
    print(len(files))


def cmd_prepare(args):
    """Load the annotation GTF, cache target gene coordinates to JSON."""
    gene = load_gene_from_annotation(args.annotation, args.gene_name)
    os.makedirs(args.output_dir, exist_ok=True)
    out = gene_cache_path(args.output_dir, args.gene_name)
    with open(out, "w") as fh:
        json.dump(gene, fh, indent=2)
    print(f"Gene '{args.gene_name}' ({gene['seqname']}:{gene['start']}-{gene['end']} "
          f"{gene['strand']}) cached → {out}", flush=True)


def cmd_merge(args):
    """Merge all per-RBP density CSVs into window × RBP matrices (raw and
    normalized), and write a single aggregated per-gene score per CLIP."""
    from functools import reduce

    pattern = os.path.join(args.output_dir, f"*_w{args.window_size}.csv")
    exclude_prefix = f"{args.gene_name}_density_matrix_w{args.window_size}"
    files = sorted(
        fp for fp in glob.glob(pattern)
        if not os.path.basename(fp).startswith(exclude_prefix)
    )
    if not files:
        sys.exit(f"ERROR: no *_w{args.window_size}.csv files found in {args.output_dir}")

    norm_totals = {}
    for value_col, dtype, suffix in [
        ("n_sites",      int,   ""),
        ("n_sites_norm", float, "_norm"),
    ]:
        dfs = []
        for fp in files:
            clip_name = os.path.basename(fp).replace(f"_w{args.window_size}.csv", "")
            df = pd.read_csv(fp)[["window_start", "window_end", value_col]]
            dfs.append(df.rename(columns={value_col: clip_name}))

        matrix = reduce(
            lambda left, right: pd.merge(left, right, on=["window_start", "window_end"], how="outer"),
            dfs,
        )
        matrix = matrix.sort_values("window_start").reset_index(drop=True)
        rbp_cols = matrix.columns[2:]
        matrix[rbp_cols] = matrix[rbp_cols].fillna(0).astype(dtype)

        out = os.path.join(args.output_dir,
                           f"{args.gene_name}_density_matrix_w{args.window_size}{suffix}.csv")
        matrix.to_csv(out, index=False)
        print(f"Merged {len(files)} files → {out}")
        print(f"Matrix: {len(matrix)} windows × {len(rbp_cols)} RBPs")

        if value_col == "n_sites_norm":
            norm_totals = matrix[rbp_cols].sum().to_dict()

    # Aggregated per-gene score: one normalized cross-link density value per
    # CLIP, summed over all windows. Written to the gene folder directly
    # (the parent of --output-dir), not the per-window/CLIP output dir.
    summary = pd.DataFrame(
        sorted(norm_totals.items()), columns=["clip", "n_cross_links"]
    )
    gene_dir = os.path.dirname(os.path.normpath(args.output_dir))
    summary_out = os.path.join(gene_dir, f"{args.gene_name}_cross_links.csv")
    summary.to_csv(summary_out, index=False)
    print(f"Aggregated per-gene cross-link scores → {summary_out}")


def cmd_run(args):
    """Process one CLIP file (one array task) and write a density CSV."""
    task_id = args.task_id
    if task_id is None:
        env_id = os.environ.get("SLURM_ARRAY_TASK_ID")
        if env_id is None:
            sys.exit("ERROR: --task-id not given and $SLURM_ARRAY_TASK_ID is not set")
        task_id = int(env_id)

    files = get_all_clip_files(args.iclip_dir, args.eclip_dir)
    if task_id >= len(files):
        sys.exit(f"ERROR: task_id {task_id} out of range (n_files={len(files)})")
    filepath, source = files[task_id]
    clip_name = clip_display_name(filepath, source)

    gene = load_gene_cached(args.output_dir, args.gene_name, args.annotation)

    print(f"[{task_id}] Processing: {filepath}", flush=True)
    sites = load_sites(filepath)
    print(f"  → {len(sites):,} sites loaded", flush=True)
    density = compute_density(gene, sites, args.window_size)

    os.makedirs(args.output_dir, exist_ok=True)
    out = os.path.join(args.output_dir, f"{clip_name}_w{args.window_size}.csv")
    density.to_csv(out, index=False)
    print(f"[{task_id}] {clip_name}: {len(density)} windows → {out}", flush=True)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Count genomic sites per window along a target gene.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # single
    ps = sub.add_parser("single",
        help="Process one BED/GFF sites file → one density CSV.")
    ps.add_argument("--gene-name",   required=True,
                    help="Gene name to look up in the annotation (e.g. NEAT1).")
    ps.add_argument("--annotation",  required=True,
                    help="Path to the full hg38 GTF annotation file.")
    ps.add_argument("--sites-file",  required=True,
                    help="BED or GFF/GTF sites file (format auto-detected by extension).")
    ps.add_argument("--window-size", type=int, default=200,
                    help="Window size in nucleotides.")
    ps.add_argument("--output",      required=True,
                    help="Path for the output CSV file.")

    # count
    pc = sub.add_parser("count", help="Print number of CLIP files.")
    pc.add_argument("--iclip-dir", default=None,
                    help="Directory containing iCLIP .gff files.")
    pc.add_argument("--eclip-dir", default=None,
                    help="Parent eCLIP directory with HepG2/ and K652/ subdirs.")

    # prepare
    pp = sub.add_parser("prepare",
        help="Cache target gene coordinates from GTF annotation.")
    pp.add_argument("--gene-name",  required=True)
    pp.add_argument("--annotation", required=True,
                    help="Path to the full hg38 GTF annotation file.")
    pp.add_argument("--output-dir", required=True)

    # run (array task)
    pr = sub.add_parser("run", help="Process one CLIP file (array task).")
    pr.add_argument("--gene-name",   required=True)
    pr.add_argument("--iclip-dir",   default=None)
    pr.add_argument("--eclip-dir",   default=None)
    pr.add_argument("--output-dir",  required=True)
    pr.add_argument("--window-size", type=int, default=200)
    pr.add_argument("--annotation",  default=None,
                    help="Fallback annotation if gene cache is missing.")
    pr.add_argument("--task-id",     type=int, default=None,
                    help="File index. Defaults to $SLURM_ARRAY_TASK_ID.")

    # merge
    pm = sub.add_parser("merge",
        help="Merge per-RBP density CSVs into a window × RBP matrix.")
    pm.add_argument("--gene-name",   required=True)
    pm.add_argument("--output-dir",  required=True)
    pm.add_argument("--window-size", type=int, default=200)

    args = parser.parse_args()
    {"single": cmd_single, "count": cmd_count,
     "prepare": cmd_prepare, "run": cmd_run,
     "merge": cmd_merge}[args.cmd](args)


if __name__ == "__main__":
    main()
