#!/usr/bin/env python
"""
m6A-switch competence and CLIP preference along whole lncRNAs (MALAT1, NEAT1_2, NORAD, XIST, TUG1).

1. --prepare
   - sequence: the gene's Ensembl-canonical transcript in GENCODE v44, as MATURE RNA (exons joined,
     introns removed; for single-exon genes this is simply the genomic sequence). NEAT1's canonical
     transcript is NEAT1-202 = NEAT1_2. Every transcript nucleotide is mapped to its genomic position.
   - m6A sites: every RMBase site on an exonic A of that transcript (gene strand).
   - each site is folded in ONE window centred on it in transcript coordinates (+/- FLANK nt, cut at
     the transcript ends; windows can span exon junctions), base pairs <= MAX_BP_SPAN nt apart
     (Lange et al. 2012 W = 200, L = 150) -- with
     fold_bpp(); partner j = most probable partner; calls:
       U-rich  : p(i,j) >= 0.5, j is a U in a run of >= 3 U's (within the window)
       UG-rich : p(i,j) >= 0.5, j is a U inside a UGU/UUG/GUG/GUU 3-mer
   - element catalogs over the whole transcript: every U-run >= 3 nt and every UGU/UUG/GUG/GUU
     occurrence (overlapping 3-mers kept separately, as in the 3'UTR analysis). An element is
     "switch" if it contains the partner U of a U-rich (resp. UG-rich) competent site, "other"
     otherwise. Elements are stored as genomic segments (one per exon they touch) with an
     element_id, so crosslinks are counted on the genome.
   -> genes.csv, per_site.csv, u_tract_catalog.csv, motif_catalog.csv, datasets.csv
1b. --features -> transcript_features.csv: exons, CDS and 5'/3' UTRs of each canonical transcript in
   transcript coordinates (1-based, inclusive), for the gene tracks of the figures (only TUG1-217
   has CDS/UTR annotation; the others are annotated as non-coding).
2. --chunk K --n-chunks N : every CLIP dataset in the chunk is loaded once and tested on each
   gene x element type: crosslinked POSITIONS (presence of crosslink evidence; read counts are not
   used, since the CLIP data come from other cell lines) in switch elements / crosslinked positions
   in all elements of that type, one-sided binomial vs. p0 = switch share of the element type
   (U-tract nucleotides; motif occurrences) -- the same test as the 3'UTR CLIP scripts.
   -> parts/part_K.csv
3. --merge -> clip_enrichment_all_datasets.csv (fold_change = observed / expected; BH q, rank and
   percentile of FC within gene x element) and rbp_summary.csv (per RBP: best dataset per gene x element).
   A dataset is TESTED on a gene x element only with >= MIN_CROSSLINKS crosslinked positions in
   that element type (column `tested`); BH correction, ranks and percentiles are computed over the
   tested datasets only, so lowly covered lncRNAs (where most datasets have a handful of
   crosslinks) don't produce significant calls from 2-3 crosslinks. Untested rows keep their raw
   counts and p-value but have no q.
   Also writes fc_q_matrix_by_rbp.csv (rows = RBPs, 2 columns per gene x element: "<gene>-<UUU|GU> FC" and
   "<gene>-<UUU|GU> q"; per RBP the dataset with the smallest q among its tested datasets, ties -> higher FC;
   empty = not tested) and fc_q_matrix_by_dataset.csv (same, one row per CLIP dataset).

Not included: SLERT, SPA1/SPA2, PNCTR (not GENCODE genes); the repeat-derived HSATIII, IGS lncRNAs
and TERRA (not resolvable in GRCh38 / with uniquely mapped CLIP reads); MIAT and FIRRE (essentially
no CLIP coverage in the K562/HepG2/HeLa/HEK293 datasets: no / almost no dataset reaches the
minimum number of crosslinked positions).

Output dir: M6A_RESULTS_ROOT/lncRNA_switch/
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import pysam
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from .competence import (
    GENOME_FASTA, M6A_BED, M6A_RESULTS_ROOT, P_THRESHOLD, U_TRACT_MIN,
    fold_bpp, has_u_tract, revcomp, write_params,
)
from .competence_gu_motif import MOTIF_LEN, MOTIFS, motif_at_position
from .clip_enrichment_u_tract import load_crosslinks, test_rbp
from .clip_enrichment_all_datasets import list_datasets
from ..paths import ANNOT_DIR  # data locations, see paths.env

GTF = f"{ANNOT_DIR}/GENCODE/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"
GENES = ["NEAT1", "MALAT1", "NORAD", "XIST", "TUG1"]
FLANK, MAX_BP_SPAN = 100, 150
OUT_DIR = os.path.join(M6A_RESULTS_ROOT, "lncRNA_switch")
ELEMENTS = {"U-tracts": ("u_tract_catalog.csv", "nt"), "GU motifs": ("motif_catalog.csv", "count")}
MIN_CROSSLINKS = 20   # minimum crosslinked positions in the element type for a dataset to be tested


# ---------------------------------------------------------------------------
# prepare
# ---------------------------------------------------------------------------
def canonical_transcripts(genes):
    """{gene: dict(transcript_id, transcript_name, chrom, strand, exons=[(start0, end), ...] 5'->3')}"""
    tx, exons = {}, {}
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t", 8)
            if x[2] not in ("transcript", "exon") or "_" in x[0] or "Ensembl_canonical" not in x[8]:
                continue
            g = re.search(r'gene_name "([^"]+)"', x[8]).group(1)
            if g not in genes:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
            if x[2] == "transcript":
                tx[g] = {"transcript_id": tid, "transcript_name": re.search(r'transcript_name "([^"]+)"', x[8]).group(1),
                         "chrom": x[0], "strand": x[6]}
            else:
                exons.setdefault(g, []).append((int(x[3]) - 1, int(x[4])))
    missing = set(genes) - set(tx)
    if missing:
        sys.exit(f"no Ensembl_canonical transcript for {missing} in {GTF}")
    for g in tx:
        tx[g]["exons"] = sorted(exons[g], reverse=tx[g]["strand"] == "-")
    return tx


def rna_seq(fasta, chrom, lo, hi, strand):
    s = fasta.fetch(chrom, lo, hi).upper()
    return (revcomp(s) if strand == "-" else s).replace("T", "U")


def mature(fasta, t):
    """-> (RNA sequence 5'->3', array: transcript index (0-based) -> genomic position (0-based))"""
    seq, gpos = [], []
    for lo, hi in t["exons"]:
        seq.append(rna_seq(fasta, t["chrom"], lo, hi, t["strand"]))
        gpos.append(np.arange(lo, hi) if t["strand"] == "+" else np.arange(hi - 1, lo - 1, -1))
    return "".join(seq), np.concatenate(gpos)


def fold_site(seq, k0):
    """Fold the window centred on transcript index k0 (0-based); j returned as a transcript index."""
    w_lo, w_hi = max(0, k0 - FLANK), min(len(seq), k0 + FLANK + 1)
    win = seq[w_lo:w_hi]
    i = k0 - w_lo + 1                                              # 1-based in the window
    best_j, best_p = fold_bpp(win, max_bp_span=MAX_BP_SPAN)
    j, p = int(best_j[i]), float(best_p[i])
    n, s1 = len(win), [""] + list(win)
    paired_u = j > 0 and p >= P_THRESHOLD and s1[j] == "U"
    motifs = motif_at_position(s1, j, n) if paired_u else set()
    return {"window_start_tx": w_lo, "window_end_tx": w_hi, "p_ij": p,
            "j_tx": (w_lo + j - 1) if j > 0 else np.nan,
            "u_rich": bool(paired_u and has_u_tract(s1, j, n)), "ug_rich": bool(paired_u and motifs),
            "ug_motifs": ",".join(sorted(motifs))}


def genomic_segments(gpos, k_start, k_end):
    """Transcript index span [k_start, k_end) -> list of genomic half-open segments (split at junctions)."""
    g = np.sort(gpos[k_start:k_end])
    breaks = np.where(np.diff(g) != 1)[0] + 1
    return [(int(part[0]), int(part[-1]) + 1) for part in np.split(g, breaks)]


def catalogs(gene, t, seq, gpos, sites):
    u_partners = set(sites.loc[sites["u_rich"], "j_tx"].astype(int))
    gu_partners = set(sites.loc[sites["ug_rich"], "j_tx"].astype(int))

    def rows(spans, partners, extra):
        out = []
        for eid, (k_s, k_e) in enumerate(spans):
            cat = "switch" if any(k_s <= j < k_e for j in partners) else "other"
            for s, e in genomic_segments(gpos, k_s, k_e):
                out.append({"gene": gene, "element_id": eid, "chrom": t["chrom"], "strand": t["strand"],
                            "start": s, "end": e, "tx_start": k_s, "tx_end": k_e, "category": cat,
                            **extra(k_s, k_e)})
        return pd.DataFrame(out)

    u_spans = [(m.start(), m.end()) for m in re.finditer(rf"U{{{U_TRACT_MIN},}}", seq)]
    gu_spans = [(k, k + MOTIF_LEN) for k in range(len(seq) - MOTIF_LEN + 1) if seq[k:k + MOTIF_LEN] in MOTIFS]
    return (rows(u_spans, u_partners, lambda a, b: {"length": b - a}),
            rows(gu_spans, gu_partners, lambda a, b: {"motif": seq[a:b]}))


def transcript_features(out_dir):
    """Exons, CDS and UTRs of each canonical transcript, mapped to mature-transcript coordinates."""
    fasta = pysam.FastaFile(GENOME_FASTA)
    txs = canonical_transcripts(GENES)
    feats = {}
    tid2gene = {t["transcript_id"]: g for g, t in txs.items()}
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t", 8)
            if x[2] not in ("CDS", "UTR") or "_" in x[0]:
                continue
            m = re.search(r'transcript_id "([^"]+)"', x[8])
            if m and m.group(1) in tid2gene:
                feats.setdefault(tid2gene[m.group(1)], []).append((x[2], int(x[3]) - 1, int(x[4])))
    rows = []
    for gene in GENES:
        t = txs[gene]
        _, gpos = mature(fasta, t)
        k_of = {int(g): k for k, g in enumerate(gpos)}

        def tx_span(lo, hi):
            ks = [k_of[g] for g in range(lo, hi) if g in k_of]
            return min(ks) + 1, max(ks) + 1

        for n, (lo, hi) in enumerate(t["exons"], start=1):
            a, b = tx_span(lo, hi)
            rows.append({"gene": gene, "feature": "exon", "number": n, "tx_start": a, "tx_end": b})
        cds = [tx_span(lo, hi) for kind, lo, hi in feats.get(gene, []) if kind == "CDS"]
        if cds:
            c0, c1 = min(a for a, _ in cds), max(b for _, b in cds)
            rows.append({"gene": gene, "feature": "CDS", "number": np.nan, "tx_start": c0, "tx_end": c1})
            for kind, lo, hi in feats[gene]:
                if kind == "UTR":
                    a, b = tx_span(lo, hi)
                    rows.append({"gene": gene, "feature": "5'UTR" if b < c0 else "3'UTR", "number": np.nan,
                                 "tx_start": a, "tx_end": b})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, "transcript_features.csv"), index=False)
    print(df.to_string(index=False))


def prepare(out_dir):
    fasta = pysam.FastaFile(GENOME_FASTA)
    txs = canonical_transcripts(GENES)
    bed = pd.read_csv(M6A_BED, sep="\t", header=None, usecols=[0, 1, 2, 3, 5],
                      names=["chrom", "start", "end", "site_id", "strand"])
    bed = bed[bed["end"] - bed["start"] == 1]
    per_site, u_cat, gu_cat, genes_rows = [], [], [], []
    for gene in GENES:
        t = txs[gene]
        seq, gpos = mature(fasta, t)
        k_of = {int(g): k for k, g in enumerate(gpos)}
        b = bed[(bed["chrom"] == t["chrom"]) & (bed["strand"] == t["strand"])
                & bed["start"].isin(k_of.keys())].drop_duplicates("start")
        rows = []
        for r in b.itertuples(index=False):
            k0 = k_of[r.start]
            if seq[k0] != "A":
                continue
            f = fold_site(seq, k0)
            rows.append({"gene": gene, "site_id": r.site_id, "chrom": t["chrom"], "pos": r.start, "strand": t["strand"],
                         "pos_in_gene": k0 + 1, **f,
                         "j_genomic": int(gpos[int(f["j_tx"])]) if not np.isnan(f["j_tx"]) else np.nan})
        sites = pd.DataFrame(rows)
        u, gu = catalogs(gene, t, seq, gpos, sites)
        n_u, n_gu = u["element_id"].nunique(), gu["element_id"].nunique()
        sw_u = u.loc[u["category"] == "switch", "element_id"].nunique()
        sw_gu = gu.loc[gu["category"] == "switch", "element_id"].nunique()
        print(f"{gene} {t['transcript_name']} {len(seq)} nt, {len(t['exons'])} exon(s): {len(b)} RMBase sites in exons, "
              f"{len(sites)} on an A; U-rich {sites['u_rich'].sum()}, UG-rich {sites['ug_rich'].sum()}; "
              f"switch U-tracts {sw_u}/{n_u}, switch GU motifs {sw_gu}/{n_gu}", flush=True)
        per_site.append(sites); u_cat.append(u); gu_cat.append(gu)
        genes_rows.append({"gene": gene, "transcript_id": t["transcript_id"], "transcript_name": t["transcript_name"],
                           "chrom": t["chrom"], "start": int(gpos.min()), "end": int(gpos.max()) + 1,
                           "strand": t["strand"], "n_exons": len(t["exons"]), "length": len(seq),
                           "n_rmbase_sites": len(b), "n_sites_on_A": len(sites),
                           "n_u_rich": int(sites["u_rich"].sum()), "n_ug_rich": int(sites["ug_rich"].sum()),
                           "n_u_tracts": n_u, "n_switch_u_tracts": sw_u, "n_motifs": n_gu, "n_switch_motifs": sw_gu})
    pd.concat(per_site).to_csv(os.path.join(out_dir, "per_site.csv"), index=False)
    pd.concat(u_cat).to_csv(os.path.join(out_dir, "u_tract_catalog.csv"), index=False)
    pd.concat(gu_cat).to_csv(os.path.join(out_dir, "motif_catalog.csv"), index=False)
    genes_df = pd.DataFrame(genes_rows)
    genes_df.to_csv(os.path.join(out_dir, "genes.csv"), index=False)
    print(genes_df.to_string(index=False))
    list_datasets().to_csv(os.path.join(out_dir, "datasets.csv"), index=False)


# ---------------------------------------------------------------------------
# CLIP
# ---------------------------------------------------------------------------
def load_targets(out_dir):
    targets = []
    for element, (cat_file, weight) in ELEMENTS.items():
        cat = pd.read_csv(os.path.join(out_dir, cat_file))
        for gene, c in cat.groupby("gene", sort=False):
            c = c.sort_values("start")
            if weight == "nt":   # nucleotides (segments of one element are disjoint)
                w = c["end"] - c["start"]
                p0 = w[c["category"] == "switch"].sum() / w.sum()
            else:                # occurrences: count elements, not their genomic segments
                e = c.drop_duplicates("element_id")
                p0 = (e["category"] == "switch").sum() / len(e)
            grouped = {(chrom, strand): (g["start"].values, g["end"].values, g["category"].values)
                       for (chrom, strand), g in c.groupby(["chrom", "strand"])}
            targets.append((gene, element, p0, grouped))
    return targets


def run_chunk(out_dir, k, n_chunks):
    ds = pd.read_csv(os.path.join(out_dir, "datasets.csv"))
    order = ds.sort_values("size_bytes", ascending=False).index
    mine = ds.loc[[i for j, i in enumerate(order) if j % n_chunks == k]]
    targets = load_targets(out_dir)
    rows = []
    for d in mine.itertuples(index=False):
        print(f"[{k}] {d.dataset}", flush=True)
        clip_by_key = load_crosslinks(d.path)
        for gene, element, p0, grouped in targets:
            n, kk = test_rbp(clip_by_key, grouped)
            obs = kk / n if n else np.nan
            rows.append({"dataset": d.dataset, "RBP": d.RBP, "assay": d.assay, "gene": gene, "element": element,
                         "n_crosslinks_in_elements": n, "n_in_switch_elements": kk,
                         "observed_rate": obs, "expected_rate": p0,
                         "odds_ratio": (obs / (1 - obs)) / (p0 / (1 - p0)) if 0 < obs < 1 else np.nan,
                         "p_value": binomtest(kk, n, p0, alternative="greater").pvalue if n else np.nan})
        del clip_by_key
    os.makedirs(os.path.join(out_dir, "parts"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "parts", f"part_{k:03d}.csv"), index=False)


def merge(out_dir, n_chunks):
    parts = sorted(glob.glob(os.path.join(out_dir, "parts", "part_*.csv")))
    if len(parts) != n_chunks:
        sys.exit(f"ERROR: {len(parts)}/{n_chunks} parts present")
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    ds = pd.read_csv(os.path.join(out_dir, "datasets.csv"))
    assert set(df["dataset"]) == set(ds["dataset"]) and len(df) == len(ds) * len(GENES) * len(ELEMENTS)
    with np.errstate(divide="ignore"):
        df["log2_odds_ratio"] = np.log2(df["odds_ratio"])
    df["fold_change"] = df["observed_rate"] / df["expected_rate"]
    df["tested"] = df["n_crosslinks_in_elements"] >= MIN_CROSSLINKS
    for _, g in df.groupby(["gene", "element"]):
        g = g[g["tested"]]
        if g.empty:
            continue
        df.loc[g.index, "q_value_bh"] = multipletests(g["p_value"], method="fdr_bh")[1]
        df.loc[g.index, "rank_fc"] = g["fold_change"].rank(ascending=False, method="min")
        df.loc[g.index, "percentile_fc"] = 100 * g["fold_change"].rank(pct=True)
    df = df.sort_values(["gene", "element", "tested", "rank_fc"], ascending=[True, True, False, True])
    df.to_csv(os.path.join(out_dir, "clip_enrichment_all_datasets.csv"), index=False)

    # per RBP and gene x element: its best dataset = highest FC among its significant (q < 0.05)
    # datasets, or among all its datasets if none is significant
    sig = df[df["tested"]].assign(sig=lambda d: d["q_value_bh"] < 0.05)
    best = (sig.sort_values(["sig", "fold_change"], ascending=[False, False])
            .groupby(["gene", "element", "RBP"], sort=False).first().reset_index())
    n_sig = sig.groupby(["gene", "element", "RBP"])["sig"].agg(n_datasets="size", n_datasets_significant="sum").reset_index()
    rbp = best.merge(n_sig, on=["gene", "element", "RBP"])[
        ["gene", "element", "RBP", "n_datasets", "n_datasets_significant", "dataset", "assay",
         "n_crosslinks_in_elements", "n_in_switch_elements", "observed_rate", "expected_rate",
         "fold_change", "p_value", "q_value_bh"]].rename(columns={"dataset": "best_dataset"})
    rbp.to_csv(os.path.join(out_dir, "rbp_summary.csv"), index=False)
    write_fc_q_matrices(df, out_dir)
    print(df.groupby(["gene", "element"]).agg(expected=("expected_rate", "first"),
                                               n_tested=("tested", "sum"),
                                               n_sig=("q_value_bh", lambda q: int((q < 0.05).sum())),
                                               median_fc=("fold_change", "median")).to_string())
    write_params(out_dir, run_params={"step": "lncrna_m6a_switch_clip", "genes": GENES,
                                      "sequence": "Ensembl-canonical transcript, mature (exons joined)", "flank": FLANK,
                                      "max_bp_span": MAX_BP_SPAN, "p_threshold": P_THRESHOLD,
                                      "u_tract_min": U_TRACT_MIN, "motifs": MOTIFS, "n_chunks": n_chunks,
                                      "test": "crosslinked positions, one-sided binomial",
                                      "min_crosslinks_tested": MIN_CROSSLINKS},
                 script_paths=[os.path.abspath(__file__)] + [
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), s)
                     for s in ("competence.py", "competence_gu_motif.py",
                               "clip_enrichment_u_tract.py", "clip_enrichment_all_datasets.py")],
                 input_paths={"gtf": GTF, "m6a_bed": M6A_BED, "genome_fasta": GENOME_FASTA})


def write_fc_q_matrices(df, out_dir):
    """Wide FC / q tables: one row per RBP (its most significant tested dataset per column) and per dataset."""
    short = {"U-tracts": "UUU", "GU motifs": "GU"}
    cols = [f"{g}-{short[e]}" for g in GENES for e in ELEMENTS]
    t = df[df["tested"]].assign(col=lambda d: d["gene"] + "-" + d["element"].map(short),
                                FC=lambda d: d["observed_rate"] / d["expected_rate"], q=lambda d: d["q_value_bh"])
    best = (t.sort_values(["q", "FC"], ascending=[True, False]).groupby(["RBP", "col"], sort=False).first().reset_index())
    for name, data, idx in [("fc_q_matrix_by_rbp.csv", best, "RBP"), ("fc_q_matrix_by_dataset.csv", t, "dataset")]:
        wide = data.pivot(index=idx, columns="col", values=["FC", "q"])
        out = pd.DataFrame(index=sorted(df[idx].unique()))
        out.index.name = idx
        for c in cols:
            for v in ("FC", "q"):
                out[f"{c} {v}"] = wide[(v, c)] if (v, c) in wide.columns else np.nan
        if idx == "dataset":
            out.insert(0, "RBP", df.drop_duplicates("dataset").set_index("dataset")["RBP"].reindex(out.index))
        else:
            best_ds = best.pivot(index="RBP", columns="col", values="dataset")
            for c in cols:   # which dataset each RBP's values come from
                out[f"{c} dataset"] = best_ds[c] if c in best_ds.columns else np.nan
        out.to_csv(os.path.join(out_dir, name))
        print(f"Wrote {os.path.join(out_dir, name)} ({len(out)} rows)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--features", action="store_true")
    ap.add_argument("--chunk", type=int)
    ap.add_argument("--n-chunks", type=int, default=40)
    ap.add_argument("--merge", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    if args.prepare:
        prepare(args.out_dir)
    elif args.features:
        transcript_features(args.out_dir)
    elif args.chunk is not None:
        run_chunk(args.out_dir, args.chunk, args.n_chunks)
    elif args.merge:
        merge(args.out_dir, args.n_chunks)
    else:
        ap.error("one of --prepare, --features, --chunk, --merge")


if __name__ == "__main__":
    main()
