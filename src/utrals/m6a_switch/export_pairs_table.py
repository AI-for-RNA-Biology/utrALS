#!/usr/bin/env python
"""
One row per pair of table_of_pairs.csv (same order), with proximal_id, distal_id, gene and, for each
region x switch definition, whether the pair's isoform has >= 1 switch-competent m6A site in that
region's window (n_switch_competent > 0 in the region's per_transcript.csv):

  R1 = [-350:0] upstream of the proximal PAS (proximal isoform)   -> regions/proximal_m350_0
  R2 = [0:200] downstream of the proximal PAS (proximal isoform)  -> regions/proximal_0_200
  R3 = [-350:0] upstream of the distal PAS (distal isoform)       -> regions/distal_m350_0
  UUU = partner U in a run of >= 3 U's; GU = partner U in a UGU/UUG/GUG/GUU 3-mer

Output: <root>/comparison/table_of_pairs_m6a_switch.csv (+ .params.json)
Usage: python export_pairs_switch_table.py [--root <regions dir>]
"""

import argparse
import os
import sys

import pandas as pd

from .competence import M6A_RESULTS_ROOT, TABLE_OF_PAIRS, write_params

DEFAULT_ROOT = os.path.join(M6A_RESULTS_ROOT, "regions")
REGIONS = {"R1": ("proximal_m350_0", "proximal_id"), "R2": ("proximal_0_200", "proximal_id"),
           "R3": ("distal_m350_0", "distal_id")}
DEFINITIONS = {"UUU": "switch_competence_u_tract", "GU": "switch_competence_gu_motif"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    args = ap.parse_args(argv)

    out = pd.read_csv(TABLE_OF_PAIRS, usecols=["proximal_id", "distal_id", "gene"])[["proximal_id", "distal_id", "gene"]]
    for r, (sub, id_col) in REGIONS.items():
        for name, ddir in DEFINITIONS.items():
            per_tx = pd.read_csv(os.path.join(args.root, sub, ddir, "per_transcript.csv"), usecols=[id_col, "n_switch_competent"])
            has_switch = per_tx.set_index(id_col)["n_switch_competent"].gt(0)
            col = f"is {name}-m6a-switch in {r}"
            out[col] = out[id_col].map(has_switch)
            if out[col].isna().any():
                sys.exit(f"ERROR: {out[col].isna().sum()} pairs have no {r} {name} result")
            out[col] = out[col].astype(bool)

    out_dir = os.path.join(args.root, "comparison")
    path = os.path.join(out_dir, "table_of_pairs_m6a_switch.csv")
    out.to_csv(path, index=False)
    print(f"Wrote {path} ({len(out)} pairs)")
    print(out.iloc[:, 3:].sum().to_string())
    write_params(out_dir, run_params={"step": "export_pairs_switch_table", "regions": REGIONS,
                                      "definition": "n_switch_competent > 0"},
                 script_paths=[os.path.abspath(__file__)],
                 input_paths={"table_of_pairs": TABLE_OF_PAIRS, "root": os.path.abspath(args.root)},
                 filename="table_of_pairs_m6a_switch.csv.params.json")


if __name__ == "__main__":
    main()
