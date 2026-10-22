#!/usr/bin/env python3
"""Welch profiles (and Fisher plots) of the integrated CLIP / RNA-modification results of one folder.

Usage: python clip_analysis.py CONFIG --layout {time_course,knockdown,treatment_grid} [--welch-only]

CONFIG is an analysis config of pipeline/datasets/<dataset>/configs/ (path_to_processed_RBPs, output_directory,
treatments, palettes, ...). The layout depends on the design of the dataset:
- time_course (MN_diff): profiles per day, genotype and fraction (utrals.welch.time_course);
- knockdown (SFPQ_TDP_FUS_KD): rows = fractions, one line per treatment (utrals.welch.knockdown);
- treatment_grid (VCP-KI-Ml240): rows = fractions, one line per genotype x drug (utrals.welch.treatment_grid).
--welch-only: Welch profiles only, no by-region profiles or Fisher outputs.
"""
import argparse
import warnings

import matplotlib

from utrals.welch import knockdown, time_course, treatment_grid
from utrals.welch.common import load_config

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
warnings.filterwarnings('ignore')

LAYOUTS = {"time_course": time_course, "knockdown": knockdown, "treatment_grid": treatment_grid}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('config', help='analysis config (YAML)')
    parser.add_argument('--layout', required=True, choices=sorted(LAYOUTS), help='figure layout of the dataset')
    parser.add_argument('--welch-only', action='store_true',
                        help='only the Welch profiles: skip the by-region profiles and all Fisher outputs')
    args = parser.parse_args()

    print(f"Loading configuration from {args.config}...")
    config = load_config(args.config)
    LAYOUTS[args.layout].run(config, welch_only=args.welch_only)
    print(f"\nCLIP analysis completed. Output saved to: {config['output_directory']}")


if __name__ == '__main__':
    main()
