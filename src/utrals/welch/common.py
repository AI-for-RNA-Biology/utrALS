"""Shared parts of the Welch profile analyses: configuration, integrated CLIP results, the per-item save loop."""
import traceback
from pathlib import Path

import yaml
import matplotlib.pyplot as plt
from cliptools.analysis import load_and_combine_integrated_CLIP

from .fisher import plot_fisher_heatmap


def load_config(config_path):
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def setup_output_directories(base_dir, subdirs):
    """Create `base_dir`/<subdir> for every {key: subdir} and return {key: Path}."""
    paths = {key: Path(base_dir) / subdir for key, subdir in subdirs.items()}
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def load_integrated(config, **kwargs):
    """Proximal, distal and stop codon results ({CLIP: DataFrame}) of `config['path_to_processed_RBPs']`
    (`kwargs`: cols, RBP_list of load_and_combine_integrated_CLIP)."""
    path = config['path_to_processed_RBPs']
    print(f"\nLoading CLIP data from {path}...")
    results = [load_and_combine_integrated_CLIP(path, prefix=prefix, **kwargs)
               for prefix in ('proximal', 'distal', 'stop_codon')]
    for prefix, res in zip(('proximal', 'distal', 'stop codon'), results):
        print(f"  Loaded {len(res)} RBPs for {prefix} analysis")
    return results


def save(path, fig=None):
    """Save the current (or given) figure, cropped to its content and without a creation date (identical file on
    re-runs), and close it."""
    (fig or plt.gcf()).savefig(path, bbox_inches='tight',
                               metadata={'CreationDate': None} if str(path).endswith('.pdf') else None)
    plt.close('all')


def for_each(items, what, plot_and_save):
    """Call plot_and_save(item) for every item; a failure is reported (with its traceback) and skipped."""
    print(f"\nGenerating {what}...")
    for i, item in enumerate(items, 1):
        print(f"  [{i}/{len(items)}] Processing {item}...")
        try:
            plot_and_save(item)
        except Exception as e:
            print(f"    Warning: Failed to process {item}\n"
                  f"    Exception type: {type(e).__name__}\n"
                  f"    Exception message: {e}\n"
                  f"    Traceback:\n{traceback.format_exc()}")
            plt.close('all')


def fisher_heatmaps(config, results_by_isoform, output_dirs, title, merge_ptb=True):
    """One Fisher heatmap per pair of `config['pairs_for_heatmaps']` and isoform type
    (results_by_isoform: {'proximal': ..., 'distal': ...}; output_dirs: {'proximal': Path, 'distal': Path})."""
    for pairs_name in config.get('pairs_for_heatmaps') or []:
        for iso, df_results in results_by_isoform.items():
            for_each([pairs_name], f"Fisher heatmaps ({iso})", lambda pairs: plot_fisher_heatmap(
                pairs_name=pairs, df_results=df_results, title=title(pairs),
                top_n=config.get('heatmap_top_n', 50),
                nb_nucleotides_upstream=config.get('nb_nucleotides_upstream', -1000),
                length_region=config.get('length_region', 50),
                save_path=str(output_dirs[iso] / f"heatmap_fisher_along_{iso}_{pairs}.pdf"),
                merge_ptb=merge_ptb))
            plt.close('all')
