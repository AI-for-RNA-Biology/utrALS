"""Figures and tables saved by the notebooks, each listed with its legend in Figures/README.md /
Tables/README.md; reproducible random draws (fix_random_seeds)."""
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .paths import FIGURES_DIR, TABLES_DIR

# ── Figure and table catalogues (Figures/README.md, Tables/README.md) ───────────
# Every file saved with save_figure() / register_figure() is listed with its legend
# in Figures/README.md, and every table saved with save_table() / register_table()
# in Tables/README.md, so they can be used without reading the notebooks.
# The legends are kept in Figures/figure_legends.json and Tables/table_legends.json,
# and each README is rebuilt from its json at each save; entries whose file no
# longer exists are dropped.
import fcntl
import json
from datetime import date
from pathlib import Path

FIGURE_LEGENDS_JSON = FIGURES_DIR / "figure_legends.json"
FIGURES_README = FIGURES_DIR / "README.md"

TABLE_LEGENDS_JSON = TABLES_DIR / "table_legends.json"
TABLES_README = TABLES_DIR / "README.md"

_SHARED_CONVENTIONS = (
    "Shared conventions: CTRL = control lines (\"Wildtype\" / \"WT\" in the tables); VCP = "
    "VCP-R155C and VCP-R191Q lines pooled; PUD = proximal usage (proximal / (proximal + distal) "
    "reads), ΔPUD > 0 = 3'UTR shortening; significant 3'UTR shift = adjusted Fisher p < 0.01, "
    "|ΔPUD| ≥ 0.15 and |ΔRUD| > 1, unless stated otherwise. Welch score = −log10(Welch p) × ΔPUD "
    "of 3'UTR pairs with vs without a CLIP (or RNA-modification) site in a 50-nt bin around the "
    "poly(A) site. Regions: R1 = [−350, −50] nt upstream of the proximal PAS (pPAS), "
    "R2 = [0, +150] nt downstream of the pPAS, R3 = [−350, −50] nt upstream of the distal PAS (dPAS)."
)
_FIGURES_INTRO = (
    "# Figures",
    "Every figure saved by the notebooks in `notebooks/`, with its legend. "
    "This file is rebuilt automatically each time a notebook saves a file with "
    "`utrals.figures.save_figure` / `register_figure` (legends are stored in "
    "`figure_legends.json`): do not edit it by hand, edit the legend in the notebook instead. "
    "`<...>` in a file name stands for a set of files saved in a loop (one per gene, CLIP, ...). "
    "The tables behind the figures are in `../Tables/` (see `../Tables/README.md`).",
)
_TABLES_INTRO = (
    "# Tables",
    "Every table (.csv / .xlsx) saved by the notebooks in `notebooks/` (one folder per notebook), "
    "with a description of its content. Unlike `Figures/`, this folder is versioned in git. "
    "This file is rebuilt automatically each time a notebook saves a table with "
    "`utrals.figures.save_table` / `register_table` (descriptions are stored in "
    "`table_legends.json`): do not edit it by hand, edit the description in the notebook instead. "
    "`<...>` in a file name stands for a set of files saved in a loop (one per region, gene, ...). "
    "Numbers are written with 12 significant digits.",
)


def _entry_exists(root, key, entry):
    if entry.get("pattern"):
        # "<CLIP>" in a pattern stands for any file name part
        return any(root.glob(re.sub(r"<[^>]+>", "*", key)))
    return (root / key).exists()


def _write_readme(readme, intro, legends):
    title, text = intro
    lines = [title, "", text, "", _SHARED_CONVENTIONS, ""]
    by_figure = {}
    for key in sorted(legends):
        parts = Path(key).parts
        top = parts[0] if len(parts) > 1 else "(top level)"
        folder = str(Path(*parts[1:-1])) if len(parts) > 2 else ""
        by_figure.setdefault(top, {}).setdefault(folder, []).append(key)
    for top, folders in by_figure.items():
        lines += [f"## {top}", ""]
        for folder, keys in sorted(folders.items(), key=lambda kv: (kv[0] != "", kv[0])):
            if folder:
                lines += [f"### `{folder}/`", ""]
            for key in keys:
                entry = legends[key]
                lines.append(f"- **`{Path(key).name}`**: {entry['legend']} "
                             f"*(`{entry['notebook']}`, saved {entry['saved']})*")
            lines.append("")
    readme.write_text("\n".join(lines))


def _register(path, legend, notebook, pattern, root, legends_json, readme, intro, caller):
    path = Path(path).resolve()
    if not path.exists():
        print(f"{caller}: {path} does not exist, not listed in {readme}")
        return
    try:
        rel = path.relative_to(root)
    except ValueError:
        print(f"{caller}: {path} is not under {root}, not listed")
        return
    key = str(rel.parent / pattern) if pattern else str(rel)

    root.mkdir(exist_ok=True)
    # Notebooks running at the same time update the same json: lock it, and replace it atomically so that
    # a reader never sees a half-written file
    with open(root / f".{legends_json.name}.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        legends = json.loads(legends_json.read_text()) if legends_json.exists() else {}
        legends[key] = {"legend": " ".join(legend.split()), "notebook": notebook,
                        "pattern": bool(pattern), "saved": date.today().isoformat()}
        legends = {k: v for k, v in legends.items() if _entry_exists(root, k, v)}
        tmp = legends_json.with_name(f".{legends_json.name}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps(legends, indent=1, ensure_ascii=False, sort_keys=True))
        os.replace(tmp, legends_json)
        _write_readme(readme, intro, legends)


def register_figure(path, legend, notebook, pattern=None):
    """Record the legend of a file saved under Figures/ and rebuild Figures/README.md.

    `legend`: what the file shows (axes, colours, statistics), as for a paper figure legend.
    `notebook`: the notebook that saves it (e.g. "Figure_3UTR_dynamics_over_time.ipynb").
    `pattern`: for files saved in a loop (one per gene, CLIP, ...), the file name to list
        instead of each file, e.g. "<CLIP>.pdf"; it replaces the file name of `path`.
    """
    _register(path, legend, notebook, pattern, FIGURES_DIR, FIGURE_LEGENDS_JSON, FIGURES_README,
              _FIGURES_INTRO, "register_figure")


def save_figure(path, legend, notebook, pattern=None, fig=None, **savefig_kwargs):
    """Save `fig` (default: the current figure) to `path` and list it with its legend in
    Figures/README.md (see register_figure). Creates the folder; bbox_inches="tight" by default.

    In a notebook, bind the notebook name once:
        save_figure = functools.partial(figures.save_figure, notebook="Figure_3UTR_dynamics_over_time.ipynb")
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    savefig_kwargs.setdefault("bbox_inches", "tight")
    if str(path).endswith(".pdf"):
        # no creation date in the file, so that re-running a notebook gives byte-identical PDFs
        savefig_kwargs.setdefault("metadata", {"CreationDate": None})
    (fig or plt.gcf()).savefig(path, **savefig_kwargs)
    register_figure(path, legend, notebook, pattern=pattern)


def register_table(path, legend, notebook, pattern=None):
    """Record the description of a table saved under Tables/ and rebuild Tables/README.md.

    `legend`: what the table holds (one row per ..., the columns, how they are computed).
    `notebook`, `pattern`: as in register_figure (e.g. pattern="<region>.csv").
    """
    _register(path, legend, notebook, pattern, TABLES_DIR, TABLE_LEGENDS_JSON, TABLES_README,
              _TABLES_INTRO, "register_table")


def _round_floats(df, digits=12):
    """Copy of `df` with its float columns rounded to `digits` significant digits."""
    df = df.copy()
    for col in df.select_dtypes("float").columns:
        df[col] = df[col].map(lambda x: float(f"{x:.{digits}g}") if pd.notna(x) else x)
    return df


def save_table(table, path, legend, notebook, pattern=None, **to_kwargs):
    """Save a table under Tables/ and list it with its description in Tables/README.md
    (see register_table). Creates the folder.

    `table`: a DataFrame, or for an .xlsx file a {sheet name: DataFrame} dict (one sheet each).
    `to_kwargs` go to DataFrame.to_csv / to_excel (e.g. index=False). Floats are written with 12 significant digits.

    In a notebook, bind the notebook name once:
        save_table = functools.partial(figures.save_table, notebook="Figure_3UTR_dynamics_over_time.ipynb")
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    # Floats get 12 significant digits: their last digits can differ between CPU types, which would make
    # re-runs on another machine give a different table for the same values
    if str(path).endswith(".xlsx"):
        sheets = table if isinstance(table, dict) else {"Sheet1": table}
        with pd.ExcelWriter(path) as writer:
            for sheet_name, df in sheets.items():
                _round_floats(df).to_excel(writer, sheet_name=sheet_name, **to_kwargs)
    else:
        to_kwargs.setdefault("float_format", "%.12g")
        table.to_csv(path, **to_kwargs)
    register_table(path, legend, notebook, pattern=pattern)


# The 3'UTR regions of the regulator tables (R1-R3, see _SHARED_CONVENTIONS), as written in their file names
REGION_FILE_NAMES = {1: "region1_upstream_pPAS", 2: "region2_downstream_pPAS", 3: "region3_upstream_dPAS"}


def lengthening_regulators_csv(region):
    """Path of the regulators of nuclear 3'UTR lengthening (D22 -> D35, CTRL) of region 1, 2 or 3 with their
    APA / nuclear-retention class, written by notebooks/Figure_regulators_lengthening.ipynb and read by the
    Figure_m6A_* notebooks, scripts/ and the m6a_regions pipeline step."""
    return (TABLES_DIR / "Figure_regulators_lengthening" /
            f"regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_{REGION_FILE_NAMES[region]}_APA_vs_nuclear_retention_class.csv")





def fix_random_seeds(seed=0):
    """Make the random parts of the figures identical on every run; call it once, at the top of a notebook.

    Seeds NumPy's global random state (used, e.g., by seaborn's strip-plot jitter), and makes seaborn's bootstrap
    (the confidence intervals of barplot, pointplot, lineplot, regplot, ...) draw from one generator seeded with `seed`
    instead of a fresh, unseeded one at every call (seaborn's default, seed=None). A call that passes its own `seed`
    keeps it. The numbers drawn depend on the order of the plots, so a figure is reproduced exactly by running the
    notebook from the top.
    """
    import functools
    import sys as sys_
    import seaborn.algorithms as sns_algorithms
    np.random.seed(seed)
    unseeded = getattr(sns_algorithms.bootstrap, "_unseeded", sns_algorithms.bootstrap)
    rng = np.random.RandomState(seed)

    @functools.wraps(unseeded)
    def bootstrap(*args, **kwargs):
        if kwargs.get("seed") is None and kwargs.get("random_seed") is None:
            kwargs.pop("random_seed", None)
            kwargs["seed"] = rng
        return unseeded(*args, **kwargs)

    bootstrap._unseeded = unseeded
    # seaborn modules import the function by name (from .algorithms import bootstrap): replace it everywhere
    for name, module in list(sys_.modules.items()):
        current = getattr(module, "bootstrap", None)
        if name.startswith("seaborn") and (current is unseeded or hasattr(current, "_unseeded")):
            module.bootstrap = bootstrap
