"""Locations of the data outside this repository, as defined in paths.env (local, git-ignored; template:
paths.env.example), and of the
repository folders written or read by the notebooks (Figures/, Tables/, cache/).

An environment variable of the same name overrides a location, e.g. PROJECT_DIR=/my/scratch; the locations defined
from it then follow. Standard library only, so that every environment can import it:

    from utrals.paths import PROJECT_DIR, ANNOT_DIR
"""
import os
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]     # src/utrals/paths.py -> repository root
PATHS_ENV = REPO_ROOT / "paths.env"


def _load(path=PATHS_ENV):
    """{KEY: value} of paths.env, ${KEY} expanded, environment variables taking precedence."""
    values = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, raw = line.split("=", 1)
        values[key] = os.environ.get(key) or re.sub(r"\$\{(\w+)\}", lambda m: values[m.group(1)], raw)
    return values


if not PATHS_ENV.exists():
    raise FileNotFoundError(f"{PATHS_ENV} not found: copy paths.env.example to paths.env and set DATA_ROOT "
                            "(the data locations of your system)")
_PATHS = _load()
globals().update(_PATHS)
FIGURES_DIR = REPO_ROOT / "Figures"     # figures saved by the notebooks (figures.save_figure)
TABLES_DIR = REPO_ROOT / "Tables"       # tables saved by the notebooks (figures.save_table)
CACHE_DIR = REPO_ROOT / "cache"         # archived web downloads (STRING, Enrichr, QuickGO)

__all__ = ["REPO_ROOT", "FIGURES_DIR", "TABLES_DIR", "CACHE_DIR", *_PATHS]
