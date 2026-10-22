"""Enrichr gene-set libraries, downloaded once and then read from cache/Enrichr/."""
import json
from datetime import date

from .paths import CACHE_DIR

# ── Cached downloads (cache/, git-tracked) ─────────────────────────────────────
# Gene-set libraries change between releases or when a server is updated; the figure notebooks read them from a
# copy saved in cache/ the first time, so that a later run gives exactly the same figures.


def get_enrichr_library(name, organism="human"):
    """Enrichr gene-set library `name` ({term: [genes]}), downloaded once with gseapy.get_library and then read from
    cache/Enrichr/<name>_<organism>.json.gz (with the download date)."""
    import gzip
    path = CACHE_DIR / "Enrichr" / f"{name}_{organism}.json.gz"
    if path.exists():
        with gzip.open(path, "rt") as f:
            return json.load(f)["library"]
    import gseapy
    library = gseapy.get_library(name=name, organism=organism)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        json.dump({"library_name": name, "organism": organism, "download_date": date.today().isoformat(),
                   "library": library}, f)
    return library
