"""Query STRING's enrichment API directly instead of downloading TSVs from the website.

Reproduces the per-category files under results/3UTR_analysis/pathways/
(<name>_GO.tsv, <name>_COMPARTMENTS.tsv, <name>_Keyword.tsv, <name>_RCTM.tsv) from a
plain gene-list file (one gene symbol per line, e.g. results/3UTR_analysis/shifts/*.csv),
computing STRING's "strength" and "signal" columns the same way as
utrals.string_analysis.get_enrichment_df() (background_size=19566 is STRING's
whole-genome human background, confirmed by back-solving the existing downloaded
pathways/*_GO.tsv files).

Usage:
    python -m utrals.string_api <gene_list.csv> <output_dir> [--species 9606]
    python -m utrals.string_api --shifts-dir <shifts_dir> --output-dir <pathways_dir>
"""

import argparse
import datetime
import gzip
import hashlib
import json
import os

import numpy as np
import pandas as pd

STRING_API = "https://string-db.org/api/json"
# Release-pinned address of STRING v12.5: the address above follows the latest STRING release, whose results change
STRING_API_V12_5 = "https://version-12-5.string-db.org/api/json"
SPECIES = 9606
# Folder where string_api_post() archives the raw STRING responses (None = no archive; the notebooks set it to
# utrals.string_analysis.STRING_CACHE_DIR so that their figures are reproduced exactly on later runs)
CACHE_DIR = None
BACKGROUND_SIZE = 19566  # STRING human whole-genome background (STRING v12.0)

# STRING "category" values to merge/rename into the pathways/ file naming convention
CATEGORY_MAP = {
    "GO": ("Process",),
    "COMPARTMENTS": ("COMPARTMENTS",),
    "Keyword": ("Keyword",),
    "RCTM": ("RCTM",),
}


def string_api_post(url, data, cache_dir=None):
    """POST `data` to the STRING API `url` and return the parsed JSON response (None on an HTTP error).

    With `cache_dir`, the raw response is archived there (gzipped JSON, with the request, STRING version and date)
    under a hash of the exact request; the same request later reads that file and never contacts STRING, so results
    stay identical even if STRING changes or goes offline. Errors are not archived. Delete a file to query again.
    """
    key = hashlib.sha256(json.dumps({"url": url, "data": data}, sort_keys=True).encode()).hexdigest()[:20]
    path = os.path.join(cache_dir, f"{url.rstrip('/').rsplit('/', 1)[-1]}_{key}.json.gz") if cache_dir else None
    if path and os.path.exists(path):
        with gzip.open(path, "rt") as f:
            return json.load(f)["response"]
    import requests  # only needed to query STRING (absent from some environments that read the archive)
    r = requests.post(url, data=data)
    if r.status_code != 200:
        print(f"STRING error {r.status_code} ({url}): {r.text[:200]}")
        return None
    response = r.json()
    if path:
        api = url.rstrip("/").rsplit("/", 1)[0]
        version = requests.get(f"{api}/version").json()[0]["string_version"]
        os.makedirs(cache_dir, exist_ok=True)
        with gzip.open(path, "wt") as f:
            json.dump({"url": url, "request": data, "string_version": version,
                       "query_date": datetime.date.today().isoformat(), "response": response}, f)
    return response


def get_string_enrichment(genes, species=SPECIES):
    """STRING enrichment of `genes` (address STRING_API; responses archived in CACHE_DIR when it is set)."""
    data = string_api_post(f"{STRING_API}/enrichment",
                           {"identifiers": "\r".join(sorted(set(genes))), "species": species,
                            "caller_identity": "ALS_UTR_RBP"},
                           cache_dir=CACHE_DIR)
    if data is None:
        raise RuntimeError("STRING enrichment request failed")
    return pd.DataFrame(data) if data else pd.DataFrame()


def get_string_enrichment_cached(genes, cache_path, species=SPECIES, query_if_missing=True):
    """Same as get_string_enrichment(), but saved to / read from `cache_path` (gzipped JSON).

    The first call queries STRING and saves its raw response together with the gene list,
    STRING version, API address and query date; later calls read that file and never contact
    STRING, so the results stay identical even if STRING's servers change or go offline
    (STRING's results for the same genes changed between 2026-10-13 and 2026-10-20).
    Raises if the saved gene list differs from `genes`; delete the file to query STRING again.
    With query_if_missing=False, a missing file raises instead of silently querying STRING again
    (use it once the response a figure was made from has been saved).
    """
    genes_sorted = sorted(set(genes))
    if os.path.exists(cache_path):
        with gzip.open(cache_path, "rt") as f:
            cache = json.load(f)
        if cache["genes"] != genes_sorted:
            raise ValueError(f"{cache_path} was saved for a different gene list; delete it to query STRING again")
        return pd.DataFrame(cache["results"])
    if not query_if_missing:
        raise FileNotFoundError(f"{cache_path} not found; it holds the STRING response the figures were made from. "
                                "Restore it (git) rather than querying STRING again, whose results change over time")

    import requests
    version = requests.get(f"{STRING_API}/version").json()[0]["string_version"]
    df = get_string_enrichment(genes_sorted, species=species)
    os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
    with gzip.open(cache_path, "wt") as f:
        json.dump({
            "string_version": version,
            "string_api": STRING_API,
            "species": species,
            "query_date": datetime.date.today().isoformat(),
            "genes": genes_sorted,
            "results": df.to_dict(orient="records"),
        }, f)
    return df


def add_strength_and_signal(df, background_size=BACKGROUND_SIZE):
    """strength = log10(observed / expected), expected = n * K / N.
    signal = harmonic mean of strength and -log10(FDR).
    n = network size inferred from the union of preferredNames across returned terms.
    """
    all_prefs = set()
    for prefs in df["preferredNames"]:
        all_prefs.update(prefs)
    n = max(len(all_prefs), 1)

    k = df["number_of_genes"].astype(float)
    K = df["number_of_genes_in_background"].astype(float)
    expected = n * K / background_size

    df = df.copy()
    df["strength"] = np.log10((k / expected).clip(lower=1e-9))
    log_fdr = -np.log10(df["fdr"].clip(lower=1e-300))
    s = df["strength"].clip(lower=1e-9)
    lf = log_fdr.clip(lower=1e-9)
    df["signal"] = 2 * s * lf / (s + lf)
    return df


def to_pathways_format(df):
    """Rename/reorder columns to match the manually-downloaded pathways/*.tsv format."""
    out = pd.DataFrame({
        "term ID": df["term"],
        "term description": df["description"],
        "observed gene count": df["number_of_genes"].astype(int),
        "background gene count": df["number_of_genes_in_background"].astype(int),
        "strength": df["strength"].round(2),
        "signal": df["signal"].round(2),
        "false discovery rate": df["fdr"],
    })
    return out.sort_values("false discovery rate")


def query_and_write(genes, output_prefix, species=SPECIES):
    """Query STRING for `genes` and write one TSV per category next to `output_prefix`."""
    df = get_string_enrichment(genes, species=species)
    if df.empty:
        print(f"  no enrichment results for {output_prefix}")
        return
    df = add_strength_and_signal(df)

    for out_name, categories in CATEGORY_MAP.items():
        sub = df[df["category"].isin(categories)]
        if sub.empty:
            continue
        out_path = f"{output_prefix}_{out_name}.tsv"
        to_pathways_format(sub).to_csv(out_path, sep="\t", index=False)
        print(f"  wrote {out_path} ({len(sub)} terms)")


def read_gene_list(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gene_list", nargs="?", help="Path to a gene-list file (one symbol per line)")
    parser.add_argument("output_dir", nargs="?", help="Directory to write <name>_<CATEGORY>.tsv files into")
    parser.add_argument("--shifts-dir", help="Batch mode: process every *.csv gene list in this directory")
    parser.add_argument("--output-dir", dest="batch_output_dir", help="Batch mode: output directory")
    parser.add_argument("--species", type=int, default=SPECIES)
    args = parser.parse_args()

    if args.shifts_dir:
        out_dir = args.batch_output_dir or args.output_dir
        os.makedirs(out_dir, exist_ok=True)
        for fname in sorted(os.listdir(args.shifts_dir)):
            if not fname.endswith(".csv"):
                continue
            name = os.path.splitext(fname)[0]
            genes = read_gene_list(os.path.join(args.shifts_dir, fname))
            print(f"Querying STRING for {name} ({len(genes)} genes)")
            query_and_write(genes, os.path.join(out_dir, name), species=args.species)
    else:
        if not args.gene_list or not args.output_dir:
            parser.error("provide gene_list and output_dir, or use --shifts-dir/--output-dir for batch mode")
        os.makedirs(args.output_dir, exist_ok=True)
        name = os.path.splitext(os.path.basename(args.gene_list))[0]
        genes = read_gene_list(args.gene_list)
        print(f"Querying STRING for {name} ({len(genes)} genes)")
        query_and_write(genes, os.path.join(args.output_dir, name), species=args.species)
