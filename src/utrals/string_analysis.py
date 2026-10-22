"""STRING protein networks and pathway enrichment of RBP lists (bubble plots, pie-node networks)."""
# requests / networkx are imported only when needed (inside the STRING functions), so that the rest of this module
# (e.g. save_figure) also works in environments without them
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np

# STRING is pinned to release 12.5, and every request goes through utrals.string_api.string_api_post, which
# archives the raw response in cache/STRING_v12.5/ (git-tracked): later runs read the archive and never contact
# STRING, so the STRING figures stay identical (string-db.org itself follows the latest release, whose results change).
import os as _os
from .paths import CACHE_DIR as _CACHE_DIR
from .string_api import STRING_API_V12_5 as STRING_API, string_api_post

STRING_CACHE_DIR = _os.path.join(str(_CACHE_DIR), "STRING_v12.5")
SPECIES = 9606  # Human


def get_string_network(genes, species=SPECIES, min_channel_score=0.4):
    """Fetch STRING interactions filtered to experiments + databases only.

    Keeps edges where (escore > 0 OR dscore > 0) AND
    max(escore, dscore) >= min_channel_score (default 0.4 = medium confidence).
    """
    data = string_api_post(f"{STRING_API}/network",
                           {"identifiers": "\r".join(sorted(set(genes))), "species": species,
                            "required_score": 0, "caller_identity": "ALS_UTR_RBP"},
                           cache_dir=STRING_CACHE_DIR)
    if not data:
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df = df[(df["escore"] > 0) | (df["dscore"] > 0)].copy()
    df["edge_weight"] = df[["escore", "dscore"]].max(axis=1)
    df = df[df["edge_weight"] >= min_channel_score]
    return df.reset_index(drop=True)


# Manual overrides for symbols STRING cannot resolve at all (neither as
# official symbol nor as a known alias), e.g. {"OLD_SYMBOL": "NEW_SYMBOL"}.
MANUAL_ALIASES = {}


def _post_get_string_ids(names, species=SPECIES):
    """Raw `get_string_ids` call -> {query: entry} (best match per query)."""
    names = sorted(set(names))  # fixed order: the request (and its archived response) must not depend on set order
    data = string_api_post(f"{STRING_API}/get_string_ids",
                           {"identifiers": "\r".join(names), "species": species, "echo_query": 1,
                            "caller_identity": "ALS_UTR_RBP"},
                           cache_dir=STRING_CACHE_DIR)
    if data is None:
        return {}
    hits = {}
    for entry in data:
        query = names[entry["queryIndex"]]
        if query not in hits:
            hits[query] = entry
    return hits


def map_to_string(names, species=SPECIES, verbose=True):
    """Map gene symbols / aliases to STRING entries: {name: entry}.

    STRING resolves aliases itself (HUR -> ELAVL1, CGI99 -> RTRAF,
    HIST2H4A -> H4C6, ...), but in a batch query it returns only ONE query
    per STRING protein: if both "HUR" and "ELAVL1" are in the list, one of
    them silently gets no entry. Names missed by the batch are therefore
    re-queried one by one, so every alias maps to its protein. Names STRING
    does not know at all can be rescued via MANUAL_ALIASES.
    """
    names_unique = list(dict.fromkeys(n for n in names if isinstance(n, str) and n))
    if not names_unique:
        return {}
    queries = {n: MANUAL_ALIASES.get(n, n) for n in names_unique}
    batch_hits = _post_get_string_ids(list(dict.fromkeys(queries.values())), species)
    mapped = {n: batch_hits[q] for n, q in queries.items() if q in batch_hits}

    for n in [n for n in names_unique if n not in mapped]:
        hit = _post_get_string_ids([queries[n]], species).get(queries[n])
        if hit is not None:
            mapped[n] = hit

    if verbose:
        aliased = {n: e["preferredName"] for n, e in mapped.items() if n != e["preferredName"]}
        if aliased:
            print(f"  STRING alias -> preferred name: {aliased}")
        unresolved = [n for n in names_unique if n not in mapped]
        if unresolved:
            print(f"  STRING could not resolve (add to MANUAL_ALIASES?): {unresolved}")
    return mapped


def resolve_string_ids(names, species=SPECIES):
    """Map gene symbols to their (best-match) STRING identifiers.

    `background_string_identifiers` requires already-resolved STRING IDs
    (unlike `identifiers`, which STRING fuzzy-matches itself) -- passing raw
    gene symbols there makes the enrichment call fail silently with a
    generic "unknown_error".
    """
    mapped = map_to_string(names, species=species)
    return list(dict.fromkeys(e["stringId"] for e in mapped.values()))


def get_string_enrichment(genes, species=SPECIES, background_genes=None):
    """Query STRING enrichment, optionally against a custom statistical background.

    When `background_genes` is given, STRING computes enrichment relative to
    that gene set (via `background_string_identifiers`) instead of the whole
    genome/proteome. Both the query and the background must be passed as
    resolved STRING IDs in that case.
    """
    genes_unique = sorted(set(genes))
    payload = {
        "identifiers": "\r".join(genes_unique),
        "species": species,
        "caller_identity": "ALS_UTR_RBP",
    }
    if background_genes is not None:
        bg_ids = resolve_string_ids(set(background_genes) | set(genes_unique), species=species)
        query_ids = resolve_string_ids(genes_unique, species=species)
        if not query_ids or not bg_ids:
            return pd.DataFrame()
        payload["identifiers"] = "\r".join(sorted(query_ids))
        payload["background_string_identifiers"] = "\r".join(sorted(set(bg_ids) | set(query_ids)))
    data = string_api_post(f"{STRING_API}/enrichment", payload, cache_dir=STRING_CACHE_DIR)
    if not data:
        return pd.DataFrame()
    if isinstance(data, list) and data and "error" in data[0]:
        print(f"STRING enrichment error: {data[0].get('message', data[0])}")
        return pd.DataFrame()
    return pd.DataFrame(data)


def _layout(G, connected, isolated):
    """Kamada-Kawai for connected nodes; ring of isolated nodes around them."""
    import networkx as nx
    pos = {}

    if connected:
        sub = G.subgraph(connected)
        try:
            pos_c = nx.kamada_kawai_layout(sub)
        except Exception:
            pos_c = nx.spring_layout(sub, seed=42)
        pos.update(pos_c)

    if isolated:
        # place isolated nodes in a circle with radius just outside the connected cluster
        coords = np.array(list(pos.values())) if pos else np.array([[0.0, 0.0]])
        radius = max(1.5, np.max(np.linalg.norm(coords, axis=1)) + 1.0)
        angles = np.linspace(0, 2 * np.pi, len(isolated), endpoint=False)
        for node, angle in zip(isolated, angles):
            pos[node] = (radius * np.cos(angle), radius * np.sin(angle))

    if not pos:
        pos = nx.circular_layout(G)

    return pos


GROUP_COLORS = {"R1": "red", "R2": "orange", "R3": "blue"}


def build_node_colors(genes, group_lists, group_colors=GROUP_COLORS):
    """Map each gene to the list of colors for every group it belongs to.

    `group_lists` is e.g. {"R1": R1_top15, "R2": R2_top15, "R3": R3_top15}.
    A gene present in several groups' lists (e.g. both R1's and R2's top 15)
    gets one color per group, so `draw_pie_nodes` renders it as a multi-color
    pie node; a gene in a single group is drawn as a solid circle.
    """
    node_colors = {}
    for gene in genes:
        colors = [group_colors[label] for label, glist in group_lists.items() if gene in glist]
        node_colors[gene] = colors if colors else ["#8936ce"]
    return node_colors


def _pie_wedge_paths(n_slices):
    """Equal-sized wedge Paths of a unit circle, `n_slices` of them."""
    from matplotlib.path import Path
    paths = []
    for i in range(n_slices):
        start, end = i / n_slices, (i + 1) / n_slices
        angles = np.linspace(2 * np.pi * start, 2 * np.pi * end, 30)
        pts = np.vstack([[0, 0], np.column_stack([np.cos(angles), np.sin(angles)]), [0, 0]])
        paths.append(Path(pts))
    return paths


def draw_pie_nodes(ax, pos, node_colors, node_size=500, alpha=0.9):
    """Draw each node as a circle, split into one equal wedge per group color.

    `node_colors`: {node: [color1, color2, ...]}.
    """
    for node, colors in node_colors.items():
        x, y = pos[node]
        for wedge, color in zip(_pie_wedge_paths(len(colors)), colors):
            ax.scatter([x], [y], marker=wedge, s=node_size, facecolor=color,
                       edgecolor="none", alpha=alpha, zorder=3)


def get_preferred_names(genes, species=SPECIES):
    """Map each input gene symbol to STRING's canonical preferredName.

    STRING resolves aliases to its own official symbol in network/enrichment
    responses (e.g. the query "HUR" comes back labeled "ELAVL1"), which can
    differ from the alias used in our RBP lists. Resolving up front keeps
    node identity -- and therefore node coloring -- consistent regardless of
    which alias was used as input.
    """
    mapped = map_to_string(genes, species=species, verbose=False)
    return {name: entry["preferredName"] for name, entry in mapped.items()}


def plot_string_network(genes, region_name, min_channel_score=0.4, ax=None, node_colors=None):
    import networkx as nx
    # sorted, not set order (which changes between Python processes): the node order sets the layout
    genes_unique = sorted(set(genes))
    df_net = get_string_network(genes_unique, min_channel_score=min_channel_score)
    name_map = get_preferred_names(genes_unique)

    G = nx.Graph()

    # Add edges using STRING preferred names
    if not df_net.empty:
        for _, row in df_net.iterrows():
            G.add_edge(row["preferredName_A"], row["preferredName_B"],
                       escore=row["escore"], dscore=row["dscore"])

    # Add genes that did not appear in any edge as isolated nodes, resolved
    # to the same preferredName space as connected nodes
    for gene in genes_unique:
        node = name_map.get(gene, gene)
        if node not in G.nodes():
            G.add_node(node)

    isolated = [n for n in G.nodes() if G.degree(n) == 0]
    connected = [n for n in G.nodes() if G.degree(n) > 0]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))

    pos = _layout(G, connected, isolated)

    node_colors = node_colors or {}
    node_colors_resolved = {name_map.get(g, g): c for g, c in node_colors.items()}
    draw_pie_nodes(
        ax, pos,
        {n: node_colors_resolved.get(n, ["#8936ce"]) for n in G.nodes()},
        node_size=500, alpha=0.9,
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold",
                            font_color="black")
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="black", width=1.5, alpha=0.7)

    ax.set_title(
        f"{region_name}\n"
        f"connected: {len(connected)}, isolated: {len(isolated)}  "
        f"[experiments + databases, score ≥ {min_channel_score}]",
        fontsize=11, weight="bold",
    )
    ax.axis("off")
    if isolated:
        print(f"  {region_name} — no qualifying interactions: {isolated}")

    return G, df_net






def dedupe_redundant_terms(df, jaccard_threshold=0.5):
    """Curate redundant pathway terms (e.g. near-duplicate GO/KEGG/Reactome
    entries covering the same genes), keeping only the most significant one.

    Terms are considered redundant when the Jaccard similarity of their
    STRING-resolved gene sets (`preferredNames`) is >= jaccard_threshold.
    Greedily walks terms from most to least significant (lowest FDR first)
    and drops any term that overlaps too much with an already-kept term.
    """
    if df.empty:
        return df
    kept_rows = []
    kept_sets = []
    for _, row in df.sort_values("fdr").iterrows():
        genes = set(row["preferredNames"])
        if any(len(genes & s) / len(genes | s) >= jaccard_threshold for s in kept_sets):
            continue
        kept_rows.append(row)
        kept_sets.append(genes)
    return pd.DataFrame(kept_rows)


def get_enrichment_df(genes, categories=("KEGG", "Reactome", "Process"),
                      species=SPECIES, background_genes=None, background_size=None,
                      jaccard_threshold=0.5):
    """Return filtered, redundancy-curated enrichment table with STRING strength
    and signal computed.

    strength = log10(observed / expected)  where expected = n * K / N
    signal   = harmonic mean of strength and -log10(FDR)
               (balances large-term bias of FDR with small-term bias of O/E ratio)
    n = network size inferred from the union of preferredNames across all returned terms.
    N = background_size, i.e. the statistical background gene set. When
        `background_genes` is provided, enrichment is computed relative to
        that custom set (e.g. all other RBPs in the analysis) rather than
        the whole genome, and N defaults to its size.
    Near-duplicate pathway terms are curated away via `dedupe_redundant_terms`.
    """
    df = get_string_enrichment(genes, species=species, background_genes=background_genes)
    if df.empty:
        return pd.DataFrame()
    df = df[df["category"].isin(categories)].copy()
    if df.empty:
        return pd.DataFrame()

    if background_size is None:
        background_size = len(set(background_genes)) if background_genes is not None else 19566

    # Network size: union of all STRING-resolved gene names in the response
    all_prefs = set()
    for prefs in df["preferredNames"]:
        all_prefs.update(prefs)
    n = max(len(all_prefs), 1)

    k = df["number_of_genes"].astype(float)
    K = df["number_of_genes_in_background"].astype(float)

    expected = n * K / background_size
    df["strength"] = np.log10((k / expected).clip(lower=1e-9))
    df["log_fdr"]  = -np.log10(df["fdr"].clip(lower=1e-300))

    # Harmonic mean of strength and -log10(FDR)
    s  = df["strength"].clip(lower=1e-9)
    lf = df["log_fdr"].clip(lower=1e-9)
    df["signal"] = 2 * s * lf / (s + lf)

    df["gene_ratio"] = k / K

    df = dedupe_redundant_terms(df, jaccard_threshold=jaccard_threshold)
    return df


def plot_enrichment_combined(region_genes, top_n=10,
                             categories=("KEGG", "Reactome", "Process"),
                             background_genes=None, jaccard_threshold=0.5, ax=None):
    """Single bubble plot: rows = pathways, columns = regions.
    Bubble size = -log10(FDR), colour = STRING signal.
    """
    region_names = list(region_genes.keys())

    enrich_per_region = {}
    for rname, genes in region_genes.items():
        df = get_enrichment_df(list(set(genes)), categories=categories,
                                background_genes=background_genes,
                                jaccard_threshold=jaccard_threshold)
        enrich_per_region[rname] = df

    # Union of top_n pathways per region
    top_terms = set()
    for df in enrich_per_region.values():
        if not df.empty:
            top_terms |= set(df.sort_values("fdr").head(top_n)["description"])

    if not top_terms:
        print("No enriched pathways found in any region.")
        return pd.DataFrame()

    # Build long-form table
    rows = []
    for rname, df in enrich_per_region.items():
        if df.empty:
            continue
        for _, row in df[df["description"].isin(top_terms)].iterrows():
            rows.append({
                "region":      rname,
                "description": row["description"],
                "log_fdr":     row["log_fdr"],
                "signal":      row["signal"],
            })
    df_plot = pd.DataFrame(rows)

    # Order pathways: best mean signal on top
    order = (df_plot.groupby("description")["signal"]
             .mean().sort_values(ascending=True).index.tolist())

    if ax is None:
        fig, ax = plt.subplots(figsize=(3 + len(region_names),
                                        1 + 0.4 * len(order)))

    norm_color = plt.Normalize(df_plot["signal"].min(), df_plot["signal"].max())
    norm_size  = plt.Normalize(0, df_plot["log_fdr"].max())
    cmap = plt.cm.YlOrRd

    x_map = {rname: i for i, rname in enumerate(region_names)}

    for _, row in df_plot.iterrows():
        x = x_map[row["region"]]
        y = order.index(row["description"])
        size  = (norm_size(row["log_fdr"]) * 400) + 20
        color = cmap(norm_color(row["signal"]))
        ax.scatter(x, y, s=size, color=color, edgecolors="grey",
                   linewidths=0.4, zorder=3)

    ax.set_xticks(range(len(region_names)))
    ax.set_xticklabels(region_names, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=9)
    ax.set_xlim(-0.5, len(region_names) - 0.5)
    ax.set_ylim(-0.5, len(order) - 0.5)
    ax.grid(linestyle="--", alpha=0.3)
    sns.despine(ax=ax, left=True, bottom=True)
    ax.tick_params(length=0)
    ax.set_title("Pathway enrichment across regions", fontsize=12, weight="bold")

    # Colourbar for signal
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm_color)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, shrink=0.4)
    cbar.set_label("Signal", fontsize=9)

    # Size legend for -log10(FDR)
    fdr_vals = [1, 2, 5]
    size_handles = [
        plt.scatter([], [], s=(norm_size(v) * 400) + 20,
                    c="grey", edgecolors="grey", alpha=0.7, label=str(v))
        for v in fdr_vals if v <= df_plot["log_fdr"].max()
    ]
    if size_handles:
        ax.legend(handles=size_handles, title="-log₁₀(FDR)",
                  loc="lower left", fontsize=8, title_fontsize=8, frameon=True)

    return df_plot
