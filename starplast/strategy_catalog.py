"""The strategies: thirty-two ways to use the combined data for inference, each testing itself.

Registered into `strategies.REGISTRY` on import. Each entry is a runner, a tester, and the prose the
Strategies tab shows -- tooltip, explanation, walkthrough and a description of its self-test. They
are grouped into eight families by MECHANISM rather than by question, because the mechanism is what
decides what a strategy can be trusted for:

    1  Search the map space          maps, walks and held-out labels             01-06
    2  Borrow from neighbours        nearest genes in the measurements           07-10
    3  Walk the networks             the measured edge layers                    11-18
    4  Learn from examples           supervised models, scored out of sample     19-23
    5  Start from a gene list        a list in, a ranking and a profile out      24-25
    6  Contrast and combine layers   where two kinds of evidence part or meet    26-28
    7  Cross species and strata      orthologs, and the genes orthology misses   29-30
    8  Combine strategies            agreement, and the understudied genes       31-32

Every tester follows one of the five patterns in `strategies`, so every verdict means the same
thing: the strategy beat the 95th percentile of the same procedure run without the information it
claims to use, by at least the stated margin.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from . import strategies as S
from .strategies import (NOISE, MIN_CLASS, Param, Strategy, StrategyResult, register,
                         default_category, default_numeric, example_set, expand, parse_grid)

MAP, NEIGHBOURS, NETWORKS, LEARN, LISTS, CONTRAST, SPECIES, COMBINE = (
    "Search the map space", "Borrow from neighbours", "Walk the networks", "Learn from examples",
    "Start from a gene list", "Contrast and combine layers", "Cross species and strata",
    "Combine strategies")

# --------------------------------------------------------------------------- shared parameters
TARGET = Param("target", "category", "Held-out label",
               "The label the strategy is scored against and never allowed to see: the column "
               "itself, anything that restates it and the experiment that produced it are removed "
               "first, by the same closure the Search tab uses.", default_category)
SAMPLE = Param("sample", "int", "Genes per map",
               "How many genes each map embeds. Zero embeds all of them; a sample of a few thousand "
               "makes a walk of many maps finish in minutes rather than an hour, at the cost of "
               "saying nothing about the genes left out.", 3000, lo=0, hi=20000, step=500)
NEIGHBORS = Param("n_neighbors", "grid", "UMAP neighbours to try",
                  "Comma-separated UMAP n_neighbors values to walk. Small values keep local detail "
                  "such as complexes; large values keep the global arrangement such as organelles. "
                  "Every value multiplies the number of maps built.", "15, 50")
MIN_DIST = Param("min_dist", "grid", "UMAP min_dist to try",
                 "Comma-separated UMAP min_dist values to walk. Zero packs similar genes tightly, "
                 "which clusters well; larger values spread them out, which reads better and "
                 "clusters worse. Every value multiplies the number of maps.", "0.0, 0.25")
MIN_CLUSTER = Param("min_cluster_size", "grid", "Cluster sizes to try",
                    "Comma-separated HDBSCAN minimum cluster sizes. Small values find small "
                    "complexes and also split organelles; large values find organelles and swallow "
                    "complexes. Cheap to widen: clustering is fast next to embedding.", "20, 50")
GENES = Param("genes", "genes", "Gene list",
              "Gene accessions, one per line or separated by commas or spaces, pasted from anywhere. "
              "Matched exactly against this organism's table, ignoring case; anything not found is "
              "reported by name rather than silently dropped.", "")
EXCLUDE = Param("exclude", "category", "Column the list came from",
                "If the gene list was made from a column of this table, name it here and that column, "
                "everything restating it and its experiment are withheld; otherwise the strategy "
                "would find the list again by reading the column that defined it.", None,
                optional=True)
SELECTIONS = Param("selection", "grid", "Cluster selection to try",
                   "HDBSCAN's cluster selection, comma-separated: 'eom' keeps the most persistent "
                   "clusters -- on the real proteome two to six very large ones -- and 'leaf' every "
                   "leaf of the cluster tree, dozens of small purer clusters with more genes left as "
                   "noise. Walking both costs little.", "eom, leaf")
SELECTION = Param("selection", "choice", "Cluster selection",
                  "HDBSCAN's cluster selection for the one map this strategy reads. 'leaf' gives "
                  "many small clusters, which is what naming a cluster needs; 'eom' gives a few "
                  "very large ones, which on this proteome no single label dominates.", "leaf",
                  choices=("leaf", "eom"))
K = Param("k", "int", "Neighbours to consult",
          "How many labelled genes vote on each call. Few neighbours follow fine local structure "
          "and are noisy; many are stable and blur small classes into large ones. Fifteen is the "
          "project's default everywhere a neighbourhood is scored.", 15, lo=3, hi=200)


def _layer_default(preferred):
    """A default that picks the first preferred layer this context actually has."""
    def pick(ctx):
        have = ctx.layers()
        for layer in preferred:
            if layer in have:
                return layer
        return have[0] if have else None
    return pick


def _need(ctx, column):
    """The column, or a ValueError a person can act on."""
    if not column or column not in ctx.nodes:
        raise ValueError(f"choose a column this table has; {column!r} is not one")
    return column


def _genes_table(ctx, positions, **columns) -> pd.DataFrame:
    """Gene id and product for some positions, plus whatever per-gene columns are passed."""
    positions = np.asarray(positions, dtype=int)
    out = pd.DataFrame({"gene_id": ctx.gene_ids[positions], "product": ctx.product(positions)})
    for k, v in columns.items():
        out[k] = np.asarray(v)
    return out


def _calls_table(ctx, pred: pd.Series, support: pd.Series, known: pd.Series,
                 limit: int = 2000, **extra) -> pd.DataFrame:
    """The genes a strategy calls that had no label of their own, strongest support first."""
    pred = pd.Series(pred).reset_index(drop=True)
    support = pd.Series(support).reset_index(drop=True)
    unknown = pd.Series(known).reset_index(drop=True).isna().to_numpy()
    rows = np.flatnonzero(pred.notna().to_numpy() & unknown)
    table = _genes_table(ctx, rows, prediction=pred.iloc[rows].to_numpy(),
                         support=support.iloc[rows].to_numpy(dtype=float),
                         **{k: np.asarray(v)[rows] for k, v in extra.items()})
    return table.sort_values("support", ascending=False, kind="stable").head(limit).reset_index(
        drop=True)


def _walk(ctx, target, blocksets: dict, nns, mds, mcss, rows, selections=("eom",)) -> list:
    """Every (feature set x n_neighbors x min_dist x min_cluster_size x selection) clustering.

    The maps are cached on the context, so a second strategy walking the same grid pays nothing.
    """
    out = []
    total = len(blocksets) * len(nns) * len(mds)
    i = 0
    for name, blocks in blocksets.items():
        for nn in nns:
            for md in mds:
                i += 1
                ctx.say(f"map {i} of {total}: {name}, n_neighbors {nn}, min_dist {md}")
                coords, pos = ctx.embed(blocks, nn, md, rows, target=target)
                for mcs in mcss:
                    for sel in selections:
                        out.append({"features": name, "n_neighbors": int(nn),
                                    "min_dist": float(md), "min_cluster_size": int(mcs),
                                    "selection": str(sel), "coords": coords, "positions": pos,
                                    "labels": ctx.cluster(coords, mcs, sel)})
    return out


def _blocksets(ctx, target, mode: str) -> dict:
    """The feature sets a walk tries: everything permitted, and optionally each family alone."""
    blocks = ctx.blocks(target)
    if not blocks:
        raise ValueError("nothing is left to build a map from once the held-out label's closure "
                         "is removed")
    sets = {"all permitted": tuple(blocks)}
    if mode == "families":
        for fam, members in ctx.families(target).items():
            if sum(len(blocks[b]) for b in members) >= 3 and len(members) < len(blocks):
                sets[fam] = tuple(members)
    return sets


def _config_label(c) -> str:
    return (f"{c['features']}, n_neighbors {c['n_neighbors']}, min_dist {c['min_dist']}, "
            f"min cluster {c['min_cluster_size']} ({c.get('selection', 'eom')})")


def _label_clusters(labels, truth_sub: pd.Series, min_label: int = 10) -> pd.DataFrame:
    """Per label, the ONE cluster that best isolates it (F1), as the Search tab scores a map."""
    from . import search
    _s, per = search.score_recovery(np.asarray(labels), pd.Series(truth_sub).reset_index(drop=True),
                                    min_label=min_label)
    return per


def _hidden_f1(labels, positions, per: pd.DataFrame, truth: pd.Series, hidden) -> float:
    """How well each label's chosen cluster recovers that label's HIDDEN genes, size-weighted.

    The cluster for each label is chosen on visible genes; this scores it on genes the choice never
    saw. A hidden label with no chosen cluster scores zero, so a mapping that covers fewer labels
    is not rewarded for it.
    """
    at = {g: i for i, g in enumerate(np.asarray(positions))}
    h = np.array([at[g] for g in hidden if g in at], dtype=int)
    if not len(h):
        return float("nan")
    hl = np.asarray(labels)[h]
    ht = pd.Series(truth).iloc[np.asarray(positions)[h]].to_numpy(dtype=object)
    chosen = dict(zip(per["label"], per["cluster"])) if len(per) else {}
    total = weight = 0.0
    for lab in {x for x in ht if isinstance(x, str)}:
        is_l = ht == lab
        n_l = int(is_l.sum())
        weight += n_l
        k = chosen.get(lab)
        if k is None:
            continue
        in_k = hl == k
        tp = int((is_l & in_k).sum())
        if tp:
            prec, rec = tp / int(in_k.sum()), tp / n_l
            total += n_l * 2 * prec * rec / (prec + rec)
    return total / weight if weight else float("nan")


# =========================================================================== 1 · map space
def _holdout_search_run(ctx, p):
    from . import search
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    rows = ctx.sample_rows(p["sample"], target)
    configs = _walk(ctx, target, _blocksets(ctx, target, p["features"]),
                    parse_grid(p["n_neighbors"], int), parse_grid(p["min_dist"]),
                    parse_grid(p["min_cluster_size"], int), rows, parse_grid(p["selection"], str))
    records, per = [], {}
    for i, c in enumerate(configs):
        sub = truth.iloc[c["positions"]].reset_index(drop=True)
        summary, table = search.score_recovery(c["labels"], sub, min_label=MIN_CLASS)
        per[i] = table
        records.append({"configuration": i, "features": c["features"],
                        "n_neighbors": c["n_neighbors"], "min_dist": c["min_dist"],
                        "min_cluster_size": c["min_cluster_size"], "selection": c["selection"],
                        "clusters": int(len(set(c["labels"].tolist()) - {NOISE})),
                        "noise": float((c["labels"] == NOISE).mean()),
                        "mean_f1": summary.get("mean_f1", 0.0), "best_f1": summary.get("best_f1", 0.0),
                        "best_label": summary.get("best_label", ""),
                        "labels_recovered": summary.get("n_labels_recovered", 0)})
    ranking = pd.DataFrame(records).sort_values("mean_f1", ascending=False, kind="stable")
    best_i = int(ranking.iloc[0]["configuration"])
    best = configs[best_i]
    sub = truth.iloc[best["positions"]].reset_index(drop=True)
    rng = ctx.rng(5)
    null = []
    for _ in range(20):
        shuffled = pd.Series(rng.permutation(sub.to_numpy(dtype=object)))
        s, _t = search.score_recovery(best["labels"], shuffled, min_label=MIN_CLASS)
        null.append(s.get("mean_f1", 0.0))
    # Each label's best cluster, if it isolates the label well enough, names its unlabelled
    # members. A cluster that is best for two labels names the one it isolates better.
    chosen = per[best_i]
    chosen = chosen[chosen["f1"] >= float(p["min_f1"])] if len(chosen) else chosen
    chosen = chosen.sort_values("f1", ascending=False).drop_duplicates("cluster") \
        if len(chosen) else chosen
    cl_full = expand(best["labels"], best["positions"], ctx.n)
    pred = pd.Series([np.nan] * ctx.n, dtype=object)
    support = pd.Series(np.nan, index=range(ctx.n))
    for r in chosen.itertuples():
        members = cl_full == r.cluster
        pred[members] = r.label
        support[members] = r.f1
    calls = _calls_table(ctx, pred, support, truth, cluster=cl_full)
    clusters = chosen.rename(columns={"f1": "f1_on_known_genes"}) if len(chosen) else chosen
    top = ranking.iloc[0]
    summary = (f"Best of {len(configs)} configurations: {_config_label(best)}. Mean F1 "
               f"{top['mean_f1']:.3f} over the {target} labels, against {np.mean(null):.3f} ± "
               f"{np.std(null):.3f} for the same clustering with the labels shuffled; "
               f"{int(top['labels_recovered'])} label(s) recovered at F1 >= 0.5, best "
               f"{top['best_label']} at {top['best_f1']:.2f}. {len(calls):,} unlabelled genes sit in "
               f"the cluster that best isolates a label at F1 >= {float(p['min_f1']):.2f}, and are "
               f"called with it; 'support' is that F1.")
    return StrategyResult("holdout_search", summary,
                          {"configurations": ranking.drop(columns="configuration"),
                           "categories": per[best_i], "calls": calls, "label clusters": clusters},
                          labels=expand(best["labels"], best["positions"], ctx.n),
                          coords=best["coords"], positions=best["positions"],
                          genes=list(calls["gene_id"][:200]) if len(calls) else [],
                          numbers={"mean_f1": float(top["mean_f1"]),
                                   "null_mean_f1": float(np.mean(null))})


def _holdout_search_test(ctx, p):
    from . import search
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    visible, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    rows = ctx.sample_rows(min(p["sample"] or ctx.n, 1500), target)
    configs = _walk(ctx, target, {"all permitted": tuple(ctx.blocks(target))},
                    parse_grid(p["n_neighbors"], int)[:2], parse_grid(p["min_dist"])[:2],
                    parse_grid(p["min_cluster_size"], int)[:2], rows,
                    parse_grid(p["selection"], str)[:2])
    placed = np.unique(np.concatenate([c["positions"] for c in configs]))
    hidden = hidden[np.isin(hidden, placed)]

    def infer(vis):
        """The configuration and each label's cluster, both chosen on visible labels only."""
        best, best_f, best_per = configs[0], -1.0, pd.DataFrame()
        for c in configs:
            s, per = search.score_recovery(c["labels"], vis.iloc[c["positions"]].reset_index(
                drop=True), min_label=10)
            if s and s["mean_f1"] > best_f:
                best, best_f, best_per = c, s["mean_f1"], per
        return _hidden_f1(best["labels"], best["positions"], best_per, truth, hidden)

    observed = infer(visible)
    rng = ctx.rng(7)
    nulls = []
    for i in range(20):
        ctx.check()
        nulls.append(infer(S.shuffled(visible, rng)))
    return S.judge("holdout_search",
                   "F1 of hidden genes in the cluster chosen for their label on known genes",
                   observed, nulls, min_effect=0.05, n_hidden=len(hidden),
                   hidden=f"25% of the {target} labels, whole orthogroups at a time",
                   null_kind="20 searches whose map and clusters were chosen on shuffled labels",
                   t0=t0, numbers={"maps": len(configs), "genes_per_map": len(rows)})


register(Strategy(
    key="holdout_search", number=1, family=MAP,
    title="Hold out a category and search for a map that finds it",
    question="Is there a combination of measurements and map settings under which a label nobody "
             "showed the map falls out as clusters -- and which unlabelled genes land in them?",
    tooltip="Hides one label entirely, walks feature sets and UMAP/HDBSCAN settings, and keeps the "
            "map whose clusters best recover that label; unlabelled genes in its pure clusters "
            "become candidates, with the shuffled-label score beside them.",
    explanation=(
        "This is the question the map was built to ask. Choose a label -- the hyperLOPIT "
        "compartment, a stage, an essentiality class -- and the strategy removes it from the "
        "inputs together with everything the leakage closure says restates it: the column, "
        "columns measurably associated with it, the experiment that produced it, and any edge "
        "layer built from it. What is left is embedded under every combination of UMAP "
        "n_neighbors and min_dist on the grid, each map is clustered at each HDBSCAN size, and "
        "each clustering is scored by how well single clusters isolate each category (the best "
        "cluster's F1 per category, weighted by category size).\n\n"
        "Why it can work: if a compartment is a biological fact, its members share an expression "
        "programme, a fitness pattern, a set of modifications -- so they should sit together in a "
        "map built only from those. Why it can fail: a walk picks the luckiest of many maps, so the "
        "best score is optimistic; every run therefore reports the same clustering's score with "
        "the labels shuffled, and the self-test chooses the map on some labels and scores it on "
        "labels it never saw. Each label is then given the ONE cluster that best isolates it, and "
        "that cluster's unlabelled members are called with it, each call carrying the cluster's "
        "F1 -- a pure-looking cluster is a lead, not a result.\n\n"
        "Choose 'families' to also walk each kind of measurement alone -- transcription only, "
        "fitness only -- which tells you WHICH evidence carries the category, not just whether "
        "the combination does."),
    walkthrough=(
        "Pick the held-out label. Compartment is the classic first question.",
        "Leave the grids at their defaults for a first pass (4 maps x 2 cluster sizes x 2 "
        "selections); widen n_neighbors to '10, 25, 50, 100' once you know the label is "
        "recoverable at all.",
        "Set 'Genes per map' to 3000 for speed, 0 for the whole proteome.",
        "Press Test first: it hides a quarter of the label, lets the search choose the map and "
        "each label's cluster on the rest, and reports how well those clusters hold the hidden "
        "genes, against searches made on shuffled labels. A FAIL means the search has nothing to "
        "find on this label.",
        "Press Run. The best map is drawn, colored by its clusters; 'configurations' ranks every "
        "map, 'categories' says which labels were recovered and at what F1.",
        "Read 'label clusters' -- the cluster chosen for each label and how well it isolates it -- "
        "then 'calls', the unlabelled genes in those clusters. Click a row to find the gene; "
        "treat a call's F1 as the most it can be trusted."),
    test_description=(
        "25% of the label is hidden (whole orthogroups together, so no gene is recovered through "
        "a visible paralog). A small walk is built on 1,500 genes; the configuration, and the one "
        "cluster that best isolates each label, are both chosen using visible labels only. "
        "Metric: for each label, the F1 of its hidden genes against its chosen cluster, weighted "
        "by size -- does the structure found on known genes hold the unknown ones? Null: the whole "
        "search repeated 20 times with the visible labels shuffled. Pass: above the null's 95th "
        "percentile by at least 0.05."),
    params=(TARGET, NEIGHBORS, MIN_DIST, MIN_CLUSTER, SELECTIONS, SAMPLE,
            Param("features", "choice", "Feature sets",
                  "'all' walks one map per setting from every permitted measurement; 'families' "
                  "also walks each kind of measurement alone, so you learn which evidence carries "
                  "the label -- at several times the cost.", "all", choices=("all", "families")),
            Param("min_f1", "float", "Call a label's cluster at F1",
                  "Each label's best cluster names its unlabelled members only if it isolates that "
                  "label at least this well on the known genes. Lower names more clusters and is "
                  "wrong more often; every call carries its cluster's F1 as its support.", 0.2,
                  lo=0.0, hi=1.0, step=0.05)),
    runner=_holdout_search_run, tester=_holdout_search_test, cost="minutes",
    needs=("a categorical column",)))


def _geneset_hunt_run(ctx, p):
    pos, missing = ctx.resolve_genes(p["genes"])
    if len(pos) < 5:
        return StrategyResult("geneset_hunt", "Paste at least five gene ids from this organism. "
                              + (f"Not found: {', '.join(missing[:20])}" if missing else ""), {})
    exclude = p.get("exclude")
    rows = ctx.sample_rows(p["sample"], must=pos)
    configs = _walk(ctx, exclude, _blocksets(ctx, exclude, p["features"]),
                    parse_grid(p["n_neighbors"], int), parse_grid(p["min_dist"]),
                    parse_grid(p["min_cluster_size"], int), rows, parse_grid(p["selection"], str))
    records = []
    for i, c in enumerate(configs):
        mask = np.isin(c["positions"], pos)
        k, prec, rec, f1 = S.set_f1(c["labels"], mask)
        records.append({"configuration": i, "features": c["features"],
                        "n_neighbors": c["n_neighbors"], "min_dist": c["min_dist"],
                        "min_cluster_size": c["min_cluster_size"], "selection": c["selection"],
                        "cluster": k,
                        "cluster_size": int((c["labels"] == k).sum()) if k != NOISE else 0,
                        "list_genes_in_it": int((mask & (c["labels"] == k)).sum()),
                        "precision": prec, "recall": rec, "f1": f1})
    ranking = pd.DataFrame(records).sort_values("f1", ascending=False, kind="stable")
    best_i = int(ranking.iloc[0]["configuration"])
    best, top = configs[best_i], ranking.iloc[0]
    # The search's own null: random sets of the same size, each given its best cluster over EVERY
    # configuration, because the observed F1 is itself a maximum over configurations.
    rng = ctx.rng(9)
    placed = best["positions"]
    size = int(np.isin(placed, pos).sum())
    null = []
    for _ in range(int(p["random_sets"])):
        ctx.check()
        r = rng.choice(placed, size=size, replace=False)
        null.append(max(S.set_f1(c["labels"], np.isin(c["positions"], r))[3] for c in configs))
    p_value = (1 + sum(v >= top["f1"] for v in null)) / (1 + len(null))
    in_cluster = best["positions"][best["labels"] == int(top["cluster"])]
    centre = best["coords"][best["labels"] == int(top["cluster"])].mean(axis=0) \
        if len(in_cluster) else np.zeros(3)
    dist = np.linalg.norm(best["coords"][best["labels"] == int(top["cluster"])] - centre, axis=1)
    members = _genes_table(ctx, in_cluster, in_list=np.isin(in_cluster, pos),
                           distance_to_centre=dist).sort_values("distance_to_centre")
    candidates = members[~members["in_list"]].reset_index(drop=True)
    missed_pos = np.setdiff1d(np.intersect1d(pos, placed), in_cluster)
    missed = _genes_table(ctx, missed_pos)
    summary = (f"{len(pos)} genes found ({len(missing)} not found). Best of {len(configs)} "
               f"configurations: {_config_label(best)} -- cluster {int(top['cluster'])} holds "
               f"{int(top['list_genes_in_it'])} of your {size} placed genes among "
               f"{int(top['cluster_size'])} (precision {top['precision']:.2f}, recall "
               f"{top['recall']:.2f}, F1 {top['f1']:.2f}). Random sets of the same size reach F1 "
               f"{np.mean(null):.2f} ± {np.std(null):.2f} at their best over the same walk "
               f"(search-corrected p = {p_value:.3g}). {len(candidates)} genes in that cluster are "
               f"not on your list.")
    return StrategyResult("geneset_hunt", summary,
                          {"configurations": ranking.drop(columns="configuration"),
                           "candidates": candidates, "cluster members": members.reset_index(
                               drop=True), "list genes not in the cluster": missed,
                           "not found": pd.DataFrame({"identifier": missing})},
                          labels=expand(best["labels"], best["positions"], ctx.n),
                          coords=best["coords"], positions=best["positions"],
                          genes=list(ctx.gene_ids[pos]),
                          numbers={"f1": float(top["f1"]), "null_f1": float(np.mean(null)),
                                   "p_value": float(p_value)})


def _list_or_example(ctx, p) -> tuple:
    """(positions, what they are, column to exclude): the user's list if long enough, else a known set."""
    pos, _missing = ctx.resolve_genes(p.get("genes", ""))
    if len(pos) >= 20:
        return pos, f"your list of {len(pos)} genes", p.get("exclude")
    # The default label first; if none of its categories is a usable size -- Plasmodium's
    # phenotype classes are all larger than a list anyone would hand over -- any other label.
    first = p.get("exclude") or default_category(ctx)
    for target in [first] + [c for c in ctx.categorical_columns() if c != first]:
        cat, members = example_set(ctx, target)
        if cat is not None:
            return members, f"the {cat!r} genes of {target}", target
    return np.array([], dtype=int), "no known set", first


def _geneset_hunt_test(ctx, p):
    t0 = time.monotonic()
    members, what, exclude = _list_or_example(ctx, p)
    if len(members) < 10:
        return S.judge("geneset_hunt", "share of hidden members in the best cluster", float("nan"),
                       [], min_effect=0.1, n_hidden=0, hidden=what, null_kind="random sets",
                       t0=t0, note="no gene set of a usable size to test on")
    rows = ctx.sample_rows(1500, must=members)
    configs = _walk(ctx, exclude, {"all permitted": tuple(ctx.blocks(exclude))},
                    parse_grid(p["n_neighbors"], int)[:2], parse_grid(p["min_dist"])[:1],
                    parse_grid(p["min_cluster_size"], int)[:2], rows,
                    parse_grid(p["selection"], str)[:2])
    rng = ctx.rng(13)

    def once(genes):
        genes = rng.permutation(genes)
        k = int(round(0.3 * len(genes)))
        hidden, query = genes[:k], genes[k:]
        best, best_f, best_k = None, -1.0, NOISE
        for c in configs:
            cl, _p, _r, f1 = S.set_f1(c["labels"], np.isin(c["positions"], query))
            if f1 > best_f:
                best, best_f, best_k = c, f1, cl
        inside = best["positions"][best["labels"] == best_k] if best_k != NOISE else np.array([])
        # The cluster's OTHER genes are its candidates: how many of them are the hidden members,
        # and how many hidden members they include. F1 of the two, so a giant cluster that holds a
        # third of any set scores its precision, not its size.
        candidates = np.setdiff1d(inside, query)
        hits = int(np.isin(candidates, hidden).sum())
        if not hits:
            return 0.0, len(hidden)
        prec, rec = hits / len(candidates), hits / len(hidden)
        return float(2 * prec * rec / (prec + rec)), len(hidden)

    observed, n_hidden = once(members)
    placed = configs[0]["positions"]
    nulls = []
    for i in range(20):
        ctx.check()
        nulls.append(once(rng.choice(placed, size=len(members), replace=False))[0])
    return S.judge("geneset_hunt", "F1 of the hidden members against the best cluster's other genes",
                   observed, nulls, min_effect=0.05, n_hidden=n_hidden,
                   hidden=f"30% of {what}", null_kind="20 random sets of the same size, same walk",
                   t0=t0, numbers={"maps": len(configs)})


register(Strategy(
    key="geneset_hunt", number=2, family=MAP,
    title="Find the map where your gene list is one cluster",
    question="Under some combination of measurements and settings, do the genes on my list fall "
             "into a single cluster -- and what else is in it?",
    tooltip="Walks maps and clusterings looking for one cluster with high precision AND recall "
            "for your list; reports that cluster's other members as candidates, and how often a "
            "random list of the same size does as well over the same walk.",
    explanation=(
        "Hand over a list -- the hits of a screen, the members of a complex, genes a paper "
        "implicates -- and the strategy walks the same space the Search tab walks, but scores each "
        "clustering by the single cluster that best captures YOUR list: precision (what share of "
        "the cluster is on the list) and recall (what share of the list is in the cluster), "
        "combined as F1. A list that is a real biological unit should, under some view of the "
        "data, gather into one cluster; the genes that gather with it are the candidates.\n\n"
        "The trap is the walk itself: over many maps and cluster sizes SOMETHING will collect a "
        "fair share of any list. So every run also scores random lists of the same size over the "
        "same maps and reports the best each achieves -- a search-corrected p-value. A list that "
        "does no better than random lists is not a unit this data can see.\n\n"
        "If the list was made from a column of this table (all genes of one compartment, say), "
        "name that column under 'Column the list came from': it is then withheld, as a held-out "
        "label would be, and the answer is not simply the column read back."),
    walkthrough=(
        "Paste your accessions into the gene list box (one per line, or comma separated).",
        "If the list came from a column of this table, choose it as the column to exclude.",
        "Press Test: with a list of 20 or more genes it hides 30% of YOUR list and asks whether "
        "the best cluster for the rest contains them; with a shorter list it tests on a known "
        "category of similar size. A FAIL on your own list means it is not coherent in this data.",
        "Press Run. The winning map is drawn; the table 'configurations' shows every map's best "
        "cluster with precision, recall and F1.",
        "Read the summary's search-corrected p-value before the candidates: it is the chance a "
        "random list of the same size does as well over the same walk.",
        "'candidates' lists the cluster's other members, nearest the cluster centre first."),
    test_description=(
        "30% of the set is hidden; the walk (1,500 genes, 2 n_neighbors x 2 cluster sizes x both "
        "selections) picks the cluster with the best F1 for the other 70%. Metric: F1 of the "
        "hidden members against that cluster's other genes -- precision is the share of the "
        "cluster's candidates that are hidden members, recall the share of hidden members among "
        "them. Null: 20 random sets of the same size through the same walk. Pass: above the "
        "null's 95th percentile by at least 0.05 (an F1 margin; random sets score about 0.04)."),
    params=(GENES, EXCLUDE, NEIGHBORS, MIN_DIST, MIN_CLUSTER, SELECTIONS, SAMPLE,
            Param("features", "choice", "Feature sets",
                  "'all' uses every permitted measurement in each map; 'families' also tries each "
                  "kind of measurement alone, which can find a list that is coherent in one kind "
                  "of evidence and invisible in the combination.", "all",
                  choices=("all", "families")),
            Param("random_sets", "int", "Random lists for the null",
                  "How many random lists of the same size are scored over the same walk to give the "
                  "search-corrected p-value. More is steadier and costs little, because the maps "
                  "are already built.", 100, lo=10, hi=2000, step=10)),
    runner=_geneset_hunt_run, tester=_geneset_hunt_test, cost="minutes",
    needs=("a gene list",)))


def _neighbour_shares(coords, labels: np.ndarray, sources, queries, k: int) -> tuple:
    """(shares queries x classes, classes): the label mix among each query's k nearest sources."""
    from sklearn.neighbors import NearestNeighbors
    classes = sorted({l for l in labels[sources]}, key=str)
    code = {c: i for i, c in enumerate(classes)}
    kk = int(min(k + 1, len(sources)))
    idx = NearestNeighbors(n_neighbors=kk).fit(coords[sources]).kneighbors(
        coords[queries], return_distance=False)
    shares = np.zeros((len(queries), len(classes)))
    for row, (q, nb) in enumerate(zip(queries, idx)):
        nb = [sources[j] for j in nb if sources[j] != q][:k]
        for j in nb:
            shares[row, code[labels[j]]] += 1
        shares[row] /= max(len(nb), 1)
    return shares, classes


def _atlas(coords, lab: np.ndarray, sources, k: int) -> pd.DataFrame:
    """Per class: AUROC and lift of neighbour share among the labelled genes, leave-one-out."""
    shares, classes = _neighbour_shares(coords, lab, sources, sources, k)
    rows = []
    for j, c in enumerate(classes):
        hit = lab[sources] == c
        if hit.sum() < 5:
            continue
        from sklearn.metrics import average_precision_score
        prevalence = float(hit.mean())
        ap = float(average_precision_score(hit, shares[:, j]))
        rows.append({"category": c, "labelled": int(hit.sum()),
                     "auroc": S.auroc(shares[:, j], hit), "auprc": ap, "prevalence": prevalence,
                     "lift": ap / prevalence if prevalence else float("nan")})
    return pd.DataFrame(rows)


def _atlas_run(ctx, p):
    from . import search
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    coords, pos, labels = ctx.blind_map(target, p["sample"], min_cluster_size=25)
    lab = truth.iloc[pos].to_numpy(dtype=object)
    sources = np.flatnonzero([isinstance(v, str) for v in lab])
    table = _atlas(coords, lab, sources, int(p["k"]))
    _s, f1 = search.score_recovery(labels, pd.Series(lab), min_label=5)
    if len(table) and len(f1):
        table = table.merge(f1[["label", "f1"]].rename(columns={"label": "category",
                                                                "f1": "best_cluster_f1"}),
                            on="category", how="left")
    if len(table):
        table["encoded"] = (table["auroc"] >= 0.7) & (table["lift"] >= 2)
        table = table.sort_values("auroc", ascending=False, kind="stable").reset_index(drop=True)
    enc = table[table["encoded"]]["category"].tolist() if len(table) else []
    summary = (f"One map built without {target} or anything restating it. {len(enc)} of "
               f"{len(table)} categories are encoded (neighbour AUROC >= 0.7 and precision lift >= "
               f"2): {', '.join(map(str, enc[:12]))}{'...' if len(enc) > 12 else ''}. The rest are "
               f"not recoverable from these measurements at this map's scale.")
    return StrategyResult("recoverability_atlas", summary, {"categories": table},
                          labels=expand(labels, pos, ctx.n), coords=coords, positions=pos)


def _atlas_test(ctx, p):
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    coords, pos, _labels = ctx.blind_map(target, min(p["sample"] or ctx.n, 2000),
                                         min_cluster_size=25)
    visible, hidden = S.hide(truth, 0.3, ctx.seed, groups=ctx.groups())
    at = {g: i for i, g in enumerate(pos)}
    hid = np.array([at[g] for g in hidden if g in at], dtype=int)
    k = int(p["k"])

    def score(vis):
        lab = vis.iloc[pos].to_numpy(dtype=object)
        sources = np.flatnonzero([isinstance(v, str) for v in lab])
        if len(sources) < 30 or len(hid) < MIN_CLASS:
            return float("nan"), float("nan")
        atlas = _atlas(coords, lab, sources, k)
        shares, classes = _neighbour_shares(coords, lab, sources, hid, k)
        true_hidden = truth.iloc[pos[hid]].to_numpy(dtype=object)
        held = {}
        for j, c in enumerate(classes):
            hit = true_hidden == c
            if hit.sum() >= 3:
                held[c] = S.auroc(shares[:, j], hit)
        atlas = atlas[atlas["category"].isin(held)].sort_values("auroc", ascending=False)
        if len(atlas) < 2:
            return float("nan"), float("nan")
        top = atlas["category"].iloc[: max(1, len(atlas) // 2)]
        rho = S.spearman(atlas["auroc"], [held[c] for c in atlas["category"]])
        return float(np.mean([held[c] for c in top])), rho

    observed, rho = score(visible)
    rng = ctx.rng(17)
    nulls = []
    for i in range(20):
        ctx.check()
        nulls.append(score(S.shuffled(visible, rng))[0])
    return S.judge("recoverability_atlas",
                   "hidden-gene AUROC of the categories the atlas ranks in its top half",
                   observed, nulls, min_effect=0.05, n_hidden=len(hid),
                   hidden=f"30% of the {target} labels", t0=t0,
                   null_kind="20 atlases and recoveries from shuffled labels",
                   numbers={"rank_agreement": rho})


register(Strategy(
    key="recoverability_atlas", number=3, family=MAP,
    title="Ask which categories the data can rediscover",
    question="Of all the categories of a label, which ones do the measurements actually encode -- "
             "and which would no map, however tuned, ever find?",
    tooltip="Builds one map blind to a label and scores every category of that label by how "
            "strongly its members' map neighbours share it; the result is an atlas of what this "
            "data knows and does not know, checked on hidden genes.",
    explanation=(
        "Before searching for a structure it is worth knowing whether there is one to find. This "
        "strategy builds a single map with the chosen label and its closure withheld, then asks, "
        "category by category, how much the map neighbourhood of a labelled gene is enriched for "
        "its own category (AUROC of the neighbour share, and its precision lift over the "
        "category's prevalence). Categories with high scores are ENCODED in the measurements -- "
        "their members behave alike -- and are worth searching for or predicting. Categories near "
        "0.5 are not: no map built from these data will recover them, and a claimed cluster for "
        "one is noise.\n\n"
        "This is knowledge about the data, not about genes: it says which biological distinctions "
        "leave a trace in expression, fitness, modification and structure, and which live only in "
        "the experiment that defined them. The self-test checks that the atlas generalises: the "
        "categories it ranks highest from known genes must be the ones best recovered for genes "
        "whose labels were hidden."),
    walkthrough=(
        "Choose the label to map (compartment, stage, phenotype class...).",
        "Keep 'Genes per map' at 3000 for a quick atlas, or 0 for the whole proteome.",
        "Press Test to confirm the atlas generalises to hidden genes on this label.",
        "Press Run and read 'categories', sorted by AUROC; 'encoded' marks the ones worth pursuing.",
        "Take an encoded category to strategy 01 or 08 for predictions; do not spend a search on "
        "a category the atlas calls not encoded."),
    test_description=(
        "30% of the label is hidden. The atlas is built from visible genes only (leave-one-out "
        "neighbour AUROC per category); hidden genes are then scored by their visible neighbours. "
        "Metric: mean hidden-gene AUROC over the categories the atlas ranks in its top half. Null: "
        "the same with the visible labels shuffled, 20 times. Pass: above the null's 95th "
        "percentile by 0.05. The rank agreement between atlas and hidden recovery is reported."),
    params=(TARGET, SAMPLE, K),
    runner=_atlas_run, tester=_atlas_test, cost="a minute", needs=("a categorical column",)))


def _coassociation(configs, rows) -> np.ndarray:
    """Share of clusterings in which each pair of `rows` was placed in the same cluster."""
    at = {g: i for i, g in enumerate(rows)}
    C = np.zeros((len(rows), len(rows)), dtype=np.float32)
    for c in configs:
        idx = np.array([at[g] for g in c["positions"]])
        for k in np.unique(c["labels"]):
            if k == NOISE:
                continue
            m = idx[c["labels"] == k]
            C[np.ix_(m, m)] += 1.0
    C /= max(len(configs), 1)
    np.fill_diagonal(C, 1.0)
    return C


def _modules(C: np.ndarray, threshold: float, min_size: int) -> np.ndarray:
    """Average-linkage modules of the co-association matrix, cut where pairs co-cluster >= threshold."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    D = 1.0 - C.astype(float)
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    lab = fcluster(Z, t=1.0 - threshold, criterion="distance")
    sizes = pd.Series(lab).value_counts()
    small = sizes[sizes < min_size].index
    lab = np.where(np.isin(lab, small), NOISE, lab)
    return lab.astype(int)


def _consensus_run(ctx, p):
    target = p.get("target")
    rows = ctx.sample_rows(min(p["sample"] or ctx.n, 4000), target)
    configs = _walk(ctx, target, {"all permitted": tuple(ctx.blocks(target))},
                    parse_grid(p["n_neighbors"], int), parse_grid(p["min_dist"]),
                    parse_grid(p["min_cluster_size"], int), rows, parse_grid(p["selection"], str))
    ctx.say("co-association over the walk")
    C = _coassociation(configs, rows)
    modules = _modules(C, float(p["threshold"]), int(p["min_module"]))
    rows_out = []
    truth = ctx.truth(target) if target else None
    for m in sorted(set(modules.tolist()) - {NOISE}):
        idx = np.flatnonzero(modules == m)
        stab = float(C[np.ix_(idx, idx)].mean())
        rec = {"module": m, "size": len(idx), "stability": stab,
               "genes": ", ".join(ctx.gene_ids[rows[idx]][:12])}
        if truth is not None:
            lab = truth.iloc[rows[idx]].dropna()
            if len(lab):
                top = lab.value_counts()
                rec.update({"labelled": len(lab), "top_label": top.index[0],
                            "top_share": float(top.iloc[0] / len(lab))})
        rows_out.append(rec)
    table = pd.DataFrame(rows_out)
    if len(table):
        table = table.sort_values("stability", ascending=False, kind="stable").reset_index(drop=True)
    members = _genes_table(ctx, rows[modules != NOISE], module=modules[modules != NOISE])
    summary = (f"{len(configs)} clusterings over {len(rows):,} genes; {len(table)} modules whose "
               f"members co-cluster in at least {float(p['threshold']):.0%} of them "
               f"({int((modules != NOISE).sum()):,} genes). A module is what survives the choice of "
               f"settings -- the part of the structure no single map's luck produced.")
    return StrategyResult("consensus_modules", summary, {"modules": table, "members": members},
                          labels=expand(modules, rows, ctx.n), coords=configs[0]["coords"],
                          positions=configs[0]["positions"])


def _consensus_test(ctx, p):
    target = _need(ctx, p.get("target") or default_category(ctx))
    rows = ctx.sample_rows(1500, target)
    configs = _walk(ctx, target, {"all permitted": tuple(ctx.blocks(target))},
                    parse_grid(p["n_neighbors"], int)[:2], parse_grid(p["min_dist"])[:1],
                    parse_grid(p["min_cluster_size"], int)[:2], rows,
                    parse_grid(p["selection"], str)[:2])
    modules = _modules(_coassociation(configs, rows), float(p["threshold"]), int(p["min_module"]))
    t0 = time.monotonic()
    truth = ctx.truth(target)
    visible, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    hidden = hidden[np.isin(hidden, rows)]
    score = lambda vis: _hidden_f1(modules, rows, _label_clusters(
        modules, vis.iloc[rows].reset_index(drop=True)), truth, hidden)
    observed = score(visible)
    rng = ctx.rng(19)
    nulls = [score(S.shuffled(visible, rng)) for _ in range(20)]
    return S.judge("consensus_modules",
                   "F1 of hidden genes in the module chosen for their label on known genes",
                   observed, nulls, min_effect=0.05, n_hidden=len(hidden),
                   hidden=f"25% of the {target} labels, whole orthogroups at a time",
                   null_kind="20 label-to-module choices made on shuffled labels", t0=t0,
                   numbers={"modules": int(len(set(modules.tolist()) - {NOISE}))})


register(Strategy(
    key="consensus_modules", number=4, family=MAP,
    title="Keep only the modules that survive the whole walk",
    question="Which groups of genes stay together whatever map settings are chosen -- the "
             "structure that is in the data rather than in one lucky configuration?",
    tooltip="Clusters many maps, counts how often each pair of genes lands in the same cluster, and "
            "keeps the groups that co-cluster in most of them; stable modules are claims about "
            "the genes, unstable ones are claims about a setting.",
    explanation=(
        "Any single map is one choice among many. Change n_neighbors from 15 to 50 and some "
        "clusters survive while others dissolve or merge. This strategy treats that as the "
        "signal: it walks the grid, clusters every map, and builds a co-association matrix -- "
        "the share of clusterings in which each pair of genes was placed together. Modules are "
        "groups whose members co-cluster in at least the chosen share of maps (average linkage "
        "on one minus co-association).\n\n"
        "A module that survives is robust to the arbitrary choices a map requires, which is the "
        "minimum a structure must be before it is interpreted. It is also label-free: no held-out "
        "label is used to build it, so any label it explains afterwards is a genuine finding. The "
        "self-test predicts hidden labels from each module's visible majority. Choosing a "
        "held-out label also removes that label's closure from the maps, so modules can then be "
        "read against it honestly."),
    walkthrough=(
        "Optionally choose a label to hold out; its closure is removed from every map, and the "
        "modules table then reports each module's dominant label.",
        "Widen the grids (three n_neighbors, two min_dist, two cluster sizes is a good start).",
        "Set the co-clustering threshold: 0.5 keeps pairs together in most maps; 0.8 keeps only "
        "the rock-solid core.",
        "Press Test, then Run. The map shows modules; 'modules' lists stability and dominant label.",
        "Treat a module with high stability and no dominant known label as a candidate new unit."),
    test_description=(
        "Modules are built from a small walk on 1,500 genes without the label. 25% of the label is "
        "hidden; the module that best isolates each label is chosen on the visible genes. Metric: "
        "the size-weighted F1 of each label's hidden genes against its chosen module. Null: 20 "
        "choices made on shuffled labels. Pass: above the null's 95th percentile by 0.05."),
    params=(Param("target", "category", "Held-out label (optional)",
                  "A label whose closure is removed from every map, and which the modules are then "
                  "read against. Leave empty to build modules from everything; the self-test then "
                  "uses the default label.", default_category, optional=True),
            Param("n_neighbors", "grid", "UMAP neighbours to try", NEIGHBORS.tip, "15, 30, 50"),
            MIN_DIST, MIN_CLUSTER, SELECTIONS, SAMPLE,
            Param("threshold", "float", "Co-clustering threshold",
                  "The share of clusterings in which two genes must be placed together to be in "
                  "the same module. Higher is stricter: fewer, smaller, more trustworthy modules.",
                  0.5, lo=0.2, hi=0.95, step=0.05),
            Param("min_module", "int", "Smallest module",
                  "Modules smaller than this are dropped as noise. A module of three genes that "
                  "co-cluster is a coincidence as often as it is a complex.", 10, lo=3, hi=500)),
    runner=_consensus_run, tester=_consensus_test, cost="minutes"))


def _associations(labels: np.ndarray, frame: pd.DataFrame, categorical, numeric) -> pd.DataFrame:
    """Every held-out feature against a clustering: chi-square or Kruskal-Wallis, BH-corrected."""
    from scipy.stats import chi2_contingency, kruskal
    labels = np.asarray(labels)
    ok = labels != NOISE
    rows = []
    for c in categorical:
        s = frame[c]
        m = ok & s.notna().to_numpy()
        if m.sum() < 20:
            continue
        tab = pd.crosstab(labels[m], s[m].to_numpy())
        tab = tab.loc[:, tab.sum(axis=0) >= 5]
        tab = tab.loc[tab.sum(axis=1) >= 5]
        if tab.shape[0] < 2 or tab.shape[1] < 2:
            continue
        chi2, pv, _dof, _e = chi2_contingency(tab.to_numpy())
        v = float(np.sqrt(chi2 / (tab.to_numpy().sum() * (min(tab.shape) - 1))))
        rows.append({"feature": c, "kind": "category", "effect": v, "p": float(pv),
                     "n": int(tab.to_numpy().sum())})
    for c in numeric:
        s = pd.to_numeric(frame[c], errors="coerce")
        m = ok & s.notna().to_numpy()
        groups = [s[m & (labels == k)].to_numpy() for k in np.unique(labels[m])]
        groups = [g for g in groups if len(g) >= 5]
        if len(groups) < 2 or sum(len(g) for g in groups) < 20:
            continue
        if np.ptp(np.concatenate(groups)) == 0:  # one value everywhere: nothing to rank
            continue
        h, pv = kruskal(*groups)
        n, kk = sum(len(g) for g in groups), len(groups)
        rows.append({"feature": c, "kind": "number", "effect": float(max(h - kk + 1, 0) / (n - kk)),
                     "p": float(pv), "n": int(n)})
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = S.bh(out["p"].to_numpy())
        out = out.sort_values(["q", "effect"], ascending=[True, False], kind="stable")
    return out.reset_index(drop=True)


def _held_out_features(ctx, used, limit: int = 120) -> tuple:
    """(categorical, numeric) columns the map did not use, for reading it afterwards."""
    cats = [c for c in ctx.categorical_columns() if c not in used]
    nums = [c for c in ctx.numeric_columns() if c not in used]
    nums = sorted(nums, key=lambda c: -int(ctx.values(c).notna().sum()))[:limit]
    return cats, nums


def _family_default(ctx):
    fams = ctx.families()
    for f in ("transcription", "expr", "fitness", "fit"):
        if f in fams:
            return f
    return max(fams, key=lambda f: len(fams[f])) if fams else None


def _tuned_blind(ctx, blocks, rows, target=None):
    """The best-structured of a few maps and clusterings, chosen WITHOUT looking at any label."""
    from . import search
    best, best_q = None, -1.0
    for nn in (15, 30):
        coords, pos = ctx.embed(blocks, nn, 0.1, rows, target=target)
        # Excess-of-mass only: leaf selection scores well on evenness while leaving most of the
        # proteome as noise, and a battery read off a third of the genes finds nothing to replicate.
        for mcs in (15, 30, 60):
            lab = ctx.cluster(coords, mcs, "eom")
            q = search.map_quality(lab)["score"]
            if q > best_q:
                best, best_q = (coords, pos, lab, nn, mcs), q
    return best


def _battery_run(ctx, p):
    fams = ctx.families()
    fam = p["map_from"]
    if fam not in fams:
        raise ValueError(f"no measurement family {fam!r}; there are {', '.join(fams)}")
    blocks = tuple(fams[fam])
    used = {c for b in blocks for c in ctx.blocks().get(b, [])}
    rows = ctx.sample_rows(p["sample"])
    coords, pos, lab, nn, mcs = _tuned_blind(ctx, blocks, rows)
    cats, nums = _held_out_features(ctx, used)
    frame = ctx.nodes.iloc[pos].reset_index(drop=True)
    for c in cats:
        frame[c] = ctx.truth(c).iloc[pos].to_numpy()
    table = _associations(lab, frame, cats, nums)
    sig = table[table["q"] < 0.05] if len(table) else table
    summary = (f"A map built from {fam} alone ({len(used)} columns), tuned for structure without "
               f"any label (n_neighbors {nn}, min cluster {mcs}): {len(sig)} of {len(table)} "
               f"held-out features differ between its clusters at q < 0.05. Strongest: "
               + ", ".join(f"{r.feature} ({r.effect:.2f})" for r in sig.head(8).itertuples())
               + ". What a map built from one kind of evidence separates, in evidence it never saw, "
                 "is what that evidence encodes.")
    return StrategyResult("blind_battery", summary, {"held-out features": table},
                          labels=expand(lab, pos, ctx.n), coords=coords, positions=pos)


def _battery_test(ctx, p):
    fams = ctx.families()
    fam = p["map_from"] if p["map_from"] in fams else _family_default(ctx)
    blocks = tuple(fams[fam])
    used = {c for b in blocks for c in ctx.blocks().get(b, [])}
    rows = ctx.sample_rows(min(p["sample"] or ctx.n, 2000))
    coords, pos, lab, _nn, _mcs = _tuned_blind(ctx, blocks, rows)
    cats, nums = _held_out_features(ctx, used, limit=60)
    frame = ctx.nodes.iloc[pos].reset_index(drop=True)
    for c in cats:
        frame[c] = ctx.truth(c).iloc[pos].to_numpy()

    def find(half):
        sub = np.full(len(lab), NOISE)
        sub[half] = lab[half]
        t = _associations(sub, frame, cats, nums)
        return t[t["q"] < 0.05]["feature"].tolist() if len(t) else []

    def replicate(feature, half, state):
        sub = np.full(len(lab), NOISE)
        sub[half] = lab[half] if state is None else state
        t = _associations(sub, frame, [feature] if feature in cats else [],
                          [feature] if feature in nums else [])
        return bool(len(t) and t["p"].iloc[0] < 0.05)

    scramble = lambda half, rng: rng.permutation(lab[half])
    return S.replication_test(ctx, "blind_battery", find, replicate, scramble, np.arange(len(pos)),
                              n_null=10, min_effect=0.2, label="held-out features the map separates")


register(Strategy(
    key="blind_battery", number=5, family=MAP,
    title="Tune a map without labels, then read what it encodes",
    question="If I build a map from one kind of evidence only -- expression, say -- and tune it for "
             "structure alone, which OTHER measurements do its clusters turn out to separate?",
    tooltip="Builds a map from one family of measurements, tunes it only for cluster structure, "
            "then tests every measurement it never saw against its clusters; replicated "
            "associations say what that kind of evidence encodes about the rest of biology.",
    explanation=(
        "The other strategies start from a question. This one starts from the data and lets the "
        "question come back. Choose one family of measurements -- transcription, fitness, "
        "modification -- and a map is built from that family only, tuned by the project's "
        "map-quality score (how much of the proteome clusters, how evenly) without reference to "
        "any label. Then the battery runs: every categorical and numeric column the map did NOT "
        "use is tested against its clusters (chi-square with Cramér's V, or Kruskal-Wallis with "
        "an eta-squared effect), and the p-values are corrected across the whole family.\n\n"
        "An association here is evidence of coupling between two kinds of biology: clusters of "
        "co-expressed genes that differ in fitness say that transcriptional programmes track "
        "essentiality; clusters that differ in compartment say expression encodes localization. "
        "Because nothing about the held-out features shaped the map, this cannot be circular -- "
        "but it CAN be a large-sample artefact, which is why the self-test demands the "
        "associations replicate on genes the discovery never used."),
    walkthrough=(
        "Choose the family to build the map from ('transcription' is the richest).",
        "Press Test: associations are found on half the genes and checked on the other half.",
        "Press Run and read 'held-out features' by q-value; the effect column says how strong.",
        "Click the map's clusters against the top features (color by the feature in the left "
        "panel) to see what the association looks like.",
        "Repeat with another family and compare: what fitness encodes and what expression encodes "
        "are different answers."),
    test_description=(
        "Pattern 5. One label-free map on up to 2,000 genes. Held-out features significant at "
        "q < 0.05 on a random half of the genes are the findings; each is re-tested (p < 0.05) on "
        "the other half. Metric: share that replicate. Null: the same with the second half's "
        "cluster labels permuted, 10 times. Pass: above the null's 95th percentile by 0.2, with "
        "at least three findings."),
    params=(Param("map_from", "choice", "Build the map from",
                  "The family of measurements the map is built from; everything else becomes a "
                  "held-out feature to read it with. Families are the kinds of evidence the table "
                  "carries: transcription, translation, fitness, modification and so on.",
                  _family_default, choices=lambda ctx: tuple(ctx.families())),
            SAMPLE),
    runner=_battery_run, tester=_battery_test, cost="a minute"))


def _cv_knn(X, truth: pd.Series, groups, k: int = 15, folds: int = 5, seed: int = 0) -> float:
    """Group-aware cross-validated correct-call rate of a k-nearest-neighbour vote."""
    t = pd.Series(truth).reset_index(drop=True)
    known = np.flatnonzero(t.notna().to_numpy())
    if len(known) < 2 * folds:
        return float("nan")
    rng = np.random.default_rng(seed)
    g = np.asarray(groups)[known]
    uniq = rng.permutation(np.unique(g))
    fold_of = {grp: i % folds for i, grp in enumerate(uniq)}
    fold = np.array([fold_of[x] for x in g])
    correct = 0
    for f in range(folds):
        test = known[fold == f]
        vis = t.copy()
        vis.iloc[test] = np.nan
        pred, _s = S.knn_vote(X, vis, k, query=test)
        correct += int((pred.iloc[test].to_numpy(dtype=object) ==
                        t.iloc[test].to_numpy(dtype=object)).sum())
    return correct / len(known)


def _units(ctx, target, unit):
    if unit == "block":
        return ctx.blocks(target)
    fams = ctx.families(target)
    blocks = ctx.blocks(target)
    return {f: [c for b in bs for c in blocks[b]] for f, bs in fams.items()}


def _ablation_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    units = _units(ctx, target, p["unit"])
    everything = [c for cols in units.values() for c in cols]
    groups = ctx.groups()
    k = int(p["k"])
    ctx.say("all permitted measurements")
    base = _cv_knn(ctx.matrix(everything), t, groups, k)
    rows = []
    for i, (u, cols) in enumerate(units.items(), 1):
        ctx.say(f"unit {i} of {len(units)}: {u}")
        rest = [c for c in everything if c not in set(cols)]
        rows.append({"unit": u, "columns": len(cols),
                     "alone": _cv_knn(ctx.matrix(cols), t, groups, k),
                     "without": _cv_knn(ctx.matrix(rest), t, groups, k) if rest else float("nan")})
    table = pd.DataFrame(rows)
    table["loss_when_removed"] = base - table["without"]
    table = table.sort_values("alone", ascending=False, kind="stable").reset_index(drop=True)
    chance = float(((t.value_counts() / t.notna().sum()) ** 2).sum())
    summary = (f"All permitted measurements call {base:.3f} of {target} labels correctly out of "
               f"fold (chance about {chance:.3f}). Alone, the strongest kinds of evidence are "
               + ", ".join(f"{r.unit} ({r.alone:.3f})" for r in table.head(5).itertuples())
               + ". 'loss_when_removed' is what the combination loses without each one -- small "
                 "when another kind of evidence carries the same information.")
    return StrategyResult("block_ablation", summary, {"evidence": table},
                          numbers={"all": base, "chance": chance})


def _ablation_test(ctx, p):
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    visible, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    units = _units(ctx, target, p["unit"])
    k = int(p["k"])
    hidden_acc = {}
    for u, cols in units.items():
        ctx.check()
        pred, _s = S.knn_vote(ctx.matrix(cols), visible, k, query=hidden)
        hidden_acc[u] = S.correct_rate(pred, truth, hidden)
    ranking = {u: _cv_knn(ctx.matrix(cols), visible, ctx.groups(), k)
               for u, cols in units.items()}
    top = max(ranking, key=lambda u: (np.nan_to_num(ranking[u], nan=-1.0), u))
    details = pd.DataFrame({"unit": list(units), "visible_cv": [ranking[u] for u in units],
                            "hidden": [hidden_acc[u] for u in units]}).sort_values(
        "visible_cv", ascending=False)
    return S.judge("block_ablation", f"hidden accuracy of the evidence ranked first ({top})",
                   hidden_acc[top], list(hidden_acc.values()), min_effect=0.02,
                   n_hidden=len(hidden), hidden=f"25% of the {target} labels",
                   null_kind=f"the {len(units)} kinds of evidence chosen at random",
                   t0=t0, details=details, quantile=80.0,
                   numbers={"rank_agreement": S.spearman(details["visible_cv"],
                                                         details["hidden"])})


register(Strategy(
    key="block_ablation", number=6, family=MAP,
    title="Find which kind of evidence carries a label",
    question="Which measurements actually carry the information about this label -- and which "
             "are redundant with others or irrelevant to it?",
    tooltip="Scores each kind of measurement alone and the combination without it, out of fold, "
            "for predicting a held-out label; says which experiments encode which biology and "
            "which could be dropped without loss.",
    explanation=(
        "A map is built from dozens of datasets, and a good result says nothing about which of "
        "them produced it. This strategy asks directly. For each kind of evidence (a family such "
        "as transcription or fitness, or each individual block), it measures how well a "
        "15-nearest-neighbour vote in that evidence alone predicts the held-out label, "
        "cross-validated with whole orthogroups held out together; then how much the combination "
        "loses when that evidence is removed.\n\n"
        "'Alone' measures how much a kind of evidence knows; 'loss when removed' measures how much "
        "it knows that nothing else does. A dataset strong alone but costing nothing when removed "
        "is redundant with another; one weak alone but costly to remove carries something unique. "
        "Both are knowledge about the experiments, and they decide where the next experiment "
        "should go. The self-test checks that the ranking is real: the evidence ranked first on "
        "the visible labels must predict the HIDDEN labels better than a kind of evidence picked "
        "at random."),
    walkthrough=(
        "Choose the label.",
        "Choose 'family' for a quick answer by kind of measurement, 'block' for every dataset "
        "(slower, many more rows).",
        "Press Test, then Run.",
        "Read 'evidence' sorted by 'alone'; then sort by 'loss_when_removed' to find the unique "
        "contributors.",
        "Use the answer to build a focused map in strategy 01 with 'families', or to decide which "
        "dataset a new experiment should extend."),
    test_description=(
        "25% of the label is hidden. Each kind of evidence is ranked by cross-validated accuracy "
        "on the visible labels; the top-ranked one then predicts the hidden labels. Metric: its "
        "hidden accuracy. Null: the hidden accuracy of every kind of evidence, i.e. choosing at "
        "random. Pass: above the null's 80th percentile by 0.02 (with ~15 kinds of evidence the "
        "95th would demand the single best, which asks more than a ranking must deliver)."),
    params=(TARGET, K,
            Param("unit", "choice", "Unit of evidence",
                  "'family' groups the datasets by kind of measurement -- a dozen or so rows, "
                  "quick. 'block' scores every dataset on its own -- a hundred rows, several "
                  "minutes, and the level at which an experiment can be named.", "family",
                  choices=("family", "block"))),
    runner=_ablation_run, tester=_ablation_test, cost="a minute", needs=("a categorical column",)))


# =========================================================================== 2 · neighbours
MIN_SHARE = Param("min_share", "float", "Call when the vote is at least",
                  "The share of the neighbourhood vote the winning label must hold before a gene is "
                  "called. Higher calls fewer genes and fewer of them wrongly; the self-test uses "
                  "the same threshold, so its number describes exactly these calls.", 0.3,
                  lo=0.0, hi=1.0, step=0.05)


def _unlabelled_but_known(vis: pd.Series, truth: pd.Series) -> np.ndarray:
    """The genes a test must call: labelled in truth, hidden from `vis`."""
    return np.flatnonzero(vis.isna().to_numpy() & truth.notna().to_numpy())


def _knn_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, cols = ctx.features(target)
    pred, share = S.knn_vote(X, truth, int(p["k"]), query=np.flatnonzero(truth.isna().to_numpy()))
    pred = pred.where(share >= float(p["min_share"]))
    calls = _calls_table(ctx, pred, share, truth)
    summary = (f"{len(calls):,} unlabelled genes called from their {int(p['k'])} nearest labelled "
               f"genes in {len(cols)} permitted measurements (vote >= {float(p['min_share']):.0%}). "
               f"Distance is taken over every measurement at once, so a gene is called by genes that "
               f"behave like it across the board -- not by any one experiment.")
    return StrategyResult("feature_knn", summary, {"calls": calls}, genes=list(calls["gene_id"][:200]))


def _knn_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    k, floor = int(p["k"]), float(p["min_share"])

    def predict(vis):
        pred, share = S.knn_vote(X, vis, k, query=_unlabelled_but_known(vis, truth))
        return pred.where(share >= floor)

    return S.label_transfer_test(ctx, "feature_knn", target, predict, n_null=10)


register(Strategy(
    key="feature_knn", number=7, family=NEIGHBOURS,
    title="Call a gene by the genes that behave like it",
    question="For a gene with no label, what label do the genes most similar to it across every "
             "permitted measurement carry?",
    tooltip="Finds each unlabelled gene's nearest labelled genes across all permitted measurements "
            "at once and calls it by their distance-weighted vote -- the plainest guilt-by-"
            "association, with no map and no clustering between the data and the answer.",
    explanation=(
        "Every map and every clustering is an approximation of one thing: which genes are similar. "
        "This strategy skips the approximation. Each gene is a point in the space of every "
        "permitted measurement (rank-scaled, so no screen dominates by units); its k nearest "
        "labelled genes vote, weighted by closeness, and the gene is called with the winning label "
        "when that label holds at least the chosen share of the vote.\n\n"
        "It is the baseline every other strategy has to beat: if a sophisticated method cannot do "
        "better than asking the neighbours, the sophistication is decoration. It is also "
        "transparent -- each call can be traced to the genes that made it. Its weakness is the "
        "curse of dimension: with hundreds of measurements, many of them missing and filled at the "
        "median, 'nearest' can be dominated by missingness rather than biology, which is why the "
        "self-test hides whole orthogroups and demands a margin over shuffled labels."),
    walkthrough=(
        "Choose the label.",
        "Leave k at 15; raise the vote threshold to 0.7 when you want fewer, safer calls.",
        "Press Test: the number is the share of hidden genes called correctly, abstentions "
        "counted as misses, against shuffled labels.",
        "Press Run and read 'calls', strongest vote first. Click a row to find the gene.",
        "Compare with strategies 08 and 11: a call that the measurements, the map and a network "
        "all agree on is strategy 31."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; each hidden gene called by its "
        "k nearest visible genes with the same vote threshold. Metric: hidden genes called "
        "correctly. Null: 10 runs with visible labels shuffled. Pass: above the null's 95th "
        "percentile by 0.05."),
    params=(TARGET, K, MIN_SHARE),
    runner=_knn_run, tester=_knn_test, needs=("a categorical column",)))


MAP_NN = Param("n_neighbors", "int", "UMAP neighbours",
               "The UMAP n_neighbors of the map the neighbours are read from. Larger values give a "
               "map whose distances reflect broad organisation; smaller ones, local detail.", 25,
               lo=2, hi=200)
MAP_MD = Param("min_dist", "float", "UMAP min_dist",
               "The UMAP min_dist of that map. Small values pack similar genes together, which is "
               "what a neighbourhood vote wants.", 0.1, lo=0.0, hi=1.0, step=0.05)


def _map_knn(ctx, p, size):
    target = _need(ctx, p["target"])
    coords, pos, labels = ctx.blind_map(target, size, int(p["n_neighbors"]), float(p["min_dist"]))
    return target, coords, pos, labels


def _map_nn_run(ctx, p):
    target, coords, pos, labels = _map_knn(ctx, p, p["sample"])
    truth = ctx.truth(target)
    sub = truth.iloc[pos].reset_index(drop=True)
    pr, sh = S.knn_vote(coords, sub, int(p["k"]), query=np.flatnonzero(sub.isna().to_numpy()))
    pred = pd.Series([np.nan] * ctx.n, dtype=object)
    share = pd.Series(np.nan, index=range(ctx.n))
    pred.iloc[pos] = pr.where(sh >= float(p["min_share"])).to_numpy()
    share.iloc[pos] = sh.to_numpy()
    calls = _calls_table(ctx, pred, share, truth)
    summary = (f"One map of {len(pos):,} genes built blind to {target}; {len(calls):,} unlabelled "
               f"genes called by their {int(p['k'])} nearest labelled neighbours in it. The map "
               f"denoises: distances in three dimensions are what survived UMAP's compression.")
    return StrategyResult("map_neighbours", summary, {"calls": calls},
                          labels=expand(labels, pos, ctx.n), coords=coords, positions=pos,
                          genes=list(calls["gene_id"][:200]))


def _map_nn_test(ctx, p):
    target, coords, pos, _labels = _map_knn(ctx, p, min(p["sample"] or ctx.n, 2500))
    truth = ctx.truth(target)
    k, floor = int(p["k"]), float(p["min_share"])
    inside = np.zeros(ctx.n, bool)
    inside[pos] = True
    tsub = truth.iloc[pos].reset_index(drop=True)

    def predict(vis):
        sub = vis.iloc[pos].reset_index(drop=True)
        pr, sh = S.knn_vote(coords, sub, k, query=_unlabelled_but_known(sub, tsub))
        out = pd.Series([np.nan] * ctx.n, dtype=object)
        out.iloc[pos] = pr.where(sh >= floor).to_numpy()
        return out

    return S.label_transfer_test(ctx, "map_neighbours", target, predict, restrict=inside,
                                 n_null=20)


register(Strategy(
    key="map_neighbours", number=8, family=NEIGHBOURS,
    title="Call a gene by its neighbours on the map",
    question="On a map built without the label, which label do a gene's nearest placed neighbours "
             "carry?",
    tooltip="Builds one map blind to the label and calls each unlabelled gene by the vote of its "
            "nearest labelled neighbours in the map's three dimensions -- what you do by eye when "
            "you hover over an unlabelled point, made systematic and scored.",
    explanation=(
        "When you look at the map colored by compartment and see a grey point inside a blue "
        "cloud, you are making this inference. The strategy makes it for every gene: a map is "
        "built with the label and its closure withheld, and each unlabelled gene is called by the "
        "distance-weighted vote of its k nearest labelled genes in the map.\n\n"
        "The difference from strategy 07 is the map. UMAP keeps local neighbourhoods and discards "
        "much of the noise in hundreds of partly missing measurements, so map neighbours can be "
        "more reliable than raw ones -- or less, when the projection tears a category apart. Which "
        "is true for a given label is exactly what running both self-tests tells you. The map is "
        "a sample unless 'Genes per map' is zero, and genes outside it are not called at all."),
    walkthrough=(
        "Choose the label; keep the map settings unless the atlas (03) suggested others.",
        "Press Test; compare its number with strategy 07's on the same label.",
        "Press Run: the map is drawn with its clusters, and 'calls' lists the calls.",
        "Color the map by the label to see each call in context; click a call to find it."),
    test_description=(
        "Pattern 1 on the genes the map places (up to 2,500): 25% of the label hidden by whole "
        "orthogroups; hidden genes called by their k nearest visible genes in the map. Null: 20 "
        "runs on shuffled visible labels. Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET, Param("sample", "int", SAMPLE.label, SAMPLE.tip, 4000, lo=0, hi=20000,
                          step=500), MAP_NN, MAP_MD, K, MIN_SHARE),
    runner=_map_nn_run, tester=_map_nn_test, cost="a minute", needs=("a categorical column",)))


def _enriched(labels: np.ndarray, vis: pd.Series, min_lift: float, q_max: float = 0.05,
              min_hits: int = 3) -> tuple:
    """(prediction over the clustered genes, the enrichment table) for guilt by cluster."""
    labels = np.asarray(labels)
    v = pd.Series(vis).reset_index(drop=True)
    known = v.notna().to_numpy() & (labels != NOISE)
    total = int(known.sum())
    counts = v[known].value_counts()
    rows = []
    for k in np.unique(labels[labels != NOISE]):
        in_k = labels == k
        draws = int((in_k & known).sum())
        if draws < 3:
            continue
        for c, hits in v[in_k & known].value_counts().items():
            if hits < min_hits:
                continue
            lift = (hits / draws) / (counts[c] / total)
            rows.append({"cluster": int(k), "label": c, "hits": int(hits), "labelled": draws,
                         "size": int(in_k.sum()), "lift": float(lift),
                         "p": S.hypergeom_sf(int(hits), draws, int(counts[c]), total)})
    table = pd.DataFrame(rows)
    pred = pd.Series([np.nan] * len(labels), dtype=object)
    if table.empty:
        return pred, table
    table["q"] = S.bh(table["p"].to_numpy())
    good = table[(table["q"] < q_max) & (table["lift"] >= min_lift)].sort_values(
        ["q", "lift"], ascending=[True, False], kind="stable")
    chosen = good.drop_duplicates("cluster")
    table["called"] = table.index.isin(chosen.index)
    for r in chosen.itertuples():
        pred[labels == r.cluster] = r.label
    return pred, table.sort_values("q", kind="stable").reset_index(drop=True)


def _guilt_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    coords, pos, labels = ctx.blind_map(target, p["sample"],
                                        min_cluster_size=int(p["min_cluster_size"]),
                                        selection=p["selection"])
    sub = truth.iloc[pos].reset_index(drop=True)
    pr, table = _enriched(labels, sub, float(p["min_lift"]))
    pred = pd.Series([np.nan] * ctx.n, dtype=object)
    pred.iloc[pos] = pr.to_numpy()
    cl = expand(labels, pos, ctx.n)
    calls = _calls_table(ctx, pred, pd.Series(np.nan, index=range(ctx.n)), truth, cluster=cl)
    if len(calls) and len(table):
        lift = table[table["called"]].set_index("cluster")["lift"]
        calls["support"] = calls["cluster"].map(lift).astype(float)
        calls = calls.sort_values("support", ascending=False, kind="stable").reset_index(drop=True)
    n_called = int(table["called"].sum()) if len(table) else 0
    summary = (f"{n_called} clusters are enriched for one {target} label at q < 0.05 and lift >= "
               f"{float(p['min_lift']):g} (hypergeometric, BH-corrected over every cluster x label); "
               f"their {len(calls):,} unlabelled members are called with it.")
    return StrategyResult("cluster_guilt", summary, {"calls": calls, "enrichment": table},
                          labels=cl, coords=coords, positions=pos, genes=list(calls["gene_id"][:200]))


def _guilt_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    _coords, pos, labels = ctx.blind_map(target, min(p["sample"] or ctx.n, 2500),
                                         min_cluster_size=int(p["min_cluster_size"]),
                                         selection=p["selection"])
    inside = np.zeros(ctx.n, bool)
    inside[pos] = True

    def predict(vis):
        pr, _t = _enriched(labels, vis.iloc[pos].reset_index(drop=True), float(p["min_lift"]))
        out = pd.Series([np.nan] * ctx.n, dtype=object)
        out.iloc[pos] = pr.to_numpy()
        return out

    return S.label_transfer_test(ctx, "cluster_guilt", target, predict, restrict=inside, n_null=20,
                                 precision=True, min_effect=0.1)


register(Strategy(
    key="cluster_guilt", number=9, family=NEIGHBOURS,
    title="Name a cluster by the label it is enriched for",
    question="Which clusters of a label-blind map hold one label far more often than chance, and "
             "what does that make of their unlabelled members?",
    tooltip="Clusters a label-blind map, tests every cluster against every label for enrichment "
            "(hypergeometric, BH-corrected), and calls the unlabelled members of significantly "
            "enriched clusters -- guilt by cluster, with a p-value rather than a majority.",
    explanation=(
        "A majority is not evidence: a cluster of five genes, three of them nuclear, is 60% "
        "nuclear by chance in a proteome where a fifth of labelled genes are. This strategy asks "
        "the statistical question instead. For each cluster and each label it computes the "
        "hypergeometric probability of seeing that many of the label among the cluster's labelled "
        "members, corrects across every cluster-label pair tested, and calls a cluster only when "
        "the label is both significant and enriched by at least the chosen lift over its share "
        "of the proteome.\n\n"
        "Enrichment is also how a cluster can be NAMED without being pure: a cluster that is 35% "
        "apicoplast when apicoplast is 3% of the proteome is an apicoplast-flavoured cluster, and "
        "its unlabelled members deserve a look even though most labelled members are elsewhere. "
        "That is also its risk, which is why the self-test scores exactly these calls on hidden "
        "genes."),
    walkthrough=(
        "Choose the label and a cluster size (25 is the project default).",
        "Press Test, then Run.",
        "Read 'enrichment': every cluster-label pair with its lift and q-value, called ones marked.",
        "Read 'calls' for the unlabelled members of called clusters; sort by 'support' (the lift)."),
    test_description=(
        "Pattern 1 scored by precision, on the genes a blind map places (up to 2,500): 25% of the "
        "label hidden; the enrichment is recomputed from visible labels, and the calls it makes "
        "on hidden genes are scored -- the strategy abstains on noise and unenriched clusters by "
        "design, so what matters is how often a call is right. Null: 20 runs on shuffled visible "
        "labels. Pass: above the null's 95th percentile by 0.1."),
    params=(TARGET, SAMPLE,
            Param("min_cluster_size", "int", "Cluster size",
                  "HDBSCAN's minimum cluster size for the one map this strategy reads. Enrichment "
                  "needs clusters big enough to hold several labelled genes, so very small values "
                  "call almost nothing.", 25, lo=5, hi=500),
            SELECTION,
            Param("min_lift", "float", "Minimum lift",
                  "How many times its proteome-wide share a label must be over-represented in a "
                  "cluster before the cluster is named for it. Two is conservative for large "
                  "classes and generous for rare ones.", 2.0, lo=1.0, hi=20.0, step=0.5)),
    runner=_guilt_run, tester=_guilt_test, cost="a minute", needs=("a categorical column",)))


def _surprise(ctx, X, labels: pd.Series, k: int, target) -> pd.DataFrame:
    """Per labelled gene: how much its neighbourhood, in measurements and networks, shares its label."""
    lab = pd.Series(labels).reset_index(drop=True).to_numpy(dtype=object)
    known = np.flatnonzero([isinstance(v, str) for v in lab])
    shares, classes = _neighbour_shares(X, lab, known, known, k)
    code = {c: i for i, c in enumerate(classes)}
    own = shares[np.arange(len(known)), [code[lab[g]] for g in known]]
    alt = shares.copy()
    alt[np.arange(len(known)), [code[lab[g]] for g in known]] = -1
    out = pd.DataFrame({"position": known, "label": lab[known], "own_share": own,
                        "alternative": [classes[j] for j in alt.argmax(axis=1)],
                        "alternative_share": alt.max(axis=1)})
    layers = ctx.measurement_layers(target)
    if layers:
        A = sum(ctx.adjacency(l, "binary") for l in layers)
        vis = pd.Series(lab)
        M, cls = S.onehot(vis)
        counts = np.asarray((A[known] @ M).todense())
        tot = counts.sum(axis=1)
        idx = {c: i for i, c in enumerate(cls)}
        net = np.where(tot > 0, counts[np.arange(len(known)), [idx[lab[g]] for g in known]]
                       / np.maximum(tot, 1), np.nan)
        out["network_share"] = net
    both = out[["own_share"] + (["network_share"] if "network_share" in out else [])]
    out["surprise"] = 1.0 - both.mean(axis=1, skipna=True)
    return out


def _outliers_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    table = _surprise(ctx, X, t, int(p["k"]), target)
    table.insert(0, "gene_id", ctx.gene_ids[table["position"]])
    table.insert(1, "product", ctx.product(table["position"]))
    table = table.drop(columns="position").sort_values("surprise", ascending=False,
                                                       kind="stable").reset_index(drop=True)
    top = table.head(int(p["top"]))
    summary = (f"{len(table):,} labelled genes scored by how little their neighbourhood shares their "
               f"{target} label; the {len(top)} most surprising are listed with the label their "
               f"neighbours carry instead. Each is a mislabel, a dual or moonlighting protein, or a "
               f"gene whose measurements are unusual -- all three worth knowing.")
    return StrategyResult("label_outliers", summary, {"surprising genes": top},
                          genes=list(top["gene_id"]))


def _outliers_test(ctx, p):
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    X, _cols = ctx.features(target)
    rng = ctx.rng(21)
    known = np.flatnonzero(t.notna().to_numpy())
    n_bad = max(10, int(0.05 * len(known)))
    if len(known) < 3 * n_bad or t.nunique() < 2:
        return S.judge("label_outliers", "AUROC for the corrupted labels", float("nan"), [],
                       min_effect=0.1, n_hidden=0, hidden="corrupted labels",
                       null_kind="random sets", t0=t0, note="too few labelled genes to corrupt")
    bad = rng.choice(known, size=n_bad, replace=False)
    freq = t.value_counts(normalize=True)
    corrupted = t.copy()
    for g in bad:
        others = freq.drop(t.iloc[g])
        corrupted.iloc[g] = rng.choice(others.index.to_numpy(), p=(others / others.sum()).to_numpy())
    table = _surprise(ctx, X, corrupted, int(p["k"]), target)
    is_bad = np.isin(table["position"].to_numpy(), bad)
    observed = S.auroc(table["surprise"], is_bad)
    nulls = [S.auroc(table["surprise"], np.isin(table["position"].to_numpy(),
                                                 rng.choice(known, size=n_bad, replace=False)))
             for _ in range(20)]
    return S.judge("label_outliers", "AUROC of surprise for the swapped labels", observed, nulls,
                   min_effect=0.1, n_hidden=n_bad,
                   hidden=f"{n_bad} labels (5%) swapped to a wrong class, in proportion to class size",
                   null_kind="20 random sets of the same size", t0=t0)


register(Strategy(
    key="label_outliers", number=10, family=NEIGHBOURS,
    title="Find genes whose label their neighbours contradict",
    question="Which labelled genes sit among genes that almost all carry a different label -- "
             "possible mislabels, dual-localized or moonlighting proteins?",
    tooltip="Scores every labelled gene by how little its measurement neighbours and network "
            "partners share its label; the top of the list is where an annotation is wrong, a "
            "protein has two lives, or the measurements are strange.",
    explanation=(
        "Every other strategy trusts the labels and predicts the missing ones. This one turns the "
        "question round: given everything else, which EXISTING labels look wrong? For each "
        "labelled gene it measures the share of its k nearest labelled genes (in the permitted "
        "measurements) that carry its own label, and the share of its labelled partners across the "
        "measured networks that do; surprise is one minus their mean.\n\n"
        "A high-surprise gene is one of three things, and each is worth knowing: an annotation "
        "error (hyperLOPIT assigns by a classifier and is wrong some fraction of the time); a "
        "protein with two locations or functions, which the neighbourhood sees as the other one; "
        "or a gene whose measurements are unusual, which is itself a finding. The self-test "
        "corrupts 5% of the labels on purpose and asks whether those genes rise to the top."),
    walkthrough=(
        "Choose the label to audit.",
        "Press Test: the AUROC says how well the method finds labels that were deliberately "
        "swapped; 0.5 would be no better than a random list.",
        "Press Run; read 'surprising genes' from the top, with the label their neighbours carry "
        "instead in 'alternative'.",
        "Check the top genes' evidence panel before believing either label."),
    test_description=(
        "5% of the labels (at least ten) are swapped to a wrong class, drawn in proportion to "
        "class size. Surprise is computed with the corrupted labels. Metric: AUROC of surprise for "
        "the swapped genes among all labelled genes. Null: 20 random sets of the same size. Pass: "
        "above the null's 95th percentile by 0.1."),
    params=(TARGET, K,
            Param("top", "int", "Genes to list",
                  "How many of the most surprising genes to list. The ranking is over every "
                  "labelled gene; this only limits how much of it is shown and saved.", 200,
                  lo=10, hi=5000, step=10)),
    runner=_outliers_run, tester=_outliers_test, needs=("a categorical column",)))


# =========================================================================== 3 · networks
LAYER = Param("layer", "layer", "Edge layer",
              "The measured network the label is walked across. Each layer is a different kind of "
              "evidence -- co-expression, co-fitness, crosslinks, structural similarity -- and a layer "
              "built from the held-out label itself is refused.",
              _layer_default(["coexpression", "cofitness", "cotranslation", "xlms", "struct",
                              "ip_ms"]))


def _usable_layer(ctx, layer, target):
    if not layer or layer not in ctx.layers():
        raise ValueError(f"this table has no edge layer {layer!r}; it has {', '.join(ctx.layers())}")
    if layer in ctx.banned_layers(target):
        raise ValueError(f"the {layer} layer is built from {target} or its family, so walking it "
                         f"would read the answer back; choose another layer")
    return layer


def _propagation_run(ctx, p):
    target = _need(ctx, p["target"])
    layer = _usable_layer(ctx, p["layer"], target)
    truth = ctx.truth(target)
    pred, strength = S.propagate(ctx.operator(layer), truth, np.flatnonzero(truth.isna().to_numpy()),
                                 restart=float(p["restart"]))
    calls = _calls_table(ctx, pred, strength, truth)
    reached = int(np.asarray(ctx.adjacency(layer).sum(axis=1)).ravel().astype(bool).sum())
    summary = (f"The {target} labels diffused across {layer} (random walk with restart "
               f"{float(p['restart']):g}); {len(calls):,} unlabelled genes reached. The layer touches "
               f"{reached:,} genes; the rest cannot be called from it at all, which is an absence of "
               f"evidence rather than evidence of anything.")
    return StrategyResult("layer_propagation", summary, {"calls": calls},
                          genes=list(calls["gene_id"][:200]))


def _propagation_test(ctx, p):
    target = _need(ctx, p["target"])
    layer = _usable_layer(ctx, p["layer"], target)
    truth = ctx.truth(target)
    op = ctx.operator(layer)
    predict = lambda vis: S.propagate(op, vis, _unlabelled_but_known(vis, truth),
                                      restart=float(p["restart"]))[0]
    return S.label_transfer_test(ctx, "layer_propagation", target, predict, n_null=10)


register(Strategy(
    key="layer_propagation", number=11, family=NETWORKS,
    title="Diffuse a label across one measured network",
    question="If labels flow along the edges of one kind of measured relationship, where do they "
             "end up -- and how much of a label does that relationship carry?",
    tooltip="Seeds each label on the genes that carry it and lets it diffuse along one edge layer "
            "by random walk with restart; genes are called by the label whose field reaches them "
            "most strongly. Also measures how much one kind of relationship knows.",
    explanation=(
        "A network carries a label if genes linked in it tend to share the label. Diffusion "
        "exploits that without needing a cluster to form: every class is seeded on its labelled "
        "genes, the seed spreads along the layer's edges (degree-normalised, so hubs -- often the "
        "most-studied genes -- do not swallow everything), and after the walk settles each gene "
        "is called by the strongest field. Genes the layer does not reach are not called.\n\n"
        "Run it layer by layer and the self-test numbers become a table of which relationships "
        "encode which biology: crosslinks should carry compartment strongly, co-fitness should "
        "carry pathway, structural similarity should carry function. A layer built from the label "
        "(the compartment layer, for compartment) is refused -- it would recover the label "
        "perfectly and mean nothing."),
    walkthrough=(
        "Choose the label and a layer.",
        "Press Test for each layer you are curious about and note the numbers: that is how much "
        "each kind of relationship knows about the label.",
        "Press Run on the best layer and read 'calls'.",
        "Lower the restart to let labels travel further (more calls, less certain)."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; the fields are seeded from "
        "visible genes only, so a hidden gene never seeds its own call. Metric: hidden genes "
        "called correctly (unreached genes count as misses). Null: 10 runs on shuffled labels. "
        "Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET, LAYER,
            Param("restart", "float", "Restart probability",
                  "At each step the walk returns to the seeds with this probability. High values "
                  "keep labels close to where they started; low values let them travel several "
                  "edges, which reaches more genes with weaker evidence.", 0.5,
                  lo=0.05, hi=0.95, step=0.05)),
    runner=_propagation_run, tester=_propagation_test, needs=("a categorical column",
                                                               "an edge layer")))


def _sources(ctx, target, k: int) -> dict:
    """name -> predict(visible, query): every permitted layer's vote, and the measurements' kNN."""
    out = {}
    for layer in ctx.measurement_layers(target):
        A = ctx.adjacency(layer)
        out[layer] = lambda vis, q, A=A: S.graph_vote(A, vis, q)[0]
    X, _cols = ctx.features(target)
    out["measurements (kNN)"] = lambda vis, q: S.knn_vote(X, vis, k, q)[0]
    return out


def _weighted_vote(ctx, sources: dict, vis: pd.Series, query) -> tuple:
    """(prediction, support, weights): each source weighted by its skill on an inner holdout."""
    inner_vis, inner = S.hide(vis, 0.2, ctx.seed + 1, groups=ctx.groups(), min_class=5)
    weights = {}
    for name, fn in sources.items():
        ctx.check()
        p_in = fn(inner_vis, inner)
        acc = S.correct_rate(p_in, vis, inner)
        weights[name] = (max(acc - S.analytic_null(p_in, vis, inner)[0], 0.0)
                         if np.isfinite(acc) else 0.0)
    query = np.asarray(query, dtype=int)
    tally: dict = {}
    for name, fn in sources.items():
        if weights[name] <= 0:
            continue
        pr = fn(vis, query)
        for g in query:
            v = pr.iloc[g]
            if isinstance(v, str):
                tally.setdefault(g, {}).setdefault(v, 0.0)
                tally[g][v] += weights[name]
    pred = pd.Series([np.nan] * len(vis), dtype=object)
    support = pd.Series(np.nan, index=range(len(vis)))
    for g, votes in tally.items():
        best = max(votes, key=votes.get)
        pred.iloc[g] = best
        support.iloc[g] = votes[best] / sum(votes.values())
    return pred, support, weights


def _vote_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    sources = _sources(ctx, target, int(p["k"]))
    pred, support, w = _weighted_vote(ctx, sources, truth, np.flatnonzero(truth.isna().to_numpy()))
    calls = _calls_table(ctx, pred, support, truth)
    weights = pd.DataFrame({"source": list(w), "weight": list(w.values())}).sort_values(
        "weight", ascending=False).reset_index(drop=True)
    summary = (f"{len(sources)} sources vote, each weighted by how far it beat chance on an inner "
               f"holdout of the known labels; {int((weights['weight'] > 0).sum())} earned a vote. "
               f"{len(calls):,} unlabelled genes called. Strongest: "
               + ", ".join(f"{r.source} ({r.weight:.2f})" for r in weights.head(4).itertuples()))
    return StrategyResult("layer_vote", summary, {"calls": calls, "source weights": weights},
                          genes=list(calls["gene_id"][:200]))


def _vote_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    sources = _sources(ctx, target, int(p["k"]))
    predict = lambda vis: _weighted_vote(ctx, sources, vis, _unlabelled_but_known(vis, truth))[0]
    return S.label_transfer_test(ctx, "layer_vote", target, predict, n_null=5)


register(Strategy(
    key="layer_vote", number=12, family=NETWORKS,
    title="Let every network vote, weighted by what it has earned",
    question="If every measured relationship and the measurements themselves vote on a gene's "
             "label, each weighted by how good it has proven to be, what is the verdict?",
    tooltip="Each permitted edge layer and the measurement-space neighbours vote on every gene; "
            "each source's vote is weighted by how far it beat chance on an inner holdout of the "
            "known labels, so the evidence that has earned trust counts most.",
    explanation=(
        "No single network reaches every gene, and no single network is right about every label. "
        "This strategy asks all of them. Each permitted layer contributes its neighbour vote and "
        "the measurement space contributes its k-nearest-neighbour vote; before voting, each "
        "source is scored on an inner holdout carved from the known labels, and its weight is its "
        "accuracy MINUS its own chance level -- a source that only guesses the commonest class "
        "earns nothing.\n\n"
        "The result reaches more genes than any one layer and is steered by the sources that know "
        "the label. The 'source weights' table is itself a finding: it ranks the kinds of "
        "evidence by how much they know about this label, net of chance -- and a layer with zero "
        "weight is one this label is invisible in, however many edges it has."),
    walkthrough=(
        "Choose the label.",
        "Press Test (five shuffled-label repeats; a minute or two on the full table).",
        "Press Run and read 'source weights' first: zero weight means a source knew nothing "
        "beyond chance about this label.",
        "Read 'calls'; support is the winning label's share of the weighted vote."),
    test_description=(
        "Pattern 1: 25% of the label hidden; the weights are learned on an inner holdout of the "
        "visible labels only, then hidden genes are called. Null: 5 runs with visible labels "
        "shuffled (weights relearned each time). Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET, K),
    runner=_vote_run, tester=_vote_test, cost="a minute", needs=("a categorical column",)))


def _physical(ctx, target):
    layers = [l for l in ("xlms", "ip_ms") if l in ctx.layers() and l not in ctx.banned_layers(
        target)]
    if not layers:
        return layers, None
    return layers, sum(ctx.adjacency(l, "w") for l in layers)


def _physical_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    layers, A = _physical(ctx, target)
    if A is None:
        return StrategyResult("physical_partners", "This table has no crosslink or IP-MS layer.", {})
    degree = np.asarray((A > 0).sum(axis=1)).ravel()
    query = np.flatnonzero(truth.isna().to_numpy() & (degree > 0))
    pred, support = S.graph_vote(A, truth, query)
    evidence = []
    for g in range(ctx.n):
        if not isinstance(pred.iloc[g], str):
            evidence.append("")
            continue
        row = A.getrow(g)
        partners = [(j, w) for j, w in zip(row.indices, row.data) if isinstance(truth.iloc[j], str)]
        partners.sort(key=lambda x: -x[1])
        evidence.append("; ".join(f"{ctx.gene_ids[j]} ({truth.iloc[j]}, {w:g})"
                                  for j, w in partners[:4]))
    calls = _calls_table(ctx, pred, support, truth, partners=np.array(evidence, dtype=object))
    summary = (f"{len(calls):,} unlabelled genes with a measured physical partner in "
               f"{' + '.join(layers)} are called by their partners' {target}, weighted by crosslink "
               f"or pulldown count. Each call names the partners that made it.")
    return StrategyResult("physical_partners", summary, {"calls": calls},
                          genes=list(calls["gene_id"][:200]))


def _physical_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    layers, A = _physical(ctx, target)
    if A is None:
        return S.judge("physical_partners", "correct calls per hidden gene", float("nan"), [],
                       min_effect=0.05, n_hidden=0, hidden="labels", null_kind="shuffled labels",
                       t0=time.monotonic(), note="no crosslink or IP-MS layer")
    has = np.asarray((A > 0).sum(axis=1)).ravel() > 0
    predict = lambda vis: S.graph_vote(A, vis, _unlabelled_but_known(vis, truth))[0]
    return S.label_transfer_test(ctx, "physical_partners", target, predict, restrict=has,
                                 n_null=20)


register(Strategy(
    key="physical_partners", number=13, family=NETWORKS,
    title="Place a protein by the proteins it physically touches",
    question="For a protein crosslinked to or pulled down with labelled proteins, what does its "
             "physical company say about where it lives and what it joins?",
    tooltip="Calls genes by the labels of their measured physical partners -- crosslinks and "
            "replicated pulldowns -- weighted by how many times the contact was seen, and names "
            "the partners behind every call so each can be checked.",
    explanation=(
        "Two proteins crosslinked by DSS were within about 30 Å of each other in the parasite: "
        "they share a compartment, and often a complex. That makes a crosslink the most direct "
        "evidence this table has about where an unlocalized protein lives, and unlike co-mention "
        "it is not biased toward famous genes. This strategy calls every protein with at least one "
        "measured physical partner by the weighted label vote of its partners (crosslink or "
        "pulldown counts as weights), and lists the partners so each call can be checked by "
        "hand.\n\n"
        "Its reach is limited to proteins in the interactomes -- a few thousand at most -- and it "
        "inherits the crosslinker's chemistry (lysines, soluble proteins). Within that reach it "
        "is expected to be the most precise localization evidence available, and the self-test "
        "measures exactly that reach."),
    walkthrough=(
        "Choose the label (compartment is the natural one).",
        "Press Test: only hidden genes with a physical partner are scored.",
        "Press Run and read 'calls'; the 'partners' column names each partner, its label and "
        "the number of contacts.",
        "For a call that matters, open the partners' evidence and the StarPath interaction page."),
    test_description=(
        "Pattern 1 restricted to genes with at least one physical partner: 25% of the label "
        "hidden by whole orthogroups; hidden genes called by their visible partners' weighted "
        "vote. Null: 20 runs on shuffled labels. Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET,),
    runner=_physical_run, tester=_physical_test, needs=("a crosslink or IP-MS layer",)))


def _ec_level(truth: pd.Series, level: int) -> pd.Series:
    """An annotation truncated to `level` dot-separated parts: EC 3.1.3.16 at level 1 is '3'."""
    import re

    def cut(v):
        if not isinstance(v, str):
            return np.nan
        tok = re.split(r"[;,]", v)[0].strip()
        m = re.match(r"^(\d+(?:\.[\d\-n]+)*)", tok)
        if level > 0 and m:
            return ".".join(m.group(1).split(".")[:level])
        return tok or np.nan
    return truth.map(cut).astype(object)


def _struct_setup(ctx, p):
    target = _need(ctx, p["target"])
    layer = _usable_layer(ctx, p["layer"], target)
    truth = _ec_level(ctx.truth(target), int(p["level"]))
    return target, layer, truth, ctx.adjacency(layer, "w")


def _struct_run(ctx, p):
    target, layer, truth, A = _struct_setup(ctx, p)
    degree = np.asarray((A > 0).sum(axis=1)).ravel()
    pred, support = S.graph_vote(A, truth, np.flatnonzero(truth.isna().to_numpy() & (degree > 0)))
    products = pd.Series(ctx.product(np.arange(ctx.n)))
    hypothetical = products.str.contains("hypothetical", case=False, na=False).to_numpy()
    calls = _calls_table(ctx, pred, support, truth, hypothetical=hypothetical)
    calls = calls.sort_values(["hypothetical", "support"], ascending=[False, False],
                              kind="stable").reset_index(drop=True)
    summary = (f"{len(calls):,} genes without a {target} annotation (at level {int(p['level'])}) "
               f"called from their {layer} neighbours; {int(calls['hypothetical'].sum()) if len(calls) else 0}"
               f" of them are hypothetical proteins. Structural similarity needs no sequence "
               f"homology, so it reaches the lineage-specific proteins orthology cannot.")
    return StrategyResult("structural_homology", summary, {"calls": calls},
                          genes=list(calls["gene_id"][:200]))


def _struct_test(ctx, p):
    target, layer, truth, A = _struct_setup(ctx, p)
    has = np.asarray((A > 0).sum(axis=1)).ravel() > 0
    predict = lambda vis: S.graph_vote(A, vis, _unlabelled_but_known(vis, truth))[0]
    return S.label_transfer_test(ctx, "structural_homology", target, predict, restrict=has,
                                 truth=truth, n_null=20)


def _ec_default(ctx):
    return "ec_number" if "ec_number" in ctx.nodes else default_category(ctx)


register(Strategy(
    key="structural_homology", number=14, family=NETWORKS,
    title="Annotate function through shared fold",
    question="What does a protein's fold -- its structural similarity to annotated proteins -- say "
             "about its enzymatic class or domain family, even without sequence homology?",
    tooltip="Transfers a functional annotation (EC class by default) along the structural-"
            "similarity layer from predicted models, calling unannotated and hypothetical proteins "
            "by their structural neighbours; fold is conserved long after sequence is not.",
    explanation=(
        "Sequence-based annotation fails exactly where this parasite is most interesting: in "
        "lineage-specific proteins whose sequence resembles nothing annotated. Fold outlives "
        "sequence, and the structural-similarity layer (Foldseek TM-score between predicted "
        "models) links proteins that share a fold whatever their sequence. This strategy carries a "
        "functional annotation -- the EC number truncated to a chosen depth, or any label -- "
        "along that layer, weighted by TM-score, to proteins that lack it.\n\n"
        "Depth matters: EC level 1 (oxidoreductase, transferase, hydrolase...) is what a shared "
        "fold reliably predicts; level 3 or 4 asks the fold to know the substrate, which it often "
        "does not. The self-test scores exactly the chosen depth on hidden annotated proteins "
        "that have structural neighbours, so the number tells you how deep to trust."),
    walkthrough=(
        "Leave the annotation column at ec_number and level 1 for a first pass.",
        "Press Test; then raise the level to 2 and test again to see how fast trust decays.",
        "Press Run; 'calls' puts hypothetical proteins first.",
        "Check a call against the structure: the neighbours' models are in the evidence panel."),
    test_description=(
        "Pattern 1 restricted to proteins with a structural neighbour: 25% of the annotation "
        "(at the chosen level) hidden by whole orthogroups; hidden proteins called by their "
        "visible structural neighbours. Null: 20 runs on shuffled annotations. Pass: above the "
        "null's 95th percentile by 0.05."),
    params=(Param("target", "column", "Annotation to transfer",
                  "The functional annotation carried along the structural layer: an EC number by "
                  "default, or any column of labels. It is truncated to the chosen level before "
                  "anything is compared.", _ec_default),
            Param("level", "int", "Annotation depth",
                  "How many dot-separated parts of the annotation to keep: 1 is the enzyme class, "
                  "4 the exact reaction. Zero keeps the whole label. Deeper is more specific and "
                  "less often predictable from fold.", 1, lo=0, hi=4),
            Param("layer", "layer", "Structural layer",
                  "The edge layer of structural similarity. Any layer can be chosen, but the "
                  "strategy's reasoning -- fold outlives sequence -- applies to the structural one.",
                  _layer_default(["struct", "domain"]))),
    runner=_struct_run, tester=_struct_test, needs=("a structural-similarity layer",)))


def _multiplex(ctx, p, target):
    from . import methods
    layers = [l for l in ctx.measurement_layers(target) if l != "orthogroup"]
    if len(layers) < 2:
        raise ValueError("multiplex communities need at least two permitted measurement layers")
    ctx.say(f"communities over {', '.join(layers)}")
    res = methods.multiplex_communities(layers, ctx.n, resolution=float(p["resolution"]),
                                        seed=ctx.seed, graph=ctx.graph,
                                        min_agreement=float(p["agreement"]), join="louvain")
    part = np.asarray(res["partition"], dtype=int)
    sizes = pd.Series(part[part >= 0]).value_counts()
    small = sizes[sizes < 5].index
    part = np.where(np.isin(part, small), NOISE, part)
    return layers, part


def _multiplex_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    layers, part = _multiplex(ctx, p, target)
    pos = np.flatnonzero(part != NOISE)
    pred, table = S.predict_from_clusters(part[pos], pos, truth, ctx.n, 0.0, 3)
    calls = _calls_table(ctx, pred, pd.Series(np.nan, index=range(ctx.n)), truth, module=part)
    rows = []
    for m in sorted(set(part.tolist()) - {NOISE}):
        idx = np.flatnonzero(part == m)
        lab = truth.iloc[idx].dropna()
        top = lab.value_counts()
        rows.append({"module": m, "size": len(idx), "labelled": len(lab),
                     "top_label": top.index[0] if len(top) else "",
                     "top_share": float(top.iloc[0] / len(lab)) if len(top) else float("nan"),
                     "genes": ", ".join(ctx.gene_ids[idx][:10])})
    modules = pd.DataFrame(rows).sort_values("size", ascending=False, kind="stable") \
        if rows else pd.DataFrame()
    summary = (f"{len(modules)} communities that the layers {', '.join(layers)} agree on; "
               f"{len(pos):,} genes placed in one. {len(calls):,} unlabelled members are called by "
               f"their community's majority {target}.")
    return StrategyResult("multiplex_modules", summary, {"modules": modules, "calls": calls},
                          labels=part, genes=list(calls["gene_id"][:200]))


def _multiplex_test(ctx, p):
    target = _need(ctx, p["target"])
    _layers, part = _multiplex(ctx, p, target)
    t0 = time.monotonic()
    pos = np.flatnonzero(part != NOISE)
    truth = ctx.truth(target)
    visible, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    hidden = hidden[np.isin(hidden, pos)]
    score = lambda vis: _hidden_f1(part[pos], pos, _label_clusters(
        part[pos], vis.iloc[pos].reset_index(drop=True)), truth, hidden)
    observed = score(visible)
    rng = ctx.rng(23)
    nulls = [score(S.shuffled(visible, rng)) for _ in range(20)]
    return S.judge("multiplex_modules",
                   "F1 of hidden genes in the community chosen for their label on known genes",
                   observed, nulls, min_effect=0.05, n_hidden=len(hidden),
                   hidden=f"25% of the {target} labels, whole orthogroups at a time",
                   null_kind="20 label-to-community choices made on shuffled labels", t0=t0,
                   numbers={"communities": int(len(set(part.tolist()) - {NOISE}))})


register(Strategy(
    key="multiplex_modules", number=15, family=NETWORKS,
    title="Find the communities several networks agree on",
    question="Which groups of genes are communities in more than one kind of measured relationship "
             "at once -- co-expressed AND co-fit AND crosslinked?",
    tooltip="Finds communities in each permitted edge layer and keeps the groupings the layers "
            "agree on, as a consensus rather than a merged graph; a module several independent "
            "relationships agree on is far more likely to be a real biological unit.",
    explanation=(
        "Merging networks into one graph lets the densest layer decide everything. This strategy "
        "does the opposite: communities are found in each layer separately (modularity "
        "optimisation), and two genes end up in the same multiplex module only when a sufficient "
        "share of the layers that cover both put them together. The result is a set of modules "
        "supported by several independent kinds of evidence.\n\n"
        "Modules are label-free, so a label they explain afterwards is a genuine finding; they are "
        "also a list of candidate complexes and pathways to read directly. The self-test chooses "
        "the community that best isolates each label on the known genes, and asks how well it "
        "holds that label's hidden genes -- which checks that the agreement carries biological "
        "signal rather than a shared artefact of the layers, such as a common bias toward "
        "abundant proteins."),
    walkthrough=(
        "Choose a label to read the modules against (it is also withheld from the layers).",
        "Keep resolution 1 and agreement 0.5 for a first pass.",
        "Press Test, then Run.",
        "Read 'modules': size, dominant label and share. A large module with no dominant label is "
        "worth opening gene by gene."),
    test_description=(
        "On genes placed in a community: 25% of the label hidden; the community that best isolates "
        "each label is chosen on the visible genes (the communities themselves never see labels). "
        "Metric: the size-weighted F1 of each label's hidden genes against its chosen community. "
        "Null: 20 choices made on shuffled labels. Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET,
            Param("resolution", "float", "Resolution",
                  "Modularity resolution for each layer's communities. Above 1 gives more, smaller "
                  "communities; below 1 fewer, larger ones.", 1.0, lo=0.1, hi=5.0, step=0.1),
            Param("agreement", "float", "Layer agreement",
                  "The share of the layers covering two genes that must put them in one community "
                  "before the consensus does. Higher is stricter and leaves more genes unplaced.",
                  0.5, lo=0.1, hi=1.0, step=0.05)),
    runner=_multiplex_run, tester=_multiplex_test, cost="a minute",
    needs=("two or more edge layers",)))


#: The most candidate pairs a link prediction scores. A layer as dense as co-fitness has millions
#: of two-step pairs; the strongest by raw support are kept, since the rest could not rank anyway.
LINK_CANDIDATES = 300000
#: Fewer edges than this and there is nothing to learn a contact from, nor a fifth worth hiding.
LINK_MIN_EDGES = 50


def _link_setup(ctx, layer):
    """(evidence layers, measurement matrix, the layer's edges as pairs, its nodes)."""
    import scipy.sparse as sp
    evidence = [l for l in ctx.layers() if l != layer and l not in S.DERIVED_LAYERS
                and l not in S.LITERATURE_LAYERS]
    T = ctx.adjacency(layer, "binary")
    U = sp.triu(T, 1).tocoo()
    edges = np.column_stack([U.row, U.col]).astype(int)
    nodes = np.flatnonzero(np.asarray(T.sum(axis=1)).ravel() > 0)
    X, _cols = ctx.features(None)
    return evidence, X, edges, nodes


def _pair_features(ctx, pairs, evidence, T, X, perm=None) -> np.ndarray:
    """Per pair: direct and two-step support in each evidence layer, closure in the layer itself,
    degree, and measurement similarity. Read at `perm` of the identities, for the null."""
    a, b = pairs[:, 0], pairs[:, 1]
    if perm is not None:
        a, b = perm[a], perm[b]
    cols = []
    for l in evidence:
        A = ctx.adjacency(l, "binary")
        cols.append(np.asarray(A[a, b]).ravel())
        cols.append(np.log1p(np.asarray(A[a].multiply(A[b]).sum(axis=1)).ravel()))
    cols.append(np.log1p(np.asarray(T[a].multiply(T[b]).sum(axis=1)).ravel()))
    deg = np.asarray(T.sum(axis=1)).ravel()
    cols.append(np.log1p(deg[a]) + np.log1p(deg[b]))
    xa, xb = X[a], X[b]
    cols.append((xa * xb).sum(axis=1) / (np.linalg.norm(xa, axis=1) * np.linalg.norm(xb, axis=1)
                                         + 1e-9))
    return np.column_stack(cols)


def _non_edges(rng, nodes, edges_set, n_pairs: int, n: int) -> np.ndarray:
    out = set()
    tries = 0
    while len(out) < n_pairs and tries < 50 * n_pairs:
        tries += 1
        a, b = rng.choice(nodes, size=2, replace=False)
        a, b = (a, b) if a < b else (b, a)
        if a * n + b not in edges_set:
            out.add((int(a), int(b)))
    return np.array(sorted(out), dtype=int).reshape(-1, 2)


def _link_model():
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight="balanced",
                                                              max_iter=2000))


def _visible_layer(n, pairs):
    import scipy.sparse as sp
    T = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
    return T.maximum(T.T).tocsr()


def _link_run(ctx, p):
    layer = p["layer"]
    if layer not in ctx.layers():
        raise ValueError(f"this table has no edge layer {layer!r}")
    evidence, X, edges, nodes = _link_setup(ctx, layer)
    if len(edges) < LINK_MIN_EDGES:
        return StrategyResult("link_prediction", f"The {layer} layer has {len(edges)} edges -- too "
                              f"few to learn what makes one.", {})
    rng = ctx.rng(31)
    n = ctx.n
    eset = set((edges[:, 0] * n + edges[:, 1]).tolist())
    T = ctx.adjacency(layer, "binary")
    neg = _non_edges(rng, nodes, eset, 5 * len(edges), n)
    if not len(neg):
        return StrategyResult("link_prediction", f"Every pair of genes in the {layer} layer is "
                              f"already an edge: there is no contact left to predict.", {})
    model = _link_model().fit(_pair_features(ctx, np.vstack([edges, neg]), evidence, T, X),
                              np.r_[np.ones(len(edges)), np.zeros(len(neg))])
    ctx.say("scoring candidate pairs")
    import scipy.sparse as sp
    C = (T @ T)
    for l in evidence:
        C = C + ctx.adjacency(l, "binary")
    sub = sp.triu(C[nodes][:, nodes], 1).tocoo()
    cand = np.column_stack([nodes[sub.row], nodes[sub.col]])
    new_pair = np.array([a * n + b not in eset for a, b in cand], dtype=bool)
    cand, strength = cand[new_pair], sub.data[new_pair]
    if len(cand) > LINK_CANDIDATES:
        cand = cand[np.argsort(-strength, kind="stable")[:LINK_CANDIDATES]]
    if not len(cand):
        return StrategyResult("link_prediction", "No candidate pairs.", {})
    score = model.decision_function(_pair_features(ctx, cand, evidence, T, X))
    order = np.argsort(-score)[: int(p["top"])]
    top = cand[order]
    support = [", ".join(l for l in evidence if ctx.adjacency(l, "binary")[a, b] > 0)
               for a, b in top]
    table = pd.DataFrame({"gene_a": ctx.gene_ids[top[:, 0]], "product_a": ctx.product(top[:, 0]),
                          "gene_b": ctx.gene_ids[top[:, 1]], "product_b": ctx.product(top[:, 1]),
                          "score": score[order],
                          "shared_partners": np.asarray(T[top[:, 0]].multiply(T[top[:, 1]]).sum(
                              axis=1)).ravel(),
                          "linked_in": support})
    summary = (f"A model of what makes a {layer} edge -- trained on its {len(edges):,} edges against "
               f"random non-edges, from {len(evidence)} other layers, shared partners and "
               f"measurement similarity -- scored {len(cand):,} candidate pairs. The top "
               f"{len(table)} are listed; each is a contact the experiment did not see and the rest "
               f"of the evidence predicts.")
    return StrategyResult("link_prediction", summary, {"predicted links": table},
                          genes=list(dict.fromkeys(list(table["gene_a"]) + list(table["gene_b"])))[:200])


def _link_test(ctx, p):
    t0 = time.monotonic()
    layer = p["layer"]
    if layer not in ctx.layers():
        raise ValueError(f"this table has no edge layer {layer!r}")
    evidence, X, edges, nodes = _link_setup(ctx, layer)
    if len(edges) < LINK_MIN_EDGES:
        return S.judge("link_prediction", "AUROC", float("nan"), [], min_effect=0.05,
                       n_hidden=0, hidden=f"{layer} edges", null_kind="permuted identities",
                       t0=t0, note=f"the {layer} layer has {len(edges)} edges, too few to learn from")
    n = ctx.n
    rng = ctx.rng(37)
    order = rng.permutation(len(edges))
    k = int(round(0.2 * len(edges)))
    hidden, visible = edges[order[:k]], edges[order[k:]]
    eset = set((edges[:, 0] * n + edges[:, 1]).tolist())
    T = _visible_layer(n, visible)
    train_neg = _non_edges(rng, nodes, eset, 5 * len(visible), n)
    test_neg = _non_edges(rng, nodes, eset | set((train_neg[:, 0] * n + train_neg[:, 1]).tolist()),
                          5 * len(hidden), n)
    y = np.r_[np.ones(len(visible)), np.zeros(len(train_neg))]
    train = np.vstack([visible, train_neg])

    def fit_score(pairs, perm):
        m = _link_model().fit(_pair_features(ctx, train, evidence, T, X, perm), y)
        return m.decision_function(_pair_features(ctx, pairs, evidence, T, X, perm))

    perms = [rng.permutation(n) for _ in range(5)]
    return S.pair_test("link_prediction", lambda pairs: fit_score(pairs, None), hidden, test_neg,
                       lambda pairs, i: fit_score(pairs, perms[i]), n_null=5, min_effect=0.05,
                       hidden=f"20% of the {layer} edges ({len(hidden):,})",
                       null_kind="5 models trained with every gene's evidence read from a random "
                                 "other gene", t0=t0, ctx=ctx)


register(Strategy(
    key="link_prediction", number=16, family=NETWORKS,
    title="Predict the contacts an interactome missed",
    question="Which pairs of proteins are probably in physical contact although the crosslinking "
             "or pulldown experiment never saw them together?",
    tooltip="Learns what distinguishes a measured contact from a random pair -- shared partners, "
            "support in the other layers, measurement similarity -- and ranks unmeasured pairs "
            "by it; every interactome misses contacts, and this says which ones.",
    explanation=(
        "An interactome is a sample. Crosslinking sees lysine-rich, abundant, soluble contacts and "
        "misses the rest; a pulldown sees what survives the wash. So among the pairs never "
        "observed are many true contacts, and the other evidence can say which. This strategy "
        "trains a logistic model on the layer's own edges against random pairs of the same "
        "proteins, using: support in each other measured layer (a direct edge, and shared "
        "neighbours), triadic closure within the layer itself (two proteins crosslinked to the "
        "same partners), each protein's degree, and the cosine similarity of their measurements. "
        "Every unobserved pair with any support is then scored.\n\n"
        "Literature layers and layers derived from the target (unwritten interactions are "
        "crosslinks minus co-mentions) are excluded as evidence. The self-test hides a fifth of the "
        "layer's edges, retrains without them, and asks whether they outrank random non-edges."),
    walkthrough=(
        "Choose the interaction layer to complete (crosslinks by default).",
        "Press Test: AUROC of the hidden edges against random pairs, against models whose evidence "
        "was scrambled across genes.",
        "Press Run; 'predicted links' lists the strongest unobserved pairs with the layers that "
        "support each.",
        "A predicted link between a labelled and an unlabelled protein is also a localization "
        "hypothesis -- take it to strategy 13."),
    test_description=(
        "Pattern 3: 20% of the layer's edges hidden; the model is trained on the rest against "
        "random non-edges and scores the hidden edges against fresh random non-edges. Metric: "
        "AUROC. Null: 5 models trained with each gene's evidence read from a random other gene "
        "(identities permuted). Pass: above the null's 95th percentile by 0.05."),
    params=(Param("layer", "layer", "Interaction layer to complete",
                  "The layer whose missing edges are predicted. Its own visible edges are used "
                  "(shared partners, degree), the other measured layers are evidence, and layers "
                  "derived from it are excluded.", _layer_default(["xlms", "ip_ms", "struct"])),
            Param("top", "int", "Pairs to list",
                  "How many of the highest-scoring unobserved pairs to list. Every candidate is "
                  "scored; this only limits what is shown and saved.", 300, lo=10, hi=20000,
                  step=50)),
    runner=_link_run, tester=_link_test, cost="a minute", needs=("an interaction layer",)))


def _literature_layer_default(ctx):
    have = ctx.layers()
    for l in ("comention_ft", "comention"):
        if l in have:
            return l
    return None


def _lit_pairs(ctx, layer, target):
    g = ctx.graph
    if not layer or f"{layer}__a" not in g:
        raise ValueError("this table has no literature co-mention layer")
    a, b = np.asarray(g[f"{layer}__a"], int), np.asarray(g[f"{layer}__b"], int)
    w = np.asarray(g[f"{layer}__w"], float)
    r = np.asarray(g.get(f"{layer}__r", w), float)
    truth = ctx.truth(target)
    la, lb = truth.iloc[a].to_numpy(dtype=object), truth.iloc[b].to_numpy(dtype=object)
    both = np.array([isinstance(x, str) and isinstance(y, str) for x, y in zip(la, lb)])
    same = np.array([x == y for x, y in zip(la, lb)]) & both
    return a, b, w, r, both, same, f"{layer}__r" in g


def _attention_run(ctx, p):
    target = _need(ctx, p["target"])
    a, b, w, r, both, same, has_r = _lit_pairs(ctx, p["layer"], target)
    k = int(p["top"])
    rank_r, rank_w = np.argsort(-r)[:k], np.argsort(-w)[:k]
    rate = lambda idx: float(same[idx][both[idx]].mean()) if both[idx].any() else float("nan")
    table = pd.DataFrame({"gene_a": ctx.gene_ids[a[rank_r]], "product_a": ctx.product(a[rank_r]),
                          "gene_b": ctx.gene_ids[b[rank_r]], "product_b": ctx.product(b[rank_r]),
                          "raw_count": w[rank_r], "corrected": r[rank_r],
                          f"same_{target}": np.where(both[rank_r], same[rank_r], np.nan)})
    raw = pd.DataFrame({"gene_a": ctx.gene_ids[a[rank_w]], "gene_b": ctx.gene_ids[b[rank_w]],
                        "raw_count": w[rank_w], "corrected": r[rank_w]})
    summary = (f"Top {k} {p['layer']} pairs by attention-corrected residual: {rate(rank_r):.0%} "
               f"share a {target} label (where both have one); by raw count: {rate(rank_w):.0%}; "
               f"all pairs: {float(same[both].mean()) if both.any() else float('nan'):.0%}. "
               + ("" if has_r else "This layer carries no residual, so both rankings are raw. ")
               + "The corrected ranking finds pairs the literature links more than their fame "
                 "predicts -- the ones worth reading about.")
    return StrategyResult("attention_correction", summary,
                          {"corrected ranking": table, "raw ranking": raw},
                          genes=list(dict.fromkeys(list(table["gene_a"]) + list(table["gene_b"])))[:200])


def _attention_test(ctx, p):
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    a, b, w, r, both, same, has_r = _lit_pairs(ctx, p["layer"], target)
    pool = np.flatnonzero(both)
    k = int(min(int(p["top"]), max(20, len(pool) // 5)))
    if len(pool) < 2 * k:
        return S.judge("attention_correction", "same-label share of the top pairs", float("nan"),
                       [], min_effect=0.05, n_hidden=0, hidden="pairs", null_kind="random pairs",
                       t0=t0, note="too few co-mentioned pairs with both genes labelled")
    top_r = pool[np.argsort(-r[pool], kind="stable")[:k]]
    top_w = pool[np.argsort(-w[pool], kind="stable")[:k]]
    rng = ctx.rng(41)
    nulls = [float(same[rng.choice(pool, size=k, replace=False)].mean()) for _ in range(20)]
    return S.judge("attention_correction",
                   f"share of the top {k} corrected pairs sharing a {target} label",
                   float(same[top_r].mean()), nulls, min_effect=0.05, n_hidden=k,
                   hidden=f"the {target} labels, never used to build the literature layer",
                   null_kind="20 random sets of co-mentioned pairs", t0=t0,
                   note="" if has_r else "the layer has no residual; raw counts were ranked",
                   numbers={"raw_ranking_share": float(same[top_w].mean())})


register(Strategy(
    key="attention_correction", number=17, family=NETWORKS,
    title="Read the literature for biology, not fame",
    question="Which pairs of genes are written about together more than their popularity "
             "explains -- and are those pairs biologically related?",
    tooltip="Ranks co-mentioned gene pairs by the residual over what each gene's publication count "
            "predicts, rather than by raw counts that reward fame; checks against held-out "
            "labels that the corrected ranking finds real relationships.",
    explanation=(
        "Raw co-mention reproduces the literature's popularity contest: the two most-published "
        "genes are co-mentioned most often whether or not they have anything to do with each "
        "other. The project's attention correction replaces each count with its residual over the "
        "count expected from each gene's own publication total, and it is on by default in the "
        "application. This strategy uses the correction as an inference tool: pairs with a high "
        "residual are pairs the literature treats as related beyond what fame explains.\n\n"
        "The test of that claim is independent biology. Genes the literature links for a reason "
        "should share a compartment or complex more often than random co-mentioned pairs do. The "
        "strategy reports that share for the corrected ranking, the raw ranking and random pairs; "
        "on the real full-text layer the corrected ranking returns known complexes where the raw "
        "one returns famous genes."),
    walkthrough=(
        "Choose the literature layer (full-text paragraphs by default) and the label to check "
        "against.",
        "Press Test: the share of the top corrected pairs that share a label, against random "
        "co-mentioned pairs; the raw ranking's share is reported beside it.",
        "Press Run and read 'corrected ranking'. Pairs with a high residual and no shared label "
        "are the literature's hypotheses the labels do not yet explain.",
        "Compare with 'raw ranking' to see what fame alone would have put at the top."),
    test_description=(
        "The label is never used to build the literature layer. Among co-mentioned pairs with "
        "both genes labelled, the top k by corrected residual are taken (k = the chosen number, "
        "capped at a fifth of the pool). Metric: the share of those pairs sharing a label. Null: "
        "20 random sets of k co-mentioned pairs. Pass: above the null's 95th percentile by 0.05."),
    params=(TARGET, Param("layer", "layer", "Literature layer",
                          "Which co-mention layer to read: abstracts, or full-text paragraphs. The "
                          "full-text layer is an open-access subset and is less dominated by "
                          "review-article lists.", _literature_layer_default),
            Param("top", "int", "Top pairs",
                  "How many of the highest-ranked pairs are compared. The test caps this at a "
                  "fifth of the pairs whose genes are both labelled, so the top is really a top.",
                  200, lo=10, hi=5000, step=10)),
    runner=_attention_run, tester=_attention_test, needs=("a co-mention layer",)))


UNWRITTEN_SKIP = ("orthogroup", "domain", "compartment")


def _support(ctx):
    """(summed binary support over measurement layers, the layers, literature union or None)."""
    layers = [l for l in ctx.measurement_layers(None) if l not in UNWRITTEN_SKIP]
    if not layers:
        return None, [], None
    M = sum(ctx.adjacency(l, "binary") for l in layers)
    lits = [l for l in S.LITERATURE_LAYERS if l in ctx.layers()]
    L = sum(ctx.adjacency(l, "binary") for l in lits) if lits else None
    return M.tocsr(), layers, L


def _unwritten_run(ctx, p):
    import scipy.sparse as sp
    M, layers, L = _support(ctx)
    if M is None:
        return StrategyResult("unwritten_links", "This table has no measurement layers.", {})
    U = sp.triu(M, 1).tocoo()
    keep = U.data >= int(p["min_layers"])
    a, b, s = U.row[keep], U.col[keep], U.data[keep]
    if L is not None and len(a):
        written = np.asarray(L[a, b]).ravel() > 0
        a, b, s = a[~written], b[~written], s[~written]
    order = np.argsort(-s, kind="stable")[: int(p["top"])]
    a, b, s = a[order], b[order], s[order]
    which = [", ".join(l for l in layers if ctx.adjacency(l, "binary")[x, y] > 0) for x, y in zip(a, b)]
    table = pd.DataFrame({"gene_a": ctx.gene_ids[a], "product_a": ctx.product(a),
                          "gene_b": ctx.gene_ids[b], "product_b": ctx.product(b),
                          "layers": s.astype(int), "linked_in": which})
    for c in ("n_papers_focal",):
        if c in ctx.nodes:
            table["papers_a"] = ctx.values(c).iloc[a].to_numpy()
            table["papers_b"] = ctx.values(c).iloc[b].to_numpy()
    summary = (f"{len(table):,} pairs linked in at least {int(p['min_layers'])} independent "
               f"measurement layers ({', '.join(layers)}) that no abstract or paragraph mentions "
               f"together. Each is a relationship the data asserts and the field has not written "
               f"down.")
    return StrategyResult("unwritten_links", summary, {"unwritten pairs": table},
                          genes=list(dict.fromkeys(list(table["gene_a"]) + list(table["gene_b"])))[:200])


def _unwritten_test(ctx, p):
    import scipy.sparse as sp
    t0 = time.monotonic()
    M, layers, L = _support(ctx)
    if M is None or L is None:
        return S.judge("unwritten_links", "AUROC", float("nan"), [], min_effect=0.02, n_hidden=0,
                       hidden="literature pairs", null_kind="permuted identities", t0=t0,
                       note="needs both measurement and literature layers")
    universe = np.flatnonzero(np.asarray((M > 0).sum(axis=1)).ravel() > 0)
    inside = np.zeros(ctx.n, bool)
    inside[universe] = True
    Lu = sp.triu(L, 1).tocoo()
    keep = inside[Lu.row] & inside[Lu.col]
    pos = np.column_stack([Lu.row[keep], Lu.col[keep]]).astype(int)
    rng = ctx.rng(43)
    lset = set((pos[:, 0] * ctx.n + pos[:, 1]).tolist())
    neg = _non_edges(rng, universe, lset, 5 * len(pos), ctx.n)
    score = lambda pairs, perm=None: np.asarray(M[(pairs[:, 0] if perm is None else perm[pairs[:, 0]]),
                                                  (pairs[:, 1] if perm is None else perm[pairs[:, 1]])]).ravel()
    perms = [rng.permutation(ctx.n) for _ in range(10)]
    return S.pair_test("unwritten_links", score, pos, neg, lambda pairs, i: score(pairs, perms[i]),
                       n_null=10, min_effect=0.02,
                       hidden=f"the literature's {len(pos):,} co-mentioned pairs, never used to "
                              f"count measurement support",
                       null_kind="10 runs with the measurement layers' gene identities permuted",
                       t0=t0, ctx=ctx)


register(Strategy(
    key="unwritten_links", number=18, family=NETWORKS,
    title="List what the data says and the literature has not written",
    question="Which gene pairs do several independent measurements link that no paper has ever "
             "mentioned together?",
    tooltip="Counts, for every gene pair, how many independent measurement layers link it, and "
            "lists the strongly supported pairs the literature never co-mentions -- after checking "
            "that measurement support predicts what the literature does write about.",
    explanation=(
        "The literature is where knowledge is recorded; the measurements are where it could come "
        "from. A pair linked by co-expression, co-fitness AND a crosslink but never mentioned "
        "together in any abstract or paragraph is a relationship the data asserts and nobody has "
        "written down -- a gap in the knowledge map, and a concrete hypothesis.\n\n"
        "That reading needs one thing to be true: that measurement support tracks the kind of "
        "relationships biologists do eventually write about. The self-test checks it directly: "
        "using the literature only as held-out truth, it asks whether pairs with more measurement "
        "support are more often co-mentioned than random pairs of the same genes, against the "
        "same score with gene identities scrambled. If that holds, the unwritten pairs at the top "
        "of the list are the ones most likely to be written next. Paralogs, shared domains and "
        "the compartment layer are left out of the support count: they link genes for reasons "
        "that are already written down."),
    walkthrough=(
        "Set how many independent layers must agree (2 is the minimum that means anything).",
        "Press Test to confirm that measurement support predicts co-mention on this table.",
        "Press Run and read 'unwritten pairs'; the 'linked_in' column names the layers, and the "
        "paper counts show whether both genes are obscure or one is famous and the link is not.",
        "A pair between an understudied gene and a well-studied one is the cheapest to follow up."),
    test_description=(
        "Pattern 3 with the literature as truth: co-mentioned pairs among genes in the measurement "
        "layers against 5 times as many random pairs. Score: number of measurement layers "
        "linking the pair. Metric: AUROC. Null: 10 runs with the measurement layers' gene "
        "identities permuted. Pass: above the null's 95th percentile by 0.02."),
    params=(Param("min_layers", "int", "Layers that must agree",
                  "How many independent measurement layers must link a pair before it is listed. "
                  "One layer is a single experiment's word; two or more is agreement between "
                  "different kinds of evidence.", 2, lo=1, hi=6),
            Param("top", "int", "Pairs to list",
                  "How many unwritten pairs to list, most supported first. Every pair is counted; "
                  "this only limits what is shown and saved.", 500, lo=10, hi=50000, step=50)),
    runner=_unwritten_run, tester=_unwritten_test, needs=("measurement and literature layers",)))


# =========================================================================== 4 · learn from examples
def _logistic(X, vis: pd.Series, query, C: float = 0.5, min_prob: float = 0.0) -> tuple:
    """(prediction, probability, model): a class-balanced logistic regression on the visible genes."""
    from sklearn.linear_model import LogisticRegression
    vis = pd.Series(vis).reset_index(drop=True)
    query = np.asarray(query, dtype=int)
    pred = pd.Series([np.nan] * len(vis), dtype=object)
    prob = pd.Series(np.nan, index=range(len(vis)))
    known = np.flatnonzero(vis.notna().to_numpy())
    y = vis.iloc[known].astype(str).to_numpy()
    if len(set(y)) < 2 or not len(query) or X.shape[1] == 0:
        return pred, prob, None
    model = LogisticRegression(C=float(C), max_iter=500, class_weight="balanced").fit(X[known], y)
    P = model.predict_proba(X[query])
    best, top = P.argmax(axis=1), P.max(axis=1)
    pred.iloc[query] = np.where(top >= float(min_prob), model.classes_[best], None)
    pred = pred.where(pred.notna(), np.nan)
    prob.iloc[query] = top
    return pred, prob, model


def _classifier_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    X, cols = ctx.features(target)
    pred, prob, model = _logistic(X, t, np.flatnonzero(t.isna().to_numpy()), p["C"],
                                  p["min_probability"])
    calls = _calls_table(ctx, pred, prob, truth)
    rows = []
    if model is not None:
        for i, c in enumerate(model.classes_):
            w = model.coef_[i] if model.coef_.shape[0] > 1 else model.coef_[0] * (1 if i else -1)
            for j in np.argsort(-w)[:8]:
                rows.append({"class": c, "measurement": cols[j], "weight": float(w[j])})
    summary = (f"A class-balanced logistic regression on {len(cols)} permitted measurements, trained "
               f"on {int(t.notna().sum()):,} labelled genes, calls {len(calls):,} unlabelled genes "
               f"(probability >= {float(p['min_probability']):.2f}). 'what drives each class' lists "
               f"the measurements with the largest weights -- the model's reasons, to be read as "
               f"hypotheses about what defines each class.")
    return StrategyResult("supervised_classifier", summary,
                          {"calls": calls, "what drives each class": pd.DataFrame(rows)},
                          genes=list(calls["gene_id"][:200]))


def _classifier_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    predict = lambda vis: _logistic(X, vis, _unlabelled_but_known(vis, truth), p["C"],
                                    p["min_probability"])[0]
    return S.label_transfer_test(ctx, "supervised_classifier", target, predict, null="analytic")


register(Strategy(
    key="supervised_classifier", number=19, family=LEARN,
    title="Train a classifier on the known genes and call the rest",
    question="Given every permitted measurement, which label does a model trained on the labelled "
             "genes assign to each unlabelled one -- and which measurements does it rely on?",
    tooltip="Fits a class-balanced logistic regression from the permitted measurements to the "
            "label, scores it on whole orthogroups it never trained on, and calls unlabelled "
            "genes with a probability; also lists the measurements each class is recognised by.",
    explanation=(
        "Neighbour votes treat every measurement as equally relevant. A classifier learns which "
        "measurements matter for which class: a dense-granule protein might be recognised by its "
        "expression in the tachyzoite and its secretion signal, a ribosomal protein by fitness "
        "and codon usage. This strategy fits a multinomial logistic regression, class-balanced so "
        "the commonest compartment does not win by default, with an L2 penalty whose strength you "
        "choose.\n\n"
        "The weights are an interpretable by-product: for each class, the measurements with the "
        "largest positive weights are what the model thinks defines it. The self-test hides whole "
        "orthogroups -- a paralog left in training makes any classifier look better than it is -- "
        "and compares the correct-call rate with the chance level given the predicted and true "
        "class frequencies (Cohen's expectation), which is what an uninformed classifier with the "
        "same output mix would score."),
    walkthrough=(
        "Choose the label; keep C at 0.5 (smaller is a stronger penalty and a simpler model).",
        "Press Test (one fit; a few seconds to a minute on the full table).",
        "Press Run; read 'what drives each class' before 'calls' -- a class recognised by "
        "sensible measurements earns more trust than one recognised by missingness.",
        "Raise the probability threshold to call fewer, surer genes."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; the model is trained on the "
        "visible genes and calls the hidden ones. Metric: hidden genes called correctly. Null: the "
        "analytic chance level for the same predicted and true class mixes (re-running a "
        "multinomial fit on shuffled labels ten times would take longer than it tells). Pass: "
        "above that level's 95% bound by 0.05."),
    params=(TARGET,
            Param("C", "float", "Regularisation C",
                  "Inverse penalty strength. Small values force a simple model that uses few "
                  "measurements strongly; large values let it fit detail, including noise. 0.5 is "
                  "a middle setting for a few hundred rank-scaled measurements.", 0.5,
                  lo=0.001, hi=100.0, step=0.1),
            Param("min_probability", "float", "Call when probability is at least",
                  "A gene is called only when the model's top probability reaches this. Zero calls "
                  "every gene; 0.6 keeps the calls the model is sure of.", 0.0, lo=0.0, hi=1.0,
                  step=0.05)),
    runner=_classifier_run, tester=_classifier_test, cost="a minute",
    needs=("a categorical column",)))


def _pu_scores(ctx, X, query, bags: int, salt: int = 0) -> np.ndarray:
    """Bagged positive-unlabelled scores: how much each gene resembles the positives, out of bag."""
    from sklearn.linear_model import LogisticRegression
    rng = ctx.rng(500 + salt)
    n = X.shape[0]
    P = np.unique(np.asarray(query, dtype=int))
    U = np.setdiff1d(np.arange(n), P)
    size = int(min(len(U), max(3 * len(P), 50)))
    total, count = np.zeros(n), np.zeros(n)
    last = np.zeros(n)
    for _ in range(int(bags)):
        ctx.check()
        ub = rng.choice(U, size=size, replace=False)
        m = LogisticRegression(C=0.5, max_iter=300, class_weight="balanced").fit(
            X[np.r_[P, ub]], np.r_[np.ones(len(P)), np.zeros(size)])
        s = m.decision_function(X)
        out = np.ones(n, bool)
        out[ub] = out[P] = False
        total[out] += s[out]
        count[out] += 1
        last = s
    scores = np.where(count > 0, total / np.maximum(count, 1), last)
    scores[P] = scores.max() + 1.0
    return scores


def _pu_run(ctx, p):
    pos, missing = ctx.resolve_genes(p["genes"])
    if len(pos) < 5:
        return StrategyResult("positive_unlabeled", "Paste at least five gene ids from this organism.",
                              {})
    X, cols = ctx.features(p.get("exclude"))
    scores = _pu_scores(ctx, X, pos, int(p["bags"]))
    order = [g for g in np.argsort(-scores) if g not in set(pos)][: int(p["top"])]
    table = _genes_table(ctx, order, score=scores[order])
    summary = (f"{len(pos)} positives against the rest of the proteome, {int(p['bags'])} bags; the "
               f"{len(table)} genes that most resemble your list are listed, scored only by models "
               f"that did not train on them. Nothing on the list is assumed about any other gene: "
               f"unlabelled is treated as unknown, not as negative.")
    return StrategyResult("positive_unlabeled", summary,
                          {"candidates": table, "not found": pd.DataFrame({"identifier": missing})},
                          genes=list(table["gene_id"]))


def _pu_test(ctx, p):
    members, what, exclude = _list_or_example(ctx, p)
    X, _cols = ctx.features(exclude)
    calls = iter(range(10 ** 6))
    return S.set_expansion_test(ctx, "positive_unlabeled", members,
                                lambda q: _pu_scores(ctx, X, q, int(p["bags"]), next(calls)),
                                n_null=10, label=what)


register(Strategy(
    key="positive_unlabeled", number=20, family=LEARN,
    title="Learn what makes your list special, from positives alone",
    question="Given only genes that ARE something -- no list of genes that are not -- which other "
             "genes look most like them?",
    tooltip="Trains many classifiers, each separating your list from a random draw of the rest of "
            "the proteome, and ranks every other gene by out-of-bag score; the right tool when "
            "all you have is positives and the unlabelled genes include unknown positives.",
    explanation=(
        "Most real gene lists are positives only: hits of a screen, members of a complex, proteins "
        "someone localized. Training 'list versus everything else' treats every unknown member as a "
        "negative and teaches the model to reject exactly the genes you want to find. "
        "Positive-unlabelled bagging avoids that: each of many models sees the list and a small "
        "random draw of other genes; because the draw is small, a hidden positive is rarely in "
        "it, and every gene is scored only by the models that did not train on it.\n\n"
        "Unlike strategy 02, no cluster has to form: the model can combine weak signals from many "
        "measurements that no map would separate on. Unlike strategy 25, it uses the "
        "measurements rather than the networks. If the list came from a column of this table, "
        "name it so its closure is withheld."),
    walkthrough=(
        "Paste the list; name the column it came from, if any.",
        "Press Test: with 20 or more genes it hides 30% of YOUR list; otherwise it tests on a "
        "known category of similar size.",
        "Press Run and read 'candidates', highest score first.",
        "Take the top candidates to strategy 24 to see what they share with the list."),
    test_description=(
        "Pattern 2: 30% of the set hidden; the other 70% are the positives. Metric: AUROC of the "
        "hidden members against every other non-query gene. Null: 10 random sets of the same size. "
        "Pass: above the null's 95th percentile by 0.1."),
    params=(GENES, EXCLUDE,
            Param("bags", "int", "Bags",
                  "How many positive-versus-random-draw models are averaged. More bags make the "
                  "ranking steadier; 15 is usually enough, and each bag is a small fit.", 15,
                  lo=3, hi=200),
            Param("top", "int", "Candidates to list",
                  "How many of the highest-scoring genes not on your list to show. Every gene is "
                  "scored; this only limits the table.", 300, lo=10, hi=10000, step=50)),
    runner=_pu_run, tester=_pu_test, needs=("a gene list",)))


def _regressor(kind: str):
    if kind == "ridge":
        from sklearn.linear_model import Ridge
        return Ridge(alpha=1.0)
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, random_state=0)


def _fit_predict(X, y: pd.Series, kind: str) -> np.ndarray:
    """Fit on the measured genes of `y`, predict every gene."""
    y = pd.Series(y).reset_index(drop=True)
    known = np.flatnonzero(y.notna().to_numpy())
    if len(known) < 20:
        return np.full(X.shape[0], np.nan)
    model = _regressor(kind).fit(X[known], y.iloc[known].to_numpy(dtype=float))
    return model.predict(X)


def _regression_matrix(ctx, target, kind, own_kind: str = "leave out"):
    banned = set(ctx.banned(target))
    if own_kind == "leave out":
        for t in ([target] if isinstance(target, str) else list(target)):
            banned |= ctx.same_kind(t)
    cols = ctx.numeric_columns(banned)
    return ctx.matrix(cols, impute=(kind == "ridge")), cols


def _out_of_fold(ctx, X, y: pd.Series, kind: str, folds: int = 5) -> np.ndarray:
    y = pd.Series(y).reset_index(drop=True)
    known = np.flatnonzero(y.notna().to_numpy())
    rng = ctx.rng(61)
    g = ctx.groups()[known]
    fold_of = {grp: i % folds for i, grp in enumerate(rng.permutation(np.unique(g)))}
    fold = np.array([fold_of[x] for x in g])
    out = np.full(len(y), np.nan)
    for f in range(folds):
        ctx.check()
        train = y.copy()
        train.iloc[known[fold == f]] = np.nan
        out[known[fold == f]] = _fit_predict(X, train, kind)[known[fold == f]]
    return out


def _regression_run(ctx, p):
    target = _need(ctx, p["target"])
    y = ctx.values(target)
    X, cols = _regression_matrix(ctx, target, p["model"], p["own_kind"])
    oof = _out_of_fold(ctx, X, y, p["model"])
    rho = S.spearman(oof, y)
    full = _fit_predict(X, y, p["model"])
    measured = y.notna().to_numpy()
    resid = (y.to_numpy() - oof)
    z = (resid - np.nanmean(resid)) / (np.nanstd(resid) or 1.0)
    surprising = _genes_table(ctx, np.flatnonzero(measured), measured=y[measured].to_numpy(),
                              expected=oof[measured], z=z[measured])
    surprising = surprising.reindex(surprising["z"].abs().sort_values(ascending=False).index).head(
        int(p["top"])).reset_index(drop=True)
    unmeasured = np.flatnonzero(~measured)
    predicted = _genes_table(ctx, unmeasured, predicted=full[unmeasured]).sort_values(
        "predicted", kind="stable").reset_index(drop=True)
    summary = (f"{target} predicted from {len(cols)} permitted measurements ({p['model']}); "
               f"out-of-fold rank correlation {rho:.3f} over {int(measured.sum()):,} measured "
               f"genes. 'surprising genes' are the measured ones furthest from what everything "
               f"else predicts -- where this gene does something its profile does not explain. "
               f"'predicted' gives a value for the {len(unmeasured):,} genes never measured.")
    return StrategyResult("trait_regression", summary,
                          {"surprising genes": surprising, "predicted": predicted},
                          genes=list(surprising["gene_id"]), numbers={"oof_spearman": rho})


def _regression_test(ctx, p):
    target = _need(ctx, p["target"])
    X, _cols = _regression_matrix(ctx, target, p["model"], p["own_kind"])
    return S.value_test(ctx, "trait_regression", ctx.values(target),
                        lambda vis: _fit_predict(X, vis, p["model"]), n_null=3, min_effect=0.1,
                        label=f"the measured {target} values")


NUMBER = Param("target", "number", "Measurement to predict",
               "The numeric measurement to predict from everything else. Its own closure -- the "
               "same experiment's other outputs and anything restating it -- is withheld first.",
               default_numeric)
OWN_KIND = Param("own_kind", "choice", "Its own kind of measurement",
                 "'leave out' removes every column measured the same way as the target -- all the "
                 "fitness screens, when the target is a fitness screen -- so the answer says what "
                 "OTHER evidence knows. 'include' keeps them, which mostly measures how well the "
                 "screens agree with each other.", "leave out", choices=("leave out", "include"))
MODEL = Param("model", "choice", "Model",
              "'boosted' is gradient-boosted trees: nonlinear, tolerant of missing values, slower. "
              "'ridge' is a penalised straight-line fit: fast, and its weights can be read.",
              "boosted", choices=("boosted", "ridge"))


register(Strategy(
    key="trait_regression", number=21, family=LEARN,
    title="Predict a measurement, and find the genes that defy the prediction",
    question="How well does everything else predict this measurement -- and which genes are far "
             "from what their profile says they should be?",
    tooltip="Predicts a numeric measurement such as a fitness score from every other permitted "
            "measurement, out of fold; reports how predictable it is, fills it in for genes never "
            "measured, and ranks the measured genes that deviate most from expectation.",
    explanation=(
        "Two kinds of knowledge come out of one regression. The first is predictability: if "
        "fitness in fibroblasts can be predicted from expression, modification and structure, "
        "then essentiality has a signature, and genes the screen missed can be given a value. "
        "The second is surprise: a gene whose measured fitness is far from what its profile "
        "predicts is doing something its profile does not explain -- an essential protein with the "
        "profile of a dispensable one is a candidate for an unusual function, or a screen "
        "artefact worth rechecking.\n\n"
        "By default the target's own kind of measurement is left out -- every knockout screen, "
        "when the target is a knockout screen -- because predicting one screen from another "
        "mostly shows that the screens agree. Predictions for measured genes are made out of "
        "fold, whole orthogroups at a time, so a gene is never predicted by a model that saw it or "
        "its paralog. The test hides a fifth of the measured values and scores the rank "
        "correlation on them against the chance distribution of a rank correlation."),
    walkthrough=(
        "Choose the measurement (the fibroblast fitness screen by default).",
        "Press Test: the rank correlation on hidden values, against shuffled training values.",
        "Press Run; the summary gives the out-of-fold correlation over every measured gene.",
        "Read 'surprising genes' by |z|: positive z is higher than predicted, negative lower.",
        "Read 'predicted' for genes the experiment never measured -- with the test's correlation "
        "as the measure of how far to trust it."),
    test_description=(
        "Pattern 4: 20% of the measured values hidden; the model is trained on the rest. Metric: "
        "rank correlation between predicted and hidden values. Null: 3 models trained on the "
        "visible values shuffled among the measured genes. Pass: above the null's 95th percentile "
        "by 0.1."),
    params=(NUMBER, MODEL, OWN_KIND,
            Param("top", "int", "Surprising genes to list",
                  "How many of the measured genes furthest from their prediction to list, by "
                  "absolute z-score of the out-of-fold residual.", 200, lo=10, hi=5000, step=10)),
    runner=_regression_run, tester=_regression_test, cost="a minute",
    needs=("a numeric column",)))


def _soft_impute(M: np.ndarray, rank: int, iterations: int = 30) -> np.ndarray:
    """Iterative low-rank completion: fill, factorise, refill the missing entries, repeat."""
    from sklearn.utils.extmath import randomized_svd
    miss = ~np.isfinite(M)
    Z = np.where(miss, 0.0, M)
    r = int(max(1, min(rank, min(M.shape) - 1)))
    R = Z
    for _ in range(int(iterations)):
        U, s, Vt = randomized_svd(Z, r, random_state=0)
        R = (U * s) @ Vt
        Z = np.where(miss, R, M)
    return R


def _column_reliability(ctx, M, cols, rank, frac=0.1, salt=0) -> pd.DataFrame:
    rng = ctx.rng(70 + salt)
    obs = np.argwhere(np.isfinite(M))
    pick = obs[rng.random(len(obs)) < frac]
    Mh = M.copy()
    Mh[pick[:, 0], pick[:, 1]] = np.nan
    R = _soft_impute(Mh, rank)
    rows = []
    for j, c in enumerate(cols):
        sel = pick[pick[:, 1] == j][:, 0]
        if len(sel) >= 20:
            rows.append({"column": c, "hidden": len(sel), "reliability": S.spearman(R[sel, j],
                                                                                    M[sel, j])})
    return pd.DataFrame(rows)


def _impute_run(ctx, p):
    cols = ctx.numeric_columns(ctx.banned(p.get("target")))
    M = ctx.matrix(cols, impute=False)
    ctx.say("measuring how well each column can be reconstructed")
    rel = _column_reliability(ctx, M, cols, int(p["rank"]))
    rel = rel.sort_values("reliability", ascending=False, kind="stable").reset_index(drop=True)
    rel.insert(1, "missing", [int(ctx.values(c).isna().sum()) for c in rel["column"]])
    ctx.say("completing the matrix")
    R = _soft_impute(M, int(p["rank"]))
    tables = {"column reliability": rel}
    column = p.get("column")
    if column and column in cols:
        j = cols.index(column)
        miss = np.flatnonzero(~np.isfinite(M[:, j]))
        r = rel.set_index("column")["reliability"].get(column, np.nan)
        tables["filled"] = _genes_table(ctx, miss, imputed_percentile=R[miss, j] + 0.5).sort_values(
            "imputed_percentile", ascending=False, kind="stable").reset_index(drop=True)
        tables["filled"]["column_reliability"] = r
    good = int((rel["reliability"] >= 0.5).sum()) if len(rel) else 0
    summary = (f"{len(cols)} measurements completed at rank {int(p['rank'])}. {good} of {len(rel)} "
               f"columns reconstruct hidden values at rank correlation >= 0.5: for those, a missing "
               f"value can be estimated from the rest of the table; for the others it cannot, and "
               f"the absence should stay an absence.")
    return StrategyResult("masked_imputation", summary, tables)


def _impute_test(ctx, p):
    t0 = time.monotonic()
    cols = ctx.numeric_columns(ctx.banned(p.get("target")))
    M = ctx.matrix(cols, impute=False)
    rel = _column_reliability(ctx, M, cols, int(p["rank"]))
    observed = float(np.nanmedian(rel["reliability"])) if len(rel) else float("nan")
    rng = ctx.rng(71)
    nulls = []
    for i in range(3):
        ctx.check()
        Mp = np.column_stack([M[rng.permutation(M.shape[0]), j] for j in range(M.shape[1])])
        r = _column_reliability(ctx, Mp, cols, int(p["rank"]), salt=i + 1)
        nulls.append(float(np.nanmedian(r["reliability"])) if len(r) else float("nan"))
    return S.judge("masked_imputation", "median per-column rank correlation on hidden entries",
                   observed, nulls, min_effect=0.1, n_hidden=int(rel["hidden"].sum()) if len(rel)
                   else 0, hidden="10% of every column's measured entries",
                   null_kind="3 completions of the table with each column shuffled independently",
                   t0=t0, details=rel)


register(Strategy(
    key="masked_imputation", number=22, family=LEARN,
    title="Fill in what was never measured, and say where that is honest",
    question="For each measurement, can its missing values be estimated from the rest of the table "
             "-- and for which measurements is that impossible?",
    tooltip="Completes the whole measurement table with a low-rank model, after first hiding a "
            "tenth of every column to measure how well each one can be reconstructed; imputes "
            "only where the reconstruction has been shown to work.",
    explanation=(
        "Most genes are missing most measurements, and every map and model here fills the gaps "
        "somehow -- usually at the median, which says 'average' about a gene nobody measured. This "
        "strategy asks whether better is possible: if the measurements are correlated (expression "
        "across stages, fitness across hosts), a low-rank model of the table can estimate a "
        "missing entry from the gene's other entries.\n\n"
        "The answer is column-specific, and the strategy measures it before using it: a tenth of "
        "each column's measured values is hidden, the table is completed, and each column gets a "
        "reliability -- the rank correlation between reconstructed and hidden values. Columns that "
        "reconstruct well can be filled; for the rest, an absence must stay an absence, which is "
        "the project's rule that unknown is grey, never zero."),
    walkthrough=(
        "Optionally choose a label whose closure should be left out of the table.",
        "Choose the column you want filled (optional) and a rank (20 is a good default).",
        "Press Test: the median reliability across columns, against a table whose columns were "
        "shuffled so nothing links them.",
        "Press Run and read 'column reliability'; then 'filled' for your chosen column, only if "
        "its reliability is high."),
    test_description=(
        "Pattern 4 across the whole table: 10% of every column's measured entries hidden, the "
        "table completed at the chosen rank. Metric: median over columns of the rank correlation "
        "on hidden entries. Null: 3 completions of a table whose columns were each shuffled "
        "independently. Pass: above the null's 95th percentile by 0.1."),
    params=(Param("target", "category", "Leave out the closure of",
                  "Optionally, a label whose closure is removed from the table before completion, "
                  "so imputed values can be used against that label without leakage.", None,
                  optional=True),
            Param("column", "number", "Column to fill",
                  "The measurement whose missing values are listed after completion. Every column "
                  "is completed; this only chooses which one to show.", default_numeric,
                  optional=True),
            Param("rank", "int", "Rank",
                  "How many underlying patterns the low-rank model may use. Too few blur distinct "
                  "programmes together; too many fit noise and reconstruct nothing.", 20,
                  lo=2, hi=200)),
    runner=_impute_run, tester=_impute_test, cost="a minute"))


def _condition_default(which):
    def pick(ctx):
        nums = set(ctx.numeric_columns())
        for cond, base in (("fit_invivo_PE", "fit_invitro_hff"), ("fit_ifng", "fit_naive_bmdm"),
                           ("expr_gametocyte_v", "expr_asexual_blood"),
                           ("piggybac_mfs", "piggybac_mis")):
            if cond in nums and base in nums:
                return cond if which == "condition" else base
        targets = ctx.numeric_targets()
        return targets[0 if which == "condition" else min(1, len(targets) - 1)] if targets else None
    return pick


def _shift(ctx, condition, baseline) -> pd.Series:
    """Condition minus what the baseline predicts, on rank scale: the condition-specific part."""
    c, b = ctx.values(condition).rank(pct=True), ctx.values(baseline).rank(pct=True)
    both = c.notna() & b.notna()
    if both.sum() < 20:
        raise ValueError(f"{condition} and {baseline} share fewer than 20 measured genes")
    slope, intercept = np.polyfit(b[both], c[both], 1)
    return (c - (slope * b + intercept)).where(both)


def _shift_run(ctx, p):
    cond, base = _need(ctx, p["condition"]), _need(ctx, p["baseline"])
    y = _shift(ctx, cond, base)
    X, cols = _regression_matrix(ctx, [cond, base], p["model"], p["own_kind"])
    oof = _out_of_fold(ctx, X, y, p["model"])
    rho = S.spearman(oof, y)
    measured = y.notna().to_numpy()
    shifted = _genes_table(ctx, np.flatnonzero(measured), shift=y[measured].to_numpy(),
                           explained=oof[measured]).sort_values("shift", kind="stable")
    full = _fit_predict(X, y, p["model"])
    no_cond = np.flatnonzero(ctx.values(cond).isna().to_numpy())
    predicted = _genes_table(ctx, no_cond, predicted_shift=full[no_cond]).sort_values(
        "predicted_shift", kind="stable").reset_index(drop=True)
    summary = (f"The part of {cond} that {base} does not explain -- the condition-specific effect -- "
               f"for {int(measured.sum()):,} genes. Other measurements predict it at out-of-fold "
               f"rank correlation {rho:.3f}: " + ("the condition-specific need has a signature, and "
                                                 "'predicted' extends it to genes never screened in "
                                                 "the condition." if rho >= 0.1 else
                                                 "little of it is predictable from the rest of the "
                                                 "table.") + " Negative shift: needed more in the "
               f"condition than the baseline predicts.")
    return StrategyResult("condition_shift", summary,
                          {"shifted genes": shifted.reset_index(drop=True), "predicted": predicted},
                          genes=list(shifted["gene_id"][:100]), numbers={"oof_spearman": rho})


def _shift_test(ctx, p):
    cond, base = _need(ctx, p["condition"]), _need(ctx, p["baseline"])
    y = _shift(ctx, cond, base)
    X, _cols = _regression_matrix(ctx, [cond, base], p["model"], p["own_kind"])
    return S.value_test(ctx, "condition_shift", y, lambda vis: _fit_predict(X, vis, p["model"]),
                        n_null=3, min_effect=0.1, label=f"the {cond}-specific shift")


register(Strategy(
    key="condition_shift", number=23, family=LEARN,
    title="Find what matters more in one condition, and why",
    question="Which genes matter more (or less) in one condition than a baseline predicts -- in "
             "the mouse rather than the dish, say -- and can the rest of the data explain which?",
    tooltip="Takes a condition screen and its baseline, keeps the part of the condition the "
            "baseline cannot explain, and asks whether the other measurements predict that "
            "condition-specific effect -- the signature of what the parasite needs only there.",
    explanation=(
        "Most genes essential in the mouse are essential in the dish too; the interesting ones "
        "are the exceptions. This strategy isolates them: the condition measurement is regressed "
        "on its baseline (both rank-scaled) and the residual -- the condition-specific shift -- "
        "becomes the target. A gene with a strongly negative shift is needed in the condition "
        "beyond what its baseline need predicts: a candidate for host interaction, immune "
        "evasion, nutrient acquisition in vivo.\n\n"
        "Then it asks whether the shift is predictable from the other measurements, with both "
        "screens' closures withheld. If it is, the condition-specific need has a signature -- "
        "secretion, host-facing localization, a stage of expression -- and genes never screened in "
        "the condition can be ranked by it. If it is not, the shift is idiosyncratic or noise, "
        "and the self-test says which. The other screens of the same kind are left out by "
        "default: an in vivo shift predicted from another in vivo screen is one experiment "
        "confirming another, not an explanation."),
    walkthrough=(
        "Choose the condition (in vivo peritoneum by default) and its baseline (fibroblasts).",
        "Press Test: can hidden genes' shifts be predicted?",
        "Press Run; read 'shifted genes' from the top (most negative: needed more in the "
        "condition).",
        "If the test passed, read 'predicted' for genes not screened in the condition."),
    test_description=(
        "Pattern 4 on the shift: 20% of genes with both measurements hidden; a model trained on "
        "the rest predicts their shift. Metric: rank correlation on hidden genes. Null: 3 models "
        "trained on shuffled shifts. Pass: above the null's 95th percentile by 0.1."),
    params=(Param("condition", "number", "Condition",
                  "The measurement taken in the condition of interest, for example fitness in the "
                  "mouse. Its closure and the baseline's are withheld from the predictors.",
                  _condition_default("condition")),
            Param("baseline", "number", "Baseline",
                  "The same measurement in the reference condition, for example fitness in "
                  "fibroblasts. What the baseline explains is removed; what remains is "
                  "condition-specific.", _condition_default("baseline")),
            MODEL, OWN_KIND),
    runner=_shift_run, tester=_shift_test, cost="a minute", needs=("two numeric columns",)))


# =========================================================================== 5 · gene lists
def _profile(ctx, query, exclude) -> pd.DataFrame:
    """Everything that distinguishes a gene set: category enrichments, shifted measurements, and
    measured networks denser inside the set than chance -- one table, one correction."""
    from scipy.stats import mannwhitneyu, poisson
    query = np.unique(np.asarray(query, dtype=int))
    banned = ctx.banned(exclude)
    inq = np.zeros(ctx.n, bool)
    inq[query] = True
    rows = []
    for c in ctx.categorical_columns():
        if c in banned:
            continue
        t = ctx.truth(c)
        known = t.notna().to_numpy()
        draws = int((inq & known).sum())
        if draws < 3:
            continue
        total = int(known.sum())
        counts = t[known].value_counts()
        for cls, h in t[inq & known].value_counts().items():
            if h < 2:
                continue
            lift = (h / draws) / (counts[cls] / total)
            rows.append({"kind": "category", "feature": c, "value": str(cls), "n": int(h),
                         "effect": float(np.log2(lift)),
                         "p": S.hypergeom_sf(int(h), draws, int(counts[cls]), total)})
    for c in ctx.numeric_columns(banned):
        v = ctx.values(c)
        m = v.notna().to_numpy()
        a, b = v[inq & m].to_numpy(), v[~inq & m].to_numpy()
        if len(a) < 5 or len(b) < 5:
            continue
        u, pv = mannwhitneyu(a, b)
        auc = u / (len(a) * len(b))
        rows.append({"kind": "number", "feature": c, "value": "higher" if auc > 0.5 else "lower",
                     "n": len(a), "effect": float(auc - 0.5), "p": float(pv)})
    for l in ctx.measurement_layers(exclude):
        A = ctx.adjacency(l, "binary")
        internal = float(A[query][:, query].sum()) / 2.0
        possible = len(query) * (len(query) - 1) / 2.0
        density = (A.nnz / 2.0) / max(ctx.n * (ctx.n - 1) / 2.0, 1)
        expected = possible * density
        if expected <= 0:
            continue
        rows.append({"kind": "layer", "feature": l, "value": "denser inside",
                     "n": int(internal), "effect": float(np.log2((internal + 0.5) / (expected + 0.5))),
                     "p": float(poisson.sf(internal - 1, expected))})
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = S.bh(out["p"].to_numpy())
        out = out.sort_values(["q", "effect"], ascending=[True, False], kind="stable")
    return out.reset_index(drop=True)


def _profile_scores(ctx, profile: pd.DataFrame, query, q_max: float = 0.05) -> np.ndarray:
    """Score every gene by how well it matches a profile: its significant measurements, enriched
    categories and dense networks, each part standardised so none dominates by count."""
    sig = profile[profile["q"] < q_max] if len(profile) else profile
    parts = []
    num = sig[sig["kind"] == "number"] if len(sig) else sig
    if len(num):
        Z = ctx.matrix(list(num["feature"]))
        parts.append(Z @ (2 * num["effect"].to_numpy()) / len(num))
    cat = sig[(sig["kind"] == "category") & (sig["effect"] > 0)] if len(sig) else sig
    if len(cat):
        s = np.zeros(ctx.n)
        for r in cat.itertuples():
            s += r.effect * (ctx.truth(r.feature) == r.value).to_numpy()
        parts.append(s / len(cat))
    lay = sig[(sig["kind"] == "layer") & (sig["effect"] > 0)] if len(sig) else sig
    if len(lay):
        inq = np.zeros(ctx.n)
        inq[np.asarray(query, dtype=int)] = 1.0
        s = np.zeros(ctx.n)
        for r in lay.itertuples():
            A = ctx.adjacency(r.feature, "binary")
            deg = np.asarray(A.sum(axis=1)).ravel()
            s += (A @ inq) / np.maximum(deg, 1)
        parts.append(s / len(lay))
    if not parts:
        return np.zeros(ctx.n)
    return sum((p_ - p_.mean()) / (p_.std() or 1.0) for p_ in parts)


def _enrichment_run(ctx, p):
    pos, missing = ctx.resolve_genes(p["genes"])
    if len(pos) < 3:
        return StrategyResult("set_enrichment", "Paste at least three gene ids from this organism.",
                              {})
    prof = _profile(ctx, pos, p.get("exclude"))
    scores = _profile_scores(ctx, prof, pos)
    order = [g for g in np.argsort(-scores) if g not in set(pos)][: int(p["top"])]
    sig = prof[prof["q"] < 0.05] if len(prof) else prof
    summary = (f"{len(pos)} genes ({len(missing)} not found) tested against {len(prof):,} features "
               f"-- categories, measurements and networks, corrected together. {len(sig)} "
               f"distinguish the list at q < 0.05; strongest: "
               + ", ".join(f"{r.feature} {r.value}" for r in sig.head(6).itertuples())
               + ". 'resembles the profile' ranks every other gene by the same features.")
    return StrategyResult("set_enrichment", summary,
                          {"profile": prof, "resembles the profile": _genes_table(
                              ctx, order, score=scores[order]),
                           "not found": pd.DataFrame({"identifier": missing})},
                          genes=list(ctx.gene_ids[pos]))


def _enrichment_test(ctx, p):
    members, what, exclude = _list_or_example(ctx, p)
    rank = lambda q: _profile_scores(ctx, _profile(ctx, q, exclude), q)
    return S.set_expansion_test(ctx, "set_enrichment", members, rank, frac=0.4, n_null=10,
                                label=what)


register(Strategy(
    key="set_enrichment", number=24, family=LISTS,
    title="Describe what your gene list has in common",
    question="What distinguishes the genes on my list from the rest -- which categories are they "
             "enriched in, which measurements are shifted, which networks are dense among them?",
    tooltip="Tests a gene list against every category, every measurement and every measured "
            "network at once, corrected together; the significant features form a profile, and "
            "the profile is then used to rank every other gene -- which is how it is tested.",
    explanation=(
        "Enrichment analysis is usually run against one annotation at a time. This strategy runs "
        "it against everything the table holds, in one family of tests with one correction: "
        "hypergeometric enrichment for every class of every categorical column, a rank-sum test "
        "(reported as AUROC minus one half) for every measurement, and, for every measured network, "
        "whether the list has more internal edges than a random set of the same size would.\n\n"
        "A profile is only a description until it predicts something. So the significant features "
        "are combined into a score -- shifted measurements weighted by their effect, enriched "
        "categories by their log lift, networks by the share of a gene's partners on the list -- "
        "and every other gene is ranked by it. The self-test builds the profile from 60% of a "
        "known set and asks whether the other 40% rank near the top: a profile that recognises "
        "members it never saw has captured what the set has in common."),
    walkthrough=(
        "Paste the list; name the column it came from, if any, so it is not simply read back.",
        "Press Test.",
        "Press Run and read 'profile' by q-value: kind, feature, direction and effect.",
        "Read 'resembles the profile' for genes that share the list's characteristics."),
    test_description=(
        "Pattern 2: 40% of the set hidden; the profile is built from the other 60% and scores "
        "every gene. Metric: AUROC of the hidden members against every other non-query gene. "
        "Null: 10 random sets of the same size. Pass: above the null's 95th percentile by 0.1."),
    params=(GENES, EXCLUDE,
            Param("top", "int", "Genes to rank",
                  "How many of the genes that best match the profile, and are not on the list, to "
                  "show.", 300, lo=10, hi=10000, step=50)),
    runner=_enrichment_run, tester=_enrichment_test, needs=("a gene list",)))


def _expansion_operator(ctx, exclude, mode: str, k: int = 15):
    key = ("expansion", exclude, mode, k)
    if key not in ctx._cache:
        ops = [ctx.operator(l) for l in ctx.measurement_layers(exclude)]
        if mode == "networks + measurements" or not ops:
            X, _cols = ctx.features(exclude)
            ops.append(S.knn_operator(X, k))
        ctx._cache[key] = sum(ops) / len(ops)
    return ctx._cache[key]


def _expansion_run(ctx, p):
    pos, missing = ctx.resolve_genes(p["genes"])
    if len(pos) < 2:
        return StrategyResult("seed_expansion", "Paste at least two gene ids from this organism.", {})
    op = _expansion_operator(ctx, p.get("exclude"), p["mode"])
    scores = S.diffuse(op, pos, restart=float(p["restart"]))
    order = [g for g in np.argsort(-scores) if g not in set(pos)][: int(p["top"])]
    links = []
    for g in order:
        parts = []
        for l in ctx.measurement_layers(p.get("exclude")):
            hits = int(ctx.adjacency(l, "binary")[g, pos].sum())
            if hits:
                parts.append(f"{l} {hits}")
        links.append(", ".join(parts))
    table = _genes_table(ctx, order, score=scores[order], direct_links_to_list=links)
    summary = (f"{len(pos)} seeds diffused across {len(ctx.measurement_layers(p.get('exclude')))} "
               f"measured networks" + (" and the measurement space" if p["mode"] ==
                                       "networks + measurements" else "")
               + f"; the {len(table)} genes the walk visits most are listed with their direct links "
                 f"to the list. A gene can rank high with no direct link at all, through shared "
                 f"neighbours -- which is the point.")
    return StrategyResult("seed_expansion", summary,
                          {"candidates": table, "not found": pd.DataFrame({"identifier": missing})},
                          genes=list(table["gene_id"]))


def _expansion_test(ctx, p):
    members, what, exclude = _list_or_example(ctx, p)
    op = _expansion_operator(ctx, exclude, p["mode"])
    return S.set_expansion_test(ctx, "seed_expansion", members,
                                lambda q: S.diffuse(op, q, restart=float(p["restart"])),
                                n_null=10, label=what)


register(Strategy(
    key="seed_expansion", number=25, family=LISTS,
    title="Grow your gene list along the networks",
    question="Starting from my genes, which others does a walk across every measured network "
             "keep returning to?",
    tooltip="Seeds a random walk on your list and lets it wander across every permitted measured "
            "network (and, optionally, a nearest-neighbour graph of the measurements); genes the "
            "walk visits most are the list's closest relatives in the data.",
    explanation=(
        "A complex, a pathway, a secretory route: the genes that belong with yours are connected "
        "to them -- not necessarily directly, but by many short paths through co-expression, "
        "co-fitness, crosslinks and structural similarity. Random walk with restart measures "
        "exactly that: the walk starts on your seeds, steps along edges, and returns to the seeds "
        "with a fixed probability, so its long-run visiting frequency is high for genes close to "
        "the list by many routes and low for genes reached by a single long path.\n\n"
        "Each layer is degree-normalised, so hubs do not attract every walk, and the layers are "
        "averaged rather than merged, so a dense layer does not drown a sparse, precise one. With "
        "'networks + measurements' a k-nearest-neighbour graph of the measurements joins the "
        "average, which lets the walk reach genes no network covers."),
    walkthrough=(
        "Paste the list; name the column it came from, if any.",
        "Choose 'networks' for relational evidence only, or 'networks + measurements' to reach "
        "genes no network covers.",
        "Press Test, then Run.",
        "Read 'candidates': the score, and which layers link each one directly to the list."),
    test_description=(
        "Pattern 2: 30% of the set hidden; the walk is seeded from the other 70%. Metric: AUROC "
        "of the hidden members against every other non-seed gene. Null: 10 random seed sets of "
        "the same size. Pass: above the null's 95th percentile by 0.1."),
    params=(GENES, EXCLUDE,
            Param("mode", "choice", "Walk across",
                  "'networks' walks only the measured edge layers. 'networks + measurements' adds "
                  "a nearest-neighbour graph of the measurement table, so genes absent from every "
                  "network can still be reached.", "networks + measurements",
                  choices=("networks", "networks + measurements")),
            Param("restart", "float", "Restart probability",
                  "How often the walk jumps back to your seeds. High keeps it close (few, near "
                  "relatives); low lets it wander (more candidates, further away).", 0.3,
                  lo=0.05, hi=0.95, step=0.05),
            Param("top", "int", "Candidates to list",
                  "How many of the most-visited genes not on your list to show. Every gene is "
                  "scored by the walk; this only limits the table that is shown and saved.", 300,
                  lo=10, hi=10000, step=50)),
    runner=_expansion_run, tester=_expansion_test, needs=("a gene list",)))


# =========================================================================== 6 · contrast layers
def _second_label(ctx, first=None, numeric_ok=True):
    """A second label to contrast with the first: the stage, the phase, or the next one there is."""
    cats = ctx.categorical_columns()
    # Topology first where it exists: it is known for every gene, so half a sample still carries
    # enough of it for a split to be significant, which a label covering a tenth of genes does not.
    for c in ("dtm_class", "stage_enriched_derived", "cellcycle_phase", "export_pred_tier",
              "pb_transferred_phenotype", "has_signal_peptide"):
        if c in cats and c != first:
            return c
    return next((c for c in cats if c != first), None)


def _coarse_category(ctx):
    """The shared label for a split or a conjunction: the coarse localization where there is one.

    Twenty-six compartments are too fine for a cluster to be homogeneous in one; the eleven of
    `lopit_unified` are what a cluster can actually agree on.
    """
    return "lopit_unified" if "lopit_unified" in ctx.categorical_columns() else default_category(ctx)


def _contrast_setup(ctx, p):
    """(a, b, coords, positions, {selection: labels}, sub-table, b is numeric) for 26 and 27.

    Both selections by default, because the two claims need different grain: a SPLIT needs a whole
    category held in one cluster, which excess-of-mass gives on a small table and leaf gives on the
    real proteome, where excess-of-mass returns a handful of clusters too large to be about anything.
    """
    a, b = _need(ctx, p["a"]), _need(ctx, p["b"])
    if a == b:
        raise ValueError("choose two different columns to contrast")
    rows = ctx.sample_rows(p["sample"], target=a)
    coords, pos = ctx.embed(tuple(ctx.blocks((a, b))), 25, 0.1, rows, target=(a, b))
    clusterings = {sel: ctx.cluster(coords, int(p["min_cluster_size"]), sel)
                   for sel in parse_grid(p["selection"], str)}
    sub = ctx.nodes.iloc[pos].reset_index(drop=True).copy()
    sub[a] = ctx.truth(a).iloc[pos].to_numpy()
    numeric_b = pd.api.types.is_numeric_dtype(ctx.nodes[b]) and not pd.api.types.is_bool_dtype(
        ctx.nodes[b])
    sub[b] = ctx.values(b).iloc[pos].to_numpy() if numeric_b else ctx.truth(b).iloc[pos].to_numpy()
    return a, b, coords, pos, clusterings, sub, numeric_b


def _flat(table: pd.DataFrame) -> pd.DataFrame:
    """Lists in a findings table as text, so it can be shown and saved."""
    out = table.copy()
    for c in out.columns:
        if out[c].map(lambda v: isinstance(v, (list, tuple))).any():
            out[c] = out[c].map(lambda v: ", ".join(map(str, v)) if isinstance(v, (list, tuple))
                                else v)
    return out


def _pooled(fn, sub, clusterings: dict, *args, **kw) -> pd.DataFrame:
    """A discovery scan run on every clustering, its findings pooled and tagged with their source."""
    parts = []
    for sel, lab in clusterings.items():
        f = fn(sub, lab, *args, **kw)
        if len(f):
            parts.append(f.assign(selection=sel))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _halves(fn, sub, clusterings, a, b, **kw):
    """`find` for pattern 5: significant findings from genes in one half, on every clustering."""
    def find(half):
        out = []
        for sel, lab in clusterings.items():
            masked = np.full(len(lab), NOISE)
            masked[half] = lab[half]
            f = fn(sub, masked, a, b, **kw)
            # A finding is a claim of significance; the ones that do not make it are not findings.
            if len(f):
                out += [dict(r, selection=sel) for r in f[f["q"] < 0.05].to_dict("records")]
        return out
    return find


def _disagreement_run(ctx, p):
    from . import discovery
    a, b, coords, pos, clusterings, sub, numeric_b = _contrast_setup(ctx, p)
    tables = {f"same {a}, split {b}": _flat(_pooled(discovery.disagreement, sub, clusterings, a, b))}
    if not numeric_b:
        tables[f"same {b}, split {a}"] = _flat(_pooled(discovery.disagreement, sub, clusterings,
                                                       b, a))
    n = sum(len(t) for t in tables.values())
    first = next(iter(clusterings.values()))
    summary = (f"A map built without {a} or {b}, clustered {' and '.join(clusterings)}; {n} clusters "
               f"agree about one and split on the other. Each is a claim that a category contains "
               f"two kinds of gene, with the measurement that separates them.")
    return StrategyResult("split_clusters", summary, tables, labels=expand(first, pos, ctx.n),
                          coords=coords, positions=pos)


def split_replicates(finding: dict, labels, a_vals, b_vals, half, numeric_b: bool) -> bool:
    """Does one split finding hold among the genes of `half`?

    For a categorical B, the finding's minority group must be enriched again among the cluster's
    A-matching genes of this half (hypergeometric p < 0.05); for a measurement, those genes must
    again fall into two modes at least `discovery.MIN_SPLIT` apart.
    """
    from . import discovery
    inh = np.zeros(len(labels), bool)
    inh[np.asarray(half, dtype=int)] = True
    members = inh & (np.asarray(labels) == finding["cluster"]) & (a_vals == finding["category"])
    if numeric_b:
        b = np.asarray(b_vals, dtype=float)
        v = b[members]
        v = v[np.isfinite(v)]
        pool = b[inh]
        pool = pool[np.isfinite(pool)]
        if len(v) < 8 or len(pool) <= len(v):
            return False
        scale = float(np.std(pool))
        gap = discovery._split(v, scale)[0]
        # Two-means cuts ANY sample in two: plain normal noise lands about 1.6 standard deviations
        # apart, over the discovery threshold. So a split replicates only if it is also wider than
        # random sets of the same size from the same half manage.
        rng = np.random.default_rng(len(v))
        null = [discovery._split(rng.choice(pool, size=len(v), replace=False), scale)[0]
                for _ in range(50)]
        return gap >= discovery.MIN_SPLIT and gap > float(np.percentile(null, 95))
    g = finding["groups"][1]
    known = inh & np.array([isinstance(x, str) for x in b_vals])
    hits = int((members & known & (b_vals == g)).sum())
    draws = int((members & known).sum())
    return hits >= 2 and S.hypergeom_sf(hits, draws, int((known & (b_vals == g)).sum()),
                                       int(known.sum())) < 0.05


def _disagreement_test(ctx, p):
    from . import discovery
    a, b, _coords, _pos, clusterings, sub, numeric_b = _contrast_setup(ctx, p)
    a_vals = sub[a].to_numpy(dtype=object)
    b_vals = sub[b].to_numpy(dtype=float if numeric_b else object)
    n = len(sub)

    def replicate(f, half, state):
        return split_replicates(f, clusterings[f["selection"]], a_vals,
                                b_vals if state is None else state, half, numeric_b)

    def scramble(half, rng):
        out = b_vals.copy()
        out[half] = out[rng.permutation(half)]
        return out

    find = _halves(discovery.disagreement, sub, clusterings, a, b, min_cluster=8, min_category=4)
    return S.replication_test(ctx, "split_clusters", find, replicate, scramble, np.arange(n),
                              n_null=20, min_effect=0.2,
                              label=f"clusters that share a {a} and split on {b}")


register(Strategy(
    key="split_clusters", number=26, family=CONTRAST,
    title="Find categories that split in two on another measurement",
    question="Which clusters agree about one thing -- a compartment -- and split cleanly on "
             "another -- a stage, a phase, a fitness level?",
    tooltip="Finds clusters that are homogeneous for one label and divided on a second label or "
            "measurement; each is a claim that a category contains two kinds of gene, checked by "
            "asking whether the split reappears in genes the discovery never saw.",
    explanation=(
        "A compartment is not one thing: the nucleus holds genes of every cell-cycle phase, the "
        "apicoplast holds essential and dispensable proteins. Where a cluster is homogeneous for "
        "one label (A) and divides on another (B) -- a categorical label into two enriched groups, "
        "or a measurement into two modes -- the data is saying the category has internal "
        "structure, and naming the measurement that carries it.\n\n"
        "Both labels are withheld from the map, so the cluster is found without either and the "
        "claim is not circular. The map is clustered both ways HDBSCAN can select clusters, and the "
        "findings of both are pooled: a split needs a whole category held in one cluster, and "
        "which selection gives that depends on the table. Run in both directions for two "
        "categorical labels: 'same compartment, split phase' and 'same phase, split compartment' "
        "are different biology. The self-test makes the findings on half the genes and asks "
        "whether each split reappears in the other half."),
    walkthrough=(
        "Choose the label that should be shared (A) and the one that should split (B).",
        "Press Test: the share of first-half findings that replicate in the second half.",
        "Press Run and read both tables; 'groups' and 'counts' (or the two means) describe the "
        "split, 'genes' lists the minority side, 'selection' which clustering found it.",
        "Color the map by B and click into the cluster to see the split."),
    test_description=(
        "Pattern 5: significant findings (q < 0.05) on a random half of the mapped genes; each is "
        "checked on the other half -- the minority group enriched again among the cluster's "
        "A-matching genes (hypergeometric p < 0.05), or the measurement bimodal again. Metric: "
        "share replicating. Null: 20 runs with B shuffled within the second half. Pass: above the "
        "null's 95th percentile by 0.2 with at least three findings."),
    params=(Param("a", "category", "Shared label (A)",
                  "The label the cluster must agree on -- what its genes have in common. Withheld "
                  "from the map, with its closure. A coarse label works best: a cluster can agree "
                  "on an organelle far more often than on one of two dozen sub-compartments.",
                  _coarse_category),
            Param("b", "column", "Splitting label or measurement (B)",
                  "The label or measurement the cluster must divide on. A category splits into "
                  "enriched groups; a measurement into two modes. Withheld from the map too.",
                  lambda ctx: _second_label(ctx, _coarse_category(ctx))),
            SAMPLE,
            Param("min_cluster_size", "int", "Cluster size",
                  "HDBSCAN minimum cluster size. A split needs enough genes on both sides, so "
                  "larger clusters find more splits and smaller ones find purer clusters.", 40,
                  lo=5, hi=500),
            SELECTIONS),
    runner=_disagreement_run, tester=_disagreement_test, cost="a minute",
    needs=("two categorical columns",)))


def _conjunction_run(ctx, p):
    from . import discovery
    a, b, coords, pos, clusterings, sub, numeric_b = _contrast_setup(ctx, p)
    if numeric_b:
        raise ValueError("a conjunction needs two categorical labels")
    found = _pooled(discovery.conjunction, sub, clusterings, a, b)
    first = next(iter(clusterings.values()))
    summary = (f"A map built without {a} or {b}; {len(found)} clusters are enriched for a "
               f"COMBINATION of the two beyond what either label alone explains (interaction ratio "
               f">= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted "
               f"AND fitness-conferring in vivo.")
    return StrategyResult("conjunctions", summary, {"conjunctions": _flat(found)},
                          labels=expand(first, pos, ctx.n), coords=coords, positions=pos)


def _conjunction_test(ctx, p):
    from . import discovery
    a, b, _coords, _pos, clusterings, sub, numeric_b = _contrast_setup(ctx, p)
    if numeric_b:
        raise ValueError("a conjunction needs two categorical labels")
    a_vals = sub[a].to_numpy(dtype=object)
    b_vals = sub[b].to_numpy(dtype=object)
    n = len(sub)

    def replicate(f, half, state):
        bv = b_vals if state is None else state
        labels = clusterings[f["selection"]]
        inh = np.zeros(n, bool)
        inh[half] = True
        base = inh & (a_vals == f["category"]) & np.array([isinstance(x, str) for x in bv])
        here = base & (labels == f["cluster"])
        hits = int((here & (bv == f["other_category"])).sum())
        return hits >= 2 and S.hypergeom_sf(hits, int(here.sum()),
                                           int((base & (bv == f["other_category"])).sum()),
                                           int(base.sum())) < 0.05

    def scramble(half, rng):
        out = b_vals.copy()
        out[half] = out[rng.permutation(half)]
        return out

    find = _halves(discovery.conjunction, sub, clusterings, a, b, min_cluster=8, min_category=4)
    return S.replication_test(ctx, "conjunctions", find, replicate, scramble, np.arange(n),
                              n_null=20, min_effect=0.2, label=f"{a} x {b} combinations")


register(Strategy(
    key="conjunctions", number=27, family=CONTRAST,
    title="Find kinds of gene defined by two labels at once",
    question="Which clusters are enriched for a COMBINATION of two labels -- more than either "
             "label alone would make them?",
    tooltip="Finds clusters enriched for a pair of labels beyond what each label's own enrichment "
            "explains; each is a kind of gene defined by both at once, checked by asking whether "
            "the second label stays enriched within the first in genes the discovery never saw.",
    explanation=(
        "Some biology lives in combinations. Secreted proteins are common and fitness-conferring "
        "genes are common, but a cluster where the secreted genes are ALSO the fitness-conferring "
        "ones is a finding neither label shows alone. The strategy withholds both labels, clusters "
        "a map, and for every cluster and every combination measures the joint enrichment against "
        "the stronger of the two single enrichments; a ratio above 1.25 means the combination says "
        "more than either margin.\n\n"
        "The replication test checks the interaction itself rather than the enrichment: within the "
        "second half's genes of the first label, is the second label still concentrated in this "
        "cluster? A cluster enriched for the first label alone cannot pass that, because the second "
        "label is then no more common inside the cluster than among the first label's genes "
        "anywhere."),
    walkthrough=(
        "Choose two categorical labels.",
        "Use a smaller cluster size than for strategy 26: a combination is a finer unit.",
        "Press Test, then Run; read 'conjunctions' by q with the interaction ratio.",
        "The 'genes' column lists unmeasured members -- candidates for the combination."),
    test_description=(
        "Pattern 5: conjunctions found on a random half; each is checked on the other half as the "
        "second label's enrichment in the cluster among genes carrying the first label "
        "(hypergeometric p < 0.05). Metric: share replicating. Null: 20 runs with the second "
        "label shuffled within the second half. Pass: above the null's 95th percentile by 0.2 "
        "with at least three findings."),
    params=(Param("a", "category", "First label",
                  "One of the two labels whose combination is sought. Withheld from the map with "
                  "its closure; a coarse label such as the unified localization works best.",
                  _coarse_category),
            Param("b", "category", "Second label",
                  "The other label. Withheld from the map too, so the cluster is found without "
                  "either and the combination cannot be read back from the input.",
                  lambda ctx: _second_label(ctx, _coarse_category(ctx))),
            SAMPLE,
            Param("min_cluster_size", "int", "Cluster size",
                  "HDBSCAN minimum cluster size. Combinations are finer than single labels, so "
                  "smaller clusters than for a split are usually right.", 15, lo=5, hi=500),
            SELECTIONS),
    runner=_conjunction_run, tester=_conjunction_test, cost="a minute",
    needs=("two categorical columns",)))


def _paralog_pairs(ctx) -> np.ndarray:
    import scipy.sparse as sp
    if "orthogroup" in ctx.layers():
        U = sp.triu(ctx.adjacency("orthogroup", "binary"), 1).tocoo()
        return np.column_stack([U.row, U.col]).astype(int)
    from .embedding import as_text
    og = as_text(ctx.nodes["orthogroup"]) if "orthogroup" in ctx.nodes else pd.Series([""] * ctx.n)
    pairs = []
    for _g, idx in og[og != ""].groupby(og[og != ""]).groups.items():
        idx = list(idx)
        if 2 <= len(idx) <= 10:
            pairs += [(i, j) for k, i in enumerate(idx) for j in idx[k + 1:]]
    return np.array(pairs, dtype=int).reshape(-1, 2)


def _divergence(ctx, pairs, target, min_shared: int) -> np.ndarray:
    cols = ctx.numeric_columns(ctx.banned(target))
    M = ctx.matrix(cols, impute=False)
    out = np.full(len(pairs), np.nan)
    for i, (a, b) in enumerate(pairs):
        m = np.isfinite(M[a]) & np.isfinite(M[b])
        if m.sum() >= min_shared:
            x, y = M[a, m], M[b, m]
            if np.ptp(x) > 0 and np.ptp(y) > 0:
                out[i] = 1.0 - float(np.corrcoef(x, y)[0, 1])
    return out


def _paralog_run(ctx, p):
    target = p.get("target")
    pairs = _paralog_pairs(ctx)
    if not len(pairs):
        return StrategyResult("paralog_divergence", "This table has no paralog pairs.", {})
    div = _divergence(ctx, pairs, target, int(p["min_shared"]))
    table = pd.DataFrame({"gene_a": ctx.gene_ids[pairs[:, 0]], "product_a": ctx.product(pairs[:, 0]),
                          "gene_b": ctx.gene_ids[pairs[:, 1]], "product_b": ctx.product(pairs[:, 1]),
                          "divergence": div})
    if target and target in ctx.nodes:
        t = ctx.truth(target)
        table[f"{target}_a"] = t.iloc[pairs[:, 0]].to_numpy()
        table[f"{target}_b"] = t.iloc[pairs[:, 1]].to_numpy()
    table = table.dropna(subset=["divergence"]).sort_values("divergence", ascending=False,
                                                            kind="stable").reset_index(drop=True)
    summary = (f"{len(table):,} paralog pairs scored by how differently they behave across the "
               f"permitted measurements (one minus the correlation of their profiles). The top of "
               f"the list is where gene duplication was followed by a change of job -- "
               f"sub- or neo-functionalisation -- and the bottom is redundancy.")
    return StrategyResult("paralog_divergence", summary, {"paralog pairs": table},
                          genes=list(dict.fromkeys(list(table["gene_a"][:100])
                                                   + list(table["gene_b"][:100]))))


def _paralog_test(ctx, p):
    t0 = time.monotonic()
    target = _need(ctx, p.get("target") or default_category(ctx))
    pairs = _paralog_pairs(ctx)
    div = _divergence(ctx, pairs, target, int(p["min_shared"])) if len(pairs) else np.array([])
    t = ctx.truth(target)
    la = t.iloc[pairs[:, 0]].to_numpy(dtype=object) if len(pairs) else np.array([])
    lb = t.iloc[pairs[:, 1]].to_numpy(dtype=object) if len(pairs) else np.array([])
    ok = np.array([isinstance(x, str) and isinstance(y, str) for x, y in zip(la, lb)], bool) \
        & np.isfinite(div) if len(pairs) else np.array([], bool)
    differ = la[ok] != lb[ok] if ok.any() else np.array([], bool)
    if ok.sum() < 20 or differ.all() or not differ.any():
        return S.judge("paralog_divergence", "AUROC", float("nan"), [], min_effect=0.05,
                       n_hidden=int(ok.sum()), hidden="paralog labels", null_kind="permuted",
                       t0=t0, note="too few paralog pairs with both genes labelled, or no contrast")
    rng = ctx.rng(81)
    d = div[ok]
    observed = S.auroc(d, differ)
    nulls = [S.auroc(rng.permutation(d), differ) for _ in range(20)]
    return S.judge("paralog_divergence",
                   f"AUROC of profile divergence for paralogs with different {target}",
                   observed, nulls, min_effect=0.05, n_hidden=int(ok.sum()),
                   hidden=f"the {target} labels, withheld from the profiles",
                   null_kind="20 random reassignments of divergence to pairs", t0=t0,
                   numbers={"pairs_that_differ": int(differ.sum())})


register(Strategy(
    key="paralog_divergence", number=28, family=CONTRAST,
    title="Find paralogs that changed jobs",
    question="Which duplicated genes behave differently across the measurements -- evidence that "
             "one copy took on a new role?",
    tooltip="Compares every paralog pair's measurement profiles and ranks them by divergence; "
            "diverged pairs are candidates for sub- or neo-functionalisation, and the test checks "
            "that divergence predicts a change of compartment the profiles never saw.",
    explanation=(
        "Gene duplication is where new functions come from, and apicomplexan genomes are full of "
        "expanded families. A pair of paralogs that behave alike -- same expression across stages, "
        "same fitness, same modifications -- is probably redundant. A pair that behave differently "
        "has probably partitioned or changed its job, and that is visible in the data long before "
        "anyone characterises either copy.\n\n"
        "Divergence is one minus the correlation of the pair's rank-scaled measurements over the "
        "columns both have (at least the chosen number). If a label is chosen, its closure is "
        "withheld from the profiles and the test asks the question that makes divergence "
        "meaningful: do paralogs that sit in DIFFERENT compartments have more divergent profiles "
        "than paralogs in the same one? If so, divergence in unlocalised pairs is evidence of a "
        "change of place or job."),
    walkthrough=(
        "Optionally choose a label (compartment); it is withheld from the profiles and shown "
        "beside each pair.",
        "Press Test: AUROC of divergence for pairs whose labels differ.",
        "Press Run and read 'paralog pairs' from the top: most diverged first.",
        "For a diverged pair, compare the two genes' evidence panels side by side."),
    test_description=(
        "Paralog pairs with both genes labelled; the labels are withheld from the profiles. "
        "Metric: AUROC of divergence for pairs whose labels differ against pairs whose labels "
        "match. Null: 20 random reassignments of the divergence values to pairs. Pass: above the "
        "null's 95th percentile by 0.05."),
    params=(Param("target", "category", "Label to compare (optional)",
                  "A label withheld from the profiles and shown beside each pair; the self-test "
                  "asks whether divergence predicts a difference in it. The default label is used "
                  "for the test when this is empty.", default_category, optional=True),
            Param("min_shared", "int", "Measurements both must have",
                  "A pair is scored only over measurements present for both genes, and only if "
                  "there are at least this many. Fewer shared measurements give noisier "
                  "divergence.", 10, lo=3, hi=500)),
    runner=_paralog_run, tester=_paralog_test, needs=("orthogroups",)))


# =========================================================================== 7 · species and strata
def _source_default(ctx):
    other = ctx.other()
    if other is None:
        return None
    nums = other.numeric_columns()
    for c in ("piggybac_mis", "fit_invitro_hff"):
        if c in nums:
            return c
    return nums[0] if nums else None


def _transfer_target_default(ctx):
    nums = ctx.numeric_columns()
    for c in ("fit_invitro_hff", "piggybac_mis"):
        if c in nums:
            return c
    return default_numeric(ctx)


def _through_orthologs(ctx, source: str) -> tuple:
    """(values mapped onto this organism's genes, source is numeric): through shared orthogroups."""
    from .embedding import as_text
    other = ctx.other()
    if other is None:
        raise ValueError("the other organism's table is not available here")
    if source not in other.nodes:
        raise ValueError(f"the other organism has no column {source!r}")
    if "orthogroup" not in ctx.nodes or "orthogroup" not in other.nodes:
        raise ValueError("both tables need an orthogroup column")
    og_this, og_other = as_text(ctx.nodes["orthogroup"]), as_text(other.nodes["orthogroup"])
    s = other.nodes[source]
    numeric = pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)
    if numeric:
        v = other.values(source)
        agg = v.groupby(og_other).median()
    else:
        v = other.truth(source)
        known = v.notna()
        agg = v[known].groupby(og_other[known]).agg(lambda x: x.value_counts().index[0])
    agg = agg[agg.index != ""]
    mapped = og_this.map(agg)
    return (mapped.astype(float) if numeric else mapped.astype(object)), numeric


def _target_numeric(ctx, target):
    s = ctx.nodes[target]
    return pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)


def _transfer_map(src: pd.Series, y: pd.Series, numeric_src: bool, numeric_y: bool):
    """Learn source -> target on the genes carrying both; return a function applying it."""
    both = src.notna() & y.notna()
    if numeric_y:
        if numeric_src:
            sign = np.sign(S.spearman(src[both], y[both])) or 1.0
            return lambda s: sign * s.astype(float)
        means = y[both].groupby(src[both]).mean()
        return lambda s: s.map(means).astype(float)
    if numeric_src:
        edges = np.nanquantile(src[both], [0.2, 0.4, 0.6, 0.8]) if both.any() else []
        bins = lambda s: pd.Series(np.digitize(s.astype(float), edges), index=s.index).where(
            s.notna()).astype(object)
        mode = y[both].groupby(bins(src[both])).agg(lambda x: x.value_counts().index[0])
        return lambda s: bins(s).map(mode)
    mode = y[both].groupby(src[both]).agg(lambda x: x.value_counts().index[0])
    return lambda s: s.map(mode)


def _transfer_run(ctx, p):
    target, source = _need(ctx, p["target"]), p["source"]
    src, numeric_src = _through_orthologs(ctx, source)
    numeric_y = _target_numeric(ctx, target)
    y = ctx.values(target) if numeric_y else ctx.truth(target)
    f = _transfer_map(src, y, numeric_src, numeric_y)
    pred = f(src)
    both = src.notna() & y.notna()
    agree = (S.spearman(pred[both], y[both]) if numeric_y
             else float((pred[both] == y[both]).mean()) if both.any() else float("nan"))
    rows = np.flatnonzero((src.notna() & y.isna()).to_numpy())
    table = _genes_table(ctx, rows, ortholog_value=src.iloc[rows].to_numpy(),
                         transferred=pred.iloc[rows].to_numpy())
    other = ctx.other()
    summary = (f"{source} carried from {other.organism if other else 'the other organism'} through "
               f"shared orthogroups to {int(src.notna().sum()):,} genes. Where both species are "
               f"measured ({int(both.sum()):,} genes), the transferred value "
               + (f"rank-correlates with {target} at {agree:.3f}" if numeric_y else
                  f"matches {target} for {agree:.0%}")
               + f". {len(table):,} genes lacking {target} receive a value from their orthologs.")
    return StrategyResult("ortholog_transfer", summary, {"transferred": table},
                          genes=list(table["gene_id"][:200]), numbers={"agreement": agree})


def _transfer_test(ctx, p):
    t0 = time.monotonic()
    target, source = _need(ctx, p["target"]), p["source"]
    src, numeric_src = _through_orthologs(ctx, source)
    numeric_y = _target_numeric(ctx, target)
    if not numeric_y:
        has = src.notna().to_numpy()
        predict = lambda vis: _transfer_map(src, vis, numeric_src, False)(src).where(src.notna())
        return S.label_transfer_test(ctx, "ortholog_transfer", target, predict, restrict=has,
                                     n_null=20)
    y = ctx.values(target)
    both = np.flatnonzero((src.notna() & y.notna()).to_numpy())
    rng = ctx.rng(91)
    hidden = rng.choice(both, size=int(round(0.25 * len(both))), replace=False) if len(both) \
        else np.array([], dtype=int)
    visible = y.copy()
    visible.iloc[hidden] = np.nan

    def score(s):
        pred = _transfer_map(s, visible, numeric_src, True)(s)
        return S.spearman(pred.iloc[hidden], y.iloc[hidden])

    observed = score(src)
    carried = np.flatnonzero(src.notna().to_numpy())
    nulls = []
    for _ in range(20):
        s = src.copy()
        s.iloc[carried] = s.iloc[rng.permutation(carried)].to_numpy()
        nulls.append(score(s))
    return S.judge("ortholog_transfer", "rank correlation of transferred and hidden values",
                   observed, nulls, min_effect=0.1, n_hidden=len(hidden),
                   hidden=f"25% of the genes measured in both species",
                   null_kind="20 runs with the orthologs' values dealt to random genes", t0=t0)


register(Strategy(
    key="ortholog_transfer", number=29, family=SPECIES,
    title="Carry what one parasite shows to the other",
    question="What does a gene's ortholog in the other parasite say about it -- its essentiality, "
             "its stage, its localization?",
    tooltip="Maps a measurement or label from the other species onto this one through shared "
            "orthogroups, learning how the two relate on genes measured in both; a Plasmodium "
            "knockout screen becomes a prediction for untested Toxoplasma genes, and back.",
    explanation=(
        "Toxoplasma and Plasmodium diverged hundreds of millions of years ago, but conserved genes "
        "often keep their jobs: a ribosomal protein essential in one is essential in the other. "
        "Each species has experiments the other lacks -- saturation mutagenesis in P. falciparum, "
        "hyperLOPIT in T. gondii -- so an ortholog is a second, independent measurement of the "
        "same gene.\n\n"
        "The two tables stay separate, as the project requires: orthology is a bridge, not a "
        "merge. The strategy aggregates the source column over each orthogroup in the other "
        "species (median, or the commonest label), learns how it relates to the target on genes "
        "measured in both -- a sign for two measurements, a lookup for labels -- and applies that to "
        "genes measured only on the other side. The test hides a quarter of the genes measured in "
        "both and compares the transfer with one whose ortholog values were dealt to random "
        "genes."),
    walkthrough=(
        "Choose the target in this organism and the source column in the other.",
        "Press Test: how well the ortholog's value predicts the hidden genes, against shuffled "
        "orthologs.",
        "Press Run; 'transferred' lists genes lacking the target that have an ortholog value.",
        "Remember what orthology cannot reach: lineage-specific genes have no ortholog -- see 30."),
    test_description=(
        "Numeric target: 25% of the genes measured in both species hidden; the relation is learned "
        "on the rest. Metric: rank correlation of transferred and hidden values. Null: 20 runs with "
        "the ortholog values permuted among genes. Categorical target: Pattern 1 on genes with an "
        "ortholog value. Pass: above the null's 95th percentile by 0.1 (numeric) or 0.05 "
        "(categorical)."),
    params=(Param("target", "column", "Target in this organism",
                  "The measurement or label to predict here. Numeric targets are scored by rank "
                  "correlation, labels by correct calls.", _transfer_target_default),
            Param("source", "column", "Source in the other organism",
                  "The measurement or label taken from the other species' orthologs. It is "
                  "aggregated over each orthogroup before it is compared with anything.",
                  _source_default, space="other")),
    runner=_transfer_run, tester=_transfer_test, needs=("the other organism's table",)))


STRATA = ("lineage-specific", "hypothetical protein", "understudied", "conserved")


def stratum_mask(ctx, name: str) -> np.ndarray:
    """Which genes belong to a named stratum of the proteome. Unknown membership is False."""
    n = ctx.n
    if name in ("lineage-specific", "conserved"):
        if "lineage_specific" in ctx.nodes:
            lin = ctx.nodes["lineage_specific"].fillna(False).astype(bool).to_numpy()
        elif "ortholog_number" in ctx.nodes:
            lin = (ctx.values("ortholog_number") == 0).to_numpy()
        else:
            lin = np.zeros(n, bool)
        return lin if name == "lineage-specific" else ~lin
    if name == "hypothetical protein":
        return pd.Series(ctx.product(np.arange(n))).str.contains("hypothetical", case=False,
                                                                 na=False).to_numpy()
    if name == "understudied":
        cols = [c for c in ("n_papers_focal", "n_papers_substantive") if c in ctx.nodes]
        if cols:
            return (sum(ctx.values(c).fillna(0) for c in cols) == 0).to_numpy()
        if "attention_depth" in ctx.nodes:
            return ctx.truth("attention_depth").isna().to_numpy()
        return np.ones(n, bool)
    raise ValueError(f"unknown stratum {name!r}; choose one of {', '.join(STRATA)}")


CONSERVATION = ("ortho", "lineage", "paralog", "conserv", "has_pf", "has_cp", "has_tg")


def _stratum_features(ctx, target):
    cols = [c for c in ctx.numeric_columns(ctx.banned(target))
            if not any(k in c for k in CONSERVATION)]
    return ctx.matrix(cols), cols


def _stratum_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    mask = stratum_mask(ctx, p["stratum"])
    X, cols = _stratum_features(ctx, target)
    query = np.flatnonzero(truth.isna().to_numpy() & mask)
    pred, share = S.knn_vote(X, truth, int(p["k"]), query=query)
    pred = pred.where(share >= float(p["min_share"]))
    calls = _calls_table(ctx, pred, share, truth)
    vis, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    hp, hs = S.knn_vote(X, vis, int(p["k"]), query=hidden)
    hp = hp.where(hs >= float(p["min_share"]))
    rows = [{"stratum": name, "hidden": int(m[hidden].sum()),
             "correct_rate": S.correct_rate(hp, truth, hidden[m[hidden]])}
            for name, m in (("in " + p["stratum"], mask), ("outside it", ~mask))]
    summary = (f"{int(mask.sum()):,} genes are {p['stratum']}; {len(calls):,} of those without a "
               f"{target} label are called from measurement neighbours, with conservation and "
               f"orthology columns removed so the call cannot lean on them. 'accuracy by stratum' "
               f"says whether these genes can be called as well as the rest.")
    return StrategyResult("stratum_focus", summary,
                          {"calls": calls, "accuracy by stratum": pd.DataFrame(rows)},
                          genes=list(calls["gene_id"][:200]))


def _stratum_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    mask = stratum_mask(ctx, p["stratum"])
    X, _cols = _stratum_features(ctx, target)
    k, floor = int(p["k"]), float(p["min_share"])

    def predict(vis):
        query = _unlabelled_but_known(vis, truth)
        pred, share = S.knn_vote(X, vis, k, query=query[mask[query]])
        return pred.where(share >= floor)

    return S.label_transfer_test(ctx, "stratum_focus", target, predict, restrict=mask, n_null=10)


register(Strategy(
    key="stratum_focus", number=30, family=SPECIES,
    title="Test inference on the genes orthology cannot reach",
    question="Can lineage-specific, hypothetical or understudied genes be called as reliably as "
             "the rest -- and what are they?",
    tooltip="Calls genes in one stratum -- lineage-specific, hypothetical, understudied -- from "
            "measurement neighbours with every conservation and orthology column removed, and "
            "measures the error rate on that stratum alone rather than borrowing the proteome's.",
    explanation=(
        "An accuracy measured over the whole proteome is dominated by conserved, well-studied "
        "genes, where every method does best. The genes that most need inference -- parasite-"
        "specific proteins with no ortholog, 'hypothetical proteins', genes nobody has published "
        "on -- are exactly where methods are weakest and where an overall accuracy flatters them "
        "most.\n\n"
        "This strategy calls genes within one stratum from their measurement neighbours, with "
        "conservation and orthology columns removed so a call cannot simply reflect how conserved "
        "a gene is, and the self-test scores ONLY hidden genes of that stratum. The result is an "
        "error rate that belongs to the genes being predicted. 'accuracy by stratum' puts it beside "
        "the rest of the proteome, which shows how much a global number would have overstated it."),
    walkthrough=(
        "Choose the label and the stratum (lineage-specific by default).",
        "Press Test: correct calls among hidden genes of the stratum only.",
        "Press Run; compare the two rows of 'accuracy by stratum'.",
        "Read 'calls' for the stratum's unlabelled genes -- trusted at the stratum's own rate."),
    test_description=(
        "Pattern 1 restricted to the stratum: 25% of the label hidden by whole orthogroups; only "
        "hidden genes inside the stratum are scored. Null: 10 runs on shuffled labels. Pass: "
        "above the null's 95th percentile by 0.05."),
    params=(TARGET,
            Param("stratum", "choice", "Stratum",
                  "Which genes to focus on: lineage-specific (no ortholog outside the lineage), "
                  "hypothetical proteins, understudied genes (no focal or substantive paper), or "
                  "conserved genes for comparison.", "lineage-specific", choices=STRATA),
            K, MIN_SHARE),
    runner=_stratum_run, tester=_stratum_test, needs=("a categorical column",)))


# =========================================================================== 8 · combine strategies
def _tri_sources(ctx, target, k: int) -> dict:
    X, _cols = ctx.features(target)
    out = {"measurements (kNN)": lambda vis, q: S.knn_vote(X, vis, k, q)[0],
           "classifier (logistic)": lambda vis, q: _logistic(X, vis, q)[0]}
    layers = ctx.measurement_layers(target)
    if layers:
        A = sum(ctx.adjacency(l, "binary") for l in layers)
        out["networks (vote)"] = lambda vis, q: S.graph_vote(A, vis, q)[0]
    return out


def _triangulate(sources: dict, vis: pd.Series, query, min_agree: int) -> tuple:
    query = np.asarray(query, dtype=int)
    preds = {name: fn(vis, query) for name, fn in sources.items()}
    pred = pd.Series([np.nan] * len(vis), dtype=object)
    agree = pd.Series(np.nan, index=range(len(vis)))
    for g in query:
        votes = pd.Series([p_.iloc[g] for p_ in preds.values() if isinstance(p_.iloc[g], str)])
        if len(votes):
            top = votes.value_counts()
            if top.iloc[0] >= min_agree:
                pred.iloc[g] = top.index[0]
                agree.iloc[g] = int(top.iloc[0])
    return pred, agree, preds


def _tri_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    sources = _tri_sources(ctx, target, int(p["k"]))
    query = np.flatnonzero(t.isna().to_numpy())
    pred, agree, preds = _triangulate(sources, t, query, int(p["min_agree"]))
    extra = {f"by {name}": pr.to_numpy(dtype=object) for name, pr in preds.items()}
    calls = _calls_table(ctx, pred, agree, truth, **extra)
    summary = (f"{len(sources)} independent strategies -- {', '.join(sources)} -- each called every "
               f"unlabelled gene; {len(calls):,} genes are called only where at least "
               f"{int(p['min_agree'])} agree. Agreement trades reach for trust: fewer calls, each "
               f"supported by different kinds of evidence.")
    return StrategyResult("triangulation", summary, {"calls": calls},
                          genes=list(calls["gene_id"][:200]))


def _tri_test(ctx, p, restrict=None, key="triangulation"):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    sources = _tri_sources(ctx, target, int(p["k"]))
    last = {}

    def predict(vis):
        pred, _agree, preds = _triangulate(sources, vis, _unlabelled_but_known(vis, truth),
                                           int(p["min_agree"]))
        last.setdefault("preds", preds)
        return pred

    result = S.label_transfer_test(ctx, key, target, predict, precision=True, n_null=5,
                                   min_effect=0.1, restrict=restrict)
    if "preds" in last:
        _vis, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
        if restrict is not None:
            hidden = hidden[np.asarray(restrict, bool)[hidden]]
        for name, pr in last["preds"].items():
            prec, n_calls = S.call_precision(pr, truth, hidden)
            result.numbers[f"{name} precision"] = prec
            result.numbers[f"{name} calls"] = n_calls
    return result


register(Strategy(
    key="triangulation", number=31, family=COMBINE,
    title="Call a gene only when independent strategies agree",
    question="Where do measurement neighbours, a trained classifier and the networks give the same "
             "answer -- and how much more often is that answer right?",
    tooltip="Runs three strategies built on different evidence -- measurement neighbours, a "
            "classifier, network partners -- and calls a gene only where enough of them agree; "
            "the test measures the precision of agreed calls against agreement between guesses.",
    explanation=(
        "Every strategy has its own failure mode: neighbours are misled by missingness, a "
        "classifier by a spurious weight, a network by a hub. Failures of different methods on "
        "different evidence are largely independent, so a call on which several agree is much "
        "less likely to be wrong than any one of them -- the logic of triangulation in "
        "measurement.\n\n"
        "The price is reach: an agreed call exists only where every method can speak, and the "
        "self-test therefore scores PRECISION (the share of agreed calls that are right) rather "
        "than calls per hidden gene, and reports each single method's precision and number of "
        "calls beside it. The null is agreement between the same methods trained on shuffled "
        "labels, which is how often chance alone produces a consensus."),
    walkthrough=(
        "Choose the label; keep 'must agree' at 2 of 3.",
        "Press Test: compare the agreed precision with each single method's (in the numbers).",
        "Press Run; 'calls' shows every method's own call beside the agreed one.",
        "Raise 'must agree' to 3 for the smallest, surest set."),
    test_description=(
        "Pattern 1 scored by precision: 25% of the label hidden; each method is trained on the "
        "visible genes and the agreed calls on hidden genes are scored. Metric: share of agreed "
        "calls that are correct. Null: 5 runs with visible labels shuffled. Pass: above the "
        "null's 95th percentile by 0.1. Single-method precisions are reported."),
    params=(TARGET, K,
            Param("min_agree", "int", "Methods that must agree",
                  "How many of the independent methods must give the same label before a gene is "
                  "called. Two of three is the usual balance; three of three is the surest and "
                  "reaches the fewest genes.", 2, lo=1, hi=3)),
    runner=_tri_run, tester=_tri_test, cost="a minute", needs=("a categorical column",)))


def _understudied_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    mask = stratum_mask(ctx, "understudied")
    sources = _tri_sources(ctx, target, int(p["k"]))
    query = np.flatnonzero(t.isna().to_numpy() & mask)
    pred, agree, _preds = _triangulate(sources, t, query, int(p["min_agree"]))
    cols = [c for c in ("n_papers_focal", "n_papers_substantive") if c in ctx.nodes]
    attention = sum(ctx.values(c).fillna(0) for c in cols) if cols else pd.Series(0.0,
                                                                                  index=range(ctx.n))
    novelty = 1.0 - attention.rank(pct=True)
    priority = (agree / len(sources)) * novelty
    calls = _calls_table(ctx, pred, priority, truth, sources_agreeing=agree.to_numpy())
    summary = (f"{int(mask.sum()):,} genes have no focal or substantive publication; "
               f"{len(calls):,} of those without a {target} label receive an agreed call, ranked by "
               f"agreement times novelty. These are the cheapest discoveries: a hypothesis the data "
               f"supports about a gene nobody has written about.")
    return StrategyResult("understudied_first", summary, {"candidates": calls},
                          genes=list(calls["gene_id"][:200]))


def _understudied_test(ctx, p):
    return _tri_test(ctx, p, restrict=stratum_mask(ctx, "understudied"), key="understudied_first")


register(Strategy(
    key="understudied_first", number=32, family=COMBINE,
    title="Put the understudied genes first",
    question="Which genes nobody has written about can the data say something trustworthy about?",
    tooltip="Makes agreed calls for genes with no focal or substantive publication and ranks them "
            "by agreement times novelty; its test scores precision on understudied genes alone, "
            "so the trust it reports is earned where the attention is not.",
    explanation=(
        "Research attention is concentrated: a few hundred genes carry most of the literature, "
        "and thousands have never been the subject of a paper. Those are where a knowledge map "
        "adds the most -- and where it must be most careful, because an inference about an "
        "unstudied gene will not be contradicted by anything anyone has written.\n\n"
        "This strategy restricts the triangulated calls of strategy 31 to understudied genes (no "
        "focal or substantive paper), ranks them by how many methods agree times how little the "
        "gene has been written about, and measures precision on hidden understudied genes only -- "
        "never borrowing trust from the well-studied genes, where every method looks better "
        "because the measurements themselves were designed around them. A call here is a first "
        "hypothesis about a gene with no literature, and the precision beside it is the only "
        "evidence for it there is."),
    walkthrough=(
        "Choose the label.",
        "Press Test: precision of agreed calls on hidden understudied genes only.",
        "Press Run; 'candidates' is ranked by agreement times novelty.",
        "Take the top candidates to their evidence panels; each is a first hypothesis about a "
        "gene with no literature."),
    test_description=(
        "Pattern 1 scored by precision and restricted to understudied genes: 25% of the label "
        "hidden; agreed calls are scored only on hidden understudied genes. Null: 5 runs with "
        "visible labels shuffled. Pass: above the null's 95th percentile by 0.1."),
    params=(TARGET, K,
            Param("min_agree", "int", "Methods that must agree",
                  "How many of the three methods must agree before an understudied gene is called; "
                  "these genes have no literature to catch a wrong call, so agreement matters "
                  "more here than anywhere.", 2, lo=1, hi=3)),
    runner=_understudied_run, tester=_understudied_test, cost="a minute",
    needs=("a categorical column", "a literature-attention column")))
