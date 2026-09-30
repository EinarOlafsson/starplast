"""Gene-gene links with their provenance: what the star map draws, and where each link came from.

The 3D map places genes; it does not say WHY two genes are related. The strategies do, but each one
says it in its own table -- a list of predicted pairs, a module assignment, a set of candidates grown
from a seed list, a call backed by partners. This module turns all of those into one kind of thing,
a link between two genes, and keeps with every link the answer to "who says so":

    kind      "measured"  an edge layer of the space's graph (`graph.npz`: crosslinks, co-fitness,
                          co-expression, shared domain, co-mention ...)
              "inferred"  a link a strategy derived from those measurements
    source    the layer key, or the strategy key
    how       what kind of statement the link is (see `HOW`)
    setting   the settings the strategy ran with, as text
    score     the strategy's own number for the link (probability, logit, count, residual ...)
    strength  that score as a percentile within its run, 0..1, so links from different strategies
              can share one width scale without pretending their scores are comparable
    run       which run produced it: "" for a measured layer, ``starplast:<key>`` for the runs shipped
              with the package, ``you:<key>:<timestamp>`` for a run made in the application
    origin    "data", "starplast" or "you"

Three stores feed the map. MEASURED links are read from the space's own graph on demand and never
copied. SHIPPED links are the default-setting runs of the strategies below on the packaged tables,
precomputed by `scripts/build_star_edges.py` into ``data/star_edges.parquet`` (recorded as an
executed notebook). USER links are written by :class:`UserEdgeStore` whenever a strategy is run in
the application, one file per run under the user cache, and appear in the map at once.

**Everything here is bounded, and every bound is a named constant.** A strategy that says "these
1,800 genes form a community" is not saying there are 1.6 million pairwise links; drawing it as a
clique would bury every other source. A module becomes, for each member, links to its `MODULE_K`
nearest co-members -- nearest first by how many measured layers already join them, then by distance
on the map the strategy built (or in the permitted measurements when it built none). A set expansion
links each candidate to its `SEED_K` nearest seeds by the same rule. A neighbour vote (07) links a
called gene to its `LABEL_K` nearest labelled genes that carry the label it was called with, found
in the same measurement space the vote used. A partner vote (13) links a call to the partners it
names, at most `PARTNER_MAX`. No run contributes more than `MAX_EDGES_PER_RUN` links.

Qt-free, so the build script and the tests use it without a display.
"""
from __future__ import annotations

import datetime
import json
import os
import re

import numpy as np
import pandas as pd

MEASURED, INFERRED = "measured", "inferred"
ORIGIN_DATA, ORIGIN_SHIPPED, ORIGIN_USER = "data", "starplast", "you"
#: The group every run made in the application falls in, whatever its strategy: one toggle and one
#: colour, so "what did MY runs add" is a single switch.
YOUR_RUNS = "your runs"

#: The shipped file, in the package's data directory.
SHIPPED_FILE = "star_edges.parquet"

#: Links per module member, per set-expansion candidate, per neighbour-vote call, per partner call.
MODULE_K = 3
SEED_K = 3
LABEL_K = 3
PARTNER_MAX = 5
#: The most links one run may add, strongest first.
MAX_EDGES_PER_RUN = 20000
#: Longest provenance note kept per link.
NOTE_MAX = 90

#: What kind of statement a link is.
HOW = {
    "layer": "measured in this layer",
    "pair": "predicted or ranked pair",
    "module": "same module or cluster, nearest co-members",
    "seed": "grown from a seed list, nearest seeds",
    "label": "called by its nearest labelled genes",
    "partner": "called through a named partner",
}

#: The strategies whose default-setting runs are shipped, in catalogue order. Every one relates
#: genes to genes; strategies that only rank or call genes one at a time add nothing to a network.
SHIPPED_STRATEGIES = (
    "holdout_search",        # 01 clusters of its best map
    "consensus_modules",     # 04 modules that survive the walk
    "feature_knn",           # 07 called gene -> nearest labelled supporters
    "cluster_guilt",         # 09 enriched clusters
    "physical_partners",     # 13 called gene -> named partners
    "multiplex_modules",     # 15 communities several layers agree on
    "link_prediction",       # 16 predicted contacts
    "attention_correction",  # 17 literature pairs, fame removed
    "unwritten_links",       # 18 multi-layer pairs never written about
    "positive_unlabeled",    # 20 seed list -> resembling genes
    "seed_expansion",        # 25 seed list -> genes the walk visits
    "neighbour_space",       # 33 pairs of the integrated space
    "network_training",      # 34 edges the networks are missing
)
#: How many example seed lists the set-expansion strategies are shipped with (the categories of the
#: default label whose size is closest to what `strategies.example_set` aims for).
SHIPPED_SEED_SETS = 3

#: Tables that restate measured edges or duplicate another table of the same run.
PAIR_TABLE_SKIP = frozenset({"measured edges, ranked", "raw ranking"})
#: The first of these a pair table carries is its score.
PAIR_SCORE_COLUMNS = ("probability", "score", "corrected", "layers", "raw_count", "support")
#: The first of these a pair table carries becomes the note.
PAIR_NOTE_COLUMNS = ("measured_in", "linked_in", "top_evidence")

#: What each measured layer is, for the legend and the hover text. Layers not named here show
#: their key.
LAYER_TEXT = {
    "comention": "co-mention in abstracts (literature)",
    "comention_ft": "co-mention in open-access full texts (literature)",
    "orthogroup": "shared orthogroup",
    "coexpression": "co-expression (stage series)",
    "cotranslation": "co-translation",
    "compartment": "shared compartment (hyperLOPIT)",
    "cofitness": "co-fitness (CRISPR screens)",
    "domain": "shared InterPro domain",
    "xlms": "crosslink MS: measured physical proximity",
    "ip_ms": "IP-MS: replicated pulldown of a tagged bait",
    "struct": "structural similarity (Foldseek)",
    "structural_hole": "structural hole (derived: biology links them, literature does not)",
    "unwritten_interaction": "measured to bind, never written about (derived)",
}

#: Columns of an edge frame in memory (`a`, `b` are positions in the node table).
COLUMNS = ("a", "b", "kind", "group", "source", "how", "setting", "score", "strength", "run",
           "origin", "created", "note")
#: Columns of a stored edge file (gene ids instead of positions, and the organism).
STORED = ("organism", "gene_a", "gene_b", "kind", "group", "source", "how", "setting", "score",
          "strength", "run", "origin", "created", "note")


# --------------------------------------------------------------------------- small helpers
def empty() -> pd.DataFrame:
    """An edge frame with no rows and every column."""
    return pd.DataFrame({c: pd.Series(dtype=("int32" if c in ("a", "b") else
                                              "float32" if c in ("score", "strength") else object))
                         for c in COLUMNS})


def strength_of(scores) -> np.ndarray:
    """Scores as percentiles within their own set, 0..1 (missing scores sit at the middle).

    A logit, a probability, a crosslink count and a residual cannot share a width scale; their
    ranks can. One link alone gets 1.0.
    """
    s = pd.Series(np.asarray(scores, dtype=float))
    if len(s) == 0:
        return np.zeros(0, dtype=np.float32)
    r = s.rank(pct=True, method="average")
    return r.fillna(0.5).to_numpy(dtype=np.float32)


def settings_text(settings: dict, limit: int = 120) -> str:
    """A run's settings as one short line; a pasted gene list is shown as its size."""
    parts = []
    for k, v in (settings or {}).items():
        if k in ("genes", "exclude") and v:
            n = len([t for t in re.split(r"[^A-Za-z0-9_.\-]+", str(v)) if t])
            if n:
                parts.append(f"{k}={n} genes")
            continue
        if v is None or v == "":
            continue
        parts.append(f"{k}={v}")
    text = ", ".join(parts)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def layer_text(layer: str) -> str:
    """What a measured layer is, in words."""
    return LAYER_TEXT.get(layer, layer)


def strategy_title(key: str) -> str:
    """'16 link_prediction: Predict the contacts ...' for a strategy key, or the key alone."""
    try:
        from . import strategies as S
        for mod in ("strategy_catalog", "strategy_graph", "strategy_learning"):
            __import__(f"{__package__}.{mod}")
        s = S.get(key)
        return f"{s.number:02d} {key}: {s.title}"
    except Exception:
        return key


def _dedupe(frame: pd.DataFrame) -> pd.DataFrame:
    """One undirected link per (pair, run, how), keeping the higher score; no self links."""
    if frame.empty:
        return frame
    a = frame["a"].to_numpy()
    b = frame["b"].to_numpy()
    frame = frame.assign(a=np.minimum(a, b).astype(np.int32), b=np.maximum(a, b).astype(np.int32))
    frame = frame[frame["a"] != frame["b"]]
    frame = frame.sort_values("score", ascending=False, kind="stable", na_position="last")
    return frame.drop_duplicates(["a", "b", "run", "how"]).reset_index(drop=True)


def _frame(a, b, score, how, note=None) -> pd.DataFrame:
    n = len(a)
    return pd.DataFrame({"a": np.asarray(a, dtype=np.int32), "b": np.asarray(b, dtype=np.int32),
                         "score": np.asarray(score, dtype=np.float32), "how": [how] * n,
                         "note": list(note) if note is not None else [""] * n})


# --------------------------------------------------------------------------- measured layers
def load_graph(code: str, gene_ids) -> dict:
    """A space's edge layers as `layer__a/b/w` arrays, or {} if the graph is over another table."""
    from . import organisms
    path = organisms.graph_path(code)
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=False)
    ids = z["gene_ids"].astype(str) if "gene_ids" in z.files else None
    if ids is None or len(ids) != len(gene_ids) or not (ids == np.asarray(gene_ids, str)).all():
        return {}
    return {k: z[k] for k in z.files if k.rsplit("__", 1)[-1] in ("a", "b", "w")}


def measured_edges(graph: dict) -> pd.DataFrame:
    """Every measured layer as links: the layer is the source and its weight the score."""
    parts = []
    for key in sorted(k for k in graph if k.endswith("__a")):
        layer = key[:-3]
        a, b = np.asarray(graph[key]), np.asarray(graph[f"{layer}__b"])
        if not len(a):
            continue
        w = np.asarray(graph.get(f"{layer}__w", np.ones(len(a))), dtype=np.float32)
        f = _frame(a, b, w, "layer")
        f["strength"] = strength_of(w)
        f["kind"], f["group"], f["source"] = MEASURED, layer, layer
        f["setting"], f["run"], f["origin"], f["created"] = "", "", ORIGIN_DATA, ""
        parts.append(f)
    if not parts:
        return empty()
    return pd.concat(parts, ignore_index=True)[list(COLUMNS)]


# --------------------------------------------------------------------------- inferred links
def _multiplex(ctx):
    """How many measurement layers join each pair (sparse, symmetric), or None without a graph."""
    cache = getattr(ctx, "_cache", {})
    if "star_multiplex" not in cache:
        M = None
        try:
            layers = ctx.measurement_layers()
        except Exception:
            layers = []
        for layer in layers:
            A = ctx.adjacency(layer, "binary")
            M = A if M is None else M + A
        cache["star_multiplex"] = M.tocsr() if M is not None else None
    return cache["star_multiplex"]


def _space(ctx, result=None, target=None) -> np.ndarray:
    """Where distances are measured: the result's own map where it has one, else the measurements.

    Rows are genes. Genes the map did not place are NaN, and a group containing one falls back to the
    measurements (`_space_for`).
    """
    if result is not None and result.coords is not None and result.positions is not None:
        coords = np.asarray(result.coords, dtype=float)
        pos = np.asarray(result.positions, dtype=int)
        if coords.ndim == 2 and len(coords) == len(pos):
            out = np.full((ctx.n, coords.shape[1]), np.nan)
            out[pos] = coords
            return out
    X, _ = ctx.features(target)
    return np.asarray(X, dtype=float)


def _space_for(ctx, space, members, target=None) -> np.ndarray:
    """The rows of `space` for these genes, or the measurement rows if any of them is unplaced."""
    X = space[members]
    if X.shape[1] and np.isfinite(X).all():
        return X
    F, _ = ctx.features(target)
    return np.asarray(F, dtype=float)[members]


def _nearest(ctx, M, space, query: np.ndarray, pool: np.ndarray, k: int, target=None) -> list:
    """For each query gene, up to `k` genes of `pool`: most shared measured layers first, then nearest.

    Returns (query, partner, layers, distance) tuples. A gene is never its own neighbour.
    """
    from sklearn.neighbors import NearestNeighbors
    pool = np.asarray(pool, dtype=int)
    query = np.asarray(query, dtype=int)
    if len(pool) == 0 or len(query) == 0 or k <= 0:
        return []
    both = np.union1d(query, pool)
    X = _space_for(ctx, space, both, target)
    row = {g: i for i, g in enumerate(both)}
    Xp = X[[row[g] for g in pool]]
    Xq = X[[row[g] for g in query]]
    nk = int(min(len(pool), k + 1))
    if Xp.shape[1]:
        dist, idx = NearestNeighbors(n_neighbors=nk).fit(Xp).kneighbors(Xq)
    else:
        dist = np.zeros((len(query), nk))
        idx = np.tile(np.arange(nk), (len(query), 1))
    sub = M[query][:, pool].tocsr() if M is not None else None
    out = []
    for r, g in enumerate(query):
        chosen, seen = [], {g}
        if sub is not None:
            s, e = sub.indptr[r], sub.indptr[r + 1]
            cols, vals = sub.indices[s:e], sub.data[s:e]
            for o in np.argsort(-vals, kind="stable"):
                p = int(pool[cols[o]])
                if p in seen:
                    continue
                d = float(np.linalg.norm(X[row[g]] - X[row[p]])) if X.shape[1] else 0.0
                chosen.append((g, p, float(vals[o]), d))
                seen.add(p)
                if len(chosen) >= k:
                    break
        for j, d in zip(idx[r], dist[r]):
            if len(chosen) >= k:
                break
            p = int(pool[j])
            if p in seen:
                continue
            chosen.append((g, p, 0.0, float(d)))
            seen.add(p)
        out.extend(chosen)
    return out


def _closeness(layers: float, dist: float) -> float:
    """Shared layers first, then nearness: `layers + 1/(1+distance)`."""
    return float(layers) + 1.0 / (1.0 + float(dist))


def _why(layers: float, where: str) -> str:
    return (f"joined in {int(layers)} measured layer{'s' if layers != 1 else ''}" if layers
            else f"nearest {where}")


def module_edges(labels, ctx, space, k: int = MODULE_K, target=None) -> pd.DataFrame:
    """Co-membership as bounded links: each member to its `k` nearest co-members (-1 = no module)."""
    labels = np.asarray(labels)
    M = _multiplex(ctx)
    a, b, score, note = [], [], [], []
    for lab in np.unique(labels[labels != -1]):
        members = np.flatnonzero(labels == lab)
        if len(members) < 2:
            continue
        for g, p, layers, d in _nearest(ctx, M, space, members, members, k, target):
            a.append(g)
            b.append(p)
            score.append(_closeness(layers, d))
            note.append(f"module {lab} ({len(members)} genes); {_why(layers, 'on its map')}")
    return _frame(a, b, score, "module", note)


def seed_edges(candidates, seeds, cand_scores, ctx, space, k: int = SEED_K) -> pd.DataFrame:
    """Set expansion as links: each candidate to its `k` nearest seeds, scored by its own rank."""
    M = _multiplex(ctx)
    by = dict(zip(np.asarray(candidates, dtype=int), np.asarray(cand_scores, dtype=float)))
    a, b, score, note = [], [], [], []
    for g, p, layers, d in _nearest(ctx, M, space, np.asarray(candidates, int),
                                    np.asarray(seeds, int), k):
        a.append(g)
        b.append(p)
        score.append(by.get(g, np.nan))
        note.append(f"seed of {len(seeds)}; {_why(layers, 'seed in the measurements')}")
    return _frame(a, b, score, "seed", note)


def label_edges(calls: pd.DataFrame, ctx, target: str, k: int = LABEL_K) -> pd.DataFrame:
    """A neighbour-vote call as links to its `k` nearest labelled genes that carry the called label.

    Found in the same space the vote used (`Context.features(target)`), so the links are the genes
    that made the call, not a new guess.
    """
    from sklearn.neighbors import NearestNeighbors
    truth = ctx.truth(target)
    X, _ = ctx.features(target)
    X = np.asarray(X, dtype=float)
    idx = ctx.index
    a, b, score, note = [], [], [], []
    if X.shape[1] == 0:
        return _frame(a, b, score, "label", note)
    for label, rows in calls.groupby("prediction", sort=False):
        pool = np.flatnonzero((truth == label).to_numpy())
        if not len(pool):
            continue
        q = np.array([idx.get(str(g).upper(), -1) for g in rows["gene_id"]])
        keep = q >= 0
        q, sup = q[keep], pd.to_numeric(rows["support"], errors="coerce").to_numpy()[keep]
        if not len(q):
            continue
        nn = NearestNeighbors(n_neighbors=int(min(k, len(pool)))).fit(X[pool])
        dist, nb = nn.kneighbors(X[q])
        for g, s, js in zip(q, sup, nb):
            for j in js:
                a.append(g)
                b.append(pool[j])
                score.append(s)
                note.append(f"called {label!r}; labelled {label!r}")
    return _frame(a, b, score, "label", note)


_PARTNER = re.compile(r"([A-Za-z0-9_.\-]+) \(([^()]*), ([0-9.]+)\)")


def partner_edges(calls: pd.DataFrame, ctx, k: int = PARTNER_MAX) -> pd.DataFrame:
    """A partner-vote call as links to the partners it names, weighted by their count."""
    idx = ctx.index
    a, b, score, note = [], [], [], []
    for g, text, pred in zip(calls["gene_id"], calls["partners"], calls.get(
            "prediction", pd.Series([""] * len(calls)))):
        i = idx.get(str(g).upper())
        if i is None:
            continue
        for m in _PARTNER.findall(str(text or ""))[:k]:
            j = idx.get(m[0].upper())
            if j is None:
                continue
            a.append(i)
            b.append(j)
            score.append(float(m[2]))
            note.append(f"called {pred!r}; partner labelled {m[1]!r}, count {m[2]}")
    return _frame(a, b, score, "partner", note)


def pair_edges(table: pd.DataFrame, ctx) -> pd.DataFrame:
    """A table of gene pairs (`gene_a`/`gene_b`, or `gene_id`/`neighbour`) as links."""
    idx = ctx.index
    if {"gene_a", "gene_b"} <= set(table.columns):
        ga, gb = table["gene_a"], table["gene_b"]
    elif {"gene_id", "neighbour"} <= set(table.columns):
        ga, gb = table["gene_id"], table["neighbour"]
    else:
        return _frame([], [], [], "pair")
    col = next((c for c in PAIR_SCORE_COLUMNS if c in table.columns), None)
    sc = (pd.to_numeric(table[col], errors="coerce").to_numpy() if col
          else np.linspace(1.0, 0.0, len(table)))           # a ranked list: its order is its score
    nc = next((c for c in PAIR_NOTE_COLUMNS if c in table.columns), None)
    notes = table[nc].fillna("").astype(str).tolist() if nc else [""] * len(table)
    a, b, score, note = [], [], [], []
    for x, y, s, t in zip(ga, gb, sc, notes):
        i, j = idx.get(str(x).upper()), idx.get(str(y).upper())
        if i is None or j is None:
            continue
        a.append(i)
        b.append(j)
        score.append(s)
        note.append((f"{nc.replace('_', ' ')}: {t}" if t else "") if nc else "")
    return _frame(a, b, score, "pair", note)


def edges_from_result(result, ctx, run: str = "", origin: str = ORIGIN_USER,
                      created: str = "", group: str | None = None) -> pd.DataFrame:
    """Every gene-gene link a strategy's result states, bounded, with its provenance.

    Reads what the result exposes: pair tables, a module or cluster assignment (`labels`), a
    candidate list grown from `settings['genes']`, a neighbour-vote's calls (07) and a partner
    vote's named partners. A result that relates no genes to each other gives an empty frame.

    :param run: the run id every link carries.
    :param origin: ``"starplast"`` for a shipped run, ``"you"`` for one made in the application.
    :param group: the toggle it is drawn under; defaults to the strategy key for a shipped run and to
        :data:`YOUR_RUNS` for the user's.
    """
    key = result.strategy
    parts = []
    tables = {n: pd.DataFrame(t) for n, t in (result.tables or {}).items()}
    for name, t in tables.items():
        if name in PAIR_TABLE_SKIP or t.empty:
            continue
        if {"gene_a", "gene_b"} <= set(t.columns) or {"gene_id", "neighbour"} <= set(t.columns):
            parts.append(pair_edges(t, ctx))
        elif {"gene_id", "partners"} <= set(t.columns):
            parts.append(partner_edges(t, ctx))
    settings = dict(result.settings or {})
    if key == "feature_knn" and "calls" in tables and settings.get("target"):
        calls = tables["calls"].dropna(subset=["prediction"])
        parts.append(label_edges(calls, ctx, settings["target"], k=LABEL_K))
    if settings.get("genes") and "candidates" in tables and not tables["candidates"].empty:
        seeds, _ = ctx.resolve_genes(settings["genes"])
        cand = tables["candidates"]
        pos = np.array([ctx.index.get(str(g).upper(), -1) for g in cand["gene_id"]])
        sc = pd.to_numeric(cand.get("score", pd.Series(np.linspace(1, 0, len(cand)))),
                           errors="coerce").to_numpy()
        keep = pos >= 0
        if len(seeds) and keep.any():
            parts.append(seed_edges(pos[keep], seeds, sc[keep], ctx, _space(ctx), k=SEED_K))
    if result.labels is not None and len(result.labels) == ctx.n:
        labels = np.asarray(result.labels)
        if (labels != -1).any():
            parts.append(module_edges(labels, ctx, _space(ctx, result), k=MODULE_K))
    parts = [p for p in parts if len(p)]
    if not parts:
        return empty()
    f = pd.concat(parts, ignore_index=True)
    f["run"] = run
    f = _dedupe(f)
    f["strength"] = np.zeros(len(f), dtype=np.float32)
    for how, rows in f.groupby("how").groups.items():
        f.loc[rows, "strength"] = strength_of(f.loc[rows, "score"])
    if len(f) > MAX_EDGES_PER_RUN:
        f = f.sort_values("strength", ascending=False, kind="stable").head(MAX_EDGES_PER_RUN)
    f["kind"], f["source"] = INFERRED, key
    f["group"] = group or (YOUR_RUNS if origin == ORIGIN_USER else key)
    f["setting"] = settings_text(settings)
    f["origin"], f["created"] = origin, created
    f["note"] = f["note"].astype(str).str.slice(0, NOTE_MAX)
    return f[list(COLUMNS)].reset_index(drop=True)


# --------------------------------------------------------------------------- storage
def to_stored(frame: pd.DataFrame, gene_ids, organism: str) -> pd.DataFrame:
    """Positions to gene ids, for writing: a stored file must survive a re-ordered node table."""
    ids = np.asarray(gene_ids, dtype=str)
    out = frame.drop(columns=["a", "b"]).copy()
    out.insert(0, "gene_b", ids[frame["b"].to_numpy(dtype=int)])
    out.insert(0, "gene_a", ids[frame["a"].to_numpy(dtype=int)])
    out.insert(0, "organism", organism)
    for c in ("score", "strength"):
        out[c] = out[c].astype(np.float32)
    return out[list(STORED)]


def from_stored(stored: pd.DataFrame, gene_ids, organism: str | None = None) -> pd.DataFrame:
    """A stored file back to positions in this node table; links to genes it lacks are dropped."""
    if stored is None or stored.empty:
        return empty()
    if organism is not None and "organism" in stored:
        stored = stored[stored["organism"].astype(str) == organism]
    index = {g: i for i, g in enumerate(np.asarray(gene_ids, dtype=str))}
    a = stored["gene_a"].astype(str).map(index)
    b = stored["gene_b"].astype(str).map(index)
    ok = (a.notna() & b.notna()).to_numpy()
    out = stored.loc[ok].drop(columns=["organism", "gene_a", "gene_b"], errors="ignore").copy()
    out.insert(0, "b", b[ok].astype(np.int32).to_numpy())
    out.insert(0, "a", a[ok].astype(np.int32).to_numpy())
    for c in COLUMNS:
        if c not in out:
            out[c] = ""
    return out[list(COLUMNS)].reset_index(drop=True)


def shipped_path() -> str:
    """Where the shipped links live (inside the package)."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", SHIPPED_FILE)


def shipped_edges(organism: str, gene_ids, path: str | None = None) -> pd.DataFrame:
    """The links precomputed from the shipped strategy runs, for one organism's table."""
    path = path or shipped_path()
    if not os.path.exists(path):
        return empty()
    return from_stored(pd.read_parquet(path), gene_ids, organism)


def run_id(key: str, when=None) -> str:
    """``you:<strategy>:<YYYYmmdd_HHMMSS>``, the id a user run's links carry."""
    when = when or datetime.datetime.now()
    return f"you:{key}:{when.strftime('%Y%m%d_%H%M%S')}"


class UserEdgeStore:
    """The links of the user's own strategy runs: one parquet and one JSON per run, kept on disk.

    Beside the other user state (`paths.user_cache_dir()`), never inside the package, so a rebuild or
    an upgrade does not lose them. The JSON carries what the hover text shows: the strategy, its
    settings, when it ran and its summary.
    """

    def __init__(self, root: str | None):
        """Keep runs under `root` (created on first write); `None` keeps them in memory only."""
        self.root = root
        self.memory: dict = {}

    def _path(self, run: str, ext: str) -> str:
        return os.path.join(self.root, re.sub(r"[^A-Za-z0-9._-]", "_", run) + ext)

    def add(self, result, ctx, organism: str, when=None) -> pd.DataFrame:
        """Derive a run's links, keep them, and return them (positions in `ctx`'s table)."""
        when = when or datetime.datetime.now()
        run = run_id(result.strategy, when)
        while run in self.memory or (self.root and os.path.exists(self._path(run, ".parquet"))):
            run += "b"                          # two runs in one second stay two runs
        created = when.isoformat(timespec="seconds")
        frame = edges_from_result(result, ctx, run=run, origin=ORIGIN_USER, created=created)
        meta = {"run": run, "strategy": result.strategy, "organism": organism, "created": created,
                "settings": settings_text(result.settings), "summary": str(result.summary)[:400],
                "links": int(len(frame))}
        self.memory[run] = (meta, frame)
        if self.root and len(frame):
            os.makedirs(self.root, exist_ok=True)
            to_stored(frame, ctx.gene_ids, organism).to_parquet(self._path(run, ".parquet"),
                                                                index=False)
            with open(self._path(run, ".json"), "w") as fh:
                json.dump(meta, fh, indent=1)
        return frame

    def runs(self, organism: str | None = None) -> list:
        """Every kept run's description, oldest first."""
        metas = {m["run"]: m for m, _ in self.memory.values()}
        if self.root and os.path.isdir(self.root):
            for f in sorted(os.listdir(self.root)):
                if f.endswith(".json"):
                    try:
                        with open(os.path.join(self.root, f)) as fh:
                            m = json.load(fh)
                        metas.setdefault(m["run"], m)
                    except (OSError, ValueError, KeyError):
                        continue
        out = [m for m in metas.values() if organism is None or m.get("organism") == organism]
        return sorted(out, key=lambda m: (m.get("created", ""), m["run"]))

    def load(self, organism: str, gene_ids) -> pd.DataFrame:
        """Every kept run's links for one organism's table (disk and this session)."""
        parts, seen = [], set()
        for run, (meta, frame) in self.memory.items():
            if meta.get("organism") == organism and len(frame):
                parts.append(frame)
                seen.add(run)
        for meta in self.runs(organism):
            if meta["run"] in seen or not self.root:
                continue
            p = self._path(meta["run"], ".parquet")
            if os.path.exists(p):
                try:
                    parts.append(from_stored(pd.read_parquet(p), gene_ids, organism))
                except Exception as exc:        # one unreadable run must not cost the others
                    print(f"starplast: could not read star-map run {p}: {exc}")
        parts = [p for p in parts if len(p)]
        return pd.concat(parts, ignore_index=True) if parts else empty()

    def remove(self, run: str) -> bool:
        """Forget one run, on disk and in memory."""
        found = self.memory.pop(run, None) is not None
        if self.root:
            for ext in (".parquet", ".json"):
                p = self._path(run, ext)
                if os.path.exists(p):
                    os.remove(p)
                    found = True
        return found


# --------------------------------------------------------------------------- the index
class EdgeIndex:
    """All links of one table, from any number of sources, looked up by gene in constant time.

    Blocks are added whole (measured, shipped, each user run) and the lookup is rebuilt lazily, so
    adding a run costs one concatenation, not a copy per gene.
    """

    def __init__(self, n: int):
        """An empty index over a table of `n` genes."""
        self.n = int(n)
        self.blocks: list = []
        self._edges = None
        self._order = None
        self._starts = None

    def add(self, frame: pd.DataFrame) -> int:
        """Add a block of links; returns how many."""
        if frame is not None and len(frame):
            self.blocks.append(frame[list(COLUMNS)])
            self._edges = None
        return 0 if frame is None else len(frame)

    @property
    def edges(self) -> pd.DataFrame:
        """Every link, one row each."""
        if self._edges is None:
            self._edges = (pd.concat(self.blocks, ignore_index=True) if self.blocks else empty())
            a = self._edges["a"].to_numpy(dtype=np.int64)
            b = self._edges["b"].to_numpy(dtype=np.int64)
            ends = np.concatenate([a, b])
            rows = np.concatenate([np.arange(len(a)), np.arange(len(a))])
            order = np.argsort(ends, kind="stable")
            self._order = rows[order]
            self._starts = np.searchsorted(ends[order], np.arange(self.n + 1))
        return self._edges

    def groups(self) -> pd.DataFrame:
        """Each group (layer, strategy, or your runs) with its kind and number of links."""
        e = self.edges
        if e.empty:
            return pd.DataFrame(columns=["group", "kind", "links"])
        g = e.groupby(["group", "kind"], sort=False).size().reset_index(name="links")
        return g

    def touching(self, gene: int, groups=None) -> pd.DataFrame:
        """The links of one gene, with `other` = the gene at the far end."""
        e = self.edges
        if e.empty or not 0 <= int(gene) < self.n:
            return e.assign(other=pd.Series(dtype="int32"))
        rows = self._order[self._starts[gene]:self._starts[gene + 1]]
        t = e.iloc[rows]
        if groups is not None:
            t = t[t["group"].isin(list(groups))]
        other = np.where(t["a"].to_numpy() == gene, t["b"].to_numpy(), t["a"].to_numpy())
        return t.assign(other=other.astype(np.int32))

    def counts(self, gene: int) -> dict:
        """group -> how many links this gene has in it."""
        t = self.touching(gene)
        return t.groupby("group").size().to_dict() if len(t) else {}


def star(index: EdgeIndex, centre: int, groups=None, depth: int = 1, per_group: int = 8,
         max_nodes: int = 160) -> tuple:
    """The neighbourhood the star map draws: (nodes, edges).

    `nodes` is a frame of `gene`, `hop` (0 centre, 1, 2) and `parent`; `edges` the links used, from
    enabled `groups` only. Around each expanded gene every group keeps its `per_group` strongest
    links -- so a dense layer (shared compartment joins hundreds of genes) cannot crowd out a sparse
    one (crosslinks) -- and the second hop keeps a third as many. At most `max_nodes` genes in all.
    """
    centre = int(centre)
    nodes = {centre: (0, -1)}
    used = []
    frontier = [centre]
    for hop in (1, 2)[:max(1, min(2, int(depth)))]:
        keep = per_group if hop == 1 else max(1, per_group // 3)
        nxt = []
        for g in frontier:
            t = index.touching(g, groups)
            if t.empty:
                continue
            if hop == 2:                    # the second hop reaches outward, not back inward
                t = t[~t["other"].isin([n for n, (h, _) in nodes.items() if h < 2])]
            t = (t.sort_values("strength", ascending=False, kind="stable")
                 .groupby("group", sort=False).head(keep))
            for row in t.itertuples(index=False):
                o = int(row.other)
                if o not in nodes:
                    if len(nodes) >= max_nodes:
                        continue
                    nodes[o] = (hop, g)
                    nxt.append(o)
                used.append(row)
        frontier = nxt
    nf = pd.DataFrame([(g, h, p) for g, (h, p) in nodes.items()], columns=["gene", "hop", "parent"])
    ef = pd.DataFrame(used) if used else empty().assign(other=pd.Series(dtype="int32"))
    if len(ef):
        ef = ef[ef["a"].isin(nf["gene"]) & ef["b"].isin(nf["gene"])]
        ef = ef.drop_duplicates(["a", "b", "run", "how", "group"]).reset_index(drop=True)
    return nf, ef
