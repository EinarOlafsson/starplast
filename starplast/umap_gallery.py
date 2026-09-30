#!/usr/bin/env python3
"""Pregenerated maps, and how well every label maps onto each of them.

The central map is one answer to "which measurements should place a gene": all of them, balanced by
block. It is not the only reasonable answer, and a label that fails to separate on it may separate
cleanly on a map built from one kind of evidence -- stage transcription, knockout fitness, protein
abundance. So a small gallery of maps is built ahead of time for each organism and shipped:

* **all measurements** at three `n_neighbors` settings (10, 25, 60), because the neighbourhood size
  decides whether a map shows many small islands or a few continents, and which one a label prefers
  is itself informative;
* **all measurements except localization**, so a localization label can be scored on a map that
  never saw a localization measurement;
* **one map per evidence family** -- the slot catalogue's axis (`strategies.Context.family_of`):
  transcription, translation, protein abundance, fitness, sequence, and so on -- where the family has
  enough columns and enough measured genes to make a map worth drawing.

Every map is built from measurements only. No categorical label column is an input (labels are what
is scored), and the recipe of each map is stored beside it. Each map is clustered with HDBSCAN through
`clustering.cluster`, as the Clusters tab does, at min_cluster_size 25, min_samples 5 and "leaf"
selection (see `CLUSTERING` for why not the tab's opening 25/25 "eom").

**The scores.** For a label (a categorical column, e.g. `compartment`) and a map:

* *Categories map to clusters* -- for each category of the label, the cluster that best recovers it
  is found and its F1 taken (precision: share of the cluster's labelled genes in the category;
  recall: share of the category's genes in that cluster). The label's score is the mean of these F1
  values weighted by category size, as `search.score_recovery` does. One difference, stated: genes
  HDBSCAN leaves unclustered still count in a category's size, so a map that calls most genes noise
  cannot score high recall on the few it clustered. Mean precision and mean recall are reported
  beside it.
* *Best single category* -- the one category that maps best, scored on the **F1 of the Wilson 95%
  lower bounds** of its precision and recall. A 2-gene category that sits alone in one cluster has
  precision 2/2 and recall 2/2, but lower bounds of 0.34 each, so it scores 0.34 rather than 1.0; a
  200-gene category at 0.9/0.9 keeps about 0.86. Size therefore counts through evidence, not through
  an arbitrary weight. The category chosen is the one whose bound beats its OWN shuffled-label
  chance by the most: on a map with one giant cluster, the majority value of a yes/no label has a
  high bound by sheer size and exactly as high a bound when shuffled, and it must not win for that.
* *Chance* -- both scores recomputed with the label's values shuffled over the same genes (five
  shuffles, fixed seeds) and averaged. "Skill" is (score - chance) / (1 - chance): 0 is what the map's
  cluster sizes give any label, 1 is perfect.
* *Circular* -- whether the map was built from a column in the label's leakage closure
  (`search.excluded_for`): a map built from hyperLOPIT probabilities recovering `compartment` is
  expected, not discovered.

The same scoring runs live for a map built in the application (`score_map`).
"""
from __future__ import annotations

import hashlib
import json
import os
import re

import numpy as np
import pandas as pd

from . import organisms

#: The files this module ships, in the package data directory.
COORDS_FILE = "umap_gallery.npz"
MANIFEST_FILE = "umap_gallery.json"
SCORES_FILE = "umap_gallery_scores.tsv"

#: The "all measurements" maps are built at these neighbourhood sizes. 25 is the application's own.
ALL_NEIGHBORS = (10, 25, 60)
#: The family left out of the "all but" map, so localization labels have a map that never saw them.
LEFT_OUT_FAMILY = "localization"
#: A family map needs this many source columns ...
MIN_FAMILY_COLUMNS = 4
#: ... and this many genes measured in at least `MIN_GENE_COVERAGE` of those columns.
MIN_FAMILY_GENES = 400
MIN_GENE_COVERAGE = 0.5
#: HDBSCAN as the Clusters tab runs it (`clustering.cluster`), with its "leaf" selection. The tab's
#: opening setting (min_cluster_size = min_samples = 25, "eom") was measured first and rejected for
#: the gallery: on these maps it returns a median of 4 clusters holding ~100% of genes (the Tg
#: transcription, translation, fitness and relation maps: 2 each), so the organisms' headline labels
#: scored at chance (median skill 0.01 for Tg compartment and Pf lopit_pf_location). Leaf selection
#: at 25/5 returns a median of 45-62 clusters, leaves 40-52% of genes unclustered, and raises the
#: median skill to 0.04-0.05. `build_umap_gallery.compare_clusterings` measures this in the notebook.
MIN_CLUSTER_SIZE = 25
CLUSTERING = {"min_cluster_size": MIN_CLUSTER_SIZE, "min_samples": 5,
              "cluster_selection_method": "leaf"}
#: Categories smaller than this are not scored at all: a one-gene category has no recall to speak of.
MIN_CATEGORY = 2
#: Label shuffles averaged for the chance level.
PERMUTATIONS = 5
#: Two-sided 95% normal quantile for the Wilson interval.
Z95 = 1.959963984540054
NOISE = -1

#: column -> (plain header, tooltip). The panel and the docs read these, so the two agree.
SCORE_COLUMNS = {
    "map": ("Map", "The pregenerated map (or the map on screen) the label was scored on."),
    "label": ("Label", "The categorical column whose categories are compared with the map's "
                       "clusters."),
    "category_f1": ("Categories → clusters",
                    "For each category, the F1 of the cluster that best recovers it; averaged over "
                    "categories weighted by their size. 1 means every category is its own cluster."),
    "category_skill": ("Skill",
                       "Categories → clusters, rescaled against the same score with the label "
                       "shuffled over the same genes: 0 is chance, 1 is perfect."),
    "category_precision": ("Precision",
                           "Size-weighted mean share of each best cluster's labelled genes that "
                           "belong to the category."),
    "category_recall": ("Recall",
                        "Size-weighted mean share of each category's genes that fall in its best "
                        "cluster. Unclustered genes count as missed."),
    "best_category": ("Best category", "The single category that maps onto one cluster furthest "
                                       "above chance, judged on the size-aware score."),
    "best_n": ("n", "Genes in the best category on this map. Small categories are scored on "
                    "lower bounds, so they need to map very cleanly to lead."),
    "best_f1_lower": ("Best category (size-aware)",
                      "F1 of the Wilson 95% lower bounds of the best category's precision and "
                      "recall. A 2-gene category alone in a cluster scores 0.34, not 1."),
    "best_skill": ("Best skill",
                   "The size-aware best-category score, rescaled against shuffled labels: 0 is "
                   "chance, 1 is perfect."),
    "best_precision": ("Best precision", "Raw precision of the best category's cluster."),
    "best_recall": ("Best recall", "Raw recall of the best category's cluster."),
    "coverage": ("Coverage", "Share of the label's labelled genes that are placed on this map. A "
                             "family map covers only the genes that family measured."),
    "genes_scored": ("Genes", "Labelled genes placed on the map: the genes the scores are over."),
    "categories": ("Categories", "Categories scored (at least 2 genes on the map)."),
    "circular": ("Circular", "Yes when the map was built from a column that restates this label "
                             "(its leakage closure): recovering it is then expected, not found."),
}


# --------------------------------------------------------------------------- arithmetic
def wilson_lower(k, n, z: float = Z95):
    """The Wilson score interval's lower bound for k successes in n trials (0 where n is 0)."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(n > 0, k / n, 0.0)
        z2 = z * z
        denom = 1.0 + z2 / n
        centre = (p + z2 / (2 * n)) / denom
        margin = z * np.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / denom
        out = np.where(n > 0, centre - margin, 0.0)
    return np.clip(out, 0.0, 1.0)


def _f1(p, r):
    p = np.asarray(p, dtype=float)
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p + r > 0, 2 * p * r / (p + r), 0.0)


def _codes(values) -> tuple:
    """Categories (as strings, absence removed) and an integer code per gene (-1 = unlabelled)."""
    from .search import ABSENCE_LABELS
    s = pd.Series(values, dtype=object)
    txt = s.astype(str).str.strip()
    ok = s.notna().to_numpy() & ~txt.str.lower().isin(ABSENCE_LABELS).to_numpy()
    cats = sorted(pd.unique(txt[ok]))
    lookup = {c: i for i, c in enumerate(cats)}
    code = np.full(len(s), -1, dtype=int)
    code[ok] = [lookup[t] for t in txt[ok]]
    return cats, code


def _best(t_code: np.ndarray, c_code: np.ndarray, n_cat: int, min_category: int) -> dict:
    """Per category: the best cluster by raw F1 and the best by lower-bound F1. Vectorised."""
    clusters = np.unique(c_code[c_code != NOISE])
    size = np.bincount(t_code, minlength=n_cat).astype(float)          # noise genes included
    if len(clusters) == 0:
        M = np.zeros((n_cat, 1))
        n_c = np.zeros(1)
        cl_ids = np.array([NOISE])
    else:
        idx = np.searchsorted(clusters, c_code)
        inc = c_code != NOISE
        M = np.zeros((n_cat, len(clusters)))
        np.add.at(M, (t_code[inc], idx[inc]), 1.0)
        n_c = M.sum(axis=0)                                            # labelled genes per cluster
        cl_ids = clusters
    with np.errstate(divide="ignore", invalid="ignore"):
        prec = np.where(n_c > 0, M / n_c, 0.0)
        rec = np.where(size[:, None] > 0, M / size[:, None], 0.0)
    f1 = _f1(prec, rec)
    lo = _f1(wilson_lower(M, np.broadcast_to(n_c, M.shape)),
             wilson_lower(M, np.broadcast_to(size[:, None], M.shape)))
    j = f1.argmax(axis=1)
    jl = lo.argmax(axis=1)
    rows = np.arange(n_cat)
    keep = size >= min_category
    return {"size": size, "keep": keep, "f1": f1[rows, j], "precision": prec[rows, j],
            "recall": rec[rows, j], "cluster": cl_ids[j], "lower": lo[rows, jl],
            "lower_precision": prec[rows, jl], "lower_recall": rec[rows, jl],
            "lower_cluster": cl_ids[jl]}


def _weighted(b: dict) -> dict:
    """The size-weighted category means of one `_best` result."""
    keep = b["keep"]
    if not keep.any():
        return {"category_f1": np.nan, "category_precision": np.nan, "category_recall": np.nan}
    w = b["size"][keep] / b["size"][keep].sum()
    return {"category_f1": float((b["f1"][keep] * w).sum()),
            "category_precision": float((b["precision"][keep] * w).sum()),
            "category_recall": float((b["recall"][keep] * w).sum())}


def skill(score, chance):
    """(score - chance) / (1 - chance): 0 at chance, 1 at perfect, negative below chance.

    Works elementwise on arrays; NaN where either is missing or chance is already perfect.
    """
    score = np.asarray(score, dtype=float)
    chance = np.asarray(chance, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(np.isfinite(score) & np.isfinite(chance) & (chance < 1.0),
                       (score - chance) / (1.0 - chance), np.nan)
    return float(out) if out.ndim == 0 else out


def _scored(truth, clusters, min_category: int, permutations: int, seed: int):
    """Categories, the observed `_best`, and the per-category / headline chance levels."""
    cats, code = _codes(truth)
    clusters = np.asarray(clusters)
    ok = code >= 0
    if ok.sum() < 2:
        return cats, None, None, None, int(ok.sum())
    t, c = code[ok], clusters[ok]
    b = _best(t, c, len(cats), min_category)
    rng = np.random.default_rng(seed)
    lows, f1s = [], []
    for _ in range(int(permutations)):
        bp = _best(rng.permutation(t), c, len(cats), min_category)
        lows.append(bp["lower"])
        f1s.append(_weighted(bp)["category_f1"])
    chance_lower = np.mean(lows, axis=0) if lows else np.full(len(cats), np.nan)
    chance_f1 = float(np.mean(f1s)) if f1s else np.nan
    return cats, b, chance_lower, chance_f1, int(ok.sum())


def per_category(truth, clusters, min_category: int = MIN_CATEGORY,
                 permutations: int = PERMUTATIONS, seed: int = 0) -> pd.DataFrame:
    """Every category of a label against one clustering: best cluster, P, R, F1, bounded F1, skill.

    `truth` and `clusters` are aligned per gene (the genes placed on one map); -1 in `clusters` is
    HDBSCAN noise. Unlabelled genes are ignored; unclustered labelled genes count as missed. Sorted
    by `skill`, the size-aware score against its own shuffled-label chance level -- the order the
    best single category is chosen in.
    """
    cats, b, chance_lower, _, _ = _scored(truth, clusters, min_category, permutations, seed)
    cols = ["category", "n", "cluster", "precision", "recall", "f1", "f1_lower", "chance",
            "skill", "scored"]
    if b is None:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame({"category": cats, "n": b["size"].astype(int),
                         "cluster": b["lower_cluster"], "precision": b["lower_precision"],
                         "recall": b["lower_recall"], "f1": _f1(b["lower_precision"],
                                                                b["lower_recall"]),
                         "f1_lower": b["lower"], "chance": chance_lower,
                         "skill": skill(b["lower"], chance_lower), "scored": b["keep"]},
                        columns=cols).sort_values(["scored", "skill"], ascending=False,
                                                  ignore_index=True)


def label_scores(truth, clusters, min_category: int = MIN_CATEGORY,
                 permutations: int = PERMUTATIONS, seed: int = 0) -> dict:
    """The scores of one label on one clustered map, with the shuffled-label chance level.

    Returns category_f1 / category_precision / category_recall (size-weighted over categories),
    the chance level and skill of category_f1, and the best single category: the one whose
    size-aware score (`best_f1_lower`, the F1 of the Wilson lower bounds) beats its OWN chance
    level by the most. Chosen on skill rather than on the raw bound, because on a map with one
    giant cluster the majority value of a yes/no label scores a high bound by sheer size, and
    exactly as high with the labels shuffled.
    """
    cats, b, chance_lower, chance_f1, n = _scored(truth, clusters, min_category, permutations,
                                                  seed)
    empty = {"genes_scored": n, "categories": 0, "category_f1": np.nan,
             "category_precision": np.nan, "category_recall": np.nan, "category_f1_chance": np.nan,
             "category_skill": np.nan, "best_category": "", "best_n": 0, "best_f1_lower": np.nan,
             "best_precision": np.nan, "best_recall": np.nan, "best_cluster": NOISE,
             "best_f1_lower_chance": np.nan, "best_skill": np.nan}
    if b is None or not b["keep"].any():
        return empty
    w = _weighted(b)
    sk = skill(b["lower"], chance_lower)
    order = np.where(b["keep"], np.nan_to_num(sk, nan=-np.inf), -np.inf)
    k = int(np.argmax(order))
    return {"genes_scored": n, "categories": int(b["keep"].sum()), **w,
            "category_f1_chance": chance_f1, "category_skill": skill(w["category_f1"], chance_f1),
            "best_category": str(cats[k]), "best_n": int(b["size"][k]),
            "best_f1_lower": float(b["lower"][k]), "best_precision": float(b["lower_precision"][k]),
            "best_recall": float(b["lower_recall"][k]), "best_cluster": int(b["lower_cluster"][k]),
            "best_f1_lower_chance": float(chance_lower[k]), "best_skill": float(sk[k])}


# --------------------------------------------------------------------------- recipes
def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")


def recipes(ctx) -> list:
    """The gallery for one organism's context: a list of dicts (id, title, blocks, n_neighbors, ...).

    Blocks are the display recipe's slot blocks (`embedding.default_spec`), grouped into families by
    the slot catalogue's axis. Only measurements are inputs: the slot catalogue's feature role
    already excludes label columns, and any categorical column that slipped through is refused.
    """
    from .embedding import default_spec, columns_for, EmbeddingSpec
    nodes = ctx.nodes
    spec = default_spec(nodes)
    cols = columns_for(nodes, EmbeddingSpec(blocks=spec.blocks))
    labels = set(ctx.categorical_columns())
    cols = {b: [c for c in cc if c not in labels] for b, cc in cols.items()}
    cols = {b: cc for b, cc in cols.items() if cc}
    family = {b: ctx.family_of(b) for b in cols}
    measured = {c: nodes[c].notna().to_numpy() for cc in cols.values() for c in cc}

    def entry(map_id, title, blocks, description, n_neighbors=25, fam=""):
        source = sorted({c for b in blocks for c in cols[b]})
        return {"id": map_id, "title": title, "family": fam, "description": description,
                "blocks": list(blocks), "columns": source, "n_neighbors": int(n_neighbors)}

    everything = list(cols)
    out = [entry(f"all_nn{k}", f"All measurements · n_neighbors {k}", everything,
                 f"Every measurement block ({len(everything)}), each scaled to equal total "
                 f"variance; neighbourhood size {k}"
                 + (" (the application's default)." if k == 25 else "."), k, "all")
           for k in ALL_NEIGHBORS]
    blind = [b for b in everything if family[b] != LEFT_OUT_FAMILY]
    if len(blind) < len(everything):
        out.append(entry("all_but_localization", "All but localization", blind,
                         "Every measurement block except the localization family, so a "
                         "localization label can be scored on a map that never saw one.", 25,
                         "all"))
    fams = {}
    for b in everything:
        fams.setdefault(family[b], []).append(b)
    for fam, blocks in sorted(fams.items(), key=lambda kv: kv[0].lower()):
        e = entry(f"family_{_slug(fam)}", f"{fam[:1].upper()}{fam[1:]} only", blocks,
                  f"The {fam} family alone: {len(blocks)} block(s).", 25, fam)
        n_cols = len(e["columns"])
        if n_cols < MIN_FAMILY_COLUMNS:
            continue
        cov = np.mean([measured[c] for c in e["columns"]], axis=0)
        if int((cov >= MIN_GENE_COVERAGE).sum()) < MIN_FAMILY_GENES:
            continue
        e["min_coverage"] = MIN_GENE_COVERAGE
        out.append(e)
    return out


def genes_for(ctx, recipe: dict) -> np.ndarray:
    """Positions of the genes a recipe places: all of them, or those measured often enough."""
    cov = recipe.get("min_coverage")
    if not cov:
        return np.arange(ctx.n)
    m = np.mean([ctx.nodes[c].notna().to_numpy() for c in recipe["columns"]], axis=0)
    return np.flatnonzero(m >= float(cov))


def build_map(ctx, recipe: dict, log=print) -> dict:
    """Embed and cluster one recipe. Returns rows, xyz (float32), clusters (int16) and provenance."""
    from .embedding import EmbeddingSpec, embed
    from .clustering import cluster
    rows = genes_for(ctx, recipe)
    # Label columns are removed from the table the map is built from, not only from the recipe: a
    # slot block resolves its own columns, and a 0/1 label such as `is_exported` is numeric enough
    # to be one of them. A map built from a label cannot then be scored against it.
    labels = [c for c in ctx.categorical_columns() if c != "gene_id"]
    sub = ctx.nodes.iloc[rows].drop(columns=labels).reset_index(drop=True)
    spec = EmbeddingSpec(name=f"gallery-{recipe['id']}", blocks=tuple(recipe["blocks"]),
                         n_neighbors=int(recipe["n_neighbors"]), n_components=3)
    coords, names, kept, meta = embed(sub, spec, log=log, strict=True, return_metadata=True)
    rows = rows[np.asarray(kept, dtype=bool)]
    leaked = sorted(set(names) & set(labels))
    if leaked:
        raise ValueError(f"{recipe['id']}: label column(s) reached the matrix: {leaked}")
    xyz = np.asarray(coords, dtype=np.float32)
    lab = np.asarray(cluster(xyz, "hdbscan", **CLUSTERING), dtype=np.int16)
    n_cl = int(len(set(lab[lab >= 0].tolist())))
    return {"rows": rows.astype(np.int32), "xyz": xyz, "clusters": lab,
            "info": {**recipe, "recipe": spec.to_dict(), "executed_method": meta["executed_method"],
                     "backend": meta["backend"], "n_features": len(names),
                     "n_genes": int(len(rows)), "n_clusters": n_cl,
                     "noise_fraction": float((lab < 0).mean()),
                     "clustering": {"algorithm": "hdbscan", **CLUSTERING},
                     "coordinates_sha256": hashlib.sha256(xyz.tobytes()).hexdigest(),
                     "versions": meta.get("versions", {})}}


# --------------------------------------------------------------------------- scoring many
def score_map(ctx, rows, clusters, labels=None, map_id: str = "on screen", columns=None,
              permutations: int = PERMUTATIONS) -> pd.DataFrame:
    """Every label's scores on one clustered map (the gallery's, or one built in the app).

    `rows` are the positions in `ctx.nodes` the map places, `clusters` one label per placed gene.
    `columns`, when known, are the map's source columns, used for the circularity flag; None leaves
    it unknown.
    """
    rows = np.asarray(rows, dtype=int)
    clusters = np.asarray(clusters)
    out = []
    for i, label in enumerate(labels if labels is not None else ctx.categorical_columns()):
        truth = ctx.truth(label)
        labelled = int(truth.notna().sum())
        sc = label_scores(truth.iloc[rows].to_numpy(), clusters, permutations=permutations,
                          seed=i)
        if columns is None:
            circ, why = None, ""
        else:
            hit = sorted(set(columns) & ctx.banned(label))
            circ, why = bool(hit), ", ".join(hit[:6]) + (" …" if len(hit) > 6 else "")
        out.append({"organism": ctx.organism, "label": label, "map": map_id,
                    "coverage": sc["genes_scored"] / labelled if labelled else np.nan,
                    **sc, "circular": circ, "circular_columns": why})
    return pd.DataFrame(out)


def score_gallery(ctx, maps: dict, labels=None, log=print) -> pd.DataFrame:
    """The label x map table for one organism. `maps` is map id -> {rows, clusters, info}."""
    frames = []
    for map_id, m in maps.items():
        log(f"scoring {ctx.organism} {map_id}")
        frames.append(score_map(ctx, m["rows"], m["clusters"], labels=labels, map_id=map_id,
                                columns=m["info"].get("columns")))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# --------------------------------------------------------------------------- files
def _data_dir() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def write(built: dict, scores: pd.DataFrame, gene_ids: dict, directory: str | None = None,
          extra: dict | None = None) -> dict:
    """Write the gallery: coordinates (npz), manifest (json), scores (tsv). Returns the paths.

    `built` is organism -> map id -> build_map() result; `gene_ids` organism -> the node table's
    gene ids, stored so a loader can refuse a gallery built over a different table.
    """
    directory = directory or _data_dir()
    arrays, manifest = {}, {"version": 1, "organisms": {}, **(extra or {})}
    for code, maps in built.items():
        ids = np.asarray(gene_ids[code]).astype(str)
        arrays[f"{code}__gene_ids"] = ids
        manifest["organisms"][code] = {
            "n_genes": int(len(ids)),
            "gene_ids_sha256": hashlib.sha256("\n".join(ids).encode()).hexdigest(),
            "maps": []}
        for map_id, m in maps.items():
            arrays[f"{code}__{map_id}__rows"] = np.asarray(m["rows"], dtype=np.int32)
            arrays[f"{code}__{map_id}__xyz"] = np.asarray(m["xyz"], dtype=np.float32)
            arrays[f"{code}__{map_id}__clusters"] = np.asarray(m["clusters"], dtype=np.int16)
            manifest["organisms"][code]["maps"].append(m["info"])
    paths = {k: os.path.join(directory, f) for k, f in
             (("coords", COORDS_FILE), ("manifest", MANIFEST_FILE), ("scores", SCORES_FILE))}
    np.savez_compressed(paths["coords"], **arrays)
    with open(paths["manifest"], "w", encoding="utf8") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=False, default=_jsonable)
        fh.write("\n")
    scores.to_csv(paths["scores"], sep="\t", index=False, float_format="%.4f")
    return paths


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray, tuple, set)):
        return list(o)
    raise TypeError(type(o).__name__)


class Gallery:
    """The shipped gallery, loaded lazily: manifest, per-map coordinates and clusters, scores."""

    def __init__(self, directory: str | None = None):
        """Point at a directory holding the three gallery files (the package data by default)."""
        self.directory = directory or _data_dir()
        self._manifest = None
        self._npz = None
        self._scores = None

    def available(self) -> bool:
        """Whether all three files are present."""
        return all(os.path.exists(os.path.join(self.directory, f))
                   for f in (COORDS_FILE, MANIFEST_FILE, SCORES_FILE))

    @property
    def manifest(self) -> dict:
        """The recipes and build facts of every map."""
        if self._manifest is None:
            with open(os.path.join(self.directory, MANIFEST_FILE), encoding="utf8") as fh:
                self._manifest = json.load(fh)
        return self._manifest

    def organisms(self) -> list:
        """Organism codes the gallery covers."""
        return list(self.manifest.get("organisms", {}))

    def maps(self, organism: str) -> list:
        """The map records (recipes, sizes, cluster counts) for one organism, in gallery order."""
        return list(self.manifest.get("organisms", {}).get(organism, {}).get("maps", []))

    def info(self, organism: str, map_id: str) -> dict:
        """One map's record."""
        for m in self.maps(organism):
            if m["id"] == map_id:
                return m
        raise KeyError(f"no gallery map {map_id!r} for {organism}")

    def _arrays(self):
        if self._npz is None:
            self._npz = np.load(os.path.join(self.directory, COORDS_FILE), allow_pickle=False)
        return self._npz

    def gene_ids(self, organism: str) -> np.ndarray:
        """The node table's gene order the gallery was built over."""
        return self._arrays()[f"{organism}__gene_ids"].astype(str)

    def load(self, organism: str, map_id: str) -> dict:
        """rows (positions in the build's node table), xyz (float32) and clusters of one map."""
        z = self._arrays()
        key = f"{organism}__{map_id}"
        return {"rows": z[key + "__rows"].astype(int), "xyz": z[key + "__xyz"],
                "clusters": z[key + "__clusters"].astype(int),
                "info": self.info(organism, map_id)}

    def aligned(self, organism: str, map_id: str, gene_ids) -> dict:
        """One map re-indexed onto another table's gene order (genes it lacks are dropped).

        The window's table normally IS the build's table; this is what keeps a gallery usable when
        it is not (a rebuilt cache with genes added or reordered) instead of drawing genes in the
        wrong places.
        """
        m = self.load(organism, map_id)
        built = self.gene_ids(organism)
        target = np.asarray(gene_ids).astype(str)
        if len(built) == len(target) and (built == target).all():
            return m
        where = {g: i for i, g in enumerate(target)}
        pos = np.array([where.get(g, -1) for g in built[m["rows"]]])
        ok = pos >= 0
        return {**m, "rows": pos[ok], "xyz": m["xyz"][ok], "clusters": m["clusters"][ok]}

    def scores(self, organism: str | None = None) -> pd.DataFrame:
        """The label x map table (for one organism, or all)."""
        if self._scores is None:
            self._scores = pd.read_csv(os.path.join(self.directory, SCORES_FILE), sep="\t",
                                       keep_default_na=False, na_values=[""])
        s = self._scores
        return s[s.organism == organism].reset_index(drop=True) if organism else s.copy()


_SHIPPED = None


def shipped() -> Gallery:
    """The packaged gallery, shared."""
    global _SHIPPED
    if _SHIPPED is None:
        _SHIPPED = Gallery()
    return _SHIPPED


def space_codes() -> list:
    """The organisms a gallery is built for: every parasite space whose table ships."""
    return organisms.codes(organisms.PARASITE, available=True)


def describe(scores: pd.DataFrame, label: str) -> str:
    """One line per map for a label, best first: the text form of the panel's table."""
    s = scores[scores.label == label].sort_values("category_skill", ascending=False)
    lines = []
    for r in s.itertuples():
        lines.append(f"{r.map:24s} categories→clusters {r.category_f1:.3f} (skill "
                     f"{r.category_skill:+.3f})  best {r.best_category!s} n={r.best_n} "
                     f"{r.best_f1_lower:.3f} (skill {r.best_skill:+.3f})"
                     + ("  [circular]" if str(r.circular) == "True" else ""))
    return "\n".join(lines)
