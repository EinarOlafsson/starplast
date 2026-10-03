"""Build the pregenerated map gallery and its label x map score table.

For every parasite space whose table ships (`starplast.organisms`), this builds the maps listed by
`starplast.umap_gallery.recipes`: all measurements at three n_neighbors settings, an "all but one
family" map per evidence family, one map per family, one per individual experiment block, and the
curated pairs and triples of families. Each is a 3D UMAP of measurements only (label columns are
removed from the table first), clustered with HDBSCAN. Every map but the three fixed reference maps
is SEARCHED over `umap_gallery.MAP_GRID` x `umap_gallery.GALLERY_CLUSTER_GRID` and the combination
with the best label-free structure score (`umap_gallery.structure`) is the one that ships. Every
categorical label is then scored on every map (`umap_gallery.score_gallery`).

Output, in `starplast/data/`:

* `umap_gallery.npz`     -- per map: placed gene positions (int32), coordinates (float32), clusters (int16)
* `umap_gallery.json`    -- per map: recipe, blocks, source columns, sizes, structure, settings searched
* `umap_gallery_scores.tsv` -- the label x map table

and the executed notebook `notebooks/umap_gallery_2026_09_30.ipynb` that records the run.

The embedding runs on the CPU (`STARPLAST_GPU=0`) so the shipped maps are the umap-learn maps a
user without CUDA would reproduce from the stored recipe. Maps are built in parallel processes and
**checkpointed one map at a time**, so a stopped run resumes instead of starting over. Run the whole
thing under a memory cap:

    systemd-run --user --scope -p MemoryMax=25G env STARPLAST_GPU=0 \\
        ~/anaconda3/envs/starplast/bin/python scripts/build_umap_gallery.py --workers 6

    python scripts/build_umap_gallery.py --no-notebook --organism <code>   # quick, no notebook
    python scripts/build_umap_gallery.py --no-checkpoint                   # rebuild everything

A user who wants more maps than ship -- every block rather than the ones with enough columns, say --
edits the thresholds in `starplast/umap_gallery.py` and runs this script with `--out` pointing
somewhere else; `umap_gallery.Gallery(directory)` loads any such gallery, and the panel takes one.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("STARPLAST_GPU", "0")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from starplast import umap_gallery as G  # noqa: E402

NOTEBOOK = os.path.join(ROOT, "notebooks", "umap_gallery_2026_09_30.ipynb")
#: One file per built map, so a stopped run resumes. Outside the package: these are not shipped.
CHECKPOINT = os.path.join(ROOT, "results", "umap_gallery_checkpoint")

_CTX: dict = {}


def context(code: str):
    """The shipped context of one organism, once per process."""
    if code not in _CTX:
        from starplast import strategies as S
        _CTX[code] = S.Context.shipped(code)
    return _CTX[code]


def _path(directory: str, code: str, map_id: str) -> str:
    return os.path.join(directory, f"{code}__{map_id}.npz")


def _save(directory: str, code: str, map_id: str, m: dict):
    """Checkpoint one built map: its arrays and its info, in one file."""
    os.makedirs(directory, exist_ok=True)
    tmp = _path(directory, code, map_id) + ".tmp.npz"
    np.savez_compressed(tmp, rows=m["rows"], xyz=m["xyz"], clusters=m["clusters"],
                        info=np.asarray(json.dumps(m["info"], default=G._jsonable)))
    os.replace(tmp, _path(directory, code, map_id))


def _load(directory: str, code: str, map_id: str):
    """One checkpointed map, or None when it is absent or unreadable."""
    path = _path(directory, code, map_id)
    if not os.path.exists(path):
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            return {"rows": z["rows"], "xyz": z["xyz"], "clusters": z["clusters"],
                    "info": json.loads(str(z["info"]))}
    except Exception:
        return None


def _one(job):
    code, recipe = job
    t0 = time.monotonic()
    m = G.build_map(context(code), recipe, log=lambda *a: None)
    m["info"]["seconds"] = round(time.monotonic() - t0, 1)
    return code, recipe["id"], m


def build(codes=None, workers: int = 6, log=print, checkpoint: str | None = CHECKPOINT) -> dict:
    """organism -> map id -> built map, in recipe order.

    Each finished map is written to `checkpoint` as it arrives and read back on the next run, so a
    build that is stopped -- or frozen by the memory guard -- resumes at the map it reached instead
    of repeating the hours before it. Pass `checkpoint=None` to rebuild everything.
    """
    codes = list(codes or G.space_codes())
    jobs = [(code, r) for code in codes for r in G.recipes(context(code))]
    done, todo = {}, []
    for code, r in jobs:
        was = _load(checkpoint, code, r["id"]) if checkpoint else None
        if was is not None and was["info"].get("columns") == r["columns"]:
            done[(code, r["id"])] = was
        else:
            todo.append((code, r))
    log(f"{len(jobs)} maps over {', '.join(codes)} on {workers} worker(s); "
        f"{len(done)} already built, {len(todo)} to build")
    if todo:
        with ProcessPoolExecutor(max_workers=max(1, int(workers))) as pool:
            for i, (code, map_id, m) in enumerate(pool.map(_one, todo), 1):
                log(f"  [{i}/{len(todo)}] {code} {map_id}: {m['info']['n_genes']:,} genes, "
                    f"{m['info']['n_features']} features, {m['info']['n_clusters']} clusters, "
                    f"{m['info']['noise_fraction']:.0%} unclustered, structure "
                    f"{m['info']['structure']:.3f} ({m['info']['seconds']} s)")
                done[(code, map_id)] = m
                if checkpoint:
                    _save(checkpoint, code, map_id, m)
    return {code: {r["id"]: done[(code, r["id"])] for c2, r in jobs if c2 == code}
            for code in codes}


def score(built: dict, log=print) -> pd.DataFrame:
    """The label x map table over every organism."""
    return pd.concat([G.score_gallery(context(code), maps, log=log)
                      for code, maps in built.items()], ignore_index=True)


def headline_label(code: str) -> str:
    """The organism's first categorical calibration target that the table carries."""
    from starplast import organisms
    cats = context(code).categorical_columns()
    return next((t for t in organisms.get(code).targets if t in cats), cats[0])


def compare_clusterings(built: dict) -> pd.DataFrame:
    """The Clusters tab's opening HDBSCAN setting against the gallery's, on every shipped map.

    Scored on each organism's headline label, so the choice of `CLUSTERING` is a measured one.
    """
    from starplast.clustering import cluster
    opening = {"min_cluster_size": 25, "min_samples": 25, "cluster_selection_method": "eom"}
    rows = []
    for code, maps in built.items():
        label = headline_label(code)
        truth = context(code).truth(label)
        for map_id, m in maps.items():
            t = truth.iloc[m["rows"]].to_numpy()
            for name, kw, lab in (("tab opening 25/25 eom", opening,
                                   cluster(m["xyz"], "hdbscan", **opening)),
                                  ("gallery", G.CLUSTERING, m["clusters"])):
                s = G.label_scores(t, lab)
                rows.append({"organism": code, "label": label, "map": map_id, "setting": name,
                             "clusters": int(len(set(lab[lab >= 0].tolist()))),
                             "unclustered": round(float((lab < 0).mean()), 3),
                             "category_f1": round(s["category_f1"], 3),
                             "category_skill": round(s["category_skill"], 3)})
    return pd.DataFrame(rows)


def maps_table(built: dict) -> pd.DataFrame:
    """One row per map: what it is built from, how it was settled, and how much structure it has."""
    rows = []
    for code, maps in built.items():
        for map_id, m in maps.items():
            i = m["info"]
            rows.append({"organism": code, "map": map_id, "group": i.get("group", ""),
                         "family": i["family"], "blocks": len(i["blocks"]),
                         "columns": len(i["columns"]), "features": i["n_features"],
                         "coverage": i.get("min_coverage", 0.0), "tuned": i.get("tuned", False),
                         "n_neighbors": i["n_neighbors"], "min_dist": i.get("min_dist"),
                         "min_cluster_size": (i.get("clustering") or {}).get("min_cluster_size"),
                         "min_samples": (i.get("clustering") or {}).get("min_samples"),
                         "genes": i["n_genes"], "clusters": i["n_clusters"],
                         "unclustered": round(i["noise_fraction"], 3),
                         "clustered": i.get("clustered"), "evenness": i.get("evenness"),
                         "separation": i.get("separation"), "structure": i.get("structure"),
                         "searched": len(i.get("settings_searched") or []),
                         "seconds": i.get("seconds")})
    return pd.DataFrame(rows).sort_values(["organism", "structure"], ascending=[True, False],
                                          ignore_index=True)


def search_table(built: dict) -> pd.DataFrame:
    """Every setting tried for every searched map, with the structure score it reached.

    This is the record of HOW each shipped map was chosen: the whole grid, the sample stage and the
    full stage, and which row won. It is also what shows the search was worth running -- if the best
    setting were always the same one, the grid would be a waste of hours.
    """
    rows = []
    for code, maps in built.items():
        for map_id, m in maps.items():
            i = m["info"]
            for r in (i.get("settings_searched") or []):
                rows.append({"organism": code, "map": map_id, **r,
                             "won": bool(r.get("stage") == "full"
                                         and r.get("n_neighbors") == i["n_neighbors"]
                                         and r.get("min_dist") == i.get("min_dist")
                                         and r.get("structure") == i.get("structure"))})
    return pd.DataFrame(rows)


def search_gain(built: dict) -> pd.DataFrame:
    """Per searched map: what the search chose, and what it was worth.

    The gain is measured at the SAMPLE stage, where every UMAP setting in the grid was scored on the
    same 1,500 genes: the best of them against the neighbourhood the gallery used to fix for every
    map (n_neighbors 25, min_dist 0.25). It is therefore the UMAP half of the gain; the clustering
    half is `cluster_spread`, the range the clustering grid covered on the winning embedding. The
    full stage cannot answer the UMAP question, because only the surviving settings are rebuilt at
    full size -- that is the point of successive halving.

    `spread_full` is the remaining disagreement between the survivors once they were rebuilt whole,
    which is how much the sample stage's ranking could still be wrong by.
    """
    rows = []
    for code, maps in built.items():
        for map_id, m in maps.items():
            i = m["info"]
            tried = i.get("settings_searched") or []
            sample = [r for r in tried if r.get("stage") == "sample"]
            full = [r for r in tried if r.get("stage") == "full"]
            if not sample:
                continue
            old = [r for r in sample if r["n_neighbors"] == 25 and r["min_dist"] == 0.25]
            best = max(r["structure"] for r in sample)
            rows.append({"organism": code, "map": map_id, "group": i.get("group", ""),
                         "shipped": i.get("structure"),
                         "sample_best": best,
                         "sample_at_old_default": old[0]["structure"] if old else None,
                         "spread_full": (round(max(r["structure"] for r in full)
                                               - min(r["structure"] for r in full), 4)
                                         if len(full) > 1 else 0.0),
                         "n_neighbors": i["n_neighbors"], "min_dist": i.get("min_dist"),
                         "min_cluster_size": (i.get("clustering") or {}).get("min_cluster_size"),
                         "min_samples": (i.get("clustering") or {}).get("min_samples")})
    t = pd.DataFrame(rows)
    if len(t):
        t["gain_on_sample"] = (t["sample_best"] - t["sample_at_old_default"]).round(4)
        t["chose_old_default"] = (t["n_neighbors"] == 25) & (t["min_dist"] == 0.25)
    return t


def write(built: dict, scores: pd.DataFrame, directory: str | None = None) -> dict:
    """Write the three shipped files and report their sizes."""
    ids = {code: context(code).gene_ids for code in built}
    paths = G.write(built, scores, ids, directory=directory, extra={
        "built": _dt.date.today().isoformat(),
        "script": "scripts/build_umap_gallery.py",
        "notebook": os.path.relpath(NOTEBOOK, ROOT),
        "groups": list(G.GROUPS),
        "search": {"umap_grid": {k: list(v) for k, v in G.MAP_GRID.items()},
                   "clustering_grid": {k: list(v) for k, v in G.GALLERY_CLUSTER_GRID.items()},
                   "sample": G.SEARCH_SAMPLE, "keep": G.SEARCH_KEEP,
                   "ranked_on": "umap_gallery.structure: sqrt(search.map_quality score x mean "
                                "silhouette of the clustered points, rescaled to 0..1). No label "
                                "is read, so the search cannot tune a map onto what it is scored on.",
                   "fixed": [f"all_nn{k}" for k in G.ALL_NEIGHBORS]},
        "scoring": {"min_category": G.MIN_CATEGORY, "permutations": G.PERMUTATIONS,
                    "bound": "Wilson 95% lower bounds of precision and recall, combined as F1",
                    "noise": "HDBSCAN-unclustered genes count in a category's size (missed)"}})
    for k, p in paths.items():
        print(f"{k:9s} {os.path.relpath(p, ROOT)}  {os.path.getsize(p) / 1e6:.2f} MB")
    return paths


def figure(built: dict, code: str, which=None, cols: int = 6, limit: int = 24):
    """Maps of one organism, first two components, colored by their clusters.

    `which` is a list of map ids; by default the `limit` best-structured maps, since a gallery of
    dozens does not fit on a page and the ones worth looking at first are the ones that found
    structure. Titled with the structure score so the figure and the table agree.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    maps = built[code]
    if which is None:
        which = [r["map"] for _, r in
                 maps_table({code: maps}).sort_values("structure", ascending=False).iterrows()]
        which = which[:limit]
    which = [m for m in which if m in maps]
    rows = max(1, (len(which) + cols - 1) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(2.4 * cols, 2.5 * rows), squeeze=False)
    for ax in np.ravel(axes):
        ax.axis("off")
    for ax, map_id in zip(np.ravel(axes), which):
        m = maps[map_id]
        lab, xyz = m["clusters"], m["xyz"]
        ax.scatter(xyz[lab < 0, 0], xyz[lab < 0, 1], s=0.6, c="#bbbbbb", linewidths=0)
        ax.scatter(xyz[lab >= 0, 0], xyz[lab >= 0, 1], s=0.6, c=lab[lab >= 0], cmap="tab20",
                   linewidths=0)
        ax.set_title(f"{map_id}\n{m['info']['n_genes']:,} genes, {m['info']['n_clusters']} cl, "
                     f"structure {m['info'].get('structure', float('nan')):.2f}", fontsize=6)
    fig.tight_layout()
    return fig


def notebook(codes=None, workers: int = 6, path: str = NOTEBOOK) -> str:
    """Run the whole build inside an executed notebook and write it."""
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from notebook_runner import ExecutedNotebook
    codes = list(codes or G.space_codes())
    nb = ExecutedNotebook("Pregenerated map gallery and label scores", {"__name__": "__notebook__"})
    nb.md("The application's central map is one choice of inputs. This notebook builds a gallery of "
          "alternatives for each organism and scores how well every categorical label maps onto each. "
          "The code is `scripts/build_umap_gallery.py` and `starplast/umap_gallery.py`; this notebook "
          "calls them, so what is shown is what ships.",
          "Six groups of maps: all measurements at three neighbourhood sizes, an 'all but one "
          "family' map per family, one map per evidence family, one per individual experiment block, "
          "and the curated pairs and triples of families.",
          "Every map is built from measurements only: label columns are removed from the table "
          "before embedding, and the build refuses a map whose matrix contains one.",
          "**Maximising structure.** Every map except the three fixed reference maps is searched over "
          f"`MAP_GRID` ({G.MAP_GRID}) x `GALLERY_CLUSTER_GRID` ({G.GALLERY_CLUSTER_GRID}) and the "
          "combination with the best LABEL-FREE structure score ships. The score is "
          "`umap_gallery.structure`: the geometric mean of `search.map_quality`'s score (share of "
          "genes clustered x evenness of the cluster sizes, zero if the clustering is degenerate) "
          "and the mean silhouette of the clustered points rescaled to 0..1. Either alone is "
          "gameable -- an even partition of an unseparated cloud scores well on the first, two crisp "
          "clusters holding 5% of the proteome on the second -- so a map has to be good at both. No "
          "label is read by the search, which is what lets the labels be scored honestly afterwards.")
    nb.code("import warnings; warnings.filterwarnings('ignore')",
            f"import sys; sys.path.insert(0, {ROOT!r}); "
            f"sys.path.insert(0, {os.path.join(ROOT, 'scripts')!r})",
            "import build_umap_gallery as B",
            "from starplast import umap_gallery as G",
            "import pandas as pd",
            "_ = pd.set_option('display.width', 220)")
    nb.md("## 1. The recipes",
          f"A family map needs at least {G.MIN_FAMILY_COLUMNS} source columns, a single-experiment "
          f"map at least {G.MIN_BLOCK_COLUMNS}, and every subset map at least {G.MIN_SET_GENES} "
          f"genes. A gene is placed when it is measured in at least `min_coverage` of the recipe's "
          f"columns; the thresholds {G.COVERAGE_STEPS} are tried in order and the strictest that "
          "still places enough genes wins, so a dense family is built over well-measured genes and a "
          "sparse one is still buildable. The 'all' maps place every gene, as the application does.")
    nb.code("rows = []",
            f"for code in {codes!r}:",
            "    ctx = B.context(code)",
            "    for r in G.recipes(ctx):",
            "        rows.append({'organism': code, 'group': r['group'], 'map': r['id'],",
            "                     'blocks': len(r['blocks']), 'columns': len(r['columns']),",
            "                     'coverage': r.get('min_coverage', 0.0),",
            "                     'genes': len(G.genes_for(ctx, r)), 'tuned': r['tune']})",
            "recipes = pd.DataFrame(rows)",
            "print(recipes.groupby(['organism', 'group']).size().to_string())",
            "print(recipes.to_string())")
    nb.md("## 2. Build every map",
          "3D UMAP (umap-learn on the CPU, seed 42, robust scaling, median imputation, each block "
          "scaled to equal total variance), then HDBSCAN with 'leaf' selection.",
          "Why 'leaf' and not the Clusters tab's opening setting (25/25, 'eom'): on these maps "
          "excess-of-mass returns one to a few clusters holding almost every gene, so every label "
          "scores at chance. The comparison cell below measures both on the shipped maps.",
          "The build is checkpointed one map at a time into `results/umap_gallery_checkpoint`, so a "
          "run stopped by the memory guard resumes where it left off.")
    nb.code(f"built = B.build({codes!r}, workers={int(workers)})",
            "maps = B.maps_table(built)",
            "print(maps.groupby('organism')[['genes', 'clusters', 'structure']].describe().round(3)"
            ".to_string())",
            "print(maps.to_string())")
    nb.md("### What the structure search was worth",
          "`gain_on_sample` is the best UMAP setting against the neighbourhood the gallery used to "
          "fix for every map (n_neighbors 25, min_dist 0.25), measured at the sample stage where all "
          "of them were scored on the same genes. `cluster_spread` is the range the clustering grid "
          "covered on the winning embedding. `chose_old_default` is how often the search came back "
          "to the setting that used to be hard-coded -- if that were always true, the grid would be "
          "hours wasted.")
    nb.code("gain = B.search_gain(built)",
            "print(gain.groupby('organism')[['shipped', 'sample_best', 'sample_at_old_default', "
            "'gain_on_sample', 'spread_full']].describe().round(3).to_string())",
            "print(gain.groupby('organism')['chose_old_default'].mean().round(3).to_string())",
            "print(gain.sort_values('gain_on_sample', ascending=False).head(15).to_string())")
    nb.code("searched = B.search_table(built)",
            "print(len(searched), 'settings tried in total')",
            "print(searched.groupby(['organism', 'stage']).size().to_string())",
            "print(searched.groupby(['n_neighbors', 'min_dist'])['structure'].mean().round(3)"
            ".to_string())",
            "print(searched.groupby(['min_cluster_size', 'min_samples'])['structure'].mean()"
            ".round(3).to_string())")
    nb.md("### The best and the worst maps by structure")
    nb.code("for code, t in maps.groupby('organism'):",
            "    print(code, 'BEST'); print(t.head(8)[['map', 'group', 'genes', 'clusters', "
            "'clustered', 'evenness', 'separation', 'structure']].to_string(index=False))",
            "    print(code, 'WORST'); print(t.tail(8)[['map', 'group', 'genes', 'clusters', "
            "'clustered', 'evenness', 'separation', 'structure']].to_string(index=False))")
    nb.code("compare = B.compare_clusterings(built)",
            "print(compare.groupby(['organism', 'setting'])[['clusters', 'unclustered', "
            "'category_skill']].median().to_string())")
    for code in codes:
        nb.md(f"### The 24 best-structured {code} maps, first two components, colored by cluster")
        nb.code(f"B.figure(built, {code!r})")
    nb.md("## 3. Score every label on every map",
          "Categories → clusters: size-weighted mean over categories of the F1 of each category's "
          "best cluster (unclustered genes count as missed). Best category: the category whose F1 of "
          "the Wilson 95% lower bounds of precision and recall beats its own shuffled-label chance "
          f"by the most. Chance: {G.PERMUTATIONS} shuffles of the label over the same genes. "
          "Circular: the map was built from a column in the label's leakage closure.")
    nb.code("scores = B.score(built)",
            "print(len(scores), 'label x map rows')",
            "scores.groupby('organism').size()")
    for code in codes:
        nb.md(f"### {code}: its first calibration target on every map",
              "Best first by 'categories → clusters' skill. A map marked circular was built from "
              "a column that restates the label.")
        nb.code(f"label = B.headline_label({code!r})",
                "print(label)",
                f"print(G.describe(scores[scores.organism == {code!r}], label))")
    nb.md("## 4. The size-aware bound on real data",
          "Across every label and map: how often the best category by raw F1 is a small category, "
          "and how often it still is once the Wilson bound and the chance correction decide.")
    nb.code("small = scores[scores.best_n < 30]",
            "print(f'best categories under 30 genes: {len(small)} of {len(scores)}')",
            "print(small[['organism', 'label', 'map', 'best_category', 'best_n', 'best_precision', "
            "'best_recall', 'best_f1_lower', 'best_skill']].round(3).head(15).to_string())")
    nb.md("## 5. Which map wins for each label",
          "The point of a gallery this size: for every label, the map its categories fall into "
          "clusters on best, and whether that map is one a user would have thought of. A map marked "
          "circular was built from a column in the label's leakage closure.")
    nb.code("best = (scores[scores.circular.astype(str) != 'True']",
            "        .sort_values('category_skill', ascending=False)",
            "        .groupby(['organism', 'label']).head(1))",
            "print(best[['organism', 'label', 'map', 'category_f1', 'category_skill', 'genes_scored'"
            "]].round(3).to_string(index=False))")
    nb.md("## 6. Write the shipped files")
    nb.code("paths = B.write(built, scores)",
            "import os",
            "print({k: round(os.path.getsize(p) / 1e6, 2) for k, p in paths.items()}, 'MB')")
    return nb.write(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", action="append", help="space code (repeatable); default all")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-notebook", action="store_true", help="build and write without a notebook")
    ap.add_argument("--no-checkpoint", action="store_true",
                    help="ignore and do not write the per-map checkpoint; rebuild everything")
    ap.add_argument("--out", help="write the gallery here instead of into the package")
    a = ap.parse_args(argv)
    if a.no_notebook:
        built = build(a.organism, a.workers, checkpoint=None if a.no_checkpoint else CHECKPOINT)
        print(maps_table(built).to_string())
        print(search_gain(built).to_string())
        write(built, score(built), directory=a.out)
    else:
        print(notebook(a.organism, a.workers))


if __name__ == "__main__":
    main()
