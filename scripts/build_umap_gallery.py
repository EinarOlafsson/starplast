"""Build the pregenerated map gallery and its label x map score table.

For every parasite space whose table ships (`starplast.organisms`), this builds the maps listed by
`starplast.umap_gallery.recipes`: all measurements at three n_neighbors settings, all but
localization, and one map per evidence family with enough columns and genes. Each is a 3D UMAP of
measurements only (label columns are removed from the table first), clustered with HDBSCAN
(`umap_gallery.CLUSTERING`). Every categorical label is then scored on every map
(`umap_gallery.score_gallery`).

Output, in `starplast/data/`:

* `umap_gallery.npz`     -- per map: placed gene positions (int32), coordinates (float32), clusters (int16)
* `umap_gallery.json`    -- per map: recipe, blocks, source columns, sizes, backend, versions
* `umap_gallery_scores.tsv` -- the label x map table

and the executed notebook `notebooks/umap_gallery_2026_09_29.ipynb` that records the run.

The embedding runs on the CPU (`STARPLAST_GPU=0`) so the shipped maps are the umap-learn maps a
user without CUDA would reproduce from the stored recipe. Maps are built in parallel processes;
run the whole thing under a memory cap:

    systemd-run --user --scope -p MemoryMax=30G env STARPLAST_GPU=0 \\
        ~/anaconda3/envs/starplast/bin/python scripts/build_umap_gallery.py --workers 8

    python scripts/build_umap_gallery.py --no-notebook --organism <code>   # quick, no notebook
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("STARPLAST_GPU", "0")

import pandas as pd  # noqa: E402

from starplast import umap_gallery as G  # noqa: E402

NOTEBOOK = os.path.join(ROOT, "notebooks", "umap_gallery_2026_09_29.ipynb")

_CTX: dict = {}


def context(code: str):
    """The shipped context of one organism, once per process."""
    if code not in _CTX:
        from starplast import strategies as S
        _CTX[code] = S.Context.shipped(code)
    return _CTX[code]


def _one(job):
    code, recipe = job
    t0 = time.monotonic()
    m = G.build_map(context(code), recipe, log=lambda *a: None)
    m["info"]["seconds"] = round(time.monotonic() - t0, 1)
    return code, recipe["id"], m


def build(codes=None, workers: int = 6, log=print) -> dict:
    """organism -> map id -> built map, in recipe order."""
    codes = list(codes or G.space_codes())
    jobs = [(code, r) for code in codes for r in G.recipes(context(code))]
    log(f"{len(jobs)} maps over {', '.join(codes)} on {workers} worker(s)")
    done = {}
    with ProcessPoolExecutor(max_workers=max(1, int(workers))) as pool:
        for code, map_id, m in pool.map(_one, jobs):
            log(f"  {code} {map_id}: {m['info']['n_genes']:,} genes, "
                f"{m['info']['n_features']} features, {m['info']['n_clusters']} clusters, "
                f"{m['info']['noise_fraction']:.0%} unclustered ({m['info']['seconds']} s)")
            done[(code, map_id)] = m
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
    """One row per map: what it is built from and how it clustered."""
    rows = []
    for code, maps in built.items():
        for map_id, m in maps.items():
            i = m["info"]
            rows.append({"organism": code, "map": map_id, "family": i["family"],
                         "blocks": len(i["blocks"]), "columns": len(i["columns"]),
                         "features": i["n_features"], "n_neighbors": i["n_neighbors"],
                         "genes": i["n_genes"], "clusters": i["n_clusters"],
                         "unclustered": round(i["noise_fraction"], 3), "seconds": i["seconds"]})
    return pd.DataFrame(rows)


def write(built: dict, scores: pd.DataFrame) -> dict:
    """Write the three shipped files and report their sizes."""
    ids = {code: context(code).gene_ids for code in built}
    paths = G.write(built, scores, ids, extra={
        "built": _dt.date.today().isoformat(),
        "script": "scripts/build_umap_gallery.py",
        "notebook": os.path.relpath(NOTEBOOK, ROOT),
        "scoring": {"min_category": G.MIN_CATEGORY, "permutations": G.PERMUTATIONS,
                    "bound": "Wilson 95% lower bounds of precision and recall, combined as F1",
                    "noise": "HDBSCAN-unclustered genes count in a category's size (missed)"}})
    for k, p in paths.items():
        print(f"{k:9s} {os.path.relpath(p, ROOT)}  {os.path.getsize(p) / 1e6:.2f} MB")
    return paths


def figure(built: dict, code: str):
    """Every map of one organism, first two components, colored by its clusters."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    maps = built[code]
    n = len(maps)
    cols = 4
    fig, axes = plt.subplots((n + cols - 1) // cols, cols, figsize=(3.2 * cols, 3.0 * ((n + cols - 1) // cols)))
    for ax in np.ravel(axes):
        ax.axis("off")
    for ax, (map_id, m) in zip(np.ravel(axes), maps.items()):
        lab = m["clusters"]
        xyz = m["xyz"]
        ax.scatter(xyz[lab < 0, 0], xyz[lab < 0, 1], s=1, c="#bbbbbb", linewidths=0)
        ax.scatter(xyz[lab >= 0, 0], xyz[lab >= 0, 1], s=1, c=lab[lab >= 0], cmap="tab20",
                   linewidths=0)
        ax.set_title(f"{map_id}\n{m['info']['n_genes']:,} genes, {m['info']['n_clusters']} clusters",
                     fontsize=8)
    fig.tight_layout()
    return fig


def notebook(codes=None, workers: int = 6, path: str = NOTEBOOK) -> str:
    """Run the whole build inside an executed notebook and write it."""
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    from notebook_runner import ExecutedNotebook
    codes = list(codes or G.space_codes())
    nb = ExecutedNotebook("Pregenerated map gallery and label scores", {"__name__": "__notebook__"})
    nb.md("The application's central map is one choice of inputs. This notebook builds a small "
          "gallery of alternatives for each organism -- all measurements at three neighbourhood "
          "sizes, all but localization, and one map per evidence family -- and scores how well every "
          "categorical label maps onto each. The code is `scripts/build_umap_gallery.py` and "
          "`starplast/umap_gallery.py`; this notebook calls them, so what is shown is what ships.",
          "Every map is built from measurements only: label columns are removed from the table "
          "before embedding, and the build refuses a map whose matrix contains one.")
    nb.code("import warnings; warnings.filterwarnings('ignore')",
            f"import sys; sys.path.insert(0, {ROOT!r}); "
            f"sys.path.insert(0, {os.path.join(ROOT, 'scripts')!r})",
            "import build_umap_gallery as B",
            "from starplast import umap_gallery as G",
            "import pandas as pd",
            "_ = pd.set_option('display.width', 200)")
    nb.md("## 1. The recipes",
          f"A family map needs at least {G.MIN_FAMILY_COLUMNS} source columns and "
          f"{G.MIN_FAMILY_GENES} genes measured in at least {G.MIN_GENE_COVERAGE:.0%} of them; only "
          "those genes are placed. The 'all' maps place every gene, as the application does.")
    nb.code("rows = []",
            f"for code in {codes!r}:",
            "    ctx = B.context(code)",
            "    for r in G.recipes(ctx):",
            "        rows.append({'organism': code, 'map': r['id'], 'blocks': len(r['blocks']),",
            "                     'columns': len(r['columns']), 'genes': len(G.genes_for(ctx, r)),",
            "                     'n_neighbors': r['n_neighbors']})",
            "recipes = pd.DataFrame(rows)",
            "print(recipes.to_string())")
    nb.md("## 2. Build every map",
          "3D UMAP (umap-learn on the CPU, seed 42, min_dist 0.25, robust scaling, median "
          "imputation, each block scaled to equal total variance), then HDBSCAN with "
          f"{G.CLUSTERING}.",
          "Why not the Clusters tab's opening setting (25/25, 'eom'): on these maps it returns one "
          "to a few clusters holding almost every gene, so every label scores at chance. The next "
          "cell measures both on the shipped maps before the gallery is written.")
    nb.code(f"built = B.build({codes!r}, workers={int(workers)})",
            "maps = B.maps_table(built)",
            "print(maps.to_string())")
    nb.code("compare = B.compare_clusterings(built)",
            "print(compare.to_string())",
            "print(compare.groupby(['organism', 'setting'])[['clusters', 'unclustered', "
            "'category_skill']].median())")
    for code in codes:
        nb.md(f"### The {code} maps, first two components, colored by cluster")
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
    nb.md("## 5. Write the shipped files")
    nb.code("paths = B.write(built, scores)")
    return nb.write(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--organism", action="append", help="space code (repeatable); default all")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-notebook", action="store_true", help="build and write without a notebook")
    a = ap.parse_args(argv)
    if a.no_notebook:
        built = build(a.organism, a.workers)
        print(maps_table(built).to_string())
        write(built, score(built))
    else:
        print(notebook(a.organism, a.workers))


if __name__ == "__main__":
    main()
