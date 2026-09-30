#!/usr/bin/env python3
"""Precompute the star map's shipped links: the default-setting strategy runs on the shipped tables.

For every parasite space whose table ships, each strategy of `star_edges.SHIPPED_STRATEGIES` is run
once with its defaults, and the gene-gene links its result states are kept (bounded as described in
`starplast/star_edges.py`). The set-expansion strategies (20, 25) need a seed list; they are run on
the `SHIPPED_SEED_SETS` categories of the default label whose size is closest to what
`strategies.example_set` aims for -- the same kind of list the panel's "example" button fills.

Nothing here changes a strategy or a calibration number: the strategies are run, not altered.

    python scripts/build_star_edges.py              # build, write, and record the notebook
    python scripts/build_star_edges.py --no-notebook

Heavy (UMAP walks); run it under a memory cap:
    systemd-run --user --scope -p MemoryMax=30G python scripts/build_star_edges.py
"""
from __future__ import annotations

import argparse
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

NOTEBOOK = os.path.join(ROOT, "notebooks", "star_edges_2026_09_29.ipynb")


def seed_sets(ctx, n: int, lo: int = 30, hi: int = 300, aim: int = 120) -> list:
    """(category, positions) for the `n` categories of the default label sized closest to `aim`."""
    from starplast import strategies as S
    target = S.default_category(ctx)
    if not target:
        return []
    t = ctx.truth(target)
    counts = t.value_counts()
    ok = counts[(counts >= lo) & (counts <= hi)]
    order = sorted(ok.index, key=lambda c: (abs(int(ok[c]) - aim), str(c)))[:n]
    return [(str(c), np.flatnonzero((t == c).to_numpy())) for c in order]


def build_space(code: str, keys=None, log=print) -> tuple:
    """(links, log frame) for one space: every shipped strategy run once, its links kept."""
    from starplast import star_edges as E
    from starplast import strategies as S
    import starplast.strategy_catalog  # noqa: F401  (registers the strategies)
    import starplast.strategy_graph  # noqa: F401
    import starplast.strategy_learning  # noqa: F401
    ctx = S.shipped(code)
    parts, rows = [], []
    for key in keys or E.SHIPPED_STRATEGIES:
        strategy = S.get(key)
        needs_genes = any(p.kind == "genes" and p.name == "genes" for p in strategy.params)
        runs = ([(f"starplast:{key}:{cat}", {"genes": " ".join(ctx.gene_ids[pos])}, cat)
                 for cat, pos in seed_sets(ctx, E.SHIPPED_SEED_SETS)] if needs_genes
                else [(f"starplast:{key}", {}, "")])
        for run, kw, cat in runs:
            t0 = time.monotonic()
            try:
                result = strategy.run(ctx, **kw)
            except Exception as exc:          # a strategy this space cannot run is recorded
                rows.append({"organism": code, "strategy": key, "run": run, "links": 0,
                             "seconds": round(time.monotonic() - t0, 1),
                             "note": f"{type(exc).__name__}: {exc}"[:160]})
                continue
            links = E.edges_from_result(result, ctx, run=run, origin=E.ORIGIN_SHIPPED)
            if cat:
                links["setting"] = (f"seed list = {cat!r} ({len(kw['genes'].split())} genes); "
                                    + links["setting"]).str.slice(0, 160)
            parts.append(E.to_stored(links, ctx.gene_ids, code))
            rows.append({"organism": code, "strategy": key, "run": run, "links": len(links),
                         "seconds": round(time.monotonic() - t0, 1),
                         "note": str(result.summary)[:160]})
            log(f"{code} {run}: {len(links):,} links in {rows[-1]['seconds']}s")
    frame = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=E.STORED)
    return frame, pd.DataFrame(rows)


def build(codes=None, keys=None, log=print) -> tuple:
    """(links, log) for every shipped parasite space."""
    from starplast import organisms
    codes = codes or organisms.codes(organisms.PARASITE, available=True)
    frames, logs = [], []
    for code in codes:
        f, l = build_space(code, keys, log)
        frames.append(f)
        logs.append(l)
    return pd.concat(frames, ignore_index=True), pd.concat(logs, ignore_index=True)


def write(frame: pd.DataFrame, path: str | None = None) -> str:
    """Write the shipped file, compactly (categorical columns, float32 scores, zstd)."""
    from starplast import star_edges as E
    path = path or E.shipped_path()
    out = frame.copy()
    for c in ("organism", "kind", "group", "source", "how", "setting", "run", "origin", "created",
              "note", "gene_a", "gene_b"):
        out[c] = out[c].astype(str).astype("category")
    out.to_parquet(path, index=False, compression="zstd")
    return path


def notebook(path: str = NOTEBOOK) -> str:
    """Build, write and record everything as an executed notebook."""
    from notebook_runner import ExecutedNotebook
    nb = ExecutedNotebook("Star map: the shipped strategy links", {"__name__": "__notebook__"})
    nb.md("The star map draws gene-gene links with their provenance. Measured links are read from "
          "each space's graph at run time. This notebook makes the other shipped part: the links "
          "stated by the default-setting runs of the strategies that relate genes to genes, on the "
          "shipped tables. The code is `scripts/build_star_edges.py` and the rules that bound each "
          "kind of link are in `starplast/star_edges.py`; this notebook calls them, so what is "
          "shown is what ships in `starplast/data/star_edges.parquet`.",
          "No strategy or calibration number is changed: each strategy is run with its defaults.")
    nb.code("import warnings; warnings.filterwarnings('ignore')",
            f"import sys; sys.path.insert(0, {ROOT!r}); "
            f"sys.path.insert(0, {os.path.join(ROOT, 'scripts')!r})",
            "import pandas as pd",
            "import build_star_edges as B",
            "from starplast import star_edges as E",
            "pd.DataFrame({'bound': ['MODULE_K', 'SEED_K', 'LABEL_K', 'PARTNER_MAX', "
            "'MAX_EDGES_PER_RUN', 'SHIPPED_SEED_SETS'], 'value': [E.MODULE_K, E.SEED_K, "
            "E.LABEL_K, E.PARTNER_MAX, E.MAX_EDGES_PER_RUN, E.SHIPPED_SEED_SETS]})")
    nb.md("## 1. Run every strategy once per space",
          "Seed-list strategies (20, 25) run once per example seed list.")
    nb.code("links, log = B.build(log=lambda m: None)",
            "log")
    nb.md("## 2. Links per source")
    nb.code("links.groupby(['organism', 'source', 'how']).size().rename('links').reset_index()")
    nb.md("## 3. What a link carries", "Three examples: a predicted pair, a module link and a "
          "seed link.")
    nb.code("pd.concat([links[links['how'] == h].head(1) for h in ('pair', 'module', 'seed')]).T")
    nb.md("## 4. Write the shipped file")
    nb.code("import os",
            "path = B.write(links)",
            "print(path, f'{os.path.getsize(path) / 1e6:.2f} MB', f'{len(links):,} links')")
    return nb.write(path)


def main(argv=None) -> int:
    """Command line: build and write, recording the notebook unless told not to."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-notebook", action="store_true", help="build and write only")
    ap.add_argument("--only", default="", help="comma-separated strategy keys (testing)")
    args = ap.parse_args(argv)
    if args.only:                               # a partial build is never written over the full one
        links, log = build(keys=[k for k in args.only.split(",") if k])
        print(log.to_string())
    elif args.no_notebook:
        links, log = build()
        print(log.to_string())
        print(write(links))
    else:
        print(notebook())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
