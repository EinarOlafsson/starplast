#!/usr/bin/env python3
"""The discovery search from a terminal, with no window anywhere in it.

The Discover tab is the wrong place to start a search that takes an hour. A window has to stay open,
the machine has to stay awake, and the run dies with the session -- which is exactly what happened
the first time this was run at scale. A search is a batch job that happens to have a GUI, not a GUI
feature, so it gets a command.

    starplast-discover --task guilt:compartment_best --budget 100
    starplast-discover --task disagreement:compartment_best:fit_invitro_hff --task knn:lopit_unified
    starplast-discover --list
    starplast-discover --read big_00_localisation

Nothing here imports Qt. That is the point and it is worth stating: a headless option that pulls in
a GUI toolkit fails on the machine it exists for, which is a server with no display.

## Why one process per task rather than a loop inside one

Each `--task` is saved the moment it finishes, and a task already saved under the same name is
skipped. So a run that is interrupted -- killed, timed out, or force-quit by somebody who saw an
unfamiliar process eating a GPU -- costs the task in flight and nothing else, and the same command
run again continues where it stopped. That property matters more than tidiness here: the failure
this is built around is a long job dying two thirds of the way through.
"""
from __future__ import annotations

import argparse
import sys

#: A task is `mode:layer` or `mode:layer:against`, which is enough to name every search this
#: application can run and short enough to type twice.
TASK_HELP = ("what to search for, as mode:layer or mode:layer:against -- "
             "e.g. guilt:compartment_best, conjunction:compartment_best:cellcycle_phase. "
             "Repeatable; each is saved separately as it finishes.")


def parse_task(text: str) -> dict:
    """`mode:layer[:against]` as a dict, or raise with what was wrong."""
    from .optimize import MODES
    bits = [b.strip() for b in str(text).split(":")]
    if len(bits) < 2 or not all(bits[:2]):
        raise ValueError(f"task {text!r} is not mode:layer or mode:layer:against")
    if bits[0] not in MODES:
        raise ValueError(f"unknown mode {bits[0]!r}; expected one of {', '.join(MODES)}")
    return {"mode": bits[0], "layer": bits[1], "against": bits[2] if len(bits) > 2 else ""}


def load_nodes(path: str = None):
    """The node table, without importing the window.

    `app.load_graph` is the usual way in and it imports Qt on the way past, which is exactly what a
    machine with no display cannot do. The parquet is read directly here, and the error names the
    command that builds it rather than letting pandas raise from inside a constructor.
    """
    import os
    import pandas as pd
    from . import paths
    where = path or os.path.join(paths.data_dir(), "nodes.parquet")
    if not os.path.exists(where):
        raise SystemExit(f"No node table at {where}. Run:  python -m starplast.build_graph")
    return pd.read_parquet(where)


def store_for(root: str = None):
    """Where searches are kept: the same place the window looks, unless told otherwise."""
    import os
    from . import paths
    from .searches import SearchStore
    return SearchStore(root or os.path.join(paths.user_cache_dir(), "searches"))


def run_task(nodes, task: dict, budget: int = 100, restarts: int = 3, seed: int = 42,
             blocks=None, exclude=(), name: str = "", store=None, log=print) -> object:
    """One search, saved under `name`. Returns the `Search`, or None if it was already done."""
    from . import gpu, optimize, searches
    store = store if store is not None else store_for()
    name = name or f"{task['mode']}_{task['layer']}"
    if any(existing == name for existing, _manifest in store.list()):
        log(f"{name}: already saved, skipping")
        return None
    # A layer that is not in the table would run the whole search and find nothing, which reads as
    # "there is nothing there" rather than as "you typed the column name wrong". Caught here, with
    # the near misses offered, because the cost of the mistake is an hour of GPU time.
    for field in ("layer", "against"):
        want = task.get(field)
        if want and want not in nodes.columns:
            near = [c for c in nodes.columns if want.lower() in str(c).lower()][:6]
            raise ValueError(f"{field} {want!r} is not a column in this table"
                             + (f" -- did you mean {', '.join(near)}?" if near else ""))
    pool = list(blocks) if blocks else [b for b in optimize.block_pool(nodes)
                                        if b not in set(exclude)]
    if not pool:
        raise ValueError("no feature blocks left to search over")
    start = {**{k: v[len(v) // 2] for k, v in optimize.SPACE.items()},
             "method": "umap", "algorithm": "hdbscan", "blocks": tuple(pool[:3])}
    log(f"{name}: {task['mode']} on {task['layer']}"
        + (f" against {task['against']}" if task["against"] else "")
        + f", up to {budget} configurations over {len(pool)} blocks ({', '.join(pool)})")
    with gpu.pinned() as backend:
        evaluate = optimize.evaluator(nodes, mode=task["mode"], layers=(task["layer"],),
                                      against=(task["against"],) if task["against"] else (),
                                      seed=seed, log=log)
        result = optimize.climb(evaluate, start, block_pool=pool, restarts=restarts,
                                max_evaluations=budget, seed=seed, log=log)
    run = searches.from_climb(result, nodes, name=name, configs=int(len(result)), seed=seed,
                              backend=backend, **task)
    where = store.save(run)
    log(f"{name}: saved to {where}")
    return run


def _summary(run) -> str:
    """One line a batch log can be read from: what won, and what it found."""
    if run is None or run.configs.empty:
        return "nothing evaluated"
    b = run.best
    winner_count = len(run.for_config(0)[1])
    total_count = len(run.findings)
    finding_text = f"{winner_count} findings"
    if total_count != winner_count:
        finding_text += f" ({total_count} across all configs)"
    bits = [f"best {float(b.score):.2f}", f"{len(run.configs)} configs", finding_text]
    for column, label in (("n_clusters", "clusters"), ("mean_auprc", "auprc"),
                          ("mean_lift", "lift"), ("algorithm", ""), ("method", "")):
        if column in run.configs.columns:
            value = b.get(column)
            if value is not None and str(value) not in ("nan", "None"):
                bits.append(f"{label} {value}".strip() if label else str(value))
    return ", ".join(bits)


def _say(*bits) -> None:
    """Print and FLUSH.

    Python block-buffers stdout when it is a pipe, which is how every long headless run is watched:
    `starplast-discover ... | tee run.log`, or a log file read from another window. Without the
    flush the per-task summaries sit in the buffer for the length of the batch and the run looks
    stalled -- and the one thing a batch job has to do is say where it has got to.
    """
    print(*bits, flush=True)


def main(argv=None) -> int:
    """`starplast-discover`. Returns a process exit code."""
    p = argparse.ArgumentParser(
        prog="starplast-discover",
        description="Search the space of maps for structure that predicts, without a window.")
    p.add_argument("--task", action="append", default=[], metavar="MODE:LAYER[:AGAINST]",
                   help=TASK_HELP)
    p.add_argument("--budget", type=int, default=100,
                   help="configurations per task (default 100)")
    p.add_argument("--restarts", type=int, default=3, help="hill-climb restarts (default 3)")
    p.add_argument("--seed", type=int, default=42, help="random seed (default 42)")
    p.add_argument("--blocks", default="", help="only these feature blocks, comma separated")
    p.add_argument("--exclude", default="literature",
                   help="feature blocks to keep out (default: literature, which counts study "
                        "effort rather than biology)")
    p.add_argument("--name", default="", help="name for a single task's saved search")
    p.add_argument("--prefix", default="", help="name prefix when several tasks are given")
    p.add_argument("--out", default="", help="where to keep searches (default: the app's cache)")
    p.add_argument("--nodes", default="", metavar="PATH",
                   help="a node table to search, if not the cached one")
    p.add_argument("--list", action="store_true", help="list saved searches and stop")
    p.add_argument("--read", default="", metavar="NAME",
                   help="print the reading of a saved search and stop")
    p.add_argument("--config", type=int, default=0,
                   help="which configuration to read, 0 being the winner (default 0)")
    p.add_argument("--top", type=int, default=8, help="findings to write out (default 8)")
    p.add_argument("--quiet", action="store_true", help="only the per-task summary lines")
    args = p.parse_args(argv)

    store = store_for(args.out or None)
    if args.list:
        rows = store.list()
        if not rows:
            _say("no saved searches")
        for name, manifest in rows:
            _say(f"{name}  ·  {manifest.get('mode', '?')} / {manifest.get('layer', '?')}"
                  f"  ·  {manifest.get('configs', '?')} configs"
                  f"  ·  {manifest.get('created', '')}")
        return 0

    if args.read:
        run = store.load(args.read)
        if run.configs.empty:
            print(f"no search called {args.read!r}", file=sys.stderr)
            return 1
        _say(run.report(load_nodes(args.nodes or None), index=args.config, top=args.top))
        return 0

    if not args.task:
        p.error("nothing to do: give at least one --task, or --list, or --read")
    try:
        tasks = [parse_task(t) for t in args.task]
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    nodes = load_nodes(args.nodes or None)
    log = (lambda *_a, **_k: None) if args.quiet else _say
    _say(f"{len(nodes):,} genes, {nodes.shape[1]} columns")
    failed = 0
    for i, task in enumerate(tasks):
        name = args.name if (args.name and len(tasks) == 1) else \
            f"{args.prefix}{i:02d}_{task['mode']}_{task['layer']}" if args.prefix else ""
        try:
            run = run_task(nodes, task, budget=args.budget, restarts=args.restarts,
                           seed=args.seed,
                           blocks=[b for b in args.blocks.split(",") if b] or None,
                           exclude=[b for b in args.exclude.split(",") if b],
                           name=name, store=store, log=log)
        except Exception as exc:                  # noqa: BLE001 -- one task must not end the batch
            # A batch that dies on its third task and takes the other seven with it is the failure
            # this command exists to avoid.
            print(f"{task['mode']}:{task['layer']}: FAILED {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            failed += 1
            continue
        if run is not None:
            _say(f"{run.name}: {_summary(run)}")
    return 1 if failed else 0


if __name__ == "__main__":                        # pragma: no cover -- exercised through main()
    raise SystemExit(main())
