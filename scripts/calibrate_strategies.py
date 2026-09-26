#!/usr/bin/env python3
"""Calibrate every strategy over a large space of settings, targets and seeds.

A strategy's self-test at its default settings says whether it works THERE. This asks where it
works: each strategy is run over a grid of its own parameters -- neighbourhood sizes, thresholds,
edge layers, models, cluster selections -- and over several held-out labels, each configuration
under several seeds (a seed changes which genes are hidden and how every map is drawn). Every
configuration therefore gets a mean and a 95% interval, not a single draw.

Results stream to ``results/calibration_<date>/runs.jsonl`` one line per run, so a stopped sweep
resumes where it stopped: a configuration already in the file is not run again. `--summarize`
turns the file into the tables the README and the Strategies tab read:

    summary.csv          one row per (organism, strategy, configuration): metric mean and 95% CI
                         across seeds, chance level, skill, pass rate
    best.csv             per strategy, the configuration with the highest lower confidence bound
    defaults.csv         per strategy, its default configuration, for the README table
    sensitivity.csv      per strategy and parameter, the mean skill at each value

"Skill" puts every metric on one scale: (observed - chance) / (1 - chance), so 0 is chance and 1
is perfect, whether the metric is a correct-call rate, an AUROC, an F1 or a replication share.

Run:  python scripts/calibrate_strategies.py --workers 12
      python scripts/calibrate_strategies.py --summarize
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import itertools
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

#: Seeds per configuration. Five gives a usable interval without a sweep of days.
SEEDS = (1, 2, 3, 4, 5)

#: Held-out labels per organism, for the strategies that take one. The coarse and fine
#: localization, a structural class, a phenotype screen and a derived stage label on T. gondii;
#: the transferred knockout phenotype, a derived stage and the export call on P. falciparum.
TARGETS = {"Tg": ("compartment", "lopit_unified", "dtm_class", "screen_any_phenotype",
                  "stage_enriched_derived"),
           "Pf": ("pb_transferred_phenotype", "stage_enriched_derived", "is_exported")}
NUMBERS = {"Tg": ("fit_invitro_hff", "fit_invivo_PE", "expr_tachy"),
           "Pf": ("piggybac_mis", "expr_schizont", "mean_plddt")}

#: The grid per strategy: parameter -> values. A parameter absent here keeps its default. "target"
#: and "number" expand to the organism's lists above. Kept to what changes the answer: every value
#: multiplies the sweep by itself times five seeds.
GRID = {
    "holdout_search": {"target": "TARGETS", "n_neighbors": ["15, 50", "10, 30", "25, 100"],
                       "min_cluster_size": ["20, 50", "10, 25"], "selection": ["eom, leaf"]},
    "geneset_hunt": {"exclude": "TARGETS", "n_neighbors": ["15, 50", "10, 30"],
                     "min_cluster_size": ["20, 50", "10, 25"]},
    "recoverability_atlas": {"target": "TARGETS", "k": [5, 15, 50]},
    "consensus_modules": {"target": "TARGETS", "threshold": [0.3, 0.5, 0.8]},
    "blind_battery": {"map_from": "FAMILIES"},
    "block_ablation": {"target": "TARGETS", "k": [5, 15, 50]},
    "feature_knn": {"target": "TARGETS", "k": [5, 15, 50], "min_share": [0.0, 0.3, 0.6]},
    "map_neighbours": {"target": "TARGETS", "n_neighbors": [10, 25, 60], "k": [5, 15, 50]},
    "cluster_guilt": {"target": "TARGETS", "min_cluster_size": [10, 25, 60],
                      "selection": ["leaf", "eom"], "min_lift": [1.5, 3.0]},
    "label_outliers": {"target": "TARGETS", "k": [5, 15, 50]},
    "layer_propagation": {"target": "TARGETS", "layer": "LAYERS", "restart": [0.2, 0.5, 0.8]},
    "layer_vote": {"target": "TARGETS", "k": [5, 15, 50]},
    "physical_partners": {"target": "TARGETS"},
    "structural_homology": {"level": [1, 2, 3]},
    "multiplex_modules": {"target": "TARGETS", "resolution": [0.5, 1.0, 2.0],
                          "agreement": [0.3, 0.5]},
    "link_prediction": {"layer": "LAYERS"},
    "attention_correction": {"target": "TARGETS", "top": [50, 200, 1000]},
    "unwritten_links": {"min_layers": [1, 2, 3]},
    "supervised_classifier": {"target": "TARGETS", "C": [0.01, 0.1, 1.0, 10.0]},
    "positive_unlabeled": {"exclude": "TARGETS", "bags": [5, 15, 40]},
    "trait_regression": {"target": "NUMBERS", "model": ["boosted", "ridge"],
                         "own_kind": ["leave out", "include"]},
    "masked_imputation": {"rank": [5, 20, 60]},
    "condition_shift": {"model": ["boosted", "ridge"], "own_kind": ["leave out", "include"]},
    "set_enrichment": {"exclude": "TARGETS"},
    "seed_expansion": {"exclude": "TARGETS", "mode": ["networks", "networks + measurements"],
                       "restart": [0.1, 0.3, 0.6]},
    "split_clusters": {"min_cluster_size": [15, 40, 80]},
    "conjunctions": {"min_cluster_size": [8, 15, 30]},
    "paralog_divergence": {"target": "TARGETS", "min_shared": [5, 10, 30]},
    "ortholog_transfer": {},
    "stratum_focus": {"target": "TARGETS", "stratum": ["lineage-specific", "hypothetical protein",
                                                       "understudied", "conserved"]},
    "triangulation": {"target": "TARGETS", "min_agree": [1, 2, 3]},
    "understudied_first": {"target": "TARGETS", "min_agree": [1, 2, 3]},
}


def _expand(values, organism, ctx):
    if values == "TARGETS":
        return [t for t in TARGETS[organism] if t in ctx.categorical_columns()]
    if values == "NUMBERS":
        return [t for t in NUMBERS[organism] if t in ctx.numeric_columns()]
    if values == "LAYERS":
        return [l for l in ctx.layers() if l not in ("compartment",)]
    if values == "FAMILIES":
        return sorted(ctx.families())
    return list(values)


def configurations(organism: str, ctx, only=None) -> list:
    """Every (strategy, settings, seed) of the sweep for one organism, in a stable order."""
    from starplast import strategies as S
    out = []
    for s in S.catalog():
        if only and s.key not in only:
            continue
        grid = GRID.get(s.key, {})
        names = [p for p in grid if any(q.name == p for q in s.params)]
        pools = [_expand(grid[p], organism, ctx) for p in names]
        combos = list(itertools.product(*pools)) if names else [()]
        for combo in combos:
            settings = dict(zip(names, combo))
            for seed in SEEDS:
                out.append({"organism": organism, "strategy": s.key, "settings": settings,
                            "seed": seed})
    return out


def config_id(cfg: dict) -> str:
    """A stable identifier for one run, so a resumed sweep can skip what it has done."""
    text = json.dumps({k: cfg[k] for k in ("organism", "strategy", "settings", "seed")},
                      sort_keys=True, default=str)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


_CTX: dict = {}


def _context(organism: str, seed: int):
    """One shipped context per organism per worker, re-seeded per run, caches kept."""
    from starplast import strategies as S
    base = _CTX.get(organism)
    if base is None:
        base = S.Context.shipped(organism)
        _CTX[organism] = base
    ctx = base.bound()
    ctx.seed = int(seed)
    return ctx


def run_one(cfg: dict) -> dict:
    """Run one configuration's self-test and return its record. Failures are records too."""
    import warnings
    warnings.filterwarnings("ignore")
    from starplast import strategies as S
    t0 = time.monotonic()
    rec = {"id": config_id(cfg), **cfg}
    try:
        ctx = _context(cfg["organism"], cfg["seed"])
        r = S.get(cfg["strategy"]).test(ctx, **cfg["settings"])
        rec.update(r.to_dict())
    except Exception as exc:                       # a configuration that cannot run is a finding
        rec.update({"verdict": "NOT RUN", "note": f"{type(exc).__name__}: {exc}"})
    rec["wall_seconds"] = round(time.monotonic() - t0, 2)
    return rec


def _init_worker():
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMBA_NUM_THREADS"):
        os.environ[var] = "2"


# --------------------------------------------------------------------------- memory
#: A sweep of this size took the machine down once: twelve pooled workers kept every map they had
#: ever drawn and grew to 12 GB each, 106 GB in all, and the OOM killer took the editor with them
#: (2026-09-25). So the sweep now schedules against MEASURED memory: every chunk of work runs in a
#: fresh process that reports its own peak RSS, the peak of a strategy becomes its estimate, and a
#: chunk starts only when its estimate fits under both budgets below.
JOB_GB = 55.0          # the sweep's own processes, summed estimates
SYSTEM_GB = 75.0       # the whole machine, measured, plus what running chunks may still grow by
EMERGENCY_GB = 92.0    # above this the youngest chunk is stopped and its work requeued
PROBE_GB = 8.0         # estimate for a strategy nobody has measured yet
CHUNK = 6              # configurations per process: enough to amortize imports, few enough to
                       # bound what a process accumulates in its caches
TIMEOUT = 2400.0       # seconds one configuration may take before its process is stopped


def used_gb() -> float:
    """System memory in use (MemTotal - MemAvailable), GiB."""
    info = {}
    with open("/proc/meminfo") as fh:
        for line in fh:
            key, value = line.split(":", 1)
            info[key] = int(value.split()[0])
    return (info["MemTotal"] - info["MemAvailable"]) / 1024 ** 2


def rss_gb(pid: int) -> float:
    try:
        with open(f"/proc/{pid}/status") as fh:
            for line in fh:
                if line.startswith("VmRSS"):
                    return int(line.split()[1]) / 1024 ** 2
    except OSError:
        pass
    return 0.0


def _serve(chunk: list, queue) -> None:
    """A worker process: run each configuration of one chunk, report each with its peak RSS."""
    import resource
    _init_worker()
    for cfg in chunk:
        rec = run_one(cfg)
        rec["peak_rss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                                   / 1024 ** 2, 3)
        rec["worker_pid"] = os.getpid()
        queue.put(rec)


def plan_chunks(todo: list, chunk: int = CHUNK) -> list:
    """(group, [cfgs]) chunks; the first chunk of every group is a single probe configuration."""
    groups = {}
    for cfg in todo:
        groups.setdefault((cfg["organism"], cfg["strategy"]), []).append(cfg)
    out = []
    for group, cfgs in groups.items():
        out.append((group, cfgs[:1]))
        rest = cfgs[1:]
        out += [(group, rest[i:i + chunk]) for i in range(0, len(rest), chunk)]
    return out


def sweep(out_dir: str, workers: int, organisms, only=None, log=print, job_gb: float = JOB_GB,
          system_gb: float = SYSTEM_GB, chunk: int = CHUNK) -> int:
    """Run every configuration not yet in runs.jsonl, appending as each finishes.

    Memory-scheduled rather than pooled: see JOB_GB. Probes come first (one configuration of every
    strategy), so every estimate is measured before the bulk of the work is placed against it.
    """
    import collections
    import multiprocessing as mp
    import queue as _queue
    from starplast import strategies as S
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "runs.jsonl")
    done, est = set(), {}
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                rec = json.loads(line)
                done.add(rec["id"])
                if rec.get("peak_rss_gb"):
                    g = (rec["organism"], rec["strategy"])
                    est[g] = max(est.get(g, 0.0), 1.2 * rec["peak_rss_gb"] + 0.5)
    todo = []
    for org in organisms:
        ctx = S.Context.shipped(org)
        todo += [c for c in configurations(org, ctx, only) if config_id(c) not in done]
    chunks = plan_chunks(todo, chunk)
    # Probes first, then the rest; within each, the slowest strategies first so the long tail
    # does not land at the end on one process.
    slow = {"multiplex_modules": 0, "masked_imputation": 1, "holdout_search": 2,
            "geneset_hunt": 3, "triangulation": 4, "understudied_first": 5,
            "positive_unlabeled": 6, "consensus_modules": 7, "blind_battery": 8}
    chunks.sort(key=lambda c: (c[0] in est, len(c[1]) > 1, slow.get(c[0][1], 99)))
    pending = collections.deque(chunks)
    attempts = collections.Counter()
    log(f"{len(done):,} runs already done; {len(todo):,} to run in {len(chunks):,} chunks, "
        f"at most {workers} processes, job budget {job_gb:.0f} GB, system limit "
        f"{system_gb:.0f} GB")
    mpx = mp.get_context("spawn")
    results = mpx.Queue()
    running = {}                  # pid -> dict(proc, group, left, est, started, last)
    t0, n_done = time.monotonic(), 0

    def estimate(group):
        return est.get(group, PROBE_GB)

    def requeue(info, why):
        group, left = info["group"], info["left"]
        if not left:
            return
        attempts[group] += 1
        est[group] = max(estimate(group) * 1.5, info["est"] * 1.5)
        if attempts[group] > 2:
            with open(path, "a") as fh:
                for cfg in left:
                    rec = {"id": config_id(cfg), **cfg, "verdict": "NOT RUN",
                           "note": f"worker stopped ({why}) three times", "wall_seconds": 0}
                    fh.write(json.dumps(rec, default=str) + "\n")
            log(f"  gave up on {len(left)} configurations of {group}: {why}")
            return
        for i in range(0, len(left), max(1, chunk // 2)):
            pending.appendleft((group, left[i:i + max(1, chunk // 2)]))
        log(f"  requeued {len(left)} configurations of {group} ({why}); estimate now "
            f"{est[group]:.1f} GB")

    with open(path, "a") as fh:
        while pending or running:
            # 1. results
            while True:
                try:
                    rec = results.get_nowait()
                except _queue.Empty:
                    break
                fh.write(json.dumps(rec, default=str) + "\n")
                fh.flush()
                n_done += 1
                g = (rec["organism"], rec["strategy"])
                if rec.get("peak_rss_gb"):
                    est[g] = max(est.get(g, 0.0), 1.2 * rec["peak_rss_gb"] + 0.5)
                info = running.get(rec.get("worker_pid"))
                if info is not None:
                    info["left"] = [c for c in info["left"] if config_id(c) != rec["id"]]
                    info["last"] = time.monotonic()
                if n_done % 25 == 0:
                    rate = n_done / max(time.monotonic() - t0, 1e-9)
                    left = sum(len(c[1]) for c in pending) + sum(len(i["left"])
                                                                   for i in running.values())
                    log(f"{n_done:,} runs this session, {rate * 3600:.0f}/h, {left:,} left "
                        f"(~{left / max(rate, 1e-9) / 3600:.1f} h); {len(running)} processes, "
                        f"system {used_gb():.1f} GB")
            # 2. finished, stuck or killed processes
            now = time.monotonic()
            for pid, info in list(running.items()):
                proc = info["proc"]
                if proc.exitcode is None and now - info["last"] > TIMEOUT:
                    proc.kill()
                    proc.join(5)
                    running.pop(pid)
                    requeue(info, "timeout")
                    continue
                if proc.exitcode is not None:
                    running.pop(pid)
                    proc.join(1)
                    # Results can still be in flight; give them a moment before calling it lost.
                    time.sleep(0.2)
                    while True:
                        try:
                            rec = results.get_nowait()
                        except _queue.Empty:
                            break
                        fh.write(json.dumps(rec, default=str) + "\n")
                        n_done += 1
                        info["left"] = [c for c in info["left"] if config_id(c) != rec["id"]]
                        g = (rec["organism"], rec["strategy"])
                        if rec.get("peak_rss_gb"):
                            est[g] = max(est.get(g, 0.0), 1.2 * rec["peak_rss_gb"] + 0.5)
                        other = running.get(rec.get("worker_pid"))
                        if other is not None:
                            other["left"] = [c for c in other["left"]
                                             if config_id(c) != rec["id"]]
                            other["last"] = time.monotonic()
                    fh.flush()
                    if info["left"]:
                        requeue(info, f"exit code {proc.exitcode}")
            # 3. emergency: the machine, not the budget
            sys_used = used_gb()
            if sys_used >= EMERGENCY_GB and running:
                pid = max(running, key=lambda p: running[p]["started"])
                info = running.pop(pid)
                info["proc"].kill()
                info["proc"].join(5)
                log(f"EMERGENCY: system at {sys_used:.1f} GB; stopped the youngest process")
                requeue(info, "memory emergency")
                time.sleep(3)
                continue
            # 4. launch what fits
            launched = True
            while pending and len(running) < workers and launched:
                launched = False
                committed = sum(i["est"] for i in running.values())
                growth = sum(max(0.0, i["est"] - rss_gb(p)) for p, i in running.items())
                sys_used = used_gb()
                for k in range(min(len(pending), 50)):
                    group, cfgs = pending[k]
                    need = estimate(group)
                    if (committed + need <= job_gb and sys_used + growth + need <= system_gb):
                        del pending[k]
                        proc = mpx.Process(target=_serve, args=(cfgs, results), daemon=True)
                        proc.start()
                        running[proc.pid] = {"proc": proc, "group": group, "left": list(cfgs),
                                             "est": need, "started": time.monotonic(),
                                             "last": time.monotonic()}
                        launched = True
                        break
            time.sleep(0.5)
    with open(os.path.join(out_dir, "memory_estimates.json"), "w") as fh:
        json.dump({f"{o}:{k}": round(v, 2) for (o, k), v in sorted(est.items())}, fh, indent=1)
    log(f"sweep finished: {n_done:,} runs this session in "
        f"{(time.monotonic() - t0) / 3600:.2f} h")
    return n_done


def _ci(values, level: float = 0.95) -> tuple:
    """Mean and a t interval across seeds; (mean, mean, mean) for a single value."""
    import numpy as np
    from scipy import stats
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], dtype=float)
    if not len(v):
        return float("nan"), float("nan"), float("nan")
    if len(v) == 1:
        return float(v[0]), float(v[0]), float(v[0])
    m, se = float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))
    h = se * stats.t.ppf(0.5 + level / 2, len(v) - 1)
    return m, m - h, m + h


def summarize(out_dir: str, log=print) -> dict:
    """runs.jsonl -> summary, best, defaults and sensitivity tables."""
    import numpy as np
    import pandas as pd
    from starplast import strategies as S
    rows = [json.loads(line) for line in open(os.path.join(out_dir, "runs.jsonl")) if line.strip()]
    df = pd.DataFrame(rows)
    df["settings_key"] = df["settings"].map(lambda d: json.dumps(d, sort_keys=True))
    df["skill"] = (df["observed"] - df["null_mean"]) / (1 - df["null_mean"]).where(
        (1 - df["null_mean"]).abs() > 1e-9)
    groups = []
    for (org, key, sk), g in df.groupby(["organism", "strategy", "settings_key"]):
        conclusive = g[g["verdict"].isin(["PASS", "FAIL"])]
        m, lo, hi = _ci(conclusive["observed"].tolist())
        sm, slo, shi = _ci(conclusive["skill"].tolist())
        groups.append({"organism": org, "strategy": key, "settings": sk,
                       "metric": g["metric"].dropna().iloc[0] if g["metric"].notna().any() else "",
                       "runs": len(g), "conclusive": len(conclusive),
                       "pass_rate": float((g["verdict"] == "PASS").mean()),
                       "observed": m, "observed_low": lo, "observed_high": hi,
                       "chance": float(conclusive["null_mean"].mean()) if len(conclusive) else
                       float("nan"),
                       "skill": sm, "skill_low": slo, "skill_high": shi,
                       "seconds": float(g["wall_seconds"].mean())})
    summary = pd.DataFrame(groups)
    summary.to_csv(os.path.join(out_dir, "summary.csv"), index=False)
    ok = summary[summary["conclusive"] >= 3]
    best = (ok.sort_values("skill_low", ascending=False).groupby(["organism", "strategy"]).head(1)
            .sort_values(["organism", "strategy"]))
    best.to_csv(os.path.join(out_dir, "best.csv"), index=False)
    defaults = []
    for org in summary["organism"].unique():
        ctx_cats = None
        for s in S.catalog():
            sub = summary[(summary["organism"] == org) & (summary["strategy"] == s.key)]
            if sub.empty:
                continue
            if ctx_cats is None:
                ctx_cats = S.Context.shipped(org)
            d = s.defaults(ctx_cats)
            match = sub[sub["settings"].map(
                lambda k: all(d.get(p) == v for p, v in json.loads(k).items()))]
            if len(match):
                defaults.append(match.iloc[0].to_dict())
    pd.DataFrame(defaults).to_csv(os.path.join(out_dir, "defaults.csv"), index=False)
    sens = []
    for (org, key), g in df.groupby(["organism", "strategy"]):
        params = sorted({p for d in g["settings"] for p in d})
        for p in params:
            vals = g["settings"].map(lambda d: d.get(p))
            for v, gg in g.groupby(vals.astype(str)):
                m, lo, hi = _ci(gg["skill"].tolist())
                sens.append({"organism": org, "strategy": key, "parameter": p, "value": v,
                             "runs": len(gg), "skill": m, "skill_low": lo, "skill_high": hi,
                             "pass_rate": float((gg["verdict"] == "PASS").mean())})
    pd.DataFrame(sens).to_csv(os.path.join(out_dir, "sensitivity.csv"), index=False)
    log(f"{len(df):,} runs -> {len(summary):,} configurations; wrote summary, best, defaults, "
        f"sensitivity to {out_dir}")
    return {"runs": len(df), "configurations": len(summary)}


README_START = "<!-- calibration:start -->"
README_END = "<!-- calibration:end -->"


def _default_settings(org: str) -> dict:
    """(org, strategy) -> the JSON of its default setting over the swept parameters."""
    from starplast import calibration as C
    from starplast import strategies as S
    ctx = S.Context.shipped(org)
    out = {}
    for s in S.catalog():
        grid = [p for p in GRID.get(s.key, {}) if p not in C.TARGET_PARAMS]
        if not grid:
            out[(org, s.key)] = json.dumps({})
            continue
        d = s.defaults(ctx)
        out[(org, s.key)] = {p: d.get(p) for p in grid}
    return out


def _match_default(runs, defaults: dict) -> dict:
    """The swept setting closest to each strategy's defaults (most parameters equal)."""
    out = {}
    for (org, key), g in runs.groupby(["organism", "strategy"]):
        d = defaults.get((org, key))
        if isinstance(d, str) or d is None:
            out[(org, key)] = d
            continue
        best, score = None, -1
        for sk in g["setting"].unique():
            cand = json.loads(sk)
            n = sum(str(cand.get(p)) == str(v) for p, v in d.items())
            if n > score:
                best, score = sk, n
        out[(org, key)] = best
    return out


def publish(out_dir: str, log=print, readme: str | None = None,
            doc: str | None = None) -> dict:
    """runs.jsonl -> starplast/data/strategy_calibration.json, the README table, the doc page."""
    import pandas as pd
    from starplast import calibration as C
    from starplast import strategies as S
    rows = [json.loads(l) for l in open(os.path.join(out_dir, "runs.jsonl")) if l.strip()]
    runs = C.runs_frame(rows)
    defaults = {}
    for org in sorted(runs["organism"].unique()):
        defaults.update(_default_settings(org))
    summary = C.summarise(runs, _match_default(runs, defaults))
    meta = {"date": _dt.date.today().isoformat(), "runs": int(len(runs)),
            "run_dir": os.path.relpath(out_dir, ROOT), "seeds": list(SEEDS),
            "targets": {o: list(t) for o, t in TARGETS.items()},
            "numbers": {o: list(t) for o, t in NUMBERS.items()},
            "method": "starplast.calibration: skill = (observed - chance)/(1 - chance); "
                      "tuned setting chosen on seeds 1-3, reported on seeds 4-5; 95% CI by "
                      "two-stage bootstrap over targets then runs"}
    path = C.write(summary, meta=meta)
    log(f"wrote {path}")
    readme = readme or os.path.join(ROOT, "README.md")
    text = open(readme).read()
    block = [README_START, "",
             f"Measured {meta['date']} over {meta['runs']:,} self-tests: every strategy's "
             f"self-test run over a grid of its settings, several held-out known labels and five "
             f"seeds. **Skill** puts every metric on one scale -- 0 is the same procedure on "
             f"shuffled data, 1 is perfect. *Tuned* is the best setting chosen on seeds 1-3 and "
             f"reported on seeds 4-5, so it is not the luckiest of many settings. Intervals are 95%, "
             f"by bootstrap over held-out targets and runs. "
             f"[Every number, per target and per setting](docs/calibration.md).", ""]
    for org, name in (("Tg", "*T. gondii*"), ("Pf", "*P. falciparum*")):
        if org in summary:
            counts = pd.Series([e["grade"] for e in summary[org].values()]).value_counts()
            block += [f"**{name}** -- " + ", ".join(f"{n} {g}" for g, n in counts.items()), "",
                      C.markdown_table(org)]
    block.append(README_END)
    if README_START in text:
        a, rest = text.split(README_START, 1)
        text = a + "\n".join(block) + rest.split(README_END, 1)[1]
    else:
        anchor = "[All 32 strategies, with their tests and measured verdicts]"
        i = text.index(anchor)
        text = text[:i] + "\n".join(block) + "\n\n" + text[i:]
    with open(readme, "w") as fh:
        fh.write(text)
    log(f"updated {readme}")
    doc = doc or os.path.join(ROOT, "docs", "calibration.md")
    with open(doc, "w") as fh:
        fh.write(calibration_doc(summary, meta))
    log(f"wrote {doc}")
    return summary


def calibration_doc(summary: dict, meta: dict) -> str:
    """The long page: per strategy, what was hidden, every target, every setting."""
    from starplast import strategies as S

    def cell(c):
        if not c or c.get("skill") is None:
            return "--"
        return (f"{c['skill']:.3f} [{c['skill_low']:.3f}, {c['skill_high']:.3f}]")
    out = ["# Strategy calibration", "",
           f"Generated by `scripts/calibrate_strategies.py --publish` on {meta['date']} from "
           f"{meta['runs']:,} self-test runs (`{meta['run_dir']}`).", "",
           "Every strategy carries a self-test that hides information already known -- labels, "
           "set members, edges or values -- asks the strategy for it back, and compares the answer "
           "with the same procedure on shuffled data. Calibration runs that test over a grid of "
           "the strategy's settings, over several held-out labels and five seeds, so each number "
           "below is a mean with an interval rather than one draw.", "",
           "* **Skill** = (observed - chance) / (1 - chance): 0 is the shuffled-data null, 1 is "
           "perfect. Negative means worse than shuffled.",
           "* **At defaults**: the strategy as it opens in the application.",
           "* **Tuned**: the setting with the best lower confidence bound on seeds 1-3, "
           "reported on seeds 4-5 only.",
           "* **95% CI**: two-stage bootstrap (held-out targets, then runs within a target).",
           "* **Pass**: share of conclusive runs beating the null's 95th percentile by the "
           "strategy's minimum effect.", "",
           "Grades: *reliable* (defaults beat chance with the interval above 0.05 and pass at "
           "least 60%), *works when tuned* (only the tuned setting does), *weak* (above chance "
           "on average but not reliably), *no skill*, *untestable* (fewer than five conclusive "
           "runs).", ""]
    for org, name in (("Tg", "Toxoplasma gondii"), ("Pf", "Plasmodium falciparum")):
        entries = summary.get(org) or {}
        if not entries:
            continue
        out += [f"## {name}", ""]
        for s in S.catalog():
            e = entries.get(s.key)
            if not e:
                continue
            d, t = e.get("default") or {}, e.get("tuned") or {}
            out += [f"### {s.number:02d} · {s.title} -- {e['grade']}", "",
                    f"{s.test_description}", "",
                    f"Metric: `{e['metric']}`. {e['runs']:,} runs, {e['settings_tested']} "
                    f"settings, targets: {', '.join(e.get('targets') or []) or '(the strategy '
                    f'chooses its own)'}.", "",
                    "| | setting | skill [95% CI] | pass [95% CI] | observed | chance | runs |",
                    "|---|---|---|---|---|---|---|"]
            for label, c in (("at defaults", d), ("tuned", t)):
                if not c:
                    continue
                setting = ", ".join(f"{k}={v}" for k, v in (c.get("setting") or {}).items())
                pr = c.get("pass_rate")
                out.append(f"| {label} | {setting or '--'} | {cell(c)} | "
                           + (f"{100 * pr:.0f}% [{100 * c['pass_low']:.0f}, "
                              f"{100 * c['pass_high']:.0f}]" if pr is not None else "--")
                           + f" | {c.get('observed') if c.get('observed') is None else round(c['observed'], 3)}"
                           f" | {c.get('chance') if c.get('chance') is None else round(c['chance'], 3)}"
                           f" | {c.get('conclusive', 0)} |")
            per = e.get("per_target") or {}
            if len(per) > 1:
                out += ["", "Tuned setting, per held-out target:", "",
                        "| target | skill [95% CI] | pass | runs |", "|---|---|---|---|"]
                for tgt, c in per.items():
                    pr = c.get("pass_rate")
                    out.append(f"| {tgt} | {cell(c)} | "
                               f"{'--' if pr is None else f'{100 * pr:.0f}%'} | "
                               f"{c.get('conclusive', 0)} |")
            sens = e.get("sensitivity") or {}
            if sens:
                out += ["", "Sensitivity (mean skill at each value, all other settings pooled):",
                        ""]
                for p, vals in sens.items():
                    parts = [f"{v}: {c['skill']:.2f}" for v, c in vals.items()
                             if c.get("skill") is not None]
                    out.append(f"* `{p}` -- " + "; ".join(parts))
            out.append("")
    return "\n".join(out) + "\n"



def main(argv=None, log=print) -> int:
    """Sweep, or summarize what a sweep has written."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=16,
                    help="most processes at once; memory decides how many actually run")
    ap.add_argument("--job-gb", type=float, default=JOB_GB)
    ap.add_argument("--system-gb", type=float, default=SYSTEM_GB)
    ap.add_argument("--organism", choices=("Tg", "Pf"), action="append")
    ap.add_argument("--only", default="", help="comma-separated strategy keys")
    ap.add_argument("--out", default=os.path.join(
        ROOT, "results", f"calibration_{_dt.date.today().isoformat()}"))
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--publish", action="store_true",
                    help="write the shipped calibration, the README table and docs/calibration.md")
    ap.add_argument("--count", action="store_true", help="only count the configurations")
    args = ap.parse_args(argv)
    only = {k for k in args.only.split(",") if k} or None
    if args.summarize:
        summarize(args.out, log=log)
        return 0
    if args.publish:
        publish(args.out, log=log)
        return 0
    orgs = args.organism or ["Tg", "Pf"]
    if args.count:
        from starplast import strategies as S
        for org in orgs:
            cfgs = configurations(org, S.Context.shipped(org), only)
            log(f"{org}: {len(cfgs):,} runs ({len(cfgs) // len(SEEDS):,} configurations x "
                f"{len(SEEDS)} seeds)")
        return 0
    sweep(args.out, args.workers, orgs, only, log=log, job_gb=args.job_gb,
          system_gb=args.system_gb)
    summarize(args.out, log=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
