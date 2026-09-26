"""How good each strategy is, measured: the calibration sweep turned into verdicts.

`scripts/calibrate_strategies.py` runs every strategy's self-test over a grid of its settings, over
several held-out labels ("targets") and five seeds. This module turns those runs into what a person
choosing a strategy needs, and what the Strategies tab and the README show:

    skill        (observed - chance) / (1 - chance): 0 is the shuffled-data null, 1 is perfect,
                 so a correct-call rate, an AUROC, an F1 and a correlation read on one scale.
    at defaults  the strategy as it opens, averaged over every target and seed.
    tuned        the best setting -- CHOSEN on seeds 1-3 and REPORTED on seeds 4-5. Choosing and
                 reporting on the same runs would report the luckiest of many settings; splitting
                 by seed keeps the reported number an honest estimate of what the choice buys.
    95% CI       a two-stage bootstrap: targets resampled, then runs within each target. Runs on
                 one target share its labels and are not independent; treating them as such gives
                 intervals several times too narrow.
    pass rate    the share of conclusive runs that beat the null's 95th percentile by the
                 strategy's minimum effect, with a Wilson interval.

A grade summarises the lot (`grade`): *reliable*, *works when tuned*, *weak*, *no skill*, or
*untestable* when too few runs were conclusive to say.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(__file__), "data")
CALIBRATION = os.path.join(DATA, "strategy_calibration.json")

#: Settings that choose WHAT is held out rather than HOW the strategy works. A setting is the rest.
TARGET_PARAMS = ("target", "exclude")
#: Seeds used to choose the best setting; the others report it.
CHOOSE_SEEDS = (1, 2, 3)
REPORT_SEEDS = (4, 5)
GRADES = ("reliable", "works when tuned", "weak", "no skill", "untestable")

_CACHE: dict = {}


# --------------------------------------------------------------------------- statistics
def wilson(k: int, n: int, z: float = 1.96) -> tuple:
    """Wilson score interval for k successes in n."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def cluster_bootstrap(values, clusters, n_boot: int = 2000, seed: int = 0) -> tuple:
    """(mean, low, high): a two-stage bootstrap -- clusters, then values within each cluster."""
    v = np.asarray(values, dtype=float)
    c = np.asarray(clusters)
    ok = np.isfinite(v)
    v, c = v[ok], c[ok]
    if not len(v):
        return (float("nan"),) * 3
    groups = [v[c == g] for g in pd.unique(c)]
    mean = float(v.mean())
    if len(v) == 1:
        return (mean, mean, mean)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for b in range(n_boot):
        picks = rng.integers(0, len(groups), len(groups))
        sample = [groups[i][rng.integers(0, len(groups[i]), len(groups[i]))] for i in picks]
        boots[b] = np.concatenate(sample).mean()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return (mean, float(lo), float(hi))


# --------------------------------------------------------------------------- from runs to verdicts
def runs_frame(rows) -> pd.DataFrame:
    """runs.jsonl records -> one row per run, with skill, target and setting split apart."""
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df[df["verdict"].notna()].copy()
    settings = df["settings"].map(lambda d: d if isinstance(d, dict) else {})
    df["target"] = settings.map(lambda d: next((str(d[p]) for p in TARGET_PARAMS if p in d), ""))
    df["setting"] = settings.map(lambda d: json.dumps(
        {k: v for k, v in sorted(d.items()) if k not in TARGET_PARAMS}, default=str))
    obs = pd.to_numeric(df.get("observed"), errors="coerce")
    chance = pd.to_numeric(df.get("null_mean"), errors="coerce")
    room = (1 - chance).where((1 - chance).abs() > 1e-9)
    df["skill"] = (obs - chance) / room
    df["conclusive"] = df["verdict"].isin(["PASS", "FAIL"])
    return df


def _cell(g: pd.DataFrame) -> dict:
    """Skill, pass rate and their intervals for one group of runs."""
    c = g[g["conclusive"]]
    mean, lo, hi = cluster_bootstrap(c["skill"], c["target"]) if len(c) else (float("nan"),) * 3
    k, n = int((c["verdict"] == "PASS").sum()), len(c)
    plo, phi = wilson(k, n)
    obs = pd.to_numeric(c.get("observed"), errors="coerce")
    chance = pd.to_numeric(c.get("null_mean"), errors="coerce")
    return {"runs": int(len(g)), "conclusive": n, "skill": mean, "skill_low": lo,
            "skill_high": hi, "pass_rate": k / n if n else float("nan"), "pass_low": plo,
            "pass_high": phi, "observed": float(obs.mean()) if n else float("nan"),
            "chance": float(chance.mean()) if n else float("nan"),
            "targets": sorted(set(c["target"]) - {""}),
            "not_run": int((g["verdict"] == "NOT RUN").sum())}


def grade(entry: dict) -> str:
    """One word for how far to trust a strategy, from its default and tuned calibration."""
    d, t = entry.get("default") or {}, entry.get("tuned") or {}
    if max(d.get("conclusive", 0), t.get("conclusive", 0)) < 5:
        return "untestable"
    if d.get("skill_low", -1) > 0.05 and d.get("pass_rate", 0) >= 0.6:
        return "reliable"
    if t.get("skill_low", -1) > 0.05 and t.get("pass_rate", 0) >= 0.5:
        return "works when tuned"
    if max(d.get("skill", 0), t.get("skill", 0)) > 0.02:
        return "weak"
    return "no skill"


def summarise(runs: pd.DataFrame, defaults: dict | None = None) -> dict:
    """{organism: {strategy: entry}} -- the calibration the application ships.

    `defaults` maps (organism, strategy) -> the default setting's JSON string, so "at defaults"
    means the strategy as it opens; without it the most common setting stands in.
    """
    out: dict = {}
    if runs.empty:
        return out
    for (org, key), g in runs.groupby(["organism", "strategy"]):
        entry = {"metric": str(g["metric"].dropna().iloc[0]) if g["metric"].notna().any() else "",
                 "runs": int(len(g)), "settings_tested": int(g["setting"].nunique()),
                 "targets": sorted(set(g["target"]) - {""})}
        dkey = (defaults or {}).get((org, key))
        if dkey is None or dkey not in set(g["setting"]):
            dkey = g["setting"].value_counts().index[0]
        entry["default"] = {"setting": json.loads(dkey), **_cell(g[g["setting"] == dkey])}
        # Choose on some seeds, report on the others.
        choose = g[g["seed"].isin(CHOOSE_SEEDS) & g["conclusive"]]
        best_key, best_low = None, -np.inf
        for sk, sg in choose.groupby("setting"):
            if len(sg) < 3:
                continue
            _m, low, _h = cluster_bootstrap(sg["skill"], sg["target"], n_boot=500)
            if low > best_low:
                best_key, best_low = sk, low
        if best_key is not None:
            report = g[(g["setting"] == best_key) & g["seed"].isin(REPORT_SEEDS)]
            entry["tuned"] = {"setting": json.loads(best_key), **_cell(report),
                              "chosen_on_seeds": list(CHOOSE_SEEDS),
                              "reported_on_seeds": list(REPORT_SEEDS)}
            per_target = {}
            for tgt, tg in g[g["setting"] == best_key].groupby("target"):
                cell = _cell(tg)
                per_target[tgt or "(none)"] = {k: cell[k] for k in
                                               ("skill", "skill_low", "skill_high", "pass_rate",
                                                "conclusive")}
            entry["per_target"] = per_target
        # Which way each setting pushes skill, over everything tested.
        sens = {}
        params = sorted({p for s in g["setting"] for p in json.loads(s)})
        for p in params:
            vals = g["setting"].map(lambda s: str(json.loads(s).get(p)))
            rows = {}
            for v, vg in g[g["conclusive"]].groupby(vals[g["conclusive"]]):
                m, lo, hi = cluster_bootstrap(vg["skill"], vg["target"], n_boot=300)
                rows[v] = {"skill": m, "skill_low": lo, "skill_high": hi, "runs": int(len(vg))}
            sens[p] = rows
        entry["sensitivity"] = sens
        entry["grade"] = grade(entry)
        out.setdefault(org, {})[key] = entry
    return out


def _clean(x):
    if isinstance(x, float):
        return None if not math.isfinite(x) else round(x, 4)
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.floating,)):
        return _clean(float(x))
    if isinstance(x, (np.integer,)):
        return int(x)
    return x


def write(calibration: dict, path: str = CALIBRATION, meta: dict | None = None) -> str:
    """Ship a calibration: the summary plus the sweep that produced it, with NaN written as null.

    The meta block is not decoration -- a calibration without the date, the run directory and the
    method that made it is a set of numbers nobody can reproduce or supersede.
    """
    with open(path, "w") as fh:
        json.dump(_clean({"meta": meta or {}, "organisms": calibration}), fh, indent=1)
    _CACHE.clear()
    return path


# --------------------------------------------------------------------------- reading it back
def load(path: str = CALIBRATION) -> dict:
    """The shipped calibration, {"meta": ..., "organisms": {org: {key: entry}}}; empty if none."""
    if path not in _CACHE:
        try:
            with open(path) as fh:
                _CACHE[path] = json.load(fh)
        except (OSError, ValueError):
            _CACHE[path] = {"meta": {}, "organisms": {}}
    return _CACHE[path]


def entry(key: str, organism: str = "Tg", path: str = CALIBRATION) -> dict | None:
    """What calibration measured for one strategy on one organism, or None if it was not measured."""
    return (load(path).get("organisms", {}).get(organism) or {}).get(key)


def tuned_settings(key: str, organism: str = "Tg", path: str = CALIBRATION) -> dict:
    """The setting calibration found best, without the held-out target (the user picks that)."""
    e = entry(key, organism, path) or {}
    return dict((e.get("tuned") or {}).get("setting") or {})


def setting_text(setting: dict) -> str:
    """A setting as text a reader cannot misparse: a grid value is bracketed, parameters are split
    by semicolons. `min_cluster_size=20, 50, n_neighbors=10, 30` read as four parameters."""
    if not setting:
        return "defaults"
    parts = []
    for k, v in setting.items():
        v = str(v)
        parts.append(f"{k}=[{v}]" if "," in v else f"{k}={v}")
    return "; ".join(parts)


def _fmt(cell: dict) -> str:
    if not cell or cell.get("skill") is None:
        return "not measured"
    return (f"skill {cell['skill']:.2f} [{cell['skill_low']:.2f}, {cell['skill_high']:.2f}], "
            f"passes {100 * (cell.get('pass_rate') or 0):.0f}% of {cell.get('conclusive', 0)}")


def sentence(key: str, organism: str = "Tg", path: str = CALIBRATION) -> str:
    """One paragraph of plain text: what calibration measured for this strategy."""
    e = entry(key, organism, path)
    if not e:
        return ""
    d, t = e.get("default") or {}, e.get("tuned") or {}
    text = (f"Calibrated: {e['grade'].upper()}. {e['runs']:,} held-out tests over "
            f"{e['settings_tested']} settings and {len(e.get('targets') or []) or 1} target(s). "
            f"At defaults: {_fmt(d)}.")
    if t:
        shown = setting_text(t.get("setting") or {})
        text += (f" Tuned ({shown}; chosen on seeds 1-3, reported on 4-5): {_fmt(t)}.")
    return text + " Skill: 0 = the shuffled-data null, 1 = perfect."


def overview(organism: str = "Tg", path: str = CALIBRATION) -> pd.DataFrame:
    """One row per calibrated strategy -- the README table as a frame."""
    from . import strategies as S
    rows = []
    for s in S.catalog():
        e = entry(s.key, organism, path)
        if not e:
            continue
        d, t = e.get("default") or {}, e.get("tuned") or {}
        rows.append({"number": s.number, "strategy": s.title, "key": s.key,
                     "metric": e.get("metric", ""), "grade": e.get("grade", ""),
                     "runs": e.get("runs"), "settings": e.get("settings_tested"),
                     "default_skill": d.get("skill"), "default_low": d.get("skill_low"),
                     "default_high": d.get("skill_high"), "default_pass": d.get("pass_rate"),
                     "tuned_skill": t.get("skill"), "tuned_low": t.get("skill_low"),
                     "tuned_high": t.get("skill_high"), "tuned_pass": t.get("pass_rate"),
                     "tuned_setting": json.dumps(t.get("setting") or {})})
    return pd.DataFrame(rows)


def markdown_table(organism: str = "Tg", path: str = CALIBRATION) -> str:
    """The README table: every strategy's measured skill on held-out known labels."""
    df = overview(organism, path)
    if df.empty:
        return "_No calibration shipped yet._\n"

    def ci(m, lo, hi):
        if m is None or (isinstance(m, float) and not math.isfinite(m)):
            return "--"
        return f"{m:.2f} [{lo:.2f}, {hi:.2f}]"

    def pct(p):
        return "--" if p is None or (isinstance(p, float) and not math.isfinite(p)) else \
            f"{100 * p:.0f}%"
    lines = ["| # | Strategy | Metric | Grade | Skill at defaults [95% CI] | Pass | "
             "Skill tuned [95% CI] | Pass | Tuned setting | Tests |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in df.itertuples():
        setting = setting_text(json.loads(r.tuned_setting)) if r.tuned_setting != "{}" else "--"
        lines.append(f"| {r.number:02d} | {r.strategy} | {r.metric} | {r.grade} | "
                     f"{ci(r.default_skill, r.default_low, r.default_high)} | "
                     f"{pct(r.default_pass)} | {ci(r.tuned_skill, r.tuned_low, r.tuned_high)} | "
                     f"{pct(r.tuned_pass)} | {setting} | {r.runs:,} |")
    return "\n".join(lines) + "\n"
