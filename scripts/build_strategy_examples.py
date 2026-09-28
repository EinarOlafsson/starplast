"""Worked examples for the strategy cards: one real failure and one real success per strategy.

Every strategy card in the Strategies tab shows two small scorecards side by side, "Fails when..."
and "Works when...". Both are real runs, chosen from the calibration sweep already on disk
(`results/calibration_2026-09-26b/runs.jsonl`, 7,640 self-tests with their settings, targets,
observed values, chance levels and scorecards). Nothing is re-measured and no calibration number
changes; this script only chooses, explains and, for the success, looks once at what the strategy
says about genes nobody has labeled.

How each example is chosen, per organism and strategy:

* **Success** -- the PASS run at the strategy's default setting with the highest skill. Where no run
  at the default passed, the best PASS at any setting, and the card says so. Where no run passed at
  all, there is no success to show and the card says why (the runs' own note).
* **Failure** -- a FAIL run on real data, preferring the default setting and, within it, the lowest
  skill: the clearest case of the strategy having nothing to say. Its failure MODE is written from
  the numbers (the label is not encoded, too few genes reachable, a signal real but below the
  required margin, sets too wide, ...). Where a strategy never failed on real data, the failure is
  its self-test on the noise table, `strategies.planted_context(null=True)`, where every label and
  edge was dealt out at random, and the card says that is what it is.
* **New information** -- the success's setting is run once (`Strategy.run`) on the shipped table,
  and its five best-supported calls, ranked genes or findings for genes WITHOUT a known label are
  kept with their support. A run is stopped after `RUN_BUDGET` seconds and recorded as skipped.

Output: `starplast/data/strategy_examples.json`, and the executed notebook
`notebooks/strategy_examples_2026_09_28.ipynb` that shows every choice.

Run:  python scripts/build_strategy_examples.py                 # both organisms, every strategy
      python scripts/build_strategy_examples.py --only feature_knn --organism Tg --no-notebook
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import json
import math
import os
import re
import signal
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RUNS = os.path.join(ROOT, "results", "calibration_2026-09-26b", "runs.jsonl")
OUT = os.path.join(ROOT, "starplast", "data", "strategy_examples.json")
NOTEBOOK = os.path.join(ROOT, "notebooks", "strategy_examples_2026_09_28.ipynb")
def _calibrated_organisms() -> tuple:
    """The registered spaces the shipped calibration covers, in registry order."""
    from starplast import calibration as CAL
    from starplast import organisms
    have = set(CAL.load().get("organisms") or {})
    return tuple(c for c in organisms.codes() if c in have)


ORGANISMS = _calibrated_organisms()
#: Seconds one `run` may take before it is given up and recorded as skipped.
RUN_BUDGET = 300
#: How many new calls or findings each success keeps.
TOP = 5
#: The scorecard metrics an example keeps: its task's four headline bars and what their chance
#: ticks are computed from.
KEEP_EXTRA = ("fold_enrichment", "null_rate", "prevalence")


# --------------------------------------------------------------------------- choosing
def load_runs(path: str = RUNS) -> pd.DataFrame:
    """The sweep's runs as `calibration.runs_frame` reads them: skill, target and setting split."""
    from starplast import calibration as CAL
    with open(path) as fh:
        rows = [json.loads(line) for line in fh if line.strip()]
    return CAL.runs_frame(rows)


def default_key(key: str, organism: str, runs: pd.DataFrame) -> str | None:
    """The default setting's JSON string, as the calibration recorded it."""
    from starplast import calibration as CAL
    e = CAL.entry(key, organism) or {}
    setting = (e.get("default") or {}).get("setting")
    if setting is None:
        return None
    return json.dumps({k: v for k, v in sorted(setting.items())}, default=str)


def default_target(key: str, organism: str) -> str:
    """The label or measurement the strategy opens on for this organism ("" if it takes none)."""
    from starplast import strategies as S
    s = S.get(key)
    ctx = _light_context(organism)
    if any(p.name == "genes" for p in s.params):
        return str(S.default_category(ctx) or "")      # a list drawn from the default label
    for p in s.params:
        if p.name in ("target", "exclude"):
            try:
                return str(p.resolve(ctx) or "")
            except Exception:
                return ""
    return ""


_LIGHT: dict = {}


def _light_context(organism: str):
    from starplast import strategies as S
    if organism not in _LIGHT:
        _LIGHT[organism] = S.Context.shipped(organism)
    return _LIGHT[organism]


def choose(runs: pd.DataFrame, key: str, organism: str) -> dict:
    """{"success": row or None, "failure": row or None, "note": ...} for one strategy.

    Success: among PASS runs at the default setting, those on the strategy's own default target
    when any passed there (the question the card opens on), else every target; the highest skill.
    Failure: among FAIL runs at the default setting (else at any), the target that failed most
    often -- a consistent failure, not a fluke -- and within it the run of median skill.
    """
    g = runs[(runs["strategy"] == key) & (runs["organism"] == organism)]
    dkey = default_key(key, organism, runs)
    target = default_target(key, organism)
    out = {"runs": int(len(g)), "passes": int((g["verdict"] == "PASS").sum()),
           "fails": int((g["verdict"] == "FAIL").sum()), "default": dkey}
    ok = g[(g["verdict"] == "PASS") & g["skill"].notna()]
    at = ok[ok["setting"] == dkey]
    pick = at if len(at) else ok
    if target:
        own = pick[(pick["target"] == target)
                   | pick["hidden"].fillna("").str.endswith(f"genes of {target}")
                   | pick["hidden"].fillna("").str.contains(f"genes of {target} ", regex=False)]
    else:
        own = pick.iloc[0:0]
    pick = own if len(own) else pick
    out["success"] = (pick.sort_values(["skill", "seed"], ascending=[False, True]).iloc[0]
                      if len(pick) else None)
    out["success_at_default"] = bool(len(at))
    bad = g[(g["verdict"] == "FAIL") & g["skill"].notna()]
    at = bad[bad["setting"] == dkey]
    pick = at if len(at) else bad
    if len(pick):
        counts = pick.groupby("target").agg(n=("skill", "size"), mean=("skill", "mean"))
        worst = counts.sort_values(["n", "mean"], ascending=[False, True]).index[0]
        cand = pick[pick["target"] == worst].sort_values(["skill", "seed"])
        out["failure"] = cand.iloc[(len(cand) - 1) // 2]
    else:
        out["failure"] = None
    out["failure_at_default"] = bool(len(at))
    notes = g.loc[g["verdict"].isin(["INCONCLUSIVE", "NOT RUN"]), "note"].dropna()
    out["note"] = str(notes.value_counts().index[0]) if len(notes) else ""
    return out


# --------------------------------------------------------------------------- one example
def _f(x) -> float:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return x if math.isfinite(x) else float("nan")


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return None if not math.isfinite(float(x)) else round(float(x), 4)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def _card(task: str, scorecard: dict) -> dict:
    from starplast import scorecard as SC
    keys = {h.key for h in SC.headline(task)} | {h.chance_key for h in SC.headline(task)
                                                  if h.chance_key} | set(KEEP_EXTRA)
    return {k: v for k, v in (scorecard or {}).items() if k in keys}


def example_from_run(row, task: str) -> dict:
    """The fields a card needs from one runs.jsonl record."""
    settings = row["settings"] if isinstance(row["settings"], dict) else {}
    target = row.get("target") or ""
    m = re.search(r"'(.+?)' genes of (\S+)", row.get("hidden") or "")
    if m:                                       # a gene-list test: the list is the question
        target = f"{m.group(1)} ({m.group(2)})"
    for name in ("layer", "map_from", "condition", "column"):
        if not target and settings.get(name):  # the question is a layer, family or column
            target = str(settings[name])
    return {"source": "calibration", "run_id": row["id"], "seed": int(row["seed"]),
            "target": target, "settings": settings,
            "verdict": row["verdict"], "metric": row.get("metric") or "",
            "observed": _f(row.get("observed")), "chance": _f(row.get("null_mean")),
            "bar": _f(row.get("null_high")), "effect": _f(row.get("effect")),
            "min_effect": _f(row.get("min_effect")), "p_value": _f(row.get("p_value")),
            "n_hidden": int(_f(row.get("n_hidden")) or 0) if math.isfinite(
                _f(row.get("n_hidden"))) else 0,
            "hidden": row.get("hidden") or "", "null_kind": row.get("null_kind") or "",
            "skill": _f(row.get("skill")), "card": _card(task, row.get("scorecard"))}


def example_from_test(test, task: str, settings: dict) -> dict:
    """The same fields from a `TestResult` (the noise-table self-test)."""
    d = test.to_dict()
    return {"source": "noise table", "run_id": "", "seed": 42, "target": "",
            "settings": {k: v for k, v in settings.items() if not isinstance(v, str) or len(v) < 80},
            "verdict": d.get("verdict"), "metric": d.get("metric") or "",
            "observed": _f(d.get("observed")), "chance": _f(d.get("null_mean")),
            "bar": _f(d.get("null_high")), "effect": _f(d.get("effect")),
            "min_effect": _f(d.get("min_effect")), "p_value": _f(d.get("p_value")),
            "n_hidden": int(d.get("n_hidden") or 0), "hidden": d.get("hidden") or "",
            "null_kind": d.get("null_kind") or "", "skill": _f(d.get("skill")),
            "note": d.get("note") or "", "card": _card(task, d.get("scorecard"))}


def _num(x, digits: int = 2) -> str:
    x = _f(x)
    return "--" if not math.isfinite(x) else f"{x:.{digits}f}"


def failure_mode(ex: dict, task: str) -> str:
    """Why a failure failed, in one or two sentences grounded in its own numbers."""
    obs, chance, bar = ex["observed"], ex["chance"], ex["bar"]
    skill, eff, need = ex["skill"], ex["effect"], ex["min_effect"]
    card = ex["card"]
    target = ex.get("target") or ""
    what = f"'{target}'" if target else "the hidden information"
    if ex["source"] == "noise table":
        head = ("On the noise table every label and edge is dealt out at random, so there is "
                "nothing to find")
        if math.isfinite(obs):
            head += (f": {ex['metric']} was {_num(obs)} against {_num(chance)} on shuffled data "
                     f"(skill {_num(skill)})")
        tail = {"FAIL": ". The test calls that a FAIL, which is the failure it exists to catch.",
                "INCONCLUSIVE": ". It made too little to score, and the test declined to judge"
                                + (f" ({ex['note']})." if ex.get("note") else "."),
                "PASS": ". It still passed: a false alarm on random data, so read its passes "
                        "with care."}.get(ex["verdict"], ".")
        return head + tail
    reasons = []
    cov = _f(card.get("coverage", card.get("value_coverage")))
    if task == "label calls" and math.isfinite(cov) and cov < 0.5:
        reasons.append(f"it could reach only {cov:.0%} of the hidden genes, so most were never "
                       f"called")
    noise = _f(card.get("noise_share"))
    if task == "cluster recovery" and math.isfinite(noise) and noise > 0.4:
        reasons.append(f"{noise:.0%} of the hidden genes were left in no cluster at all")
    if task == "set retrieval":
        ret, prec = _f(card.get("returned")), _f(card.get("precision"))
        if math.isfinite(ret) and math.isfinite(prec) and prec < 0.2:
            reasons.append(f"the set was too wide: {ret:,.0f} genes returned, only {prec:.0%} of "
                           f"them members")
    if ex["n_hidden"] and ex["n_hidden"] < 30:
        reasons.append(f"only {ex['n_hidden']} hidden items could be scored")
    if (not math.isfinite(skill) or skill <= 0.02) and math.isfinite(chance) and chance >= 0.8 \
            and task == "label calls":
        core = (f"one class dominates {what}, so guessing it on shuffled data already scores "
                f"{_num(chance, 3)}; the strategy's {_num(obs, 3)} is no better than that (skill "
                f"{_num(skill)}), so what it reads does not separate the classes")
    elif not math.isfinite(skill) or skill <= 0.02:
        core = (f"{what} is not encoded in what this strategy reads: {ex['metric']} was "
                f"{_num(obs, 3)} against {_num(chance, 3)} on shuffled data (skill {_num(skill)})")
    elif math.isfinite(bar) and obs <= bar:
        core = (f"the signal is too weak to tell from luck: {_num(obs, 3)} did not clear "
                f"{_num(bar, 3)}, what shuffled data reaches one time in twenty (skill {_num(skill)})")
    elif math.isfinite(need) and math.isfinite(eff) and eff < need:
        core = (f"the signal is real but small: {_num(obs, 3)} beat shuffled data ({_num(chance, 3)}), "
                f"but by {_num(eff, 3)}, short of the {_num(need, 3)} margin a PASS requires")
    else:
        core = (f"{ex['metric']} was {_num(obs, 3)} against {_num(chance, 3)} on shuffled data "
                f"(skill {_num(skill)})")
    text = f"{core}."
    if reasons:
        text += " Also, " + "; ".join(reasons) + "."
    return text[0].upper() + text[1:]


def success_text(ex: dict, new: dict) -> str:
    """What the success showed, and what the run then added."""
    text = (f"{ex['metric'][0].upper() + ex['metric'][1:]} reached {_num(ex['observed'])} "
            f"against {_num(ex['chance'])} on shuffled data (skill {_num(ex['skill'])}, "
            f"{ex['n_hidden']:,} scored).")
    if new.get("skipped"):
        return text + f" Not re-run here: {new['skipped']}."
    if new.get("rows"):
        return text + f" Run once on every gene, it {new['what']}; the top {len(new['rows'])} are listed."
    return text + f" Run once on every gene, it {new.get('what') or 'added nothing new'}."


# --------------------------------------------------------------------------- new information
def _fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}"
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if not math.isfinite(v):
            return "--"
        return f"{v:.2e}" if (v != 0 and abs(v) < 0.001) else f"{v:.3g}"
    return str(v)


def _label_column(key: str, settings: dict, ctx) -> str | None:
    from starplast import strategies as S
    for name in ("target", "exclude", "a"):
        if settings.get(name):
            return str(settings[name])
    s = S.get(key)
    for p in s.params:
        if p.name in ("target", "exclude"):
            v = p.resolve(ctx)
            if v:
                return str(v)
    return None


def _unknown_ids(ctx, column: str | None) -> set | None:
    """Gene ids with no value in `column`, or None when there is no column to judge by."""
    if not column or column not in ctx.nodes:
        return None
    s = ctx.nodes[column]
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        missing = ctx.values(column).isna().to_numpy()
    else:
        missing = ctx.truth(column).isna().to_numpy()
    return set(ctx.gene_ids[missing])


def _gene_rows(frame, unknown, call, support, ascending=False, sort=None, label="") -> tuple:
    """Top rows of a per-gene table for genes without a known label."""
    f = frame
    if unknown is not None and "gene_id" in f:
        f = f[f["gene_id"].astype(str).isin(unknown)]
    n = len(f)
    by = sort or support
    if by in f:
        key = f[by]
        if callable(ascending):
            f = f.iloc[np.argsort(ascending(key.to_numpy()), kind="stable")]
        else:
            f = f.sort_values(by, ascending=ascending, kind="stable")
    rows = [{"gene_id": str(r["gene_id"]), "product": str(r.get("product", "")),
             "call": call(r) if callable(call) else _fmt(r.get(call, "")),
             "support": support_text(r, support, label)} for _, r in f.head(TOP).iterrows()]
    return rows, n


def support_text(r, support, label="") -> str:
    if callable(support):
        return support(r)
    return f"{label or support} {_fmt(r.get(support))}"


def _pair_rows(frame, score, keep=None, label="", call="linked") -> tuple:
    f = frame if keep is None else frame[keep(frame)]
    n = len(f)
    f = f.sort_values(score, ascending=False, kind="stable") if score in f else f
    rows = [{"gene_id": f"{r['gene_a']} + {r['gene_b']}",
             "product": f"{r.get('product_a', '')} / {r.get('product_b', '')}",
             "call": call, "support": f"{label or score} {_fmt(r.get(score))}"}
            for _, r in f.head(TOP).iterrows()]
    return rows, n


def _finding_rows(frame, call, support, sort=None, ascending=False) -> tuple:
    f = frame
    if sort and sort in f:
        f = f.sort_values(sort, ascending=ascending, kind="stable")
    rows = [{"gene_id": "", "product": "", "call": call(r), "support": support(r)}
            for _, r in f.head(TOP).iterrows()]
    return rows, len(frame)


def _gene_list(ctx, hidden: str, column: str | None) -> str:
    """The gene list a gene-list strategy's self-test was built from: the hidden category's genes."""
    m = re.search(r"'(.+?)' genes of (\S+)", hidden or "")
    if not m:
        return ""
    cat, col = m.group(1), m.group(2)
    col = col if col in ctx.nodes else column
    t = ctx.truth(col)
    return "\n".join(ctx.gene_ids[(t == cat).to_numpy()])


def new_information(key: str, ctx, settings: dict, hidden: str = "") -> dict:
    """Run the strategy once and keep its top new calls, ranked genes or findings."""
    from starplast import strategies as S
    s = S.get(key)
    settings = {k: v for k, v in settings.items() if any(p.name == k for p in s.params)}
    column = _label_column(key, settings, ctx)
    if any(p.name == "genes" for p in s.params) and not settings.get("genes"):
        settings["genes"] = _gene_list(ctx, hidden, column)
        m = re.search(r"'(.+?)' genes of (\S+)", hidden or "")
        if m and m.group(2) in ctx.nodes:
            column = m.group(2)            # the list's own label: new means not labeled there
    t0 = time.monotonic()
    result = _bounded(lambda: s.run(ctx, **settings))
    seconds = round(time.monotonic() - t0, 1)
    if isinstance(result, str):
        return {"skipped": result, "seconds": seconds}
    unknown = _unknown_ids(ctx, column)
    listed = set(str(g).strip() for g in str(settings.get("genes") or "").split())
    what_col = f"no known {column}" if column else "no label"
    if listed:
        # A gene-list strategy's new information is a gene NOT on the list: a new member.
        unknown = set(ctx.gene_ids) - listed
        what_col = "no place on the list"
    rows, n, what = _extract(key, result.tables, unknown, what_col, settings)
    return {"rows": _clean(rows), "n_new": int(n), "what": what, "seconds": seconds,
            "column": column or "", "summary": result.summary[:600]}


def _extract(key: str, tables: dict, unknown, what_col: str, settings: dict) -> tuple:
    """(rows, how many, a phrase) from one strategy's result tables."""
    t = {k: v for k, v in tables.items() if isinstance(v, pd.DataFrame)}

    def first(*names):
        for n in names:
            if n in t and len(t[n]):
                return t[n]
        return None

    calls = first("calls", "candidates")
    if key in ("positive_unlabeled", "seed_expansion") and calls is not None:
        rows, n = _gene_rows(calls, unknown, lambda r: "likely member", "score")
        return rows, n, f"ranked {n:,} genes with {what_col} as likely members of the list"
    if key == "set_enrichment":
        f = first("resembles the profile")
        if f is not None:
            rows, n = _gene_rows(f, unknown, lambda r: "resembles the list", "score")
            return rows, n, f"ranked {n:,} genes with {what_col} by how much they resemble the list"
    if key == "geneset_hunt":
        f = first("candidates")
        if f is not None:
            f = f[~f["in_list"].astype(bool)] if "in_list" in f else f
            rows, n = _gene_rows(f, unknown, lambda r: "in the list's cluster",
                                 "distance_to_centre", ascending=True, label="distance to centre")
            return rows, n, f"placed {n:,} genes with {what_col} in the cluster that holds the list"
    if key == "understudied_first" and calls is not None:
        rows, n = _gene_rows(calls, unknown, "prediction",
                             lambda r: f"support {_fmt(r.get('support'))}, "
                                       f"{_fmt(r.get('sources_agreeing'))} sources agree")
        return rows, n, f"called {n:,} little-studied genes with {what_col}"
    if key == "conformal_calls" and calls is not None:
        rows, n = _gene_rows(calls, unknown, "prediction",
                             lambda r: f"set {r.get('prediction_set')} (size {_fmt(r.get('set_size'))})",
                             sort="set_size", ascending=True)
        return rows, n, f"gave prediction sets to {n:,} genes with {what_col}"
    if calls is not None and "prediction" in calls:
        rows, n = _gene_rows(calls, unknown, "prediction", "support")
        return rows, n, f"called {n:,} genes with {what_col}"
    if key == "consensus_modules":
        mods, mem = first("modules"), first("members")
        if mods is not None and mem is not None:
            top = mods.set_index("module")[["top_label", "top_share", "stability"]]
            f = mem.join(top, on="module")
            f = f[f["top_label"].notna()]
            rows, n = _gene_rows(f, unknown,
                                 lambda r: f"module {r['module']} ({r['top_label']})",
                                 lambda r: f"{_fmt(r['top_share'])} of its labeled genes, "
                                           f"stability {_fmt(r['stability'])}",
                                 sort="top_share")
            return rows, n, f"placed {n:,} genes with {what_col} in stable modules"
    if key == "label_outliers":
        f = first("surprising genes")
        if f is not None:
            rows, n = _gene_rows(f, None, lambda r: f"labeled {r['label']}, looks {r['alternative']}",
                                 "surprise")
            return rows, n, f"flagged {n:,} labeled genes whose label looks wrong"
    if key in ("trait_regression", "conformal_values", "condition_shift", "masked_imputation",
               "ortholog_transfer"):
        spec = {"trait_regression": ("predicted", "predicted"),
                "condition_shift": ("predicted", "predicted_shift")
                if first("predicted") is not None else ("shifted genes", "shift"),
                "conformal_values": ("predicted intervals", "predicted"),
                "masked_imputation": ("filled", "imputed_percentile"),
                "ortholog_transfer": ("transferred", "transferred")}[key]
        f = first(spec[0])
        if f is not None:
            col = spec[1]
            if key == "conformal_values":
                f = f.assign(width=f["upper"] - f["lower"])
                sup = (lambda r: f"interval {_fmt(r['lower'])} to {_fmt(r['upper'])}")
                rows, n = _gene_rows(f, unknown, lambda r: f"predicted {_fmt(r[col])}", sup,
                                     sort="width", ascending=True)
                return rows, n, f"predicted {n:,} unmeasured genes, each with a 90% interval"
            if key == "masked_imputation":
                sup = (lambda r: f"column reliability {_fmt(r.get('column_reliability'))}")
                rows, n = _gene_rows(f, None, lambda r: f"percentile {_fmt(r[col])}", sup,
                                     sort=col, ascending=lambda v: -np.abs(v - 50))
                return rows, n, f"filled {n:,} missing measurements"
            med = float(np.nanmedian(pd.to_numeric(f[col], errors="coerce"))) if len(f) else 0.0
            ext = (lambda v: -np.abs(v - med))
            if col == "shift":                  # every gene measured: the shift itself is new
                rows, n = _gene_rows(f, None, lambda r: f"condition-specific shift {_fmt(r[col])}",
                                     lambda r: f"explained {_fmt(r.get('explained'))}", sort=col,
                                     ascending=lambda v: -np.abs(v))
                return rows, n, f"measured the condition-specific shift of {n:,} genes"
            rows, n = _gene_rows(f, unknown if key != "ortholog_transfer" else None,
                                 lambda r: f"{col.replace('_', ' ')} {_fmt(r[col])}",
                                 lambda r: f"furthest from the median, {_fmt(med)}", sort=col,
                                 ascending=ext)
            return rows, n, f"predicted {n:,} genes with {what_col.replace('known', 'measured')}"
    if key in ("link_prediction", "attention_correction", "unwritten_links", "paralog_divergence",
               "neighbour_space", "network_training"):
        spec = {"link_prediction": ("predicted links", "score", None, "predicted link"),
                "attention_correction": ("corrected ranking", "corrected", None,
                                         "linked beyond fame"),
                "unwritten_links": ("unwritten pairs", "layers", None, "measured, unwritten"),
                "paralog_divergence": ("paralog pairs", "divergence", None, "diverged paralogs"),
                "neighbour_space": ("the space", "probability",
                                    lambda f: f["kind"].astype(str) != "measured",
                                    "predicted pair"),
                "network_training": ("gaps", "probability", None, "predicted pair")}[key]
        f = first(spec[0])
        if f is not None:
            rows, n = _pair_rows(f, spec[1], spec[2], call=spec[3])
            phrase = {"link_prediction": "proposed {n:,} links not in the layer",
                      "attention_correction": "re-ranked {n:,} pairs once fame is taken out",
                      "unwritten_links": "found {n:,} measured pairs nobody has written about",
                      "paralog_divergence": "ranked {n:,} paralog pairs by how far they diverged",
                      "neighbour_space": "proposed {n:,} unmeasured pairs",
                      "network_training": "proposed {n:,} unmeasured pairs"}[key]
            if n == 0 and key == "neighbour_space":
                return rows, n, (f"ranked already-measured pairs at the top: all {len(f):,} "
                                 f"pairs it listed are recorded by some layer, so this run "
                                 f"proposes no unmeasured pair")
            return rows, n, phrase.format(n=n)
    if key == "recoverability_atlas":
        f = first("categories")
        if f is not None:
            rows, n = _finding_rows(f, lambda r: f"{r['category']} is recoverable",
                                    lambda r: f"AUROC {_fmt(r['auroc'])}, lift {_fmt(r['lift'])}",
                                    sort="auroc")
            return rows, n, f"ranked {n:,} categories by how well the data recovers them"
    if key == "blind_battery":
        f = first("held-out features")
        if f is not None:
            rows, n = _finding_rows(f, lambda r: f"map separates {r['feature']}",
                                    lambda r: f"effect {_fmt(r['effect'])}, q {_fmt(r['q'])}",
                                    sort="q", ascending=True)
            return rows, n, f"tested {n:,} held-out features the map was never shown"
    if key == "block_ablation":
        f = first("evidence")
        if f is not None:
            rows, n = _finding_rows(f, lambda r: f"{r['unit']} carries the label",
                                    lambda r: f"alone {_fmt(r['alone'])}, "
                                              f"loss when removed {_fmt(r['loss_when_removed'])}",
                                    sort="alone")
            return rows, n, f"ranked {n:,} kinds of evidence by what each carries alone"
    if key in ("split_clusters", "conjunctions"):
        frames = [v for v in t.values() if len(v) and "category" in v]
        if frames:
            f = pd.concat(frames, ignore_index=True)
            other = "other_category" if "other_category" in f else "other"
            rows, n = _finding_rows(
                f, lambda r: f"{r['category']} x {r.get(other, '')} ({r.get('layer', '')})",
                lambda r: f"purity {_fmt(r.get('purity'))}, {_fmt(r.get('n_cluster'))} genes",
                sort="purity")
            return rows, n, f"found {n:,} structures neither label shows alone"
    if calls is not None:
        rows, n = _gene_rows(calls, unknown, "prediction" if "prediction" in calls else "gene_id",
                             "support" if "support" in calls else "gene_id")
        return rows, n, f"listed {n:,} genes"
    return [], 0, "produced no new calls on this setting"


class _Budget(Exception):
    pass


def _bounded(fn):
    """fn() or, past `RUN_BUDGET` seconds or on an error, a string saying why not."""
    def alarm(_sig, _frame):
        raise _Budget()
    old = signal.signal(signal.SIGALRM, alarm)
    signal.alarm(RUN_BUDGET)
    try:
        return fn()
    except _Budget:
        return f"the run took longer than {RUN_BUDGET // 60} minutes, the budget for this page"
    except Exception as exc:                     # recorded, never raised: one strategy, one entry
        return f"the run stopped with {type(exc).__name__}: {str(exc)[:160]}"
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


# --------------------------------------------------------------------------- building
def noise_failure(key: str, task: str) -> dict:
    """The strategy's self-test on the noise table, as a failure example."""
    from starplast import strategies as S
    ctx = S.planted_context(null=True)
    s = S.get(key)
    settings = s.defaults(ctx)
    test = _bounded(lambda: s.test(ctx, **settings))
    if isinstance(test, str):
        return {"source": "noise table", "verdict": "NOT RUN", "note": test, "card": {},
                "explanation": f"The noise-table self-test could not run: {test}."}
    ex = example_from_test(test, task, settings)
    ex["explanation"] = failure_mode(ex, task)
    return ex


def build_one(key: str, organism: str, runs: pd.DataFrame, ctx, run_new: bool = True) -> dict:
    """The two worked examples for one strategy on one organism."""
    from starplast import strategies as S
    s = S.get(key)
    pick = choose(runs, key, organism)
    out = {"task": s.task, "runs": pick["runs"], "passes": pick["passes"], "fails": pick["fails"]}
    if pick["failure"] is not None:
        f = example_from_run(pick["failure"], s.task)
        f["at_default"] = pick["failure_at_default"]
        f["explanation"] = failure_mode(f, s.task)
    else:
        f = noise_failure(key, s.task)
        f["why_noise"] = (f"none of its {pick['runs']:,} calibration runs on real data failed "
                          f"({pick['passes']:,} passed)")
    out["failure"] = f
    if pick["success"] is not None:
        ok = example_from_run(pick["success"], s.task)
        ok["at_default"] = pick["success_at_default"]
        if run_new:
            ctx.seed = ok["seed"]
            new = new_information(key, ctx, dict(ok["settings"]), ok["hidden"])
        else:
            new = {"skipped": "not re-run in this build"}
        ok["new"] = new
        ok["explanation"] = success_text(ok, new)
    else:
        ok = {"source": "none", "verdict": "none",
              "explanation": (f"None of its {pick['runs']:,} calibration runs on real data passed"
                              + (f"; the usual reason: {pick['note']}" if pick["note"] else "")
                              + ". There is no success to show.")}
    out["success"] = ok
    return _clean(out)


def build(organisms=ORGANISMS, only=None, run_new: bool = True, log=print) -> dict:
    """Every example for every strategy and organism, as the JSON the application ships."""
    import warnings
    warnings.filterwarnings("ignore")
    from starplast import strategies as S
    runs = load_runs()
    data = {"meta": {"built": _dt.date.today().isoformat(),
                     "runs": os.path.relpath(RUNS, ROOT), "run_budget_seconds": RUN_BUDGET}}
    for org in organisms:
        ctx = S.Context.shipped(org)
        _LIGHT[org] = ctx
        data[org] = {}
        for s in S.catalog():
            if only and s.key not in only:
                continue
            t0 = time.monotonic()
            with contextlib.redirect_stderr(open(os.devnull, "w")):
                data[org][s.key] = build_one(s.key, org, runs, ctx, run_new=run_new)
            e = data[org][s.key]
            log(f"{org} {s.number:02d} {s.key}: failure {e['failure'].get('source')} "
                f"{e['failure'].get('verdict')}, success {e['success'].get('verdict')} "
                f"({time.monotonic() - t0:.0f}s)")
    return data


def write(data: dict, path: str = OUT, merge: bool = True) -> str:
    """Write the JSON, keeping entries for strategies and organisms this build did not touch."""
    old = {}
    if merge and os.path.exists(path):
        with open(path) as fh:
            old = json.load(fh)
    for org, entries in data.items():
        if org == "meta":
            old["meta"] = entries
        else:
            old.setdefault(org, {}).update(entries)
    with open(path, "w") as fh:
        json.dump(old, fh, indent=1, sort_keys=False)
        fh.write("\n")
    return path


def summary(data: dict) -> pd.DataFrame:
    """One row per organism and strategy: where each example came from, and its skill."""
    rows = []
    for org in ORGANISMS:
        for key, e in (data.get(org) or {}).items():
            f, ok = e["failure"], e["success"]
            rows.append({"organism": org, "strategy": key, "task": e["task"],
                         "failure_from": f.get("source"), "failure_verdict": f.get("verdict"),
                         "failure_target": f.get("target"), "failure_skill": f.get("skill"),
                         "success_target": ok.get("target"), "success_skill": ok.get("skill"),
                         "new_rows": len((ok.get("new") or {}).get("rows") or []),
                         "new_skipped": (ok.get("new") or {}).get("skipped", "")})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- the notebook
def notebook(organisms=ORGANISMS, only=None, run_new: bool = True, path: str = NOTEBOOK) -> str:
    from notebook_runner import ExecutedNotebook
    nb = ExecutedNotebook("Worked examples for the strategy cards", {"__name__": "__notebook__"})
    nb.md("Each strategy card shows one real failure and one real success. This notebook chooses "
          "them from the calibration sweep on disk and, for each success, runs the strategy once "
          "to list what it says about genes without a known label. The code is "
          "`scripts/build_strategy_examples.py`; this notebook calls it, so what is shown is what "
          "ships in `starplast/data/strategy_examples.json`.",
          "Nothing here changes a calibration number: the failure and the success are records "
          "already in `runs.jsonl`, and the explanations are written from their own numbers.")
    nb.code("import warnings; warnings.filterwarnings('ignore')",
            f"import sys; sys.path.insert(0, {ROOT!r}); sys.path.insert(0, {os.path.join(ROOT, 'scripts')!r})",
            "import build_strategy_examples as B",
            "runs = B.load_runs()",
            "print(len(runs), 'runs;', runs['verdict'].value_counts().to_dict())")
    nb.md("## 1. The choice rule, strategy by strategy",
          "Success: the PASS at the default setting with the highest skill. Failure: the FAIL at "
          "the default setting with the lowest skill, else any FAIL, else the self-test on the "
          "noise table.")
    nb.code("from starplast import strategies as S",
            "rows = []",
            f"for org in {tuple(organisms)!r}:",
            "    for s in S.catalog():",
            f"        if {only!r} and s.key not in {only!r}: continue",
            "        c = B.choose(runs, s.key, org)",
            "        rows.append({'organism': org, 'strategy': s.key, 'runs': c['runs'], "
            "'passes': c['passes'], 'fails': c['fails'], "
            "'success_at_default': c['success_at_default'], "
            "'failure_at_default': c['failure_at_default']})",
            "import pandas as pd",
            "choice = pd.DataFrame(rows)",
            "choice")
    nb.md("## 2. Build every example",
          f"Each success is run once on the shipped table; a run longer than {RUN_BUDGET} seconds "
          "is stopped and recorded as skipped.")
    nb.code(f"data = B.build({tuple(organisms)!r}, only={only!r}, run_new={run_new!r})",
            "table = B.summary(data)", "table")
    nb.md("## 3. Where the failures came from",
          "A real-data failure where the sweep has one; the noise table where it has none.")
    nb.code("table.groupby(['organism', 'failure_from']).size()")
    nb.md("## 4. Two examples in full")
    nb.code("import json",
            "org = next(o for o in data if o != 'meta')",
            "for k in list(data[org])[:2]:",
            "    e = data[org][k]",
            "    print(org, k)",
            "    print('FAILS WHEN:', e['failure']['explanation'])",
            "    print('WORKS WHEN:', e['success']['explanation'])",
            "    print(json.dumps((e['success'].get('new') or {}).get('rows'), indent=1))")
    nb.md("## 5. Every explanation",
          "The sentence each card shows under its two mini scorecards, for every strategy.")
    nb.code("for o in data:",
            "    if o == 'meta': continue",
            "    for k, e in data[o].items():",
            "        print(f'{o} {k}\\n  FAILS WHEN: {e[\"failure\"][\"explanation\"]}'",
            "              f'\\n  WORKS WHEN: {e[\"success\"][\"explanation\"]}')")
    nb.md("## 6. Write the shipped file")
    nb.code("B.write(data)")
    return nb.write(path)


# --------------------------------------------------------------------------- the docs page
DOC = os.path.join(ROOT, "docs", "strategy_cards.md")
def _italic(code: str) -> str:
    from starplast import organisms
    genus, _, rest = organisms.get(code).species.partition(" ")
    return f"*{genus[0]}. {rest}*"


SPECIES = {c: _italic(c) for c in ORGANISMS}


def _bar(position: float, width: int = 12) -> str:
    if not isinstance(position, (int, float)) or not math.isfinite(position):
        return "`" + "·" * width + "`"
    n = int(round(position * width))
    return "`" + "█" * n + "░" * (width - n) + "`"


def _cell(bar: dict) -> str:
    from starplast import scorecard as SC
    v = SC.fmt_value(bar["scale"], bar["value"])
    if v == "--":
        return "--"
    text = f"{_bar(bar['position'])} **{v}**"
    if math.isfinite(_f(bar["low"])) and math.isfinite(_f(bar["high"])):
        text += (f" [{SC.fmt_value(bar['scale'], bar['low'])}, "
                 f"{SC.fmt_value(bar['scale'], bar['high'])}]")
    if math.isfinite(_f(bar["chance"])):
        text += f" · chance {SC.fmt_value(bar['scale'], bar['chance'])}"
    return text


def card_markdown(key: str, data: dict | None = None) -> str:
    """One strategy's card as Markdown: bars for both organisms, the explainer, the examples."""
    from starplast import calibration as CAL
    from starplast import scorecard as SC
    from starplast import strategies as S
    from starplast import strategy_explainers as EX
    s = S.get(key)
    grades = " · ".join(f"{SPECIES[o]} **{(CAL.entry(key, o) or {}).get('grade') or 'not calibrated'}**"
                        for o in ORGANISMS)
    lines = [f"## {s.number:02d} · {s.title}", "",
             f"`{s.method}` · {s.task} · {grades}", "", f"*{s.question}*", "",
             f"| | {' | '.join(SPECIES[o] for o in ORGANISMS)} |", "|---|" + "---|" * len(ORGANISMS)]
    bars = {}
    for o in ORGANISMS:
        e = CAL.entry(key, o) or {}
        d = e.get("default") or {}
        card = d.get("scorecard") or {}
        bars[o] = SC.headline_bars(s.task, card, d)
    for i, h in enumerate(SC.headline(s.task)):
        lines.append(f"| **{h.label}** <br><sub>{h.technical}</sub> | "
                     + " | ".join(_cell(bars[o][i]) for o in ORGANISMS) + " |")
    lines += ["", "**About this test**", "", EX.markdown(key), ""]
    for o in ORGANISMS:
        ex = ((data or {}).get(o) or {}).get(key) or EX.examples(key, o)
        if not ex:
            continue
        f, ok = ex.get("failure") or {}, ex.get("success") or {}
        tag = " (noise table)" if f.get("source") == "noise table" else ""
        lines.append(f"**{SPECIES[o]}.** *Fails when* {f.get('target') or 'random data'}{tag}: "
                     f"{f.get('explanation', '')}  ")
        lines.append(f"*Works when* {ok.get('target') or '--'}: {ok.get('explanation', '')}")
        rows = (ok.get("new") or {}).get("rows") or []
        if rows:
            lines.append("")
            for r in rows:
                who = f"`{r['gene_id']}` {r.get('product', '')} -- " if r.get("gene_id") else ""
                lines.append(f"  - {who}{r.get('call', '')} ({r.get('support', '')})")
        lines.append("")
    return "\n".join(lines)


def write_docs(path: str = DOC, data: dict | None = None) -> str:
    """docs/strategy_cards.md: every strategy's card, as the Strategies tab shows it."""
    from starplast import scorecard as SC
    from starplast import strategies as S
    head = ["# Strategy cards", "",
            "Generated by `scripts/build_strategy_examples.py --docs-only`; do not edit by hand.", "",
            "Each strategy shows the same four bars, in the same places: **Better than chance** "
            "(skill: 0 is the same procedure on shuffled data, 1 is perfect), **Reach** (how much "
            "of the question it can speak to), and the two metrics of its task a biologist asks "
            "about first, in plain words. Values are the mean over the calibration's held-out tests "
            "at default settings, with the 95% interval and the chance level.", "",
            "| Task | Reach means | Bar 3 | Bar 4 |", "|---|---|---|---|"]
    for t in SC.TASKS:
        r, (a, b) = SC.REACH[t], SC.HEADLINE[t]
        head.append(f"| {t} | {r.technical}{' -- ' + r.note if r.note else ''} | "
                    f"{a.label} ({a.technical}) | {b.label} ({b.technical}) |")
    head += ["", "Below the bars, **About this test** says what the strategy does, how it is "
             "evaluated, what failure looks like and what success looks like. The worked examples "
             "are real calibration runs: the failure is the target it failed on most often at "
             "default settings (or, where it never failed on real data, its self-test on a table "
             "of random labels and edges), the success its best pass, re-run once to list what it "
             "says about genes without a known label.", ""]
    body = [card_markdown(s.key, data) for s in S.catalog()]
    with open(path, "w") as fh:
        fh.write("\n".join(head) + "\n" + "\n".join(body))
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--organism", choices=ORGANISMS, action="append")
    ap.add_argument("--only", default="", help="comma-separated strategy keys")
    ap.add_argument("--no-run", action="store_true", help="skip the new-information runs")
    ap.add_argument("--no-notebook", action="store_true",
                    help="write the JSON without recording the notebook (for development)")
    ap.add_argument("--docs-only", action="store_true",
                    help="regenerate docs/strategy_cards.md from the shipped JSON")
    args = ap.parse_args(argv)
    organisms = tuple(args.organism or ORGANISMS)
    only = sorted({k for k in args.only.split(",") if k}) or None
    if args.docs_only:
        print(write_docs())
        return 0
    if args.no_notebook:
        print(write(build(organisms, only, run_new=not args.no_run)))
    else:
        print(notebook(organisms, only, run_new=not args.no_run))
    print(write_docs())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
