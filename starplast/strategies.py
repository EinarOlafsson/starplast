"""Strategies: named ways of turning the combined data into a claim, each with a test of itself.

A *strategy* is a way of reasoning from the measurements to something nobody has measured -- hold a
category out and look for a map that finds it again, hand over a gene list and look for the one
cluster that holds it, carry a label along a crosslink, predict a screen from the other screens. The
analysis tabs are the instruments; a strategy is a way of playing them, written down so that it can
be repeated, explained and -- the part that matters -- checked.

Every strategy carries four things, and the Strategies tab shows all four:

* a **tooltip** that says in one breath what it infers and why that is worth doing;
* an **explanation** of why the inference works, what it assumes, and how it fails;
* a **walkthrough**: what to choose, what to press, and what to read in the result;
* a **self-test**. Known information is hidden, the strategy is asked to recover it from everything
  else, and the answer is scored against the SAME procedure run on shuffled labels (or random gene
  sets, or permuted node identities). A strategy passes only when it beats the 95th percentile of
  its own null by a stated margin. The null is re-run rather than assumed, because every procedure
  here has a chance level of its own -- a clustering calls some genes right by luck, a random gene
  list sits in some cluster -- and a raw accuracy without its chance level is the number this
  project exists not to publish.

The leakage guard is applied before anything else. A strategy scored against `compartment` never
sees `compartment`, anything measurably restating it, anything the same experiment produced, or an
edge layer built from it: :meth:`Context.banned` and :meth:`Context.banned_layers` are the project's
own closure (`search.excluded_detail`, `search.excluded_edges`), not a second copy of it.

The strategies themselves live in `strategy_catalog`; this module is the machinery they share:
the :class:`Context` they run against, the :class:`StrategyResult` and :class:`TestResult` they
return, the inference primitives several of them use, and the five holdout test patterns.

Headless and Qt-free, so everything here can be driven from a script or a notebook::

    from starplast import strategies
    ctx = strategies.Context.shipped("Tg")
    s = strategies.get("holdout_search")
    print(s.test(ctx).summary())
    result = s.run(ctx, target="compartment")
"""
from __future__ import annotations

import math
import os
import re
import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd

from . import scorecard as SC
from .stopping import Stopped

#: Label values that mean "not measured" rather than naming a class. Shared with `search`, so a
#: strategy and a search agree on what counts as a label.
ABSENT = frozenset({"unassigned", "unknown", "", "nan", "none", "unlabelled", "unlabeled", "<na>"})
#: The smallest class a strategy is scored on. Below this a class cannot be split into a visible
#: part and a hidden part that both mean anything.
MIN_CLASS = 15
#: A self-test with fewer hidden items than this is reported as inconclusive rather than scored.
MIN_HIDDEN = 10
#: Clustering noise, as every clustering in this project writes it.
NOISE = -1
#: Columns that are never features. Identifiers and free text are not measurements, and literature
#: attention is excluded everywhere in this project: a map organised by how often a gene is written
#: about is a map of the field, not of the parasite.
NEVER_FEATURES = re.compile(r"^(gene_id|sequence|product|alphafold_accession|gene_name|symbol)$"
                            r"|^n_papers_|^attention|^lit_|^n_abstracts|^mention")
#: Edge layers that are computed FROM other layers, and what they are computed from. A link
#: predictor scored against `xlms` must not read `unwritten_interaction`, which is `xlms` minus the
#: pairs the literature mentions -- the declared-family guard cannot see a derivation like that.
DERIVED_LAYERS = {"unwritten_interaction": ("xlms", "ip_ms", "comention", "comention_ft"),
                  "structural_hole": ("coexpression", "cofitness", "comention", "comention_ft")}
#: Literature layers. Kept apart from the measurement layers wherever a strategy asks what the
#: MEASUREMENTS say, because co-mention follows attention.
LITERATURE_LAYERS = ("comention", "comention_ft")
#: Parameter kinds the panel knows how to draw.
PARAM_KINDS = ("category", "number", "column", "layer", "genes", "int", "float", "grid", "choice")


def _quiet(_message):
    """A log that says nothing, for the library calls a strategy makes on its own behalf."""


# --------------------------------------------------------------------------- small numerics
def auroc(score, positive) -> float:
    """Area under the ROC curve: the chance a random positive outranks a random negative.

    Missing scores rank below everything, which is what a missing score means here: the method had
    nothing to say about that gene. Ties are averaged. NaN when either side is empty.
    """
    from scipy.stats import rankdata
    s = np.asarray(score, dtype=float).copy()
    y = np.asarray(positive, dtype=bool)
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    finite = np.isfinite(s)
    if not finite.all():
        floor = s[finite].min() - 1.0 if finite.any() else 0.0
        s[~finite] = floor
    r = rankdata(s)
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def spearman(a, b) -> float:
    """Rank correlation over the pairs where both are finite. NaN below five pairs or on a constant."""
    from scipy.stats import spearmanr
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5 or np.ptp(a[ok]) == 0 or np.ptp(b[ok]) == 0:
        return float("nan")
    return float(spearmanr(a[ok], b[ok]).statistic)


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values. Every scan here tests many things at once."""
    p = np.asarray(p, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def hypergeom_sf(hits: int, draws: int, successes: int, total: int) -> float:
    """P(at least `hits` successes in `draws` draws), exactly. 1.0 when nothing was drawn."""
    from scipy.stats import hypergeom
    if draws <= 0 or hits <= 0:
        return 1.0
    return float(hypergeom.sf(hits - 1, total, successes, draws))


def parse_grid(text, kind=float) -> tuple:
    """A comma-separated list of values, as the panel's grid boxes hold them. Blanks are skipped."""
    if isinstance(text, (list, tuple)):
        return tuple(kind(v) for v in text)
    return tuple(kind(v.strip()) for v in str(text).split(",") if v.strip())


def expand(values, positions, n: int, fill=NOISE) -> np.ndarray:
    """Values defined on some genes, written into an array over all `n` genes."""
    out = np.full(n, fill, dtype=np.asarray(values).dtype if len(values) else int)
    out[np.asarray(positions, dtype=int)] = values
    return out


# --------------------------------------------------------------------------- results
@dataclass
class StrategyResult:
    """What one run of a strategy produced, including when that is nothing.

    `tables` are the outputs a person reads, first one first. `labels` (one per gene, -1 where the
    gene has none) and `coords` over `positions` are what the map can show; both are optional because
    most strategies produce a list, not a map. `genes` are the ones worth highlighting.
    """
    strategy: str
    summary: str
    tables: dict = field(default_factory=dict)
    labels: np.ndarray | None = None
    coords: np.ndarray | None = None
    positions: np.ndarray | None = None
    genes: list = field(default_factory=list)
    numbers: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)
    seconds: float = 0.0

    @property
    def ok(self) -> bool:
        """Whether the run produced anything to read beyond its summary."""
        return any(len(t) for t in self.tables.values())

    def save(self, folder: str) -> list:
        """Write every table as CSV and the summary and settings as JSON. Returns the paths."""
        import json
        os.makedirs(folder, exist_ok=True)
        paths = []
        for name, table in self.tables.items():
            path = os.path.join(folder, f"{self.strategy}_{re.sub(r'[^A-Za-z0-9]+', '_', name)}.csv")
            pd.DataFrame(table).to_csv(path, index=False)
            paths.append(path)
        path = os.path.join(folder, f"{self.strategy}_summary.json")
        with open(path, "w") as fh:
            json.dump({"strategy": self.strategy, "summary": self.summary,
                       "settings": self.settings, "numbers": self.numbers,
                       "seconds": self.seconds}, fh, indent=1, default=str)
        return paths + [path]

    def plot(self, ctx: "Context | None" = None, color=None, ax=None, size: float = 3.0):
        """The map this result was computed on, colored by its clusters or by a column of `ctx`.

        Only for results that carry a map (`coords`); a list-producing strategy has none. Gray is
        a gene with no cluster or no value -- the application's rule, kept here.
        """
        import matplotlib.pyplot as plt
        if self.coords is None or self.positions is None:
            raise ValueError(f"{self.strategy} produced no map to plot")
        xyz = np.asarray(self.coords)
        pos = np.asarray(self.positions, dtype=int)
        if color is not None and ctx is not None and color in ctx.nodes:
            values = ctx.truth(color).iloc[pos].to_numpy(dtype=object)
        elif self.labels is not None:
            values = np.asarray(self.labels)[pos].astype(object)
            values = np.where(values == NOISE, None, values)
        else:
            values = np.array([None] * len(pos), dtype=object)
        if ax is None:
            fig = plt.figure(figsize=(7, 6))
            ax = fig.add_subplot(projection="3d")
        known = np.array([v is not None and v == v for v in values])
        ax.scatter(*xyz[~known].T, s=size, c="#8a8a8a", alpha=0.35, linewidths=0)
        cats = sorted({str(v) for v in values[known]})
        cmap = plt.get_cmap("tab20", max(len(cats), 1))
        for i, c in enumerate(cats):
            m = known & (values.astype(str) == c)
            ax.scatter(*xyz[m].T, s=size * 1.6, color=cmap(i % 20), label=c if len(cats) <= 20
                       else None, linewidths=0)
        if 0 < len(cats) <= 20:
            ax.legend(fontsize=7, markerscale=3, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        ax.set_axis_off()
        return ax

    def _repr_html_(self) -> str:
        from html import escape
        first = next(iter(self.tables.items()), None)
        table = (f"<p><b>{escape(first[0])}</b> ({len(first[1]):,} rows; first 10)</p>"
                 + pd.DataFrame(first[1]).head(10).to_html(index=False)) if first else ""
        return (f"<p><b>{escape(self.strategy)}</b> -- {escape(self.summary)}</p>"
                f"<p>tables: {escape(', '.join(f'{k} ({len(v):,})' for k, v in self.tables.items()))}"
                f"</p>{table}")


@dataclass
class TestResult:
    """A strategy's self-test: what was hidden, how much came back, and what luck alone gives.

    `null` holds the metric under the null procedure -- the same strategy on shuffled labels, random
    gene sets or permuted identities -- so the chance level is measured rather than assumed. Where a
    re-run is too costly the analytic expectation is used instead and `null_kind` says so.
    """
    strategy: str
    metric: str
    observed: float
    null: list
    min_effect: float
    n_hidden: int
    hidden: str
    null_kind: str
    details: pd.DataFrame = field(default_factory=pd.DataFrame)
    seconds: float = 0.0
    note: str = ""
    quantile: float = 95.0
    analytic: tuple | None = None
    numbers: dict = field(default_factory=dict)
    #: How many hidden items make a verdict conclusive. Hidden genes for most patterns; findings to
    #: replicate for pattern 5, where three findings are already a real test.
    min_hidden: int = MIN_HIDDEN
    #: Which of the scorecard's tasks this test measured ("label calls", "ranking", ...), and that
    #: task's standard metrics on the same hidden genes -- see :mod:`starplast.scorecard`.
    task: str = ""
    scorecard: dict = field(default_factory=dict)

    @property
    def null_mean(self) -> float:
        """The chance level: the mean of the null runs, or the analytic expectation."""
        if self.analytic is not None:
            return float(self.analytic[0])
        vals = [v for v in self.null if np.isfinite(v)]
        return float(np.mean(vals)) if vals else float("nan")

    @property
    def null_sd(self) -> float:
        """How much the null runs vary, which is what makes a gap over them mean something."""
        if self.analytic is not None:
            return float(self.analytic[1])
        vals = [v for v in self.null if np.isfinite(v)]
        return float(np.std(vals)) if len(vals) > 1 else 0.0

    @property
    def null_high(self) -> float:
        """The bar to clear: the `quantile`-th percentile of the null, analytic or re-run."""
        if self.analytic is not None:
            z = {95.0: 1.645, 99.0: 2.326, 90.0: 1.282, 80.0: 0.842}.get(float(self.quantile), 1.645)
            return float(self.analytic[0] + z * self.analytic[1])
        vals = [v for v in self.null if np.isfinite(v)]
        return float(np.percentile(vals, self.quantile)) if vals else float("nan")

    @property
    def effect(self) -> float:
        """Observed minus chance: the part of the answer that luck does not explain."""
        return float(self.observed - self.null_mean)

    @property
    def skill(self) -> float:
        """(observed - chance) / (1 - chance): 0 is chance, 1 is perfect -- one scale for all."""
        chance = self.null_mean
        if not (np.isfinite(self.observed) and np.isfinite(chance)) or chance >= 1:
            return float("nan")
        return float((self.observed - chance) / (1 - chance))

    def card(self) -> pd.DataFrame:
        """The scorecard as one table: the verdict block every test shares, then the task's metrics
        in their standard order, each with its one-line reading."""
        from . import scorecard as SC
        head = {"verdict": self.verdict, "metric": self.metric, "observed": self.observed,
                "chance": self.null_mean, "bar": self.null_high, "p_value": self.p_value,
                "skill": self.skill, "n_hidden": self.n_hidden}
        rows = [{"section": "verdict", "metric": label, "key": k, "value": head[k],
                 "reading": text} for k, label, text in SC.VERDICT]
        if self.task in SC.TASKS:
            for _, r in SC.frame(self.scorecard, self.task).iterrows():
                rows.append({"section": self.task, "metric": r["metric"], "key": r["key"],
                             "value": r["value"], "reading": r["reading"]})
        return pd.DataFrame(rows)

    @property
    def p_value(self) -> float:
        """Share of null runs at least as good as the observed, with the usual +1 correction."""
        if self.analytic is not None:
            from scipy.stats import norm
            sd = self.analytic[1]
            return float(norm.sf((self.observed - self.analytic[0]) / sd)) if sd > 0 else float(
                self.observed <= self.analytic[0])
        vals = [v for v in self.null if np.isfinite(v)]
        return float((1 + sum(v >= self.observed for v in vals)) / (1 + len(vals)))

    @property
    def conclusive(self) -> bool:
        """Enough was hidden, and a number came back, for the verdict to mean anything."""
        return self.n_hidden >= self.min_hidden and np.isfinite(self.observed) and np.isfinite(
            self.null_high)

    @property
    def passed(self) -> bool:
        """Beats its own null's high percentile AND by at least the stated margin."""
        return bool(self.conclusive and self.observed > self.null_high
                    and self.effect >= self.min_effect)

    @property
    def verdict(self) -> str:
        """PASS, FAIL or INCONCLUSIVE -- the last when too little could be hidden to score."""
        if not self.conclusive:
            return "INCONCLUSIVE"
        return "PASS" if self.passed else "FAIL"

    def summary(self) -> str:
        """One line, leading with the verdict and the number it rests on."""
        if not self.conclusive:
            return (f"INCONCLUSIVE -- {self.note or 'too little could be hidden to score'} "
                    f"({self.n_hidden} hidden)")
        return (f"{self.verdict} -- {self.metric} {self.observed:.3f} against {self.null_mean:.3f} "
                f"± {self.null_sd:.3f} under {self.null_kind} (bar {self.null_high:.3f}, "
                f"p {self.p_value:.3g}); effect {self.effect:+.3f}, needed {self.min_effect:+.3f}; "
                f"{self.n_hidden:,} hidden")

    def save(self, folder: str) -> list:
        """Write the verdict as JSON and the per-class details as CSV. Returns the paths."""
        import json
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, f"{self.strategy}_test.json")
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=1)
        paths = [path]
        if len(self.details):
            paths.append(os.path.join(folder, f"{self.strategy}_test_details.csv"))
            self.details.to_csv(paths[-1], index=False)
        return paths

    def _repr_html_(self) -> str:
        from html import escape
        colour = {"PASS": "#2e8b57", "FAIL": "#c0392b"}.get(self.verdict, "#888")
        rest = escape(self.summary().split(" -- ", 1)[-1])
        rows = ""
        if self.task:
            from . import scorecard as SC
            for _, r in SC.frame(self.scorecard, self.task).iterrows():
                v = r["value"]
                shown = "--" if not np.isfinite(v) else (f"{v:,.0f}" if abs(v) >= 100 or float(
                    v).is_integer() and abs(v) >= 2 else f"{v:.3f}")
                rows += (f"<tr><td title='{escape(r['reading'])}'>{escape(r['metric'])}</td>"
                         f"<td style='text-align:right'>{shown}</td></tr>")
            rows = (f"<p><b>Scorecard ({escape(self.task)})</b></p><table>{rows}</table>"
                    if rows else "")
        return (f"<p><b style='color:{colour}'>{self.verdict}</b> -- {rest}</p>"
                f"<p><i>Hidden: {escape(self.hidden)}. Null: {escape(self.null_kind)}.</i></p>"
                + rows)

    def to_dict(self) -> dict:
        """JSON-safe, for the shipped record of what each strategy measured on the real data."""
        clean = lambda v: None if v is None or (isinstance(v, float) and not math.isfinite(v)) else v
        return {"strategy": self.strategy, "verdict": self.verdict, "metric": self.metric,
                "observed": clean(float(self.observed)), "null_mean": clean(self.null_mean),
                "null_sd": clean(self.null_sd), "null_high": clean(self.null_high),
                "p_value": clean(self.p_value), "effect": clean(self.effect),
                "min_effect": self.min_effect, "n_hidden": int(self.n_hidden),
                "hidden": self.hidden, "null_kind": self.null_kind, "note": self.note,
                "seconds": round(float(self.seconds), 1), "skill": clean(self.skill),
                "task": self.task,
                "scorecard": {k: clean(float(v)) for k, v in self.scorecard.items()},
                "numbers": {k: clean(v) for k, v in self.numbers.items()}}


def judge(strategy: str, metric: str, observed: float, null, *, min_effect: float,
          n_hidden: int, hidden: str, null_kind: str, t0: float, details=None, note: str = "",
          quantile: float = 95.0, analytic=None, numbers=None,
          min_hidden: int = MIN_HIDDEN, task: str = "", scorecard=None) -> TestResult:
    """Build a :class:`TestResult` -- the one place a verdict's arithmetic lives.

    `task` names the scorecard task and `scorecard` carries its metrics; an inconclusive test may
    leave both empty.
    """
    return TestResult(strategy=strategy, metric=metric, observed=float(observed),
                      null=[float(v) for v in (null or [])], min_effect=float(min_effect),
                      n_hidden=int(n_hidden), hidden=hidden, null_kind=null_kind,
                      details=details if details is not None else pd.DataFrame(),
                      seconds=time.monotonic() - t0, note=note, quantile=quantile,
                      analytic=analytic, numbers=dict(numbers or {}), min_hidden=int(min_hidden),
                      task=task, scorecard=dict(scorecard or {}))


# --------------------------------------------------------------------------- parameters and strategies
@dataclass(frozen=True)
class Param:
    """One setting a strategy takes, with the reason it exists.

    `default` may be a callable taking the :class:`Context`, because the right default depends on the
    organism -- `compartment` is the natural held-out label for Toxoplasma and does not exist for
    Plasmodium.
    """
    name: str
    kind: str
    label: str
    tip: str
    default: object = None
    lo: float = 0
    hi: float = 0
    step: float = 1
    choices: tuple = ()
    #: May be left empty -- a gene list with no column to exclude, say. The panel offers "(none)".
    optional: bool = False
    #: Whose columns a column parameter names: this organism's, or the `other` one's.
    space: str = "this"

    def resolve(self, ctx) -> object:
        """The default for this context."""
        return self.default(ctx) if callable(self.default) else self.default

    def options(self, ctx) -> tuple:
        """The fixed choices of a `choice` parameter; a callable is asked, since some depend on data."""
        return tuple(self.choices(ctx)) if callable(self.choices) else tuple(self.choices)


@dataclass
class Strategy:
    """A named way of inferring something, with its explanation, walkthrough and self-test."""
    key: str
    number: int
    title: str
    family: str
    question: str
    tooltip: str
    explanation: str
    walkthrough: tuple
    test_description: str
    params: tuple
    runner: Callable
    tester: Callable
    cost: str = "seconds"
    needs: tuple = ()
    method: str = ""
    #: The scorecard task its self-test reports (see :mod:`starplast.scorecard`), and the
    #: techniques its method is built from (see :mod:`starplast.techniques`).
    task: str = ""
    techniques: tuple = ()

    @property
    def name(self) -> str:
        """The title with the method it runs in brackets, as every list and heading shows it."""
        return f"{self.title} ({self.method})"

    def defaults(self, ctx) -> dict:
        """Every parameter at its default for this context."""
        return {p.name: p.resolve(ctx) for p in self.params}

    def settings(self, ctx, **overrides) -> dict:
        """Defaults with `overrides` applied. An unknown name is an error, not silently ignored."""
        unknown = set(overrides) - {p.name for p in self.params}
        if unknown:
            raise ValueError(f"{self.key} has no parameter {', '.join(sorted(unknown))}")
        out = self.defaults(ctx)
        out.update({k: v for k, v in overrides.items() if v is not None})
        return out

    def run(self, ctx, **params) -> StrategyResult:
        """Run the strategy on `ctx` and return what it found."""
        p = self.settings(ctx, **params)
        t0 = time.monotonic()
        result = self.runner(ctx, p)
        result.settings, result.seconds = p, time.monotonic() - t0
        return result

    def test(self, ctx, **params) -> TestResult:
        """Hide known information, ask the strategy for it back, and score that against its null."""
        p = self.settings(ctx, **params)
        return self.tester(ctx, p)

    def parameters(self, ctx: "Context | None" = None) -> pd.DataFrame:
        """The settings as a table: name, kind, default (for `ctx`), and why each exists."""
        return pd.DataFrame([{"name": p.name, "label": p.label, "kind": p.kind,
                              "default": p.resolve(ctx) if ctx is not None else (
                                  None if callable(p.default) else p.default),
                              "why": p.tip} for p in self.params])

    @property
    def metrics(self) -> tuple:
        """The scorecard metrics its self-test reports, in the task's standard order."""
        return SC.TASKS[self.task].metrics

    def techniques_table(self) -> pd.DataFrame:
        """The techniques its method is built from: name, kind, what each does and why."""
        from . import techniques
        return techniques.glossary(self.techniques)

    def scorecard_table(self) -> pd.DataFrame:
        """Its scorecard's metrics, each with definition, range, chance level and how to read it."""
        return SC.glossary(self.task)

    def _repr_html_(self) -> str:
        from html import escape
        from .techniques import TECHNIQUES
        steps = "".join(f"<li>{escape(x)}</li>" for x in self.walkthrough)
        tech = ", ".join(escape(TECHNIQUES[t].name) for t in self.techniques)
        mets = ", ".join(escape(SC.METRICS[m].label) for m in self.metrics)
        return (f"<p><b>{self.number:02d} · {escape(self.name)}</b> ({escape(self.family)})</p>"
                f"<p><i>{escape(self.question)}</i></p><p>{escape(self.tooltip)}</p>"
                f"<p><b>Techniques:</b> {tech}.<br><b>Scorecard ({escape(self.task)}):</b> "
                f"{mets}.</p>"
                f"<ol>{steps}</ol><p><b>How it is tested:</b> {escape(self.test_description)}</p>")

    def help_text(self) -> str:
        """Everything the panel shows about this strategy, as plain text."""
        from .techniques import TECHNIQUES
        steps = "\n".join(f"  {i}. {s}" for i, s in enumerate(self.walkthrough, 1))
        tech = "\n".join(f"  {TECHNIQUES[t].name}: {TECHNIQUES[t].what}" for t in self.techniques)
        mets = "\n".join(f"  {SC.METRICS[m].label}: {SC.METRICS[m].reading}" for m in self.metrics)
        return (f"{self.number:02d} · {self.name}\n\n{self.question}\n\n"
                f"Method: {self.method}\n{tech}\n\n{self.explanation}\n\n"
                f"Scorecard ({self.task})\n{mets}\n\n"
                f"Walkthrough\n{steps}\n\nHow it is tested\n{self.test_description}")


REGISTRY: dict = {}


def register(strategy: Strategy) -> Strategy:
    """Add a strategy to the catalogue. A second strategy under one key is refused."""
    if strategy.key in REGISTRY:
        raise ValueError(f"strategy {strategy.key!r} registered twice")
    if not strategy.method.strip():
        raise ValueError(f"{strategy.key}: name the method it runs, e.g. method=\"UMAP + HDBSCAN\"")
    from .techniques import TECHNIQUES
    if strategy.task not in SC.TASKS:
        raise ValueError(f"{strategy.key}: task must be one of {', '.join(SC.TASKS)}")
    unknown = [t for t in strategy.techniques if t not in TECHNIQUES]
    if not strategy.techniques or unknown:
        raise ValueError(f"{strategy.key}: list its techniques from starplast.techniques"
                         + (f"; unknown: {', '.join(unknown)}" if unknown else ""))
    bad = [p.name for p in strategy.params if p.kind not in PARAM_KINDS]
    if bad:
        raise ValueError(f"{strategy.key}: unknown parameter kind for {', '.join(bad)}")
    REGISTRY[strategy.key] = strategy
    return strategy


def catalog() -> list:
    """Every strategy, in catalogue order."""
    from . import strategy_catalog  # noqa: F401 -- registers on import
    from . import strategy_graph  # noqa: F401 -- registers on import
    from . import strategy_learning  # noqa: F401 -- registers on import
    return sorted(REGISTRY.values(), key=lambda s: s.number)


def get(key: str) -> Strategy:
    """One strategy by key, or a KeyError naming the ones that exist."""
    catalog()
    if key not in REGISTRY:
        raise KeyError(f"no strategy {key!r}; there are {', '.join(sorted(REGISTRY))}")
    return REGISTRY[key]


def families() -> list:
    """The families, in the order their first strategy appears."""
    seen = []
    for s in catalog():
        if s.family not in seen:
            seen.append(s.family)
    return seen


def overview(organism: str = "Tg") -> pd.DataFrame:
    """Every strategy in one table: number, key, family, name, method, the question it answers,
    and -- where the calibration sweep measured it -- its grade and skill at defaults and tuned."""
    from . import calibration as C
    rows = []
    for s in catalog():
        e = C.entry(s.key, organism) or {}
        d, t = e.get("default") or {}, e.get("tuned") or {}
        rows.append({"number": s.number, "key": s.key, "family": s.family, "name": s.name,
                     "method": s.method, "task": s.task, "question": s.question, "cost": s.cost,
                     "grade": e.get("grade", ""), "skill_default": d.get("skill"),
                     "skill_tuned": t.get("skill")})
    return pd.DataFrame(rows)


def metrics(task: str | None = None) -> pd.DataFrame:
    """Every scorecard metric (or one task's), with its definition, range, chance level and how to
    read it -- the glossary the Strategies tab shows beside each number."""
    return SC.glossary(task)


def techniques() -> pd.DataFrame:
    """Every technique the strategies are built from: what it does and why a strategy uses it."""
    from . import techniques as TQ
    return TQ.glossary()


def calibration(key: str, organism: str = "Tg") -> dict:
    """What the calibration sweep measured for `key`: grade, skill at defaults and tuned with 95%
    intervals, per held-out target and per setting. Empty if it has not been measured."""
    from . import calibration as C
    return dict(C.entry(key, organism) or {})


def tuned(key: str, organism: str = "Tg") -> dict:
    """The setting calibration found best for `key` -- pass it on: ``run(key, **tuned(key))``."""
    from . import calibration as C
    return C.tuned_settings(key, organism)


_SHIPPED: dict = {}


def shipped(organism: str = "Tg") -> "Context":
    """The packaged table of one organism, loaded once and kept, so repeated calls share maps."""
    if organism not in _SHIPPED:
        _SHIPPED[organism] = Context.shipped(organism)
    return _SHIPPED[organism]


def run(key: str, organism: str = "Tg", ctx: "Context | None" = None, **settings) -> "StrategyResult":
    """Run strategy `key` on the shipped table of `organism` (or on `ctx`) with `settings`.

    The one-line entry point: ``strategies.run("geneset_hunt", genes=my_list)``.
    """
    return get(key).run(ctx or shipped(organism), **settings)


def test(key: str, organism: str = "Tg", ctx: "Context | None" = None, **settings) -> "TestResult":
    """Self-test strategy `key` with `settings`: hide what is known, ask for it back, score it."""
    return get(key).test(ctx or shipped(organism), **settings)


# --------------------------------------------------------------------------- the context
def _guess_organism(gene_ids) -> str:
    ids = [str(g) for g in list(gene_ids)[:200]]
    return "Pf" if ids and sum(g.upper().startswith(("PF3D7", "PF")) for g in ids) > len(ids) / 2 \
        else "Tg"


class Context:
    """One organism's table and graph, with the leakage guard and the caches strategies share.

    Built once and handed to every strategy, so a map one strategy embedded is not embedded again by
    the next, and a closure computed for `compartment` is computed once. `log` is called with
    progress; `should_stop` is asked between expensive steps and a stop unwinds with `jobs.Stopped`,
    which the job runner reports as stopped rather than failed.

    `graph` is a mapping of `layer__a`, `layer__b`, `layer__w` (and optionally `layer__r`) arrays,
    positions into `nodes`; ``None`` loads the shipped graph for the organism, and only if it was
    built over this exact table. `other` is the other organism's context (or table), for strategies
    that cross species; ``None`` loads the shipped one on first use.
    """

    def __init__(self, nodes: pd.DataFrame, graph=None, other=None, organism: str | None = None,
                 seed: int = 42, log=None, should_stop=None):
        """Wrap a gene table; nothing expensive happens until a strategy asks for it."""
        self.nodes = nodes.reset_index(drop=True)
        self.n = len(self.nodes)
        self.gene_ids = (self.nodes["gene_id"].astype(str).to_numpy() if "gene_id" in self.nodes
                         else np.array([str(i) for i in range(self.n)]))
        self.organism = organism or _guess_organism(self.gene_ids)
        self.seed = int(seed)
        self.log = log or _quiet
        self.should_stop = should_stop
        self._graph = graph
        self._other = other
        self._cache: dict = {}

    @classmethod
    def shipped(cls, organism: str = "Tg", **kw) -> "Context":
        """The packaged table for one organism, as the application loads it."""
        from . import paths
        name = "pf_nodes.parquet" if organism == "Pf" else "nodes.parquet"
        return cls(pd.read_parquet(paths.cache_file(name)), organism=organism, **kw)

    def bound(self, log=None, should_stop=None) -> "Context":
        """The same table, graph and caches, with a job's own progress log and stop flag.

        Two jobs running at once must not report through each other's log or stop each other, and
        they should still share the maps either one has already built -- so the caches are shared
        and only the two callbacks are the job's own.
        """
        import copy
        _ = self.graph                      # loaded once, here, so every bound copy shares it
        out = copy.copy(self)
        out.log = log or _quiet
        out.should_stop = should_stop
        return out

    # ---------------------------------------------------------------- control
    def check(self):
        """Unwind with `jobs.Stopped` if a stop was asked for."""
        if self.should_stop is not None and self.should_stop():
            raise Stopped("strategy stopped")

    def say(self, message: str):
        """Report progress, honouring a stop first."""
        self.check()
        self.log(message)

    def rng(self, salt: int = 0) -> np.random.Generator:
        """A generator seeded from the context's seed, so every run is repeatable."""
        return np.random.default_rng(self.seed + int(salt))

    # ---------------------------------------------------------------- identity
    @property
    def index(self) -> dict:
        """gene id -> position, case-insensitive."""
        if "index" not in self._cache:
            self._cache["index"] = {g.upper(): i for i, g in enumerate(self.gene_ids)}
        return self._cache["index"]

    def resolve_genes(self, text) -> tuple:
        """(positions, unresolved): a gene list, as typed or pasted, matched to this table.

        Separators are anything that is not part of an identifier, so a column pasted from a
        spreadsheet, a comma list and a space list all work. Matching is exact on the accession,
        ignoring case -- a partial match would silently answer for the wrong gene.
        """
        if isinstance(text, (list, tuple, np.ndarray, pd.Series)):
            tokens = [str(t).strip() for t in text]
        else:
            tokens = re.split(r"[^A-Za-z0-9_.\-]+", str(text or ""))
        found, missing = [], []
        for t in tokens:
            if not t:
                continue
            i = self.index.get(t.upper())
            (found if i is not None else missing).append(i if i is not None else t)
        return np.array(sorted(set(found)), dtype=int), sorted(set(missing))

    def product(self, positions) -> list:
        """Product descriptions for some genes, blank where the table has none."""
        if "product" not in self.nodes:
            return [""] * len(positions)
        from .embedding import as_text
        return as_text(self.nodes["product"]).iloc[np.asarray(positions, dtype=int)].tolist()

    def groups(self) -> np.ndarray:
        """A group id per gene -- its orthogroup, so paralogs are hidden together, never split.

        A gene whose paralog stays visible is recovered by resemblance to its own copy, which is not
        the inference being tested. Genes with no orthogroup are their own group.
        """
        if "groups" not in self._cache:
            if "orthogroup" in self.nodes:
                from .embedding import as_text
                og = as_text(self.nodes["orthogroup"]).to_numpy()
                own = np.array([f"__{i}" for i in range(self.n)])
                self._cache["groups"] = np.where(og != "", og, own)
            else:
                self._cache["groups"] = np.array([f"__{i}" for i in range(self.n)])
        return self._cache["groups"]

    # ---------------------------------------------------------------- columns
    def truth(self, column: str) -> pd.Series:
        """A column as labels: plain strings, with every way of saying 'not measured' as NaN."""
        from .embedding import as_text
        if column not in self.nodes:
            raise ValueError(f"{column!r} is not a column of this table")
        s = as_text(self.nodes[column]).str.strip()
        s = s.where(~s.str.lower().isin(ABSENT))
        return s.astype(object).where(s.notna(), np.nan)

    def values(self, column: str) -> pd.Series:
        """A column as numbers, NaN where unmeasured. Booleans are 0/1."""
        s = self.nodes[column]
        if pd.api.types.is_bool_dtype(s):
            return s.astype("float64")
        return pd.to_numeric(s, errors="coerce").astype("float64")

    def categorical_columns(self, min_labelled: int = 30, max_classes: int = 60) -> list:
        """Columns that can be held out as a category: 2..`max_classes` classes, enough labels."""
        if "categorical" not in self._cache:
            out = []
            for c in self.nodes.columns:
                if NEVER_FEATURES.search(c):
                    continue
                s = self.nodes[c]
                if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
                    continue
                t = self.truth(c).dropna()
                if len(t) >= min_labelled and 2 <= t.nunique() <= max_classes:
                    out.append(c)
            self._cache["categorical"] = out
        return list(self._cache["categorical"])

    def numeric_columns(self, exclude=(), min_values: int = 30) -> list:
        """Measurements a strategy may use as features: numeric, varying, measured often enough."""
        if "numeric" not in self._cache:
            out = []
            floor = max(min_values, int(0.02 * self.n))
            for c in self.nodes.columns:
                if NEVER_FEATURES.search(c):
                    continue
                s = self.nodes[c]
                if not (pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s)):
                    continue
                v = self.values(c)
                if v.notna().sum() >= floor and v.nunique() >= 2:
                    out.append(c)
            self._cache["numeric"] = out
        banned = set(exclude)
        return [c for c in self._cache["numeric"] if c not in banned]

    def numeric_targets(self) -> list:
        """Numeric columns worth predicting: at least 50 values and at least 5 distinct ones."""
        return [c for c in self.numeric_columns() if self.values(c).notna().sum() >= 50
                and self.values(c).nunique() >= 5]

    # ---------------------------------------------------------------- the leakage guard
    def banned(self, target) -> set:
        """Every column a strategy scored against `target` must not see, by the project's closure."""
        if not target:
            return set()
        targets = [target] if isinstance(target, str) else list(target)
        out = set()
        for t in targets:
            key = ("banned", t)
            if key not in self._cache:
                from . import search
                if t in self.nodes:
                    try:
                        self._cache[key] = set(search.excluded_for(self.nodes, t))
                    except Exception as exc:          # a closure that cannot be computed bans more
                        self.log(f"closure for {t}: {type(exc).__name__}: {exc}; banning {t} only")
                        self._cache[key] = {t}
                else:
                    self._cache[key] = {t}
            out |= self._cache[key]
        return out

    def banned_layers(self, target) -> set:
        """Every edge layer a strategy scored against `target` must not traverse."""
        if not target:
            return set()
        from . import search
        targets = [target] if isinstance(target, str) else list(target)
        out = set()
        for t in targets:
            try:
                # The declared family's layers AND every layer built from a column the closure
                # removes: a co-expression layer built from the stage series is the stage series.
                out |= set(search.excluded_layers(self.nodes, t))
            except Exception:                          # an unknown target declares no family
                pass
            # A layer is also banned when the table says it was built from this very column: the
            # `compartment` layer joins genes sharing a compartment, whatever the catalogue says.
            if t in self.layers():
                out.add(t)
        return out

    def features(self, target=None, exclude=(), columns=None) -> tuple:
        """(matrix, column names): the permitted measurements, rank-scaled and median-filled."""
        cols = list(columns) if columns is not None else self.numeric_columns(
            set(exclude) | self.banned(target))
        return self.matrix(cols), cols

    def matrix(self, columns, impute: bool = True) -> np.ndarray:
        """Columns rank-scaled to [-0.5, 0.5]; missing values at the median (0) unless not `impute`.

        Ranks rather than z-scores because the screens here have inverted sign conventions, spreads
        that differ 64-fold and heavy tails -- the same reason `embedding` defaults to rank.
        """
        key = ("matrix", tuple(columns), impute)
        if key not in self._cache:
            if not columns:
                X = np.zeros((self.n, 0))
            else:
                frame = pd.DataFrame({c: self.values(c) for c in columns})
                X = (frame.rank(pct=True) - 0.5).to_numpy(dtype=float)
                # Copied before the fill: a pandas buffer can arrive read-only, and writing into one
                # in place is a bug that has shipped here three times already.
                X = np.array(X, dtype=float, copy=True)
                if impute:
                    X[~np.isfinite(X)] = 0.0
            self._cache[key] = X
        return self._cache[key]

    def blocks(self, target=None) -> dict:
        """block name -> columns, for the blocks a map scored against `target` may be built from.

        The project's own blocks where the table carries them, dropped whole when any column is
        banned -- the rule `search` applies, since a block's other columns came from the same
        experiment. A table the slot catalogue does not describe is grouped by column prefix.
        """
        key = ("blocks", target)
        if key not in self._cache:
            banned = self.banned(target)
            out = {}
            try:
                from .embedding import EmbeddingSpec, columns_for, default_spec
                spec = default_spec(self.nodes)
                for b, cols in columns_for(self.nodes, EmbeddingSpec(blocks=spec.blocks)).items():
                    cols = [c for c in cols if not NEVER_FEATURES.search(c)]
                    if cols and not set(cols) & banned:
                        out[b] = cols
            except Exception:                          # a table no slot describes
                out = {}
            if sum(len(v) for v in out.values()) < 0.5 * len(self.numeric_columns(banned)):
                # Grouped by prefix, and a group with ANY banned column goes whole -- the same rule
                # as above. Dropping only the banned columns once let `cycle_t2`, associated with
                # the held-out phase at 0.7999, carry the phase into a map built to be blind to it.
                groups = {}
                for c in self.numeric_columns():
                    groups.setdefault(c.split("_")[0], []).append(c)
                out = {g: cols for g, cols in groups.items() if not set(cols) & banned}
            self._cache[key] = out
        return dict(self._cache[key])

    def family_of(self, block: str) -> str:
        """The kind of measurement a block is: the slot catalogue's axis, where it has one.

        The axis, not the block's name: `Tg_essentiality_in_a_second_background` and
        `Tg_genetic_interaction_delta_gra17` are knockout screens like `Tg_fitness_hff_in_vitro`,
        and grouping by the word in the name let them stand in for the fitness screens a strategy
        had just left out -- predicting fibroblast fitness at 0.89 from "other" evidence that was a
        second knockout screen. A table the catalogue does not describe falls back to the prefix.
        """
        axes = self._cache.get("axes")
        if axes is None:
            try:
                from . import slots
                axes = {sl.key: sl.axis for sl in slots.all_slots(self.organism)
                        if getattr(sl, "axis", None)}
            except Exception:                       # a table no slot describes
                axes = {}
            self._cache["axes"] = axes
        if block in axes:
            return axes[block]
        parts = block.split("_")
        return parts[1] if len(parts) > 1 and parts[0] in ("Tg", "Pf") else parts[0]

    def families(self, target=None) -> dict:
        """family -> blocks, so a walk can ask which KIND of measurement carries a signal."""
        out = {}
        for b in self.blocks(target):
            out.setdefault(self.family_of(b), []).append(b)
        return out

    # ---------------------------------------------------------------- the graph
    @property
    def graph(self) -> dict:
        """The edge layers, as arrays keyed `layer__a` and so on. Empty when there is none."""
        if self._graph is None:
            from . import paths
            name = "pf_graph.npz" if self.organism == "Pf" else "graph.npz"
            path = paths.cache_file(name)
            g = {}
            if os.path.exists(path):
                z = np.load(path, allow_pickle=True)
                ids = z["gene_ids"].astype(str) if "gene_ids" in z.files else None
                if ids is not None and len(ids) == self.n and (ids == self.gene_ids).all():
                    g = {k: z[k] for k in z.files if k.rsplit("__", 1)[-1] in ("a", "b", "w", "r")}
                else:
                    # Edge endpoints are positions, so a graph over a different table is a graph
                    # about different genes. Refused, loudly, rather than used.
                    self.log(f"{name} was built over a different gene table; no edge layer is used")
            self._graph = g
        return self._graph

    def layers(self) -> list:
        """The edge layers that exist and carry at least one edge."""
        g = self.graph
        return sorted(k[:-3] for k in g if k.endswith("__a") and len(g[k]))

    def measurement_layers(self, target=None) -> list:
        """Layers that are measurements rather than literature or derived, minus the banned ones."""
        banned = self.banned_layers(target)
        return [l for l in self.layers() if l not in LITERATURE_LAYERS
                and l not in DERIVED_LAYERS and l not in banned]

    def adjacency(self, layer: str, weight: str = "w"):
        """One layer as a symmetric sparse matrix of raw weights (`w`, `r`, or `binary`)."""
        key = ("adj", layer, weight)
        if key not in self._cache:
            import scipy.sparse as sp
            g = self.graph
            a, b = np.asarray(g[f"{layer}__a"], int), np.asarray(g[f"{layer}__b"], int)
            if weight == "binary":
                w = np.ones(len(a))
            else:
                w = np.asarray(g.get(f"{layer}__{weight}", g.get(f"{layer}__w", np.ones(len(a)))),
                               dtype=float)
            W = sp.coo_matrix((w, (a, b)), shape=(self.n, self.n)).tocsr()
            W = W.maximum(W.T).tocsr()
            W.setdiag(0)
            W.eliminate_zeros()
            self._cache[key] = W
        return self._cache[key]

    def operator(self, layer: str):
        """One layer degree-normalised for a random walk, as `methods.layer_matrix` builds it."""
        key = ("op", layer)
        if key not in self._cache:
            from . import methods
            self._cache[key] = methods.layer_matrix(layer, self.n, graph=self.graph)
        return self._cache[key]

    def other(self) -> "Context | None":
        """The other organism, for strategies that cross species. None when it is not available."""
        if isinstance(self._other, pd.DataFrame):
            self._other = Context(self._other, graph={}, seed=self.seed, log=self.log)
        if self._other is None:
            try:
                self._other = Context.shipped("Pf" if self.organism == "Tg" else "Tg",
                                              graph={}, seed=self.seed, log=self.log)
            except Exception as exc:                   # no second arm in this checkout
                self.log(f"no second organism: {type(exc).__name__}: {exc}")
                self._other = False
        return self._other or None

    # ---------------------------------------------------------------- maps
    def sample_rows(self, size: int | None, target=None, must=()) -> np.ndarray:
        """Positions to embed: all of them, or `size` of them favouring the ones that can be scored.

        `must` are always included -- a gene list's members, say. Labelled genes of `target` fill up
        to 70% of the rest, so a sample always carries something to score against; the remainder is
        unlabelled genes, because a map of labelled genes only is a map of the studied proteome.
        """
        if not size or size >= self.n:
            return np.arange(self.n)
        rng = self.rng(11)
        must = np.unique(np.asarray(must, dtype=int))
        rest = np.setdiff1d(np.arange(self.n), must)
        labelled = np.array([], dtype=int)
        if target:
            known = self.truth(target).notna().to_numpy() if target in self.nodes else np.zeros(
                self.n, bool)
            labelled = rng.permutation(np.intersect1d(rest, np.flatnonzero(known)))
        room = max(int(size) - len(must), 0)
        take = list(labelled[:int(0.7 * room)])
        others = rng.permutation(np.setdiff1d(rest, take))
        take += list(others[:room - len(take)])
        return np.sort(np.concatenate([must, np.asarray(take, dtype=int)])).astype(int)

    def embed(self, blocks, n_neighbors: int = 25, min_dist: float = 0.1, rows=None,
              target=None) -> tuple:
        """(coords, positions): a 3-D map of `rows` built from `blocks`, by the project's builder.

        Cached on everything that changes the answer. `positions` is which genes the map places --
        all of `rows` unless the builder dropped one with no usable measurement.
        """
        from .embedding import EmbeddingSpec, SLOT_BLOCKS, BLOCKS, embed
        rows = np.arange(self.n) if rows is None else np.asarray(rows, dtype=int)
        blocks = tuple(blocks)
        key = ("embed", blocks, int(n_neighbors), float(min_dist), rows.tobytes(), target)
        if key not in self._cache:
            self.check()
            known = self.blocks(target)
            real = tuple(b for b in blocks if b in SLOT_BLOCKS or b in BLOCKS)
            extra = tuple(c for b in blocks if b not in real for c in known.get(b, [b]))
            sub = self.nodes.iloc[rows].reset_index(drop=True)
            spec = EmbeddingSpec(name="strategy", blocks=real, extra_columns=extra,
                                 na_policy="median", scaling="rank",
                                 n_neighbors=int(max(2, min(n_neighbors, len(rows) - 1))),
                                 min_dist=float(min_dist), random_state=self.seed)
            out = embed(sub, spec, log=_quiet)
            coords, kept = np.asarray(out[0], dtype=float), np.asarray(out[2])
            if kept.dtype == bool:
                kept = np.flatnonzero(kept)
            self._cache[key] = (coords, rows[kept.astype(int)])
        return self._cache[key]

    def cluster(self, coords, min_cluster_size: int = 25, selection: str = "eom") -> np.ndarray:
        """HDBSCAN labels for a map, -1 for noise, as every clustering here is made.

        `selection` is HDBSCAN's: "eom" keeps the most persistent clusters, "leaf" every leaf of
        the cluster tree. On the real proteome "eom" returns two to six clusters, the largest
        holding most of the map -- a structure no label can dominate -- while "leaf" returns dozens
        of small, purer ones and calls the rest noise. Strategies that name clusters need "leaf";
        a walk tries both.
        """
        from .clustering import cluster
        self.check()
        mcs = int(max(5, min(min_cluster_size, max(len(coords) // 4, 5))))
        return np.asarray(cluster(np.asarray(coords), "hdbscan", min_cluster_size=mcs,
                                  cluster_selection_method=str(selection)), dtype=int)

    def blind_map(self, target=None, size: int | None = 3000, n_neighbors: int = 25,
                  min_dist: float = 0.1, min_cluster_size: int = 25, must=(),
                  selection: str = "leaf") -> tuple:
        """(coords, positions, labels): one map built without looking at `target` at all."""
        rows = self.sample_rows(size, target=target, must=must)
        blocks = tuple(self.blocks(target))
        coords, pos = self.embed(blocks, n_neighbors, min_dist, rows, target=target)
        return coords, pos, self.cluster(coords, min_cluster_size, selection)

    def same_kind(self, column: str) -> set:
        """Every column measured the same WAY as `column` -- its whole family of blocks.

        Not a leak in the closure's sense: fibroblast fitness and macrophage fitness are different
        experiments. But predicting one knockout screen from another mostly says the screens agree,
        and strategies that ask what OTHER evidence says leave the whole family out.
        """
        blocks = self.blocks()
        fam = next((self.family_of(b) for b, cols in blocks.items() if column in cols), None)
        if fam is None:
            return set()
        return {c for b, cols in blocks.items() if self.family_of(b) == fam for c in cols}


# --------------------------------------------------------------------------- hiding things
def hide(truth: pd.Series, frac: float = 0.25, seed: int = 0, groups=None,
         min_class: int = MIN_CLASS) -> tuple:
    """(visible labels, hidden positions): a share of the labelled genes lose their label.

    Classes smaller than `min_class` are dropped from both sides -- they cannot be split into a part
    to learn from and a part to test on. With `groups`, whole groups are hidden together, so a gene
    is never recovered by resemblance to a paralog left visible.
    """
    t = pd.Series(truth).reset_index(drop=True).astype(object)
    counts = t.value_counts()
    keep = set(counts[counts >= min_class].index)
    t = t.where(t.isin(keep), np.nan)
    labelled = np.flatnonzero(t.notna().to_numpy())
    rng = np.random.default_rng(seed)
    if groups is None:
        hidden = []
        for c in sorted(keep, key=str):
            members = np.flatnonzero((t == c).to_numpy())
            k = int(round(frac * len(members)))
            hidden += list(rng.choice(members, size=k, replace=False)) if k else []
        hidden = np.array(sorted(hidden), dtype=int)
    else:
        g = np.asarray(groups)[labelled]
        order = rng.permutation(np.unique(g))
        target = int(round(frac * len(labelled)))
        chosen, total = set(), 0
        sizes = pd.Series(g).value_counts()
        for grp in order:
            if total >= target:
                break
            chosen.add(grp)
            total += int(sizes[grp])
        hidden = np.sort(labelled[np.isin(g, list(chosen))]).astype(int)
    visible = t.copy()
    visible.iloc[hidden] = np.nan
    return visible, hidden


def shuffled(visible: pd.Series, rng: np.random.Generator) -> pd.Series:
    """The same labels on the same genes, dealt out at random: the null every label test uses."""
    out = visible.copy()
    known = np.flatnonzero(out.notna().to_numpy())
    out.iloc[known] = out.iloc[rng.permutation(known)].to_numpy()
    return out


def correct_rate(pred: pd.Series, truth: pd.Series, positions) -> float:
    """Share of `positions` called correctly. An abstention (NaN) counts as not correct.

    Deliberately not accuracy-over-calls, which a strategy can raise by abstaining on anything hard;
    counting abstentions as misses makes coverage part of the score.
    """
    positions = np.asarray(positions, dtype=int)
    if not len(positions):
        return float("nan")
    p = pd.Series(pred).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    t = pd.Series(truth).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    return float(np.mean([a == b for a, b in zip(p, t)]))


def call_precision(pred: pd.Series, truth: pd.Series, positions) -> tuple:
    """(precision over calls, calls made) on `positions` -- for strategies that trade reach for trust."""
    positions = np.asarray(positions, dtype=int)
    p = pd.Series(pred).reset_index(drop=True).iloc[positions]
    t = pd.Series(truth).reset_index(drop=True).iloc[positions]
    called = p.notna().to_numpy()
    if not called.any():
        return float("nan"), 0
    return float((p[called].to_numpy(dtype=object) == t[called].to_numpy(dtype=object)).mean()), \
        int(called.sum())


def per_class(pred: pd.Series, truth: pd.Series, positions) -> pd.DataFrame:
    """Per hidden class: how many, how many recovered, and how precise the calls of that class were."""
    positions = np.asarray(positions, dtype=int)
    p = pd.Series(pred).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    t = pd.Series(truth).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    rows = []
    for c in sorted({x for x in t if isinstance(x, str)}):
        is_c = t == c
        called_c = p == c
        rows.append({"class": c, "hidden": int(is_c.sum()),
                     "recovered": int((is_c & called_c).sum()),
                     "recall": float((is_c & called_c).sum() / max(is_c.sum(), 1)),
                     "called": int(called_c.sum()),
                     "precision": float((is_c & called_c).sum() / called_c.sum())
                     if called_c.sum() else float("nan")})
    return pd.DataFrame(rows)


def analytic_null(pred: pd.Series, truth: pd.Series, positions) -> tuple:
    """(mean, sd) of `correct_rate` if predictions were independent of the truth.

    Cohen's chance term: the product of the predicted and true class shares, summed. Used where
    re-running a costly strategy on shuffled labels ten times would take longer than it tells.
    """
    positions = np.asarray(positions, dtype=int)
    p = pd.Series(pred).reset_index(drop=True).iloc[positions]
    t = pd.Series(truth).reset_index(drop=True).iloc[positions]
    n = len(positions)
    if not n:
        return float("nan"), float("nan")
    ps = p.value_counts(dropna=True) / n
    ts = t.value_counts(dropna=True) / n
    mean = float(sum(ps.get(c, 0.0) * ts.get(c, 0.0) for c in ps.index))
    return mean, float(math.sqrt(max(mean * (1 - mean), 1e-12) / n))


# --------------------------------------------------------------------------- inference primitives
class _ClassScores:
    """Per-class scores riding on a prediction. pandas deep-copies `attrs` on every operation, so
    the holder copies as itself: the scores are never modified after they are attached."""
    __slots__ = ("frame",)

    def __init__(self, frame: pd.DataFrame):
        self.frame = frame

    def __deepcopy__(self, memo):
        return self


def with_class_scores(pred: pd.Series, scores: pd.DataFrame) -> pd.Series:
    """Attach per-class scores (genes x classes, higher = more likely) to a label prediction.

    The label test reads them back for the scorecard's macro AUROC and AUPRC. A prediction without
    them is still scored on everything else; those two are reported missing, not estimated.
    """
    pred.attrs["class_scores"] = _ClassScores(scores)
    return pred


def class_scores_of(pred) -> pd.DataFrame | None:
    """The per-class scores attached by :func:`with_class_scores`, or None."""
    holder = getattr(pred, "attrs", {}).get("class_scores")
    return holder.frame if isinstance(holder, _ClassScores) else None


def expand_prediction(pred, positions, n: int) -> pd.Series:
    """A prediction over a subset of genes (a map's sample) spread over all `n`, scores included."""
    positions = np.asarray(positions, dtype=int)
    out = pd.Series([np.nan] * n, dtype=object)
    out.iloc[positions] = pd.Series(pred).to_numpy()
    scores = class_scores_of(pred)
    if scores is not None:
        full = pd.DataFrame(np.nan, index=range(n), columns=list(scores.columns))
        full.iloc[positions] = scores.to_numpy()
        with_class_scores(out, full)
    return out


def _score_frame(n: int, query, matrix, classes) -> pd.DataFrame:
    out = np.full((n, len(classes)), np.nan)
    if len(query):
        out[np.asarray(query, dtype=int)] = matrix
    return pd.DataFrame(out, columns=list(classes))


def knn_vote(X: np.ndarray, visible: pd.Series, k: int = 15, query=None) -> tuple:
    """(prediction, vote share) from the `k` nearest labelled genes, distance-weighted.

    A labelled gene is never its own neighbour. Genes are compared on the rows of `X` given, which is
    whatever space the strategy chose -- measurements, or a map.
    """
    from sklearn.neighbors import NearestNeighbors
    vis = pd.Series(visible).reset_index(drop=True)
    known = np.flatnonzero(vis.notna().to_numpy())
    n = len(vis)
    query = np.arange(n) if query is None else np.asarray(query, dtype=int)
    pred = pd.Series([np.nan] * n, dtype=object)
    share = pd.Series(np.nan, index=range(n))
    if len(known) < 2 or not len(query) or X.shape[1] == 0:
        return pred, share
    kk = int(min(k + 1, len(known)))
    nn = NearestNeighbors(n_neighbors=kk).fit(X[known])
    dist, idx = nn.kneighbors(X[query])
    labels = vis.iloc[known].to_numpy(dtype=object)
    classes = sorted(set(labels), key=str)
    col = {c: j for j, c in enumerate(classes)}
    matrix = np.zeros((len(query), len(classes)))
    for row, (q, d, i) in enumerate(zip(query, dist, idx)):
        # At least two labelled genes are consulted, so dropping the query itself never empties this.
        keep = known[i] != q
        d, i = d[keep][:k], i[keep][:k]
        w = 1.0 / (d + 1e-6)
        votes: dict = {}
        for lab, wt in zip(labels[i], w):
            votes[lab] = votes.get(lab, 0.0) + wt
        best = max(votes, key=votes.get)
        total = sum(votes.values())
        pred.iloc[q] = best
        share.iloc[q] = votes[best] / total
        for lab, wt in votes.items():
            matrix[row, col[lab]] = wt / total
    return with_class_scores(pred, _score_frame(n, query, matrix, classes)), share


def onehot(visible: pd.Series) -> tuple:
    """(matrix genes x classes, class names) for the labelled genes; unlabelled rows are zero."""
    import scipy.sparse as sp
    vis = pd.Series(visible).reset_index(drop=True)
    classes = sorted(vis.dropna().unique(), key=str)
    code = {c: i for i, c in enumerate(classes)}
    known = np.flatnonzero(vis.notna().to_numpy())
    cols = np.array([code[v] for v in vis.iloc[known]], dtype=int)
    M = sp.csr_matrix((np.ones(len(known)), (known, cols)), shape=(len(vis), len(classes)))
    return M, classes


def graph_vote(A, visible: pd.Series, query=None, min_support: float = 0.0) -> tuple:
    """(prediction, support) from labelled neighbours in one weighted graph, one hop.

    Abstains where a gene has no labelled neighbour: no edge is no evidence, and a guess from class
    frequencies would be a guess dressed as a network result.
    """
    vis = pd.Series(visible).reset_index(drop=True)
    n = len(vis)
    query = np.arange(n) if query is None else np.asarray(query, dtype=int)
    pred = pd.Series([np.nan] * n, dtype=object)
    support = pd.Series(np.nan, index=range(n))
    M, classes = onehot(vis)
    if not classes or not len(query):
        return pred, support
    counts = np.asarray((A[query] @ M).todense())
    total = counts.sum(axis=1)
    best = counts.argmax(axis=1)
    for row, q in enumerate(query):
        if total[row] > min_support and counts[row, best[row]] > 0:
            pred.iloc[q] = classes[best[row]]
            support.iloc[q] = float(counts[row, best[row]] / total[row])
    shares = np.where(total[:, None] > 0, counts / np.maximum(total[:, None], 1e-12), np.nan)
    return with_class_scores(pred, _score_frame(n, query, shares, classes)), support


def propagate(operator, visible: pd.Series, query=None, restart: float = 0.5,
              iterations: int = 30) -> tuple:
    """(prediction, field strength): random walk with restart, one field per class.

    The fields are seeded from visible genes only, so a hidden gene's own label never seeds its own
    score. Abstains where no class reaches the gene at all.
    """
    vis = pd.Series(visible).reset_index(drop=True)
    n = len(vis)
    query = np.arange(n) if query is None else np.asarray(query, dtype=int)
    M, classes = onehot(vis)
    pred = pd.Series([np.nan] * n, dtype=object)
    strength = pd.Series(np.nan, index=range(n))
    if not classes:
        return pred, strength
    seed = np.asarray(M.todense(), dtype=float)
    seed = seed / np.maximum(seed.sum(axis=0, keepdims=True), 1.0)
    P = seed.copy()
    for _ in range(int(iterations)):
        P = (1 - restart) * (operator @ P) + restart * seed
    F = P[query]
    best = F.argmax(axis=1)
    top = F[np.arange(len(query)), best]
    for row, q in enumerate(query):
        if top[row] > 0:
            pred.iloc[q] = classes[best[row]]
            strength.iloc[q] = float(top[row] / max(F[row].sum(), 1e-12))
    rows = F.sum(axis=1, keepdims=True)
    fields = np.where(rows > 0, F / np.maximum(rows, 1e-12), np.nan)
    return with_class_scores(pred, _score_frame(n, query, fields, classes)), strength


def diffuse(operator, seeds, restart: float = 0.3, iterations: int = 30) -> np.ndarray:
    """Random walk with restart from a set of seed genes: one score per gene, seeds included."""
    n = operator.shape[0]
    s = np.zeros(n)
    s[np.asarray(seeds, dtype=int)] = 1.0 / max(len(seeds), 1)
    p = s.copy()
    for _ in range(int(iterations)):
        p = (1 - restart) * (operator @ p) + restart * s
    return p


def knn_operator(X: np.ndarray, k: int = 15):
    """A symmetric, degree-normalised k-nearest-neighbour graph over the rows of `X`.

    What lets a measurement matrix be walked like an edge layer, so a diffusion can use the table and
    the measured networks in one operator.
    """
    import scipy.sparse as sp
    from sklearn.neighbors import NearestNeighbors
    n = X.shape[0]
    kk = int(min(k + 1, n))
    dist, idx = NearestNeighbors(n_neighbors=kk).fit(X).kneighbors(X)
    rows = np.repeat(np.arange(n), kk)
    W = sp.csr_matrix((np.ones(n * kk), (rows, idx.ravel())), shape=(n, n))
    W = W.maximum(W.T).tocsr()
    W.setdiag(0)
    W.eliminate_zeros()
    d = np.asarray(W.sum(axis=1)).ravel()
    inv = np.divide(1.0, np.sqrt(d), out=np.zeros_like(d), where=d > 0)
    D = sp.diags(inv)
    return D @ W @ D


def cluster_majority(labels: np.ndarray, visible: pd.Series, min_share: float = 0.5,
                     min_labelled: int = 5) -> dict:
    """cluster -> (label, share, labelled count) where one visible label dominates the cluster."""
    labels = np.asarray(labels)
    vis = pd.Series(visible).reset_index(drop=True).to_numpy(dtype=object)
    out = {}
    for k in np.unique(labels):
        if k == NOISE:
            continue
        here = [v for v in vis[labels == k] if isinstance(v, str)]
        if len(here) < min_labelled:
            continue
        s = pd.Series(here).value_counts()
        share = float(s.iloc[0] / len(here))
        if share >= min_share:
            out[int(k)] = (str(s.index[0]), share, len(here))
    return out


def predict_from_clusters(labels: np.ndarray, positions, visible: pd.Series, n: int,
                          min_share: float = 0.5, min_labelled: int = 5) -> tuple:
    """(prediction over all genes, the cluster table): every member of a dominated cluster is called."""
    positions = np.asarray(positions, dtype=int)
    vis = pd.Series(visible).reset_index(drop=True)
    sub = vis.iloc[positions].reset_index(drop=True)
    major = cluster_majority(labels, sub, min_share, min_labelled)
    pred = pd.Series([np.nan] * n, dtype=object)
    for row, k in enumerate(np.asarray(labels)):
        if int(k) in major:
            pred.iloc[positions[row]] = major[int(k)][0]
    table = pd.DataFrame([{"cluster": k, "label": v[0], "share": v[1], "labelled": v[2],
                           "size": int((np.asarray(labels) == k).sum())}
                          for k, v in sorted(major.items())])
    return pred, table


def set_f1(labels: np.ndarray, members: np.ndarray) -> tuple:
    """(best cluster, precision, recall, F1) of a gene set against one clustering, noise excluded.

    `members` is a boolean mask over the clustered genes. The best cluster is the one with the
    highest F1 -- the question is whether ONE cluster holds the set, not how many it touches.
    """
    labels = np.asarray(labels)
    members = np.asarray(members, dtype=bool)
    n_set = int(members.sum())
    if n_set == 0:
        return NOISE, 0.0, 0.0, 0.0
    ok = labels != NOISE
    if not ok.any():
        return NOISE, 0.0, 0.0, 0.0
    lab = labels[ok]
    ids, inv = np.unique(lab, return_inverse=True)
    size = np.bincount(inv)
    hits = np.bincount(inv, weights=members[ok].astype(float))
    prec = hits / size
    rec = hits / n_set
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0.0)
    b = int(np.argmax(f1))
    return int(ids[b]), float(prec[b]), float(rec[b]), float(f1[b])


def default_category(ctx: Context) -> str | None:
    """The held-out label a strategy starts on: the measured localization where one exists."""
    cats = ctx.categorical_columns()
    for c in ("compartment", "lopit_unified", "pb_transferred_phenotype", "stage_enriched_derived",
              "dtm_class", "export_pred_tier"):
        if c in cats:
            return c
    return cats[0] if cats else None


def default_numeric(ctx: Context) -> str | None:
    """The numeric trait a regression starts on: the canonical fitness screen where one exists."""
    nums = ctx.numeric_targets()
    for c in ("fit_invitro_hff", "piggybac_mis", "piggybac_mfs"):
        if c in nums:
            return c
    return nums[0] if nums else None


def example_set(ctx: Context, target: str | None, lo: int = 30, hi: int = 300,
                aim: int = 120) -> tuple:
    """(category, positions): a known gene set to test a gene-list strategy on.

    The category of `target` whose size is closest to `aim` within [lo, hi] -- big enough to split
    into a query and a hidden part, small enough to be a set someone would actually hand over.
    """
    if not target or target not in ctx.nodes:
        return None, np.array([], dtype=int)
    t = ctx.truth(target)
    counts = t.value_counts()
    ok = counts[(counts >= lo) & (counts <= hi)]
    if ok.empty:
        return None, np.array([], dtype=int)
    cat = str((ok.index.to_series().map(lambda c: abs(ok[c] - aim))).sort_values(
        kind="stable").index[0])
    return cat, np.flatnonzero((t == cat).to_numpy())


# --------------------------------------------------------------------------- the five test patterns
def label_transfer_test(ctx: Context, key: str, target: str, predict: Callable, *,
                        frac: float = 0.25, n_null: int = 10, min_effect: float = 0.05,
                        restrict=None, null: str = "shuffle", metric: str = None,
                        precision: bool = False, note: str = "", truth=None) -> TestResult:
    """Pattern 1: hide a share of a label, predict it back, and compare with shuffled labels.

    `predict(visible)` returns a prediction over all genes (NaN to abstain) from the visible labels
    alone. `truth` replaces the column's own labels when a strategy scores a derived label -- the
    first digit of an EC number, say. `restrict` is a boolean mask over genes: only hidden genes inside it are scored -- how a
    strategy that can only speak about some genes (those with a crosslink partner, say) is tested
    on exactly those. With `precision`, the metric is precision over calls instead of correct calls
    per hidden gene -- for strategies whose point is to abstain.
    """
    t0 = time.monotonic()
    truth = ctx.truth(target) if truth is None else pd.Series(truth).reset_index(drop=True)
    visible, hidden = hide(truth, frac, ctx.seed, groups=ctx.groups())
    if restrict is not None:
        hidden = hidden[np.asarray(restrict, dtype=bool)[hidden]]
    ctx.say(f"{key}: {len(hidden):,} labelled genes hidden, "
            f"{int(visible.notna().sum()):,} left to learn from")
    if len(hidden) < MIN_HIDDEN:
        return judge(key, metric or "correct calls per hidden gene", float("nan"), [],
                     min_effect=min_effect, n_hidden=len(hidden), hidden=f"{target} labels",
                     null_kind="shuffled labels", t0=t0, task=SC.T_LABEL,
                     note=note or "too few labelled genes this strategy can speak about")

    def score(pred):
        if precision:
            return call_precision(pred, truth, hidden)[0]
        return correct_rate(pred, truth, hidden)

    pred = predict(visible)
    observed = score(pred)
    calls = int(pd.Series(pred).iloc[hidden].notna().sum())
    details = per_class(pred, truth, hidden)
    card = SC.label_calls(pred, truth, hidden, class_scores_of(pred))
    nulls, analytic = [], None
    if null == "analytic" and not precision:
        analytic = analytic_null(pred, truth, hidden)
    else:
        rng = ctx.rng(101)
        for i in range(int(n_null)):
            ctx.say(f"{key}: null {i + 1} of {n_null}")
            v = score(predict(shuffled(visible, rng)))
            # A null that calls nothing has no precision to report, and a strategy that on
            # shuffled labels never finds anything confident enough to call is doing exactly what
            # it should: counted as zero, not as missing, or an ideal null would make every test
            # inconclusive.
            nulls.append(0.0 if precision and not np.isfinite(v) else v)
    # Precision rests on the calls made, so that is what has to be numerous enough to judge.
    scored = calls if precision else len(hidden)
    return judge(key, metric or ("precision of calls on hidden genes" if precision
                                 else "correct calls per hidden gene"),
                 observed, nulls, min_effect=min_effect, n_hidden=scored,
                 hidden=(f"{frac:.0%} of the {target} labels, whole orthogroups at a time"
                         + (f"; {calls:,} calls made on {len(hidden):,} hidden genes"
                            if precision else "")),
                 null_kind=("the analytic chance level" if analytic is not None
                            else f"{n_null} runs on shuffled labels"),
                 t0=t0, details=details, note=note, analytic=analytic,
                 numbers={"calls": calls, "coverage": calls / max(len(hidden), 1)},
                 task=SC.T_LABEL, scorecard=card)


def set_expansion_test(ctx: Context, key: str, members, rank: Callable, *, frac: float = 0.3,
                       n_null: int = 10, min_effect: float = 0.1, universe=None,
                       label: str = "the gene set", note: str = "") -> TestResult:
    """Pattern 2: hand over part of a gene set, and ask whether the rest ranks near the top.

    `rank(query_positions)` scores every gene. The metric is the AUROC of the hidden members against
    every other gene not in the query; the null is the same procedure on random sets of the same
    size drawn from the same `universe` (all genes by default).
    """
    t0 = time.monotonic()
    members = np.unique(np.asarray(members, dtype=int))
    universe = np.arange(ctx.n) if universe is None else np.asarray(universe, dtype=int)
    rng = ctx.rng(202)

    cards = []

    def once(pos):
        pos = rng.permutation(pos)
        k = int(round(frac * len(pos)))
        hidden, query = pos[:k], pos[k:]
        scores = np.asarray(rank(query), dtype=float)
        cand = np.setdiff1d(universe, query)
        cards.append((scores[cand], np.isin(cand, hidden)))
        return auroc(scores[cand], np.isin(cand, hidden)), len(hidden)

    if len(members) < 10:
        return judge(key, "AUROC of hidden members", float("nan"), [], min_effect=min_effect,
                     n_hidden=0, hidden=label, null_kind="random gene sets", t0=t0,
                     note=note or "fewer than ten genes to split", task=SC.T_RANK)
    observed, n_hidden = once(members)
    card = SC.ranking(*cards[0])
    nulls = []
    for i in range(int(n_null)):
        ctx.say(f"{key}: random set {i + 1} of {n_null}")
        nulls.append(once(rng.choice(universe, size=len(members), replace=False))[0])
    return judge(key, "AUROC of hidden members against every other gene", observed, nulls,
                 min_effect=min_effect, n_hidden=n_hidden,
                 hidden=f"{frac:.0%} of {label} ({len(members)} genes)",
                 null_kind=f"{n_null} random sets of the same size", t0=t0, note=note,
                 task=SC.T_RANK, scorecard=card)


def pair_test(key: str, score_fn: Callable, positives: np.ndarray, negatives: np.ndarray,
              null_fn: Callable, *, n_null: int = 5, min_effect: float = 0.05, hidden: str,
              null_kind: str, t0: float, note: str = "", ctx: Context | None = None,
              metric: str = "AUROC of hidden pairs against random non-pairs") -> TestResult:
    """Pattern 3: hidden pairs against non-pairs, scored by AUROC, against a permuted null.

    `score_fn(pairs)` and `null_fn(pairs, i)` return one score per pair (rows of `[a, b]`). `metric`
    names how the non-pairs were drawn, because that decides what the number means: against random
    non-pairs an AUROC mostly measures degree, against degree-matched ones it measures the pair.
    """
    pairs = np.vstack([positives, negatives]) if len(negatives) else positives
    y = np.r_[np.ones(len(positives), bool), np.zeros(len(negatives), bool)]
    scored = np.asarray(score_fn(pairs), dtype=float)
    observed = auroc(scored, y)
    nulls = []
    for i in range(int(n_null)):
        if ctx is not None:
            ctx.say(f"{key}: null {i + 1} of {n_null}")
        nulls.append(auroc(null_fn(pairs, i), y))
    return judge(key, metric, observed, nulls,
                 min_effect=min_effect, n_hidden=len(positives), hidden=hidden,
                 null_kind=null_kind, t0=t0, note=note, task=SC.T_RANK,
                 scorecard=SC.ranking(scored, y))


def value_test(ctx: Context, key: str, y: pd.Series, predict: Callable, *, frac: float = 0.2,
               n_null: int = 3, min_effect: float = 0.1, label: str = "the values",
               note: str = "") -> TestResult:
    """Pattern 4: hide a share of measured values, predict them, and rank-correlate.

    `predict(visible_y)` returns a prediction for every gene from the visible values only. The bar
    is the exact chance distribution of a rank correlation between independent variables -- centred
    on zero with standard deviation 1/sqrt(n - 1) -- because three or five refits on shuffled values
    give a 95th percentile too unstable to judge by: on a table with nothing in it, three refits
    once averaged -0.10 and let a correlation of 0.10 through. The refits are still run and
    reported, as a check that a model trained on nothing predicts nothing.
    """
    t0 = time.monotonic()
    y = pd.Series(y).reset_index(drop=True).astype(float)
    measured = np.flatnonzero(y.notna().to_numpy())
    rng = ctx.rng(303)
    hidden = np.sort(rng.choice(measured, size=int(round(frac * len(measured))), replace=False)) \
        if len(measured) else np.array([], dtype=int)
    if len(hidden) < MIN_HIDDEN:
        return judge(key, "rank correlation on hidden values", float("nan"), [],
                     min_effect=min_effect, n_hidden=len(hidden), hidden=label,
                     null_kind="shuffled values", t0=t0, note=note or "too few measured values",
                     task=SC.T_VALUES)
    visible = y.copy()
    visible.iloc[hidden] = np.nan
    pred = np.asarray(predict(visible), dtype=float)
    observed = spearman(pred[hidden], y.iloc[hidden])
    nulls = []
    for i in range(int(n_null)):
        ctx.say(f"{key}: null {i + 1} of {n_null}")
        v = visible.copy()
        known = np.flatnonzero(v.notna().to_numpy())
        v.iloc[known] = v.iloc[rng.permutation(known)].to_numpy()
        nulls.append(spearman(np.asarray(predict(v), dtype=float)[hidden], y.iloc[hidden]))
    return judge(key, "rank correlation of predicted and hidden values", observed, nulls,
                 min_effect=min_effect, n_hidden=len(hidden), hidden=f"{frac:.0%} of {label}",
                 null_kind="the chance distribution of a rank correlation", t0=t0, note=note,
                 analytic=(0.0, 1.0 / math.sqrt(max(len(hidden) - 1, 1))),
                 numbers={"shuffled_refits": float(np.nanmean(nulls)) if nulls else float("nan")},
                 task=SC.T_VALUES, scorecard=SC.values(pred[hidden], y.iloc[hidden].to_numpy()))


def replication_test(ctx: Context, key: str, find: Callable, replicate: Callable,
                     scramble: Callable, positions, *, n_null: int = 20, min_effect: float = 0.2,
                     min_findings: int = 3, label: str = "findings", note: str = "") -> TestResult:
    """Pattern 5: make findings on half the genes, and ask whether they hold on the other half.

    `find(half_a)` returns findings from genes in `half_a` only; `replicate(finding, half_b, state)`
    says whether one holds on `half_b`; `scramble(half_b, rng)` returns a `state` in which the second
    half's evidence is permuted, and ``None`` means the real one. The metric is the share of
    findings that replicate; the null is that share with the second half scrambled.
    """
    t0 = time.monotonic()
    positions = np.asarray(positions, dtype=int)
    rng = ctx.rng(404)
    perm = rng.permutation(positions)
    half_a, half_b = np.sort(perm[: len(perm) // 2]), np.sort(perm[len(perm) // 2:])
    findings = find(half_a)
    n = len(findings)
    if n < min_findings:
        return judge(key, "share of findings that replicate", float("nan"), [],
                     min_effect=min_effect, n_hidden=n, hidden=f"half of the genes' {label}",
                     null_kind="the second half scrambled", t0=t0, min_hidden=min_findings,
                     note=note or f"only {n} findings on the first half, {min_findings} needed",
                     task=SC.T_REPL)
    rate = lambda state: float(np.mean([bool(replicate(f, half_b, state)) for f in findings]))
    observed = rate(None)
    nulls = []
    for i in range(int(n_null)):
        ctx.check()
        nulls.append(rate(scramble(half_b, rng)))
    # Scored on the findings, so what counts as enough is how many findings there were to check.
    return judge(key, "share of first-half findings that replicate on the second half", observed,
                 nulls, min_effect=min_effect, n_hidden=n, min_hidden=min_findings,
                 hidden=f"{label} made on half the genes, checked on the other half",
                 null_kind=f"{n_null} runs with the second half scrambled", t0=t0,
                 note=note, numbers={"findings": n}, task=SC.T_REPL,
                 scorecard=SC.replication(int(round(observed * n)), n, nulls))


# --------------------------------------------------------------------------- a planted table
def planted_context(n: int = 480, seed: int = 7, null: bool = False) -> Context:
    """A synthetic organism whose answers are known, for testing that each strategy can find them.

    Every signal the strategies look for is planted: localization carried by expression and by a
    crosslink layer, a second label carried by fitness, protein families shared through a structural
    layer, paralogs, a second species with orthologs and correlated fitness, attention that tracks
    fame rather than biology, and a condition screen that depends on a host-facing trait. With
    `null`, the same table is built with every label and every edge dealt out at random, so nothing
    is there to find -- which is what every self-test must then report.

    Column names follow the real tables where a strategy's default reaches for them, so the same
    defaults work here and on the shipped data.
    """
    rng = np.random.default_rng(seed)
    comps = ["nucleus", "cytosol", "mitochondrion", "apicoplast", "micronemes", "dense granules"]
    phases = ["G1", "S", "M", "C"]
    comp = rng.choice(len(comps), size=n)
    # Paralog pairs. Most kept their partner's job and profile; the rest diverged, taking their own
    # compartment and their own measurements -- the difference strategy 28 looks for.
    pairs = [(i, i + 1) for i in range(0, n - 1, 9)]
    kept = rng.random(len(pairs)) < 0.6
    for (i, j), same in zip(pairs, kept):
        if same:
            comp[j] = comp[i]
    # Phase is independent of compartment except in two compartments that subdivide by it: the
    # structure strategies 26 and 27 exist to find.
    phase = rng.choice(len(phases), size=n)
    phase[comp == 0] = rng.choice([0, 1], size=int((comp == 0).sum()))
    phase[comp == 1] = rng.choice([2, 3], size=int((comp == 1).sum()))
    phase[comp == 2] = rng.choice([0, 2], size=int((comp == 2).sum()))
    # A stage label whose features separate particular COMBINATIONS of compartment and stage: the
    # kind of gene strategy 27 looks for, defined by both labels at once.
    stages = ["tachyzoite", "bradyzoite", "oocyst"]
    stage = rng.choice(len(stages), size=n)
    combos = {(5, 1): 0, (3, 2): 1, (4, 0): 2, (2, 1): 3}
    combo_centers = rng.normal(scale=4.0, size=(len(combos), 10))
    stage_centers = rng.normal(scale=0.5, size=(len(stages), 10))
    family = rng.choice(18, size=n)
    essential = rng.normal(size=n)
    host_facing = (comp == 5).astype(float) + 0.3 * rng.normal(size=n)
    frame: dict = {"gene_id": [f"TGSYN_{i:06d}" for i in range(n)]}
    # Thirty moderately informative columns rather than a few strong ones: each stays under the
    # closure's 0.8 association threshold (a column that restates the label is removed as a
    # near-copy, which is right), while together they separate compartments -- as real expression
    # series do.
    n_expr = 30
    centers = rng.normal(scale=1.1, size=(len(comps), n_expr))
    expr = centers[comp] + rng.normal(size=(n, n_expr))
    fam_centers = rng.normal(size=(18, 4))
    prot = fam_centers[family] + 0.8 * rng.normal(size=(n, 4)) + 0.3 * centers[comp][:, :4]
    for (i, j), same in zip(pairs, kept):
        if same:
            expr[j] = expr[i] + 0.4 * rng.normal(size=n_expr)
            prot[j] = prot[i] + 0.4 * rng.normal(size=4)
    for j in range(n_expr):
        frame[f"expr_stage{j}"] = expr[:, j]
    phase_centers = rng.normal(scale=3.0, size=(len(phases), 6))
    cyc = phase_centers[phase] + rng.normal(size=(n, 6))
    for j in range(6):
        frame[f"cycle_t{j}"] = cyc[:, j]
    stg = stage_centers[stage] + rng.normal(size=(n, 10))
    for i in range(n):
        k = combos.get((int(comp[i]), int(stage[i])))
        if k is not None:
            stg[i] += combo_centers[k]
    for j in range(10):
        frame[f"stagefeat_{j}"] = stg[:, j]
    for j, name in enumerate(["fit_invitro_hff", "fit_naive_bmdm", "fit_ifng", "fit_extra"]):
        frame[name] = essential + 0.4 * rng.normal(size=n) + 0.3 * (phase == j % 4)
    frame["fit_invivo_PE"] = frame["fit_invitro_hff"] + 1.2 * host_facing + 0.3 * rng.normal(size=n)
    # Essentiality leaves a trace outside the fitness screens -- abundance, conservation of fold --
    # so a fitness score is predictable from measurements its closure does not remove.
    prot[:, 3] += 0.9 * essential
    for j in range(4):
        frame[f"prot_feat{j}"] = prot[:, j]
    for j in range(3):
        frame[f"noise_n{j}"] = rng.normal(size=n)
    # Missingness, as every real screen has.
    for c in list(frame)[1:]:
        miss = rng.random(n) < 0.08
        frame[c] = np.where(miss, np.nan, frame[c])
    labelled = rng.random(n) < 0.75
    frame["compartment"] = np.where(labelled, np.array(comps)[comp], "unassigned")
    frame["cellcycle_phase"] = np.where(rng.random(n) < 0.8, np.array(phases)[phase], "")
    frame["stage_enriched_derived"] = np.where(rng.random(n) < 0.85, np.array(stages)[stage], "")
    frame["ec_number"] = [f"{1 + f % 6}.{f}.1.1 (synthetic)" if rng.random() < 0.7 else ""
                          for f in family]
    ogs = [f"OG_{i}" for i in range(n)]
    for i, j in pairs:
        ogs[j] = ogs[i]
    frame["orthogroup"] = ogs
    frame["lineage_specific"] = rng.random(n) < 0.35
    frame["product"] = np.where(rng.random(n) < 0.4, "hypothetical protein", "synthetic protein")
    fame = rng.pareto(1.5, size=n)
    frame["n_papers_focal"] = np.floor(fame).astype(int)
    frame["n_papers_substantive"] = np.floor(fame * 0.5).astype(int)
    frame["attention_depth"] = np.where(fame > 1, "focal", np.where(fame > 0.3, "incidental", ""))
    nodes = pd.DataFrame(frame)

    def pairs_within(groups, per_gene: int, weight=None):
        a, b = [], []
        for g in np.unique(groups):
            members = np.flatnonzero(groups == g)
            if len(members) < 2:
                continue
            for i in members:
                for j in rng.choice(members, size=min(per_gene, len(members) - 1), replace=False):
                    if i < j:
                        a.append(i); b.append(j)
        a, b = np.array(a, int), np.array(b, int)
        w = np.ones(len(a)) if weight is None else weight(len(a))
        return a, b, w

    graph: dict = {}
    # Crosslinks inside compartments, weighted by crosslink count.
    a, b, w = pairs_within(comp * 10 + (family % 3), 3, lambda k: rng.integers(1, 5, size=k))
    graph.update({"xlms__a": a, "xlms__b": b, "xlms__w": w.astype(float)})
    # Co-expression within compartments; co-fitness within phases.
    a, b, w = pairs_within(comp, 4, lambda k: rng.uniform(0.6, 1.0, size=k))
    graph.update({"coexpression__a": a, "coexpression__b": b, "coexpression__w": w})
    a, b, w = pairs_within(phase, 3, lambda k: rng.uniform(0.6, 1.0, size=k))
    graph.update({"cofitness__a": a, "cofitness__b": b, "cofitness__w": w})
    # Structural similarity within families.
    a, b, w = pairs_within(family, 3, lambda k: rng.uniform(0.7, 0.95, size=k))
    graph.update({"struct__a": a, "struct__b": b, "struct__w": w})
    # Paralogs.
    og = np.asarray(ogs)
    a, b, w = pairs_within(og, 1)
    graph.update({"orthogroup__a": a, "orthogroup__b": b, "orthogroup__w": w})
    # Literature: raw counts follow fame, the residual follows shared compartment.
    k = 3 * n
    a = rng.integers(0, n, size=k)
    b = rng.integers(0, n, size=k)
    keep = a != b
    a, b = a[keep], b[keep]
    same = comp[a] == comp[b]
    raw = 1 + np.floor(3 * (fame[a] + fame[b]))
    resid = np.where(same, 2.0, 0.2) + 0.2 * rng.normal(size=len(a))
    graph.update({"comention_ft__a": a, "comention_ft__b": b, "comention_ft__w": raw,
                  "comention_ft__r": resid})
    # And the literature mostly writes about pairs the measurements also link.
    lit_a, lit_b = graph["coexpression__a"][::3], graph["coexpression__b"][::3]
    graph.update({"comention__a": lit_a, "comention__b": lit_b,
                  "comention__w": np.ones(len(lit_a)), "comention__r": np.ones(len(lit_a))})
    # The target layer built from the label itself, which the guard must ban.
    a, b, w = pairs_within(np.where(labelled, comp, -1 - np.arange(n)), 5)
    graph.update({"compartment__a": a, "compartment__b": b, "compartment__w": w})

    # A second species: orthologs of half the genes, with fitness that tracks essentiality.
    half = rng.choice(n, size=n // 2, replace=False)
    other = pd.DataFrame({
        "gene_id": [f"PFSYN_{i:06d}" for i in range(len(half))],
        "orthogroup": np.asarray(ogs)[half],
        "piggybac_mis": essential[half] + 0.5 * rng.normal(size=len(half)),
        "localization_other": np.array(comps)[comp[half]],
    })

    if null:
        # Every column dealt out at random, independently, and every edge rewired: no label is
        # carried by any measurement, and no measurement by any other, so nothing is left to find.
        for c in [c for c in nodes.columns if c != "gene_id"]:
            nodes[c] = rng.permutation(nodes[c].to_numpy())
        for key in [k for k in graph if k.endswith("__a")]:
            layer = key[:-3]
            graph[f"{layer}__a"] = rng.integers(0, n, size=len(graph[key]))
            graph[f"{layer}__b"] = rng.integers(0, n, size=len(graph[key]))
        other["piggybac_mis"] = rng.permutation(other["piggybac_mis"].to_numpy())
        other["localization_other"] = rng.permutation(other["localization_other"].to_numpy())
    other_ctx = Context(other, graph={}, organism="Pf", seed=seed)
    return Context(nodes, graph=graph, other=other_ctx, organism="Tg", seed=seed)
