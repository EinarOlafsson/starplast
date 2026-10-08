"""The track record: hold every labelled gene out in turn, and keep what each strategy said.

A self-test hides a quarter of a label once and reports one number. That answers "does this strategy
beat chance on this label", which is what a verdict needs, and it leaves the question a biologist
actually asks unanswered: *could this software have told me what I already know about THIS gene?*

So this module runs the hold-out to completion. The labelled genes are split into folds by
ORTHOGROUP, every fold is hidden in turn, and the strategy is asked for the genes it cannot see --
so each labelled gene is held out exactly once, is never predicted by a model that saw it or its
paralog, and ends with a row of its own saying what was predicted, whether that was right, and how
sure the strategy was. Those rows are the ledger; `summary` pools them at the level being asked
about, from one gene up to one strategy, and every rate carries the count it rests on and a Wilson
interval, because a rate without those is a number nobody can act on.

    ledger = evaluate(ctx, "feature_knn", "compartment")      # one row per labelled gene
    summary(ledger, "gene")                                   # how each gene fared
    summary(ledger, "class")                                  # and each class, with what it is confused with

Abstentions are kept apart from errors throughout: a strategy that declines to answer is not wrong,
and averaging the two together would hide exactly the difference the scorecard exists to show.
"""
from __future__ import annotations

import os
import time
from urllib.parse import quote

import numpy as np
import pandas as pd

from . import organisms
from . import strategies as S

#: Folds the labelled genes are split into. Five is the project's cross-validation default: with
#: fewer, each model sees too little; with more, the run costs more than the extra precision is
#: worth.
FOLDS = 5

#: Below this many scored genes a rate is noise. Reported as a count instead, never as a percentage.
MIN_FOR_RATE = 5

#: What a ledger row holds. `prediction` and `correct` are null where the strategy abstained, which
#: is a third outcome and not a wrong answer.
COLUMNS = ("organism", "strategy", "target", "setting_key", "seed", "fold", "mode",
           "set_name", "set_size", "gene", "gene_id", "truth", "prediction", "correct",
           "abstained", "support")


def folds_by_group(ctx, positions, folds: int = FOLDS, seed: int = 0) -> np.ndarray:
    """A fold number per position, whole orthogroups kept together.

    Groups are dealt round-robin from a shuffled list rather than cut into blocks, so one enormous
    family cannot put a quarter of the genes into one fold.
    """
    positions = np.asarray(positions, dtype=int)
    groups = np.asarray(ctx.groups())[positions]
    rng = np.random.default_rng(seed)
    order = rng.permutation(np.unique(groups))
    sizes = pd.Series(groups).value_counts()
    load = np.zeros(int(folds))
    where = {}
    for group in order:                       # biggest first, into whichever fold is emptiest
        pick = int(np.argmin(load))
        where[group] = pick
        load[pick] += float(sizes[group])
    return np.array([where[g] for g in groups], dtype=int)


def _codes(values) -> pd.Series:
    return pd.Series(values).astype("object").where(pd.notna(pd.Series(values)), None)


def evaluate(ctx, key: str, target: str, settings: dict | None = None, folds: int = FOLDS,
             seed: int = 0, log=None) -> pd.DataFrame:
    """Hold out every labelled gene in turn and record what `key` said about each.

    One row per labelled gene: the fold it was hidden in, the label it carries, what the strategy
    predicted without seeing it or its orthogroup, whether that was right, and the strategy's own
    support for the call. Returns an empty frame, rather than raising, when the strategy cannot
    speak about this target at all -- a refusal is a fact about the pairing, not an error.
    """
    say = log or (lambda _m: None)
    strategy = S.get(key)
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    labelled = np.flatnonzero(kept.notna().to_numpy())
    if len(labelled) < folds * 2:
        return pd.DataFrame(columns=list(COLUMNS))
    fold_of = folds_by_group(ctx, labelled, folds, seed)
    chosen = _permitted(ctx, dict(strategy.settings(ctx, target=target, **(settings or {}))), target)
    setting_key = ", ".join(f"{k}={v}" for k, v in sorted(chosen.items()) if k != "target")
    rows = []
    for fold in range(int(folds)):
        ctx.check()
        hidden = labelled[fold_of == fold]
        if not len(hidden):
            continue
        say(f"{key} on {target}: fold {fold + 1} of {folds}, {len(hidden):,} genes held out")
        visible = kept.copy()
        visible.iloc[hidden] = np.nan
        t0 = time.monotonic()
        pred, support = _predict(ctx, strategy, visible, hidden, chosen, target)
        say(f"  {time.monotonic() - t0:.1f}s")
        if pred is None:
            return pd.DataFrame(columns=list(COLUMNS))
        p = pd.Series(pred).reset_index(drop=True)
        sup = pd.Series(support).reset_index(drop=True) if support is not None else None
        for gene in hidden:
            call = p.iloc[gene]
            said = isinstance(call, str)
            rows.append({
                "organism": ctx.organism, "strategy": key, "target": target,
                "setting_key": setting_key, "seed": int(seed), "fold": int(fold), "mode": "together",
                "set_name": None, "set_size": len(hidden),
                "gene": int(gene), "gene_id": str(ctx.gene_ids[gene]),
                "truth": str(kept.iloc[gene]), "prediction": call if said else None,
                "correct": bool(call == kept.iloc[gene]) if said else None,
                "abstained": not said,
                "support": float(sup.iloc[gene]) if sup is not None and said
                and np.isfinite(sup.iloc[gene]) else float("nan")})
    return pd.DataFrame(rows, columns=list(COLUMNS))


def _permitted(ctx, chosen: dict, target: str) -> dict:
    """The settings with any banned edge layer replaced by the first permitted one.

    A strategy's default layer can be one the target bans -- label diffusion defaults to
    coexpression, which IS the stage label's source. Run in the Strategies tab, the strategy refuses
    it; here, silently walking it once "recovered" the stage label at 100%. The replacement goes into
    the settings, so the ledger's `setting_key` says which layer was actually walked.
    """
    layer = chosen.get("layer")
    if layer and layer in ctx.banned_layers(target):
        allowed = [l for l in ctx.measurement_layers(target) if l != "orthogroup"]
        chosen = {**chosen, "layer": allowed[0] if allowed else None}
    return chosen

def _predict(ctx, strategy, visible: pd.Series, query, settings: dict, target: str):
    """(prediction, support) from one strategy, given the visible labels only.

    The strategies do not share one prediction entry point -- each `runner` returns tables meant for
    a person to read -- so this asks the same primitives their own self-tests ask, chosen by what the
    strategy is built from (`Strategy.techniques`). A strategy whose method cannot be asked about
    specific genes returns `(None, None)` and is skipped rather than guessed at.
    """
    from . import strategy_catalog as C
    techniques = set(strategy.techniques)
    key = strategy.key
    X = None
    if {"knn", "logistic_regression", "sgc", "random_forest", "pu_bagging"} & techniques:
        X, _cols = ctx.features(target)
    if key in ("feature_knn", "stratum_focus"):
        pred, share = S.knn_vote(X, visible, int(settings.get("k", 15)), query=query)
        floor = float(settings.get("min_share", 0.0) or 0.0)
        return pred.where(share >= floor), share
    if key == "supervised_classifier":
        pred, prob, _m = C._logistic(X, visible, query, float(settings.get("C", 0.5)),
                                     float(settings.get("min_probability", 0.0)))
        return pred, prob
    if key == "random_forest":
        from . import strategy_learning as L
        pred, prob, _m = L._forest(X, visible, query, settings.get("trees", 300),
                                   settings.get("min_leaf", 2), ctx.seed)
        return pred, prob
    if key == "graph_convolution":
        from . import strategy_learning as L
        H, _c, _l = L._smoothed(ctx, target, int(settings.get("hops", 2)))
        pred, prob, _m = C._logistic(H, visible, query, float(settings.get("C", 0.5)),
                                     float(settings.get("min_probability", 0.0)))
        return pred, prob
    if key == "stacking":
        from . import strategy_learning as L
        pred, support, _trust = L._stack(ctx, L._stack_bases(ctx, target, int(settings.get("k", 15))),
                                         visible, query, int(settings.get("folds", 5)))
        return pred, support
    if key == "layer_propagation":
        layer = settings.get("layer") or (ctx.measurement_layers(target) or [None])[0]
        if layer is None or layer in ctx.banned_layers(target):
            return None, None
        pred, strength = S.propagate(ctx.operator(layer), visible, query,
                                     float(settings.get("restart", 0.5)))
        return pred, strength
    if key == "physical_partners":
        layers, A = C._physical(ctx, target)
        if A is None:
            return None, None
        pred, support = S.graph_vote(A, visible, query)
        return pred, support
    if key == "layer_vote":
        pred, support, _w = C._weighted_vote(ctx, C._sources(ctx, target, int(settings.get("k", 15))),
                                             visible, query)
        return pred, support
    if key in ("map_neighbours", "cluster_guilt"):
        # Both read a map built BLIND to the target, so the map itself does not depend on which
        # genes are hidden and is built once for every fold; only the vote or the enrichment is
        # redone. `blind_map` caches on the context, so the second fold pays nothing.
        size = int(settings.get("sample") or 0) or None
        if key == "map_neighbours":
            coords, pos, _labels = ctx.blind_map(target, min(size or ctx.n, 2500),
                                                 int(settings.get("n_neighbors", 25)),
                                                 float(settings.get("min_dist", 0.1)))
        else:
            coords, pos, labels = ctx.blind_map(
                target, min(size or ctx.n, 2500),
                min_cluster_size=int(settings.get("min_cluster_size", 25)),
                selection=settings.get("selection", "leaf"))
        inside = np.intersect1d(np.asarray(query, dtype=int), pos)
        if not len(inside):
            return None, None
        at = {g: i for i, g in enumerate(pos)}
        sub = pd.Series(visible).reset_index(drop=True).iloc[pos].reset_index(drop=True)
        if key == "map_neighbours":
            pr, share = S.knn_vote(coords, sub, int(settings.get("k", 15)),
                                   query=[at[g] for g in inside])
            floor = float(settings.get("min_share", 0.0) or 0.0)
            return (S.expand_prediction(pr.where(share >= floor), pos, ctx.n),
                    S.expand_prediction(share, pos, ctx.n))
        pr, _table = C._enriched(labels, sub, float(settings.get("min_lift", 2.0)))
        return S.expand_prediction(pr, pos, ctx.n), None
    if key == "structural_homology":
        layer = C._usable_layer(ctx, settings.get("layer"), target)
        if layer is None:
            return None, None
        A = ctx.adjacency(layer, "w")
        truth_at_level = C._ec_level(ctx.truth(target), int(settings.get("level", 1)))
        visible = pd.Series(visible).reset_index(drop=True).where(truth_at_level.notna())
        return S.graph_vote(A, visible, query)
    if key in ("triangulation", "understudied_first"):
        pred, agree, _p = C._triangulate(C._tri_sources(ctx, target, int(settings.get("k", 15))),
                                         visible, query, int(settings.get("min_agree", 2)))
        return pred, agree
    return None, None


def supported(organism: str | None = None) -> list:
    """The strategies this can hold genes out for, in catalogue order.

    Those that call a label for named genes, and whose prediction can therefore be attributed to the
    gene it is about. A map walk or a module search answers about structure, not about one gene, and
    keeps its existing self-test instead.
    """
    keys = ("feature_knn", "map_neighbours", "cluster_guilt", "layer_propagation", "layer_vote",
            "physical_partners", "structural_homology", "supervised_classifier", "stratum_focus",
            "triangulation", "understudied_first", "graph_convolution", "random_forest", "stacking")
    known = {s.key for s in S.catalog()}
    return [k for k in keys if k in known]


# --------------------------------------------------------------------------- reading the ledger
def _rate(frame: pd.DataFrame) -> dict:
    """right/wrong/abstained counts, the rate over ANSWERED genes, and its Wilson interval."""
    from .calibration import wilson
    n = int(len(frame))
    abstained = int(frame["abstained"].sum()) if n else 0
    answered = n - abstained
    right = int(frame["correct"].fillna(False).sum()) if n else 0
    low, high = wilson(right, answered) if answered else (float("nan"), float("nan"))
    return {"genes": n, "answered": answered, "right": right, "wrong": answered - right,
            "abstained": abstained,
            "rate": right / answered if answered else float("nan"),
            "rate_low": low, "rate_high": high,
            "enough": answered >= MIN_FOR_RATE}


def summary(ledger: pd.DataFrame, level: str = "strategy") -> pd.DataFrame:
    """Pool the ledger at one level: "gene", "class", "target" or "strategy".

    Every row carries the counts it rests on, the rate over ANSWERED genes, a Wilson 95% interval and
    `enough`, which is False where too few genes were answered for a percentage to mean anything.
    Rows are never pooled across settings: a gene called right at one setting and wrong at another is
    not half right, it is a gene whose answer depends on the setting, and that is worth seeing.
    """
    if not len(ledger):
        return pd.DataFrame()
    by = {"gene": ["organism", "gene", "gene_id", "truth"],
          "class": ["organism", "target", "truth", "strategy", "setting_key"],
          "target": ["organism", "target", "strategy", "setting_key"],
          "strategy": ["organism", "strategy", "setting_key"]}[level]
    rows = []
    for keys, part in ledger.groupby(by, dropna=False, observed=True):
        row = dict(zip(by, keys if isinstance(keys, tuple) else (keys,)))
        row.update(_rate(part))
        if level == "gene":
            row["strategies"] = int(part["strategy"].nunique())
            wrong = part[(part["correct"] == False)]            # noqa: E712 -- a nullable column
            row["called_instead"] = ", ".join(sorted(set(wrong["prediction"].dropna()))[:3])
        if level == "class":
            wrong = part[(part["correct"] == False)]            # noqa: E712
            row["confused_with"] = ", ".join(
                wrong["prediction"].dropna().astype(str).value_counts().index[:3])
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values([c for c in ("organism", "target", "strategy", "truth", "gene_id")
                            if c in out.columns], kind="stable").reset_index(drop=True)


def for_gene(ledger: pd.DataFrame, gene_id: str) -> pd.DataFrame:
    """Every strategy's verdict on one gene: what it said, whether it was right, how sure it was."""
    part = ledger[ledger["gene_id"].astype(str) == str(gene_id)]
    keep = ["strategy", "setting_key", "truth", "prediction", "correct", "abstained", "support"]
    return part[keep].sort_values(["correct", "support"], ascending=[False, False],
                                  kind="stable").reset_index(drop=True)


def sentence(ledger: pd.DataFrame, gene_id: str) -> str:
    """One line for the gene card: the condensed form, with the detail a click away."""
    part = for_gene(ledger, gene_id)
    if not len(part):
        return ""
    answered = part[~part["abstained"].astype(bool)]
    right = int(answered["correct"].fillna(False).sum())
    truth = str(part["truth"].iloc[0])
    if not len(answered):
        return f"Known {truth}: no strategy would call it (all {len(part)} abstained)."
    return (f"Known {truth}: recovered by {right} of {len(answered)} strategies that answered"
            + (f", {len(part) - len(answered)} abstained." if len(part) > len(answered) else "."))


# --------------------------------------------------------------------------- sets held out together
def evaluate_sets(ctx, key: str, target: str, sets: dict, settings: dict | None = None,
                  seed: int = 0, log=None) -> pd.DataFrame:
    """Hide each named set of genes TOGETHER, and record what the strategy said about its members.

    A different question from the fold-by-fold run. There, a gene is hidden while its neighbours
    keep their labels, so the answer may rest on one close relative; here the whole set goes at once,
    which asks whether the evidence still reaches those genes when the obvious witnesses are gone
    too. Hiding a whole compartment is the strongest form: nothing of that class is left to copy
    from, so a label-caller cannot name it at all, and the honest result is that it abstains or
    places the genes somewhere else -- which is worth seeing stated rather than assumed.

    `sets` maps a name to gene positions. Returns the same rows as `evaluate`, with `mode` "set",
    the set's name and its size, so the two can be concatenated and read side by side.
    """
    say = log or (lambda _m: None)
    strategy = S.get(key)
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    chosen = _permitted(ctx, dict(strategy.settings(ctx, target=target, **(settings or {}))), target)
    setting_key = ", ".join(f"{k}={v}" for k, v in sorted(chosen.items()) if k != "target")
    rows = []
    for name, members in sets.items():
        ctx.check()
        hidden = np.asarray([g for g in np.asarray(members, dtype=int)
                             if isinstance(kept.iloc[g], str)], dtype=int)
        if not len(hidden):
            continue
        say(f"{key} on {target}: holding out {name} ({len(hidden):,} genes) together")
        visible = kept.copy()
        visible.iloc[hidden] = np.nan
        pred, support = _predict(ctx, strategy, visible, hidden, chosen, target)
        if pred is None:
            return pd.DataFrame(columns=list(COLUMNS))
        p = pd.Series(pred).reset_index(drop=True)
        sup = pd.Series(support).reset_index(drop=True) if support is not None else None
        for gene in hidden:
            call = p.iloc[gene]
            said = isinstance(call, str)
            rows.append({
                "organism": ctx.organism, "strategy": key, "target": target,
                "setting_key": setting_key, "seed": int(seed), "fold": -1, "mode": "set",
                "set_name": str(name), "set_size": int(len(hidden)),
                "gene": int(gene), "gene_id": str(ctx.gene_ids[gene]),
                "truth": str(kept.iloc[gene]), "prediction": call if said else None,
                "correct": bool(call == kept.iloc[gene]) if said else None,
                "abstained": not said,
                "support": float(sup.iloc[gene]) if sup is not None and said
                and np.isfinite(sup.iloc[gene]) else float("nan")})
    return pd.DataFrame(rows, columns=list(COLUMNS))


def class_sets(ctx, target: str) -> dict:
    """One set per label value: every gene of that class, hidden at once.

    The hardest hold-out there is. Nothing of the class remains to copy, so this measures whether a
    class is recoverable from the evidence itself or only by resemblance to its own members.
    """
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    return {str(value): np.flatnonzero((kept == value).to_numpy())
            for value in sorted(kept.dropna().unique(), key=str)}


def random_sets(ctx, target: str, sizes=(1, 2, 5, 20, 100), repeats: int = 5,
                seed: int = 0) -> dict:
    """Random sets of several sizes: how fast does a strategy lose a gene as its company is hidden?

    The degradation curve. One gene hidden alone is the easiest case there is; a hundred hidden with
    it is the case a real screen presents, where a whole list of genes is unknown at once.
    """
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    labelled = np.flatnonzero(kept.notna().to_numpy())
    rng = np.random.default_rng(seed)
    out = {}
    for size in sizes:
        if size > len(labelled):
            continue
        for repeat in range(int(repeats)):
            out[f"random {size} #{repeat + 1}"] = rng.choice(labelled, size=int(size),
                                                             replace=False)
    return out


def set_summary(ledger: pd.DataFrame) -> pd.DataFrame:
    """Per held-out set: how many members came back, and whether they stayed together.

    `together` is the share of the set given the SAME label as each other, whatever that label was:
    a set that is placed consistently but wrongly has still been recognised as one thing, and that
    is a different finding from a set scattered across the map.
    """
    part = ledger[ledger["mode"] == "set"]
    if not len(part):
        return pd.DataFrame()
    rows = []
    for keys, group in part.groupby(["organism", "strategy", "target", "setting_key", "set_name"],
                                    dropna=False, observed=True):
        row = dict(zip(("organism", "strategy", "target", "setting_key", "set_name"), keys))
        row["set_size"] = int(group["set_size"].iloc[0])
        row.update(_rate(group))
        # As strings: the shipped file stores labels as categories, and a categorical value_counts
        # lists every category, zero counts included -- it once named classes nothing was called.
        called = group["prediction"].dropna().astype(str)
        row["together"] = (float(called.value_counts().iloc[0] / len(called)) if len(called)
                           else float("nan"))
        row["placed_at"] = ", ".join(called.value_counts().index[:2]) if len(called) else ""
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["strategy", "set_name"],
                                          kind="stable").reset_index(drop=True)


# --------------------------------------------------------------------------- the shipped record
_SHIPPED: dict = {}


def shipped_path() -> str:
    """Where the built record lives in the package."""
    import os
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                        "track_record.parquet")


def shipped(organism: str | None = None) -> pd.DataFrame:
    """The built record, loaded once and kept. Empty when it has not been built here.

    Read rather than recomputed: the folds are a one-off build (`scripts/build_track_record.py`),
    so asking what happened to a gene costs a lookup and not a model fit.
    """
    import os
    if "all" not in _SHIPPED:
        path = shipped_path()
        try:
            _SHIPPED["all"] = pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame(
                columns=list(COLUMNS))
        except Exception:                                   # a half-written or foreign file
            _SHIPPED["all"] = pd.DataFrame(columns=list(COLUMNS))
    frame = _SHIPPED["all"]
    if organism is None or not len(frame):
        return frame
    # Kept per organism: a gene card asks several times per click, and filtering a million rows each
    # time was a tenth of a second apiece.
    key = ("organism", str(organism))
    if key not in _SHIPPED:
        _SHIPPED[key] = frame[frame["organism"] == str(organism)]
    return _SHIPPED[key]


#: Biological labels recorded beyond each space's declared targets. Columns that say where a label
#: came from (`compartment_source`, `ortholopit_donors`, `screen_scorers_agree`, `chromosome`, ...)
#: are left out on purpose: "recovering" which dataset assigned a compartment means nothing.
EXTRA_LABELS = {organisms.TOXOPLASMA: ("cellcycle_phase", "screen_actin_phenotype",
                                      "screen_apicoplast_phenotype",
                                      "screen_egress_phenotype", "screen_replication_phenotype")}


def labels(ctx) -> list:
    """The labels worth recording for this table, the default first."""
    from . import organisms
    have = set(ctx.categorical_columns())
    try:
        declared = list(organisms.get(ctx.organism).targets)
    except Exception:
        declared = []
    wanted = [S.default_category(ctx)] + declared + list(EXTRA_LABELS.get(ctx.organism, ()))
    # A `_derived` label is computed from measurements the strategies themselves read (the stage
    # label is the argmax of the expression columns): recovering it is circular, and the first build
    # duly "recovered" it at 91-100%. It is not recorded.
    return [t for t in dict.fromkeys(wanted) if t in have and not t.endswith("_derived")]

def recorded_targets(organism: str, ledger: pd.DataFrame | None = None) -> list:
    """The labels the record holds for one organism, the default first."""
    if ledger is None and ("targets", str(organism)) in _SHIPPED:
        return _SHIPPED[("targets", str(organism))]
    frame = shipped(organism) if ledger is None else ledger
    if not len(frame):
        return []
    have = [str(t) for t in pd.unique(frame["target"])]
    try:
        from . import organisms
        declared = [t for t in organisms.get(organism).targets if t in have]
    except Exception:                                       # an organism the registry lacks
        declared = []
    out = list(dict.fromkeys(declared[:1] + have))
    if ledger is None:
        _SHIPPED[("targets", str(organism))] = out
    return out


def default_target(organism: str, ledger: pd.DataFrame | None = None) -> str | None:
    """The label a view means when it names none: the organism's first declared target."""
    targets = recorded_targets(organism, ledger)
    return targets[0] if targets else None


def gene_html(gene_id: str, organism: str, ledger: pd.DataFrame | None = None,
              mode: str = "together", target: str | None = None) -> str:
    """The gene card's line and its table, as HTML: the condensed form first, the detail under it.

    One sentence a reader can take in -- what is known, and how many strategies recovered it when it
    was hidden -- then a row per strategy saying what each one actually said. Empty string when the
    record says nothing about this gene, so the card simply does not grow a section. `mode` "alone"
    reads the rows `alone` computed on demand, for a label the shipped record does not cover.
    """
    from html import escape
    frame = shipped(organism) if ledger is None else ledger
    folds = frame[frame["mode"] == mode] if len(frame) else frame
    if len(folds) and mode == "together":
        # One label at a time: the default unless asked, so the card never mixes two labels' rows.
        target = target or default_target(organism, frame)
        folds = folds[folds["target"].astype(str) == str(target)]
    one = for_gene(folds, gene_id) if len(folds) else pd.DataFrame()
    if not len(one):
        return ""
    line = sentence(folds, gene_id)
    # The class is a link: one click from this gene to how every strategy fares on its whole class.
    target = str(folds[folds["gene_id"].astype(str) == str(gene_id)]["target"].iloc[0])
    truth = str(one["truth"].iloc[0])
    link = (f"starplast://class/{quote(target, safe='')}/{quote(truth, safe='')}")
    # Only the shipped label has a class page to go to; an on-demand answer is about one gene.
    line_html = (escape(line).replace(escape(truth), f"<a href='{link}'>{escape(truth)}</a>", 1)
                 if mode == "together" else escape(line))
    heading = ("If this gene were unknown"
               if mode == "together" and target == (default_target(organism) or target)
               else f"If its {escape(target.replace('_', ' '))} were unknown")
    rows = ""
    for r in one.itertuples():
        said = "—" if r.abstained else escape(str(r.prediction))
        mark = "·" if r.abstained else ("✓" if r.correct else "✗")
        colour = "#888888" if r.abstained else ("#2e8b57" if r.correct else "#c0392b")
        support = "" if not np.isfinite(r.support) else f"{r.support:.2f}"
        rows += (f"<tr><td>{escape(S.get(r.strategy).title)}</td>"
                 f"<td style='color:{colour}'>{mark} {said}</td>"
                 f"<td align='right' style='color:#888'>{support}</td></tr>")
    return (f"<h4>{heading}</h4><p>{line_html}</p>"
            f"<p style='color:#888'>Each strategy below was asked about this gene with its label "
            f"hidden, and its orthogroup hidden with it, so nothing could answer by copying a "
            f"paralog. ✓ right · ✗ wrong · · declined to answer; the last column is how sure it "
            f"was.</p><table cellspacing='0' cellpadding='3'>{rows}</table>")


def class_html(target: str, label: str, organism: str,
               ledger: pd.DataFrame | None = None) -> str:
    """One class: how well each strategy recovers it, and what it is mistaken for.

    The level above the gene. A class every strategy misses is not a failure of any one of them --
    it is a statement that this distinction is not in the measurements, which is worth knowing
    before a screen is designed around it.
    """
    from html import escape
    frame = shipped(organism) if ledger is None else ledger
    folds = frame[(frame["mode"] == "together") & (frame["target"].astype(str) == str(target))]
    here = folds[folds["truth"].astype(str) == str(label)]
    if not len(here):
        return ""
    rows = ""
    for r in summary(here, "class").sort_values("rate", ascending=False,
                                                na_position="last").itertuples():
        rate = ("too few" if not r.enough else
                f"{r.rate:.0%} <span style='color:#888'>[{r.rate_low:.0%}, {r.rate_high:.0%}]</span>")
        rows += (f"<tr><td>{escape(S.get(r.strategy).title)}</td>"
                 f"<td align='right'>{r.right}/{r.answered}</td><td>{rate}</td>"
                 f"<td style='color:#888'>{escape(str(r.confused_with or ''))}</td></tr>")
    # Distinct genes, not rows: every strategy contributes a row per gene, so the row count is the
    # gene count times the number of strategies -- 10,766 for a class of 769, which it once said.
    genes = int(here["gene"].nunique())
    best = summary(here, "class").sort_values("rate", ascending=False, na_position="last")
    lead = (f"{escape(str(label))}: {genes:,} genes, asked of {here['strategy'].nunique()} "
            f"strategies. "
            + (f"Best recovered by {escape(S.get(best.iloc[0]['strategy']).title)} "
               f"({best.iloc[0]['right']} of {best.iloc[0]['answered']})."
               if len(best) and best.iloc[0]["enough"] else
               "No strategy answers enough of them to judge."))
    up = (f"<p><a href='starplast://target/{quote(str(target), safe='')}'>"
          f"all {escape(str(target))} classes ▸</a></p>")
    return (f"<h4>{escape(str(label))}</h4><p>{lead}</p>{up}"
            f"<table cellspacing='0' cellpadding='3'><tr><th align='left'>strategy</th>"
            f"<th>right</th><th align='left'>rate [95%]</th>"
            f"<th align='left'>called instead</th></tr>{rows}</table>")


def target_html(target: str, organism: str, ledger: pd.DataFrame | None = None) -> str:
    """Every class of one information category on a line: how recoverable each one is.

    The level above the class. Each class names its best strategy and that strategy's rate, and
    links down to the class page; classes no strategy answers often enough to judge are said to be
    so rather than ranked, and classes no strategy recovers are listed as such -- they are the
    distinctions the measurements do not carry.
    """
    from html import escape
    frame = shipped(organism) if ledger is None else ledger
    folds = frame[(frame["mode"] == "together") & (frame["target"].astype(str) == str(target))]
    if not len(folds):
        return ""
    by_class = summary(folds, "class")
    by_class = by_class[by_class["answered"] > 0]
    genes = folds.groupby(folds["truth"].astype(str), observed=True)["gene"].nunique()
    lines = []
    for label, part in by_class.groupby(by_class["truth"].astype(str)):
        judged = part[part["enough"].astype(bool)].sort_values("rate", ascending=False)
        best = judged.iloc[0] if len(judged) else None
        lines.append((-1.0 if best is None else float(best["rate"]), label, best))
    rows = ""
    for rate, label, best in sorted(lines, key=lambda t: (-t[0], t[1])):
        link = (f"<a href='starplast://class/{quote(str(target), safe='')}/"
                f"{quote(label, safe='')}'>{escape(label)}</a>")
        if best is None:
            verdict = "<span style='color:#888'>too few answered to judge</span>"
        elif best["right"] == 0:
            verdict = "<span style='color:#888'>never recovered</span>"
        else:
            # The short name: a table of 24 classes cannot carry 24 full titles. The class page,
            # one click down, names each strategy in full.
            verdict = (f"{best['rate']:.0%} by {escape(str(best['strategy']).replace('_', ' '))} "
                       f"<span style='color:#888'>({best['right']}/{best['answered']})</span>")
        rows += (f"<tr><td>{link}</td><td align='right'>{int(genes.get(label, 0)):,}</td>"
                 f"<td>{verdict}</td></tr>")
    lead = (f"{len(lines)} classes, {int(folds['gene'].nunique()):,} genes held out. "
            f"Best strategy per class; click a class for all of them.")
    # The same category from the strategies' side: two measures, because a strategy built for small
    # classes loses plain accuracy while winning per class, and either alone would misjudge it.
    strat_rows = ""
    scored = []
    for key, part in folds.groupby(folds["strategy"].astype(str)):
        b = _baselines(part)
        if b is not None:
            scored.append((b["balanced"], key, b))
    # The guess is judged on every held-out gene once, not on any one strategy's answered subset.
    every = folds.drop_duplicates("gene")["truth"].astype(str)
    ref = {"commonest": float(every.value_counts(normalize=True).iloc[0]),
           "chance": 1.0 / max(every.nunique(), 1)}
    for _bal, key, b in sorted(scored, reverse=True):
        strat_rows += (f"<tr><td>{escape(key.replace('_', ' '))}</td>"
                       f"<td align='right'>{b['accuracy']:.0%}</td>"
                       f"<td align='right'>{b['balanced']:.0%}</td></tr>")
    strategies = ""
    if strat_rows:
        strategies = (f"<p style='color:#888'>By strategy, over the genes each answered: right overall, "
                      f"and right per class averaged so small classes count equally. Guessing the "
                      f"commonest class scores about {ref['commonest']:.0%} and {ref['chance']:.0%}."
                      f"</p><table cellspacing='0' cellpadding='3'><tr><th align='left'>strategy</th>"
                      f"<th>overall</th><th>per class</th></tr>{strat_rows}</table>")
    return (f"<h4>{escape(str(target))}</h4><p>{lead}</p>"
            f"<table cellspacing='0' cellpadding='3'><tr><th align='left'>class</th>"
            f"<th>genes</th><th align='left'>best recovered</th></tr>{rows}</table>{strategies}")


def weakest(strategy: str, organism: str, target: str | None = None, limit: int = 3,
            ledger: pd.DataFrame | None = None) -> str:
    """One line for a strategy card: where this strategy is weakest, by class.

    Only classes with enough answered genes to judge; a class it answered twice says nothing about
    it, and saying so would be worse than saying nothing.
    """
    frame = shipped(organism) if ledger is None else ledger
    folds = frame[(frame["mode"] == "together") & (frame["strategy"].astype(str) == str(strategy))]
    target = target or default_target(organism, frame)
    folds = folds[folds["target"].astype(str) == str(target)]
    if not len(folds):
        return ""
    rows = summary(folds, "class")
    rows = rows[rows["enough"].astype(bool)].sort_values("rate", kind="stable")
    pooled = _rate(folds)
    if not len(rows):
        return (f"Right on {pooled['right']:,} of {pooled['answered']:,} held-out genes "
                f"({pooled['rate']:.0%}).")
    worst = ", ".join(f"{r.truth} ({r.right} of {r.answered})" for r in rows.head(limit).itertuples())
    return (f"Right on {pooled['right']:,} of {pooled['answered']:,} held-out genes "
            f"({pooled['rate']:.0%}); weakest on {worst}.")


def my_list(ctx, genes, target: str | None = None, strategies=None, log=None) -> pd.DataFrame:
    """Hide YOUR genes together and ask every strategy what they are: would it have found them?

    `genes` is any list of identifiers -- a screen's hits, a complex, a pull-down -- resolved the way
    the gene-list strategies resolve them. Each strategy is asked about all of them at once, with
    their labels hidden, so the answer says whether the evidence reaches your genes or only reaches
    genes whose neighbours are already known. Returns the per-gene rows; `set_summary` pools them.

        rows = track_record.my_list(ctx, ["TGME49_294550", "TGME49_244470"], "compartment")
        track_record.set_summary(rows)[["strategy", "right", "answered", "together"]]
    """
    found, _missing = ctx.resolve_genes(genes) if not isinstance(genes, np.ndarray) else (genes, [])
    positions = np.asarray(found, dtype=int)
    target = target or S.default_category(ctx)
    parts = []
    for key in (strategies or supported()):
        rows = evaluate_sets(ctx, key, target, {"your list": positions}, log=log)
        if len(rows):
            parts.append(rows)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=list(COLUMNS))


def _baselines(folds: pd.DataFrame) -> dict | None:
    """Accuracy and mean per-class recall, beside what always naming the commonest class scores.

    Two measures because strategies trade one for the other: a strategy built to find small classes
    (label diffusion gives every class the same seed mass) loses plain accuracy to the commonest-class
    guess while recovering small classes far above chance. Judged on accuracy alone it would look
    worse than not looking at the data, which it is not.
    """
    answered = folds[~folds["abstained"].astype(bool)]
    if len(answered) < MIN_FOR_RATE:
        return None
    truth = answered["truth"].astype(str)
    correct = answered["correct"].fillna(False).astype(bool)
    return {"accuracy": float(correct.mean()),
            "commonest": float(truth.value_counts(normalize=True).iloc[0]),
            "balanced": float(correct.groupby(truth).mean().mean()),
            "chance": 1.0 / truth.nunique()}


def _record_folds(strategy: str, organism: str, target: str | None) -> pd.DataFrame:
    try:
        frame = shipped(organism)
    except Exception:                                     # no record built, or a foreign one
        return pd.DataFrame()
    if not len(frame):
        return frame
    folds = frame[(frame["mode"] == "together") & (frame["strategy"].astype(str) == str(strategy))]
    target = target or default_target(organism, frame)
    return folds[folds["target"].astype(str) == str(target)]


def beats_baseline(strategy: str, organism: str, target: str | None = None) -> bool | None:
    """False only when the strategy loses to always naming the commonest class on BOTH accuracy and
    mean per-class recall; None where there is no record to judge from."""
    b = _baselines(_record_folds(strategy, organism, target))
    if b is None:
        return None
    return bool(b["accuracy"] > b["commonest"] or b["balanced"] > b["chance"])


def record_phrase(strategy: str, organism: str, target: str | None = None) -> str:
    """A clause for a recommendation: how often this strategy was right with labels hidden.

    Only from the shipped record and only for its label: if `target` is given and is not the label
    the record holds, nothing is said rather than quoting a different label's rate as this one's.
    """
    folds = _record_folds(strategy, organism, target)
    b = _baselines(folds) if len(folds) else None
    if b is None:
        return ""
    r = _rate(folds)
    label = str(folds["target"].iloc[0]).replace("_", " ")
    head = (f"with {label} hidden, right on {r['right']:,} of the {r['answered']:,} genes it "
            f"answered ({b['accuracy']:.0%}")
    if b["accuracy"] > b["commonest"]:
        return f"{head}; {b['commonest']:.0%} by always naming the commonest class)"
    if b["balanced"] > b["chance"]:
        return (f"{head}, below the {b['commonest']:.0%} of always naming the commonest class, "
                f"because it favours small classes: {b['balanced']:.0%} averaged over classes, "
                f"against {b['chance']:.0%} by chance)")
    return f"{head}, no better than always naming the commonest class)"

def list_html(rows: pd.DataFrame, rows_of: dict | None = None) -> str:
    """What `my_list` found, condensed first: a line per strategy, then a gene-by-strategy grid.

    `rows_of` maps a gene id to its row in the window's table, so each gene links to its card.
    """
    from html import escape
    if not len(rows):
        return ("<h4>Would they have found your genes?</h4><p style='color:#888'>None of your genes "
                "carries a label these strategies can be tested on.</p>")
    rows = rows.assign(strategy=rows["strategy"].astype(str), gene_id=rows["gene_id"].astype(str))
    target = str(rows["target"].iloc[0]).replace("_", " ")
    genes = list(dict.fromkeys(rows["gene_id"]))
    keys = list(dict.fromkeys(rows["strategy"]))
    summary_rows = ""
    for key in keys:
        part = rows[rows["strategy"] == key]
        r = _rate(part)
        called = part["prediction"].dropna().astype(str)
        placed = ", ".join(called.value_counts().index[:2]) if len(called) else "—"
        summary_rows += (f"<tr><td>{escape(key.replace('_', ' '))}</td>"
                         f"<td align='right'>{r['right']}/{r['answered']}</td>"
                         f"<td style='color:#888'>{escape(placed)}</td></tr>")
    best = max(keys, key=lambda k: _rate(rows[rows["strategy"] == k])["right"])
    b = _rate(rows[rows["strategy"] == best])
    lead = (f"Your {len(genes)} labelled gene{'s' if len(genes) != 1 else ''}, hidden together "
            f"and asked their {escape(target)}. "
            f"Best: {escape(best.replace('_', ' '))}, right on {b['right']} of {b['answered']}.")
    head = "".join(f"<th title='{escape(k)}'>{i + 1}</th>" for i, k in enumerate(keys))
    grid = ""
    marks = {(g, k): ("·", "#888888") for g in genes for k in keys}
    for r in rows.itertuples():
        marks[(r.gene_id, r.strategy)] = (("·", "#888888") if r.abstained else
                                          ("✓", "#2e8b57") if r.correct else ("✗", "#c0392b"))
    for g in genes:
        truth = str(rows[rows["gene_id"] == g]["truth"].iloc[0])
        name = (f"<a href='starplast://gene/{int(rows_of[g])}'>{escape(g)}</a>"
                if rows_of and g in rows_of else escape(g))
        cells = "".join(f"<td align='center' style='color:{marks[(g, k)][1]}'>{marks[(g, k)][0]}</td>"
                        for k in keys)
        grid += f"<tr><td>{name}</td><td style='color:#888'>{escape(truth)}</td>{cells}</tr>"
    legend = " · ".join(f"{i + 1} {escape(k.replace('_', ' '))}" for i, k in enumerate(keys))
    return (f"<h4>Would they have found your genes?</h4><p>{lead}</p>"
            f"<table cellspacing='0' cellpadding='3'><tr><th align='left'>strategy</th>"
            f"<th>right</th><th align='left'>placed them at</th></tr>{summary_rows}</table>"
            f"<p style='color:#888'>Gene by gene (✓ right · ✗ wrong · · declined): {legend}</p>"
            f"<table cellspacing='0' cellpadding='2'><tr><th align='left'>gene</th>"
            f"<th align='left'>known</th>{head}</tr>{grid}</table>")


#: Fast enough to answer for one gene while someone waits: a pass of each takes seconds on the full
#: Toxoplasma table. The slow learners (graph convolution, stacking) are in the shipped record.
ALONE = ("feature_knn", "layer_vote", "physical_partners", "structural_homology", "random_forest")


def alone(ctx, gene, target: str | None = None, strategies=ALONE, cache: bool = True,
          log=None) -> pd.DataFrame:
    """Hide ONE gene, with its orthogroup, and ask each fast strategy what it is.

    For any label, not only the shipped one: the shipped record covers each organism's default
    label, and this answers the same question for the others on demand. The orthogroup goes too,
    exactly as in the folds, so a paralogue cannot give its sibling away. `gene` is a position or an
    identifier. Results are cached per release, organism, label and gene.

        T.alone(ctx, "TGME49_294550", "function")
    """
    from . import __version__, paths
    position = int(gene) if isinstance(gene, (int, np.integer)) else ctx.index.get(str(gene).upper())
    if position is None:
        return pd.DataFrame(columns=list(COLUMNS))
    target = target or S.default_category(ctx)
    gene_id = str(ctx.gene_ids[position])
    path = os.path.join(paths.user_cache_dir(), "track_record",
                        f"alone_{__version__}_{ctx.organism}_{_safe(target)}_{_safe(gene_id)}.parquet")
    if cache and os.path.exists(path):
        try:
            done = pd.read_parquet(path)
            if set(strategies) <= set(done["strategy"].astype(str)):
                return done[done["strategy"].astype(str).isin(strategies)].reset_index(drop=True)
        except Exception:                                 # a damaged cache is rebuilt, not trusted
            pass
    groups = np.asarray(ctx.groups())
    siblings = np.flatnonzero(groups == groups[position])
    parts = []
    for key in strategies:
        rows = evaluate_sets(ctx, key, target, {gene_id: siblings}, log=log)
        rows = rows[rows["gene"] == position]
        if len(rows):
            parts.append(rows.assign(mode="alone"))
    out = (pd.concat(parts, ignore_index=True) if parts
           else pd.DataFrame(columns=list(COLUMNS)))
    if cache and len(out):
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            out.to_parquet(path, index=False)
        except OSError:                                   # a read-only cache costs time, not answers
            pass
    return out


def _safe(text: str) -> str:
    """A string as a file-name component."""
    return "".join(c if c.isalnum() or c in "-." else "_" for c in str(text))[:80]
