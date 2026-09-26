"""Strategies 35-39: learners that say how sure they are, and learners that combine evidence.

    35  Call genes with a stated error rate                (split conformal)
    36  Smooth the measurements along the networks          (graph convolution + logistic)
    37  Let a random forest find what defines a label        (random forest + permutation importance)
    38  Learn how much to trust each kind of evidence        (stacked logistic regression)
    39  Predict a value with an interval that holds          (gradient boosting + split conformal)

35 and 39 add what no other strategy offers: an error rate stated BEFORE the answer is seen.
Conformal prediction turns any model's scores into sets (or intervals) that contain the truth for at
least 1 - alpha of new genes, whatever the model, as long as the calibration genes resemble the new
ones -- so a call is made only where the set is a single label, and the self-test checks the stated
coverage on genes the calibration never saw. 36 lets a linear model see a gene's network
neighbourhood; 37 lets a non-linear one find thresholds and interactions; 38 learns, out of fold,
how far to trust each kind of evidence for this label.

Every calibration and every meta-model split is by whole orthogroups, as every test here is.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from . import scorecard as SC
from . import strategies as S
from .strategies import MIN_CLASS, Param, Strategy, StrategyResult, register
from .strategy_catalog import (K, MODEL, NUMBER, OWN_KIND, TARGET, _calls_table, _genes_table,
                               _logistic, _need, _regression_matrix, _unlabelled_but_known)

#: The ninth family: models that state their own error rate, see the networks, find thresholds, or
#: learn how to combine the others. Numbered 35-39, after the families they build on.
ADVANCED = "Advanced models"

ALPHA = Param("alpha", "float", "Error rate allowed (alpha)",
              "The share of new genes whose true answer may fall outside what is returned. 0.1 "
              "promises that at least 90% of prediction sets contain the true label (or 90% of "
              "intervals the true value); smaller alpha gives bigger sets and wider intervals.",
              0.1, lo=0.01, hi=0.5, step=0.01)


def _group_split(ctx, positions, frac: float, salt: int) -> tuple:
    """(keep, held): `positions` split so whole orthogroups land on one side, about `frac` held."""
    positions = np.asarray(positions, dtype=int)
    g = ctx.groups()[positions]
    rng = ctx.rng(900 + salt)
    order = rng.permutation(np.unique(g))
    sizes = pd.Series(g).value_counts()
    target, total, chosen = frac * len(positions), 0, set()
    for grp in order:
        if total >= target:
            break
        chosen.add(grp)
        total += int(sizes[grp])
    held = np.isin(g, list(chosen))
    return positions[~held], positions[held]


def _conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """The split-conformal threshold: the ceil((n + 1)(1 - alpha))-th smallest score."""
    scores = np.sort(np.asarray(scores, dtype=float))
    n = len(scores)
    if not n:
        return float("inf")
    k = int(math.ceil((n + 1) * (1 - alpha)))
    return float("inf") if k > n else float(scores[k - 1])


# =========================================================================== 35 · conformal calls
def _probabilities(ctx, X, vis: pd.Series, query, p) -> pd.DataFrame:
    """Class scores for `query` from a model trained on `vis` only (genes x classes, NaN elsewhere)."""
    if p["model"] == "kNN":
        pred, _s = S.knn_vote(X, vis, int(p["k"]), query=query)
    else:
        pred, _p, _m = _logistic(X, vis, query, float(p["C"]))
    scores = S.class_scores_of(pred)
    return scores if scores is not None else pd.DataFrame(index=range(len(vis)))


def _conformal(ctx, X, vis: pd.Series, query, p) -> tuple:
    """(prediction, set size, set text, thresholds): split conformal label sets for `query`.

    The visible labelled genes are split by orthogroup into a training part and a calibration
    part. The model learns on the first; on the second, each gene's nonconformity is one minus the
    score its TRUE label received. The threshold is that score's conformal quantile -- per class
    when `per_class`, so rare classes get their own guarantee -- and a gene's set is every label
    scoring within it. A call is made only where the set holds exactly one label.
    """
    vis = pd.Series(vis).reset_index(drop=True)
    query = np.asarray(query, dtype=int)
    labelled = np.flatnonzero(vis.notna().to_numpy())
    train, cal = _group_split(ctx, labelled, 0.25, 35)
    fit_vis = vis.copy()
    fit_vis.iloc[cal] = np.nan
    P = _probabilities(ctx, X, fit_vis, np.union1d(cal, query), p)
    classes = list(P.columns)
    alpha, per_class = float(p["alpha"]), p["thresholds"] == "per class"
    truth_cal = vis.iloc[cal].to_numpy(dtype=object)
    score_of = {c: j for j, c in enumerate(classes)}
    Pm = P.to_numpy(dtype=float)
    nonconf = np.array([1.0 - (Pm[g, score_of[t]] if t in score_of and np.isfinite(Pm[g, score_of[t]])
                               else 0.0) for g, t in zip(cal, truth_cal)])
    overall = _conformal_quantile(nonconf, alpha)
    q = {c: overall for c in classes}
    if per_class:
        # A class needs at least ceil(1 / alpha) - 1 calibration genes before its own quantile is
        # finite; below that it would enter every set and no gene could ever be called. Such rare
        # classes use the overall threshold -- their guarantee is the marginal one, and the
        # thresholds table says which classes those are.
        need = int(math.ceil(1 / alpha)) - 1
        for c in classes:
            here = nonconf[truth_cal == c]
            if len(here) >= need:
                q[c] = _conformal_quantile(here, alpha)
    n = len(vis)
    pred = pd.Series([np.nan] * n, dtype=object)
    size = pd.Series(np.nan, index=range(n))
    sets = pd.Series([""] * n, dtype=object)
    for g in query:
        row = Pm[g]
        members = [c for j, c in enumerate(classes)
                   if np.isfinite(row[j]) and 1.0 - row[j] <= q[c]]
        size.iloc[g] = len(members)
        sets.iloc[g] = " | ".join(members)
        if len(members) == 1:
            pred.iloc[g] = members[0]
    scores = P.copy()
    rest = np.setdiff1d(np.arange(n), query)
    scores.iloc[rest] = np.nan
    return S.with_class_scores(pred, scores), size, sets, q


def _conformal_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    X, _cols = ctx.features(target)
    query = np.flatnonzero(t.isna().to_numpy())
    pred, size, sets, q = _conformal(ctx, X, t, query, p)
    calls = _calls_table(ctx, pred, 1.0 / size.where(size > 0), truth,
                         prediction_set=sets.to_numpy(), set_size=size.to_numpy())
    reach = size.iloc[query]
    labelled = t.notna().to_numpy()
    thresholds = pd.DataFrame({"class": list(q), "threshold": list(q.values()),
                               "labelled_genes": [int((t[labelled] == c).sum()) for c in q]})
    summary = (f"Split conformal prediction at alpha {float(p['alpha']):.2f}: for at least "
               f"{1 - float(p['alpha']):.0%} of new genes the prediction set contains the true "
               f"{target} label. Of {len(query):,} unlabelled genes, {int((reach == 1).sum()):,} get a "
               f"single label (called), {int((reach > 1).sum()):,} a set of several (the data cannot "
               f"decide between them) and {int((reach == 0).sum()):,} an empty set (unlike any "
               f"calibration gene).")
    sets_table = _genes_table(ctx, query, prediction_set=sets.iloc[query].to_numpy(),
                              set_size=size.iloc[query].to_numpy())
    return StrategyResult("conformal_calls", summary,
                          {"calls": calls, "prediction sets": sets_table.sort_values(
                              "set_size", kind="stable").reset_index(drop=True),
                           "thresholds": thresholds}, genes=list(calls["gene_id"][:200]))


def _efficiency(size: pd.Series, classes: int) -> float:
    """1 - (mean set size - 1) / (classes - 1): 1 when every set is one label, 0 when every set holds
    them all. An empty set counts as one label here; the coverage check counts it as a miss."""
    if classes < 2 or not len(size):
        return float("nan")
    s = np.clip(size.to_numpy(dtype=float), 1, classes)
    return float(1.0 - (s.mean() - 1.0) / (classes - 1.0))


def _conformal_test(ctx, p):
    """Validity and efficiency: do the sets keep their promise, and how much do they narrow it?

    A conformal set that holds every label always keeps its promise and says nothing, so coverage
    alone cannot be the test. The verdict is efficiency -- how far the sets narrow the possibilities
    -- against the same procedure on shuffled labels, where the model learns nothing and the sets
    must grow to keep the promise. Coverage is reported beside it and should sit at or above
    1 - alpha; the scorecard scores the single-label calls.
    """
    import time
    t0 = time.monotonic()
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    visible, hidden = S.hide(truth, 0.25, ctx.seed, groups=ctx.groups())
    classes = int(visible.dropna().nunique())
    metric = "set efficiency: 1 - (mean set size - 1) / (classes - 1)"
    if len(hidden) < S.MIN_HIDDEN or classes < 2:
        return S.judge("conformal_calls", metric, float("nan"), [], min_effect=0.05,
                       n_hidden=len(hidden), hidden=f"{target} labels", null_kind="shuffled labels",
                       t0=t0, note="too few labelled genes to hide", task=SC.T_LABEL)
    pred, size, sets, _q = _conformal(ctx, X, visible, hidden, p)
    observed = _efficiency(size.iloc[hidden], classes)
    t = truth.iloc[hidden].to_numpy(dtype=object)
    inside = [str(a) in str(x).split(" | ") for a, x in zip(t, sets.iloc[hidden])]
    rng = ctx.rng(351)
    nulls, null_cover = [], []
    for i in range(10):
        ctx.say(f"conformal_calls: null {i + 1} of 10")
        _pn, sz, st, _qn = _conformal(ctx, X, S.shuffled(visible, rng), hidden, p)
        nulls.append(_efficiency(sz.iloc[hidden], classes))
        null_cover.append(float(np.mean([str(a) in str(x).split(" | ")
                                         for a, x in zip(t, st.iloc[hidden])])))
    card = SC.label_calls(pred, truth, hidden, S.class_scores_of(pred))
    return S.judge("conformal_calls", metric, observed, nulls, min_effect=0.05,
                   n_hidden=len(hidden),
                   hidden=f"25% of the {target} labels, whole orthogroups at a time",
                   null_kind="10 runs on shuffled labels, calibrated the same way", t0=t0,
                   task=SC.T_LABEL, scorecard=card,
                   numbers={"set_coverage": float(np.mean(inside)),
                            "promised_coverage": 1 - float(p["alpha"]),
                            "null_set_coverage": float(np.mean(null_cover)),
                            "mean_set_size": float(size.iloc[hidden].mean()),
                            "classes": classes,
                            "singleton_share": float((size.iloc[hidden] == 1).mean())})


register(Strategy(
    key="conformal_calls", number=35, family=ADVANCED,
    title="Call genes with a stated error rate",
    method="split conformal prediction",
    task="label calls", techniques=("split_conformal", "logistic_regression", "knn"),
    question="Which genes can be given a label with a guaranteed error rate -- and for which does "
             "the data leave two or more labels equally possible?",
    tooltip="Turns a classifier's scores into prediction sets that contain the true label for at "
            "least 1 - alpha of new genes, calibrated on orthogroups the model never saw; calls a "
            "gene only where its set is a single label, and lists the genes the data cannot decide.",
    explanation=(
        "Every other label strategy returns a call and, at best, a score. A score is not an error "
        "rate: a probability of 0.8 from a logistic regression is right far more or far less than "
        "80% of the time depending on the class and the data. Conformal prediction fixes that "
        "without assumptions about the model. The labelled genes are split by orthogroup; the "
        "model learns on one part, and on the other the strategy records how low a score the TRUE "
        "label can get. The threshold that covers all but alpha of those genes then defines, for "
        "every new gene, the set of labels that score within it -- and at least 1 - alpha of such "
        "sets contain the truth, for any model, as long as new genes resemble the calibration "
        "genes.\n\n"
        "A set of one label is a call made with that guarantee behind it. A set of two labels is a "
        "finding in itself -- the data cannot tell those compartments apart for this gene -- and an "
        "empty set says the gene is unlike every calibration gene. With per-class thresholds "
        "(Mondrian conformal) each class gets its own guarantee, so a rare compartment is not "
        "covered by borrowing the common ones' accuracy. The self-test checks both the calls and "
        "the promise: the precision of single-label calls against shuffled labels, and the share "
        "of hidden genes whose set contains their true label, against the 1 - alpha it promised."),
    walkthrough=(
        "Choose the label; keep alpha at 0.1 (90% of sets hold the truth).",
        "Press Test: read 'set_coverage' in numbers -- it should be close to or above the promised "
        "coverage -- and 'mean_set_size' out of 'classes': how far the data narrows the answer.",
        "Press Run: 'calls' are single-label sets; 'prediction sets' lists every unlabelled gene "
        "with its set, smallest first.",
        "Lower alpha for a stronger promise (bigger sets, fewer calls); switch the model to kNN "
        "for a faster, model-free base score."),
    test_description=(
        "25% of the label hidden by whole orthogroups; the model is trained and calibrated on the "
        "rest (itself split by orthogroup) and builds a set for every hidden gene. Metric: set "
        "efficiency, 1 - (mean set size - 1) / (classes - 1) -- how far the sets narrow the "
        "possibilities. Null: 10 runs on shuffled labels, where the model learns nothing and the "
        "sets must grow to keep the promise. Pass: above the null's 95th percentile by 0.05. Also "
        "reported: set coverage on the hidden genes beside the 1 - alpha promised, and the "
        "scorecard of the single-label calls."),
    params=(TARGET, ALPHA,
            Param("model", "choice", "Base model",
                  "Where the scores the sets are built from come from: a class-balanced logistic "
                  "regression on every permitted measurement, or the distance-weighted vote of the "
                  "nearest labelled genes. The guarantee holds for either.", "logistic",
                  choices=("logistic", "kNN")),
            Param("thresholds", "choice", "Guarantee",
                  "'per class' calibrates a threshold for each label separately (Mondrian "
                  "conformal), so the promise holds within every class; 'overall' uses one "
                  "threshold, which can over-cover common classes and under-cover rare ones.",
                  "per class", choices=("per class", "overall")),
            Param("C", "float", "Regularisation C",
                  "Inverse penalty strength of the logistic base model: small values force a "
                  "simple model, large values let it fit detail. Unused by the kNN base model.",
                  0.5, lo=0.001, hi=100.0, step=0.1),
            K),
    runner=_conformal_run, tester=_conformal_test, cost="minutes",
    needs=("a categorical column",)))


# =========================================================================== 36 · graph convolution
def _smoothed(ctx, target, hops: int) -> tuple:
    """(features, columns, layers): the measurements, then their average over 1..hops network steps.

    Simplified graph convolution: every permitted measured layer is merged into one graph with
    self-loops, symmetrically normalised, and applied to the measurement matrix `hops` times. Each
    power is kept as its own block, so the model can weigh the gene against its neighbourhood.
    """
    import scipy.sparse as sp
    key = ("sgc", target, int(hops))
    if key in ctx._cache:
        return ctx._cache[key]
    X, cols = ctx.features(target)
    layers = ctx.measurement_layers(target)
    blocks = [X]
    if layers and hops > 0:
        A = sum(ctx.adjacency(l, "binary") for l in layers)
        A = (A > 0).astype(float) + sp.eye(ctx.n, format="csr")
        d = np.asarray(A.sum(axis=1)).ravel()
        D = sp.diags(1.0 / np.sqrt(d))
        Ahat = D @ A @ D
        H = X
        for _ in range(int(hops)):
            H = np.asarray(Ahat @ H)
            blocks.append(H)
    out = (np.hstack(blocks), cols, layers)
    ctx._cache[key] = out
    return out


def _sgc_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    H, cols, layers = _smoothed(ctx, target, int(p["hops"]))
    pred, prob, model = _logistic(H, t, np.flatnonzero(t.isna().to_numpy()), float(p["C"]),
                                  float(p["min_probability"]))
    calls = _calls_table(ctx, pred, prob, truth)
    rows = []
    if model is not None:
        w = np.abs(model.coef_)
        per = len(cols)
        total = w.sum()
        for h in range(H.shape[1] // per):
            share = float(w[:, h * per:(h + 1) * per].sum() / total) if total else float("nan")
            rows.append({"block": "the gene itself" if h == 0 else f"{h}-step neighbourhood",
                         "share of model weight": share})
    summary = (f"{len(cols)} permitted measurements smoothed over {len(layers)} measured layer(s) "
               f"for {int(p['hops'])} step(s), and a class-balanced logistic regression trained on "
               f"{int(t.notna().sum()):,} labelled genes calls {len(calls):,} unlabelled ones. "
               f"'where the model looks' says how much weight it puts on a gene's own measurements "
               f"against its neighbourhood's.")
    return StrategyResult("graph_convolution", summary,
                          {"calls": calls, "where the model looks": pd.DataFrame(rows)},
                          genes=list(calls["gene_id"][:200]))


def _sgc_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    H, _cols, _layers = _smoothed(ctx, target, int(p["hops"]))
    predict = lambda vis: _logistic(H, vis, _unlabelled_but_known(vis, truth), float(p["C"]),
                                    float(p["min_probability"]))[0]
    return S.label_transfer_test(ctx, "graph_convolution", target, predict, null="analytic")


register(Strategy(
    key="graph_convolution", number=36, family=ADVANCED,
    title="Smooth the measurements along the networks, then classify",
    method="graph convolution + logistic regression",
    task="label calls", techniques=("sgc", "logistic_regression"),
    question="Does a gene's label follow from its own measurements together with those of its "
             "network neighbours -- and how much does the model lean on each?",
    tooltip="Averages every permitted measurement over each gene's partners in the measured "
            "networks, one and two steps out, and trains a class-balanced logistic regression on "
            "the gene's own profile beside its neighbourhood's -- a graph neural network without "
            "the black box.",
    explanation=(
        "Strategy 19 sees a gene's own measurements; strategies 11-13 see its partners' labels. "
        "Neither sees its partners' MEASUREMENTS, which is what a graph neural network learns "
        "from. This strategy does it in the simplest form that works well in practice, a "
        "simplified graph convolution: the permitted measured layers are merged into one graph "
        "with self-loops, normalised by degree, and multiplied into the measurement matrix once "
        "and twice. The result is three blocks per gene -- itself, its neighbourhood, and its "
        "neighbourhood's neighbourhood -- and a logistic regression learns from all three.\n\n"
        "Only measurements travel along the edges, never labels, so a hidden gene's label cannot "
        "leak through its neighbours; the held-out label's closure is removed from both the "
        "columns and the layers first. The 'where the model looks' table is a finding: a "
        "compartment whose weight sits on the neighbourhood blocks is one the networks encode "
        "better than the gene's own profile does. For genes with no edges the neighbourhood "
        "blocks equal the gene itself, so the strategy degrades to strategy 19 rather than "
        "failing."),
    walkthrough=(
        "Choose the label; keep two steps and C at 0.5.",
        "Press Test and compare its scorecard with strategy 19's: the difference is what the "
        "networks add.",
        "Press Run; read 'where the model looks' before the calls.",
        "Try one step if the networks are dense -- two steps on a dense graph average most of the "
        "proteome together."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; the smoothed features are "
        "computed from measurements only, and the model is trained on the visible genes. Metric: "
        "hidden genes called correctly. Null: the analytic chance level for the same predicted and "
        "true class mixes. Pass: above that level's 95% bound by 0.05."),
    params=(TARGET,
            Param("hops", "int", "Network steps",
                  "How many steps of neighbourhood averaging are added as feature blocks. Zero is "
                  "strategy 19; one adds direct partners; two adds partners of partners, which "
                  "helps on sparse networks and blurs everything on dense ones.", 2, lo=0, hi=4),
            Param("C", "float", "Regularisation C",
                  "Inverse penalty strength of the logistic regression. Smoothing multiplies the "
                  "number of features, so a stronger penalty (smaller C) than strategy 19's is "
                  "often right.", 0.5, lo=0.001, hi=100.0, step=0.1),
            Param("min_probability", "float", "Call when probability is at least",
                  "A gene is called only when the model's top probability reaches this. Zero calls "
                  "every gene; 0.6 keeps the calls the model is sure of.", 0.0, lo=0.0, hi=1.0,
                  step=0.05)),
    runner=_sgc_run, tester=_sgc_test, cost="a minute",
    needs=("a categorical column", "measured edge layers")))


# =========================================================================== 37 · random forest
def _forest(X, vis: pd.Series, query, trees: int, leaf: int, seed: int) -> tuple:
    """(prediction, probability, model): a class-balanced random forest on the visible genes."""
    from sklearn.ensemble import RandomForestClassifier
    vis = pd.Series(vis).reset_index(drop=True)
    query = np.asarray(query, dtype=int)
    pred = pd.Series([np.nan] * len(vis), dtype=object)
    prob = pd.Series(np.nan, index=range(len(vis)))
    known = np.flatnonzero(vis.notna().to_numpy())
    y = vis.iloc[known].astype(str).to_numpy()
    if len(set(y)) < 2 or not len(query) or X.shape[1] == 0:
        return pred, prob, None
    model = RandomForestClassifier(n_estimators=int(trees), min_samples_leaf=int(leaf),
                                   max_features="sqrt", class_weight="balanced_subsample",
                                   n_jobs=4, random_state=int(seed)).fit(X[known], y)
    P = model.predict_proba(X[query])
    pred.iloc[query] = model.classes_[P.argmax(axis=1)]
    prob.iloc[query] = P.max(axis=1)
    pred = S.with_class_scores(pred, S._score_frame(len(vis), query, P, model.classes_))
    return pred, prob, model


def _permutation_importance(model, X, y, columns, repeats: int, rng) -> tuple:
    """(mean, sd) loss of balanced accuracy when each of `columns` is shuffled among these genes."""
    from sklearn.metrics import balanced_accuracy_score
    base = balanced_accuracy_score(y, model.predict(X))
    mean, sd = [], []
    for j in columns:
        drops = []
        for _ in range(int(repeats)):
            Xp = X.copy()
            Xp[:, j] = rng.permutation(Xp[:, j])
            drops.append(base - balanced_accuracy_score(y, model.predict(Xp)))
        mean.append(float(np.mean(drops)))
        sd.append(float(np.std(drops)))
    return np.array(mean), np.array(sd)


def _forest_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    X, cols = ctx.features(target)
    pred, prob, model = _forest(X, t, np.flatnonzero(t.isna().to_numpy()), p["trees"],
                                p["min_leaf"], ctx.seed)
    calls = _calls_table(ctx, pred, prob, truth)
    importance = pd.DataFrame()
    if model is not None:
        # Importance is measured on orthogroups the forest did not train on: impurity importance
        # alone rewards measurements with many distinct values, whatever they predict.
        ctx.say("measuring permutation importance on held-out orthogroups")
        labelled = np.flatnonzero(t.notna().to_numpy())
        _train, held = _group_split(ctx, labelled, 0.25, 37)
        fit = t.copy()
        fit.iloc[held] = np.nan
        _p, _pr, m2 = _forest(X, fit, held, p["trees"], p["min_leaf"], ctx.seed)
        if m2 is not None and len(held):
            top = np.argsort(-m2.feature_importances_)[:30]
            mean, sd = _permutation_importance(m2, X[held], t.iloc[held].astype(str).to_numpy(),
                                               top, 3, ctx.rng(370))
            importance = pd.DataFrame({"measurement": [cols[j] for j in top], "importance": mean,
                                       "sd": sd,
                                       "impurity_importance": m2.feature_importances_[top]}
                                      ).sort_values("importance", ascending=False,
                                                    kind="stable").reset_index(drop=True)
    summary = (f"A random forest of {int(p['trees'])} trees on {len(cols)} permitted measurements, "
               f"trained on {int(t.notna().sum()):,} labelled genes, calls {len(calls):,} unlabelled "
               f"ones. 'what defines the label' ranks measurements by how much balanced accuracy is "
               f"lost on held-out orthogroups when each is shuffled.")
    return StrategyResult("random_forest", summary,
                          {"calls": calls, "what defines the label": importance},
                          genes=list(calls["gene_id"][:200]))


def _forest_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    predict = lambda vis: _forest(X, vis, _unlabelled_but_known(vis, truth), p["trees"],
                                  p["min_leaf"], ctx.seed)[0]
    return S.label_transfer_test(ctx, "random_forest", target, predict, null="analytic")


register(Strategy(
    key="random_forest", number=37, family=ADVANCED,
    title="Let a random forest find what defines a label",
    method="random forest + permutation importance",
    task="label calls", techniques=("random_forest", "permutation_importance"),
    question="Which measurements, in which combinations and past which thresholds, define a label "
             "-- and which unlabelled genes carry that definition?",
    tooltip="Trains a class-balanced random forest on every permitted measurement, calls the "
            "unlabelled genes, and ranks the measurements by how much accuracy is lost on held-out "
            "orthogroups when each is shuffled -- the non-linear counterpart of strategy 19.",
    explanation=(
        "A logistic regression draws straight boundaries: more of this measurement, more of that "
        "class. Biology is often not like that. A secreted protein might be recognised by high "
        "tachyzoite expression AND a signal peptide, but not by either alone; an essential gene "
        "by fitness below a threshold whatever its other values. A random forest finds such "
        "interactions and thresholds by growing hundreds of decision trees on random subsets of "
        "genes and measurements and letting them vote, with each class weighted so the commonest "
        "compartment does not win by default.\n\n"
        "What it gains in flexibility it loses in transparency, so the strategy reports "
        "permutation importance rather than the forest's own impurity importance: each of the "
        "forest's top thirty measurements is shuffled among genes of orthogroups the forest never "
        "trained on, and the loss of balanced accuracy is recorded. A measurement the forest merely "
        "split on often scores near zero here, which is the point. Compare its scorecard with "
        "strategy 19's: if the forest is not better, the label's signal is linear and the simpler "
        "model's weights are the better explanation."),
    walkthrough=(
        "Choose the label; keep 300 trees and a minimum leaf of 2.",
        "Press Test and compare the scorecard with strategies 19 and 36.",
        "Press Run; read 'what defines the label' -- measurements whose importance is within two "
        "standard deviations of zero do not matter, whatever their impurity importance.",
        "Raise the minimum leaf size if the labels are noisy; more trees only make the answer "
        "steadier, never more flexible."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; the forest is trained on the "
        "visible genes and calls the hidden ones. Metric: hidden genes called correctly. Null: the "
        "analytic chance level for the same predicted and true class mixes. Pass: above that "
        "level's 95% bound by 0.05."),
    params=(TARGET,
            Param("trees", "int", "Trees",
                  "How many decision trees vote. More trees give a steadier answer and a smoother "
                  "probability, never a more flexible model; 300 is past the point of change for a "
                  "few hundred measurements.", 300, lo=20, hi=2000, step=20),
            Param("min_leaf", "int", "Genes per leaf at least",
                  "The smallest group of genes a tree may end on. Larger leaves smooth the model "
                  "and resist noisy labels; 1 lets every tree memorise its training genes.", 2,
                  lo=1, hi=50)),
    runner=_forest_run, tester=_forest_test, cost="a minute",
    needs=("a categorical column",)))


# =========================================================================== 38 · stacking
def _stack_bases(ctx, target, k: int) -> dict:
    """The evidence stacked: measurement neighbours, a linear model, and the measured networks."""
    X, _cols = ctx.features(target)
    out = {"measurements (kNN)": lambda vis, q: S.knn_vote(X, vis, k, q)[0],
           "measurements (logistic)": lambda vis, q: _logistic(X, vis, q)[0]}
    layers = ctx.measurement_layers(target)
    if layers:
        A = sum(ctx.adjacency(l, "binary") for l in layers)
        out["networks (vote)"] = lambda vis, q: S.graph_vote(A, vis, q)[0]
    return out


def _meta_features(preds: dict, rows, classes) -> np.ndarray:
    """Per base: its score for every class (0 where it abstained) and whether it spoke at all."""
    blocks = []
    for pred in preds.values():
        sc = S.class_scores_of(pred)
        M = np.zeros((len(rows), len(classes)))
        if sc is not None:
            M = sc.reindex(columns=classes).iloc[rows].to_numpy(dtype=float)
        spoke = np.isfinite(M).any(axis=1).astype(float)[:, None]
        blocks += [np.nan_to_num(M, nan=0.0), spoke]
    return np.hstack(blocks)


def _stack(ctx, bases: dict, vis: pd.Series, query, folds: int = 5) -> tuple:
    """(prediction, support, trust): a logistic meta-model over out-of-fold base predictions."""
    from sklearn.linear_model import LogisticRegression
    vis = pd.Series(vis).reset_index(drop=True)
    query = np.asarray(query, dtype=int)
    classes = sorted(vis.dropna().unique(), key=str)
    labelled = np.flatnonzero(vis.notna().to_numpy())
    n = len(vis)
    pred = pd.Series([np.nan] * n, dtype=object)
    support = pd.Series(np.nan, index=range(n))
    if len(classes) < 2 or len(labelled) < 5 * folds or not len(query):
        return pred, support, pd.DataFrame()
    g = ctx.groups()[labelled]
    rng = ctx.rng(380)
    fold_of = {grp: i % folds for i, grp in enumerate(rng.permutation(np.unique(g)))}
    fold = np.array([fold_of[x] for x in g])
    Z = np.zeros((len(labelled), (len(classes) + 1) * len(bases)))
    for f in range(folds):
        ctx.check()
        rows = labelled[fold == f]
        train = vis.copy()
        train.iloc[rows] = np.nan
        preds = {name: fn(train, rows) for name, fn in bases.items()}
        Z[fold == f] = _meta_features(preds, rows, classes)
    meta = LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced").fit(
        Z, vis.iloc[labelled].astype(str).to_numpy())
    preds = {name: fn(vis, query) for name, fn in bases.items()}
    P = meta.predict_proba(_meta_features(preds, query, classes))
    pred.iloc[query] = meta.classes_[P.argmax(axis=1)]
    support.iloc[query] = P.max(axis=1)
    pred = S.with_class_scores(pred, S._score_frame(n, query, P, meta.classes_))
    w = np.abs(meta.coef_).sum(axis=0)
    width = len(classes) + 1
    trust = pd.DataFrame({"evidence": list(bases),
                          "share of meta-model weight": [float(w[i * width:(i + 1) * width].sum()
                                                               / max(w.sum(), 1e-12))
                                                         for i in range(len(bases))]}
                         ).sort_values("share of meta-model weight", ascending=False,
                                       kind="stable").reset_index(drop=True)
    return pred, support, trust


def _stack_run(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    t = truth.where(truth.map(truth.value_counts()) >= MIN_CLASS)
    bases = _stack_bases(ctx, target, int(p["k"]))
    pred, support, trust = _stack(ctx, bases, t, np.flatnonzero(t.isna().to_numpy()),
                                  int(p["folds"]))
    calls = _calls_table(ctx, pred, support, truth)
    summary = (f"{len(bases)} kinds of evidence, each run out of fold on the labelled genes "
               f"({int(p['folds'])} folds by orthogroup), and a logistic meta-model learns how far "
               f"to trust each for {target}; it calls {len(calls):,} unlabelled genes. 'trust by "
               f"evidence' is what it learned.")
    return StrategyResult("stacking", summary, {"calls": calls, "trust by evidence": trust},
                          genes=list(calls["gene_id"][:200]))


def _stack_test(ctx, p):
    target = _need(ctx, p["target"])
    truth = ctx.truth(target)
    bases = _stack_bases(ctx, target, int(p["k"]))
    predict = lambda vis: _stack(ctx, bases, vis, _unlabelled_but_known(vis, truth),
                                 int(p["folds"]))[0]
    return S.label_transfer_test(ctx, "stacking", target, predict, null="analytic")


register(Strategy(
    key="stacking", number=38, family=ADVANCED,
    title="Learn how much to trust each kind of evidence",
    method="stacked logistic regression",
    task="label calls", techniques=("stacking", "knn", "logistic_regression", "network_vote",
                                    "grouped_cv"),
    question="Given measurement neighbours, a linear model and the measured networks, how should "
             "their answers be combined for THIS label -- and what does the combination call?",
    tooltip="Runs three kinds of evidence out of fold on the labelled genes and trains a logistic "
            "meta-model on their class scores, so the weight each earns is learned per label and "
            "per class rather than fixed; reports that trust beside the calls.",
    explanation=(
        "Strategy 31 calls a gene when enough independent strategies agree; strategy 12 weights "
        "each source by one accuracy figure. Stacking learns the combination instead. Each kind "
        "of evidence -- the nearest genes in the measurements, a logistic regression on them, and "
        "the vote of network partners -- predicts every labelled gene out of fold, from orthogroups "
        "it was not trained on, and a second logistic regression learns from those predictions "
        "which evidence to believe for which class. The networks might be decisive for complexes "
        "and useless for secreted proteins; the meta-model can learn exactly that, and whether a "
        "network's silence (no labelled partner) is itself informative.\n\n"
        "Because the meta-model only ever sees out-of-fold predictions, it cannot learn to trust a "
        "base model for having memorised its training genes -- the classic failure of naive "
        "blending. 'trust by evidence' reports each source's share of the meta-model's weight, a "
        "per-label answer to which experiments matter. Stacking is usually the most accurate "
        "strategy here and the slowest; its scorecard beside strategies 07, 13 and 19 says how much "
        "combining is worth for your label."),
    walkthrough=(
        "Choose the label; keep 5 folds and 15 neighbours.",
        "Press Test (several model fits; a few minutes) and compare the scorecard with 07, 19 and "
        "31.",
        "Press Run; read 'trust by evidence' before the calls.",
        "Fewer folds are faster and give the meta-model noisier training data; more are slower "
        "and steadier."),
    test_description=(
        "Pattern 1: 25% of the label hidden by whole orthogroups; base predictions for the "
        "meta-model are out of fold within the visible genes only, so no hidden label reaches "
        "either level. Metric: hidden genes called correctly. Null: the analytic chance level for "
        "the same predicted and true class mixes. Pass: above that level's 95% bound by 0.05."),
    params=(TARGET, K,
            Param("folds", "int", "Out-of-fold splits",
                  "How many orthogroup folds the base predictions are made in. Each fold's genes "
                  "are predicted by base models that never saw them, which is what keeps the "
                  "meta-model honest.", 5, lo=3, hi=10)),
    runner=_stack_run, tester=_stack_test, cost="minutes",
    needs=("a categorical column",)))


# =========================================================================== 39 · conformal values
def _intervals(ctx, X, y_vis: pd.Series, kind: str, alpha: float) -> tuple:
    """(prediction, lower, upper, half-width): split conformal intervals around a regression.

    The measured genes are split by orthogroup; the model learns on one part, and the conformal
    quantile of the absolute errors on the other is the half-width that covers at least 1 - alpha
    of new genes.
    """
    y_vis = pd.Series(y_vis).reset_index(drop=True).astype(float)
    known = np.flatnonzero(y_vis.notna().to_numpy())
    train, cal = _group_split(ctx, known, 0.25, 39)
    fit = y_vis.copy()
    fit.iloc[cal] = np.nan
    from .strategy_catalog import _fit_predict
    pred = _fit_predict(X, fit, kind)
    half = _conformal_quantile(np.abs(pred[cal] - y_vis.iloc[cal].to_numpy()), alpha)
    return pred, pred - half, pred + half, half


def _interval_run(ctx, p):
    target = _need(ctx, p["target"])
    y = ctx.values(target)
    X, cols = _regression_matrix(ctx, target, p["model"], p["own_kind"])
    pred, lo, hi, half = _intervals(ctx, X, y, p["model"], float(p["alpha"]))
    missing = np.flatnonzero(y.isna().to_numpy())
    filled = _genes_table(ctx, missing, predicted=pred[missing], lower=lo[missing],
                          upper=hi[missing]).sort_values("predicted", kind="stable")
    measured = np.flatnonzero(y.notna().to_numpy())
    yv = y.iloc[measured].to_numpy()
    outside = (yv < lo[measured]) | (yv > hi[measured])
    surprises = _genes_table(ctx, measured[outside], measured=yv[outside],
                             predicted=pred[measured][outside], lower=lo[measured][outside],
                             upper=hi[measured][outside])
    surprises["distance_outside"] = np.maximum(surprises["lower"] - surprises["measured"],
                                               surprises["measured"] - surprises["upper"])
    surprises = surprises.sort_values("distance_outside", ascending=False,
                                      kind="stable").reset_index(drop=True)
    summary = (f"{target} predicted from {len(cols)} permitted measurements ({p['model']}), with "
               f"intervals of +/- {half:.3g} that contain the measured value for at least "
               f"{1 - float(p['alpha']):.0%} of new genes. {len(missing):,} genes without a "
               f"measurement receive an interval; {int(outside.sum()):,} measured genes fall outside "
               f"theirs (expected by chance: about {float(p['alpha']) * len(measured):,.0f}), listed "
               f"as surprises.")
    return StrategyResult("conformal_values", summary,
                          {"predicted intervals": filled.reset_index(drop=True),
                           "surprises": surprises},
                          genes=list(surprises["gene_id"][:200]))


def _interval_test(ctx, p):
    target = _need(ctx, p["target"])
    y = ctx.values(target)
    X, _cols = _regression_matrix(ctx, target, p["model"], p["own_kind"])
    seen = {}

    def predict(vis):
        pred, lo, hi, half = _intervals(ctx, X, vis, p["model"], float(p["alpha"]))
        if not seen:
            hidden = np.flatnonzero(vis.isna().to_numpy() & y.notna().to_numpy())
            yh = y.iloc[hidden].to_numpy()
            seen.update({"interval_coverage": float(np.mean((yh >= lo[hidden]) & (yh <= hi[hidden])))
                         if len(hidden) else float("nan"),
                         "promised_coverage": 1 - float(p["alpha"]),
                         "interval_width_in_sd": float(2 * half / np.nanstd(yh))
                         if len(hidden) and np.nanstd(yh) > 0 else float("nan")})
        return pred

    result = S.value_test(ctx, "conformal_values", y, predict, n_null=3, min_effect=0.1,
                          label=f"the {target} values")
    result.numbers.update(seen)
    return result


register(Strategy(
    key="conformal_values", number=39, family=ADVANCED,
    title="Predict a value with an interval that holds",
    method="gradient boosting / ridge + split conformal",
    task="values", techniques=("gradient_boosting", "ridge", "split_conformal", "grouped_cv"),
    question="For a gene never measured, what value is expected -- and within what range, with a "
             "guaranteed chance of containing the truth?",
    tooltip="Predicts a measurement from every other permitted one and wraps each prediction in an "
            "interval calibrated on orthogroups the model never saw, so at least 1 - alpha of new "
            "genes fall inside; measured genes outside their interval are listed as surprises.",
    explanation=(
        "Strategy 21 predicts a value; this one says how far to trust it, in the measurement's own "
        "units. The measured genes are split by orthogroup: a model learns on one part, and on "
        "the other the strategy records how large its errors are. The conformal quantile of those "
        "errors becomes the half-width of every interval, and at least 1 - alpha of new genes' "
        "true values fall inside -- a guarantee that holds for any model, provided new genes "
        "resemble the calibration genes.\n\n"
        "Two readings follow. For unmeasured genes, an interval narrower than the measurement's "
        "spread is a prediction worth acting on, and one as wide as the spread says the other "
        "measurements know little about this one. For measured genes, a value outside its interval "
        "is a surprise the rest of the data cannot explain -- by construction about alpha of genes "
        "will be, so the list is ranked by how far outside, and the count is shown beside what "
        "chance predicts. The self-test scores the predictions against hidden values and checks "
        "the coverage promised on them."),
    walkthrough=(
        "Choose a measurement (a fitness score is the classic case); keep alpha at 0.1.",
        "Press Test: read 'interval_coverage' against 'promised_coverage', and the interval width "
        "in standard deviations.",
        "Press Run; 'predicted intervals' fills unmeasured genes, 'surprises' lists measured ones "
        "outside their range.",
        "Switch to ridge for a faster, linear model; leave its own kind of measurement out so the "
        "answer is not one screen predicting another."),
    test_description=(
        "Pattern 4: 20% of the measured values hidden; the model is trained and calibrated on the "
        "rest (split by orthogroup) and predicts the hidden ones. Metric: rank correlation of "
        "predicted and hidden values. Null: the chance distribution of a rank correlation, with "
        "refits on shuffled values reported. Pass: above its 95th percentile by 0.1. Also "
        "reported: the share of hidden values inside their intervals, beside the coverage "
        "promised."),
    params=(NUMBER, ALPHA, MODEL, OWN_KIND),
    runner=_interval_run, tester=_interval_test, cost="a minute",
    needs=("a numeric column",)))
