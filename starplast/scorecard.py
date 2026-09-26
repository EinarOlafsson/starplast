"""Scorecards: the same metrics, in the same order, for every strategy that does the same kind of task.

A self-test's verdict rests on one number judged against its null. That number is chosen to be the
honest one for the strategy, but on its own it cannot say *how* a strategy is good or bad: whether it
is precise but reaches few genes, whether it finds the common classes and misses the rare ones,
whether a ranking is good at the top or only on average. The scorecard adds that, and standardises
it: every strategy states which of six tasks it performs, and every strategy performing a task
reports that task's metrics in the same order, computed the same way on the same hidden genes.

The six tasks, and what is hidden in each:

* **label calls** -- a share of a label is hidden and called back (localization, stage, ...);
* **ranking** -- hidden positives (set members, edges, corrupted labels) are ranked among negatives;
* **set retrieval** -- a set of genes is returned and compared with the hidden members;
* **cluster recovery** -- clusters chosen on visible labels are scored on hidden ones;
* **values** -- hidden measurements (fitness, abundance) are predicted;
* **replication** -- findings made on half the genes are checked on the other half.

Every metric is defined once, in :data:`METRICS`, with its range, its chance level and how to read
it; :func:`glossary` returns that as a table, and the application shows it wherever a metric
appears. Metrics that need something a strategy does not produce -- per-class scores for a macro
AUROC, say -- are reported as missing rather than guessed, and the glossary says why.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

NAN = float("nan")


@dataclass(frozen=True)
class Metric:
    """One metric: what it measures, its range, what chance gives, and how to read it."""
    key: str
    label: str
    task: str
    definition: str
    range: str
    chance: str
    reading: str
    higher_is_better: bool = True


@dataclass(frozen=True)
class Task:
    """A kind of inference, with the metrics every strategy performing it reports, in order."""
    key: str
    description: str
    hidden: str
    metrics: tuple


_M = []


def _metric(key, label, task, definition, range_, chance, reading, higher=True):
    _M.append(Metric(key, label, task, definition, range_, chance, reading, higher))


# --------------------------------------------------------------------------- label calls
T_LABEL = "label calls"
_metric("accuracy", "Accuracy", T_LABEL,
        "Hidden genes called with their true label, divided by ALL hidden genes. A gene the "
        "strategy declined to call counts as wrong.",
        "0 to 1", "about the sum of squared class shares (Cohen's chance term)",
        "The headline for label calls. Counting abstentions as errors stops a strategy looking "
        "accurate by calling only the easy genes; coverage and precision of calls split it apart.")
_metric("coverage", "Coverage", T_LABEL,
        "Share of hidden genes that received any call.",
        "0 to 1", "not applicable (a property of the strategy, not of luck)",
        "How far the strategy reaches. A network strategy cannot call a gene with no edges; a "
        "low coverage with a high precision of calls is a precise but narrow tool.")
_metric("precision_of_calls", "Precision of calls", T_LABEL,
        "Correct calls divided by calls made (also called selective accuracy).",
        "0 to 1", "as accuracy, among the genes called",
        "How far to trust one call. Equals accuracy when coverage is 1.")
_metric("macro_precision", "Macro precision", T_LABEL,
        "For each true class: of the genes called that class, the share that truly are; averaged "
        "over classes with equal weight. A class never called scores 0.",
        "0 to 1", "about the average class share",
        "Whether calls of the RARE classes can be trusted too; a strategy that only ever calls the "
        "commonest class scores low.")
_metric("macro_recall", "Macro recall (balanced accuracy)", T_LABEL,
        "For each true class: the share of its hidden genes called correctly; averaged over classes "
        "with equal weight. Identical to balanced accuracy.",
        "0 to 1", "1 / number of classes, for a caller that ignores the data",
        "Whether every class is found, not just the large ones. Compare with accuracy: a large gap "
        "means the strategy lives on the big classes.")
_metric("macro_f1", "Macro F1", T_LABEL,
        "Per class, the harmonic mean of precision and recall; averaged over classes with equal "
        "weight.",
        "0 to 1", "low; roughly the average class share",
        "One number balancing finding each class and being right when calling it, with rare classes "
        "counted as much as common ones.")
_metric("weighted_f1", "Weighted F1", T_LABEL,
        "Per-class F1 averaged with each class weighted by its number of hidden genes.",
        "0 to 1", "about the sum of squared class shares",
        "Macro F1's counterpart that follows the class sizes; close to accuracy when coverage is "
        "high.")
_metric("kappa", "Cohen's kappa", T_LABEL,
        "Agreement between calls and truth corrected for the agreement their class frequencies "
        "alone would give: (observed - expected) / (1 - expected). No call is its own category.",
        "-1 to 1", "0",
        "Accuracy with chance removed: 0 is no better than matching class frequencies, 1 is "
        "perfect. Comparable across labels with different numbers and sizes of classes.")
_metric("mcc", "Matthews correlation (MCC)", T_LABEL,
        "The multiclass Matthews correlation coefficient (Gorodkin's R_K) between calls and truth, "
        "with no call as its own category.",
        "-1 to 1", "0",
        "A correlation between the call and the truth that stays honest under strong class "
        "imbalance; often the single most informative number for an unbalanced label.")
_metric("macro_auroc", "Macro AUROC", T_LABEL,
        "For each class, the probability that a random hidden gene of that class gets a higher "
        "score for the class than a random hidden gene of another class; averaged over classes. "
        "Needs per-class scores.",
        "0 to 1", "0.5",
        "How well the strategy's scores separate each class from the rest, before any threshold "
        "is chosen. Missing for strategies that call labels without scoring every class.")
_metric("macro_auprc", "Macro AUPRC", T_LABEL,
        "For each class, the area under the precision-recall curve of its one-vs-rest scores "
        "(average precision); averaged over classes. Needs per-class scores.",
        "0 to 1", "the class's share among hidden genes, averaged",
        "Like macro AUROC but dominated by the top of each ranking, which is where calls are made; "
        "far below AUROC means good separation overall but a noisy top.")

# --------------------------------------------------------------------------- ranking
T_RANK = "ranking"
_metric("auroc", "AUROC", T_RANK,
        "Probability that a random hidden positive is ranked above a random negative (area under "
        "the ROC curve); ties count half.",
        "0 to 1", "0.5",
        "Separation over the whole ranking. Insensitive to how rare positives are, so a high AUROC "
        "can coexist with a poor top of the list -- read it with AUPRC.")
_metric("auprc", "AUPRC", T_RANK,
        "Average precision: the mean, over the hidden positives, of the precision of the ranking "
        "down to that positive (area under the precision-recall curve).",
        "0 to 1", "the prevalence (share of positives among everything ranked)",
        "How clean the top of the ranking is. Its chance level is the prevalence, so compare it "
        "with that (AUPRC lift), never with 0.5.")
_metric("auprc_lift", "AUPRC lift", T_RANK,
        "AUPRC divided by the prevalence, its value for a random ordering.",
        "0 to 1/prevalence", "1",
        "How many times better than a random ordering the ranking is, where it matters. Comparable "
        "across tests with different prevalences, unlike AUPRC itself.")
_metric("prevalence", "Prevalence", T_RANK,
        "Hidden positives divided by everything ranked.",
        "0 to 1", "not applicable (a property of the test)",
        "The chance level of AUPRC and of precision at any depth; context for every other number.")
_metric("pauroc_10", "Partial AUROC (FPR <= 10%)", T_RANK,
        "AUROC restricted to false-positive rates up to 10%, rescaled (McClish) so 0.5 is chance "
        "and 1 perfect.",
        "0.5 to 1 (after rescaling)", "0.5",
        "Separation in the part of the ranking anyone would act on; what matters for screening.")
_metric("r_precision", "R-precision", T_RANK,
        "Precision among the top R items, where R is the number of hidden positives. Equal to "
        "recall at that depth.",
        "0 to 1", "the prevalence",
        "If you took as many candidates as there are true positives, the share that would be "
        "right.")
_metric("precision_at_1pct", "Precision @ top 1%", T_RANK,
        "Share of positives among the top 1% of the ranking (at least one item).",
        "0 to 1", "the prevalence",
        "How good the very first candidates are -- the list anyone would test first.")
_metric("enrichment_at_1pct", "Enrichment @ top 1%", T_RANK,
        "Precision at the top 1% divided by the prevalence.",
        "0 to 1/prevalence", "1",
        "How many times more positives the first 1% holds than a random 1% would.")
_metric("recall_at_10pct", "Recall @ top 10%", T_RANK,
        "Share of all hidden positives found in the top 10% of the ranking.",
        "0 to 1", "0.1",
        "How much of the answer a short list contains.")
_metric("max_f1", "Best F1", T_RANK,
        "The highest F1 (harmonic mean of precision and recall) over every cut-off of the ranking.",
        "0 to 1", "about 2 x prevalence / (1 + prevalence) (taking everything)",
        "The best single-threshold trade-off the ranking allows; optimistic, since the cut-off is "
        "chosen on the answer.")
_metric("ndcg", "nDCG", T_RANK,
        "Normalised discounted cumulative gain: positives count 1 / log2(rank + 1), divided by the "
        "value of a perfect ranking.",
        "0 to 1", "depends on prevalence; low when positives are rare",
        "A whole-ranking score that weights the top most, smoothly rather than at one cut-off.")

# --------------------------------------------------------------------------- set retrieval
T_SET = "set retrieval"
_metric("precision", "Precision", T_SET,
        "Share of the returned genes that are hidden members.",
        "0 to 1", "the members' share of the candidate genes",
        "How many of the candidates are real.")
_metric("recall", "Recall", T_SET,
        "Share of the hidden members that were returned.",
        "0 to 1", "the returned share of the candidate genes",
        "How much of the set was found.")
_metric("f1", "F1", T_SET,
        "Harmonic mean of precision and recall.",
        "0 to 1", "low; about the members' share for a random return of the same size",
        "One number that is only high when both are.")
_metric("jaccard", "Jaccard index", T_SET,
        "Overlap of returned and hidden sets divided by their union.",
        "0 to 1", "near 0",
        "How close the returned set is to being exactly the hidden set; stricter than F1.")
_metric("set_mcc", "Matthews correlation (MCC)", T_SET,
        "Correlation between 'returned' and 'is a hidden member' over all candidate genes.",
        "-1 to 1", "0",
        "Accounts for the genes correctly NOT returned as well; honest when the set is small.")
_metric("fold_enrichment", "Fold enrichment", T_SET,
        "Precision divided by the members' share of the candidate genes.",
        "0 to 1/share", "1",
        "How many times more members the returned set holds than a random set of its size.")
_metric("returned", "Genes returned", T_SET,
        "How many genes the strategy returned as the set.",
        "0 to all", "not applicable",
        "Context: a high recall from returning half the genome is not a finding.")

# --------------------------------------------------------------------------- cluster recovery
T_CLUSTER = "cluster recovery"
_metric("weighted_f1_clusters", "Weighted F1", T_CLUSTER,
        "For each label, F1 of its hidden genes against the ONE cluster chosen for it on visible "
        "genes; averaged with labels weighted by their hidden genes. A label with no cluster scores "
        "0.",
        "0 to 1", "set by permuting hidden labels over the same clusters (the verdict's null)",
        "Whether a label that nobody showed the clustering falls out as a cluster, scored on genes "
        "the choice never saw.")
_metric("weighted_precision_clusters", "Weighted precision", T_CLUSTER,
        "For each label, the share of hidden genes in its chosen cluster that carry it; "
        "label-size weighted.",
        "0 to 1", "about the label's share",
        "How pure the chosen clusters are.")
_metric("weighted_recall_clusters", "Weighted recall", T_CLUSTER,
        "For each label, the share of its hidden genes that landed in its chosen cluster; "
        "label-size weighted.",
        "0 to 1", "about the cluster's share of the genes",
        "How completely each label is gathered into one cluster.")
_metric("ari", "Adjusted Rand index (ARI)", T_CLUSTER,
        "Agreement between the clustering and the hidden labels over all pairs of hidden genes, "
        "corrected so random partitions of the same sizes score 0.",
        "-0.5 to 1", "0",
        "Whole-partition agreement, independent of which cluster was chosen for which label.")
_metric("nmi", "Normalised mutual information (NMI)", T_CLUSTER,
        "Mutual information between clusters and hidden labels, divided by the mean of their "
        "entropies.",
        "0 to 1", "small but above 0 for many small clusters",
        "How much knowing a gene's cluster tells you about its label.")
_metric("homogeneity", "Homogeneity", T_CLUSTER,
        "1 minus the uncertainty about the label left once the cluster is known, relative to the "
        "label's own entropy.",
        "0 to 1", "near 0",
        "Whether each cluster holds one label.")
_metric("completeness", "Completeness", T_CLUSTER,
        "1 minus the uncertainty about the cluster left once the label is known, relative to the "
        "clustering's entropy.",
        "0 to 1", "near 0",
        "Whether each label sits in one cluster. Its harmonic mean with homogeneity is the "
        "V-measure, which equals NMI.")
_metric("noise_share", "Unclustered share", T_CLUSTER,
        "Share of hidden genes HDBSCAN left as noise (in no cluster).",
        "0 to 1", "not applicable", "Context: noise genes cannot be recovered by any cluster.",
        higher=False)

# --------------------------------------------------------------------------- values
T_VALUES = "values"
_metric("spearman", "Spearman rho", T_VALUES,
        "Rank correlation between predicted and hidden measured values.",
        "-1 to 1", "0 (standard deviation 1/sqrt(n - 1))",
        "Whether the ORDER of the genes is predicted, whatever the scale. The headline for values, "
        "robust to outliers and to a prediction on a different scale.")
_metric("pearson", "Pearson r", T_VALUES,
        "Linear correlation between predicted and hidden values.",
        "-1 to 1", "0",
        "Like Spearman but on the values themselves, so a few extreme genes can dominate it.")
_metric("kendall", "Kendall tau-b", T_VALUES,
        "Share of gene pairs ordered the same way by prediction and truth, minus the share ordered "
        "oppositely, corrected for ties.",
        "-1 to 1", "0",
        "The most directly interpretable rank agreement: (1 + tau) / 2 is the chance a random "
        "pair is ordered correctly.")
_metric("r2", "R-squared (out of sample)", T_VALUES,
        "1 minus the squared prediction error over the variance of the hidden values around their "
        "own mean.",
        "below 0 to 1", "0 or below (predicting the mean scores 0)",
        "The share of the hidden values' variance the prediction explains. Negative when the "
        "prediction is further off than the mean -- common for a good ranking on the wrong scale.")
_metric("nrmse", "Normalised RMSE", T_VALUES,
        "Root-mean-square error divided by the standard deviation of the hidden values.",
        "0 upward", "1 (predicting the mean)",
        "Typical error in units of the measurement's own spread: below 1 beats the mean.",
        higher=False)
_metric("mae", "Mean absolute error", T_VALUES,
        "Mean absolute difference between predicted and hidden values, in the measurement's "
        "units.",
        "0 upward", "the mean absolute deviation of the hidden values",
        "Typical size of an error, in the same units as the data.", higher=False)
_metric("top_decile_recall", "Top-decile recall", T_VALUES,
        "Share of the genes in the true top 10% of hidden values that are also in the predicted "
        "top 10%.",
        "0 to 1", "0.1",
        "Whether the extreme genes -- usually the interesting ones -- are predicted extreme.")
_metric("bottom_decile_recall", "Bottom-decile recall", T_VALUES,
        "The same for the bottom 10% -- for fitness scores, the most essential genes.",
        "0 to 1", "0.1", "Whether the genes at the other extreme are found.")
_metric("value_coverage", "Coverage", T_VALUES,
        "Share of hidden values that received a prediction.",
        "0 to 1", "not applicable", "Every other value metric is computed on these genes.")

# --------------------------------------------------------------------------- replication
T_REPL = "replication"
_metric("replication_rate", "Replication rate", T_REPL,
        "Share of the findings made on one half of the genes that hold on the other half.",
        "0 to 1", "the rate with the second half's evidence scrambled",
        "Whether the findings are properties of the genes or of the sample they were found in.")
_metric("findings", "Findings made", T_REPL,
        "Number of findings made on the first half.",
        "0 upward", "not applicable", "Context: how much there was to replicate.")
_metric("replicated", "Findings replicated", T_REPL,
        "Number of those findings that held on the second half.",
        "0 to findings", "findings x the null rate",
        "The findings worth reading first; each is a claim that held twice.")
_metric("null_rate", "Replication by chance", T_REPL,
        "The replication rate with the second half's evidence scrambled, averaged over runs.",
        "0 to 1", "is itself the chance level",
        "What replication looks like for findings with nothing behind them.", higher=False)
_metric("replication_lift", "Replication lift", T_REPL,
        "Replication rate divided by the replication rate by chance.",
        "0 upward", "1", "How many times more often real findings replicate than chance ones.")

METRICS: dict = {m.key: m for m in _M}

TASKS: dict = {
    T_LABEL: Task(T_LABEL, "Call a label for genes that lack it.",
                  "a share of a label's genes, whole orthogroups at a time",
                  tuple(m.key for m in _M if m.task == T_LABEL)),
    T_RANK: Task(T_RANK, "Rank candidates so the true ones come first.",
                 "set members, edges or corrupted labels, ranked among negatives",
                 tuple(m.key for m in _M if m.task == T_RANK)),
    T_SET: Task(T_SET, "Return a set of genes that belong with a query.",
                "part of a gene set, to be returned among all other genes",
                tuple(m.key for m in _M if m.task == T_SET)),
    T_CLUSTER: Task(T_CLUSTER, "Find clusters that correspond to a label nobody showed them.",
                    "a share of a label, scored against clusters chosen on the rest",
                    tuple(m.key for m in _M if m.task == T_CLUSTER)),
    T_VALUES: Task(T_VALUES, "Predict a measured value for genes without one.",
                   "a share of a measurement's values",
                   tuple(m.key for m in _M if m.task == T_VALUES)),
    T_REPL: Task(T_REPL, "Make findings that hold beyond the genes they were made on.",
                 "half of the genes, on which first-half findings are checked",
                 tuple(m.key for m in _M if m.task == T_REPL)),
}

#: The verdict block every scorecard starts with, whatever the task.
VERDICT = (
    ("verdict", "Verdict", "PASS: above the 95th percentile of the null AND by the stated margin. "
     "FAIL otherwise. INCONCLUSIVE when too little could be hidden to score."),
    ("metric", "Judged on", "The one metric the verdict rests on, chosen per strategy as the honest "
     "test of its claim."),
    ("observed", "Observed", "That metric on the hidden genes."),
    ("chance", "Chance", "The same metric for the same procedure on shuffled labels, random sets or "
     "permuted identities -- measured, not assumed."),
    ("bar", "Bar", "The 95th percentile of the null runs: what luck reaches one time in twenty."),
    ("p_value", "p", "Share of null runs at least as good as the observed, (1 + k) / (1 + runs)."),
    ("skill", "Skill", "(observed - chance) / (1 - chance): 0 is chance, 1 is perfect, negative is "
     "worse than chance. Puts every metric on one scale."),
    ("n_hidden", "Hidden", "How many genes, pairs or findings the test scored."),
)


def glossary(task: str | None = None) -> pd.DataFrame:
    """Every metric (or one task's), in scorecard order, with its definition and how to read it."""
    rows = [{"task": "every test", "key": k, "metric": label, "definition": text,
             "range": "", "chance": "", "reading": "", "higher_is_better": True}
            for k, label, text in VERDICT] if task is None else []
    for t in ([task] if task else list(TASKS)):
        for k in TASKS[t].metrics:
            m = METRICS[k]
            rows.append({"task": m.task, "key": m.key, "metric": m.label,
                         "definition": m.definition, "range": m.range, "chance": m.chance,
                         "reading": m.reading, "higher_is_better": m.higher_is_better})
    return pd.DataFrame(rows)


def explain(key: str) -> str:
    """One metric's definition, range, chance level and reading, as a paragraph."""
    for k, label, text in VERDICT:
        if k == key:
            return f"{label}: {text}"
    m = METRICS[key]
    return (f"{m.label}: {m.definition} Range {m.range}; chance {m.chance}. {m.reading}")


def frame(card: dict, task: str) -> pd.DataFrame:
    """A scorecard as a table in the task's order: metric, value, and its one-line reading."""
    rows = []
    for k in TASKS[task].metrics:
        m = METRICS[k]
        rows.append({"metric": m.label, "key": k, "value": card.get(k, NAN),
                     "chance": m.chance, "reading": m.reading})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- the arithmetic
def _finite(x) -> float:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return NAN
    return x if math.isfinite(x) else NAN


def _ordered(task: str, values: dict) -> dict:
    return {k: _finite(values.get(k, NAN)) for k in TASKS[task].metrics}


def _auroc(score, positive) -> float:
    from scipy.stats import rankdata
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    ok = np.isfinite(score)
    score, positive = score[ok], positive[ok]
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if not n_pos or not n_neg:
        return NAN
    r = rankdata(score)
    return float((r[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _average_precision(score, positive) -> float:
    """Average precision with ties broken pessimistically-neutral (stable by score, then input)."""
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    ok = np.isfinite(score)
    score, positive = score[ok], positive[ok]
    if not positive.any():
        return NAN
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(positive, score))


def ranking(score, positive) -> dict:
    """The ranking scorecard for `score` (higher = more likely positive) against `positive`.

    Genes or pairs with a non-finite score are ranked last, tied: a strategy that cannot score an
    item has not ranked it high.
    """
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    if not len(score) or not positive.any() or positive.all():
        return _ordered(T_RANK, {"prevalence": positive.mean() if len(positive) else NAN})
    low = np.nanmin(score[np.isfinite(score)]) - 1.0 if np.isfinite(score).any() else 0.0
    s = np.where(np.isfinite(score), score, low)
    n, n_pos = len(s), int(positive.sum())
    prev = n_pos / n
    order = np.argsort(-s, kind="stable")
    hits = positive[order]
    cum = np.cumsum(hits)
    depth = np.arange(1, n + 1)
    prec, rec = cum / depth, cum / n_pos
    f1 = np.where(prec + rec > 0, 2 * prec * rec / np.maximum(prec + rec, 1e-12), 0.0)
    top1 = max(1, int(math.ceil(0.01 * n)))
    top10 = max(1, int(math.ceil(0.10 * n)))
    ap = _average_precision(s, positive)
    gains = hits / np.log2(depth + 1)
    ideal = (1.0 / np.log2(np.arange(1, n_pos + 1) + 1)).sum()
    try:
        from sklearn.metrics import roc_auc_score
        pauc = float(roc_auc_score(positive, s, max_fpr=0.1))
    except ValueError:
        pauc = NAN
    return _ordered(T_RANK, {
        "auroc": _auroc(s, positive), "auprc": ap, "auprc_lift": ap / prev,
        "prevalence": prev, "pauroc_10": pauc, "r_precision": float(prec[n_pos - 1]),
        "precision_at_1pct": float(prec[top1 - 1]),
        "enrichment_at_1pct": float(prec[top1 - 1] / prev),
        "recall_at_10pct": float(rec[top10 - 1]), "max_f1": float(f1.max()),
        "ndcg": float(gains.sum() / ideal)})


def label_calls(pred, truth, positions, class_scores: pd.DataFrame | None = None) -> dict:
    """The label-call scorecard on `positions`: `pred` may abstain (NaN); `truth` must be labelled.

    `class_scores` (genes x classes, higher = more likely) adds macro AUROC and AUPRC; without it
    those two are missing, not estimated.
    """
    positions = np.asarray(positions, dtype=int)
    p = pd.Series(pred).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    t = pd.Series(truth).reset_index(drop=True).iloc[positions].to_numpy(dtype=object)
    keep = np.array([isinstance(x, str) for x in t])
    p, t, positions = p[keep], t[keep], positions[keep]
    n = len(t)
    if not n:
        return _ordered(T_LABEL, {})
    called = np.array([isinstance(x, str) for x in p])
    correct = called & (p == t)
    NONE = "\x00no call"
    pc = np.where(called, p, NONE).astype(object)
    classes = sorted(set(t), key=str)
    precs, recs, f1s, sup = [], [], [], []
    for c in classes:
        is_c, call_c = t == c, pc == c
        tp = int((is_c & call_c).sum())
        prec = tp / int(call_c.sum()) if call_c.any() else 0.0
        rec = tp / int(is_c.sum())
        precs.append(prec)
        recs.append(rec)
        f1s.append(2 * prec * rec / (prec + rec) if prec + rec else 0.0)
        sup.append(int(is_c.sum()))
    sup = np.asarray(sup, dtype=float)
    # Kappa and MCC over the full confusion matrix, with "no call" as a predicted category.
    cats = sorted(set(classes) | set(pc), key=str)
    idx = {c: i for i, c in enumerate(cats)}
    C = np.zeros((len(cats), len(cats)))
    for a, b in zip(t, pc):
        C[idx[a], idx[b]] += 1
    po = np.trace(C) / n
    row, col = C.sum(axis=1), C.sum(axis=0)
    pe = float((row * col).sum()) / n ** 2
    kappa = (po - pe) / (1 - pe) if pe < 1 else NAN
    cov_ytyp = np.trace(C) * n - (row * col).sum()
    cov_ypyp = n ** 2 - (col * col).sum()
    cov_ytyt = n ** 2 - (row * row).sum()
    mcc = cov_ytyp / math.sqrt(cov_ypyp * cov_ytyt) if cov_ypyp > 0 and cov_ytyt > 0 else 0.0
    out = {"accuracy": correct.mean(), "coverage": called.mean(),
           "precision_of_calls": correct[called].mean() if called.any() else NAN,
           "macro_precision": float(np.mean(precs)), "macro_recall": float(np.mean(recs)),
           "macro_f1": float(np.mean(f1s)),
           "weighted_f1": float((np.asarray(f1s) * sup).sum() / sup.sum()),
           "kappa": kappa, "mcc": mcc}
    if class_scores is not None and len(class_scores):
        S = pd.DataFrame(class_scores).reset_index(drop=True)
        au, ap = [], []
        for c in classes:
            if c not in S.columns:
                continue
            s = S[c].iloc[positions].to_numpy(dtype=float)
            y = t == c
            if y.all() or not y.any():
                continue
            au.append(_auroc(np.nan_to_num(s, nan=-np.inf), y))
            ap.append(_average_precision(np.nan_to_num(s, nan=-1e300), y))
        au = [a for a in au if np.isfinite(a)]
        ap = [a for a in ap if np.isfinite(a)]
        out["macro_auroc"] = float(np.mean(au)) if au else NAN
        out["macro_auprc"] = float(np.mean(ap)) if ap else NAN
    return _ordered(T_LABEL, out)


def set_retrieval(returned, members, universe) -> dict:
    """Returned genes against hidden members, over a candidate universe (both subsets of it)."""
    universe = np.unique(np.asarray(universe, dtype=int))
    r = np.isin(universe, np.asarray(returned, dtype=int))
    m = np.isin(universe, np.asarray(members, dtype=int))
    n, tp = len(universe), int((r & m).sum())
    fp, fn = int((r & ~m).sum()), int((~r & m).sum())
    tn = n - tp - fp - fn
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else NAN
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    denom = math.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    share = m.mean() if n else NAN
    return _ordered(T_SET, {
        "precision": prec, "recall": rec, "f1": f1,
        "jaccard": tp / (tp + fp + fn) if tp + fp + fn else NAN,
        "set_mcc": (tp * tn - fp * fn) / denom if denom else 0.0,
        "fold_enrichment": prec / share if share else NAN, "returned": int(r.sum())})


def cluster_recovery(cluster_of_hidden, truth_of_hidden, chosen: dict, noise: int = -1) -> dict:
    """Clusters (one per hidden gene) against hidden labels, with each label's chosen cluster."""
    cl = np.asarray(cluster_of_hidden)
    t = np.asarray(truth_of_hidden, dtype=object)
    keep = np.array([isinstance(x, str) for x in t])
    cl, t = cl[keep], t[keep]
    if not len(t):
        return _ordered(T_CLUSTER, {})
    wf = wp = wr = w = 0.0
    for lab in sorted(set(t), key=str):
        is_l = t == lab
        n_l = int(is_l.sum())
        w += n_l
        k = chosen.get(lab)
        if k is None:
            continue
        in_k = cl == k
        tp = int((is_l & in_k).sum())
        prec = tp / int(in_k.sum()) if in_k.any() else 0.0
        rec = tp / n_l
        wp += n_l * prec
        wr += n_l * rec
        wf += n_l * (2 * prec * rec / (prec + rec) if prec + rec else 0.0)
    from sklearn import metrics as skm
    # Noise genes are each their own singleton: HDBSCAN did not group them, so they must not be
    # scored as one large cluster that happens to agree with itself.
    lab_cl = np.where(cl == noise, -(np.arange(len(cl)) + 2), cl)
    h, c, _v = skm.homogeneity_completeness_v_measure(t.astype(str), lab_cl)
    return _ordered(T_CLUSTER, {
        "weighted_f1_clusters": wf / w, "weighted_precision_clusters": wp / w,
        "weighted_recall_clusters": wr / w,
        "ari": skm.adjusted_rand_score(t.astype(str), lab_cl),
        "nmi": skm.normalized_mutual_info_score(t.astype(str), lab_cl),
        "homogeneity": h, "completeness": c, "noise_share": float((cl == noise).mean())})


def values(pred, truth) -> dict:
    """Predicted against hidden measured values; NaN predictions lower coverage, nothing else."""
    from scipy import stats
    pred = np.asarray(pred, dtype=float)
    truth = np.asarray(truth, dtype=float)
    have = np.isfinite(truth)
    pred, truth = pred[have], truth[have]
    if not len(truth):
        return _ordered(T_VALUES, {})
    ok = np.isfinite(pred)
    out = {"value_coverage": ok.mean()}
    p, y = pred[ok], truth[ok]
    if len(y) >= 3 and np.std(y) > 0 and np.std(p) > 0:
        err = p - y
        k = max(1, int(round(0.1 * len(y))))
        top_t, top_p = set(np.argsort(-y)[:k]), set(np.argsort(-p)[:k])
        bot_t, bot_p = set(np.argsort(y)[:k]), set(np.argsort(p)[:k])
        out.update({"spearman": stats.spearmanr(p, y).statistic,
                    "pearson": stats.pearsonr(p, y).statistic,
                    "kendall": stats.kendalltau(p, y).statistic,
                    "r2": 1.0 - float((err ** 2).sum()) / float(((y - y.mean()) ** 2).sum()),
                    "nrmse": math.sqrt(float((err ** 2).mean())) / float(np.std(y)),
                    "mae": float(np.abs(err).mean()),
                    "top_decile_recall": len(top_t & top_p) / k,
                    "bottom_decile_recall": len(bot_t & bot_p) / k})
    return _ordered(T_VALUES, out)


def replication(replicated: int, findings: int, null_rates) -> dict:
    """Findings made on one half, checked on the other, against the rate with the evidence scrambled."""
    nulls = [v for v in (null_rates or []) if np.isfinite(v)]
    null = float(np.mean(nulls)) if nulls else NAN
    rate = replicated / findings if findings else NAN
    return _ordered(T_REPL, {"replication_rate": rate, "findings": findings,
                             "replicated": replicated, "null_rate": null,
                             "replication_lift": rate / null if null and null > 0 else NAN})
