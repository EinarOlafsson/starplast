"""Nested, label-locked validation for cluster-to-annotation inference.

Embeddings and clusters are unsupervised, but choosing their hyperparameters by the same labels later
reported as validation is still leakage.  This module separates those uses: an inner validation split
chooses a structure, an untouched outer fold estimates its performance, and a cluster-to-category map
is learned only from the training labels.  Genes without a sufficiently supported/pure cluster abstain.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

ABSTAIN = "__abstain__"


def _known(values) -> np.ndarray:
    """Measured labels only; explicit absence words are unknown, not biological categories."""
    from .search import ABSENCE_LABELS
    series = pd.Series(values, dtype="object")
    text = series.astype("object").where(series.notna(), "").astype(str).str.lower()
    return series.notna().to_numpy() & ~text.isin(ABSENCE_LABELS).to_numpy()


@dataclass(frozen=True)
class ClusterCall:
    """One cluster's training-only category call and its empirical support."""
    category: str
    support: int
    probability: float


def fit_cluster_calls(labels, truth, train_index, *, min_support=5,
                      min_probability=0.55) -> dict[int, ClusterCall]:
    """Learn cluster-to-category calls from training indices only; noise never receives a call."""
    labels = np.asarray(labels, dtype=int)
    truth = np.asarray(truth, dtype=object)
    train = np.asarray(train_index, dtype=int)
    out = {}
    for cluster in np.unique(labels[train]):
        if cluster < 0:
            continue
        values = truth[train[labels[train] == cluster]]
        values = values[pd.notna(values)]
        if not len(values):
            continue
        categories, counts = np.unique(values.astype(str), return_counts=True)
        best = int(np.argmax(counts))
        support = int(counts[best])
        probability = float(support / counts.sum())
        if support >= min_support and probability >= min_probability:
            out[int(cluster)] = ClusterCall(str(categories[best]), support, probability)
    return out


def predict_cluster_calls(labels, calls, index) -> np.ndarray:
    """Predict selected indices, using an explicit abstention where the training map has no call."""
    labels = np.asarray(labels, dtype=int)
    index = np.asarray(index, dtype=int)
    return np.asarray([calls.get(int(labels[i]), ClusterCall(ABSTAIN, 0, 0.0)).category
                       for i in index], dtype=object)


def score_predictions(truth, predicted, categories) -> tuple[dict, pd.DataFrame]:
    """Macro F1 plus per-category precision/recall; abstentions count as false negatives."""
    from sklearn.metrics import precision_recall_fscore_support

    y = np.asarray(truth, dtype=str)
    p = np.asarray(predicted, dtype=str)
    categories = np.asarray(sorted(set(map(str, categories))), dtype=object)
    precision, recall, f1, support = precision_recall_fscore_support(
        y, p, labels=categories, zero_division=0)
    per = pd.DataFrame({"category": categories, "precision": precision, "recall": recall,
                        "f1": f1, "support": support.astype(int)})
    covered = p != ABSTAIN
    summary = {
        "macro_f1": float(np.mean(f1)) if len(f1) else 0.0,
        "macro_precision": float(np.mean(precision)) if len(precision) else 0.0,
        "macro_recall": float(np.mean(recall)) if len(recall) else 0.0,
        "coverage": float(covered.mean()) if len(covered) else 0.0,
        "covered_accuracy": float((p[covered] == y[covered]).mean()) if covered.any() else 0.0,
    }
    return summary, per


def nested_structure_cv(structures: dict[str, np.ndarray], truth: pd.Series, *, folds=5,
                        seed=42, min_support=5, min_probability=0.55,
                        split_labels: list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None = None):
    """Select structures on inner validation labels and score once on each untouched outer fold.

    Returns ``(outer, per_category, validation_scores, splits)``.  ``splits`` can be reused with
    permuted labels so the null repeats the complete selection procedure on identical partitions.
    """
    from sklearn.model_selection import StratifiedKFold, train_test_split

    values = truth.astype("object").to_numpy()
    labeled = np.flatnonzero(_known(values))
    y = values[labeled].astype(str)
    categories = sorted(set(y))
    if len(categories) < 2:
        raise ValueError("nested validation needs at least two target categories")
    if min(pd.Series(y).value_counts()) < folds:
        raise ValueError(f"every category needs at least {folds} labels for stratified folds")

    if split_labels is None:
        split_labels = []
        outer = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
        for fold, (train_val_pos, test_pos) in enumerate(outer.split(labeled, y)):
            train_pos, val_pos = train_test_split(
                train_val_pos, test_size=0.2, random_state=seed + fold,
                stratify=y[train_val_pos])
            split_labels.append((labeled[train_pos], labeled[val_pos], labeled[test_pos]))

    outer_rows, per_rows, validation_rows = [], [], []
    for fold, (train, validation, test) in enumerate(split_labels):
        ranked = []
        for structure_id, labels in structures.items():
            calls = fit_cluster_calls(labels, values, train, min_support=min_support,
                                      min_probability=min_probability)
            predicted = predict_cluster_calls(labels, calls, validation)
            summary, _ = score_predictions(values[validation], predicted, categories)
            row = {"fold": fold, "structure_id": structure_id, **summary}
            validation_rows.append(row)
            ranked.append(row)
        best = max(ranked, key=lambda row: (row["macro_f1"], row["coverage"],
                                            row["covered_accuracy"], row["structure_id"]))
        labels = structures[best["structure_id"]]
        calls = fit_cluster_calls(labels, values, np.r_[train, validation],
                                  min_support=min_support, min_probability=min_probability)
        predicted = predict_cluster_calls(labels, calls, test)
        summary, per = score_predictions(values[test], predicted, categories)
        outer_rows.append({"fold": fold, "selected_structure": best["structure_id"],
                           "n_train": len(train) + len(validation), "n_test": len(test), **summary})
        per.insert(0, "fold", fold)
        per.insert(1, "selected_structure", best["structure_id"])
        per_rows.append(per)
    return (pd.DataFrame(outer_rows), pd.concat(per_rows, ignore_index=True),
            pd.DataFrame(validation_rows), split_labels)


def choose_final_structure(validation_scores: pd.DataFrame) -> str:
    """Choose the structure with the best mean inner-fold score, never the outer test score."""
    mean = validation_scores.groupby("structure_id", as_index=False)[
        ["macro_f1", "coverage", "covered_accuracy"]].mean()
    mean = mean.sort_values(["macro_f1", "coverage", "covered_accuracy", "structure_id"],
                            ascending=[False, False, False, True])
    if mean.empty:
        raise ValueError("no validation scores")
    return str(mean.iloc[0].structure_id)


def pooled_category_scores(per_category: pd.DataFrame) -> pd.DataFrame:
    """Mean and spread of locked outer-fold category metrics."""
    return (per_category.groupby("category", as_index=False)
            .agg(precision=("precision", "mean"), recall=("recall", "mean"),
                 f1=("f1", "mean"), f1_sd=("f1", "std"), support=("support", "sum"))
            .sort_values("f1", ascending=False).reset_index(drop=True))


def inference_table(labels, truth: pd.Series, gene_ids, category_cv: pd.DataFrame, *,
                    min_support=5, min_probability=0.55) -> pd.DataFrame:
    """Call unlabeled genes from the final structure, carrying cluster and locked-CV confidence."""
    values = truth.astype("object").to_numpy()
    known_mask = _known(values)
    known = np.flatnonzero(known_mask)
    unknown = np.flatnonzero(~known_mask)
    calls = fit_cluster_calls(labels, values, known, min_support=min_support,
                              min_probability=min_probability)
    cv = category_cv.set_index("category") if not category_cv.empty else pd.DataFrame()
    rows = []
    for index in unknown:
        cluster = int(np.asarray(labels)[index])
        call = calls.get(cluster)
        if call is None:
            continue
        row = {"gene_id": str(np.asarray(gene_ids)[index]), "proposed": call.category,
               "cluster": cluster, "cluster_training_support": call.support,
               "cluster_training_probability": call.probability,
               "evidence_status": "inferred_from_held_out_structure"}
        if not cv.empty and call.category in cv.index:
            row.update({"nested_cv_precision": float(cv.loc[call.category, "precision"]),
                        "nested_cv_recall": float(cv.loc[call.category, "recall"]),
                        "nested_cv_f1": float(cv.loc[call.category, "f1"])})
        rows.append(row)
    return pd.DataFrame(rows)
