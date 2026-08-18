#!/usr/bin/env python3
"""Ways of answering a recipe other than clustering a UMAP, judged by the same control.

Instruction 47. Until now every claim this program made came out of one procedure: embed a feature
matrix, cluster it, and see whether a held-out label comes back. That is a deliberately conservative
design and it stays. But it has no BASELINE -- "mean F1 0.20" cannot be called good or bad by anyone
-- and it flattens away 302,756 edges that never inform a prediction.

**The unification that makes the comparison honest: every method produces a PARTITION.** For
UMAP + HDBSCAN the partition is the clustering. For a classifier it is the out-of-fold predicted
class. For propagation it is the label whose diffusion score is highest. Once a method has produced
one, everything downstream is unchanged -- `score_recovery` measures it per label, `infer` names the
unlabelled genes in concentrated groups, and `dominant_by_cluster` asks whether the independent
control corroborates them. One question, several methods, one control judging all of them.

**Supervised methods must be scored OUT OF FOLD, and this is the whole difficulty.** A classifier's
partition is derived from the label it is being scored against, so scoring its training predictions
would report memorisation as recovery -- a mistake far easier to make here than in the clustering,
which never sees the label at all. Every supervised partition below is assembled from cross-validated
predictions, each gene classified by a model that never saw it.

The leakage closure matters MORE here, not less. A UMAP blurs a 0.7-associated column; a penalised
regression puts a coefficient on it and a boosted tree splits on it twice.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: Folds for every cross-validated method, and the same folds for all of them so two methods on one
#: question are compared on identical splits rather than on their own luck.
FOLDS = 5

#: Below this many labelled genes a class cannot be learned OR scored, matching the floor the
#: clustering path already uses. A class with eight members is fitted by luck.
MIN_CLASS = 15


def _usable(truth: pd.Series, min_class: int = MIN_CLASS):
    """The genes and classes a supervised method may be trained and scored on.

    Rare classes are dropped rather than kept, and this is a real choice: keeping them lets a model
    score well by never predicting them, and the per-label table then shows a row of zeros that reads
    as a failure of the biology rather than of the design.
    """
    values = truth.astype("object").where(truth.notna(), None).to_numpy()
    known = np.array([v is not None for v in values], dtype=bool)
    counts = pd.Series([str(v) for v in values[known]]).value_counts()
    keep = set(counts[counts >= min_class].index)
    mask = np.array([v is not None and str(v) in keep for v in values], dtype=bool)
    return mask, np.array([str(v) for v in values], dtype=object)


def out_of_fold(X: np.ndarray, truth: pd.Series, model, seed: int = 42,
                folds: int = FOLDS, min_class: int = MIN_CLASS) -> tuple:
    """Cross-validated predictions for the labelled genes, and a full-data model for the rest.

    Returns ``(partition, fitted, classes)``. `partition` carries the out-of-fold predicted class for
    every labelled gene and the full model's prediction for every unlabelled one, as integer codes --
    the same shape a clustering produces, so the scoring path does not know which method made it.

    The unlabelled genes are predicted by a model trained on ALL the labelled ones, which is correct
    and is worth saying plainly: they are the answer, they were never in any fold, and holding a
    fold out of the model that names them would only make that answer worse.
    """
    from sklearn.model_selection import StratifiedKFold
    mask, values = _usable(truth, min_class)
    if mask.sum() < folds * 2:
        return np.full(len(X), -1), None, []
    y = values[mask]
    classes = sorted(set(y))
    codes = {c: i for i, c in enumerate(classes)}
    yi = np.array([codes[v] for v in y])
    partition = np.full(len(X), -1)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    index = np.arange(len(X))[mask]
    for train, test in splitter.split(X[mask], yi):
        model.fit(X[mask][train], yi[train])
        partition[index[test]] = model.predict(X[mask][test])
    model.fit(X[mask], yi)
    rest = ~mask
    if rest.any():
        partition[rest] = model.predict(X[rest])
    return partition, model, classes


def logistic(X: np.ndarray, truth: pd.Series, seed: int = 42, folds: int = FOLDS,
             C: float = 0.1) -> dict:
    """L1-penalised logistic regression: the baseline this project has never had.

    Chosen as the first alternative method for three reasons, and the third is the one that matters.
    It yields a COEFFICIENT PER COLUMN, so an answer reads "this label is predicted by these six
    measurements" -- something a cluster can never say. It has no cluster-size floor, so labels too
    rare to cluster become answerable. And it is a yardstick: if a penalised linear model recovers a
    label better than the tuned UMAP and its tuned clustering, then the map is a picture rather than
    an inference engine, and everything more elaborate should be judged against this number rather
    than against nothing.

    L1 rather than L2 because the answer wanted is WHICH measurements, not a weight on all 361 of
    them; `class_weight="balanced"` because these labels are wildly unequal and an unweighted fit
    predicts the common class and calls it accuracy.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    # One-vs-rest, stated rather than inherited. `liblinear` used to apply this scheme implicitly for
    # multiclass problems and scikit-learn is removing that -- it warns now and raises in 1.8. Wrapping
    # it keeps the same fit and the same per-class coefficients, which is what makes this method worth
    # having; switching to a multinomial solver instead would silently change every number here.
    model = make_pipeline(
        StandardScaler(),
        OneVsRestClassifier(
            LogisticRegression(penalty="l1", C=C, solver="liblinear", class_weight="balanced",
                               max_iter=2000, random_state=seed)))
    partition, fitted, classes = out_of_fold(X, truth, model, seed=seed, folds=folds)
    return {"partition": partition, "classes": classes, "model": fitted,
            "settings": {"method": "logistic", "penalty": "l1", "C": C, "folds": folds}}


def coefficients(fitted, classes, feature_names) -> pd.DataFrame:
    """Which measurements the model actually used, per class, strongest first.

    The reason a linear baseline earns its place. A cluster says "these genes go together"; this says
    "this label is predicted by these columns, in this direction", which is a sentence a biologist
    can disagree with -- and disagreement is what makes a claim testable.
    """
    if fitted is None:
        return pd.DataFrame()
    clf = fitted[-1] if hasattr(fitted, "__getitem__") else fitted
    # OneVsRestClassifier keeps a model per class rather than one `coef_`, so the per-class weights
    # are stacked back into the shape the rest of this function expects.
    coef = (np.vstack([np.ravel(e[-1].coef_ if hasattr(e, "__getitem__") else e.coef_)
                       for e in clf.estimators_])
            if hasattr(clf, "estimators_") else np.atleast_2d(clf.coef_))
    rows = []
    for i, weights in enumerate(coef):
        label = classes[i] if len(coef) > 1 else f"{classes[1]} vs {classes[0]}"
        for j in np.argsort(-np.abs(weights))[:20]:
            if weights[j] == 0:
                continue
            rows.append({"label": label, "feature": feature_names[j],
                         "coefficient": float(weights[j])})
    return pd.DataFrame(rows)
