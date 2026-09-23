"""Gene prediction with held-out groups, train-only preprocessing and explicit support.

UMAP is an optional intermediate representation, independent of the display map.
Unknown labels are never negatives. Each outer fold reserves separate training
and calibration groups; scores and uncertainty on test genes use neither their
labels nor their fitted feature distribution. Network runs are explicitly
transductive over a fixed, target-independent graph.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import softmax
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             f1_score, log_loss, mean_absolute_error, mean_squared_error, r2_score)
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold, StratifiedKFold, KFold
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor

METHODS = ("prior", "linear", "boosted", "neighbors", "pca", "umap", "multiview", "network")


@dataclass
class TaskSpec:
    """One outcome and validation contract; all settings are saved with a result."""
    target: str
    kind: str = "classification"
    method: str = "linear"
    features: tuple = ()
    group_column: str | None = "orthogroup"
    folds: int = 3
    seed: int = 42
    dimensions: int = 20
    neighbors: int = 15
    min_feature_coverage: float = .05
    min_gene_coverage: float = .1
    min_probability: float = 0.0
    calibrate: bool = True
    interval_alpha: float = .1
    exclude: tuple = ()
    blocks: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.kind not in {"classification", "regression"} or self.method not in METHODS:
            raise ValueError("unsupported prediction task or method")
        if self.folds < 2 or self.dimensions < 1 or self.neighbors < 1:
            raise ValueError("folds >= 2, dimensions >= 1 and neighbors >= 1 required")
        if not 0 < self.interval_alpha < 1:
            raise ValueError("interval_alpha must be between zero and one")
        for name in ('min_feature_coverage', 'min_gene_coverage', 'min_probability'):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be between zero and one")
        if self.method == "network" and self.kind != "classification":
            raise ValueError("network propagation currently supports classification")
        if len(self.features) != len(set(self.features)):
            raise ValueError("feature names must be unique")


@dataclass
class PredictionResult:
    """Held-out evaluations and unknown-gene hypotheses with reproducible provenance."""
    spec: TaskSpec
    predictions: pd.DataFrame
    metrics: dict
    per_class: pd.DataFrame
    folds: list
    provenance: dict
    model: object = field(default=None, repr=False)

    def save(self, directory):
        """Write ordinary CSV/JSON artifacts; fitted Python objects are not pickled."""
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        self.predictions.to_csv(out / "predictions.csv", index=False)
        self.per_class.to_csv(out / "per_class.csv", index=False)
        payload = {"spec": asdict(self.spec), "metrics": self.metrics,
                   "folds": self.folds, "provenance": self.provenance}
        (out / "run.json").write_text(json.dumps(payload, indent=2, default=_json, allow_nan=False) + "\n")
        return out


def _json(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def known_labels(series):
    """Boolean mask of observed categorical labels, preserving measured 0/False."""
    from .holdout_cv import _known
    return _known(series.to_numpy(dtype=object))


def group_ids(nodes, column):
    """Keep related genes together; genes with no group receive unique singleton IDs."""
    if column is None:
        return np.array([f"gene:{gene}" for gene in nodes.gene_id])
    if column not in nodes:
        raise ValueError(f"group column {column!r} is absent; choose an available grouping or random folds")
    values = nodes[column].astype(object)
    known = known_labels(values)
    return np.array([f"group:{v}" if ok else f"gene:{gene}"
                     for v, ok, gene in zip(values, known, nodes.gene_id)])


class FoldMatrix(BaseEstimator, TransformerMixin):
    """Fit feature filtering, median imputation, robust scaling and view weights on training genes."""
    def __init__(self, min_coverage=.05, blocks=None):
        self.min_coverage = min_coverage
        self.blocks = blocks

    def fit(self, X, y=None):
        """Learn columns and all normalization statistics from this training frame only."""
        values = X.to_numpy(dtype=float, copy=True)
        finite = np.isfinite(values)
        self.columns_ = list(X.columns[(finite.mean(axis=0) >= self.min_coverage) & (finite.sum(axis=0) >= 2)])
        if not self.columns_:
            raise ValueError("no measured training features remain")
        raw = X[self.columns_].to_numpy(dtype=float, copy=True)
        raw[~np.isfinite(raw)] = np.nan
        self.median_ = np.nanmedian(raw, axis=0)
        low, high = np.nanpercentile(raw, [25, 75], axis=0)
        self.scale_ = np.where(high - low > 1e-9, high - low, 1.)
        self.weights_ = np.ones(len(self.columns_))
        scaled = np.clip(np.nan_to_num((raw-self.median_) / self.scale_), -20, 20)
        # One column belongs to one view. Explicit overlapping view assignments
        # are refused rather than counted twice.
        assignments = {}
        for block, cols in (self.blocks or {}).items():
            for col in cols:
                if col in assignments:
                    raise ValueError(f"feature {col} belongs to multiple views")
                assignments[col] = block
        labels = np.array([assignments.get(c, _view(c)) for c in self.columns_])
        self.view_columns_ = {}
        for label in sorted(set(labels)):
            indices = np.flatnonzero(labels == label)
            variance = float(scaled[:, indices].var(axis=0).sum())
            self.weights_[indices] = 1 / np.sqrt(variance) if variance > 1e-12 else 1.
            self.view_columns_[label] = indices
        return self

    def transform(self, X):
        """Apply the frozen training transform; missing values map to training medians."""
        values = X[self.columns_].to_numpy(dtype=float, copy=True)
        values[~np.isfinite(values)] = np.nan
        return np.clip(np.nan_to_num((values-self.median_) / self.scale_), -20, 20) * self.weights_

    def coverage(self, X):
        """Fraction of retained measurements actually observed for each gene."""
        return np.isfinite(X[self.columns_].to_numpy(dtype=float, copy=True)).mean(axis=1)


def _view(column):
    from .datasets import provenance
    if column.startswith("esm_"):
        return "protein_sequence"
    source = provenance(column)
    return source.key if source else column.split("_", 1)[0]


class MaskedFactors(BaseEstimator, TransformerMixin):
    """Train a low-rank representation by reconstructing only observed entries.

    This is a masked matrix-completion baseline, not MOFA+. View balancing is
    supplied by FoldMatrix. New rows infer factors from observed features only.
    """
    def __init__(self, components=20, iterations=12, ridge=.1, seed=42):
        self.components = components
        self.iterations = iterations
        self.ridge = ridge
        self.seed = seed

    def fit(self, X, y=None):
        """Alternate low-rank reconstruction and replacement of missing entries."""
        mask = np.isfinite(X)
        filled = np.nan_to_num(X)
        count = min(self.components, max(1, min(filled.shape)-1))
        for _ in range(self.iterations):
            pca = PCA(n_components=count, random_state=self.seed, svd_solver="randomized")
            factors = pca.fit_transform(filled)
            reconstructed = pca.inverse_transform(factors)
            filled = np.where(mask, X, reconstructed)
        self.mean_, self.components_ = pca.mean_, pca.components_
        return self

    def transform(self, X):
        """Solve regularized factors without using unknown feature entries."""
        result = np.zeros((len(X), len(self.components_)))
        identity = np.eye(len(self.components_)) * self.ridge
        for i, row in enumerate(X):
            observed = np.isfinite(row)
            basis = self.components_[:, observed]
            if observed.any():
                result[i] = np.linalg.solve(basis @ basis.T + identity,
                                           basis @ (row[observed] - self.mean_[observed]))
        return result


def _splits(indices, y, groups, spec, n_splits=None):
    count = min(n_splits or spec.folds, len(np.unique(groups[indices])))
    if count < 2:
        raise ValueError("at least two independent groups are required")
    if spec.kind == "classification":
        splitter = StratifiedGroupKFold(n_splits=count, shuffle=True, random_state=spec.seed)
        pairs = splitter.split(indices, y[indices], groups[indices])
    else:
        # Randomize group names before deterministic GroupKFold for old sklearn versions.
        unique = np.unique(groups[indices]); shuffled = unique.copy()
        np.random.default_rng(spec.seed).shuffle(shuffled)
        order = dict(zip(unique, shuffled))
        pairs = GroupKFold(n_splits=count).split(indices, y[indices], [order[g] for g in groups[indices]])
    return [(indices[a], indices[b]) for a,b in pairs]


def _training_columns(nodes, indices, spec):
    from .search import excluded_detail
    # Association is measured only inside the fitting subset, never test/calibration.
    banned = excluded_detail(nodes.iloc[indices], spec.target, scope="target_family")
    banned.update({name: "explicit exclusion" for name in spec.exclude})
    banned[spec.target] = "target"
    if spec.group_column:
        banned[spec.group_column] = "validation group"
    candidates = spec.features or tuple(nodes.select_dtypes(include=np.number).columns)
    columns = [c for c in candidates if c in nodes and c not in banned
               and pd.api.types.is_numeric_dtype(nodes[c])]
    # Literature attention is a sampling bias, not a default biological predictor.
    columns = [c for c in columns if not c.startswith(("n_publications", "n_fulltext", "n_papers_"))]
    if not columns:
        raise ValueError("no numeric features remain after target leakage exclusions")
    return columns, banned


class _Fitted:
    def __init__(self, nodes, y, indices, spec, classes):
        self.spec, self.classes = spec, classes
        self.columns, self.excluded = _training_columns(nodes, indices, spec)
        self.transformer = FoldMatrix(spec.min_feature_coverage, spec.blocks).fit(nodes.iloc[indices][self.columns])
        self.reducer = None
        self.model = self._model(nodes, y, indices)
        self.temperature = 1.
        self.interval_radius = None
        self.calibration_status = "not_calibrated"

    def _matrix(self, nodes):
        X = self.transformer.transform(nodes[self.columns])
        if self.spec.method == "multiview":
            raw = nodes[self.transformer.columns_].to_numpy(dtype=float, copy=True)
            X[~np.isfinite(raw)] = np.nan
        return self.reducer.transform(X) if self.reducer is not None else X

    def _model(self, nodes, y, indices):
        spec = self.spec
        X = self._matrix(nodes.iloc[indices])
        if spec.method in {"pca", "umap", "multiview"}:
            count = min(spec.dimensions, max(1, min(X.shape)-1))
            if spec.method == "pca":
                self.reducer = PCA(n_components=count, random_state=spec.seed)
            elif spec.method == "multiview":
                self.reducer = MaskedFactors(count, seed=spec.seed)
            else:
                import umap
                self.reducer = umap.UMAP(n_components=count, n_neighbors=min(spec.neighbors, len(X)-1),
                                         min_dist=0., metric="euclidean", random_state=spec.seed,
                                         init="random", n_jobs=1)
            self.reducer.fit(X)
            X = self.reducer.transform(X)
        classifier = spec.kind == "classification"
        if spec.method == "prior" or (classifier and len(np.unique(y[indices])) < 2):
            model = DummyClassifier(strategy="prior") if classifier else DummyRegressor()
        elif spec.method == "boosted":
            model = (HistGradientBoostingClassifier if classifier else HistGradientBoostingRegressor)(
                max_iter=120, max_leaf_nodes=15, l2_regularization=1., random_state=spec.seed)
        elif spec.method == "linear":
            model = LogisticRegression(C=1., max_iter=1500, random_state=spec.seed) if classifier else Ridge(alpha=10.)
        else:
            model = (KNeighborsClassifier if classifier else KNeighborsRegressor)(
                n_neighbors=min(spec.neighbors, len(indices)), weights="distance")
        model.fit(X, y[indices])
        return model

    def probabilities(self, nodes):
        local = self.model.predict_proba(self._matrix(nodes))
        probabilities = np.zeros((len(nodes), len(self.classes)))
        for j, label in enumerate(self.model.classes_):
            probabilities[:, int(label)] = local[:, j]
        # Temperature scaling preserves class ordering and uses calibration labels only.
        return softmax(np.log(np.clip(probabilities, 1e-12, 1)) / self.temperature, axis=1)

    def predict(self, nodes):
        return self.model.predict(self._matrix(nodes))

    def calibrate(self, nodes, y, indices):
        if not len(indices):
            return
        if self.spec.kind == "classification":
            from scipy.optimize import minimize_scalar
            raw = self.probabilities(nodes.iloc[indices])
            if len(indices) < max(20, 2*len(self.classes)) or len(np.unique(y[indices])) < 2:
                self.calibration_status = "insufficient_calibration_labels"
                return
            logs = np.log(np.clip(raw, 1e-12, 1))
            def objective(log_temperature):
                probabilities = softmax(logs / np.exp(log_temperature), axis=1)
                return log_loss(y[indices], probabilities, labels=np.arange(len(self.classes)))
            result = minimize_scalar(objective, bounds=(-2, 3), method="bounded")
            self.temperature = float(np.exp(result.x))
            self.calibration_status = "held_out_temperature_scaling"
        else:
            residual = np.abs(y[indices] - self.predict(nodes.iloc[indices]))
            quantile = min(1., np.ceil((len(indices)+1)*(1-self.spec.interval_alpha)) / len(indices))
            self.interval_radius = float(np.quantile(residual, quantile, method="higher"))
            self.calibration_status = "held_out_absolute_residual_interval"


def _fit(nodes, y, train, spec, groups, classes):
    fit, calibration = train, np.array([], dtype=int)
    if spec.calibrate and len(np.unique(groups[train])) >= 5:
        fit, calibration = _splits(train, y, groups, spec, n_splits=5)[0]
    fitted = _Fitted(nodes, y, fit, spec, classes)
    if spec.calibrate:
        fitted.calibrate(nodes, y, calibration)
    return fitted, fit, calibration


def _evaluate(y, predicted, probabilities, known, accepted, classes, kind):
    index = np.flatnonzero(known)
    covered = index[accepted[index]]
    metrics = {"test_genes": len(index), "supported_genes": len(covered),
               "coverage": float(len(covered)/len(index)), "evaluation": "out_of_fold"}
    per_class = []
    if kind == "classification":
        # Unsupported calls count as misses in F1/recall; selective accuracy is separate.
        calls = predicted[index].astype(int).copy()
        calls[~accepted[index]] = -1
        metrics.update(accuracy=float(accuracy_score(y[index], calls)),
                       macro_f1=float(f1_score(y[index], calls, labels=np.arange(len(classes)), average="macro", zero_division=0)),
                       balanced_accuracy=float(np.mean([np.mean(calls[y[index]==c]==c) for c in range(len(classes))])),
                       log_loss=float(log_loss(y[index], probabilities[index], labels=np.arange(len(classes)))),
                       brier=float(np.mean(np.sum((probabilities[index]-np.eye(len(classes))[y[index]])**2, axis=1))))
        if len(covered):
            metrics["selective_accuracy"] = float(accuracy_score(y[covered], predicted[covered]))
        for c, label in enumerate(classes):
            truth = y[index] == c
            selected = calls == c
            tp = int(np.sum(truth & selected))
            ap = float(average_precision_score(truth, probabilities[index, c]))
            per_class.append({"label": label, "support": int(truth.sum()), "prevalence": float(truth.mean()),
                              "average_precision": ap, "precision": tp/max(1,int(selected.sum())),
                              "recall": tp/max(1,int(truth.sum()))})
        metrics["macro_average_precision"] = float(np.mean([r["average_precision"] for r in per_class]))
    else:
        from scipy.stats import spearmanr
        if len(covered):
            metrics.update(mae=float(mean_absolute_error(y[covered],predicted[covered])),
                           rmse=float(np.sqrt(mean_squared_error(y[covered],predicted[covered]))))
        if len(covered) > 1:
            metrics["r2"] = float(r2_score(y[covered], predicted[covered]))
            correlation = spearmanr(y[covered], predicted[covered]).statistic
            metrics["spearman"] = float(correlation) if np.isfinite(correlation) else None
    return metrics, pd.DataFrame(per_class)


def run(nodes, spec, log=print):
    """Evaluate held-out groups, then generate hypotheses for genuinely unlabelled genes.

    All folds share the same frozen specification. For method selection, compare
    these development results and evaluate the selected method on a separate
    future dataset before claiming prospective performance.
    """
    if spec.method == "network":
        raise ValueError("use run_network with explicit target-independent edge layers")
    if "gene_id" not in nodes or nodes.gene_id.isna().any() or nodes.gene_id.duplicated().any():
        raise ValueError("one unique nonmissing gene_id per row is required")
    if spec.target not in nodes:
        raise ValueError(f"target {spec.target!r} is absent")
    nodes = nodes.reset_index(drop=True)
    if spec.kind == "classification":
        known = known_labels(nodes[spec.target])
        classes = sorted(set(nodes.loc[known, spec.target].map(str)))
        if len(classes) < 2:
            raise ValueError("classification requires at least two observed classes")
        mapping = {label:i for i,label in enumerate(classes)}
        y = np.array([mapping.get(str(v), -1) for v in nodes[spec.target]])
    else:
        y = pd.to_numeric(nodes[spec.target], errors="coerce").to_numpy(dtype=float, copy=True)
        known = np.isfinite(y)
        classes = []
    groups = group_ids(nodes, spec.group_column)
    index = np.flatnonzero(known)
    if len(index) < 10:
        raise ValueError("at least ten labelled genes are required for evaluation")
    predicted = np.full(len(nodes), np.nan)
    probabilities = np.full((len(nodes), len(classes)), np.nan)
    coverage = np.zeros(len(nodes)); fold_ids = np.full(len(nodes), -1)
    lower = np.full(len(nodes), np.nan); upper = lower.copy()
    reports = []
    for fold, (train, test) in enumerate(_splits(index,y,groups,spec)):
        log(f"Fold {fold+1}/{spec.folds}: {len(train)} training, {len(test)} held-out genes")
        fitted, fit, calibration = _fit(nodes,y,train,spec,groups,classes)
        if set(groups[fit]) & set(groups[test]) or set(groups[calibration]) & set(groups[test]):
            raise AssertionError("validation groups overlap")
        coverage[test] = fitted.transformer.coverage(nodes.iloc[test])
        if classes:
            probabilities[test] = fitted.probabilities(nodes.iloc[test])
            predicted[test] = probabilities[test].argmax(axis=1)
        else:
            predicted[test] = fitted.predict(nodes.iloc[test])
            if fitted.interval_radius is not None:
                lower[test] = predicted[test]-fitted.interval_radius
                upper[test] = predicted[test]+fitted.interval_radius
        fold_ids[test] = fold
        reports.append({"fold":fold,"training_ids":nodes.gene_id.iloc[fit].tolist(),
                        "calibration_ids":nodes.gene_id.iloc[calibration].tolist(),
                        "test_ids":nodes.gene_id.iloc[test].tolist(),
                        "features":fitted.transformer.columns_,"excluded":fitted.excluded,
                        "calibration":fitted.calibration_status,"temperature":fitted.temperature,
                        "interval_radius":fitted.interval_radius,
                        "estimator":type(fitted.model).__name__})
    log("Fitting unknown-gene predictions using the same training/calibration contract")
    final, fit, calibration = _fit(nodes,y,index,spec,groups,classes)
    unknown = np.flatnonzero(~known)
    if len(unknown):
        coverage[unknown] = final.transformer.coverage(nodes.iloc[unknown])
        if classes:
            probabilities[unknown] = final.probabilities(nodes.iloc[unknown])
            predicted[unknown] = probabilities[unknown].argmax(axis=1)
        else:
            predicted[unknown] = final.predict(nodes.iloc[unknown])
            if final.interval_radius is not None:
                lower[unknown] = predicted[unknown]-final.interval_radius
                upper[unknown] = predicted[unknown]+final.interval_radius
    accepted = coverage >= spec.min_gene_coverage
    if classes:
        accepted &= probabilities.max(axis=1) >= spec.min_probability
    metrics, per_class = _evaluate(y,predicted,probabilities,known,accepted,classes,spec.kind)
    table = pd.DataFrame({"gene_id":nodes.gene_id,"group":groups,"fold":fold_ids,
                          "role":np.where(known,"held_out_evaluation","unlabelled_candidate"),
                          "truth":nodes[spec.target],"feature_coverage":coverage,
                          "supported":accepted,"evidence_status":"model_prediction"})
    table["prediction"] = [classes[int(v)] for v in predicted] if classes else predicted
    table.loc[~accepted,"prediction"] = None
    table["abstention_reason"] = np.where(accepted,"",np.where(coverage < spec.min_gene_coverage,
                                                             "insufficient_measured_features","below_probability_threshold"))
    if classes:
        for c,label in enumerate(classes):
            table[f"probability::{label}"] = probabilities[:,c]
        table["score"] = probabilities.max(axis=1)
    else:
        table["interval_lower"], table["interval_upper"] = lower, upper
        valid = known & accepted & np.isfinite(lower)
        if valid.any():
            metrics["interval_coverage"] = float(np.mean((y[valid]>=lower[valid])&(y[valid]<=upper[valid])))
            metrics["mean_interval_width"] = float(np.mean(upper[valid]-lower[valid]))
    table["calibration_status"] = final.calibration_status
    for report in reports:
        table.loc[table.fold.eq(report["fold"]),"calibration_status"] = report["calibration"]
    fingerprint = hashlib.sha256(pd.util.hash_pandas_object(nodes,index=True).values.tobytes()).hexdigest()
    from importlib.metadata import version
    provenance = {"data_sha256":fingerprint,"ordered_gene_ids":nodes.gene_id.tolist(),
                  "evaluation_mode":"inductive_group_holdout" if spec.group_column else "inductive_random_holdout",
                  "final_training_ids":nodes.gene_id.iloc[fit].tolist(),
                  "final_calibration_ids":nodes.gene_id.iloc[calibration].tolist(),
                  "executed_method":spec.method,"classes":classes,
                  "final_estimator":type(final.model).__name__,
                  "versions":{name:version(name) for name in ('numpy','pandas','scikit-learn','scipy','umap-learn')},
                  "limitations":["Calibration may shift for poorly studied, unlabelled genes.",
                                 "Residual intervals have no guaranteed coverage under family/domain shift.",
                                 "Development cross-validation is not independent prospective validation."]}
    return PredictionResult(spec,table,metrics,per_class,reports,provenance,final)


def run_multilabel(nodes, targets, **kwargs):
    """Predict each measured trait independently; missing entries are never negative.

    Targets must be separate observed binary columns, with NaN for unassayed genes.
    Every target and its registered family is excluded from the other target runs
    to avoid treating an unmeasured sibling phenotype as available input.
    """
    exclude = tuple(set(targets) | set(kwargs.pop("exclude", ())))
    result = {}
    for target in targets:
        values = set(pd.to_numeric(nodes[target], errors="coerce").dropna().unique())
        if not values <= {0,1}:
            raise ValueError(f"{target} is not an observed binary trait")
        result[target] = run(nodes, TaskSpec(target, exclude=exclude, **kwargs))
    return result
