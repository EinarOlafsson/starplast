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
    from .holdout_cv import _known
    known = _known(values)
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
    if len(classes) < 2:
        return np.full(len(X), -1), None, []
    codes = {c: i for i, c in enumerate(classes)}
    yi = np.array([codes[v] for v in y])
    partition = np.full(len(X), -1)
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    index = np.arange(len(X))[mask]
    for train, test in splitter.split(X[mask], yi):
        model.fit(X[mask][train], yi[train])
        partition[index[test]] = model.predict(X[mask][test])
    model.fit(X[mask], yi)
    from .holdout_cv import _known
    rest = ~_known(truth)
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
    if not (hasattr(clf, "coef_") or hasattr(clf, "estimators_")):
        # Not every method has coefficients. A propagator's evidence is the graph, and there is no
        # per-column weight to report -- an empty frame says that, where an exception would claim the
        # method was broken.
        return pd.DataFrame()
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


# --------------------------------------------------------------------------- propagation (47.2)
#: How much of each step returns to the seed. 0.5 keeps a walk local enough that a label stays near
#: the genes that carry it; at 0.1 the scores approach the stationary distribution of the graph,
#: which is a statement about node degree rather than about the label.
RESTART = 0.5

#: Power iterations. The walk converges geometrically at this restart, and 30 puts the residual far
#: below the differences the argmax is decided on.
ITERATIONS = 30


def layer_matrix(layer: str, n: int, graph=None):
    """One edge layer as a symmetric, degree-normalised sparse operator.

    Symmetric normalisation (D^-1/2 W D^-1/2) rather than row-stochastic, because a random walk on
    the row-stochastic matrix concentrates on hubs -- and the hubs of the two co-mention layers are
    the genes people write about most. Attention bias is a correctness problem in this project, not a
    cosmetic one, and the normalisation is the first place it gets in.
    """
    import scipy.sparse as sp
    from . import paths
    if graph is None:
        graph = np.load(paths.cache_file("graph.npz"), allow_pickle=True)
    a, b = graph[f"{layer}__a"], graph[f"{layer}__b"]
    # `.files` on an npz, plain membership on a mapping: a test builds one of these by hand, and a
    # loader that only accepts the on-disk shape forces every test through a temporary file.
    keys = getattr(graph, "files", graph)
    w = graph[f"{layer}__w"] if f"{layer}__w" in keys else np.ones(len(a))
    W = sp.coo_matrix((np.asarray(w, dtype=float), (a, b)), shape=(n, n)).tocsr()
    W = W + W.T                                     # the layers are stored once per undirected edge
    degree = np.asarray(W.sum(axis=1)).ravel()
    inv = np.divide(1.0, np.sqrt(degree), out=np.zeros_like(degree), where=degree > 0)
    D = sp.diags(inv)
    return D @ W @ D


def graph_size():
    """How many nodes the shipped graph was built over, or None when there is no graph here.

    Edge endpoints are stored as positions in the node table, which makes them meaningless against
    any other table. This is what lets a caller be told that rather than shown a confident answer
    about the wrong genes.
    """
    import os
    from . import paths
    path = paths.cache_file("graph.npz")
    if not os.path.exists(path):
        return None
    return int(np.load(path, allow_pickle=True)["xyz"].shape[0])


class Propagator:
    """Random walk with restart over one edge layer, wearing a classifier's interface.

    Deliberately shaped like an estimator so it can go through `out_of_fold` unchanged: the honesty
    problem is identical to the classifier's and the solution has to be too. A gene's own seed makes
    its own score enormous, so scoring a gene the walk was seeded from would report the seeding, not
    the biology. Every gene is therefore scored by a walk seeded only from the OTHER folds.

    `X` here is a column of node positions rather than features -- the walk's evidence is the graph,
    not the table -- which is what lets the same cross-validation machinery drive both methods.
    """

    def __init__(self, matrix, restart: float = RESTART, iterations: int = ITERATIONS):
        self.matrix, self.restart, self.iterations = matrix, restart, iterations
        self.scores_ = None

    def fit(self, X, y):
        """Diffuse each class from the genes that carry it, and keep the resulting fields."""
        positions = np.asarray(X).ravel().astype(int)
        classes = sorted(set(int(v) for v in y))
        fields = []
        for c in classes:
            seed = np.zeros(self.matrix.shape[0])
            members = positions[np.asarray(y) == c]
            if len(members):
                seed[members] = 1.0 / len(members)
            p = seed.copy()
            for _ in range(self.iterations):
                p = (1 - self.restart) * (self.matrix @ p) + self.restart * seed
            fields.append(p)
        self.scores_, self.classes_ = np.vstack(fields), np.array(classes)
        return self

    def predict(self, X):
        """The class whose field is highest at each node.

        Nodes with no positive propagated support abstain with code -1. A class
        tie at zero is absence of evidence, not support for the first class.
        """
        positions = np.asarray(X).ravel().astype(int)
        scores = self.scores_[:, positions]
        best = np.argmax(scores, axis=0)
        return np.where(scores[best, np.arange(len(positions))] > 0,
                        self.classes_[best], -1)



def propagation(layer: str, index: np.ndarray, truth: pd.Series, n_nodes: int, seed: int = 42,
                folds: int = FOLDS, graph=None) -> dict:
    """Answer by diffusing the label across one edge layer.

    The reason this is worth having is the constraint it removes. A clustering needs fifteen labelled
    genes inside one cluster before it can say anything, which is why several questions in the
    shipped catalogue return nothing and why a 221-gene control cannot corroborate anything.
    Propagation has no such floor: seeded with the genes that ARE labelled, it scores all 8,140 and
    ranks them.

    One layer, named. The thirteen layers mean thirteen different things and merging them is the
    mistake this project's second design decision exists to prevent -- and two of them are
    attention-biased, so a walk over the merged graph would carry the literature's popularity contest
    into every answer.
    """
    matrix = layer_matrix(layer, n_nodes, graph=graph)
    X = np.asarray(index).reshape(-1, 1)
    partition, fitted, classes = out_of_fold(X, truth, Propagator(matrix), seed=seed, folds=folds)
    return {"partition": partition, "classes": classes, "model": fitted,
            "settings": {"method": f"propagation:{layer}", "restart": RESTART,
                         "iterations": ITERATIONS, "folds": folds}}


# --------------------------------------------------------------------------- boosting (47.4)
def raw_matrix(nodes: pd.DataFrame, columns) -> tuple:
    """The chosen columns as numbers, with missing values LEFT MISSING.

    Every other method in this project goes through `build_matrix`, whose four missing-value policies
    all resolve absence one way or another: impute the median, drop the column, drop the gene, or add
    an indicator. Each is a decision made before the model sees the data, and one of them -- dropping
    genes -- changes WHICH GENES the labels describe, a bug class this project has already been
    bitten by.

    A histogram-boosted tree needs none of that: it learns a split for the missing branch directly,
    so "not measured" becomes evidence rather than a value invented to stand in for it. That is the
    reason this method is here and not a marginal accuracy gain, so the matrix is built raw.
    """
    keep = [c for c in columns if c in nodes.columns
            and pd.api.types.is_numeric_dtype(nodes[c])]
    if not keep:
        return np.zeros((len(nodes), 0)), []
    return nodes[keep].to_numpy(dtype=float), keep


class _DropConstant:
    """Drop columns with fewer than two distinct finite values, measured on the data it is fit to.

    A column with one value cannot produce a split, so this changes no answer -- and since
    scikit-learn 1.9 it cannot even be BINNED: the histogram binner takes a sliding window of two
    over a column's distinct values and raises `window shape cannot be larger than input array
    shape`, which names neither the column nor the cause. NaN does not count as a value, because the
    boosted model treats missingness as a direction at a split rather than as a level.
    """

    def fit(self, X, y=None):
        """Record which columns vary in this matrix."""
        self.keep_ = [j for j in range(X.shape[1])
                      if np.unique(X[np.isfinite(X[:, j]), j]).size > 1]
        self.dropped_ = [j for j in range(X.shape[1]) if j not in self.keep_]
        return self

    def transform(self, X):
        """Take those columns, in order."""
        return X[:, self.keep_]

    def fit_transform(self, X, y=None):
        """Fit and transform in one pass, which is what a pipeline calls."""
        return self.fit(X, y).transform(X)

    def get_params(self, deep=True):
        """No parameters, but a pipeline clones its steps and clone() asks."""
        return {}

    def set_params(self, **_):
        """Likewise, and there is nothing to set."""
        return self


def boosted(X: np.ndarray, truth: pd.Series, seed: int = 42, folds: int = FOLDS,
            max_iter: int = 200) -> dict:
    """Histogram gradient boosting, scored out of fold like every other supervised method here.

    `HistGradientBoostingClassifier` rather than a new dependency: it is already in scikit-learn,
    which this project already requires, and it takes NaN natively. Adding xgboost or lightgbm for
    the same capability would put two gigabytes of wheels behind a method whose whole argument is
    that it needs less preprocessing, not more.
    """
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.pipeline import make_pipeline
    # The constant-column guard is INSIDE the estimator, which is the only place that can be right:
    # a column can be constant within one cross-validation fold and varying in the next, so a filter
    # applied once to the whole matrix does not protect the fold that actually breaks. A pipeline
    # step remembers the columns it kept, so `predict` and permutation importance both still take
    # the full matrix and every score stays aligned with its own column name.
    model = make_pipeline(_DropConstant(),
                          HistGradientBoostingClassifier(max_iter=max_iter, random_state=seed,
                                                         early_stopping=False))
    partition, fitted, classes = out_of_fold(X, truth, model, seed=seed, folds=folds)
    dropped = len(getattr(fitted[0], "dropped_", [])) if fitted is not None else 0
    return {"partition": partition, "classes": classes, "model": fitted,
            "settings": {"method": "boosted", "max_iter": max_iter, "folds": folds,
                         "missing": "native", "constant_features_dropped": int(dropped)}}


def importances(fitted, feature_names, X: np.ndarray = None, truth: pd.Series = None,
                seed: int = 42) -> pd.DataFrame:
    """Which measurements the tree relied on, by permutation rather than by split count.

    Split counts are the cheap answer and a misleading one: a column with many distinct values offers
    more places to split and accumulates a high count without predicting anything. Permutation
    importance asks the question that matters -- how much worse does the model get when this column
    is shuffled -- and is measured on the fitted model rather than read off its structure.
    """
    if fitted is None or X is None or truth is None or not len(feature_names):
        return pd.DataFrame()
    from sklearn.inspection import permutation_importance
    mask, values = _usable(truth)
    if mask.sum() < 20:
        return pd.DataFrame()
    classes = sorted(set(values[mask]))
    codes = {c: i for i, c in enumerate(classes)}
    y = np.array([codes[v] for v in values[mask]])
    out = permutation_importance(fitted, X[mask], y, n_repeats=3, random_state=seed, n_jobs=1)
    table = pd.DataFrame({"feature": list(feature_names), "importance": out.importances_mean,
                          "sd": out.importances_std})
    return table.sort_values("importance", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- multiplex (47.3)
def layer_communities(layer: str, n_nodes: int, resolution: float = 1.0, seed: int = 42,
                      graph=None) -> np.ndarray:
    """Communities of ONE layer, as a partition over all nodes. Unreached nodes get -1.

    Greedy modularity via networkx rather than Leiden: leidenalg and igraph are not installed here
    and neither is worth two more dependencies for a method whose purpose is to be COMPARED against
    the clustering rather than to win. Resolution is exposed because modularity has a resolution
    limit and a single value would silently decide how big a community is allowed to be.
    """
    import networkx as nx
    import scipy.sparse as sp
    matrix = sp.triu(layer_matrix(layer, n_nodes, graph=graph), k=1).tocoo()
    G = nx.Graph()
    G.add_nodes_from(range(n_nodes))
    G.add_weighted_edges_from(zip(matrix.row.tolist(), matrix.col.tolist(), matrix.data.tolist()))
    labels = np.full(n_nodes, -1)
    groups = nx.community.greedy_modularity_communities(G, weight="weight",
                                                        resolution=resolution)
    for i, members in enumerate(groups):
        if len(members) < 2:
            continue                    # a community of one is a node with no community
        labels[list(members)] = i
    # An isolated node joins no community, and saying so is the point: it carries no relational
    # evidence, and the enrichment gate downstream must not treat a bag of isolates as a group.
    degree = np.asarray(layer_matrix(layer, n_nodes, graph=graph).sum(axis=1)).ravel()
    labels[degree == 0] = -1
    return labels


def multiplex_communities(layers, n_nodes: int, resolution: float = 1.0, seed: int = 42,
                          graph=None, min_agreement: float = 0.5) -> dict:
    """Communities that several layers AGREE on, built as a consensus rather than a merged graph.

    Design decision 2 of this project says the thirteen edge types are not one graph and are never
    merged silently. Summing adjacency matrices would do exactly that, and would let the two
    attention-biased layers -- which reproduce the literature's popularity contest -- pull every
    community toward the well-published genes.

    So each layer is clustered ON ITS OWN, and the consensus is over the resulting PARTITIONS: two
    genes are joined when at least `min_agreement` of the layers that can see them both put them
    together. A merge is then a statement several independent measurements make, and the layers that
    cannot see a pair abstain rather than voting no -- which matters here, because layer sizes span
    four orders of magnitude and `ip_ms` has 64 edges against `compartment`'s 118,712.
    """
    import scipy.sparse as sp
    partitions = [layer_communities(L, n_nodes, resolution=resolution, seed=seed, graph=graph)
                  for L in layers]
    if not partitions:
        return {"partition": np.full(n_nodes, -1), "classes": [], "model": None,
                "unsupervised": True, "settings": {"method": "multiplex", "layers": []}}
    together = sp.csr_matrix((n_nodes, n_nodes), dtype=float)
    seen = sp.csr_matrix((n_nodes, n_nodes), dtype=float)
    for labels in partitions:
        for community in np.unique(labels[labels >= 0]):
            members = np.flatnonzero(labels == community)
            if len(members) < 2 or len(members) > n_nodes // 2:
                continue        # a community holding half the graph is a bisection, not a community
            block = sp.csr_matrix((np.ones(len(members)), (members, np.zeros(len(members), int))),
                                  shape=(n_nodes, 1))
            together = together + block @ block.T
        visible = np.flatnonzero(labels >= 0)
        block = sp.csr_matrix((np.ones(len(visible)), (visible, np.zeros(len(visible), int))),
                              shape=(n_nodes, 1))
        seen = seen + block @ block.T
    agreement = together.multiply(seen.power(-1.0))
    agreement.data[np.isnan(agreement.data)] = 0.0
    consensus = (agreement >= min_agreement).astype(float)
    consensus.setdiag(0)
    consensus.eliminate_zeros()
    n_groups, labels = sp.csgraph.connected_components(consensus, directed=False)
    sizes = pd.Series(labels).value_counts()
    singletons = set(sizes[sizes < 2].index)
    out = np.array([-1 if v in singletons else v for v in labels])
    return {"partition": out, "classes": [], "model": None, "unsupervised": True,
            "settings": {"method": "multiplex", "layers": list(layers), "resolution": resolution,
                         "min_agreement": min_agreement, "groups": int(n_groups)}}


def classification_recovery(codes, truth, classes):
    """Score fixed class identities on held-out predictions; abstentions are errors.

    Returns the legacy summary/table shape for recipe consumers, with an explicit
    score_kind. Class labels are never remapped by looking at evaluation answers.
    """
    from .holdout_cv import ABSTAIN, _known, score_predictions
    codes = np.asarray(codes, dtype=int)
    values = pd.Series(truth, dtype="object")
    eligible = _known(values) & values.astype(str).isin(classes).to_numpy()
    decoded = np.asarray([classes[c] if 0 <= c < len(classes) else ABSTAIN for c in codes])
    summary, per = score_predictions(values[eligible], decoded[eligible], classes)
    per = per.rename(columns={"category": "label", "support": "n_label"})
    per["cluster"] = [classes.index(label) for label in per.label]
    per["n_in_cluster"] = [int(((decoded == label) & eligible & (values.astype(str) == label)).sum())
                           for label in per.label]
    total = per.n_label.sum()
    summary.update(mean_f1=float((per.f1 * per.n_label).sum() / total) if total else 0.0,
                   best_f1=float(per.f1.max()) if len(per) else 0.0,
                   best_label=str(per.loc[per.f1.idxmax(), "label"]) if len(per) else "",
                   n_labels_scored=len(per), n_labels_recovered=int((per.f1 >= .5).sum()),
                   score_kind="out_of_fold_fixed_class_prediction")
    return summary, per


def classification_candidates(codes, truth, gene_ids, classes, validation,
                              min_precision=.3, min_enrichment=2., min_support=15):
    """Attach class-level CV evidence to unknown genes without changing model calls.

    CV class precision describes performance among labelled validation genes; it
    is not a calibrated probability for an individual, potentially shifted gene.
    """
    from .holdout_cv import _known
    columns = ["gene_id", "predicted", "cluster", "cv_class_precision", "cv_class_recall",
               "enrichment", "base_rate", "n_labelled_in_cluster", "evidence_status"]
    known = _known(truth)
    values = pd.Series(truth, dtype="object").astype(str).to_numpy()
    eligible = known & np.isin(values, classes)
    rows = []
    for code, category in enumerate(classes):
        scores = validation[validation.label == category]
        if scores.empty:
            continue
        row = scores.iloc[0]
        prior = float((values[eligible] == category).mean()) if eligible.any() else 0.
        lift = float(row.precision / prior) if prior else 0.
        if row.n_label < min_support or row.precision < min_precision or lift < min_enrichment:
            continue
        for gene in np.asarray(gene_ids)[~known & (np.asarray(codes) == code)]:
            rows.append(dict(gene_id=str(gene), predicted=category, cluster=code,
                             cv_class_precision=float(row.precision), cv_class_recall=float(row.recall),
                             enrichment=lift, base_rate=prior, n_labelled_in_cluster=int(row.n_label),
                             evidence_status="model_prediction_with_class_level_cv_support"))
    return pd.DataFrame(rows, columns=columns)
