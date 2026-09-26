"""The techniques the strategies are built from, each explained once.

A strategy's name ends with its method in brackets; its `techniques` lists the pieces that method is
made of, and every piece is defined here: what it does, and why a strategy would use it. The
Strategies tab shows these beside each strategy, and :func:`glossary` returns them as a table, so a
reader meets "HDBSCAN" or "random walk with restart" with its explanation attached.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Technique:
    """One technique: its name, what kind of step it is, what it does, and why a strategy uses it."""
    key: str
    name: str
    kind: str
    what: str
    why: str


_T = []


def _t(key, name, kind, what, why):
    _T.append(Technique(key, name, kind, what, why))


EMBED, CLUSTER, NEIGHBOUR, NETWORK, MODEL, STATS = (
    "embedding", "clustering", "neighbours", "networks", "models", "statistics")

_t("umap", "UMAP", EMBED,
   "Uniform Manifold Approximation and Projection: places every gene in three dimensions so that "
   "genes with similar measurements stay close, preserving local neighbourhoods over global "
   "distances. n_neighbors sets how local, min_dist how tightly points may pack.",
   "Makes the combined measurements into a map whose structure can be clustered and looked at; "
   "the map is only ever built from columns the held-out label's closure permits.")
_t("hdbscan", "HDBSCAN", CLUSTER,
   "Hierarchical density-based clustering: finds dense groups of any shape and leaves sparse genes "
   "unclustered (noise) instead of forcing them into a group. min_cluster_size is the smallest "
   "group it reports; 'eom' keeps the most persistent clusters, 'leaf' the finest.",
   "Clusters a map without being told how many clusters there are, and admits when a gene belongs "
   "nowhere.")
_t("settings_walk", "Settings walk", CLUSTER,
   "Builds and clusters the map under every combination of settings on a grid (and optionally "
   "each family of measurements alone), and keeps the combination that scores best on labels "
   "visible to it.",
   "No single map setting is right for every question; the walk searches them, and the self-test "
   "then scores the chosen one on labels the choice never saw.")
_t("co_association", "Co-association consensus", CLUSTER,
   "Counts, for every pair of genes, the share of clusterings that put them together, then joins "
   "genes by average linkage on one minus that share.",
   "Keeps only structure that survives the arbitrary choices a map requires.")
_t("map_quality", "Label-free map quality", CLUSTER,
   "Scores a clustering by how much of the proteome it places and how evenly, without looking at "
   "any label.",
   "Lets a map be tuned without the labels it will later be tested against.")
_t("knn", "k-nearest neighbours (kNN)", NEIGHBOUR,
   "Finds the k genes most similar to a gene (Euclidean distance on standardised measurements, or "
   "on map coordinates) and lets their labels vote, each weighted by 1 / distance.",
   "The plainest guilt by association: no model, no clustering, and every call can be traced to "
   "the genes that made it.")
_t("knn_graph", "kNN graph", NEIGHBOUR,
   "Links every gene to its k nearest neighbours in measurement space, symmetrised and "
   "degree-normalised, so a table can be walked like a network.",
   "Lets measurement similarity and measured networks be combined in one graph.")
_t("network_vote", "Network neighbour vote", NETWORK,
   "Calls a gene by the labels of its direct partners in one or more edge layers, each partner "
   "weighted by the edge's weight (contact counts, correlation strength). Abstains where a gene "
   "has no labelled partner.",
   "Uses measured relationships -- crosslinks, pulldowns, co-expression -- directly as evidence.")
_t("random_walk_restart", "Random walk with restart", NETWORK,
   "Diffuses from seed genes along the network: at each step a walker moves to a neighbour or, "
   "with the restart probability, jumps back to a seed. The stationary visiting frequency scores "
   "every gene.",
   "Reaches genes several steps from the seeds while still favouring close ones; restart sets how "
   "far the influence spreads.")
_t("greedy_modularity", "Greedy modularity communities", NETWORK,
   "Merges groups of genes greedily while the network's modularity -- edges inside groups beyond "
   "what degree alone predicts -- increases. The resolution parameter sets the community size.",
   "Finds communities in each layer separately, before any layer is allowed to dominate.")
_t("louvain", "Louvain consensus", NETWORK,
   "Joins genes that a sufficient share of the layers put together, and splits the joined graph "
   "into communities by Louvain modularity optimisation.",
   "Keeps groupings several independent kinds of evidence agree on without chaining everything "
   "into one component.")
_t("tm_score", "Structural similarity (TM-score)", NETWORK,
   "Foldseek TM-score between predicted structures: 0 to 1, above about 0.5 usually the same fold, "
   "regardless of sequence identity.",
   "Fold is conserved long after sequence is not, so it can annotate proteins no sequence search "
   "reaches.")
_t("chance_weighting", "Chance-weighted voting", NETWORK,
   "Scores each evidence source on an inner holdout of the known labels and weights its vote by "
   "its accuracy minus its own chance level; a source that only guesses the commonest class gets "
   "no vote.",
   "Lets the evidence that has earned trust for this label count most, and says which that is.")
_t("agreement", "Agreement of independent callers", NETWORK,
   "Calls a gene only where at least a set number of callers built on different evidence give the "
   "same label.",
   "Trades reach for precision: independent errors rarely agree.")
_t("logistic_regression", "Logistic regression", MODEL,
   "A linear model of the log-odds of a class (or of an edge) from standardised features, with "
   "classes balanced and L2 regularisation (C: smaller is stronger).",
   "Interpretable -- each feature has a weight -- and hard to overfit; the baseline any fancier "
   "model must beat.")
_t("gradient_boosting", "Gradient boosting", MODEL,
   "An ensemble of shallow decision trees, each fitted to the errors of the ones before "
   "(scikit-learn's histogram gradient boosting). Handles missing values and non-linear effects.",
   "Captures interactions and thresholds a linear model cannot, at the cost of interpretability.")
_t("ridge", "Ridge regression", MODEL,
   "Linear regression with an L2 penalty that shrinks every coefficient towards zero.",
   "A stable linear alternative to boosting, with missing values imputed first.")
_t("pu_bagging", "Positive-unlabelled bagging", MODEL,
   "Trains many classifiers, each on the positives against a small random draw of other genes, "
   "and scores every gene only by the models that did not train on it (out of bag).",
   "When all you have are positives, treating every other gene as a negative teaches the model to "
   "reject the hidden positives you want to find; small random draws rarely contain them.")
_t("soft_impute", "Soft-impute (low-rank SVD)", MODEL,
   "Completes a matrix with missing entries by repeatedly filling them from a truncated singular "
   "value decomposition of the current estimate, at a chosen rank.",
   "Borrows strength across all measurements at once; columns are checked on hidden entries "
   "before any imputed value is trusted.")
_t("residualisation", "Residualisation on a baseline", MODEL,
   "Regresses the condition's measurement on its baseline (both rank-scaled) and keeps the "
   "residual: the part the baseline does not explain.",
   "Separates what matters only in one condition from what matters everywhere.")
_t("spectral_embedding", "Spectral embedding", MODEL,
   "Node vectors from a truncated SVD of a layer's adjacency matrix; a pair is described by the "
   "elementwise product of its two vectors.",
   "Captures network position beyond direct evidence; kept only where it beats the interpretable "
   "baseline on a degree-matched null.")
_t("triadic_closure", "Shared partners (triadic closure)", NETWORK,
   "Counts the partners two genes share in a layer; two genes linked to the same partners are "
   "often linked themselves.",
   "The classic link-prediction feature, used alongside support from other layers.")
_t("degree_matched", "Degree-matched non-pairs", STATS,
   "Negative pairs drawn so each end has a degree similar to the positives' ends.",
   "Against random non-pairs, a model scores highly just by knowing which genes are well "
   "connected; matched non-pairs leave only what the evidence says about the pair.")
_t("attention_residual", "Attention-corrected co-mention", STATS,
   "Replaces each pair's co-mention count by its residual over what the two genes' own publication "
   "counts predict (log2 observed / expected).",
   "Stops famous genes being linked just because they are famous.")
_t("support_count", "Multi-layer support count", NETWORK,
   "Counts, for each gene pair, how many independent measurement layers link it.",
   "Independent kinds of evidence agreeing on a pair is stronger than any one of them.")
_t("profile_correlation", "Profile correlation", STATS,
   "One minus the correlation of two genes' profiles across every permitted measurement both "
   "have.",
   "Measures how differently two paralogs behave, whatever they are called.")
_t("orthogroup_mapping", "Orthogroup mapping", STATS,
   "Aggregates a measurement or label over each orthogroup in the other species (median, or the "
   "commonest label) and learns its relation to this species' values on genes measured in both.",
   "Carries evidence across species without merging the two tables.")
_t("hypergeometric", "Hypergeometric enrichment", STATS,
   "The probability of drawing at least the observed number of a label's genes into a group of "
   "that size by chance (one-sided Fisher's exact test).",
   "The standard test that a cluster or gene list holds more of a label than chance.")
_t("rank_sum", "Rank-sum test", STATS,
   "Mann-Whitney U: whether a measurement's values are higher (or lower) in one group of genes "
   "than in the rest, using ranks only.",
   "Tests a measurement against a gene list without assuming a distribution.")
_t("chi_square", "Chi-square test (Cramer's V)", STATS,
   "Tests whether a categorical label is distributed differently across clusters; Cramer's V is "
   "the effect size, 0 to 1.",
   "Asks what a map built from one kind of measurement encodes about another kind.")
_t("kruskal_wallis", "Kruskal-Wallis test (eta-squared)", STATS,
   "Rank-based test that a measurement differs between clusters; eta-squared is the effect size.",
   "The numeric counterpart of the chi-square test above.")
_t("bh_fdr", "Benjamini-Hochberg FDR", STATS,
   "Adjusts many p-values together so the expected share of false discoveries among the "
   "significant ones stays below the chosen rate (q).",
   "Every strategy that tests many clusters, labels or measurements at once corrects for it.")
_t("grouped_cv", "Orthogroup-grouped cross-validation", STATS,
   "Splits genes into folds by whole orthogroups, so a gene is never predicted by a model that saw "
   "it or its paralog.",
   "Paralogs share measurements; without grouping, a model can look accurate by recognising "
   "family members.")
_t("ablation", "Evidence ablation", STATS,
   "Scores each kind of measurement alone and the combination without it, out of fold.",
   "Says which experiments carry which biology and which could be dropped without loss.")
_t("bimodal_split", "Bimodal split", STATS,
   "Divides a cluster's values on a second measurement into two modes and requires them to be at "
   "least a set distance apart.",
   "Finds categories whose members fall into two kinds on another measurement.")

_t("split_conformal", "Split conformal prediction", MODEL,
   "Calibrates any model on genes it did not train on: the errors (or the scores of the true "
   "labels) on those genes set a threshold so that the returned sets or intervals contain the truth "
   "for at least 1 - alpha of new genes.",
   "Turns scores into a stated error rate that holds whatever the model, and says which genes the "
   "data cannot decide.")
_t("sgc", "Simplified graph convolution", NETWORK,
   "Multiplies the measurement matrix by the degree-normalised network (with self-loops) once per "
   "step, so each gene gains features that average its neighbours' measurements.",
   "Lets a linear model learn from a gene's network neighbourhood -- the core of a graph neural "
   "network -- while staying transparent and leak-free.")
_t("random_forest", "Random forest", MODEL,
   "Hundreds of decision trees, each grown on a bootstrap sample of genes with a random subset of "
   "measurements at every split, voting on the class; classes are re-weighted to balance.",
   "Finds thresholds and interactions a linear model cannot, and is robust to scale and outliers.")
_t("permutation_importance", "Permutation importance", STATS,
   "Shuffles one measurement among held-out genes and records how much balanced accuracy drops; "
   "repeated to give a mean and a spread.",
   "Measures what a model actually relies on for new genes, unlike impurity importance, which "
   "rewards measurements with many distinct values.")
_t("stacking", "Stacking", MODEL,
   "Trains a meta-model on the out-of-fold predictions of several base models, so the combination "
   "is learned on genes no base model trained on.",
   "Learns per label and per class how far to trust each kind of evidence, without rewarding a "
   "base model for memorising its training genes.")

TECHNIQUES: dict = {t.key: t for t in _T}


def glossary(keys=None) -> pd.DataFrame:
    """The techniques (all, or `keys` in that order) as a table: name, kind, what and why."""
    chosen = [TECHNIQUES[k] for k in keys] if keys is not None else _T
    return pd.DataFrame([{"key": t.key, "technique": t.name, "kind": t.kind, "what": t.what,
                          "why": t.why} for t in chosen])
