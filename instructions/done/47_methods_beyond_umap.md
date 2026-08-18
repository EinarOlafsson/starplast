# 47 — Inference by network propagation, graph learning and regression, compared honestly

**Status: open. Requested by the user 2026-08-18.**

## Why

Every claim this program makes is currently produced one way: build a UMAP of a feature matrix,
cluster it, and check whether a held-out label comes back. That is a deliberately conservative
design and it should stay. But it flattens away most of what is on this disk — **13 edge layers and
302,756 edges** contribute nothing to a single prediction, and the pipeline has no baseline, so
"mean F1 0.20" cannot presently be called good or bad by anyone.

Four families of method would each add something the clustering cannot, and the point of doing them
together is not variety. It is that a recipe already names its inputs, its holdout and its control,
so every method can be run on IDENTICAL data behind IDENTICAL guards and judged by the same
independent control. That comparison is the deliverable.

## 0. BLOCKING PREREQUISITE: leakage closure must understand edges

`excluded_detail` reasons entirely about node-table COLUMNS. The graph is a separate leak surface and
the guard cannot currently see it:

**`compartment` is both a 118,712-edge layer and the most-used holdout in this project.** A
propagation or graph-learning method fed the compartment edge layer while holding out the compartment
label would recover it perfectly, and every existing guard would pass.

The catalogue already models edges (`slots.edge_types`, the `edge:<type>` pattern), so closure has to
return banned EDGE LAYERS beside banned columns, by the same family and class scoping. This lands
FIRST. Supervised and propagation methods are far better at finding leaks than UMAP is: UMAP blurs a
0.7-associated column, a gradient-boosted tree exploits it.

## 1. Regularised linear models — first, because it is the missing yardstick

L1 logistic regression per label on the closed feature set, same holdout, same seed, same folds.

* it yields a COEFFICIENT PER COLUMN, so an answer reads "this label is predicted by these six
  measurements" — something a cluster can never say;
* no cluster-size floor, so labels too rare to cluster become answerable;
* **it is the baseline this project does not have.** If plain logistic regression recovers
  `compartment` better than UMAP + HDBSCAN, then the map is a visualisation rather than an inference
  engine, and everything below should be judged against that number rather than against nothing.

## 2. Network propagation — the cheapest real win

Random-walk-with-restart and label propagation, run PER LAYER, never on the merged graph.

This removes the constraint that has shaped every result so far: a cluster needs 15 labelled genes
before it can support a claim, which is why the immunity axis ships nothing and why
`iedb_epitope_count` (221 genes) cannot serve as a control. Propagation has no such floor — seeded
with 34 known genes it returns a ranked list over all 8,140.

Per layer and not merged, because the layers mean different things and two of them are
attention-biased. `comention` and `comention_ft` reproduce the literature's popularity contest, so
centrality computed on them measures publication history. The layers that are NOT attention-biased —
`xlms` (measured physical proximity), `ip_ms`, `struct`, `cofitness` — are where a network method
could find something the literature has not.

## 3. Multiplex community detection

Leiden with per-layer weights, or consensus across layers, as a direct alternative to UMAP +
HDBSCAN under the same scoring. It respects decision 2 — the 13 edge types are not one graph and are
never merged silently — which the current feature matrix violates by construction.

## 4. Gradient boosting, then graph learning

Gradient boosting handles missing values natively, which matters more here than usual: the
`na_policy` currently drops or imputes rows, and that changes WHICH GENES the labels describe — a bug
class this project has already been bitten by. SHAP then gives per-gene attribution, which is what
makes a prediction actionable at a bench.

A relational GNN (one relation per edge type) is the most powerful option and the LAST to build:
8,140 nodes is small, it will overfit, and it is the least interpretable thing here. Build it only if
the baseline says the simpler methods are leaving something on the table.

Separately, **link prediction on `xlms`/`ip_ms` is a different question this data supports**:
predicting missing physical interactions rather than gene labels.

## The shape of the deliverable

`Recipe` gains `method`, and every method inherits the closure, the seed, the validation holdout and
the corroboration test unchanged. One question, five methods, one control judging all of them.

## Done when

* Closure bans edge layers as well as columns, with a test that a graph method cannot eat the layer
  its holdout came from.
* A recipe can be run by any implemented method without changing anything else about it.
* The linear baseline is reported beside the clustering for the same question, on the same folds.
* Propagation answers at least one question the clustering could not, BECAUSE of the cluster floor.
* Every method carries a shuffled-label negative control, as the clustering already does.
* A comparison table exists: per question, per method — genes named, precision, and the fraction the
  independent control corroborates.
