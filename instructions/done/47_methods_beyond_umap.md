# 47 — Inference by network propagation, graph learning and regression, compared honestly

**Status: DONE 2026-08-18**, except the GNN, which is deliberately not built — see the end.

## The answer, measured across all twenty shipped questions

Every method ran the same recipes: same leakage closure, same seed, same holdout, same validation
control, same statistics. 84 runs, in `instructions/done/47_method_comparison.csv`.

| method | questions | mean F1 | genes named | corroborated | rate |
|---|---|---|---|---|---|
| propagation · struct | 20 | **0.530** | 975 | 645 | 66% |
| propagation · xlms | 20 | 0.504 | 630 | 189 | 30% |
| logistic (L1) | 20 | 0.501 | 4,105 | 2,752 | **67%** |
| umap + hdbscan | 20 | **0.126** | 942 | 433 | 46% |

**Every alternative recovers the holdout about four times better than the map does**, across twenty
questions rather than one lucky case. On the first question specifically, gradient boosting reaches
mean F1 0.552 and best-label 0.787 — the strongest single result in the project.

The plain reading: **the UMAP is a visualisation, and it is not this project's inference engine.**
That is worth knowing before more is built on top of it, and it is a question nobody could ask before
because there was nothing to compare against.

Two results worth more than the ranking:

* **A 2,842-edge layer beats the entire 361-column feature matrix.** `xlms` is measured crosslink
  proximity — not attention-biased, not predicted, not derived — and diffusing a label across it
  recovers that label better than embedding every column in the table and clustering the result.
* **`struct` names 645 corroborated genes of 975.** Foldseek structural similarity needs no
  orthology, so it reaches lineage-specific effectors that every homology-based route misses.

## What each method contributes that the clustering cannot

* **Linear (L1 logistic)** — a coefficient per column, so an answer reads "this label is predicted by
  these six measurements" rather than "these genes go together". It is also the baseline: without it
  "mean F1 0.20" could not be called good or bad.
* **Propagation, per layer** — no cluster-size floor. A clustering needs fifteen labelled genes
  INSIDE one cluster; propagation seeded with the genes that are labelled scores all 8,140 and ranks
  them, which is what makes sparse labels answerable at all.
* **Boosting** — missing values stay missing. Every other path resolves absence before the model sees
  it, and one of those policies changes WHICH GENES the labels describe. A tree learns a split for
  the missing branch, so "not measured" becomes evidence. A test pins this directly: a column whose
  values carry nothing but whose ABSENCE tracks the label is learned at >0.9, and every imputing
  policy destroys that signal by construction.
* **Multiplex consensus** — communities each layer finds ON ITS OWN, joined where at least half the
  layers that can see a pair agree. Never a summed adjacency: that would merge the thirteen types
  design decision 2 forbids merging, and would let the two attention-biased layers pull every
  community toward the well-published genes. Measured: 267 components, 6 usable groups, mean F1
  0.287, and it names nobody — the consensus is conservative, and that is reported rather than tuned
  away.

## The guards that had to come first

**Closure now bans edge layers**, and the family is closed at ANY scope. Class scope — what every
recipe uses — banned nothing for `compartment`, because `Tg_shared_compartment` sits under
`molecular relationships` while the label sits under `cell organization`: the class of the target
does not contain the layer built from the target. For columns this never arises, since association,
provenance and same-quantity catch the family by other routes; an edge layer has none of those.

**Enforced, not merely reported**: walking `compartment` while holding it out is refused, a bare
`propagation` is refused with the list of layers, and a walk over a SUBSET of the node table is
refused because edge endpoints are positions in the full table — a subset walks between the wrong
genes and stays in range once it is large enough, which is the silent version of the same bug.

**Every supervised method is scored only on genes no fold of its model ever saw**, with a test that
attacks it: forty columns of pure noise must not be "recovered". The shuffled-label negative control
the clustering carries applies to each of them.

## The GNN is deliberately not built, and this is the reason

This instruction said to build it last and only "if the baseline says the simpler methods are leaving
something on the table". The baseline now exists and says the opposite: a penalised linear model and
a random walk over a two-thousand-edge layer already quadruple the map's recovery, and boosting beats
both. A relational GNN on 8,140 nodes would overfit, would be the least interpretable thing in the
project, and would need a heavy dependency to reach numbers the simple methods have already reached.

Left open rather than closed off: if the acquisition campaign grows the graph substantially, or if a
question appears that needs to combine layers rather than compare them, the gate this instruction set
would be met and it should be built then.

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
