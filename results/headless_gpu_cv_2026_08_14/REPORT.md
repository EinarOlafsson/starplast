# Leakage-closed headless GPU validation — 2026-08-14

## Outcome

The expanded Toxoplasma cache was tested with cuML 26.08.00 on an RTX 3090. The walk evaluated 120
three-dimensional UMAPs and 2,280 UMAP/clustering structures over the fixed 8,140-gene proteome. Every
reported score comes from five locked outer folds after structure selection on inner validation labels.
Cluster-to-category calls are learned from training labels only and may abstain. Twenty-five label
permutations repeat the entire nested selection procedure for each experiment.

| held-out target | exclusion scope | structures | nested macro-F1 | coverage | permutation p | proposed rows |
|---|---|---:|---:|---:|---:|---:|
| localization (`compartment`) | **strict biology class** | 608 | 0.0456 ± 0.0052 | 7.6% | 0.0385 | 233 |
| localization (`compartment`) | same target family only | 608 | 0.0446 ± 0.0137 | 7.4% | 0.0385 | 174 |
| cell-cycle phase | **strict biology class** | 456 | 0.1357 ± 0.0207 | 23.0% | 0.0385 | 972 |
| cell-cycle phase | same target family only | 608 | 0.1832 ± 0.0405 | 30.4% | 0.0385 | 1,525 |

`0.0385 = 1 / (25 + 1)` is the smallest p-value this permutation count can resolve: none of the 25
selection-aware null runs matched the observed score. Statistical separation from shuffled labels does
not make the low macro-F1 practically strong. The strict biology-class experiments are the primary
results; target-family exclusion is a deliberately narrower sensitivity analysis.

## What can be proposed responsibly

Strict localization holdout selected an all-independent, 187-feature UMAP (`n_neighbors=60`,
`min_dist=0`) with HDBSCAN. Its useful signal is narrow:

| category | held-out precision | recall | F1 | proposed unlabeled genes |
|---|---:|---:|---:|---:|
| PM - peripheral 1 | 0.900 | 0.282 | 0.413 | 186 |
| dense granules | 0.704 | 0.300 | 0.418 | 5 |
| nucleus - chromatin | 0.574 | 0.121 | 0.199 | 42 |

Thus, 191 localization proposals belong to categories with both precision ≥ 0.60 and F1 ≥ 0.20.
They are candidates for orthogonal follow-up, not measured localizations. Most localization categories
received no defensible call.

Strict cell-cycle holdout selected a measured-only, 24-feature UMAP (`n_neighbors=30`, `min_dist=0`)
with k-means clustering (`k=80`):

| phase | held-out precision | recall | F1 | proposed unlabeled genes |
|---|---:|---:|---:|---:|
| M | 0.647 | 0.218 | 0.323 | 177 |
| C | 0.485 | 0.151 | 0.224 | 338 |
| S | 0.614 | 0.067 | 0.114 | 457 |
| G1b | 0.050 | 0.011 | 0.018 | 0 |
| G1a | 0.000 | 0.000 | 0.000 | 0 |

Only the 177 M-phase proposals meet both precision ≥ 0.60 and F1 ≥ 0.20. S-phase precision is
moderate but recall and F1 are low; G1 phases are not inferable by this result. The higher sensitivity
score when transcription features remain available quantifies why the strict gene-expression holdout
must remain the primary analysis.

## Search space

- Feature collections: all independent slots, measured-only slots, broad-coverage slots, and
  non-transcription slots where that collection survived the holdout.
- UMAP: 3 components; `n_neighbors` 15, 30, 60, 120; `min_dist` 0 or 0.2; rank scaling; median
  imputation; seed 42.
- HDBSCAN: minimum cluster size 20, 50, 100, 200; minimum samples default or 10; EOM or leaf
  selection.
- K-means: 20, 40, or 80 clusters.
- Calls: minimum five training labels and cluster purity ≥ 0.55; otherwise abstain.

## Files and interpretation

Each experiment directory contains `manifest.json`, all structure recipes, inner-validation scores,
locked outer-fold and per-category scores, the permutation null, the best coordinates/clusters, and
`inferred_genes.csv`. Inference rows explicitly say `inferred_from_held_out_structure` and carry cluster
support/purity plus nested-CV precision, recall and F1. They must not be merged into measured annotation
columns.

This is transductive unsupervised validation over a fixed proteome: all genes contribute their
non-target features to UMAP, but target labels are excluded from features and used only inside the
locked label splits. The analysis validates the complete selection procedure, not a biological causal
mechanism and not each individual candidate independently.
