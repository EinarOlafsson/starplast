# Automated walk: UMAP, then clustering, then per-category scoring

Depends on [15_umap_gallery.md](15_umap_gallery.md) for its output to be inspectable.

For every point in the UMAP hyperparameter walk:

1. build the embedding,
2. run a **HDBSCAN hyperparameter search** on it and keep the best clustering,
3. compute precision, recall and F1 for **each category** of the held-out variable.

The result is "how well does each category map onto a cluster, for each configuration" — which is
what lets a user pick a map where the GRAs are clean even if everything else is a mess.

## Three optimisation targets, offered rather than chosen

- best single-category mapping (find a map that nails one compartment),
- best average across categories (find a map that is good everywhere),
- both, ranked on each, showing the frontier.

`starplast/objectives.py` already implements the scoring, including the guards. Use it rather than
writing new metrics — see [21_wire_objectives_into_search.md](21_wire_objectives_into_search.md),
which must land first.

## What must not be reimplemented

`search.score_recovery` already excludes absence labels, and absence labels are what made the
negative control read 0.654 instead of 0.339. It also reports `n_labels_recovered`, which is the
column that matters: mean F1 is inflated by a single cluster whenever one class dominates, which is
exactly the failure a per-category table invites.

Clusters must be visualisable in this mode too, not only categories.

## Done when

A walk produces a table of (configuration x category) scores that can be sorted by any of the three
targets, and clicking a row opens that embedding in the gallery with its clustering shown.
