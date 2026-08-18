# 46 — A paper-ready PDF for every recipe run

**Status: open. Requested by the user 2026-08-18.**

## What to do

At the end of a recipe run the user must be able to (a) evaluate the result in the application, as
they can now, and (b) open a **PDF** built for a paper: the figures that would go in a manuscript,
not a screenshot of the app.

## The layout the user asked for

**One row per UMAP.** A run tunes several maps and clusters each; the interesting ones all belong in
the document so the reader can see what was tried, not only what won. Each row, left to right:

1. **A hyperparameter table for the UMAP** — n_neighbors, min_dist, metric, components, seed, and
   which library built it (a map built by cuML's UMAP is a different map of the same data).
2. **A hyperparameter table for the clustering** — algorithm, min_cluster_size, min_samples,
   selection method, and the degeneracy verdict.
3. **The UMAP coloured by cluster.**
4. **The UMAP coloured by the primary holdout label.**
5. **The UMAP coloured by the validation holdout label.**
6. **Overlap-quality graphs across ALL labels** — how well each label maps onto the clustering, for
   the primary and for the validation holdout, so the two can be read against each other.
7. **Detail panels for a few interesting labels** — the ones that recovered best, and the ones where
   the primary and the control disagree, which is where a reader learns whether to believe the run.

## What "quality of cluster vs label vs validation label" has to show

The point of the figure set is the three-way comparison, so it must be legible as one:

* per label: precision, recall, F1 against its best cluster;
* per (label, cluster) pair: the full matrix, since a label split cleanly across three clusters is a
  finding and a single best-cluster number hides it;
* the same for the validation label on the SAME clusters;
* and the agreement between them per cluster — the number that says whether the answer is
  corroborated or standing on the primary alone.

## Constraints

* **Vector output.** A figure destined for a paper is not a PNG of a screen.
* **The document must be self-describing**: the question, the inputs after leakage closure, what was
  excluded and by which mechanism, the seed, and the version. A figure a reader cannot trace back to
  a rebuildable run is not evidence.
* **A 3D embedding drawn in 2D must say which components it is showing**, rather than silently
  projecting and letting a reader assume they are seeing the map.
* Headless: it must build with no display, because the suite has none.

## Done when

* A recipe run produces a PDF with one row per retained UMAP, in the order above.
* The PDF is reachable from the application at the end of a run, and from a script.
* The three-way comparison (cluster vs label vs validation label) is legible in one figure set.
* The document names the question, the closure, the seed and the version.
* It builds headless, and a test asserts the file is a valid multi-page PDF with the expected pages.
