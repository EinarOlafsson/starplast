# Structure search, corrected scoring

The same four targets as `../search_2026_08_12/`, re-run after `score_recovery` began excluding labels
that mean *not measured* (`unassigned`, `unknown`, empty, `nan`).

The first run's table was computed under a different rule, so none of its numbers are comparable with
these. It is kept as the evidence that justified the change.

## Why a re-run was needed rather than a recomputation

The exclusion changes which configuration *wins*, not only the winner's score. Once `unassigned` stops
inflating a run's mean, a different combination of blocks and hyperparameters can come out on top — so
the honest comparison required running the whole grid again, 328 runs per target.

## Method

Identical to the first run: 41 dataset combinations × n_neighbors {15, 50} × min_dist {0.0, 0.25} ×
min_cluster_size {25, 60}, on a seeded 3,000-gene subsample (seed 42), median NA policy, rank scaling.
Every scorable run's embedding is stored with its recipe.

## Results

_(filled in below when the run completes)_
