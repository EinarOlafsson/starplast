# 34 — Two metric columns that are wrong, and one summary line that misleads

**Status: open. Not started.** Small, independent, and each one currently prints a number nobody
should trust. Written 2026-08-14.

## 1. `mean_auprc` is nonsense for a continuous layer

Seen in `bigB_00_guilt_fit_invitro_hff`: `mean_auprc = 0.0032`.

`metrics.ranking` treats the truth column as categorical — `_classes` takes `value_counts` and keeps
anything appearing at least ten times. Handed a float column like `fit_invitro_hff`, every distinct
value becomes its own "category", and whatever survives the count threshold produces curves that
mean nothing. The yield score is unaffected (it never consults these), but the column is written into
every saved search and into the results table the tab displays.

Fix: `metrics.report` and `metrics.ranking` should detect a numeric truth column and return NaN for
the ranking metrics rather than computing them, the way `curve_scores` already returns NaN for a
class with no negatives — "a question that was not asked" rather than a fabricated number.

Better, if it is wanted: for a continuous layer the analogous quantity is a **rank correlation**
between a gene's leave-one-out neighbourhood mean and its own value. That is the continuous version
of the same claim ("do neighbours predict this gene's value") and would fill the column with
something real. Optional; NaN is the minimum acceptable fix.

## 2. `trustworthiness` is never computed

`metrics.report(labels, truth, coords, X=None)` computes trustworthiness only when handed the
feature matrix `X`. `optimize.evaluator` calls it with `coords` and never with `X`, so the column is
NaN in every row of every search run so far.

The awkwardness is that `embedding.embed` returns `(coords, names, rows)` and builds the matrix
internally, so `evaluator` has no `X` to pass without building it twice. Options, in order of
preference:

1. Have `embed` optionally return the matrix it built, and cache it in `optimize.evaluator`'s
   embedding cache alongside the coordinates. One extra array per distinct embedding.
2. Call `build_matrix` in `coords_for` and pass it to `embed` — needs an `embed(X=...)` entry point.
3. Drop the column. Defensible: trustworthiness is subsampled and quadratic, and it answers "did the
   projection invent the clusters" which nothing currently optimises. If it is dropped, remove it
   from `metrics.report` too rather than leaving a parameter nobody passes.

Whichever is chosen, no results table should carry a column that is structurally always NaN.

## 3. The CLI summary line counts the wrong findings

`discover._summary` prints `len(run.findings)`, which is every finding from **every configuration**
in the run — 2,315 for a search whose winner found 42. It reads as the winner's count and is not.

Fix: report the winner's own count, and the run-wide total separately if it is wanted at all:

    bigA_00_guilt_compartment_best: best 254.67, 100 configs, 42 findings (2315 across all configs)

`searches.Search.for_config(0)[1]` is the winner's findings and is already the right accessor.

## Acceptance

Each of the three is independently testable and small. 100% coverage, no `pragma`. Worth one commit
each rather than one commit for all three — they touch different modules and only the third is
cosmetic.
