# Full proteome: the subsample was flattering

Every number in `../search_2026_08_12_corrected/` was computed on a seeded 3,000-gene subsample. This
run scores cell-cycle phase on all 8,140 genes, with all 873 measured phases in play rather than the
316 that happened to fall in the sample.

## The negative result holds, and strengthens

| | subsample (3,000 genes) | full proteome (8,140) |
|---|---|---|
| mean F1 | 0.489 | **0.396** |
| best cluster purity | 52% | **31%** |
| predictions at the 80% bar | 0 | **0** |

Across the ten best configurations, no cluster exceeds **31%** purity. The bar is 80%. So "the map
recovers cell-cycle phase but cannot predict it" is not an artefact of subsampling — it is more clearly
true on the whole proteome than on part of it.

## Why the subsample scored higher

Smaller samples produce tighter, purer clusters. With 3,000 genes a cluster holds fewer members and is
more likely to be dominated by one phase; with 8,140 the same structure absorbs more of everything.
Nothing about the subsample was wrong — it was seeded, reproducible and honestly reported — but its
numbers are **optimistic relative to the full data**, and that applies to all four targets rather than
only this one.

That is the reason the whole table was re-run at full size: see `full_summary.json` and the per-target
CSVs here.

## The cost objection was wrong

90 full-proteome runs took 269 seconds — about three seconds each, the same as the subsample. Sampling
was buying nothing. A paper should not quote subsample estimates when the real thing is this affordable.

## Files

| file | what |
|---|---|
| `full_<target>.csv` | every run's score and settings, with the excluded columns per row |
| `fulllab_<target>.csv` | per-label precision, recall and F1 |
| `full_summary.json` | the winning configuration per target |
