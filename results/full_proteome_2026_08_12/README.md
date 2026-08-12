# Held-out recovery on the full proteome

All four targets, 328 configurations each, over all 8,140 genes rather than the
3,000-gene subsample the earlier searches used. Absence labels excluded
throughout, as in `../search_2026_08_12_corrected/`.

The subsample was flattering, and the honest summary is that at full scale the
map recovers one label: the majority stage class.

## Headline

| target | role | mean F1 | clusters | noise | labels recovered |
|---|---|---|---|---|---|
| `stage_enriched_derived` | POSITIVE control | 0.675 | 2 | 0.16 | **1 of 3** |
| `cellcycle_phase` | measured | 0.398 | 3 | 0.05 | **0 of 5** |
| `attention_depth` | NEGATIVE control | 0.322 | 9 | 0.52 | **0 of 3** |
| `compartment` | measured | 0.193 | 30 | 0.11 | **1 of 24** |

Against the subsample: 0.709 → 0.675, 0.489 → 0.398, 0.339 → 0.322,
0.228 → 0.193. Every target fell, and the two measured targets fell furthest.

## Why mean F1 alone would have misled

`stage_enriched_derived` scores 0.675, the highest of the four, from a
**two-cluster** solution in which all three stage labels best-match the *same*
cluster. Its per-label recalls are 1.000, 1.000, 1.000 — the signature of one
cluster containing everything that was labelled. Precision is then just class
prevalence: oocyst is 1,348 of 1,877 labelled genes (72%), so oocyst alone
scores F1 0.838 and drags the mean up.

Mean F1 over labels is inflated by a trivial partition whenever one class
dominates. `n_labels_recovered` is not, which is why the search records it, and
it is the column to read first. On that column the positive control recovers the
majority class and nothing else, and neither measured target recovers anything.

This is not a degenerate *clustering* — HDBSCAN found 2 clusters, not 1. It is a
degenerate *match*: the labels do not distinguish the clusters that exist.

## What this does and does not say

It does not say the map is empty. The ranking is stable across both scales and
in the expected order — derived stage above measured cell cycle above the
attention control above localisation — and a ranking that survives a 2.7×
increase in genes is worth something.

It does say that **no claim of the form "the map recovers X" is supported at full
scale**, because at the threshold the search uses, cell cycle recovers 0 of 5
phases and localisation 1 of 24 compartments. The subsample result that cell
cycle (0.489) cleanly cleared the negative control (0.339) narrows to 0.398 vs
0.322 on the full proteome, and the margin is carried by labels that are not
individually recovered.

`compartment` remains below the negative control at both scales. Localisation is
recovered worse than how much a gene has been studied.

## Files

`full_<target>.csv` — one row per configuration, 320 each. `mean_f1`,
`n_clusters`, `noise_frac`, `n_labels_recovered`, `n_labels_scored`, and the
`excluded` column naming the features withheld to break circularity.

`fulllab_<target>.csv` — per-label precision, recall and F1 for the
best-matching cluster, for every configuration. This is where the recall-1.0
signature above is visible.

`run_full_proteome.py` — the script that produced all of the above, committed
beside its output. `search.search()` is a library function with no command line,
so the script is the only honest record of the grid that was run: it passes
`sample_size=None`, which is what makes this the full proteome rather than a
subsample. Run it as `python run_full_proteome.py <output-dir>`.
