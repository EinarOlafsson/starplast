# The same battery, with the circularity guard fixed

All four targets, 328 configurations each, over all 8,140 genes — the run script and grid of
`../full_proteome_2026_08_12/`, re-run after three leaks were found in the guard on 2026-08-13
(see `HANDOFF.md` 3f-2 and `instructions/done/29_circularity_guard_leaks.md`).

## The answer: nothing moved

| target | role | 2026-08-12 | this run |
|---|---|---|---|
| `stage_enriched_derived` | POSITIVE control | 0.675 | **0.675** |
| `cellcycle_phase` | measured | 0.398 | **0.398** |
| `attention_depth` | NEGATIVE control | 0.322 | **0.322** |
| `compartment` | measured | 0.193 | **0.193** |

All four identical to the last digit — winning block combination and per-label best F1 included
(0.839, 0.467, 0.383, 0.500). That is the
result, and it is worth as much as a change would have been: the published table is not affected by
the leaks.

## Why not, and where the leaks DID reach

`search.search` with `block_sets=None` sweeps a hard-coded base of six blocks —
`expression_summary`, `expression_raw`, `fitness_screens`, `published_screens`, `protein_features`,
`interactions`. It contains neither `localization` nor `literature`, so the leaking columns
(`lopit_prob_map`, `lopit_prob_mcmc`, `lopit_methods_agree`; `n_papers_focal / substantive /
incidental`) were never available to a published row. The sweep now names the blocks it used and the
ones it left out, so this is readable off the log rather than deduced from the source.

**The interface is the other case.** `AnalysisPanel.run_search` builds its combinations from every
block that has columns in the table, `localization` included. Run from the Search tab over a
600-gene subsample, that block won outright — mean F1 0.259, against 0.192 for the best measurement
block — which is how the leak was found and is exactly the configuration a user would have believed.

## What changed in the guard

| test | catches | added |
|---|---|---|
| measured association ≥ 0.8 | the undeclared copy | v1.4 |
| declared derivation | the joint function (argmax of three columns) | v1.4 |
| shared provenance | the same experiment's OTHER outputs | 2026-08-13 |
| shared quantity | the same thing estimated another way | 2026-08-13 |

Excluded for `compartment` now: 14 columns rather than 5. For `attention_depth`: 7 rather than 2.
For `cellcycle_phase`: 2 rather than 1.

## Reproducing

    python run_full_proteome.py <output directory>

~25 minutes per target on 32 cores. The saved embeddings are not committed — each one is rebuildable
from the recipe in its row, and 984 npz files are a working directory rather than a result.
