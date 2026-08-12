# Structure search, 2026-08-12

Four targets, 328 runs each: 41 dataset combinations x n_neighbors {15, 50} x
min_dist {0.0, 0.25} x min_cluster_size {25, 60}, on a seeded 3,000-gene
subsample (seed 42).

| target | role | best mean F1 |
|---|---|---|
| `stage_enriched_derived` | positive control (computed from expression) | 0.709 |
| `attention_depth` | **negative control** (study effort) | 0.654 |
| `cellcycle_phase` | measured | 0.489 |
| `compartment` | measured | 0.484 |

**The negative control outscored both measured targets**, and the per-label
tables here say why: essentially all of it comes from recovering the
never-named class at F1 0.77, while the real attention tiers score 0.14-0.35.
The map separates MEASURED from UNMEASURED genes, which is a property of
missingness rather than of the literature.

The same artefact inflated `compartment`: its best label was `unassigned` at
0.39, and excluding that pseudo-label the run scores **0.207**, with real
compartments at 0.21-0.35.

`cellcycle_phase` has no absence class, so its 0.489 stands: C 0.54, S 0.53,
M 0.37, G1b 0.25.

These tables were produced BEFORE `score_recovery` began excluding absence
labels. They are kept as the evidence for that change. Re-running under the new
rule is the first item in `instructions/done/10_search_more_targets.txt`.
