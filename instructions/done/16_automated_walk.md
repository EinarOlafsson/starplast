# Automated walk: UMAP, then clustering, then per-category scoring — DONE 2026-08-12

Built on what task 15 and task 21 already provided rather than as a new subsystem: `search.search`
already embedded, clustered at each `min_cluster_size` and computed a per-label table. What was
missing was that none of it could be watched, none of it could be looked at, and the per-label table
never reached the interface.

## What was added

- **`search` streams.** `on_run` is called with a `RunStep` as each EMBEDDING finishes, carrying its
  best clustering among the sizes tried, its per-category scores, its coordinates and its labels.
  Per embedding rather than per row: the rows for one embedding differ only in the clustering, and a
  gallery of the same map five times is not a gallery. Every row still reaches the table — what is
  streamed is what can be looked at, not a subset of what was scored.
- **The winner is chosen by the objective in force**, not by a fixed metric. "Best" is not a fixed
  thing; a step that ignored the objective would show a different map from the one the ranking names.
- **A per-category table in the Search tab**, filling as the walk runs: one row per configuration and
  category, each carrying enough of the recipe to rebuild its map. A mean over categories hides the
  case this table exists for — a map where the GRAs are clean and everything else is a mess is
  exactly what you want when you are looking for GRAs.
- **The frontier.** `search.frontier` marks the rows nothing beats on both mean F1 and best F1. The
  task asked for the two targets offered rather than chosen, and ranked on each; ranking on each is
  two sorts, and the frontier is what neither sort shows — a configuration second on both is often
  the one to use and tops neither list. Marked, not filtered: a dominated configuration is still a
  result.
- **Clusters are visualisable**, as the task required. Each configuration appears in the gallery
  coloured by the clustering that was scored, with noise grey for the same reason grey means unknown
  everywhere else. Clicking a row of either table rebuilds that exact configuration into the central
  view, and a per-category row also names the cluster its score is about.

## What was deliberately not done

**"Keep the best clustering" is not "discard the others".** Every `min_cluster_size` keeps its row;
the best is what the gallery shows and what the step carries. Discarding the rest would hide the
sensitivity to a parameter the project already knows is decisive, and the table is where that shows.

## Verified

- 1,469 tests pass headless; `search.py` at 100%.
- A real search on the fixture emits one step per embedding, whose clustering is the highest-scoring
  one under the objective in force, at the same coordinate scale as any other map.
- The frontier's arithmetic is pinned by a case with a genuinely dominated row — the first version of
  that test asserted the wrong answer, and the code was right.
