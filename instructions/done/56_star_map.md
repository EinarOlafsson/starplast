# 56 · Star map: gene-gene links by source, with provenance, including the user's own runs

**Status: complete (2026-09-29, worktree branch off nightly 0.48.0; not merged, not pushed).**

The user, 2026-09-29:

> there should be a network map the user can navigate where they can visualize how different
> strategies link together different genes, like a star map, where each gene is a node and its
> edges are connections to other genes through data or inference done in starplast; their own runs
> should show up here as well as the runs already done by starplast.

## Plan

1. `starplast/star_edges.py` (Qt-free): one link store with provenance (kind measured/inferred,
   group, source, how, setting, score, strength, run, origin, created, note). Measured links read
   from each space's graph; inferred links derived from strategy results, bounded.
2. `scripts/build_star_edges.py` -> `starplast/data/star_edges.parquet`, recorded as
   `notebooks/star_edges_2026_09_29.ipynb`: default-setting runs of 13 strategies on both shipped
   tables.
3. `UserEdgeStore` under `paths.user_cache_dir()/star_edges`: every Strategies-tab run adds its links.
4. `starplast/star_map.py`: the view (radial, 1-2 hops, per-source toggles and counts, hover
   provenance, click re-centres, follow selection, select on 3D map). `app.py` gets a 3-line hook
   (`_star_map`) and 2 lines in `on_pick`.

Rules: no strategy algorithm or calibration number changes; no new organism literals.

## Bounds (module constants in `star_edges.py`)

* Pair tables (16, 17, 18, 33, 34): taken as listed (already top-N). `measured edges, ranked` and
  `raw ranking` are skipped (restatements).
* Modules/clusters (01, 04, 09, 15, any result with `labels`): each member -> `MODULE_K`=3 nearest
  co-members, first by the number of measurement layers joining them, then by distance on the
  strategy's map (or the permitted measurements when it has none).
* Set expansion (20, 25): each candidate -> `SEED_K`=3 nearest seeds, same rule.
* Neighbour vote (07): each call -> `LABEL_K`=3 nearest labelled genes carrying the called label,
  in `Context.features(target)` (the vote's own space).
* Partner vote (13): the partners it names, at most `PARTNER_MAX`=5.
* No run adds more than `MAX_EDGES_PER_RUN`=20,000 (Tg multiplex_modules reaches this cap).

## Verified

* Shipped file 1.16 MB, 87,115 links (Tg 45,081; Pf 42,034; per source in the notebook). Wheel was
  55 MB before; limit 75 MB.
* `tests/test_star_map.py` (18 tests): provenance of measured and inferred links, skipped tables,
  module caps, the per-run cap, seed/label/partner links, per-source and node caps in `star`, user
  store round trip and per-organism isolation, shipped file size and provenance, tooltips on every
  control, click re-centres and Back, hover text, a user run appears at once, toggles, follow.
* Screenshots looked at (offscreen): window with the star map dock, 1-hop with hover on a user-run
  link, 2-hop.

## Decided not to

* Draw modules as cliques (a 1,800-gene community would be 1.6 M links).
* Ship per-strategy colours for user runs: one "your runs" group, one toggle; the hover names the
  strategy and run.
* Add strategies that call genes one at a time (19, 21, 35-39 ...): they state no gene-gene link.
