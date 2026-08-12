# The UMAP gallery: see every embedding, as it is generated — DONE 2026-08-12

Requested: grid and scroll modes, thumbnails appearing as they compute, click to expand into the
central view where a gene can be clicked like any other map. The hard part named in the task was the
contract change; that is what was done first.

## The contract change

`tuning.walk_umap_iter` is the walk now, and it **yields one `WalkStep` per configuration** as each
is computed. `walk_umap` is that collected and ranked, keeping its old signature and return, plus an
`on_step` callback — so every existing caller is unaffected and a caller that wants to watch does not
have to re-run the sweep to see it.

A step carries the scores AND the map: `coords`, and `genes` as a boolean mask over the node table.
The mask is the part that matters. A walk embeds a subsample, so without it there is no way to say
which gene a coordinate belongs to, and a map whose points cannot be named cannot be clicked into.

Two things had to be fixed to make that mask true rather than approximately true:

- **The subsample is now taken in table order.** `rng.choice` returns its picks in random order, so
  row *i* of the embedding was an arbitrary gene: the mask, the saved gene ids and the clicked point
  each named a different one. Sorting picks the same genes — same seed, same set — in the order the
  table has them. `search.py` already did this; the walk did not.
- **Each configuration is written through `EmbeddingStore` as it is computed**, with the
  hyperparameters that run actually used. A walk stopped half way now leaves behind everything it
  finished, with the recipe to rebuild each one.

## The gallery

`starplast/gallery.py`, a dock along the bottom, fed from `AnalysisPanel.walk_step` (emitted from the
worker thread, queued to the GUI thread by Qt). Grid mode is a wall of thumbnails at one size on one
background, so a difference between two tiles is a difference between two maps; scroll mode is one
large view with a slider stepping through in order. Both fill in as the walk runs, and the walk table
on the Map tab now fills a row at a time as well, ranked only when the whole sweep is in.

In scroll mode the view follows the newest map **only while the slider is already at the end** — the
log-tail rule. A walk that yanked the view to run 13 while run 12 was being read would be unusable
during exactly the period the gallery exists for.

**Thumbnails are painted, not grabbed from GL.** The task suggested `grabFramebuffer`; 288 offscreen
framebuffers need a live GL context each, are slow, and fail on the headless machines the suite runs
on. A thumbnail is a small orthographic projection of the same coordinates: deterministic, fast
enough to keep up with a walk, and testable without a display. Both axes take the same scale, so a
map that is genuinely elongated still looks elongated. Colour comes from `Window.colours`, so a
thumbnail is coloured by whatever the map is coloured by and grey means unknown in both.

Clicking a tile calls `use_embedding`, the same path as "build this map" — so the expanded map *is*
the central view, with picking, colour modes and edges. Not a picture of a map.

## The bug this uncovered

`use_embedding` left every gene the embedding did not cover **at the origin**. With a full-proteome
build that is nothing; with a walk configuration it is 7,340 of 8,140 genes stacked in the middle of
the map, drawn at 6% alpha, pickable, counted in the status bar and included in the class filter.
That is absence rendered as a value, in the application whose first constraint is that it must not
be. Unplaced genes are now recorded in `Window.placed`, drawn at alpha 0, excluded from
`visible_mask` (which is what points, centroids, edges and the gene count all read), and excluded
from picking via `Map3D.pickable` — hiding alone is not enough, since an invisible point is still the
nearest point to a click where it sits.

`embedding.normalise` was also lifted out of `embed` so a walk map arrives at the same extent as any
other, and its centring moved to float64 with the cast to float32 last: done the other way round, a
cloud whose spread is small beside its offset from the origin loses that spread to cancellation.

## Verified

- 1,290 tests pass headless, including 40 new ones. `gallery.py` and `tuning.py` are at 100%.
- Rendered under Xvfb against the real cache and **looked at**: four configurations as four
  distinguishable thumbnails with their scores, scroll mode, and an expanded configuration in the 3D
  view whose shape matches its thumbnail. Looking is what caught the mode chooser still reading
  "grid" while the scroll page was showing.
- Run end to end under the pandas 3 interpreter on the real 8,140-gene table: walk streams, table and
  gallery fill, expanding places 600 genes and hides 7,540 at alpha 0, a gene in the expanded map
  clicks through to its evidence panel.

## What this unblocks

Task 16 (automated walk: UMAP → clustering → per-category scoring) needs exactly this contract; its
per-configuration results now have somewhere to appear and something to be looked at.
