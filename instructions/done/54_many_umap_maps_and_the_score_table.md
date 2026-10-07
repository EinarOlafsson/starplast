# 54. Many more pregenerated maps, searched for structure, and a navigable maps panel

Opened 2026-09-30, from user feedback on 0.49.0, verbatim:

> there should be way more UMAPs in maps, basically every category by itself and then logical
> combinations, finding max structure should be a goal in making all these UMAPs. the table below
> only changes when the first UMAP is chosen in maps. the second choice and so on don't change the
> table.

Three parts.

## 1. The score table only updated for the first map chosen — BUG

**Reproduced** before any change, offscreen, on the real `Window` (not a stub), by clicking list
items through `QListWidget.itemClicked` and comparing the table's cell *text* between clicks:

| Scenario | Table changed on the 1st map click | on the 2nd |
|---|---|---|
| The panel's **default** view, "One label, every map" | **no** | **no** |
| "One map, every label" | yes | yes |

**Cause.** Not a signal connected once and not a build-once flag. `refresh()` filtered the table by
the view and the label only:

```python
if by_label:
    s = s[s.label == self.label_box.currentText()]   # no mention of the chosen map
```

In the label view the rows ARE the maps, one per map, for one label — so they are the same table
whichever map is current, and nothing in it moved when a map was picked: no marked row, no selected
row, no change to the sentence underneath. That view is the panel's default, which is why the user
saw a table that never responded. (The "first choice appeared to work" because the very first click
is also what first puts a gallery map in the central view and first fills the sentence line.)

**Fix.** Both views now follow the current map. The label view marks the map on screen with
`SHOWN_MARK` on its row, selects and scrolls to that row, and the line under the tree names the map
with its recipe and structure score; the map view is filtered to the current map as before. Choosing
a map through the list, through the map combo (`_map_box_changed`, which previously only refreshed
when the index moved) or by double-clicking a row all go through `select_map`.

**Regression test.** `test_choosing_a_second_map_changes_the_score_table` picks map A, then B, then A
again, **in both views**, and compares every cell's text: A ≠ B, B ≠ A-again, and A-again == A. It
also asserts exactly one row carries the marker and that it follows the chosen map.

## 2. Many more maps, tuned for structure

`umap_gallery.recipes` now emits six groups (`GROUPS`) instead of three kinds:

1. all measurements at n_neighbors 10 / 25 / 60 — the fixed reference maps, settings declared;
2. all but one family, for every substantial family (not only localization);
3. one map per evidence family;
4. one map per individual experiment block with at least `MIN_BLOCK_COLUMNS` columns;
5. pairs and 6. triples of families, from the curated `COMBINATIONS`, each carrying the biological
   question it asks.

Genes are placed at the strictest of `COVERAGE_STEPS` that still places `MIN_SET_GENES` genes, so a
dense family is built over well-measured genes and a sparse one is still buildable.

**Maximising structure.** Every map but the three reference maps is searched over `MAP_GRID` ×
`GALLERY_CLUSTER_GRID` and the best `structure` ships. `structure` is label-free and is the geometric
mean of `search.map_quality`'s score (partition: clustered share × evenness, zero on a degenerate
clustering) and the mean silhouette of the clustered points rescaled to 0–1 (geometry). Either alone
is gameable; the geometric mean requires both. The search is successive halving as `search.tune_umap`
does it. Every setting tried is recorded in the manifest beside the winner.

Nothing in `search.py`, the strategy catalogue or the calibration was changed: `map_quality` is read,
not modified.

## 3. A navigable panel

`QTreeWidget` grouped by `GROUPS`, a filter box, a sort box (gallery order / structure / how well the
chosen label maps / genes placed) and a Group checkbox that flattens it into one best-first list. Each
row shows the map's size, cluster count, structure score and the chosen label's skill; its tooltip has
the description, the full recipe and the number of settings searched. Thumbnails are drawn lazily per
expanded heading and subsampled to `THUMB_POINTS`, because `gallery.thumbnail` paints point by point.
Every control has a tooltip.

## Done when

- [x] The bug is reproduced, its cause written down, fixed, and covered by a test that compares cell
      values in both views.
- [ ] Both organisms' galleries are rebuilt with the expanded, searched catalogue; coordinates
      float32, clusters int16, the three files under 25 MB and the wheel under 75 MB
      (`scripts/check_wheel.py`).
- [ ] The build is an executed notebook (`notebooks/umap_gallery_2026_09_30.ipynb`) recording, per
      map, its recipe, the settings searched, the structure score and the genes placed.
- [ ] `docs/guide.md` documents the groups, the structure score, the navigation and the command to
      build more locally.
- [x] Offscreen screenshots of the panel with many maps and of the table changing between maps, looked
      at and iterated on.
- [ ] Full suite green.

## Not done, and why

- The gallery is **not** trimmed to the best N per category: the whole catalogue fits the size budget,
  so trimming would lose maps for nothing. The documented command for building an even larger gallery
  locally is in `docs/guide.md` instead.
- Single-column blocks get no map of their own. A 3D UMAP of one measurement is a line; those blocks
  are covered by their family's map.
