# 58 · Star map: stop the hover twitch, show the larger network, let the user define the connections

**Status: implemented (2026-09-30, worktree branch `worktree-agent-a1bd3d8ee164b2ab5` off nightly
0.49.0; not merged, not pushed).**

The user, 2026-09-30:

> the new network feature is great but when a gene is hovered the graph becomes twitchy. I'd also
> like a way to visualize the larger network and define the connections in the network view.

## 1. The twitch: what it actually was

Reproduced offscreen before any change was made, by centring a gene and sending one `QTest.mouseMove`
over a link: the view's transform went from scale 0.852823 to 0.848790 and the viewport lost two
pixels of height. Three things combined.

1. **The hover text resized the view.** `self.info` is a word-wrapping `QLabel` in the same vertical
   layout as the view. A hover wrote a longer or shorter provenance string into it, the label's
   size hint changed, the layout gave the view fewer or more pixels, `_View.resizeEvent` fired and
   called `panel.fit()`, and `fitInView` re-scaled and re-centred the whole graph. **This is the
   twitch.** It happened on every pointer move onto a different item.
2. **`fit()` framed `scene.itemsBoundingRect()`**, which is not constant: `_Edge.hoverEnterEvent`
   widened the hovered pen by two units, which grew that item's bounding rectangle, which grew the
   scene's — so even a hover that did not change the label's height changed what the next fit framed.
3. **Nothing was debounced.** Sweeping the pointer across a dense star fires hover events at the
   mouse's rate, and each one rebuilt rich text and could trigger 1 and 2.

### The fix

* `self.info` and `self.headline` have fixed heights (`INFO_HEIGHT`, `HEADLINE_HEIGHT`) and an
  `Ignored` horizontal size policy: what they say can no longer change the geometry of anything.
* `_View.resizeEvent` re-frames only when `ev.oldSize() != ev.size()`.
* `fit()` frames `self._frame_rect`, computed in `redraw` from the layout's own positions and from
  nothing else.
* `_Edge.hoverEnterEvent` changes the pen's **colour** only; `_Edge.boundingRect` adds a fixed
  `HOVER_MARGIN`, so the item's rectangle never depends on the hover state.
* Hover text is queued and rendered by a single-shot `QTimer` (`HOVER_DEBOUNCE_MS = 40`).
* Layouts are cached on `_layout_key()` = (mode, centre, definition signature, depth, per source,
  hops, gene cap), and the link-filtering result on the definition's signature. Hover is in neither.

**Done when:** `test_hovering_every_visible_item_never_moves_the_graph` hovers every `_Node` and
`_Edge` of a drawn star and requires the view transform, the visible scene rectangle, the scene's
bounding width and every node position to be identical, to the bit, before and after each one;
`test_hovering_the_large_view_moves_nothing` does the same for the canvas.

## 2. The larger network

`star_edges` (Qt-free) gains `neighbourhood`, `overview`, `clusters` and `layout`; `star_map` gains
a mode box and `_NetworkCanvas`.

| Mode | What it shows | Caps |
|---|---|---|
| one gene (star) | unchanged | `MAX_NODES` 160 |
| neighbourhood | breadth-first from the centre, strongest links first | `MAX_HOPS` 4, `NEIGHBOURHOOD_MAX` 2,000 genes |
| whole network | everything the rules leave, 12 largest clusters coloured | `OVERVIEW_MAX` 3,000 genes, `LARGE_EDGE_MAX` 30,000 links |

Degradation is by keeping the strongest links and then the best-linked genes, and `overview` returns
the sentence that says so; the view prints it. `clusters` is weighted label propagation visiting
genes in id order with ties to the smallest label, so it is deterministic. `layout` is a spring
layout whose repulsion is summed over a `LAYOUT_GRID` 12x12 grid of mass centres rather than over
every pair, so it is linear in the genes: 2,000 genes in about 0.7 s, cached per layout key and
never run from a paint handler.

`_NetworkCanvas` builds one `QPainterPath` per (source, strength bucket) once, in scene coordinates,
and paints them through the current transform into a cached pixmap. Pan, zoom and hover repaint that
pixmap plus a small overlay. Minimap, fit, click-to-recentre and double-click-to-3D-map all work.

## 3. Defining the connections

`Definition` (Qt-free, hashable, JSON round-tripping) with four rules: which sources may contribute,
`min_strength`, `max_per_gene` (strongest first, a link needs room at both ends) and `min_sources`
("at least N different sources say so"). `apply_definition` returns the kept links and a count of
what each rule removed; `counts_text` is the live headline. `DefinitionStore` keeps named
definitions in `star_edges/star_definitions.json` under `paths.user_cache_dir()`, written through a
temporary file. Every control carries a tooltip saying what it does and why.

## Standing constraints honoured

No strategy algorithm, calibration number or shipped data table changed. No new organism literals.
`app.py` is untouched — `install()` already builds the panel, and it now also hands it a
`DefinitionStore`. `umap_gallery*.py`, `glass.py` and `strategy_panel.py` are untouched.
