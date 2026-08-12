# The UMAP gallery: see every embedding, as it is generated

The prerequisite for everything else in the analysis workflow. Today "walk hyperparameters" returns
a table of numbers and no pictures, and "build this map" replaces the single central view, so there
is no way to look at 288 embeddings.

## What to build

- **Grid mode** — a wall of thumbnails, one per configuration, labelled with its parameters. Click
  one to expand it into the central view.
- **Scroll mode** — one large view with a slider stepping through configurations in order.
- **Both populate incrementally.** A UMAP appears the moment it is computed, not when the walk ends.
  This is the difference between a tool and a progress bar: a 288-run walk is half an hour, and the
  user needs to be looking at run 12 while run 13 computes.
- **An expanded UMAP behaves exactly like the main view** — click a gene for its evidence, colour by
  anything, draw edges. Not a picture of a map; a map.

## The hard part

`tuning.walk_umap` currently returns a table at the end. Incremental display requires it to emit per
configuration instead — a generator, or a callback. That contract change is the real work; the
widgets are straightforward once embeddings arrive one at a time.

Storage: the walk already writes each embedding through `EmbeddingStore`, so the gallery should read
from the store rather than holding 288 coordinate arrays in memory. Thumbnails can be rendered
offscreen with `grabFramebuffer` at low resolution.

## Done when

A walk of 20 configurations shows 20 thumbnails appearing one by one, any of which can be expanded
and clicked into, while the walk is still running.
