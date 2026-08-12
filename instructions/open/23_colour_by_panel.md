# The left panel becomes "color by", and colors by anything

Requested 2026-08-12. The panel says **filter by**, which is half of what it does and the less
important half. It should say **color by**, and be the one place a user chooses what the map is
colored by -- because looking at structure colored by many different held-out variables is what this
application is for.

## What goes in it

- **Categorical columns**, as now.
- **Clusterings**, as they are generated. Each run appears with its name and its cluster count.
  `AnalysisPanel.clusters_ready` already emits the labels and `Window.use_clusters` already colors
  by them; what is missing is that a clustering is not kept, named, or listed -- so the second one
  replaces the first with no way back.
- **Numeric columns**, binned. The user chooses the number of bins. A continuous column colored by a
  ramp is already offered elsewhere; binning it makes it behave like a category, which is what makes
  it comparable with a clustering.

## Runs

- Auto-named by timestamp to the second on creation, so two runs a minute apart are distinguishable
  without anyone typing anything.
- Renameable, with name-and-save at the bottom of the panel.
- Saved with their full recipe. `EmbeddingStore` already does this for embeddings; clusterings need
  the same, or a run cannot be rebuilt and the name is a label on nothing.

## The thing to be careful about

A clustering of a subsample has fewer labels than the map has genes. `Window.colours` already refuses
to draw a mismatched clustering rather than lining labels up by position -- keep that. A run listed
in the panel must know which genes it applies to, not just how many.

## Done when

Choosing a clustering from the panel colors the map by it; two runs can both be kept and switched
between; a numeric column can be binned into a chosen number of bins and colored as a category; runs
carry a timestamp name that can be edited and saved.
