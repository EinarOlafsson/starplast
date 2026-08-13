# Save and load results, per tab and all at once

Requested 2026-08-12.

Every tab computes something expensive and shows it in a table: the hyperparameter walk, the
clustering walk, the held-out battery and its per-category view, the recovery search and its
per-category view, the validation scores and their candidate list. A search is minutes to hours. All
of it is lost when the window closes.

## What to build

- **Per tab**: save that tab's results to a file, and load one back into it.
- **All tabs at once**: one file holding everything, loaded back in one action.
- Reachable from the File menu for the whole set, and from each results table's right-click menu for
  that table.

## What makes a saved file worth having

A row in these tables already carries its own recipe — blocks, policy, scaling, seed, subsample,
excluded columns, cluster size. That is what makes a row clickable back into a map. **So a loaded
table must be as clickable as a computed one**: load a search result from last week, click a row,
and the map rebuilds. If loading produces a table you can only read, the feature is a screenshot.

Which means the file has to record WHICH table it came from, or a validation table loaded into the
walk tab produces rows that cannot be rebuilt and an error that names the wrong thing.

## Done when

A search can be saved, the window closed, the file loaded, and a row clicked to rebuild its map;
and one file can carry every tab's results at once.
