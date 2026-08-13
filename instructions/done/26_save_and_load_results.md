# Save and load results, per tab and all at once — DONE 2026-08-12

`starplast/results.py`, wired into every results table's right-click menu and into the File menu.

## The format, and why it is this one

    one table   a CSV whose FIRST LINE is `# starplast-table: recovery_search`
    every tab   a zip of those CSVs plus a manifest

A zip of CSVs rather than a pickle: these are results. Someone will want to open one in a
spreadsheet, or read it on a machine that has never had this program installed, and a format that
needs the program to read it loses the results the day the program stops building.

## The marker line is the point

A row in these tables is a **recipe** — blocks, missing-value policy, scaling, both UMAP
hyperparameters, the clustering size, the seed, the subsample, the excluded columns — which is what
makes clicking it rebuild a map. So a loaded table has to be as clickable as a computed one, and a
test proves exactly that: save a search, empty the table, load the file, click row 0, and the rebuild
receives the same blocks.

Which is why the file names its table. A validation file dropped into the walk tab would produce rows
nobody can rebuild and an error message about the wrong thing, so a mismatch is **refused**, naming
the tab it came from. A file with no marker — someone's own CSV, or one saved before this existed —
loads, and the status line says it did not say where it came from.

## Small decisions

- **Empty tables are left out of a bundle.** One listing eight tables of which six are empty invites
  the reader to think six analyses returned nothing, when they were never run.
- **One unreadable member costs that table, not the bundle.** A bundle is a session's work.
- **A bundle from a later version names the tables this one cannot show.** "3 of 5 loaded" with no
  names is not something anyone can act on.

## Verified

`results.py` at 100%, `analysis_panel.py` at 99%; 1,755 tests pass.
