# Import data, with the preprocessing offered rather than assumed — DONE 2026-08-12

`starplast/importer.py` and `File ▸ Import data…` (Ctrl+I).

## Reads what this project already reads

CSV, TSV, `.txt`, Excel and parquet. The separator is sniffed rather than assumed — half the
supplementary tables in this field are `.txt` that are really TSV, and guessing wrongly gives one
column with everything in it, which looks like a file with no data rather than a file read badly.
Workbooks list their sheets, because published supplements routinely put the table on sheet 3 behind
a legend, and re-judge the numbers when the sheet changes.

Identifiers go through `tuning.import_table` and the identity layer, as instructed — including the
repair regexes, so `TgME49.208830` and `gene|TGME49_208830|v2` resolve instead of matching nothing,
and previous and strain accessions map forward. That is the failure that cost a published screen
every one of its rows.

## THE RANGE DECIDES, NOT THE FILENAME

The dialog leads with a table of each column's range, and the guess at what the numbers are is read
off those numbers with its reason attached: *"values reach 16,520, so this has NOT been logged"* is
something a person can agree or disagree with; "fpkm" alone is not. The same GEO series is the
example both ways — its FPKM file reaches 16,520, and the columns derived from it in the node table
are still called `_FPKM` and stop at 9.7 because they were logged upstream.

Writing that guess turned up a bug in it: a column running 0.001 to 0.05 was being called a logged
intensity, when it is a probability or a p-value, and it would then have been centred per sample as
though the columns were instrument runs. A logged intensity spans several units.

## Every choice, offered and recorded

Quantification (`sources.TRANSFORM`), scaling (`embedding.SCALINGS`), missing values
(`embedding.NA_POLICIES`), duplicate handling (mean/median/max/min/first), and a sign flip with the
warning that pooling screens without rank-normalising first is where the 64x spread bites. Each one
uses the implementation the rest of the program uses, not a second copy of it.

`indicator` is deliberately not handled at import: the embedding adds missingness as its own
weighted block, and doing it here as well would count absence twice.

Every import returns a **record** — the file, the identifier column, any repair, and each choice —
kept in `Window.imports` and logged. An imported column whose provenance is a memory of which
dropdowns were set is a column nobody can defend three weeks later.

## Never confused with a measurement

Columns are prefixed (`imported_`) and joined **in memory**: the cache on disk is untouched, and a
test asserts the node table's columns are unchanged after an import. Genes absent from the file are
missing, not zero. The Data tab gains an `imported (N columns)` block that feeds the map like any
other — and, like any other, cannot afterwards be evidence about the clusters.

## Verified

`importer.py` at 100%; 1,803 tests pass, including a round trip that imports a table, builds the
feature matrix from the imported block, and colours the map by one of its columns.
