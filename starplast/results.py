#!/usr/bin/env python3
"""Saving what the analysis tabs computed, and loading it back so it is still clickable.

A search is minutes to hours and every tab's table is lost when the window closes. What makes that
worth fixing here rather than by exporting a CSV is what a row in these tables IS: not a summary but
a recipe -- blocks, missing-value policy, scaling, both UMAP hyperparameters, the clustering size,
the seed, the subsample and the excluded columns. That is why clicking a row rebuilds its map.

**So a loaded table has to be as clickable as a computed one.** Load a search from last week, click
a row, and the map rebuilds. A loader that produced a table you could only read would be a
screenshot with extra steps.

That requires one thing the plain CSV export does not do: recording WHICH table a file came from. A
validation table dropped into the walk tab makes rows nobody can rebuild and an error message about
the wrong thing, so the file names its table and loading refuses a mismatch instead of guessing.

Two shapes, the same format underneath:

    one table   a CSV whose first line is `# starplast-table: recovery_search`
    every tab   a zip of those CSVs plus a manifest saying what is in it

A zip of CSVs rather than a pickle or an HDF5: these are results, someone will want to open one in a
spreadsheet or read it on a machine that has never had this program installed, and a format that
requires the program to read it is a format that loses the results the day the program stops
building.
"""
from __future__ import annotations

import io
import json
import os
import zipfile

import pandas as pd

from .logging_util import get_logger

_log = get_logger(__name__)

#: The marker line. Read back by `table_kind`, and deliberately a comment so that everything else
#: can open the file as an ordinary CSV.
MARKER = "# starplast-table:"

#: What a bundle carries besides the tables: the version that wrote it and when. Enough to explain a
#: file that will not load in two years, and not so much that it pretends to be provenance -- the
#: provenance is in the rows.
MANIFEST = "manifest.json"


def save_table(path: str, kind: str, df: pd.DataFrame) -> str:
    """Write one table, naming which table it is on the first line."""
    with open(path, "w", encoding="utf8", newline="") as fh:
        fh.write(f"{MARKER} {kind}\n")
        df.to_csv(fh, index=False)
    _log.info("results: wrote %d rows of %s to %s", len(df), kind, path)
    return path


def table_kind(path: str) -> str:
    """Which table a file holds, or "" if it does not say.

    A file that does not say is still loadable -- somebody's own CSV, or one saved before this
    existed -- but the caller is the one that decides whether to accept it, because only the caller
    knows which table it is being dropped into.
    """
    try:
        with open(path, encoding="utf8") as fh:
            first = fh.readline().strip()
    except (OSError, UnicodeDecodeError):
        # A file that cannot be read as text names no table -- which is the honest answer, and lets
        # the caller report "could not read this" rather than dying on the first line.
        return ""
    return first[len(MARKER):].strip() if first.startswith(MARKER) else ""


def load_table(path: str) -> pd.DataFrame:
    """Read a saved table, with or without the marker line."""
    skip = 1 if table_kind(path) else 0
    return pd.read_csv(path, skiprows=skip)


def save_bundle(path: str, tables: dict, meta: dict | None = None) -> str:
    """Write every tab's results into one file.

    Empty tables are skipped rather than written as headers: a bundle that lists eight tables of
    which six are empty invites the reader to think six analyses returned nothing, when they were
    never run.
    """
    from . import __version__
    kept = {k: df for k, df in tables.items() if df is not None and len(df)}
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(MANIFEST, json.dumps({"starplast": __version__,
                                         "tables": {k: int(len(df)) for k, df in kept.items()},
                                         **(meta or {})}, indent=1))
        for kind, df in kept.items():
            buf = io.StringIO()
            buf.write(f"{MARKER} {kind}\n")
            df.to_csv(buf, index=False)
            z.writestr(f"{kind}.csv", buf.getvalue())
    _log.info("results: wrote %d tables to %s", len(kept), path)
    return path


def load_bundle(path: str) -> tuple:
    """Read a bundle: `(tables, manifest)`.

    A member that will not parse costs that table and says so, rather than the whole file: a bundle
    is a session's work and losing seven tables to one bad row would be its own kind of failure.
    """
    tables, meta = {}, {}
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            if name == MANIFEST:
                try:
                    meta = json.loads(z.read(name).decode("utf8"))
                except ValueError:
                    meta = {}
                continue
            if not name.endswith(".csv"):
                continue
            try:
                text = z.read(name).decode("utf8")
                head, _, body = text.partition("\n")
                kind = head[len(MARKER):].strip() if head.startswith(MARKER) else name[:-4]
                tables[kind] = pd.read_csv(io.StringIO(body if head.startswith(MARKER) else text))
            except Exception as exc:
                _log.warning("results: %s in %s could not be read: %s: %s",
                             name, path, type(exc).__name__, exc)
    return tables, meta


def is_bundle(path: str) -> bool:
    """Whether a path is a bundle rather than a single table."""
    return zipfile.is_zipfile(path) if os.path.exists(path) else path.endswith(".starplast")
