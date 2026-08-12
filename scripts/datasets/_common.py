"""Shared machinery for the per-dataset scripts.

Each script in this directory covers one dataset in the registry: where it comes from, how it is
fetched, how its identifiers are standardised, and what it contributes. They exist so that a source
can be inspected, re-fetched or debugged on its own, without running a build that touches all 28.

**What these scripts are not.** They are not a second implementation of the build. The shipped
columns are assembled by `starplast.build_graph.load_nodes()` through the domain modules
(`localisation`, `expression`, `screens`, `cellcycle`), and a per-dataset copy of that logic would
drift from it and eventually lie about what produced the data. What these scripts do is the part
that is genuinely per-dataset: fetch the file, read it, resolve its accessions to current ToxoDB
ME49, and report the coverage that resolution achieves. Each script names the module that does the
authoritative normalisation for its columns.

The accession step is not a formality. Published supplements cite whatever identifier was current
when they were written, and `TGGT1_` turns out to be more common than `TGME49_`. A dataset keyed on
type I accessions joins zero rows without this, silently -- a join matching nothing looks exactly
like a dataset with no coverage.
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from starplast import datasets, identity, paths  # noqa: E402

_INDEX = None


def index():
    """The accession resolver, built once per process.

    Includes the GT1 and VEG strain tables. Leaving them out was a real bug: the cell-cycle table is
    964 rows of `TGGT1_` accessions and contributed nothing at all until they were added.
    """
    global _INDEX
    if _INDEX is None:
        data = paths.data_dir()
        gene_ids = pd.read_parquet(paths.cache_file("nodes.parquet")).gene_id
        _INDEX = identity.build_index(gene_ids, os.path.join(data, "toxodb_identity.tsv"),
                                      log=lambda *a: None)
        identity.add_strain_accessions(
            _INDEX, {"GT1": os.path.join(data, "toxodb_strain_gt1.tsv"),
                     "VEG": os.path.join(data, "toxodb_strain_veg.tsv")}, log=lambda *a: None)
    return _INDEX


def resolve(acc):
    """One accession to its current ME49 id, or None. Ambiguity is not guessed at."""
    hit = index().lookup.get(identity.norm(acc))
    return hit[0] if hit else None


def fetch(key: str, log=print) -> str | None:
    """Download the dataset if it is not already local. Returns the path, or None if unfetchable."""
    d = datasets.get(key)
    path = datasets.local_path(key)
    if path and os.path.exists(path):
        log(f"already local: {path}")
        return path
    ok, why = datasets.fetchable(key)
    if not ok:
        # Twelve of these have no downloadable URL. Saying so beats a stack trace, and beats
        # pretending the dataset is optional.
        log(f"cannot fetch automatically: {why}")
        log(f"  get it from: {d.url or d.citation}")
        log(f"  and put it at: {d.path}")
        return None
    return datasets.ensure(key, log=log)


def read(path: str, sep=None, sheet=None) -> pd.DataFrame:
    """Read a table, choosing the reader by extension rather than by the file's name.

    Extensions lie in this corpus -- the cell-cycle supplement is tab-separated with a `.csv`
    extension -- so `sep` is passed explicitly by the scripts that need it rather than sniffed.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet if sheet is not None else 0)
    if ext in (".parquet",):
        return pd.read_parquet(path)
    return pd.read_csv(path, sep=sep if sep is not None else ",", low_memory=False)


def id_column(df: pd.DataFrame) -> str | None:
    """The column most likely to hold gene accessions, by how many of its values resolve.

    Chosen by measurement rather than by name: these files label the column `Gene ID`, `gene_id`,
    `Accession`, `TGME49`, `Product` and several other things, and picking by name means a new
    spelling silently yields nothing.
    """
    best, best_hits = None, 0
    for c in df.columns:
        vals = df[c].dropna().astype(str).head(300)
        if vals.empty:
            continue
        hits = sum(resolve(v) is not None for v in vals)
        if hits > best_hits:
            best, best_hits = c, hits
    return best


def standardise(df: pd.DataFrame, col: str | None = None, log=print) -> pd.DataFrame:
    """Add a `gene_id` column of current ME49 accessions, and report what resolved.

    Rows that do not resolve are kept with a null `gene_id`. Dropping them would quietly shrink the
    dataset and make a broken identifier column look like a small one.
    """
    col = col or id_column(df)
    if col is None:
        log("no column in this file resolves to ToxoDB accessions")
        return df.assign(gene_id=pd.NA)
    out = df.copy()
    out["gene_id"] = [resolve(v) for v in out[col].astype(str)]
    n_ok = int(out.gene_id.notna().sum())
    log(f"identifier column {col!r}: {n_ok:,} of {len(out):,} rows resolved "
        f"({n_ok / max(len(out), 1) * 100:.0f}%), {out.gene_id.nunique()} distinct genes")
    return out


def describe(key: str, log=print) -> None:
    """Print the registry entry: the provenance, in the script that uses it."""
    d = datasets.get(key)
    log(f"{d.name}")
    log(f"  level/kind : {d.level} / {d.kind}")
    log(f"  provides   : {d.provides}")
    log(f"  columns    : {', '.join(d.columns) if d.columns else '(edges or inputs only)'}")
    log(f"  coverage   : {d.coverage}")
    log(f"  citation   : {d.citation}")
    if d.pmid:
        log(f"  PMID       : {d.pmid}")
    if d.accession:
        log(f"  accession  : {d.accession}")
    if d.url:
        log(f"  url        : {d.url}")
    log(f"  local path : {d.path}")
    if d.derived_from:
        log(f"  derived from: {', '.join(d.derived_from)}")
    if d.note:
        log(f"  note       : {d.note}")


def run(key: str, sep=None, sheet=None, id_col=None, normalised_by="build_graph.load_nodes()"):
    """Fetch, read, standardise and report one dataset. The body of every script here."""
    describe(key)
    print()
    path = fetch(key)
    if not path:
        return None
    df = read(path, sep=sep, sheet=sheet)
    print(f"read {len(df):,} rows x {len(df.columns)} columns from {os.path.basename(path)}")
    out = standardise(df, id_col)
    print()
    print(f"authoritative normalisation for the shipped columns: starplast.{normalised_by}")
    return out
