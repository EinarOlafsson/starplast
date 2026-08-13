#!/usr/bin/env python3
"""Bringing a user's own table into the map, with every choice made along the way visible.

Everything here already exists somewhere in this project, applied to published data:
`tuning.import_table` resolves identifiers through the identity layer, `sources.normalize` applies
the transform a quantification type implies, `embedding` supplies the scalings and the missing-value
policies. What was missing is that none of it was reachable for a file of your own, so importing
meant either editing the node table by hand or not importing.

The one rule that governs the whole module, learned the expensive way and written on the wall of
`build_graph`:

    **THE RANGE DECIDES, NOT THE FILENAME.**

`GSE108740_FPKM.xlsx` reaches 16,520 and is real FPKM; the columns derived from it in the node table
are still called `rna108740_*_FPKM` and top out at 9.7, because they were logged upstream. Take the
log twice and real variation compresses to nothing; skip it and one gene dominates every distance.
Neither is visible in the output. So `describe_columns` reports each column's range and
`suggest_quantification` reads the numbers rather than the name -- and the interface shows both and
lets a person say which is right.

Nothing here writes to the node table on disk. An imported column lives in the session, is prefixed
so it cannot be mistaken for a measurement that shipped with the cache, and carries a record of
every choice that produced it.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from .embedding import NA_POLICIES, SCALINGS, _scale
from .logging_util import get_logger
from .sources import TRANSFORM
from .tuning import IMPORT_RX, suggest_gene_column, suggest_repair

_log = get_logger(__name__)

#: What can be read. Everything this project has already imported from: delimited text, Excel, and
#: parquet -- which is what the cache itself is written as, so a table exported from here comes back.
READABLE = (".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls", ".parquet")

#: How to combine several rows for one gene. `mean` by default because that is what `import_table`
#: has always done; `max` is here because a hit table often lists one row per guide or per peptide
#: and the strongest is the claim being made.
DUPLICATES = ("mean", "median", "max", "min", "first")

#: Quantification types, from `sources.TRANSFORM`, plus the honest option of leaving values alone.
QUANTIFICATIONS = ("none",) + tuple(TRANSFORM)


def read_any(path: str, sheet=0) -> pd.DataFrame:
    """Read a table in whatever format it arrived in.

    The separator is sniffed rather than assumed: half the supplementary tables in this field are
    `.txt` that are really TSV, and a `.csv` exported from Excel in a European locale is
    semicolon-separated. Guessing wrongly gives one column with everything in it, which looks like a
    file with no data rather than a file read badly.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=sheet)
    return pd.read_csv(path, sep=None, engine="python")


def sheets_in(path: str) -> list:
    """The sheet names of a workbook, or [] for anything else.

    Published supplements routinely put the table on sheet 3 behind a legend and a blank sheet, so
    which sheet is a choice the user has to be able to make.
    """
    if os.path.splitext(path)[1].lower() not in (".xlsx", ".xls"):
        return []
    try:
        return list(pd.ExcelFile(path).sheet_names)
    except Exception as exc:
        _log.warning("import: %s could not be opened as a workbook: %s: %s",
                     path, type(exc).__name__, exc)
        return []


def numeric_columns(df: pd.DataFrame, gene_column: str = "") -> list:
    """The columns worth importing: the numbers, minus the identifier column."""
    return [c for c in df.columns
            if c != gene_column and pd.api.types.is_numeric_dtype(df[c])]


def describe_columns(df: pd.DataFrame, columns=None) -> pd.DataFrame:
    """Each column's range, missingness and shape -- what the quantification decision rests on.

    Shown before anything is transformed, because this is the table on which a person decides
    whether they are looking at FPKM or at something that has already been logged. `max` is the
    number that gives it away: a column that reaches 16,520 has not been logged, one that stops at
    9.7 has.
    """
    columns = list(columns if columns is not None else numeric_columns(df))
    rows = []
    for c in columns:
        v = pd.to_numeric(df[c], errors="coerce")
        ok = v.dropna()
        rows.append({"column": c, "n": int(len(ok)),
                     "missing": float(v.isna().mean()) if len(v) else 1.0,
                     "min": float(ok.min()) if len(ok) else np.nan,
                     "median": float(ok.median()) if len(ok) else np.nan,
                     "max": float(ok.max()) if len(ok) else np.nan,
                     "negatives": float((ok < 0).mean()) if len(ok) else 0.0})
    return pd.DataFrame(rows, columns=["column", "n", "missing", "min", "median", "max",
                                       "negatives"])


def suggest_quantification(df: pd.DataFrame, columns=None) -> tuple:
    """A guess at the quantification type, and the reason for it, read off the NUMBERS.

    Never off the filename. The guess is a starting point for a human, which is why it comes back
    with its reason attached: "values reach 16,520, so this has not been logged" is something a
    person can agree or disagree with, and "fpkm" alone is not.
    """
    d = describe_columns(df, columns)
    if d.empty:
        return "none", "no numeric columns to judge"
    hi, lo = float(d["max"].max()), float(d["min"].min())
    negative = float(d["negatives"].mean())
    if negative > 0.2 and hi < 40:
        return "lfc", (f"{negative:.0%} of values are negative and the largest is {hi:.1f} -- that "
                       f"is a log fold change, already logged and already centred")
    if hi > 100:
        return "fpkm", (f"values reach {hi:,.0f}, so this has NOT been logged; if these are counts "
                        f"or TPM rather than FPKM the transform is the same")
    if lo >= 0 and 3 <= hi <= 40:
        # The floor matters: a column running 0.001 to 0.05 is a probability or a p-value, not a
        # logged intensity, and calling it one would centre it per sample as though the columns were
        # runs of an instrument. A logged intensity spans several units.
        return "log_intensity", (f"values run {lo:.1f} to {hi:.1f}, which is the range of something "
                                 f"already logged -- logging it again would compress real variation "
                                 f"to nothing")
    return "none", f"values run {lo:.3g} to {hi:.3g}; no transform is obviously right"


def preprocess(df: pd.DataFrame, gene_column: str = "", columns=None, quantification: str = "none",
               scaling: str = "none", na_policy: str = "median", duplicates: str = "mean",
               flip: bool = False, prefix: str = "imported_", resolve=None,
               repair: bool = True, log=print) -> tuple:
    """A user's table to gene-keyed, transformed columns -- returning the record of what was done.

    Returns `(frame, record)`. The record is not decoration: an imported column whose provenance is
    a memory of which dropdowns were set is a column nobody can defend three weeks later, and this
    project's whole discipline is that a number carries what produced it.

    The order is the order these things have to happen in. Identifiers first, because a table that
    resolves to nothing should say so before anything is computed. Then duplicates, because two rows
    for one gene must become one before a rank means anything. Then the sign, then the
    quantification transform, then the missing-value policy, then the scaling -- each one operating
    on what the previous produced.
    """
    from .tuning import import_table
    record = {"rows_in": int(len(df)), "quantification": quantification, "scaling": scaling,
              "na_policy": na_policy, "duplicates": duplicates, "flip": bool(flip),
              "prefix": prefix}
    column = gene_column or suggest_gene_column(df)[0]
    if column is None:
        raise ValueError("no column looks like a gene identifier -- choose one, or add a regex to "
                         "reshape it")
    pattern = replacement = None
    if repair and not df[column].astype(str).str.contains(IMPORT_RX, regex=True, na=False).any():
        # The repair regexes from `import_table`'s own suggestions: `TgME49.208830` and
        # `gene|TGME49_208830|v2` are both real formats from published supplements, and both match
        # nothing unrepaired.
        suggestion = suggest_repair(df[column])
        if suggestion:
            pattern, replacement = suggestion
            log(f"identifiers reshaped with {pattern!r}")
    record.update(gene_column=column, repair=pattern)

    keep = list(columns) if columns is not None else numeric_columns(df, column)
    if not keep:
        raise ValueError("no numeric columns to import from that table")
    wide = import_table(df[[column] + keep], gene_column=column, pattern=pattern,
                        replacement=replacement or r"\g<0>", resolve=resolve, prefix="", log=log)
    record["rows_resolved"] = int(len(wide))
    if not len(wide):
        raise ValueError("no row resolved to a gene id -- check the identifier column")

    # `import_table` already averages duplicates; anything else has to be done on the raw rows.
    if duplicates != "mean":
        ids = df[column].astype(str).str.extract(f"({IMPORT_RX.pattern})", expand=False, flags=2)
        ids = ids.str.upper()
        if resolve is not None:
            ids = ids.map(lambda x: resolve(x) if isinstance(x, str) else x)
        raw = df.loc[ids.notna(), keep].apply(pd.to_numeric, errors="coerce")
        raw.insert(0, "gene_id", ids[ids.notna()].values)
        wide = getattr(raw.groupby("gene_id"), duplicates)()
        log(f"duplicate rows combined by {duplicates}")

    if flip:
        # Screens disagree about which sign means "worse". Offered rather than guessed, and the
        # warning travels with it: pooling screens without rank-normalizing first is how the 64x
        # spread between them bites.
        wide = -wide
        log("sign flipped -- rank-normalize before pooling this with another screen")

    if quantification != "none":
        from .sources import normalize
        wide = normalize(wide, quantification, log=log)

    if na_policy == "drop_genes":
        before = len(wide)
        wide = wide.dropna()
        log(f"drop_genes: kept {len(wide):,} of {before:,} rows with no missing value")
    elif na_policy == "drop_columns":
        keep_cols = [c for c in wide.columns if wide[c].isna().mean() <= 0.5]
        log(f"drop_columns: kept {len(keep_cols)} of {wide.shape[1]} columns")
        wide = wide[keep_cols]
    elif na_policy == "median":
        wide = wide.fillna(wide.median(numeric_only=True))
    # `indicator` is deliberately NOT handled here: the embedding adds missingness indicators itself,
    # weighted below a measurement, and doing it twice would count absence twice.

    if scaling != "none":
        wide = pd.DataFrame(_scale(wide.to_numpy(dtype=float), scaling),
                            index=wide.index, columns=wide.columns)

    wide.columns = [f"{prefix}{c}" for c in wide.columns]
    record.update(columns=list(wide.columns), genes=int(len(wide)))
    log(f"imported {wide.shape[1]} column(s) for {len(wide):,} genes")
    _log.info("import: %s", record)
    return wide, record


def merge_into(nodes: pd.DataFrame, imported: pd.DataFrame) -> pd.DataFrame:
    """Attach imported columns to the node table, in memory, keyed by gene id.

    A copy, never the file. The cache is measurement that shipped with the program; an imported
    column is the user's own and must not become indistinguishable from it by being written into the
    same table on disk. The prefix does the same job for anyone reading a column name.
    """
    out = nodes.copy()
    joined = out[["gene_id"]].merge(imported, left_on="gene_id", right_index=True, how="left")
    for c in imported.columns:
        out[c] = joined[c].to_numpy()
    return out
