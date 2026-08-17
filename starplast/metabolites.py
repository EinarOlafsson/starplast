#!/usr/bin/env python3
"""The metabolite table: rows are compounds, not genes.

Three slots on the metabolism axis ask about metabolites -- their steady-state levels, the flux
through them, the composition of membrane lipids -- and none of them can be answered by a table whose
rows are genes. A metabolite is not a gene, and the mapping between them is a metabolic model rather
than a measurement. Those slots carry `unit="metabolite"` for that reason, and this is the table they
resolve against.

It is the same shape as the host tables instruction 39 describes: a second table, keyed by its own
identifier, reached through slots that declare their unit rather than through rows bolted onto the
parasite table.

## What identifies a row

The compound name as the study reported it, normalised for case and punctuation only. There is no
Toxoplasma-wide metabolite identifier the way there is for genes -- studies report names, sometimes
HMDB or KEGG accessions, often neither -- so joining two studies means matching names, and matching
names is lossy. With one study in the table that cost is not yet paid; the moment a second arrives it
is the first problem to solve, and it should be solved with an identifier and not with fuzzy matching.
"""
from __future__ import annotations

import io
import os
import re
import zipfile

import numpy as np
import pandas as pd

TABLE = "metabolites.parquet"

#: The iron-deprivation study, which measures all three things a metabolite slot can ask for in one
#: deposit: steady-state levels, and glucose and glutamine labelling.
SOURCE = "PMC13170339"
LEVELS = ("mbio.03788-25-s0004.xlsx", "data_normalised")
LABELLING = {"glucose": ("mbio.03788-25-s0006.xlsx", "Relative abundance glucose"),
             "glutamine": ("mbio.03788-25-s0006.xlsx", "Relative abundance glutamine")}


def norm(name) -> str | None:
    """A compound name reduced to what two studies could plausibly agree on."""
    if not isinstance(name, str):
        return None
    key = re.sub(r"[^a-z0-9]+", "", name.strip().lower())
    return key or None


def _levels(archive: zipfile.ZipFile) -> pd.DataFrame:
    """Steady-state abundance: the log2 fold change under iron deprivation, per compound."""
    member, sheet = LEVELS
    d = pd.read_excel(io.BytesIO(archive.read(member)), sheet_name=sheet, header=1)
    # The sheet carries two blocks side by side under one header row; the volcano block holds the
    # summary, and its columns are unnamed because the real names sit in the first data row.
    body = d.iloc[1:].copy()
    body.columns = list(d.iloc[0].fillna(pd.Series(d.columns, index=d.columns)))
    name = body.columns[0]
    fold = next((c for c in body.columns if isinstance(c, str) and "log2" in c.lower()), None)
    padj = next((c for c in body.columns
                 if isinstance(c, str) and "p.ajusted" in c.lower().replace(" ", "")), None)
    if fold is None:
        return pd.DataFrame()
    out = pd.DataFrame({
        "metabolite": body[name].astype(str),
        "metabolite_level_log2fc_iron_depleted": pd.to_numeric(body[fold], errors="coerce")})
    if padj is not None:
        out["metabolite_level_padj"] = pd.to_numeric(body[padj], errors="coerce")
    return out.dropna(subset=["metabolite_level_log2fc_iron_depleted"])


def _labelled_fraction(archive: zipfile.ZipFile, member: str, sheet: str) -> pd.Series:
    """How much of a compound's pool came from the labelled precursor, in control conditions.

    One minus the unlabelled isotopologue. That is the standard readout and the only one that
    survives compounds of different carbon number: M1 through M55 are not comparable between a
    three-carbon and a twenty-carbon molecule, and M0 is.
    """
    # One of the two sheets carries the table title in its first row and the other does not, so the
    # header row is found rather than assumed. Assuming it cost the glucose arm silently: the sheet
    # read fine, had no column called `Metabolite`, and returned empty.
    blob = io.BytesIO(archive.read(member))
    d = pd.read_excel(blob, sheet_name=sheet)
    if "Metabolite" not in d.columns:
        blob.seek(0)
        d = pd.read_excel(blob, sheet_name=sheet, header=1)
    if not {"Metabolite", "M0"} <= set(d.columns):
        return pd.Series(dtype=float)
    control = d[d.get("ExpCondition", "").astype(str).str.lower() == "control"]
    if control.empty:
        return pd.Series(dtype=float)
    m0 = pd.to_numeric(control["M0"], errors="coerce")
    frame = pd.DataFrame({"metabolite": control["Metabolite"].astype(str), "m0": m0}).dropna()
    return 1.0 - frame.groupby("metabolite")["m0"].mean()


def build(archive_path: str, lipid_archive: str | None = None) -> pd.DataFrame:
    """Assemble the metabolite table from the polar-metabolite study, plus lipid species if given.

    The two arrive as separate archives because they are separate studies measuring disjoint
    classes of compound; `lipids` explains why appending them by name is safe here and what would
    make it unsafe.
    """
    if not os.path.exists(archive_path):
        return pd.DataFrame()
    with zipfile.ZipFile(archive_path) as archive:
        names = set(archive.namelist())
        parts = []
        if LEVELS[0] in names:
            levels = _levels(archive)
            if not levels.empty:
                levels["key"] = levels["metabolite"].map(norm)
                parts.append(levels.dropna(subset=["key"]).set_index("key"))
        for precursor, (member, sheet) in LABELLING.items():
            if member not in names:
                continue
            got = _labelled_fraction(archive, member, sheet)
            if got.empty:
                continue
            frame = got.rename(f"labelled_fraction_{precursor}").to_frame()
            frame["metabolite"] = frame.index
            frame.index = pd.Index([norm(i) for i in frame.index], name="key")
            parts.append(frame[frame.index.notna()])
    if not parts:
        return pd.DataFrame()
    # Every part carries `metabolite`, so the join drops the duplicate and the name is never
    # missing. A guard for the missing case used to sit here and could not fire; unreachable
    # defensive code tells a reader the case is handled when nothing handles it.
    out = parts[0]
    for extra in parts[1:]:
        out = out.join(extra.drop(columns=["metabolite"]), how="outer")
    out = out[~out.index.duplicated()]
    out["metabolite"] = out["metabolite"].fillna(pd.Series(out.index, index=out.index))
    out = out.reset_index()
    if lipid_archive:
        from . import lipids
        extra = lipids.build(lipid_archive)
        if not extra.empty:
            # Appended, not joined: a lipid species and a polar metabolite are never the same row,
            # and a key that did collide would mean the naming assumption in `lipids` had broken.
            extra = extra[~extra["key"].isin(set(out["key"]))]
            out = pd.concat([out, extra], ignore_index=True)
    return out.drop(columns=["key"])


def load(base: str) -> pd.DataFrame:
    """The shipped metabolite table, or an empty frame if it has not been built."""
    path = os.path.join(base, "starplast", "data", TABLE)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
