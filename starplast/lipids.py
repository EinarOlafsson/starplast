#!/usr/bin/env python3
"""Membrane lipid composition: rows are lipid species, and the species are the parasite's own.

The third metabolite slot asks what the parasite's membranes are made of. A lipidome of an infected
culture cannot answer it. The host cell is most of the lipid in the dish, and a lipid measured there
is not thereby a parasite lipid -- which is why the slot sat empty while several Toxoplasma
lipidomics papers were on the shelf. This source answers it because the vesicles were collected from
tachyzoites AFTER egress, in host-cell-free medium, so what is in them the parasite carried out.

## Why these read as parasite composition rather than host carry-over

The study grew the parasite in four host backgrounds -- fibroblast, intestinal epithelium, myotube,
Vero -- whose own lipidomes are not remotely alike: between any two of them 1,018 to 1,362 species
differ at FDR < 0.05. Between the vesicles the same parasite released in those four backgrounds, 0 to
4 do. A composition set by the host would move when the host is replaced; this one does not. Both
numbers are in the shipped archive, and `tests/test_lipids.py` recomputes them rather than trusting
the sentence -- the whole case for calling these parasite lipids rests on that contrast, so it is
checked and not cited.

## Two columns, because there are two questions

`lipid_ev_level_log2` is how much of a species is present. It is the published processed matrix,
averaged over the nineteen vesicle samples, and nothing else.

`lipid_ev_vs_host_clr` is how much more of it is present than in the host cell the parasite grew in.
That question needs the loading difference removed first: vesicles carry far less total material than
a cell, so every raw EV-minus-cell difference is negative (-7 to -14 log2) and says only that a
vesicle is smaller than a cell. Centring each sample on its own mean before differencing asks about
proportion instead of amount, which is the question the word "composition" is asking.

The published archive has its own EV-minus-cell column. It is not used: it correlates 0.92 with every
reconstruction attempted here but matches none of them to better than 3 log units, so its transform
is not known. A number whose definition cannot be reproduced cannot be verified, and this column is
computed here from the published per-sample matrices instead, by the rule stated above.

## The key, and why the name collision the metabolite table warns about does not arise

Rows join the metabolite table on the same normalised name. `metabolites` warns that matching a
second study by compound NAME is lossy and is the first thing to fix when a second study arrives.
This is that second study, and the cost is not paid here: lipid shorthand (`PI 38:4`, `SM 34:1;O2`)
and polar metabolite common names (`citrate`, `glutamine`) are different naming systems that do not
overlap, and none of the 194 species collides with any of the 1,102 compounds already in the table.
That is asserted in the tests. A third study of polar metabolites will still need an identifier.
"""
from __future__ import annotations

import io
import os
import zipfile

import pandas as pd

from .metabolites import norm

#: EVs released by post-egress tachyzoites in host-cell-free medium, in four host backgrounds.
SOURCE = "PMC12913473"

#: Abundance: the processed, imputed, log2 vesicle matrix as published.
LEVELS = ("Table2.xlsx", "Processed_EV_Log2_Imp (2)")

#: Composition against the background it was drawn from. The two sheets are the intersected species,
#: so the vesicle and the host cell matrices carry the same rows and the difference is well defined.
PAIRED = ("Table3.xlsx", "Intersection_EV_Log2_Imp", "Intersection_Cell_Log2_Imp")

#: The four host backgrounds, as they appear inside the sample column names.
HOSTS = ("Fibro", "IPEC", "Myo", "Vero")

#: The two sheets that carry the conservation argument: how many species differ between backgrounds,
#: in the host cells and in the vesicles the parasite released in them.
CONSERVATION = (("Table1.xlsx", "Summary_Pairwise"), ("Table2.xlsx", "Summary"))


def _sheet(archive: zipfile.ZipFile, member: str, sheet: str) -> pd.DataFrame:
    """One sheet of one workbook inside the supplementary archive."""
    return pd.read_excel(io.BytesIO(archive.read(member)), sheet_name=sheet)


def _samples(frame: pd.DataFrame, host: str) -> pd.DataFrame:
    """The replicate columns for one host background.

    Matched on the host name appearing in the column, because the replicate numbering is not
    contiguous -- the vesicle arm is missing Myo_2, and a range would quietly read the wrong column.
    """
    return frame[[c for c in frame.columns if host in c]]


def _clr(frame: pd.DataFrame) -> pd.DataFrame:
    """Each sample centred on its own mean, so columns compare proportion and not amount."""
    return frame.sub(frame.mean(axis=0), axis=1)


def ev_levels(archive: zipfile.ZipFile) -> pd.Series:
    """Mean log2 abundance of each lipid species across the vesicle samples."""
    member, sheet = LEVELS
    d = _sheet(archive, member, sheet)
    if "Lipid_names" not in d.columns:
        return pd.Series(dtype=float)
    values = d.drop(columns=["Lipid_names"]).apply(pd.to_numeric, errors="coerce")
    return pd.Series(values.mean(axis=1).to_numpy(), index=d["Lipid_names"].astype(str))


def ev_vs_host(archive: zipfile.ZipFile) -> pd.Series:
    """COMPUTED here: composition of the vesicle against the host cell it was grown in.

    Per host background, centre both matrices sample-wise, average the replicates, and subtract.
    The four backgrounds are then averaged, which is the quantity the study's own conservation
    result licenses: if the vesicle composition does not depend on the host, one number stands for
    all four, and averaging them is not hiding a spread.
    """
    member, ev_sheet, cell_sheet = PAIRED
    ev, cell = _sheet(archive, member, ev_sheet), _sheet(archive, member, cell_sheet)
    if "Lipid_names" not in ev.columns or "Lipid_names" not in cell.columns:
        return pd.Series(dtype=float)
    ev = ev.set_index(ev["Lipid_names"].astype(str)).drop(columns=["Lipid_names"])
    cell = cell.set_index(cell["Lipid_names"].astype(str)).drop(columns=["Lipid_names"])
    ev, cell = ev.apply(pd.to_numeric, errors="coerce"), cell.apply(pd.to_numeric, errors="coerce")
    parts = []
    for host in HOSTS:
        here, there = _samples(ev, host), _samples(cell, host)
        if here.empty or there.empty:
            continue
        parts.append(_clr(here).mean(axis=1) - _clr(there).mean(axis=1))
    if not parts:
        return pd.Series(dtype=float)
    return pd.concat(parts, axis=1).mean(axis=1)


def conservation(archive: zipfile.ZipFile) -> tuple:
    """(host-cell contrasts, vesicle contrasts): how many species differ between backgrounds.

    The evidence that these are parasite lipids, read straight out of the archive so a test can
    assert on it. Returns the two count columns, one row per pair of host backgrounds.
    """
    out = []
    for member, sheet in CONSERVATION:
        d = _sheet(archive, member, sheet)
        col = next((c for c in d.columns if "signif" in str(c).lower()), None)
        out.append(pd.to_numeric(d[col], errors="coerce").dropna() if col else pd.Series(dtype=float))
    return tuple(out)


def build(archive_path: str) -> pd.DataFrame:
    """Lipid-species rows for the metabolite table, keyed the way that table is keyed."""
    if not os.path.exists(archive_path):
        return pd.DataFrame()
    with zipfile.ZipFile(archive_path) as archive:
        names = set(archive.namelist())
        if LEVELS[0] not in names:
            return pd.DataFrame()
        levels = ev_levels(archive)
        against = ev_vs_host(archive) if PAIRED[0] in names else pd.Series(dtype=float)
    if levels.empty:
        return pd.DataFrame()
    out = pd.DataFrame({"metabolite": levels.index.astype(str),
                        "lipid_ev_level_log2": levels.to_numpy()})
    out["lipid_ev_vs_host_clr"] = out["metabolite"].map(against)
    out["key"] = out["metabolite"].map(norm)
    out = out.dropna(subset=["key"])
    return out[~out["key"].duplicated()].reset_index(drop=True)
