#!/usr/bin/env python3
"""The arrayed imaging screen: what a parasite looks like when a gene is switched off.

Every other fitness measurement in this map is a growth rate. A gene whose loss stops the parasite
growing gets a number, and the number says nothing about WHY. This screen is the other kind: 320
genes disrupted one per well, imaged, and each well given a phenotype by eye -- so a gene can fail
at egress specifically rather than merely fail.

## The category codes, which had to be worked out rather than read

The table gives each gene one or more codes -- `E3`, `F1 A2`, `R2` -- and the archive carries no
legend; it lives in a figure that ships as an image. Guessing what a letter means is exactly how a
column ends up labelled with the wrong biology, so the mapping here is verified twice over:

* The paper names its phenotype categories in three independent places and names the same four each
  time: replication, apicoplast, F-actin, and egress. Four categories, four letters, one initial
  each, one-to-one.
* Then the check that does not depend on initials at all. The paper names exactly two genes as the
  egress mutants it went on to characterise -- CGP (`TGGT1_240380`) and SLF (`TGGT1_208420`) -- and
  both of them carry an `E` code in this table and nothing else does the work. If `E` meant anything
  other than egress, those two rows would not read that way. `tests/test_phenotype_screen.py`
  asserts it, so the mapping is a claim the suite can lose.

The SUBSCRIPT is not interpreted. `E3` and `E4` differ in something the paper does not define in any
text available here -- severity, penetrance, or which replicate -- and a severity column invented out
of a digit would be a number with no measurement behind it. Only the presence of a category is read.

## Why absence means two different things, and why that had to be encoded

320 genes were screened and 99 got a phenotype. A gene in the library with no `E` was looked at and
did not have an egress defect; a gene outside the library was never looked at. Those are opposite
statements and both would arrive as an empty cell. So the screened set carries `False` and everything
else stays missing -- without that, "not tested" would read as "tested and normal" across 7,800 genes.

## What this does not cover

The screen's own figure legend calls it a screen for actin dynamics, apicoplast segregation and
egress. Invasion appears in the abstract as a property of hits characterised afterwards, not as a
category anything was scored into. So this fills the egress half of its slot and the detail says so;
a per-gene invasion phenotype for Toxoplasma is still not published.
"""
from __future__ import annotations

import os
import re

import pandas as pd

#: Nature Microbiology 2022 supplementary tables, PMID 35538310.
SOURCE = "41564_2022_1114_MOESM4_ESM.xlsx"

#: The 320 genes put into the library, and the 99 that came out with a phenotype.
LIBRARY = ("Table 2 library gRNAs", 2)
CALLS = ("Table 3B", 4)

#: Letter -> category, established as the docstring describes.
CATEGORIES = {"E": "egress", "F": "actin", "A": "apicoplast", "R": "replication"}

#: The two investigators scored independently, in their own columns.
SCORERS = ("inv1", "inv2")

#: Subscripts are written as unicode digits in some cells and ASCII in others.
_SUBSCRIPT = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def _find(base: str, name: str) -> str:
    for root, _dirs, files in os.walk(os.path.join(base, "datasets")):
        if name in files:
            return os.path.join(root, name)
    return os.path.join(base, "datasets", name)


def categories(cell) -> set:
    """Every phenotype category named in one investigator's cell.

    A cell holds one clone's calls, or several separated by `/` when the gene was picked more than
    once. The letter is what is read; the digit after it is deliberately ignored.
    """
    if not isinstance(cell, str):
        return set()
    text = cell.translate(_SUBSCRIPT)
    return {CATEGORIES[letter] for letter in re.findall(r"([EFAR])\s*\d", text)
            if letter in CATEGORIES}


def _accession(value) -> str | None:
    """`201270` and `TGGT1_201270` are the same gene written two ways in one workbook."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    match = re.search(r"(?:TGGT1_)?(\d{6})", text)
    return f"TGGT1_{match.group(1)}" if match else None


def screen(base: str, log=print, resolve=None) -> pd.DataFrame:
    """One row per screened gene, with a boolean per phenotype category."""
    path = _find(base, SOURCE)
    if not os.path.exists(path):
        log("phenotype screen: source not present")
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if LIBRARY[0] not in book.sheet_names or CALLS[0] not in book.sheet_names:
        log("phenotype screen: the archive does not carry both sheets")
        return pd.DataFrame()
    library = book.parse(LIBRARY[0], header=LIBRARY[1])
    calls = book.parse(CALLS[0], header=CALLS[1])
    calls.columns = (list(calls.columns[:1]) + ["name", "picked"]
                     + list(SCORERS[:max(0, len(calls.columns) - 3)]))
    screened = [a for a in (_accession(v) for v in library.iloc[:, 0]) if a]
    if not screened:
        log("phenotype screen: no genes in the library sheet")
        return pd.DataFrame()
    called = {}
    for _index, row in calls.iterrows():
        gene = _accession(row.iloc[0])
        if gene is None:
            continue
        seen = [categories(row.get(column)) for column in SCORERS if column in calls.columns]
        called[gene] = seen
    out = pd.DataFrame({"gene_id": sorted(set(screened))})
    for category in sorted(set(CATEGORIES.values())):
        out[f"screen_{category}_phenotype"] = [
            any(category in s for s in called.get(g, [])) for g in out["gene_id"]]
    out["screen_any_phenotype"] = [bool(called.get(g)) and any(called[g]) for g in out["gene_id"]]
    # Agreement is only defined where both investigators scored the same clone at all.
    out["screen_scorers_agree"] = [
        (called[g][0] == called[g][1]) if (g in called and len(called[g]) == 2
                                           and called[g][0] and called[g][1]) else None
        for g in out["gene_id"]]
    if resolve is not None:
        out["gene_id"] = [resolve(g) or g for g in out["gene_id"]]
    out = out[~out["gene_id"].duplicated()]
    hits = int(out["screen_any_phenotype"].sum())
    log(f"phenotype screen (PMID 35538310): {len(out):,} genes screened, {hits} with a phenotype, "
        f"{int(out['screen_egress_phenotype'].sum())} at egress")
    return out.set_index(pd.Index(out.pop("gene_id"), name="gene_id"))
