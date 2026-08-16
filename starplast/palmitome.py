#!/usr/bin/env python3
"""S-palmitoylation, from the palmitome ToxoDB serves for Foe et al. 2015.

The paper's supplementary tables are behind PMC's download interstitial and the article is not open
access, so the numbers come from ToxoDB's own query service for the same dataset -- which returns
the authors' fold differences per gene rather than a re-analysis.

## The two comparisons, and why only one of them is palmitoylation

Both label proteins with 17-ODYA, an alkyne palmitate analogue that click chemistry then pulls down.

* **against hydroxylamine** -- hydroxylamine cleaves thioester bonds, which is the bond an
  S-palmitoyl group makes. A protein enriched here lost its label to hydroxylamine, so the label was
  on a cysteine thioester. This is the S-palmitoylation-specific comparison.
* **against palmitic acid** -- competition with unlabelled palmitate. It shows the label is
  fatty-acid-dependent, but it does not distinguish a thioester from an amide, so it includes
  N-myristoylated and otherwise lipidated proteins.

Both are shipped because they answer different questions and the difference between them is
informative. The slot is filled by the hydroxylamine column.

## What the sign means

ToxoDB reports a signed fold DIFFERENCE, not a ratio: -3.12 means three-fold down, not a negative
quantity. Reading it as a ratio and taking its logarithm would produce NaN for every depleted protein
and quietly drop half the table, so it is converted explicitly.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

#: Which downloaded table fills which column. The file names are written by the fetch below; the
#: comparison each one holds is stated in the column name, because "palmitome" alone does not say
#: whether a number is thioester-specific.
COMPARISONS = {
    "toxodb_palmitome_hydroxylamine.tsv": "palmitome_odya_vs_hydroxylamine_log2",
    "toxodb_palmitome_palmitate.tsv": "palmitome_odya_vs_palmitate_log2",
}

#: Shipped inside the package rather than left in the dataset archive: two files of a few
#: hundred rows, and shipping them is what lets a fresh checkout rebuild the columns instead
#: of having to re-run a ToxoDB query whose parameters would then live only in a note.
FOLDER = os.path.join("starplast", "data")


def signed_log2(fold) -> float:
    """A signed fold difference as a log2 ratio. -4 becomes -2, not NaN."""
    value = pd.to_numeric(fold, errors="coerce")
    if pd.isna(value) or value == 0:
        return float("nan")
    return float(np.log2(value)) if value > 0 else float(-np.log2(-value))


def palmitome(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Both palmitome comparisons as log2 columns, indexed by resolved gene id.

    Absent genes are absent rather than zero. This is an enrichment experiment: a protein the
    pulldown never saw is not a protein with no palmitoylation, and 7,600 zeros would say it was.
    """
    out = pd.DataFrame()
    for filename, column in COMPARISONS.items():
        path = os.path.join(base, FOLDER, filename)
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t")
        if d.shape[1] < 2:
            continue
        genes = d[d.columns[0]].astype(str)
        if resolve is not None:
            genes = genes.map(lambda g: resolve(g) or g)
        values = pd.Series([signed_log2(v) for v in d[d.columns[1]]], index=genes.to_numpy())
        # A gene reported twice is one protein measured twice; the stronger enrichment is the one
        # the study makes a claim about, and averaging a hit with a non-hit would erase it.
        series = values.groupby(level=0).max()
        out = out.join(series.to_frame(column), how="outer") if len(out) else series.to_frame(column)
        log(f"palmitome: {column}, {int(series.notna().sum()):,} genes")
    return out
