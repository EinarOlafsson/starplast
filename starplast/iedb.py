#!/usr/bin/env python3
"""Antibody epitopes per gene, from IEDB directly rather than through ToxoDB.

ToxoDB integrates IEDB and exposes a count per gene, which is what fills `T-cell epitope content`.
That count is not split by epitope type. IEDB's own API is, and the split is the whole point here:
`seroreactivity / antigenicity` is a question about ANTIBODIES, and answering it with a number
dominated by T-cell assays would be answering a different question with a plausible column.

`bcell_search` is the antibody half. Filtering it to Toxoplasma source antigens gives 939 assay
records over 44 antigens, which reduce to 241 distinct epitope sequences.

## Distinct epitopes, not assay records

A protein studied by twenty groups accumulates twenty records for the same peptide, so counting
records would rank antigens by how fashionable they are. Counting distinct sequences asks how much
of the protein antibodies actually recognise, which is what the slot means by antigenicity.

## Why the mapping is by product description

IEDB identifies antigens by UniProt accession and a name, and the names are verbatim ToxoDB product
descriptions -- `SAG-related sequence SRS29B`, `Dense granule protein GRA6`. Matching those exactly
resolves 34 of 44. Matching by the trailing symbol instead resolves only 24 and loses SRS29B, which
is SAG1 and the single most studied antigen in the organism, so the description is the better key and
the symbol is only a fallback. The eleven that remain are generic -- `Uncharacterized protein`,
`PA14 domain-containing protein` -- and match many genes each, so they are dropped rather than
guessed.
"""
from __future__ import annotations

import os

import pandas as pd

TABLE = "iedb_bcell_epitopes.tsv"
COLUMN = "n_bcell_epitopes"


def bcell_epitopes(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Distinct antibody epitope sequences per gene.

    Absent is absent, not zero: IEDB records what somebody tested, and a protein nobody has raised
    an antibody against is not a protein with no epitopes.
    """
    path = os.path.join(base, "starplast", "data", TABLE)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    if d.shape[1] < 2:
        return pd.DataFrame()
    genes = d[d.columns[0]].astype(str)
    if resolve is not None:
        genes = genes.map(lambda g: resolve(g) or g)
    values = pd.Series(pd.to_numeric(d[d.columns[1]], errors="coerce").to_numpy(),
                       index=genes.to_numpy(), dtype=float).groupby(level=0).max()
    log(f"IEDB antibody epitopes: {int(values.notna().sum()):,} genes, "
        f"{int(values.sum()):,} distinct epitopes")
    return values.to_frame(COLUMN)
