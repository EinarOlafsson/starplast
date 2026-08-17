#!/usr/bin/env python3
"""Drug sensitivity per gene: does disrupting it change survival under a compound.

The second curated table in the map, and it exists for the same reason `resistance` does -- no
genome-wide chemogenomic screen has been published for Toxoplasma, so the answer lives in individual
knockout studies rather than in a deposit. Six catalogue sweeps looked for the screen and correctly
did not find one; what they were testing was the sentence "no screen exists", which is not the same
statement as "the question cannot be answered". A knockout with a measured shift in sensitivity to a
named compound answers it one gene at a time.

## The bar

A row exists only if a knockout or knockdown showed a **measured** change in sensitivity to a
**named** compound, read from primary open-access text. Not "an inhibitor of this protein kills the
parasite" -- that is target engagement and has its own slot. The direction is recorded, and so is
`unchanged`, because a transporter deleted with no effect on analog sensitivity is a result and
suppressing it would leave the column looking like a list of hits.

## What the background qualifier is for

TgENT3's effect was measured in a ΔTgAT1 background rather than alone, so it is a genetic-interaction
result and not a clean single-gene one. `background` says so per row. Reading those two rows as though
the gene had been deleted on its own would be the same error as reading a double mutant's phenotype
off one of its genes.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Effect:
    """One gene against one compound, and everything needed to check it."""
    gene: str
    product: str          # as the CURRENT annotation gives it, so a mistyped accession fails a test
    compound: str
    direction: str        # resistant | sensitive | unchanged
    background: str       # the strain the deletion was made in
    pmid: str
    evidence: str


DIRECTIONS = ("resistant", "sensitive", "unchanged")

#: Curated 2026-08-17 from primary open-access text, read directly.
EFFECTS = (
    Effect("TGME49_244440", "adenosine transporter AT1", "Ara-A", "resistant", "parental",
           "41025776", "plaque and vacuole size across an Ara-A dose range in the single knockout"),
    Effect("TGME49_244440", "adenosine transporter AT1", "5-FU", "resistant", "parental",
           "41025776", "resistance to the pyrimidine analog, reported as unexpected for a purine "
                       "transporter"),
    # Measured in a delta-TgAT1 background, which makes these genetic-interaction results.
    Effect("TGME49_233130", "nucleoside transporter protein", "Ara-A", "resistant", "ΔTgAT1",
           "41025776", "the double knockout formed larger vacuoles than ΔTgAT1 alone under 10 µM"),
    Effect("TGME49_233130", "nucleoside transporter protein", "5-FU", "sensitive", "ΔTgAT1",
           "41025776", "heightened sensitivity in the same double knockout"),
    # Real negatives. A transporter deleted with no effect is a measurement.
    Effect("TGME49_500147", "nucleoside transporter, putative", "Ara-A", "unchanged", "parental",
           "41025776", "indistinguishable from parental across every phenotypic assay"),
    Effect("TGME49_500147", "nucleoside transporter, putative", "5-FU", "unchanged", "parental",
           "41025776", "indistinguishable from parental across every phenotypic assay"),
)


def table() -> pd.DataFrame:
    """One row per gene, summarising every compound it was tested against."""
    if not EFFECTS:
        return pd.DataFrame()
    frame = pd.DataFrame([vars(e) for e in EFFECTS])
    grouped = frame.groupby("gene")
    shifted = frame[frame["direction"] != "unchanged"].groupby("gene")
    out = pd.DataFrame({
        "drug_compounds_tested": grouped["compound"].nunique(),
        "drug_sensitivity_shifts": shifted["compound"].nunique(),
        "drug_sensitivity_directions": grouped.apply(
            lambda g: ";".join(f"{c}:{d}" for c, d in zip(g["compound"], g["direction"])),
            include_groups=False),
    })
    out["drug_sensitivity_shifts"] = out["drug_sensitivity_shifts"].fillna(0).astype(int)
    out.index = pd.Index(out.index, name="gene_id")
    return out


def sensitivity(base: str = "", log=print, resolve=None) -> pd.DataFrame:
    """The curated table, resolved to current accessions. `base` is unused, kept for symmetry."""
    out = table()
    if out.empty:
        return out
    if resolve is not None:
        out.index = pd.Index([resolve(g) or g for g in out.index], name="gene_id")
        out = out[~out.index.duplicated()]
    log(f"drug sensitivity (curated): {len(out)} genes over "
        f"{len({e.compound for e in EFFECTS})} compounds, "
        f"{sum(1 for e in EFFECTS if e.direction == 'unchanged')} tested-and-unchanged")
    return out
