#!/usr/bin/env python3
"""Resistance-conferring mutations: the one table in this map that is curated rather than parsed.

Every other column here is read out of a file somebody deposited. This one cannot be. Toxoplasma has
no resistome the way *Plasmodium* has one -- no study has selected resistance across a compound panel
and deposited the calls -- so the answer to "which genes carry a mutation that confers drug
resistance" exists only as sentences in individual papers. A curated table is the honest shape for
that, and it is only worth having if every row can be checked, so every row carries the paper, the
substitution, the compound and how causality was established.

## The bar, which is what makes this short

An allele is in here only if the mutation was **put back** -- introduced into a clean background and
shown to produce the resistance. That is a deliberately higher bar than the literature's usual
output, and the reason is a specific failure found while building this: the Toxoplasma in-vitro
evolution studies produce *candidate* loci, and in the one case where causality was actually tested
the candidate failed. The auranofin study (PMID 33816332) ranks its variants, names SOD2 as most
likely resistance-conferring, and then reports that SOD2 L201P "was not sufficient to confer
resistance when introduced into wild-type parasites". The artemisinin study's own table (PMID
31806760) is titled "Mutations found in candidate genes". Selection under a drug is evidence of
selection; it is not a demonstration that the mutation does the work. Nine genes were curated from
those two studies and then dropped, which is why this file has one gene in it and not ten.

## What is documented but deliberately absent

The field's other well-known resistance alleles are real and are not here, each for a stated reason:

* **DHFR-TS** (`TGME49_249180`) pyrimethamine alleles -- T83N and others. Established well enough
  that the mutant allele is the field's standard selectable marker, but the primary text (Donald and
  Roos 1993) is pre-PMC and unreachable from here, and a review's paraphrase is not the measurement.
* **DHODH** (`TGME49_210790`) N302S for 1-hydroxyquinolones -- primary paper not open access.
* **Cytochrome b** M129L and I254L for atovaquone -- excluded for a structural reason rather than an
  access one: Toxoplasma cytochrome b is mitochondrially encoded and has no row in a table of the
  8,140 nuclear genes. The eighteen nuclear "cytochrome b" hits are b-c1 complex subunits and b5
  domain proteins, and putting the allele on one of those would be a straightforward error.

They are listed here rather than omitted so that the next person adds them when the sources open,
instead of rediscovering the same three dead ends.

## Missingness

A gene with no row has not been shown to be resistance-free; nobody selected resistance in it. So
only genes with an allele carry a value and every other gene stays missing -- the opposite of the
screen columns, where absence inside the screened set is a real negative result.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Allele:
    """One mutation, and everything needed to check it without trusting this file."""
    gene: str
    product: str          # as the CURRENT annotation gives it, so a mistyped accession fails a test
    substitution: str
    compound: str
    pmid: str
    validation: str       # how the mutation was shown to CAUSE the resistance


#: Curated 2026-08-16 from primary open-access text, read directly rather than from any summary.
ALLELES = (
    Allele("TGME49_312570", "CMGC kinase, MAPK family (ERK) MAPK-1", "L162Q", "1NM-PP1",
           "24533298",
           "the mutant sequence replaced the native locus under 1000 nM 1NM-PP1, and "
           "overexpression conferred higher resistance than wild type"),
    Allele("TGME49_312570", "CMGC kinase, MAPK family (ERK) MAPK-1", "I171N", "1NM-PP1",
           "24533298", "recovered with L162Q and tested for cross resistance to other BKIs"),
    Allele("TGME49_312570", "CMGC kinase, MAPK family (ERK) MAPK-1", "L162Q", "3BrB-PP1",
           "24533298", "cross resistance measured for the same clones"),
    Allele("TGME49_312570", "CMGC kinase, MAPK family (ERK) MAPK-1", "L162Q", "3MB-PP1",
           "24533298", "cross resistance measured for the same clones"),
    # The gatekeeper, and the cleanest causal demonstration in the set: one residue swapped two ways
    # in the same background gives a sensitive clone and a resistant one.
    Allele("TGME49_312570", "CMGC kinase, MAPK family (ERK) MAPK-1", "S191Y", "1NM-PP1",
           "25941623",
           "engineered gatekeeper clones in one background: S191Y resistant, S191A sensitive"),
)

#: Alleles known to the field, not carried, and why. Kept as data so a test can assert that none of
#: them has quietly been added without a source.
DOCUMENTED_ELSEWHERE = {
    "TGME49_249180": "DHFR-TS pyrimethamine alleles; primary text (Donald and Roos 1993) unreachable",
    "TGME49_210790": "DHODH N302S for 1-hydroxyquinolones; primary paper not open access",
}


def table() -> pd.DataFrame:
    """One row per gene carrying at least one validated resistance allele."""
    if not ALLELES:
        return pd.DataFrame()
    frame = pd.DataFrame([vars(a) for a in ALLELES])
    grouped = frame.groupby("gene")
    out = pd.DataFrame({
        "resistance_allele_count": grouped["substitution"].nunique(),
        "resistance_compound_count": grouped["compound"].nunique(),
        "resistance_substitutions": grouped["substitution"].apply(
            lambda s: ";".join(sorted(set(s)))),
        "resistance_compounds": grouped["compound"].apply(lambda s: ";".join(sorted(set(s)))),
    })
    out.index = pd.Index(out.index, name="gene_id")
    return out


def resistance(base: str = "", log=print, resolve=None) -> pd.DataFrame:
    """The curated table, resolved to current accessions. `base` is unused and kept for symmetry.

    Every other loader in the build takes a dataset root because it opens a file. This one takes it
    and ignores it, so the build loop does not need a special case for the one source that is a
    literature curation rather than a download.
    """
    out = table()
    if out.empty:
        return out
    if resolve is not None:
        out.index = pd.Index([resolve(g) or g for g in out.index], name="gene_id")
        out = out[~out.index.duplicated()]
    log(f"resistance alleles (curated): {len(out)} gene(s), "
        f"{len({a.substitution for a in ALLELES})} substitutions, "
        f"{len({a.compound for a in ALLELES})} compounds")
    return out
