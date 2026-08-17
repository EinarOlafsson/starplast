#!/usr/bin/env python3
"""Fitness through the enteric and sexual cycle: does disrupting a gene cost the parasite oocysts.

The third curated table in the map, and it exists for a reason worth stating precisely. The sweeps
that closed this slot looked for a *pooled screen* through the enteroepithelial stages and correctly
found none -- the sexual cycle runs only in a felid, and nobody has put a barcoded library through a
cat. But the slot asks whether disrupting a gene changes fitness in the gut, and a pooled screen is
only one instrument that answers it. Feeding one knockout to a cat and counting the oocysts it sheds
answers the same question one gene at a time, and four labs have done it.

## The bar

A row exists only if a defined gene disruption was carried **through the feline stage** and the
enteric outcome was **measured** -- oocysts counted, or the sporulation rate scored. Read from
primary open-access text. A paper that says the cat experiment "should be carried out" does not
count, and one that says oocysts were seen but "the numbers were not quantified" does not either;
both were found and both were left out.

## Why the unchanged rows are the most valuable ones here

Four of the ten rows are `unchanged`, and they come from a deletion of all four LEA genes at once.
That is not a weaker result than a single knockout -- it is a stronger one. When a joint deletion of
a paralogue cluster leaves oocyst yield alone (30 against 34 million from paired kittens), no member
of the cluster is individually required, and redundancy cannot be hiding the phenotype. A single
knockout of one LEA gene could not have established that.

## What `background` is for, again

The AAH2 row measured in a delta-aah1 background is a genetic-interaction result, and the double
mutant is indistinguishable from delta-aah1 alone -- so reading a severe defect onto AAH2 from it
would be wrong. `background` says so per row, exactly as it does in `drug_sensitivity`.

## The graded outcomes are deliberate

`abolished`, `reduced` and `unchanged` are not one scale collapsed into three bins; they are what the
assays distinguish. HAP2 sheds a few mis-shapen oocysts that never sporulate, which is not the same
event as Grx5 shedding a third as many that sporulate poorly, and neither is AAH1 shedding four logs
fewer. The magnitudes live in `evidence` because they are not comparable across cats, strains and
inocula -- one numeric column here would invent a precision the experiments do not have.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class Enteric:
    """One gene carried through the enteric cycle, and everything needed to check it."""
    gene: str
    product: str          # as the CURRENT annotation gives it, so a mistyped accession fails a test
    outcome: str          # oocyst yield: abolished | reduced | unchanged
    sporulation: str      # abolished | reduced | unchanged | not measured
    background: str       # the strain the deletion was made in
    strain: str           # the parental line, which has to be oocyst-competent to be asked at all
    pmid: str
    evidence: str


OUTCOMES = ("abolished", "reduced", "unchanged")
SPORULATION = OUTCOMES + ("not measured",)

#: Curated 2026-08-17 from primary open-access text, read directly. Accessions resolved through the
#: ToxoDB identity index (`gene_name`), never by guessing from a product string -- which mattered
#: here, because ToxoDB spells both AAH products "hydrolase" where the papers say "hydroxylase".
MEASUREMENTS = (
    Enteric("TGME49_227100", "glutaredoxin 5", "reduced", "reduced", "parental", "Pru (type II)",
            "PMC12942651",
            "cats fed tissue cysts orally; total oocysts over the shedding period 2.5e7 in Pru "
            "against 8e6 in the deletion, ~70% down, and sporulation ~80% against ~30%. A "
            "complemented line was built, so the phenotype is tied to the locus"),
    Enteric("TGME49_287510", "aromatic amino acid hydrolase AAH1", "reduced", "not measured",
            "parental", "ME49 dhxg::Luc", "28288194",
            "~1e3 total oocysts per animal against 1e6-1e7 for wild type, in two of two cats. "
            "Sporulation could not be scored because the yield was too low to quantify"),
    Enteric("TGME49_212740", "aromatic amino acid hydrolase AAH2", "reduced", "reduced",
            "parental", "ME49 dhxg::Luc", "28288194",
            "milder and variable: much lower in two of three cats, the third only ~10-fold down to "
            "~1e5. Sporulation 75-80% for wild type against ~60%"),
    # Genetic interaction, not a clean single-gene result: the double is no worse than delta-aah1.
    Enteric("TGME49_212740", "aromatic amino acid hydrolase AAH2", "reduced", "not measured",
            "dAAH1", "ME49 dhxg::Luc", "28288194",
            "the double deletion's oocyst yield is severe (~1e3 per animal) but indistinguishable "
            "from delta-aah1 alone, so this row carries no evidence that AAH2 adds to the defect"),
    Enteric("TGME49_285940", "male gamete fusion factor HAP2, putative", "reduced", "abolished",
            "parental", "CZ clone H3", "30728393",
            "oocysts fell below reliable enumeration (<1,600 per two-day faecal sample against "
            "11 +/- 6 million for the parental line), were mis-shapen, and showed no diploidy by "
            "ploidy qPCR, so fertilisation did not occur. Sporulation was 0% at every timepoint "
            "against 45/69/82% for the parental. The authors report a second mutation in an intron "
            "of TGME49_223060 in this line and argue it is unlikely to contribute"),
    # Real negatives, and the strongest rows in the table: all four deleted together changed nothing.
    Enteric("TGME49_276850", "LEA domain-containing protein LEA850", "unchanged", "unchanged",
            "dLEA cluster (all four)", "ME49 dhpt luc+", "36809045",
            "paired kittens yielded 30 million oocysts from wild type against 34 million from the "
            "cluster deletion, with no morphological difference and no excystation or invasion "
            "defect"),
    Enteric("TGME49_276860", "LEA domain-containing protein LEA860", "unchanged", "unchanged",
            "dLEA cluster (all four)", "ME49 dhpt luc+", "36809045",
            "as LEA850: the four-gene deletion left oocyst yield alone"),
    Enteric("TGME49_276870", "LEA domain-containing protein LEA870", "unchanged", "unchanged",
            "dLEA cluster (all four)", "ME49 dhpt luc+", "36809045",
            "as LEA850: the four-gene deletion left oocyst yield alone"),
    Enteric("TGME49_276880", "LEA domain-containing protein LEA880", "unchanged", "unchanged",
            "dLEA cluster (all four)", "ME49 dhpt luc+", "36809045",
            "as LEA850. Sporulation was ~65% against ~85%, which the authors place inside the "
            "60-90% range they see normally, so it is recorded as unchanged rather than reduced"),
)


def table() -> pd.DataFrame:
    """One row per gene, summarising every enteric measurement made on it."""
    if not MEASUREMENTS:
        return pd.DataFrame()
    frame = pd.DataFrame([vars(m) for m in MEASUREMENTS])
    grouped = frame.groupby("gene")
    scored = frame[frame["sporulation"] != "not measured"].groupby("gene")
    out = pd.DataFrame({
        "enteric_oocyst_yield": grouped.apply(
            lambda g: ";".join(f"{b}:{o}" for b, o in zip(g["background"], g["outcome"])),
            include_groups=False),
        "enteric_sporulation": scored.apply(
            lambda g: ";".join(f"{b}:{s}" for b, s in zip(g["background"], g["sporulation"])),
            include_groups=False),
        "enteric_measurements": grouped.size(),
    })
    out["enteric_sporulation"] = out["enteric_sporulation"].fillna("not measured")
    out.index = pd.Index(out.index, name="gene_id")
    return out


def fitness(base: str = "", log=print, resolve=None) -> pd.DataFrame:
    """The curated table, resolved to current accessions. `base` is unused, kept for symmetry."""
    out = table()
    if out.empty:
        return out
    if resolve is not None:
        out.index = pd.Index([resolve(g) or g for g in out.index], name="gene_id")
        out = out[~out.index.duplicated()]
    log(f"enteric fitness (curated): {len(out)} genes over "
        f"{len({m.pmid for m in MEASUREMENTS})} studies, "
        f"{sum(1 for m in MEASUREMENTS if m.outcome == 'unchanged')} tested-and-unchanged")
    return out
