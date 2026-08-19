#!/usr/bin/env python3
"""Codon usage bias, computed from the coding sequences and nothing else.

Three numbers per gene, all sequence-derived, which is what the slot asks for:

* **ENC** (Wright's effective number of codons) -- 20 when a gene uses one codon per amino acid, 61
  when it uses all of them evenly. Reference-free: it is computed from the gene's own codon
  homozygosity, so it needs no notion of which genes are highly expressed and cannot smuggle
  expression into a sequence column.
* **GC3** -- GC content at synonymous third positions. The compositional half of the same story, and
  the thing ENC is usually plotted against.
* **CAI** -- codon adaptation index against ribosomal proteins as the reference set. Ribosomal
  proteins are the standard choice because they are translated heavily in every organism, and using
  them rather than this map's own expression columns is what keeps the column from being a
  restatement of `expr_tachy`. That choice is the one judgement call here and it is declared in the
  registry as well as stated here.

## Why this is not circular, and where it nearly was

CAI's reference set is normally "the most highly expressed genes", which for this map would mean
building a sequence column out of an expression column and then being surprised when they correlate.
The ribosomal set is chosen by product annotation instead. The correlation with expression that
remains -- and there is one, which is the whole reason CAI is interesting -- is then a finding rather
than an artefact of construction.
"""
from __future__ import annotations

import gzip
import math
import os
from collections import defaultdict

import numpy as np
import pandas as pd

TABLE = "toxodb_cds.tsv.gz"

#: The standard genetic code, as codon -> amino acid. Stop codons map to `*` and are excluded from
#: every statistic here: a stop is not a synonymous choice.
BASES = "TCAG"
AMINO = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODE = {a + b + c: AMINO[i]
        for i, (a, b, c) in enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)}

#: Amino acids with exactly one codon have no synonymous choice to make, so they carry no
#: information about bias and are excluded from ENC and CAI alike.
SYNONYMS = defaultdict(list)
for _codon, _aa in CODE.items():
    if _aa != "*":
        SYNONYMS[_aa].append(_codon)
DEGENERATE = {aa: codons for aa, codons in SYNONYMS.items() if len(codons) > 1}


def translate(cds: str) -> str:
    """A coding sequence to its protein, in frame 1, with the terminator trimmed.

    Here because a phosphosite table can key on nothing this project's identity index carries -- the
    CDPK1 supplement reports numeric ids from a 2017 annotation -- while giving the 15-residue window
    around each site. A window is an identifier if it occurs in exactly one protein, and checking
    that requires the proteome the CDS table already ships. Anything that is not a clean triplet of
    ACGT becomes `X`, so a match is never made across a gap that was silently closed.
    """
    seq = str(cds).strip().upper().replace("U", "T")
    return "".join(CODE.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - 2, 3)).rstrip("*")


def codon_counts(cds: str) -> dict:
    """Codons in frame, skipping anything that is not a clean triplet of ACGT."""
    seq = cds.strip().upper().replace("U", "T")
    counts: dict = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in CODE:
            counts[codon] = counts.get(codon, 0) + 1
    return counts


def gc3(counts: dict) -> float:
    """GC at the third position of degenerate codons. NaN when there are none to look at."""
    gc = at = 0
    for codon, n in counts.items():
        if CODE[codon] == "*" or len(SYNONYMS[CODE[codon]]) == 1:
            continue
        if codon[2] in "GC":
            gc += n
        else:
            at += n
    return gc / (gc + at) if gc + at else float("nan")


def enc(counts: dict) -> float:
    """Wright's effective number of codons, from per-family homozygosity.

    Families with too few observations to estimate homozygosity are dropped and their class average
    is used in their place -- Wright's own prescription. A gene with no usable family at all gets
    NaN rather than 61, because "we could not tell" is not "unbiased".
    """
    by_class: dict = defaultdict(list)
    for aa, codons in DEGENERATE.items():
        observed = [counts.get(c, 0) for c in codons]
        n = sum(observed)
        if n < 2:
            continue
        homozygosity = sum((x / n) ** 2 for x in observed)
        f = (n * homozygosity - 1) / (n - 1)
        if f > 0:
            by_class[len(codons)].append(f)
    # 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6 -- the counts are how many families of each size the code has.
    total, sizes = 2.0, {2: 9, 3: 1, 4: 5, 6: 3}
    for size, weight in sizes.items():
        values = by_class.get(size)
        if not values:
            # Wright substitutes the F3 estimate when isoleucine is unusable; more generally, a class
            # with nothing observed is filled from the nearest class that has something, and a gene
            # with nothing at all anywhere is refused below.
            values = next((by_class[s] for s in (4, 2, 6, 3) if by_class.get(s)), None)
            if not values:
                return float("nan")
        total += weight / (sum(values) / len(values))
    return min(total, 61.0)


def reference_weights(counts_by_gene: dict) -> dict:
    """Relative adaptiveness per codon, from a pooled reference set of genes."""
    pooled: dict = defaultdict(int)
    for counts in counts_by_gene.values():
        for codon, n in counts.items():
            pooled[codon] += n
    weights = {}
    for aa, codons in DEGENERATE.items():
        best = max(pooled.get(c, 0) for c in codons)
        if not best:
            continue
        for codon in codons:
            # A codon the reference set never uses would give log(0); Sharp and Li's convention is a
            # small non-zero value, and 0.01 is the usual one.
            weights[codon] = max(pooled.get(codon, 0) / best, 0.01)
    return weights


def cai(counts: dict, weights: dict) -> float:
    """Geometric mean of relative adaptiveness over a gene's degenerate codons."""
    total, n = 0.0, 0
    for codon, k in counts.items():
        w = weights.get(codon)
        if w:
            total += k * math.log(w)
            n += k
    return math.exp(total / n) if n else float("nan")


def codon_usage(base: str, resolve=None, reference=None, log=print,
                table: str | None = None, nodes: str | None = None) -> pd.DataFrame:
    """ENC, GC3 and CAI per gene.

    `reference` is the gene ids of the CAI reference set. Left None it is taken from the product
    descriptions in the shipped node table -- the ribosomal proteins -- and if that table is not
    present the CAI column is left out rather than computed against an arbitrary set.

    `table` and `nodes` name the CDS file and the node table to read, so the Plasmodium arm can use
    this construction instead of writing its own. Sharing it is the point: ENC, GC3 and a CAI against
    ribosomal proteins are definitions, and two arms computing them differently would make a
    difference between the arms unreadable.
    """
    path = os.path.join(base, "starplast", "data", table or TABLE)
    if not os.path.exists(path):
        return pd.DataFrame()
    with gzip.open(path, "rt", errors="replace") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    counts_by_gene = {}
    for row in rows[1:]:
        if len(row) < 2 or not row[1]:
            continue
        gene = (resolve(row[0]) or row[0]) if resolve is not None else row[0]
        counts = codon_counts(row[1])
        if counts:
            counts_by_gene[gene] = counts
    if not counts_by_gene:
        return pd.DataFrame()

    out = pd.DataFrame(index=pd.Index(sorted(counts_by_gene), name="gene_id"))
    out["codon_enc"] = [enc(counts_by_gene[g]) for g in out.index]
    out["codon_gc3"] = [gc3(counts_by_gene[g]) for g in out.index]

    if reference is None:
        reference = _ribosomal(base, nodes or "nodes.parquet")
    chosen = {g: counts_by_gene[g] for g in (reference or ()) if g in counts_by_gene}
    if chosen:
        weights = reference_weights(chosen)
        out["codon_cai_ribosomal"] = [cai(counts_by_gene[g], weights) for g in out.index]
        log(f"codon usage: {out.shape[1]} columns, {len(out):,} genes, "
            f"CAI against {len(chosen)} ribosomal proteins")
    else:
        log(f"codon usage: {out.shape[1]} columns, {len(out):,} genes, no CAI reference set")
    return out


def _ribosomal(base: str, nodes: str = "nodes.parquet") -> list:
    """Ribosomal protein gene ids from the shipped node table, for the CAI reference set."""
    path = os.path.join(base, "starplast", "data", nodes)
    if not os.path.exists(path):
        return []
    n = pd.read_parquet(path, columns=["gene_id", "product"])
    hit = n["product"].astype(str).str.contains("ribosomal protein", case=False, na=False)
    return n.loc[hit, "gene_id"].astype(str).tolist()
