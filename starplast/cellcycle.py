#!/usr/bin/env python3
"""Cell-cycle phase, and stage enrichment — kept rigorously apart, because one is measured and one is not.

Half of what this project exists to do could not be evaluated before this module. The stated aim is to
walk dataset combinations and find embeddings where **location or cell-cycle stage** maps onto the
structure, so that a gene landing in that structure inherits a testable prediction. Location has a
target: 3,827 hyperLOPIT compartment labels to hold out. Cell-cycle stage had none — searching all 165
node columns for anything cyclical returned only continuous expression (`expr_tachy`, tachyzoite FPKM,
oocyst iTRAQ ratios). There was nothing discrete to hold out, so the question could not be asked, let
alone answered.

## The measured label

Xue et al. 2020 (eLife 9:e54129, PMID 32065584) sequenced single parasites and assigned each gene a
cell-cycle phase from where its expression peaks along pseudotime. Supplementary file 3 gives 964 genes
a phase from **G1a, G1b, S, M, C**. That is an independent measurement of cell-cycle behaviour: it comes
from single-cell expression this project does not otherwise use, and it is a label a paper committed to,
not a threshold chosen here.

    cellcycle_phase             G1a / G1b / S / M / C     964 genes, measured
    cellcycle_pseudotime        1 / 2 / 3                 8,590 genes, coarser

Those files use **TGGT1_** accessions — the RH strain, type I — while this project is built on ME49.
The Pru differentiation files in the same paper use TGME49_. Two strains inside one supplementary set,
told apart only by the accession prefix; routing through the identity layer is mandatory here, not a
precaution.

## The derived label, and why it is named that way

`stage_enriched_derived` assigns each gene the life-cycle stage its own expression is highest in. It is
useful and it is **not evidence**: it is a restatement of the expression columns it was computed from.

Holding it out against an embedding built on those same columns would be circular by construction, and
this project has already been bitten by exactly that. `compartment` fed one embedding and the cluster
battery reported its twin `lopit_map` and its derivations as the top held-out discoveries at V = 0.96;
held out honestly, the best of 164 runs reached mean F1 0.362. So the derived column carries `_derived`
in its name, the registry records which columns produced it, and `search.py`'s circularity guard —
which measures association with the inputs rather than trusting names — will mark it derived on its own.

The naming is the cheap part. The guard is what actually protects the result.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

# Where the study's files sit under the dataset root, and which column of each carries the label.
STUDY = os.path.join("transcription", "scRNAseq", "32065584")
PHASE_FILE = "cellcycle_phase_RH.csv"
PSEUDOTIME_FILE = "pseudotime_cluster_RH.csv"

PHASES = ("G1a", "G1b", "S", "M", "C")


def _read(path: str) -> pd.DataFrame:
    """Tab-separated despite the .csv extension, with an unnamed first column holding the accession.

    Read as comma-separated it parses as a single column of 5,691 strings and every downstream join
    silently produces nothing — the failure looks like "this dataset covers no genes", not like a
    parse error.
    """
    d = pd.read_csv(path, sep="\t")
    return d.rename(columns={d.columns[0]: "gene_id"})


def _clean(series: pd.Series) -> pd.Series:
    """Strip the layered quoting the export carries: a phase arrives as '\"\"\"C\"\"\"'."""
    return series.astype(str).str.replace('"', "", regex=False).str.strip()


def phase_labels(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Measured cell-cycle phase per gene, on canonical accessions.

    `resolve` is the identity layer. Without it every TGGT1_ accession fails to join and the result is
    an empty frame rather than an error, which is the failure mode this project keeps meeting.
    """
    from . import paths

    p = paths.find(STUDY, PHASE_FILE)
    if not p:
        log(f"cell cycle: {PHASE_FILE} not found under any dataset root")
        return pd.DataFrame()

    d = _read(p)
    d["cellcycle_phase"] = _clean(d["cell_cycle_phase"])
    d = d[d.cellcycle_phase.isin(PHASES)]

    ids = d.gene_id.astype(str).str.strip()
    if resolve is not None:
        mapped = ids.map(lambda a: resolve(a) or None)
        lost = int(mapped.isna().sum())
        if lost:
            log(f"cell cycle: {lost} of {len(ids)} accessions did not resolve")
        d = d.assign(gene_id=mapped).dropna(subset=["gene_id"])

    # One row per gene. A gene appearing under two probes with different phases is a genuine
    # ambiguity, so it is dropped rather than resolved by taking whichever sorted first.
    counts = d.groupby("gene_id").cellcycle_phase.nunique()
    ambiguous = set(counts[counts > 1].index)
    if ambiguous:
        log(f"cell cycle: {len(ambiguous)} genes carry conflicting phases and are left unlabelled")
    d = d[~d.gene_id.isin(ambiguous)].drop_duplicates("gene_id").set_index("gene_id")

    out = d[["cellcycle_phase"]]
    log(f"cell cycle: {len(out)} genes with a measured phase "
        f"({', '.join(f'{k} {v}' for k, v in out.cellcycle_phase.value_counts().items())})")
    return out


def pseudotime_clusters(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Coarse pseudotime cluster for the whole RH gene set. Broader coverage, weaker label."""
    from . import paths

    p = paths.find(STUDY, PSEUDOTIME_FILE)
    if not p:
        return pd.DataFrame()
    d = _read(p)
    d["cellcycle_pseudotime"] = pd.to_numeric(_clean(d["pseudotime_cluster"]), errors="coerce")
    ids = d.gene_id.astype(str).str.strip()
    if resolve is not None:
        d = d.assign(gene_id=ids.map(lambda a: resolve(a) or None)).dropna(subset=["gene_id"])
    d = d.dropna(subset=["cellcycle_pseudotime"]).drop_duplicates("gene_id").set_index("gene_id")
    log(f"cell cycle: {len(d)} genes with a pseudotime cluster")
    return d[["cellcycle_pseudotime"]]


# --------------------------------------------------------------------------- derived stage classes
# Which node columns speak for which life-cycle stage. Named explicitly rather than pattern-matched:
# a regex over column names would silently absorb any new column whose name happened to match, and the
# derivation's provenance is the whole reason this column is allowed to exist at all.
#
# One column per stage, all three from the same family, deliberately:
#
#   * The oocyst iTRAQ columns are RATIOS between isobaric labels (115:113, 116:114). A ratio is a
#     comparison between two channels, not an abundance, so putting one on the same axis as an
#     abundance and taking an argmax compares two different kinds of quantity. They also cover only
#     2,068 genes against 7,739.
#   * Giving one stage three columns and another one column makes the three-column stage's mean less
#     noisy than its rivals'. The comparison is then partly about how much data each stage happens to
#     have, which is study effort, not biology.
#
# expr_tachy / expr_cyst / expr_sporulated are one measurement family on one scale covering the same
# genes, which is the only version of this comparison that means what it appears to mean.
STAGE_COLUMNS = {
    "tachyzoite": ("expr_tachy",),
    "bradyzoite": ("expr_cyst",),
    "oocyst": ("expr_sporulated",),
}

MIN_MARGIN = 0.5     # z-units the winning stage must lead by before the call is made


def stage_enrichment(nodes: pd.DataFrame, log=print) -> pd.DataFrame:
    """Assign each gene the stage its expression is highest in. DERIVED, never evidence.

    Every contributing column is z-scored first, because they are not on one scale — FPKM, iTRAQ ratios
    and log intensities cannot be compared as raw numbers, and whichever column happened to have the
    largest units would otherwise win every gene.

    A gene is left unlabelled unless one stage leads the next by `MIN_MARGIN`. Assigning every gene a
    class would manufacture confident labels for the flat majority, and a label that is really a coin
    toss is worse than an absent one: it looks like a measurement in every table it appears in.

    Read the class counts with one caveat in mind. Tachyzoite and bradyzoite expression correlate at
    r = 0.89 while sporulated oocyst correlates with both at ~0.30, so a gene high in tachyzoite is
    usually also high in bradyzoite and neither wins by a margin, whereas oocyst wins easily. The
    result is 1,380 oocyst calls against 230 tachyzoite -- which reflects how separable the stages are
    from each other, not how many genes each stage uses. The margin rule is doing its job here; it is
    the interpretation that has to stay careful.
    """
    present = {stage: [c for c in cols if c in nodes.columns] for stage, cols in STAGE_COLUMNS.items()}
    usable = {s: c for s, c in present.items() if c}
    if len(usable) < 2:
        log(f"stage enrichment: only {list(usable)} available; not enough to compare")
        return pd.DataFrame(index=nodes.index)

    scores = {}
    for stage, cols in usable.items():
        block = nodes[cols].apply(pd.to_numeric, errors="coerce")
        z = (block - block.mean()) / block.std(ddof=0).replace(0, np.nan)
        scores[stage] = z.mean(axis=1)
    S = pd.DataFrame(scores)

    # Only genes measured in at least two stages can be compared at all; idxmax on an all-NaN row is
    # meaningless (and deprecated), so those rows are excluded before the argmax rather than after.
    comparable = S.notna().sum(axis=1).ge(2)
    best = pd.Series(pd.NA, index=S.index, dtype=object)
    best[comparable] = S[comparable].idxmax(axis=1)
    top2 = np.sort(S.to_numpy(), axis=1)[:, ::-1][:, :2]
    margin = pd.Series(top2[:, 0] - top2[:, 1], index=S.index)

    call = best.where(comparable & (margin >= MIN_MARGIN))
    out = pd.DataFrame({"stage_enriched_derived": call,
                        "stage_margin_derived": margin.where(call.notna())}, index=nodes.index)
    n = int(call.notna().sum())
    log(f"stage enrichment (DERIVED): {n} of {len(nodes)} genes called from "
        f"{sum(len(c) for c in usable.values())} columns across {len(usable)} stages; "
        f"{', '.join(f'{k} {v}' for k, v in call.value_counts().items())}")
    return out


def add_all(base: str, nodes: pd.DataFrame, resolve=None, log=print) -> pd.DataFrame:
    """Attach the measured cell-cycle columns and the derived stage column to the node table."""
    n = nodes
    for fn in (phase_labels, pseudotime_clusters):
        t = fn(base, resolve=resolve, log=log)
        if not t.empty:
            for c in t.columns:
                n[c] = n.gene_id.map(t[c])
    derived = stage_enrichment(n, log=log)
    for c in derived.columns:
        n[c] = derived[c].values
    return n
