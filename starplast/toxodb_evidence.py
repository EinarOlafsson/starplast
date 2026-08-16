#!/usr/bin/env python3
"""Per-gene evidence that ToxoDB integrates and nobody else publishes as a table.

ToxoDB does not only host datasets; it runs them through its own pipelines and exposes the result as
searches whose report can return the per-gene value. That is the route used here, and it is the only
route for some of these: the epitope mapping is ToxoDB's join of IEDB against the ME49 proteome, and
the ChIP-on-chip scores were published as array data that nobody has since re-tabulated per gene.

Every column is the value the search returns, unmodified except for the sign convention below. What
is chosen here is which comparison to ask for, and that choice is recorded per source.

## The sign convention, again

Like the palmitome, the RNA-seq searches report a signed fold DIFFERENCE: -1257.4 means
1257-fold down, not a negative expression. It is converted with the same `signed_log2`, for the same
reason -- passing it to a logarithm would turn every down-regulated gene into NaN.
"""
from __future__ import annotations

import os

import pandas as pd

from .palmitome import signed_log2

#: Each downloaded report, the column it becomes, and whether its value is a signed fold difference.
#: `query` records what was asked for, because a search that can be asked five ways produces five
#: different columns and the note beside the data must say which one this is.
#: The H4 acetylation ChIP is the counter-example that makes the refusal below firm. Same site, same
#: assay type, same pipeline, same query shape -- and it behaves the way an active mark must:
#: rho +0.36 with expression, +0.43 with promoter ATAC, -0.02 with fitness, and genes in its top
#: decile are expressed four times as highly as those in its bottom. So the H3K4me1 result is not an
#: artefact of how these reports are read here.
#:
#: NOT here, and deliberately: the Einstein H3K4me1 ChIP-on-chip. Its 231 "marked" genes have LESS
#: accessible promoters than unmarked genes (-0.51 against +0.39, Mann-Whitney p = 3e-50) and its
#: score correlates with promoter ATAC at rho = -0.37. H3K4me1 marks active and poised chromatin, so
#: that is backwards, and no reading of the assay makes it right. Whether the fault is in the
#: peak-to-gene assignment, in the array's coverage, or in what the dataset actually contains is not
#: something this project can settle -- and an unexplained backwards correlation is exactly the shape
#: of a column that would look like data and rank genes wrongly. The report is in
#: `datasets/quarantine/2026_08_16_toxodb/`, and `chromatin state - histone marks` stays empty.
SOURCES = (
    ("toxodb_epitopes.tsv", "iedb_epitope_count", False,
     "GenesWithEpitopes, organism=Toxoplasma gondii ME49, confidence High+Medium+Low"),
    ("toxodb_h4_acetylation.tsv", "h4_acetylation_chip_score", False,
     "GenesByChIPchip Hakimi/Ali genome-wide H4 K5-K8-K12-K16 acetylation, within 1 kb, no floor"),
    ("toxodb_macrophage.tsv", "macrophage_expression_percentile", False,
     "GenesByRNASeq Saeij 29 strains, ME49-infected murine macrophages, percentile, channel 1"),
    ("toxodb_enteroepithelial.tsv", "ees_vs_tachyzoite_log2", True,
     "GenesByRNASeq Ramakrishnan enteroepithelial, EES1-5 against tachyzoites, sense strand"),
)


def evidence(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Every downloaded ToxoDB search report as one column each.

    The value column is taken by POSITION -- the second -- rather than by name. ToxoDB's tabular
    report ships display names as the header (`Epitope Count`, `Score`, `Fold Change (Avg)`), and
    those change with the site's release while the report's shape does not.
    """
    out = pd.DataFrame()
    for filename, column, is_fold, _query in SOURCES:
        path = os.path.join(base, "starplast", "data", filename)
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t")
        if d.shape[1] < 2:
            continue
        genes = d[d.columns[0]].astype(str)
        if resolve is not None:
            genes = genes.map(lambda g: resolve(g) or g)
        raw = d[d.columns[1]]
        values = pd.Series([signed_log2(v) for v in raw] if is_fold
                           else pd.to_numeric(raw, errors="coerce").to_numpy(),
                           index=genes.to_numpy(), dtype=float)
        # Strongest rather than mean, for the same reason the palmitome takes the maximum: two rows
        # for one gene are one gene measured twice, and averaging a hit with a non-hit erases it.
        series = values.groupby(level=0).max()
        out = out.join(series.to_frame(column), how="outer") if len(out) else series.to_frame(column)
        log(f"toxodb evidence: {column}, {int(series.notna().sum()):,} genes")
    return out
