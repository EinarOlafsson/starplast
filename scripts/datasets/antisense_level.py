#!/usr/bin/env python3
"""Antisense transcription (via ToxoDB)

Percentile of antisense signal at this gene, across the life cycle

    level / kind : transcription / RNAseq
    provides     : antisense_expression_percentile
    coverage     : 8,140 genes (100%)
    accession    : ToxoDB full life-cycle transcriptome, Antisense
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByRNASeqtgonME49_tgme49_spor_ocyst_rnaseq_ebi_rnaSeq_RSRCPercentile/reports/attributesTabular
    local path   : starplast/data/antisense_level.tsv

Quirks that cost time once:
    The LEVEL of antisense transcription and not its change. An earlier attempt at this slot
    used the sense/antisense CHANGE across a tachyzoite time course and was refused because two
    strain time courses of it shared 12 of their top 200 genes where chance gives 28. The level
    reproduces: against the independent Gregory ME49 series it shares 197 of its top 500 where
    chance gives 31, a six-fold enrichment. How much antisense a gene has is a property of the
    gene; how much it changed was a property of the run. Correlates with sense transcription at
    only rho = +0.23, so it is not a restatement of expression.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/antisense_level.py
"""
from _common import run

KEY = "antisense_level"

if __name__ == "__main__":
    run(KEY)
