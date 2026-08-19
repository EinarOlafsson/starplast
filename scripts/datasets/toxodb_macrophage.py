#!/usr/bin/env python3
"""Expression in infected macrophages (via ToxoDB)

Expression percentile in ME49-infected murine macrophages

    level / kind : transcription / RNAseq
    provides     : macrophage_expression_percentile
    coverage     : 8,140 genes (100%)
    accession    : ToxoDB Saeij 29 strains
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByRNASeqtgonME49_Saeij_Jeroen_strains_rnaSeq_RSRCPercentile/reports/attributesTabular
    local path   : starplast/data/toxodb_macrophage.tsv

Quirks that cost time once:
    The ME49 arm of a 29-strain panel, so the strain matches the rest of the map. Correlates at
    rho +0.89 with the fibroblast transcriptome (`expr_tachy`). That is a FINDING -- the
    parasite's transcriptional programme is largely independent of which host cell it is in --
    and not a construction: it is an independent measurement in a different host context, with
    nothing to declare in derived_from. It is written down here so that nobody counts the two as
    independent evidence when they agree, which they mostly will.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_macrophage.py
"""
from _common import run

KEY = "toxodb_macrophage"

if __name__ == "__main__":
    run(KEY)
