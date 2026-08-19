#!/usr/bin/env python3
"""Histone H4 acetylation (ChIP-chip, via ToxoDB)

Genome-wide H4 K5/K8/K12/K16 acetylation score within 1 kb of the gene

    level / kind : DNA / ChIPchip
    provides     : h4_acetylation_chip_score
    coverage     : 7,515 genes (92%)
    accession    : ToxoDB Hakimi/Ali H4 acetylation
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByChIPchiptgonME49_chipChipExper_Hakimi_Ali_RSRC/reports/attributesTabular
    local path   : starplast/data/toxodb_h4_acetylation.tsv

Quirks that cost time once:
    Verified as an ACTIVE mark must behave: rho +0.36 with transcription, +0.43 with promoter
    ATAC, -0.02 with fitness, and genes in its top decile are expressed four times as highly as
    those in its bottom. The Einstein H3K4me1 report from the same site, the same assay type and
    the same query shape does the opposite -- its marked genes have LESS accessible promoters --
    and is in quarantine. This entry is the counter-example that says that refusal is about the
    data and not about the reader.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_h4_acetylation.py
"""
from _common import run

KEY = "toxodb_h4_acetylation"

if __name__ == "__main__":
    run(KEY)
