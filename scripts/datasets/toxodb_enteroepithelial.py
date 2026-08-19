#!/usr/bin/env python3
"""Enteroepithelial stage transcriptome (via ToxoDB)

Expression in the feline enteroepithelial stages against tachyzoites

    level / kind : transcription / RNAseq
    provides     : ees_vs_tachyzoite_log2
    coverage     : 7,739 genes (95%)
    accession    : ToxoDB Ramakrishnan enteroepithelial
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByRNASeqtgonME49_Ramakrishnan_enteroepithelial_stages_ebi_rnaSeq_RSRC/reports/attributesTabular
    local path   : starplast/data/toxodb_enteroepithelial.tsv

Quirks that cost time once:
    EES1-5 averaged against tachyzoites, sense strand, no fold-change floor. The only life stage
    in the map that happens inside a cat. Verified against stage markers: GRA11B, which is
    merozoite-specific, comes out at +8.99 log2 and the family A/B/C merozoite antigens at
    +2.66, against a genome median of -0.14. ToxoDB reports a SIGNED fold difference, so -1257.4
    means 1257-fold down and is converted rather than logged.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_enteroepithelial.py
"""
from _common import run

KEY = "toxodb_enteroepithelial"

if __name__ == "__main__":
    run(KEY)
