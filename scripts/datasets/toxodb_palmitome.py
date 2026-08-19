#!/usr/bin/env python3
"""S-palmitoylome (Foe 2015, via ToxoDB)

17-ODYA enrichment per gene, against hydroxylamine and against palmitate

    level / kind : post_translation / proteomics
    provides     : palmitome_odya_vs_hydroxylamine_log2, palmitome_odya_vs_palmitate_log2
    coverage     : 470 and 488 genes
    PMID         : 26468752
    accession    : ToxoDB Foe palmitome
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByProteomicsDirecttgonGT1_quantitativeMassSpec_Foe_Lipidome_Palmitoylome_RSRC/reports/attributesTabular
    local path   : starplast/data/toxodb_palmitome_hydroxylamine.tsv

Quirks that cost time once:
    The paper is not open access and PMC serves its supplementary spreadsheets only through a
    download interstitial, so the numbers come from ToxoDB's own query service for the same
    dataset -- the authors' fold differences, not a re-analysis. Two comparisons and only one is
    palmitoylation: hydroxylamine cleaves thioester bonds, which is the bond an S-palmitoyl
    group makes, so that column is thioester-specific; the palmitate competition shows only that
    the label is fatty-acid-dependent and includes N-myristoylated proteins. Verified against
    known substrates: ROP5 +2.20, GAP45 +1.34, AMA1 +0.93, MLC1 +0.71, IMC proteins +0.34,
    against a measured-gene median of -0.17. ToxoDB reports a SIGNED fold difference and not a
    ratio -- -3.12 means three-fold down -- so reading it as a ratio would have made every
    depleted protein NaN and dropped half the table.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_palmitome.py
"""
from _common import run

KEY = "toxodb_palmitome"

if __name__ == "__main__":
    run(KEY)
