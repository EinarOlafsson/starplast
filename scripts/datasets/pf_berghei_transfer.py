#!/usr/bin/env python3
"""P. berghei knockout fitness, transferred to falciparum

Relative growth of berghei knockouts, carried onto their falciparum orthologs

    level / kind : DNA / CRISPR_screen
    provides     : pb_transferred_phenotype, pb_transferred_growth_rate, pb_transfer_confidence
    coverage     : 2,448 falciparum genes of 2,578 berghei mutants
    PMID         : 28708996
    accession    : Cell 2017 Table S1 (PlasmoGEM)
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5509546/supplementaryFiles
    local path   : datasets/reference/plasmodb/pb_transfer/28708996/mmc1.xlsx

Quirks that cost time once:
    A TRANSFER, and it says so in every column name, because instruction 39 requires one to be
    visible rather than folded into the measured slot -- transfer berghei fitness onto
    falciparum, hold out falciparum fitness and recover it, and you have measured orthology.
    What makes this one safe is that the orthology is not ours: the screen's own table names a
    falciparum gene per row, and no falciparum gene is named by two berghei ones, so nothing is
    dropped for ambiguity and nothing is derived. Forty rows name a transcript rather than a
    gene and are stripped, the same suffix the phosphosite loader handles; kept whole they would
    have vanished for not matching an accession. CHECKED AGAINST THE RECEIVING ARM'S OWN SCREEN,
    which is the check a transfer has to pass: berghei-essential genes have a median piggyBac
    mutagenesis index of 0.160, slow ones 0.394 and dispensable ones 0.996 -- monotonic across
    two species and two unrelated methods, barcoded knockouts in mice against saturation
    mutagenesis in culture, p = 7e-107 -- and 65 of 71 ribosomal proteins come out essential.
    The 12 mutants the screen calls `Insufficient data` keep their confidence and lose their
    phenotype and growth rate: that phrase is the absence of a measurement, not a middle value.
    The other 3,272 falciparum genes are UNSCREENED, not dispensable, and stay missing.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_berghei_transfer.py
"""
from _common import run

KEY = "pf_berghei_transfer"

if __name__ == "__main__":
    run(KEY)
