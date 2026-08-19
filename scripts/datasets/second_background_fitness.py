#!/usr/bin/env python3
"""Fitness in the reporter strain (COMPUTED)

Guide depletion over eight passages of ordinary growth

    level / kind : DNA / CRISPR_screen
    provides     : crispr_reporter_strain_p8_log2
    coverage     : 262 genes
    accession    : GSE132237
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132237/suppl/GSE132237_RAW.tar
    local path   : datasets/quarantine/2026_08_16_pride/Tg/essentiality_in_a_second_background/

Quirks that cost time once:
    The SAME archive as the differentiation screen, answering its other question. Its passage
    arms -- p8 against the input library -- are ordinary tachyzoite growth, and the reporter
    line is not the type I RH the genome-wide screens use, which is what makes this a second
    background rather than a repeat. A targeted library of nucleic-acid binding genes, so 262
    genes and not the genome. Verified by agreeing with the RH screen where it should: rho =
    +0.62 against fit_invitro_hff over 130 shared genes, close enough that the direction and the
    join are right and far enough that it is not a copy. The deposit sat in a folder named for
    this slot all day while only its differentiation arms were read.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/second_background_fitness.py
"""
from _common import run

KEY = "second_background_fitness"

if __name__ == "__main__":
    run(KEY)
