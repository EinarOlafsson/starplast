#!/usr/bin/env python3
"""Arrayed splitCas9 imaging screen

What a parasite looks like when a gene is off: egress, actin, apicoplast, replication

    level / kind : DNA / imaging_screen
    provides     : screen_egress_phenotype, screen_actin_phenotype, screen_apicoplast_phenotype, screen_replication_phenotype, screen_any_phenotype, screen_scorers_agree
    coverage     : 319 genes screened, 99 with a phenotype, 35 at egress
    PMID         : 35538310
    accession    : Nat Microbiol 41564-2022-01114 Supplementary Tables 2 and 3B
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-022-01114-y/MediaObjects/41564_2022_1114_MOESM4_ESM.xlsx
    local path   : datasets/DNA/imaging_screen/35538310/41564_2022_1114_MOESM4_ESM.xlsx

Quirks that cost time once:
    The only per-gene invasion-or-egress phenotype table published for Toxoplasma, and it covers
    the EGRESS half: the screen's own figure legend calls it a screen for actin dynamics,
    apicoplast segregation and egress, and invasion is a property its hits were shown to have
    afterwards rather than a category anything was scored into. Two earlier passes read this
    slot as blocked because the other candidates promise invasion and egress in their titles and
    deliver it by characterising one gene. The category codes had to be EARNED: the workbook
    ships no legend, it lives in a figure that is an image, so `E` is read as egress because the
    paper names exactly four categories in three places and because the two genes it names as
    its egress mutants, CGP TGGT1_240380 and SLF TGGT1_208420, both carry an E -- which the test
    asserts, so the mapping can fail. The SUBSCRIPT is deliberately not read: E3 and E4 differ
    in something no accessible text defines, and a severity invented from a digit is a number
    with no measurement behind it. Missingness carries the other half of the meaning -- a
    screened gene with no egress call was looked at and was normal, and the 7,800 unscreened
    genes stay missing, because collapsing those would tell the map that nearly every gene has
    been checked.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/splitcas9_imaging_screen.py
"""
from _common import run

KEY = "splitcas9_imaging_screen"

if __name__ == "__main__":
    run(KEY)
