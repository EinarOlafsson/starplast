#!/usr/bin/env python3
"""Protein melting temperature (mineCETSA)

Where this protein's melting curve sits, in degrees

    level / kind : post_translation / proteomics
    provides     : melting_temperature_tm
    coverage     : 3,120 proteins (38%)
    PMID         : 35976251
    accession    : eLife 80336 supplementary file 3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9436416/supplementaryFiles
    local path   : starplast/data/melting_temperature.tsv

Quirks that cost time once:
    The same file as the thermal-shift column and a different question: that one asks how far
    the curve MOVES when calcium is added, this asks where it SITS. Median 55.1 C, which is
    where protein melting temperatures live. Curves with R2 below 0.8 or a Tm outside 30-80 C
    are dropped -- the fit reports values up to 8,563, which is a failed fit and not a
    thermophile. Verified by reproducing across independent replicates at rho = +0.78 over 1,623
    proteins; a Tm that did not reproduce would be describing the run.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/melting_temperature.py
"""
from _common import run

KEY = "melting_temperature"

if __name__ == "__main__":
    run(KEY)
