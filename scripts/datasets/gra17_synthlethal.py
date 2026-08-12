#!/usr/bin/env python3
"""GRA17 synthetic-lethal screen

RH and RH-delta-gra17 phenotype by passage; MAGeCK p-values

    level / kind : DNA / CRISPR_screen
    provides     : crispr_gra17ko_phenotype, crispr_gra17_synthlethal_delta, crispr_gra17_candidate
    coverage     : 7,553 (genome-wide)
    citation     : Genome-wide CRISPR screen identifies genes synthetically lethal with GRA17, a nutrient channel encoding gene in Toxoplasma
    PMID         : 37498952
    url          : https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1011543.s001&type=supplementary
    local path   : datasets/DNA/CRISPR_screen/37498952/

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gra17_synthlethal.py
"""
from _common import run

KEY = "gra17_synthlethal"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
