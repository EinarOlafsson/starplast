#!/usr/bin/env python3
"""Oxidative-stress CRISPR screen

Screening score per gene under oxidative challenge

    level / kind : DNA / CRISPR_screen
    provides     : oxidative_stress_screen_score
    coverage     : 7,384 genes (91%)
    PMID         : 34163449
    accession    : PMC8216390 Data Sheet 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8216390/supplementaryFiles
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/fitness_oxidative_stress/

Quirks that cost time once:
    The authors' own `Screening score` sheet, not a recomputation from the guide counts beside
    it in the same workbook -- they published the score, so it is theirs to define. Negative is
    required. Verified by the sign and the extreme: catalase comes out at -6.15, essentially the
    bottom of the whole screen, and it is the enzyme that disposes of hydrogen peroxide.
    Peroxiredoxin (-1.57), superoxide dismutase (-1.01), glutaredoxin (-0.71) and thioredoxin
    (-0.63) all sit below the genome median of -0.38.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/oxidative_stress_screen.py
"""
from _common import run

KEY = "oxidative_stress_screen"

if __name__ == "__main__":
    run(KEY)
