#!/usr/bin/env python3
"""Calcium thermal-shift proteome (mineCETSA)

How far a protein's melting curve moves when calcium is added

    level / kind : post_translation / proteomics
    provides     : cetsa_calcium_ed_score
    coverage     : 2,348 proteins
    PMID         : 35976251
    accession    : PMC9436416 Supplementary file 3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9436416/supplementaryFiles
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/thermal_shift/

Quirks that cost time once:
    The authors' own Euclidean-distance score, not recomputed from the ten temperature points
    published beside it. Verified against the calcium sensors: CAM1 and CAM2 sit at the 98th
    percentile and CAM3 at the 83rd, and a protein whose melting curve does not move when
    calcium is added is not calcium-binding. One caveat for anyone comparing against the paper
    -- its headline conclusion is about PP1, and PP1 is unremarkable in THIS column. That claim
    comes from the zaprinast time course in the same paper, a different experiment; this is the
    calcium mineCETSA sheet. PXD033642, the deposit for the same study, publishes only
    identifications and could not have filled this slot.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/thermal_shift_cetsa.py
"""
from _common import run

KEY = "thermal_shift_cetsa"

if __name__ == "__main__":
    run(KEY)
