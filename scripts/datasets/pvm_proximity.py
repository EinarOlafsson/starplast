#!/usr/bin/env python3
"""PVM proximity labelling

Whether the study placed this protein at the parasitophorous vacuole membrane

    level / kind : post_translation / proteomics
    provides     : pvm_proximity_positive
    coverage     : 1,274 genes (73 positive)
    PMID         : 34749525
    accession    : mBio 00260-21 Data Set S1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8576527/supplementaryFiles
    local path   : starplast/data/pvm_proximity.tsv

Quirks that cost time once:
    One of the few binaries in the map with REAL zeros: the study publishes a likely-negative
    list of 1,201 genes beside its 73 positives, so a zero here is a measurement and not a gap
    -- for those 1,201 and only for them. Everything else is NaN. Verified by what the positives
    are: 53 of 73 are dense granule proteins against 0 of 1,201 negatives (Fisher p = 2e-77),
    and dense granule proteins are exactly what Toxoplasma secretes into the vacuole and inserts
    into the membrane it shares with the host cytosol.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pvm_proximity.py
"""
from _common import run

KEY = "pvm_proximity"

if __name__ == "__main__":
    run(KEY)
