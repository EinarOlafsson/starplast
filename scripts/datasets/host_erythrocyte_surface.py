#!/usr/bin/env python3
"""Human red blood cell SURFACE proteome, by population

Which host proteins are reachable from outside the cell the merozoite invades

    level / kind : reference / proteomics
    provides     : rbc_surface_copies_uk, rbc_surface_copies_senegal, rbc_surface_found_uk, rbc_surface_found_senegal
    coverage     : 267 plasma-membrane proteins, 230 in both populations
    PMID         : 31552303
    accession    : Commun Biol 0596 Supplementary Data 2A
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6754445/supplementaryFiles
    local path   : datasets/host/erythrocyte/31552303/42003_2019_596_MOESM6_ESM.xlsx

Quirks that cost time once:
    Plasma membrane profiling, so this is the OUTSIDE of an intact red cell rather than
    everything in it -- a different question from the fractionated proteome beside it, and the
    one an invasion receptor slot asks. The authors' stated aim is candidate Plasmodium
    receptors. The two donor populations, nine UK and nine Senegalese, ship as their own columns
    because the difference IS the result: averaging a Duffy-positive population with a Duffy-
    negative one would erase the best-known receptor polymorphism in malaria. A zero is the
    paper's own encoding for `not identified in this population`, which is why the found flags
    ship beside the counts -- for ACKR1 the zero is the West African Duffy-negative phenotype
    and not a detection failure, and only the flag separates those. Reconciles with the paper's
    own sheets: 230 proteins in both populations, 11 in UK donors only, 26 in Senegalese donors
    only. Self-validating against numbers measured long before mass spectrometry -- band 3 at
    1.3 million copies per cell and glycophorin A at 3.3 million, with basigin present in both
    populations.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_erythrocyte_surface.py
"""
from _common import run

KEY = "host_erythrocyte_surface"

if __name__ == "__main__":
    run(KEY)
