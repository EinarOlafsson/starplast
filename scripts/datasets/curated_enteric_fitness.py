#!/usr/bin/env python3
"""Enteric / sexual-cycle fitness per gene (CURATED)

Gene disruptions carried through the feline stage with oocyst output measured

    level / kind : reference / literature
    provides     : enteric_oocyst_yield, enteric_sporulation, enteric_measurements
    coverage     : 8 genes, 4 studies
    PMID         : 28288194
    accession    : PMIDs 28288194, 30728393, 36809045 and PMC12942651
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5363998/fullTextXML

Quirks that cost time once:
    THIRD curated source, and it exists because the sweeps that closed this slot looked for a
    POOLED SCREEN through the enteroepithelial stages and correctly found none -- nobody has put
    a barcoded library through a cat. But the slot asks whether disrupting a gene costs the
    parasite oocysts, and feeding one knockout to a cat answers that one gene at a time. Two
    papers were found and REFUSED for failing the bar: one says the cat experiment 'should be
    carried out', one says oocysts were seen but the numbers were 'not quantified'. Four of the
    ten rows are unchanged, and they are the strongest rows here -- deleting all four LEA genes
    at once left oocyst yield alone (30 against 34 million from paired kittens), which a single
    knockout could not have established because redundancy could have hidden it. Yield and
    sporulation are separate columns because HAP2 sheds a few oocysts that never sporulate while
    Grx5 sheds fewer that sporulate poorly, and those are different events. Magnitudes stay in
    the evidence text: they are not comparable across cats, strains and inocula, and one numeric
    column would invent a precision the experiments do not have.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/curated_enteric_fitness.py
"""
from _common import run

KEY = "curated_enteric_fitness"

if __name__ == "__main__":
    run(KEY)
