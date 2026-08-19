#!/usr/bin/env python3
"""Drug sensitivity per gene (CURATED)

Knockouts with a measured shift in sensitivity to a named compound

    level / kind : reference / literature
    provides     : drug_compounds_tested, drug_sensitivity_shifts, drug_sensitivity_directions
    coverage     : 3 genes, 2 compounds
    PMID         : 41025776
    accession    : mBio, PMID 41025776
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12607627/fullTextXML

Quirks that cost time once:
    SECOND curated source, and it exists because six catalogue sweeps looked for a genome-wide
    chemogenomic screen and correctly found none -- while the sentence they tested, 'no screen
    exists', is not the same statement as 'the question cannot be answered'. That is the same
    error the lipid slot took three passes to notice. The bar is a MEASURED shift under a NAMED
    compound, not 'an inhibitor of this protein kills the parasite', which is target engagement
    and has its own slot; a test asserts no row's evidence reads that way. `unchanged` rows are
    KEPT -- a transporter deleted with no effect on analog sensitivity is a result, and dropping
    those would leave the column looking like a list of hits. TgENT3's two rows were measured in
    a ΔTgAT1 background and say so, because reading a double mutant's phenotype off one of its
    genes is its own error. The per-row product check earned its keep immediately: the
    annotation calls TGME49_244440 'adenosine transporter AT1', which is independent
    confirmation the accession is TgAT1.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/curated_drug_sensitivity.py
"""
from _common import run

KEY = "curated_drug_sensitivity"

if __name__ == "__main__":
    run(KEY)
