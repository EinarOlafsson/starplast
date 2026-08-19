#!/usr/bin/env python3
"""CDPK1 substrates (thiophosphate labelling)

Thiophosphorylated peptides per gene from analog-sensitive CDPK1

    level / kind : post_translation / proteomics
    provides     : cdpk1_thiophospho_peptides
    coverage     : 361 genes
    PMID         : 37933960
    accession    : eLife 85654 supplementary file 6
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10629828/supplementaryFiles
    local path   : starplast/data/cdpk1_substrates.tsv

Quirks that cost time once:
    The enriched fraction only, one gene per peptide from the master accession the search engine
    assigned -- counting every protein a shared peptide maps to would credit ambiguous peptides
    several times. The count tracks abundance, as every phosphoproteomic count does, and the top
    of it is HSP70, HSP90 and BiP. What says it is nonetheless CDPK1's substrate set is the
    enrichment: microneme proteins are 11-fold over-represented (Fisher p = 2e-05) and CDPK1 is
    the kinase that governs microneme secretion, myosin A is in it, and so is the HOOK protein
    that the paper exists to report.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/cdpk1_substrates.py
"""
from _common import run

KEY = "cdpk1_substrates"

if __name__ == "__main__":
    run(KEY)
