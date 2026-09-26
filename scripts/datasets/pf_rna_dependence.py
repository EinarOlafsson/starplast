#!/usr/bin/env python3
"""RNA-dependent proteins (R-DeeP)

Whether a protein's complex falls apart when the RNA is digested

    level / kind : post_translation / RDeeP
    provides     : rna_dependent, rna_dependence_qvalue
    coverage     : 3,671 proteins (64%)
    citation     : Hollin T et al., Proteome-wide identification of RNA-dependent proteins and an emerging role for RNAs in Plasmodium falciparum protein complexes. Nat Commun 2024;15:1365
    PMID         : 38355719
    accession    : MassIVE MSV000093488 / Nat Commun Supplementary Data 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10866993/supplementaryFiles
    local path   : datasets/post_translation/RDeeP/38355719/41467_2024_45519_MOESM4_ESM.xlsx

Quirks that cost time once:
    A question nothing in either table asked: not what a protein binds, but whether RNA is
    holding its complex together. One lysate down a sucrose gradient twice, with and without
    RNase, 25 fractions each; a protein that shifts towards the light fractions was being held
    by RNA. It is NOT RNA binding -- a protein can shift because its partner binds RNA, which is
    why the authors call them RNA-dependent. The 898 significant shifters at q < 0.05 are the
    paper's own number, reproduced from the deposit's q-values, and the flag is 0 rather than
    blank for the other 2,773 proteins the run quantified: that experiment did test them. The
    classes come out the right way round -- RNA helicases enriched (28 of 61, odds 2.7, p 2e-4),
    the proteasome 0 of 14 -- while ribosomal proteins are DEPLETED (23 of 137), which is the
    reminder that this is not a column about binding RNA.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_rna_dependence.py
"""
from _common import run

KEY = "pf_rna_dependence"

if __name__ == "__main__":
    run(KEY)
