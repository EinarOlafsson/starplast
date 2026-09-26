#!/usr/bin/env python3
"""Spatial proteome of the schizont (hyperLOPIT)

Which of 24 cellular niches each protein sits in, and the classifier's confidence

    level / kind : post_translation / LOPIT
    provides     : lopit_pf_location, lopit_pf_svm_score
    coverage     : 1,646 classified of 3,000 (29%)
    citation     : Chisholm SA et al., The spatial proteome of the Plasmodium falciparum schizont illuminates the composition and evolutionary trajectories of its organelles. Nat Commun 2026;17:6192
    PMID         : 42218142
    accession    : PXD070842 / Nat Commun Supplementary Data 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13369866/supplementaryFiles
    local path   : datasets/post_translation/LOPIT/42218142/41467_2026_73664_MOESM3_ESM.xlsx

Quirks that cost time once:
    The FIRST subcellular localization of any kind in the Plasmodium table, which had no
    compartment at all -- the question hyperLOPIT answers for Toxoplasma in this same map.
    Reproduces the paper exactly: 3,000 proteins mapped, 1,646 classified, 24 niches, and the
    markers where they must be (RAP1 rhoptries, MAHRP1 Maurer's cleft, ACP apicoplast, EXP2 and
    HSP101 at the PVM, GAPDH cytosol). `unknown` is shipped as a MISSING label, not a 25th
    niche. The two-experiment classifier (S1-S2) ships rather than the three-experiment one: it
    is the paper's headline and classifies more proteins, and shipping both would be one
    measurement twice.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_spatial_proteome.py
"""
from _common import run

KEY = "pf_spatial_proteome"

if __name__ == "__main__":
    run(KEY)
