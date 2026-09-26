#!/usr/bin/env python3
"""Population and between-species variation per gene

Non-synonymous variation in field isolates, and dN/dS against Plasmodium orthologs

    level / kind : reference / sequence
    provides     : field_pnps_adj, field_variant_fraction, dnds_laverania, dnds_plasmodium
    coverage     : 4,282-5,232 genes
    citation     : Chisholm SA et al., The spatial proteome of the Plasmodium falciparum schizont illuminates the composition and evolutionary trajectories of its organelles. Nat Commun 2026;17:6192
    PMID         : 42218142
    accession    : Nat Commun Supplementary Data 3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13369866/supplementaryFiles
    local path   : datasets/post_translation/LOPIT/42218142/41467_2026_73664_MOESM5_ESM.xlsx

Quirks that cost time once:
    The second half of the spatial-proteome paper, and a different question from its first:
    pN/pS and variant fraction are computed over FIELD isolates, so they answer the field-
    variation slot that had nothing, while dN/dS over Laverania and over Plasmodium orthologs
    measures selection since the species split. They are not restatements of the PlasmoDB lab-
    strain SNP columns already shipped: pN/pS agrees with them at rho 0.29, which is a related
    quantity measured on a different population, not a copy. Registered as its own entry because
    the columns belong to different slots from the localization ones and provenance is a
    statement about which measurements a dataset produced.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_field_variation.py
"""
from _common import run

KEY = "pf_field_variation"

if __name__ == "__main__":
    run(KEY)
