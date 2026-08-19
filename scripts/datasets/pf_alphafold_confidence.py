#!/usr/bin/env python3
"""Plasmodium model confidence and disorder (AlphaFold DB)

Mean pLDDT per protein, and the fraction of it at each confidence band

    level / kind : reference / structure
    provides     : mean_plddt, plddt_fraction_very_low, plddt_fraction_low, plddt_fraction_confident, plddt_fraction_very_high, alphafold_accession
    coverage     : 5,098 of 5,720 genes
    accession    : AlphaFold DB API, per UniProt accession
    url          : https://alphafold.ebi.ac.uk/api/prediction/
    local path   : datasets/reference/plasmodb/plasmodb_pf3d7_alphafold.tsv

Quirks that cost time once:
    Fetched per protein because the bulk proteome archive for this organism is not where the
    documented path says it is. The FRACTIONS matter as much as the mean: a protein half well-
    folded and half disordered has the same mean as one uniformly mediocre, and the slot asks
    about disorder as well as confidence. A gene can carry several UniProt accessions -- 876 do,
    mostly the variant surface families where each field isolate's allele has its own entry --
    so the fetch tries them in order and `alphafold_accession` records which one supplied the
    model, or a number could not be traced back to a structure. Validated on an ordering rather
    than a total: proteins carrying a recognised InterPro domain model at median pLDDT 73.5
    against 56.0 for those without (p = 2e-159), because a domain is a thing that folds. The
    correlation with protein length is NEGATIVE at -0.555, which is not a fault -- it is this
    proteome's low-complexity asparagine insertions, which are long and disordered.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_alphafold_confidence.py
"""
from _common import run

KEY = "pf_alphafold_confidence"

if __name__ == "__main__":
    run(KEY)
