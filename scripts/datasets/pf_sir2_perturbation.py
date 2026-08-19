#!/usr/bin/env python3
"""Plasmodium transcription under Sir2 knockout

Wild type and sir2a / sir2b knockout at ring, trophozoite and schizont

    level / kind : transcription / microarray
    provides     : sir2_wt_ring, sir2_wt_trophozoite, sir2_wt_schizont, sir2a_ko_ring, sir2a_ko_trophozoite, sir2a_ko_schizont, sir2b_ko_ring, sir2b_ko_trophozoite, sir2b_ko_schizont
    coverage     : 5,615 genes
    accession    : PlasmoDB Sir2 KO Marray
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : datasets/reference/plasmodb/plasmodb_pf3d7_sir2_perturbation.tsv

Quirks that cost time once:
    Shipped as stated CONDITIONS and not as a knockout-minus-wild-type contrast, which is the
    quantity anyone will want. The independent check on that contrast came out ambiguous and the
    column would have stated more confidence than there is: Sir2a silences subtelomeric var
    genes, so var should rise in the sir2a knockout, and it does in ring (+0.135, p = 4e-12) and
    schizont (+0.130, p = 2e-38) but FALLS in trophozoite (-0.240, p = 3e-22), with effects
    small against a spread of 0.7. That fits the canonical result being subset-specific and var
    probes cross-hybridising across sixty paralogues, but it is not a clean confirmation. The
    conditions themselves are unambiguous -- PlasmoDB names them, the values are log
    intensities, and the medians align across arrays within 0.1, which a test asserts because it
    is the precondition that makes differencing them meaningful at all.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_sir2_perturbation.py
"""
from _common import run

KEY = "pf_sir2_perturbation"

if __name__ == "__main__":
    run(KEY)
