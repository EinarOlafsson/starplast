#!/usr/bin/env python3
"""Plasmodium transcription at febrile temperature

Wild type and two mutants at 37 C and at the 41 C of a malarial fever

    level / kind : transcription / RNAseq
    provides     : febrile_wt_37c, febrile_wt_41c, febrile_lrr5ko_37c, febrile_lrr5ko_41c, febrile_dhcko_37c, febrile_dhcko_41c
    coverage     : 5,791 genes
    accession    : PlasmoDB Pfal3D7 Febrile temps RNA-Seq
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : datasets/reference/plasmodb/plasmodb_pf3d7_febrile.tsv

Quirks that cost time once:
    Shipped as CONDITIONS rather than as a 41-versus-37 contrast, the same restraint as the Sir2
    entry and for a different reason. There the check on the contrast contradicted itself; here
    it came out NULL. Heat shock proteins move by a median log2 of +0.08 against -0.07 for
    everything else (p = 0.2), so a fever does not measurably induce them -- which is consistent
    with what is known, since this organism's chaperones are constitutively high rather than
    stress-induced, and the genes that do rise are Maurer's cleft two-TM proteins and stevor at
    four to five log2, matching published fever-driven surface remodelling. But a null result on
    the one available prediction is not a validation, and a derived column would imply it had
    passed one. A test asserts the arms are on a comparable scale, which is the precondition for
    the caller making the contrast themselves.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_febrile_stress.py
"""
from _common import run

KEY = "pf_febrile_stress"

if __name__ == "__main__":
    run(KEY)
