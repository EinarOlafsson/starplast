#!/usr/bin/env python3
"""Co-translation layer (COMPUTED)

Gene pairs whose ribosome footprints covary

    level / kind : translation / RiboSeq
    provides     : edge:cotranslation
    coverage     : 6,231 edges over 7,437 genes
    derived from : rpf99395_intracellular_r1, rpf129869_confluent_r1, rpf245775_parent_tachy_r1

Quirks that cost time once:
    COMPUTED here: the same construction co-expression uses over transcripts, applied to the 22
    ribosome-footprint columns, at r >= 0.95 and the top 25 neighbours. DERIVED, so it is
    declared -- an embedding built on the RPF columns must not then be validated against this
    layer. It is not a copy of co-expression: 97% of its edges are not co-expression edges,
    Jaccard 0.003. What says it is co-TRANSLATION is that ribosomal proteins pair with each
    other 464 times where chance gives 3; they are made together stoichiometrically, which is
    the textbook case of co-translational regulation. 85 of its edges are also measured
    crosslink contacts.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/cotranslation_edges.py
"""
from _common import run

KEY = "cotranslation_edges"

if __name__ == "__main__":
    run(KEY)
