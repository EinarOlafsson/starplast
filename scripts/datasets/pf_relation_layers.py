#!/usr/bin/env python3
"""Plasmodium relation layers (COMPUTED)

Gene pairs sharing an orthogroup or a domain, and pairs whose stages covary

    level / kind : reference / graph
    provides     : edge:orthogroup, edge:domain, edge:coexpression
    coverage     : 1,741 + 24,123 + 63,158 pairs
    local path   : starplast/data/pf_graph.npz
    derived from : orthogroup, interpro_ids, expr_ring, expr_schizont, expr_sporozoite

Quirks that cost time once:
    A SECOND graph file, because an edge is a pair of indices into a table and a falciparum gene
    has no index in the Toxoplasma one. Constructions are copied from the Toxoplasma arm rather
    than re-invented, so that a difference between the arms means the biology differs and not
    that the edges were drawn by different rules. The one deliberate change is the domain
    weight. Two proteins can share more than one InterPro domain, and emitting the pair once per
    shared domain draws the same edge repeatedly, which reads as repeated evidence; the naive
    construction gave 10,764 duplicate emissions among 34,887 here, nearly all inside the var,
    rifin and stevor families that share whole multi-domain architectures. The count is now the
    weight, so a pair sharing eleven domains says so. The Toxoplasma arm emits no duplicates at
    all today -- checked rather than assumed -- and would acquire the same fault the moment its
    annotation gained a pair sharing two domains.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_relation_layers.py
"""
from _common import run

KEY = "pf_relation_layers"

if __name__ == "__main__":
    run(KEY)
