#!/usr/bin/env python3
"""Codon usage bias (COMPUTED)

Effective number of codons, GC3, and codon adaptation index

    level / kind : reference / sequence
    provides     : codon_enc, codon_gc3, codon_cai_ribosomal
    coverage     : 8,140 genes (100%)
    accession    : ToxoDB ME49
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/toxodb_cds.tsv.gz

Quirks that cost time once:
    COMPUTED here from the coding sequences fetched 2026-08-16. ENC is Wright's effective number
    of codons and GC3 the synonymous third-position GC, both reference-free. CAI's reference set
    is the 158 ribosomal proteins, chosen by product annotation and NOT by this map's expression
    columns -- the usual choice, 'the most highly expressed genes', would have built a sequence
    column out of an expression column and then found them correlated. Verified by the signs
    translational selection predicts: ribosomal proteins are more biased than the rest (ENC 46.5
    against 54.0), and CAI rises with transcription (rho +0.37) and with protein abundance (rho
    +0.24) while ENC falls with both. Those correlations are a finding here rather than a
    construction.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_codon_usage.py
"""
from _common import run

KEY = "toxodb_codon_usage"

if __name__ == "__main__":
    run(KEY)
