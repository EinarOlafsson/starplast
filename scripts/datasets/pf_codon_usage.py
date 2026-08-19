#!/usr/bin/env python3
"""Plasmodium codon usage (COMPUTED)

Effective number of codons, GC3, and CAI against the ribosomal proteins

    level / kind : reference / annotation
    provides     : codon_enc, codon_gc3, codon_cai_ribosomal
    coverage     : 5,318 genes
    accession    : PlasmoDB-68 Pfalciparum3D7 AnnotatedCDSs
    url          : https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-68_Pfalciparum3D7_AnnotatedCDSs.fasta
    local path   : starplast/data/plasmodb_cds.tsv.gz
    derived from : length

Quirks that cost time once:
    COMPUTED through the SAME code as the Toxoplasma arm rather than reimplemented --
    `codons.codon_usage` gained a table and node-table parameter for it. ENC and GC3 are
    definitions and the CAI reference set is 'the ribosomal proteins' in both arms, so two
    implementations could only differ by being wrong in one of them. Sharing it buys the first
    measurement the two arms can be COMPARED on, and the comparison is the validation: GC3
    median 0.150 here against 0.583 in Toxoplasma, and ENC 37.6 against 53.9 -- P. falciparum
    has the most AT-rich genome of any eukaryote, so extreme codon bias is what has to appear,
    and a test fails if the two arms ever converge. The sequence report API returned 422, 400
    and 500 to three different request shapes; the static release FASTA is what works.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_codon_usage.py
"""
from _common import run

KEY = "pf_codon_usage"

if __name__ == "__main__":
    run(KEY)
