#!/usr/bin/env python3
"""Novel transcript models (Nanopore, via ToxoDB)

How many novel TALON transcript models long reads support for this gene

    level / kind : transcription / LongRead
    provides     : novel_transcript_models
    coverage     : 798 genes (10%)
    accession    : ToxoDB Stuart/Ralph nanopore
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByLongReadEvidence_tgonME49_Stuart_Ralph_nanopore_rnaSeqNextflow_RSRC/reports/attributesTabular
    local path   : starplast/data/toxodb_nanopore_isoforms.tsv

Quirks that cost time once:
    Incomplete-splice-match, novel-in-collection and novel-not-in-collection models, at five
    supporting reads or more. `Known` is excluded because the annotated model says nothing about
    isoform use, and `Genomic` because it is unspliced. Its power is in PRESENCE rather than
    magnitude -- most genes that have a novel model have one. Verified on that basis: genes with
    a novel model have a median of 6 exons against 4 for genes without, Mann-Whitney p = 9e-43,
    which is the relationship alternative splicing has to produce. Absent is NOT zero: a gene
    with no novel model here may simply not have been sequenced deeply enough.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_nanopore_isoforms.py
"""
from _common import run

KEY = "toxodb_nanopore_isoforms"

if __name__ == "__main__":
    run(KEY)
