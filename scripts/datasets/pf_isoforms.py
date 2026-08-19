#!/usr/bin/env python3
"""Plasmodium long-read transcript models

Transcript models per gene, and how many the annotation does not contain

    level / kind : transcription / nanopore
    provides     : n_transcript_models, novel_transcript_models
    coverage     : 1,857 genes, 2,498 models, 238 novel
    PMID         : 40316999
    accession    : Malar J 05376 Supplementary Data 2
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12046715/supplementaryFiles
    local path   : datasets/transcription/isoforms/40316999/SupplementaryData2.xlsx

Quirks that cost time once:
    SQANTI classifications of long-read models. `full-splice_match` is the reference transcript
    recovered and is NOT counted as novel -- only novel_in_catalog, novel_not_in_catalog and
    fusion are, which is what the Toxoplasma column of the same name counts. Absence is
    sequencing depth rather than a statement that a gene has one transcript, so unseen genes
    stay missing instead of reading as 1. Found while looking for something else: this paper was
    opened for its m6A data, whose Pf arm turned out to be a 43-gene intersection with P. vivax
    rather than a methylome, and was refused for that -- the isoform table beside it is the
    usable one. The obvious correlation holds: more expressed genes yield more models (rho
    +0.36), which is detection depth and is why the count is not read as isoform diversity.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_isoforms.py
"""
from _common import run

KEY = "pf_isoforms"

if __name__ == "__main__":
    run(KEY)
