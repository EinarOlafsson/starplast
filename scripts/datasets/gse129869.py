#!/usr/bin/env python3
"""Host-context parasite ribosome profiling

Parasite ribosome footprints, RNA and translation efficiency in two HFF states

    level / kind : translation / Ribo-seq
    provides     : rpf129869_confluent_r1, rpf129869_confluent_r2, rpf129869_confluent_r3, rpf129869_subconfluent_r1, rpf129869_subconfluent_r2, rpf129869_subconfluent_r3, rna129869_confluent_r1, rna129869_confluent_r2, rna129869_confluent_r3, rna129869_subconfluent_r1, rna129869_subconfluent_r2, rna129869_subconfluent_r3, te129869_confluent_r1, te129869_confluent_r2, te129869_confluent_r3, te129869_subconfluent_r1, te129869_subconfluent_r2, te129869_subconfluent_r3
    coverage     : measured at build time
    citation     : Holmes MJ et al., Simultaneous Ribosome Profiling of Human Host Cells Infected with Toxoplasma gondii. mSphere 2019;4:e00292-19
    PMID         : 31167946
    accession    : GSE129869
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129869/suppl/GSE129869_RAW.tar
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE129869_RAW.tar

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse129869.py
"""
from _common import run

KEY = "gse129869"

if __name__ == "__main__":
    run(KEY)
