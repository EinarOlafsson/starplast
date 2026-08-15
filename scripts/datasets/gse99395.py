#!/usr/bin/env python3
"""Intracellular/extracellular ribosome profiling

Ribosome footprints, matched RNA and relative translation efficiency

    level / kind : translation / Ribo-seq
    provides     : rpf99395_extracellular_r1, rpf99395_extracellular_r2, rpf99395_intracellular_r1, rpf99395_intracellular_r2, rna99395_extracellular_r1, rna99395_extracellular_r2, rna99395_intracellular_r1, rna99395_intracellular_r2, te99395_extracellular_r1, te99395_extracellular_r2, te99395_intracellular_r1, te99395_intracellular_r2
    coverage     : measured at build time
    citation     : Hassan MA et al., Comparative ribosome profiling uncovers a dominant role for translational control in Toxoplasma gondii. BMC Genomics 2017;18:961
    PMID         : 29228904
    accession    : GSE99395
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99395/suppl/GSE99395_Raw_counts.txt.gz
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE99395_Raw_counts.txt.gz

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse99395.py
"""
from _common import run

KEY = "gse99395"

if __name__ == "__main__":
    run(KEY)
