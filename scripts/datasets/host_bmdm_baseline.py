#!/usr/bin/env python3
"""Baseline mouse bone-marrow macrophage transcriptome

How much of each gene an unstimulated macrophage transcribes (TPM)

    level / kind : reference / transcription
    provides     : bmdm_tpm
    coverage     : 15,437 mouse genes
    citation     : John SV et al., Classically activated macrophages undergo functionally significant nucleotide metabolism remodelling driven by nitric oxide. Nat Metab 2025;7:1681-1702
    PMID         : 40759751
    accession    : GSE267544
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267544/suppl/GSE267544_fpkm_BMDM.txt.gz
    local path   : datasets/host/bmdm/GSE267544/GSE267544_fpkm_BMDM.txt.gz

Quirks that cost time once:
    The three unstimulated M0 samples only; the LPS + IFN-gamma arms are a stimulus, not
    infection, and do not answer any slot. Keyed on Ensembl through reviewed UniProt, never on
    the deposit's symbols, which are old (Emr1, Irg1). Markers behave: Lyz2 27,505 TPM, Cd68
    2,791, Csf1r 1,276, Adgre1 903, Alb 0.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_bmdm_baseline.py
"""
from _common import run

KEY = "host_bmdm_baseline"

if __name__ == "__main__":
    run(KEY)
