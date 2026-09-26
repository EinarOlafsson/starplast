#!/usr/bin/env python3
"""High-resolution ribosome profiling (5'UTR study)

Translation efficiency per replicate, footprints over matched RNA

    level / kind : translation / Ribo-seq
    provides     : te302107_tachy_r1, te302107_tachy_r2
    coverage     : 5,992 (73.6%)
    citation     : Peters ML et al., 5' untranslated regions tune Toxoplasma translation. bioRxiv 2025, doi:10.1101/2025.07.14.664749 (preprint)
    accession    : GSE302107
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302107/suppl/GSE302107_RPKM_and_TE.xlsx
    local path   : datasets/translation/riboseq/GSE302107/GSE302107_RPKM_and_TE.xlsx

Quirks that cost time once:
    ME49 in HFF, two replicates of paired footprints and total RNA. The authors' TE is a LINEAR
    ratio (verified: footprint RPKM / RNA RPKM for every row); shipped as log2 to match the
    other TE columns, with a zero ratio -- no footprints counted -- as missing rather than minus
    infinity. The most reproducible TE in the organism (replicates rho 0.973, against 0.87 for
    GSE99395 and 0.67 for GSE245775) and concordant with both (rho 0.78-0.79). Ribosomal
    proteins are the high-TE class (median ratio 2.86 against 1.08). A preprint: cite the
    journal version once it exists. The host sheet (7,943 human genes) is not used yet.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse302107_riboseq.py
"""
from _common import run

KEY = "gse302107_riboseq"

if __name__ == "__main__":
    run(KEY)
