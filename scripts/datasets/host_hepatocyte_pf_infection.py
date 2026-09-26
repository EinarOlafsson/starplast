#!/usr/bin/env python3
"""Human hepatocyte response to P. falciparum

Infected against uninfected primary hepatocytes, per host gene: log2 change and FDR

    level / kind : reference / transcription
    provides     : hepatocyte_pf_infection_log2fc, hepatocyte_pf_infection_padj
    coverage     : 12,550 human genes
    citation     : Kanyal A et al., Inhibition of Wnt signaling in primary human hepatocytes promotes Plasmodium falciparum liver stage development. PLoS Pathog 2025;21:e1013800
    PMID         : 41428747
    accession    : GSE263643
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263643/suppl/GSE263643_deseq2-normalized_counts.txt.gz
    local path   : datasets/host/hepatocyte_infection/GSE263643/GSE263643_deseq2-normalized_counts.txt.gz

Quirks that cost time once:
    The paper never compares infected with uninfected wells; the contrast is made here from its
    DESeq2-normalized counts, treatment and donor blocked, freshly isolated cells excluded. WEAK
    by construction -- bulk wells in which few cells are infected, two donors -- and only 25
    genes reach padj < 0.05; what does is a type-I interferon response in both strains (IFI44L,
    CXCL10, CXCL11, RSAD2, IFIT1). Whether the uninfected wells received mock mosquito material
    is not described, so part of that signal may come from the sporozoite preparation.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_hepatocyte_pf_infection.py
"""
from _common import run

KEY = "host_hepatocyte_pf_infection"

if __name__ == "__main__":
    run(KEY)
