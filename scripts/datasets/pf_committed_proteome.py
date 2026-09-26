#!/usr/bin/env python3
"""Proteome of sexually committed parasites

How much more or less of each protein a committed parasite carries

    level / kind : translation / proteomics
    provides     : committed_vs_asexual_log2fc, committed_vs_asexual_fdr
    coverage     : 1,950 proteins (34%)
    citation     : Venugopal K et al., Defining the proteome of sexually committed parasites in Plasmodium falciparum. Mol Cell Proteomics 2026;25:101505
    PMID         : 41482054
    accession    : PXD059080 / MCP Table S3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12878696/supplementaryFiles
    local path   : datasets/translation/proteomics/41482054/mmc4.xlsx

Quirks that cost time once:
    Commitment, not the gametocyte: the decision a cycle before the stage V proteome already in
    the map, isolated by sorting on MSRP1, which this study establishes as the marker. The
    combined contrast over both reporter lines and both sorting directions ships, because the
    four arms are one comparison done four ways. MSRP1 itself is +1.68 at FDR 0 and that is a
    positive control and nothing more -- it is what the sort was done on. The check that means
    something is the paper's finding that merozoite surface proteins separate the populations:
    MSP1 +0.34 at FDR 0, MSP2 +0.67.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_committed_proteome.py
"""
from _common import run

KEY = "pf_committed_proteome"

if __name__ == "__main__":
    run(KEY)
