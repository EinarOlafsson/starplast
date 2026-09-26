#!/usr/bin/env python3
"""Proteome and transcriptome without iron

Change in each protein, and in each transcript, after 24 h of iron depletion

    level / kind : translation / proteomics
    provides     : iron_depletion_protein_log2fc, iron_depletion_protein_padj, iron_depletion_rna_log2fc
    coverage     : 5,047 protein / 3,113 RNA
    citation     : Hanna JC et al., Global translational and metabolic remodeling during iron deprivation in Toxoplasma gondii. mBio 2026;17:e0378825
    PMID         : 41925342
    accession    : mBio Tables S1 and S2
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13170339/supplementaryFiles
    local path   : datasets/translation/proteomics/41925342/mbio.03788-25-s0002.xlsx

Quirks that cost time once:
    The same paper the metabolome_iron row cites, whose proteome and RNA-seq that row never
    took. Genome-scale on the protein side (5,052 groups, 3 + 3 replicates); the RNA side is
    only the 3,113 genes of the paper's joint analysis, a significance-filtered subset, so a
    missing RNA value is 'not reported' rather than 'unchanged'. The paper's 61 iron-sulfur
    proteins shift down on balance (median -0.04 against +0.02; 12 of the 16 that change
    significantly fall), a modest effect. Despite the title there is NO ribosome profiling in
    this paper -- translation is measured by microscopy -- so it fills protein abundance under
    stress, not translation, and the earlier reading of it as a translation dataset was wrong.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/iron_depletion_proteome.py
"""
from _common import run

KEY = "iron_depletion_proteome"

if __name__ == "__main__":
    run(KEY)
