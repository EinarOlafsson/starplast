#!/usr/bin/env python3
"""m6A methylation per transcript

How many canonical methylation sites a transcript has, and how fully methylated

    level / kind : transcription / RNAseq
    provides     : m6a_n_canonical_sites, m6a_canonical_stoichiometry
    coverage     : 5,285 genes (92%)
    citation     : Levendis JM et al., m6A positions polyadenylation in Plasmodium falciparum. bioRxiv 2026, doi:10.64898/2026.05.19.726191 (preprint)
    accession    : bioRxiv 10.64898/2026.05.19.726191 Supplementary Tables 6-8
    url          : https://www.biorxiv.org/content/10.64898/2026.05.19.726191v1.supplementary-material
    local path   : datasets/transcription/m6a/pf_m6a_2026/SupplementaryTables_media-2.xlsx

Quirks that cost time once:
    Nanopore direct RNA sequencing, which reads the modification on the molecule rather than
    inferring it from an antibody pulldown -- so the numbers are per transcript and per site
    rather than per peak. Sites are the union over three time points; stoichiometry is averaged
    only over time points where the transcript was read at least 20 times, because below that a
    fraction is a coin-flip. Validated against an orthogonal chemistry: 67% of GLORI-seq sites
    fall within one nucleotide of a nanopore site. Only the control arm is used. The deposit's
    two groups are UNLABELLED; the second has about half the methylation, consistent with the
    paper's knock-sideways of the methyltransferase, which is enough to choose the control arm
    but not enough to ship a difference, so no difference is shipped.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_m6a_nanopore.py
"""
from _common import run

KEY = "pf_m6a_nanopore"

if __name__ == "__main__":
    run(KEY)
