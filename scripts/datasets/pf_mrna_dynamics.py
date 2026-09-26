#!/usr/bin/env python3
"""mRNA synthesis and decay rates through the blood-stage cycle

Transcripts made per minute, and transcripts lost per minute, at each gene's peak

    level / kind : transcription / RNAseq
    provides     : transcription_rate_4tu, mrna_decay_rate_4tu
    coverage     : 4,373 / 4,420 genes
    citation     : Painter HJ, Chung NC, Sebastian A, Albert I, Storey JD, Llinas M. Genome-wide real-time in vivo transcriptional dynamics during Plasmodium falciparum blood-stage development. Nat Commun 2018;9:2656
    PMID         : 29985403
    accession    : GSE114621 / Nat Commun Supplementary Data 2
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6037754/supplementaryFiles
    local path   : datasets/transcription/RNAseq/29985403/41467_2018_4966_MOESM5_ESM.xlsx

Quirks that cost time once:
    4-thiouracil labelling at every one of the 48 hours, which is why this is the only RNA-
    stability measurement for this organism. FLUXES, in transcripts per minute at the gene's own
    peak -- NOT half-lives, and not comparable with the Toxoplasma actinomycin columns, which
    are fractions remaining: an abundant transcript loses more transcripts per minute than a
    scarce stable one, and the decay column tracks abundance at rho 0.42 for exactly that
    reason. The paper's claim that transcription runs in every stage reproduces (616-962 genes
    peak in each of the six windows), and the timing agrees with a series this study had no part
    in: of the genes whose transcription peaks in a ring window, 87-92% also peak in the shipped
    ring expression column. The peak-stage labels in the deposit are NOT shipped -- they would
    restate the stage expression series the table already carries.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_mrna_dynamics.py
"""
from _common import run

KEY = "pf_mrna_dynamics"

if __name__ == "__main__":
    run(KEY)
