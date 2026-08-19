#!/usr/bin/env python3
"""Plasmodium intraerythrocytic cycle timing

When in the 48-hour cycle each transcript peaks, and how strongly it cycles

    level / kind : transcription / RNAseq
    provides     : idc_peak_hour, idc_cycling_amplitude
    coverage     : 5,038 genes timed of 5,499; 461 do not cycle strongly enough to place
    PMID         : 34668757
    accession    : GSE163144
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163144/suppl/
    local path   : datasets/transcription/idc_timecourse/34668757/

Quirks that cost time once:
    Sixteen timepoints three hours apart across 48 hours, two replicates, read from ONE cell of
    the deposit's square: the 3D7 line in normal (HbAA) red cells. The study's variable is
    sickle-trait haemoglobin and its second line is FUP, so the other three cells are
    perturbations rather than a reference -- and which sample is which comes from the deposit's
    own metadata file rather than from the sample titles, which encode the genotype and not the
    hour. Samples are TPM and sum to 1e6 exactly, checked rather than assumed. The peak is a
    PHASE, not the largest column: the axis wraps, since hour 48 is hour 0 of the next cycle,
    and an argmax splits the invasion peak between the last timepoint and the first -- it put
    AMA1 and PTRAMP, textbook invasion transcripts, at hour 3, where a culture synchronised at
    invasion is still carrying the merozoite's mRNA. Fitted through the first Fourier harmonic
    (`cellcycle.cyclic_phase`, shared rather than written here) the markers land where the
    biology says: MSP1 45.7 h, SERA5 42.9 h, RhopH2 36.2 h, KAHRP 25.1 h, SBP1 15.1 h -- and
    AMA1 at 2.5 h, which is five hours from MSP1 across the wrap and not the other way round.
    READ IT AS A CIRCLE. Genes whose first harmonic explains less than 40% of their variation
    are left missing, because a flat profile still has an angle.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_idc_timing.py
"""
from _common import run

KEY = "pf_idc_timing"

if __name__ == "__main__":
    run(KEY)
