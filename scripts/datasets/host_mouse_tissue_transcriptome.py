#!/usr/bin/env python3
"""FANTOM5 mouse brain and skeletal muscle

How much of each gene the mouse tissues a bradyzoite persists in transcribe

    level / kind : reference / transcription
    provides     : brain_tpm, skeletal_muscle_tpm
    coverage     : 13,739 mouse genes; brain 9,763, skeletal muscle 8,635
    PMID         : 24670764
    accession    : E-MTAB-3579
    url          : https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-3579/download/zip?fileType=rnaseq-baseline-rpkms&accessKey=
    local path   : datasets/host/fantom5/24670764/E-MTAB-3579-mouse-tissue-tpm.tsv

Quirks that cost time once:
    Two caveats that must travel with these numbers.
    
    First, a BLANK here does not mean unmeasured. This is the Expression Atlas's protein-coding,
    above-cutoff export, and both columns bottom out at exactly 0.5 -- so a gene present in the
    table with no value was measured and fell below 0.5 TPM. It is still left as missing rather
    than written as zero, because the file does not say which of the two it is for genes absent
    from it entirely.
    
    Second, SKELETAL MUSCLE is juvenile. It is the only skeletal-muscle sample in the atlas;
    biceps femoris is skeletal muscle so the tissue is right and the age is not, and it ships
    with the age said rather than quietly relabelled adult. Brain, by contrast, is the mean of
    the four ADULT regions -- cerebral cortex, cerebellum, hippocampal formation, olfactory
    brain -- an average across regions rather than replicates, which is what the slot wants
    because a cyst is not confined to one region.
    
    These are CAGE tag values labelled TPM; they order genes within a tissue and across tissues,
    but are not interchangeable with RNA-seq TPM. Validated on markers that must separate and
    do, by four to five orders of magnitude in both directions: Acta1 is 102,086 in muscle
    against 1.0 in brain and Ckm 54,636 against 0.6, while Snap25 is 1,171 in brain and below
    cutoff in muscle, as is Gfap at 34.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_mouse_tissue_transcriptome.py
"""
from _common import run

KEY = "host_mouse_tissue_transcriptome"

if __name__ == "__main__":
    run(KEY)
