#!/usr/bin/env python3
"""Foldseek structural similarity (Plasmodium)

Which parasite proteins fold alike, without asking whether they are related

    level / kind : post_translation / structure
    provides     : n_struct_similar, edge:struct
    coverage     : 4,571 pairs / 1,620 genes at TM >= 0.7
    derived from : mean_plddt

Quirks that cost time once:
    COMPUTED HERE by `scripts/run_foldseek.py`, so there is nothing to download; the models come
    from the AlphaFold reference proteome UP000001450 (5,168 of them). TM-align rather than
    foldseek's faster 3Di+AA mode, and the same TM >= 0.7 cut as the Toxoplasma layer, because
    two arms whose `structural similarity` means two different numbers cannot be compared. The
    threshold is applied by the LOADER, not the search, so the shipped pair table keeps
    everything above the e-value cut and the cut can be revisited without re-running an hour of
    alignment.
    
    PREFILTERED, not exhaustive, and that is a caveat rather than a detail. The exhaustive
    search -- every one of the 26.7 million pairs TM-aligned -- was attempted twice and finished
    neither time, so no exhaustive result exists to compare this against and it CANNOT be
    claimed that no pair was missed. Two things argue the loss is small: at TM >= 0.7 the 3Di
    k-mer prefilter is retaining exactly what it was built to retain, and the per-query cap
    provably does not bind -- the busiest query returned 159 rows against a ceiling of 300. If a
    later run needs certainty, the tell that the cap has started to bind is queries returning
    exactly `--max-seqs` rows.
    
    Validated on what it should and should not say: structural neighbours share an orthogroup
    13.9% of the time against 0.35% for random pairs, a 39-fold enrichment, which is the
    positive control; and 3,937 of the 4,571 pairs (86%) join genes in DIFFERENT orthogroups,
    which is the layer's whole reason to exist. The Toxoplasma claim that it `reaches genes
    homology cannot` does NOT transfer verbatim: PlasmoDB gives all 5,720 genes an orthogroup,
    so there are none without one to reach, and the cross-orthogroup share is the honest form of
    that statement here.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_foldseek_struct.py
"""
from _common import run

KEY = "pf_foldseek_struct"

if __name__ == "__main__":
    run(KEY)
