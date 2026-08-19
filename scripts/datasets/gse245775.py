#!/usr/bin/env python3
"""Differentiation ribosome profiling (eIF1.2)

RPF and RNA counts, and their ratio, in tachyzoites and pre-bradyzoites

    level / kind : translation / RiboSeq
    provides     : rpf245775_parent_tachy_r1, rpf245775_parent_tachy_r2, rpf245775_parent_tachy_r3, rpf245775_parent_prebrady_r1, rpf245775_parent_prebrady_r2, rpf245775_parent_prebrady_r3, rpf245775_eif12ko_tachy_r1, rpf245775_eif12ko_tachy_r2, rpf245775_eif12ko_tachy_r3, rpf245775_eif12ko_prebrady_r1, rpf245775_eif12ko_prebrady_r2, rpf245775_eif12ko_prebrady_r3, rna245775_parent_tachy_r1, rna245775_parent_tachy_r2, rna245775_parent_tachy_r3, rna245775_parent_prebrady_r1, rna245775_parent_prebrady_r2, rna245775_parent_prebrady_r3, rna245775_eif12ko_tachy_r1, rna245775_eif12ko_tachy_r2, rna245775_eif12ko_tachy_r3, rna245775_eif12ko_prebrady_r1, rna245775_eif12ko_prebrady_r2, rna245775_eif12ko_prebrady_r3, te245775_parent_tachy_r1, te245775_parent_tachy_r2, te245775_parent_tachy_r3, te245775_parent_prebrady_r1, te245775_parent_prebrady_r2, te245775_parent_prebrady_r3, te245775_eif12ko_tachy_r1, te245775_eif12ko_tachy_r2, te245775_eif12ko_tachy_r3, te245775_eif12ko_prebrady_r1, te245775_eif12ko_prebrady_r2, te245775_eif12ko_prebrady_r3
    coverage     : 7,880 genes (97%)
    PMID         : 38782906
    accession    : GSE245775
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245775/suppl/GSE245775_RAW.tar
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/stage_conversion_phenotype/

Quirks that cost time once:
    The `prebrady` arms are 48 hours in RPMI pH 8.3 at ambient CO2, which the submitters label
    `cell type: pre-bradyzoites`. Named for that and not for `bradyzoite`: these are not tissue
    cysts. Verified by reproducing the paper's own result from the counts -- BFD1 rises 3.6 log2
    on conversion in the parental line and 2.6 in the knockout, and BFD2 rises 1.4 and 0.1, so
    the knockout's failure to induce BFD2 is visible in the column itself. LDH2, BAG1 and SRS
    also rise on conversion, which is what confirms the arms are not swapped. Reached the map
    proposed for `stage-conversion phenotype`, which it is not.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse245775.py
"""
from _common import run

KEY = "gse245775"

if __name__ == "__main__":
    run(KEY)
