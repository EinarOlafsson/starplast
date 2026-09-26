#!/usr/bin/env python3
"""Protein melting temperature across the blood-stage cycle

The temperature at which each protein leaves solution in an intact cell

    level / kind : translation / proteomics
    provides     : melting_temperature_tm, melting_temperature_sd
    coverage     : 2,039 proteins (35.6%)
    citation     : Pazicky S et al., MAP-X reveals distinct protein complex dynamics across Plasmodium falciparum blood stages. Nat Microbiol 2025;10:3229-3244
    PMID         : 41315737
    accession    : PXD056075 / Nat Microbiol Supplementary Table 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12615847/supplementaryFiles
    local path   : datasets/translation/proteomics/41315737/41564_2025_2173_MOESM3_ESM.xlsx

Quirks that cost time once:
    DERIVED, and it has to be: no melting temperature is published anywhere. PRIDE holds 281 GB
    of raw mass spectrometry and nothing else, the paper tabulates per-temperature abundances,
    and the authors' own package computes the melting point inside a fitting step without
    tabulating it. `scripts/fit_meltome.py` fits it with the acceptance the Toxoplasma column
    already uses (R-squared at least 0.8, 30-80 degrees), keeping 68.9% of 22,764 curves. Self-
    validating the way a melting point should be: replicates agree at rho 0.84 with a spread of
    1.15 degrees; subunits of one complex melt together (0.19 degrees across the T-complex, 0.42
    across the proteasome alpha ring, against 5.57 for size-matched random sets, p = 9e-5);
    glycolysis and the proteasome are stable near 62.5 degrees and HSP70/90 labile at 50.4; and
    it agrees with the independent Toxoplasma measurement through orthology at +0.28. ONE
    column, not seven: a protein's melting point varies 3.7 degrees across the cycle against
    1.15 of replicate noise, so per-stage columns would claim a resolution the data lacks. The
    worst run by fit acceptance is the replicate the authors themselves excluded, found
    independently here and dropped.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_meltome.py
"""
from _common import run

KEY = "pf_meltome"

if __name__ == "__main__":
    run(KEY)
