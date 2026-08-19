#!/usr/bin/env python3
"""Plasmodium palmitome (observed only)

Proteins observed S-palmitoylated, with the motif prediction deliberately excluded

    level / kind : post_translation / palmitoylation
    provides     : is_palmitoylated
    coverage     : 503 proteins
    PMID         : 36250062
    accession    : Front Cell Infect Microbiol Table 3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9556994/supplementaryFiles
    local path   : datasets/post_translation/palmitome/36250062/Table_3.xlsx

Quirks that cost time once:
    Read from the OBSERVED column and not from the sheet named for it. The same workbook carries
    `nrPalmitoylatedProteins`, whose name says palmitoylated and whose 3,105 rows are the union
    of palmitoyl-ABLE -- a motif prediction over 2,902 proteins -- and the 503 actually
    observed. Taking that sheet at its name would have called 54% of the proteome palmitoylated,
    against published palmitomes of 400 to 500, and the first rows being PfEMP1 and rifin is
    what gave it away. Validated on substrates and on mechanism rather than on a total: GAP45
    and CDPK1, the canonical Plasmodium substrates, are both present, and membrane proteins are
    enriched 2.1-fold among the palmitoylated (44% against 27%, p = 7e-15), which is what a
    membrane-anchoring modification has to do. ARO is a known miss -- no palmitome is complete,
    and absence here means not observed.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_palmitome.py
"""
from _common import run

KEY = "pf_palmitome"

if __name__ == "__main__":
    run(KEY)
