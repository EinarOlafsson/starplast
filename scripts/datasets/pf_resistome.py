#!/usr/bin/env python3
"""In vitro evolution resistome and field variation

How often a gene mutated under compound selection, and whether the paper calls it a target

    level / kind : DNA / in_vitro_evolution
    provides     : resistance_selection_clones, resistance_selection_compounds, resistance_selection_variants, resistance_target_compounds, pf6_field_dnds, pf6_field_nonsyn_snvs
    coverage     : 732 selected / 4,941 field genes
    citation     : Luth MR et al., Systematic in vitro evolution in Plasmodium falciparum reveals key determinants of drug resistance. Science 2024;386:eadk9893
    PMID         : 39607932
    accession    : Science Supplementary Data 3, 5, 6
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11809290/supplementaryFiles
    local path   : datasets/DNA/in_vitro_evolution/39607932/SupplementaryData3_SNVs-INDELs.xlsx

Quirks that cost time once:
    724 clones, each evolved against one of 118 compounds and sequenced whole-genome; both
    numbers reproduce from the deposit. Two kinds of column, and the difference is the point.
    The COUNTS are of coding variants per gene and their top is AP2-G (13 compounds) and PfEMP1,
    genes that mutate under prolonged culture whatever the drug, so a count is not a claim of
    resistance. `resistance_target_compounds` is the paper's own classification, restricted to
    the two classes its hypergeometric test supports, and its top is the canonical set: PfATP4
    and PfMDR1 at 5 compounds each, the prodrug-activating esterase at 4, then cytochrome b,
    CARL, PI4K beta and the Niemann-Pick C1-related protein at 3, with DHODH, PfCRT, the tRNA
    ligases and DHFR-TS behind them. Cytochrome b is in the deposit under the pre-2010 name
    `mal_mito_3` and matches no accession pattern: it is mapped on the deposit's own description
    naming exactly one product in the shipped table, because dropping it would lose the
    atovaquone gene and its 32 selected clones. A gene with no selected mutation gets NO value
    rather than a zero -- 118 compounds are not a test of the other 4,600 genes -- which is the
    opposite of the choice made for the R-DeeP flag, where the run did quantify every protein it
    reports. The Pf6 columns are field variation over 5,970 isolates, with the deposit's -1 for
    'not computable' read as missing.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_resistome.py
"""
from _common import run

KEY = "pf_resistome"

if __name__ == "__main__":
    run(KEY)
