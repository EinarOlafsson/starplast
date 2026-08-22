#!/usr/bin/env python3
"""GTEx median expression in two host cell types

How much of each gene the host cell transcribes, in the cell each parasite lives in

    level / kind : reference / transcription
    provides     : fibroblast_tpm, hepatocyte_tpm
    coverage     : 19,087 human genes over 2 tissues
    PMID         : 32913098
    accession    : GTEx v10
    url          : https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
    local path   : datasets/host/gtex/32913098/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

Quirks that cost time once:
    TWO of GTEx's sixty-eight tissues, and which two is the whole judgement.
    `Cells_Cultured_fibroblasts` IS the cell Toxoplasma is grown in, and v10's
    `Liver_Hepatocyte` is laser-captured hepatocyte rather than liver tissue -- the difference
    between answering a cell-type slot and substituting the organ around it. `Liver` and
    `Skin_*` are in the same file and are deliberately NOT taken: a liver is not a hepatocyte
    and skin is not dermis, which is the substitution this campaign already measured and refused
    when it looked at the Human Protein Atlas for cell-type slots.
    
    Keyed through `host.uniprot_index`, because GTEx is Ensembl and the host table is UniProt.
    Reviewed entries only: the raw mapping is 74% ambiguous since it lists every TrEMBL fragment
    beside the canonical entry, and restricting to Swiss-Prot takes that to 0.9%, of which the
    remainder is dropped rather than guessed. 39,808 of 59,033 GTEx rows have no reviewed
    accession at all -- non-coding RNA and pseudogenes -- and are left out rather than carried
    as blanks.
    
    Validated on markers that must separate: albumin is 12,448 TPM in hepatocyte against 0.24 in
    fibroblast, APOA1 3,330 against 0.40, and in the other direction COL1A1 is 4,009 in
    fibroblast against 3.4 in hepatocyte and fibronectin 21,268 against 200. Three to four
    orders of magnitude the right way round in both directions is the check that the Ensembl
    keying is correct.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_gtex_transcriptome.py
"""
from _common import run

KEY = "host_gtex_transcriptome"

if __name__ == "__main__":
    run(KEY)
