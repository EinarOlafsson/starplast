#!/usr/bin/env python3
"""How much a gene differs between sequenced Toxoplasma strains.

ToxoDB counts SNPs per gene across every strain it has high-throughput sequencing for, split by what
the substitution does to the protein. That split is the whole value of the table: a gene with many
synonymous and few nonsynonymous changes is under purifying selection, and a gene with the reverse is
under diversifying selection -- most likely because the immune system is looking at it.

## Why counts and not a diversity statistic

A pi or a Tajima's D would be a stronger summary and would need the alignments, which are not
published per gene. Counts are what ToxoDB serves, they are comparable across genes once length is
accounted for, and they answer the question the slot asks -- is this gene the same in every strain.

## What the numbers look like when they are right

Nonsynonymous SNPs per kb, measured on the shipped table: SRS surface antigens 112, the ROP5 / ROP18
/ GRA15 virulence loci 55, all genes 30, ribosomal proteins 2.9. That ordering is the check. The SRS
family is the most polymorphic thing in the genome because it is what host immunity sees, the
virulence loci are the classic strain-typing markers, and the ribosomal core is conserved. A load
that does not reproduce it has joined the wrong column or the wrong genes.
"""
from __future__ import annotations

import os

import pandas as pd

#: ToxoDB's display names on the left, because the report ships them as the header, and what the map
#: calls each column on the right. `All Strains` is ToxoDB's phrase for every strain it has
#: high-throughput sequencing for, which is what makes this a strain-variation measurement rather
#: than a comparison of two reference genomes.
COLUMNS = {
    "Total SNPs All Strains": "snp_total_all_strains",
    "NonSynonymous SNPs All Strains": "snp_nonsynonymous",
    "Synonymous SNPs All Strains": "snp_synonymous",
    "Non-Coding SNPs All Strains": "snp_noncoding",
    "SNPs with Stop Codons All Strains": "snp_stop_codon",
}

TABLE = "toxodb_strain_snps.tsv"


def strain_snps(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Per-gene SNP counts across sequenced strains, indexed by resolved gene id.

    Zero is a real measurement here and not a missing value: 690 genes carry no SNP in any sequenced
    strain, which is a statement about the gene rather than about coverage. Only genes ToxoDB did not
    report at all are absent.
    """
    path = os.path.join(base, "starplast", "data", TABLE)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    present = [c for c in COLUMNS if c in d.columns]
    if not present or d.columns[0] not in ("Gene ID", "gene_id"):
        return pd.DataFrame()
    genes = d[d.columns[0]].astype(str)
    if resolve is not None:
        genes = genes.map(lambda g: resolve(g) or g)
    out = d[present].apply(pd.to_numeric, errors="coerce").rename(columns=COLUMNS)
    out.index = pd.Index(genes, name="gene_id")
    out = out.groupby(level=0).max()
    log(f"strain variation (ToxoDB HTS SNPs): {out.shape[1]} columns, {len(out):,} genes")
    return out
