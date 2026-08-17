#!/usr/bin/env python3
"""The *Plasmodium falciparum* table: the second species, in its own table rather than merged.

Instruction 39 settled the shape and this is its first parasite table beyond Toxoplasma. Nothing is
merged. `PF3D7_*` and `TGME49_*` do not align, and neither does the data behind them -- piggyBac
essentiality is *falciparum*, PlasmoGEM barcoded knockouts are *berghei*, relapse biology is *vivax*.
Merged into one frame the missing values would encode which species was convenient to work on, and
since missingness is a feature here the embedding would learn laboratory history and it would look
like biology. Slots carry an organism, and a `Pf_` slot resolves against this table.

## Where it comes from

One PlasmoDB report, the same WDK service that answers six of the Toxoplasma slots. That single call
is worth more here than it is there, because PlasmoDB curates into gene attributes what for
Toxoplasma had to be found one paper at a time: the piggyBac saturation screen, the cross-strain SNP
counts, InterPro and orthology all arrive as columns of one table rather than as five acquisitions.

## The transcript-versus-gene trap

The report is served by the `transcript` record type, so a gene with two transcripts is two rows --
5,791 rows for 5,720 genes. Every quantity here is per gene, so the rows are collapsed and the
collapse is checked. Taking the report at face value would have silently double-weighted 71 genes.

## What is deliberately not read

`SignalP Peptide` covers 10.5% and reads as a label rather than a score; it is kept, but as a
boolean, and the absence of a call is not read as the absence of a signal peptide -- 89% of genes
having no entry is a statement about what was run, not about the proteins.
"""
from __future__ import annotations

import os

import pandas as pd

TABLE = "pf_nodes.parquet"

#: The organism this table is built for. PlasmoDB serves many; one table holds one.
ORGANISM = "Plasmodium falciparum 3D7"

#: PlasmoDB's report header -> the column name in this table. Names mirror the Toxoplasma table
#: where the quantity is the same, so a slot pattern reads the same on both arms.
COLUMNS = {
    "Gene ID": "gene_id",
    "Gene Type": "gene_type",
    "Product Description": "product",
    "Protein Length": "length",
    "Chromosome": "chromosome",
    "Transcript Length": "transcript_length",
    "# Exons in Transcript": "exon_count",
    "Molecular Weight": "molecular_weight",
    "Isoelectric Point": "isoelectric_point",
    "# TM Domains": "n_tm",
    "Ortholog count": "ortholog_number",
    "Paralog count": "paralog_number",
    "Ortholog Group": "orthogroup",
    "Total SNPs All Strains": "snp_total_all_strains",
    "NonSynonymous SNPs All Strains": "snp_nonsynonymous",
    "Synonymous SNPs All Strains": "snp_synonymous",
    "Non-Coding SNPs All Strains": "snp_noncoding",
    "SNPs with Stop Codons All Strains": "snp_stop_codon",
    "NonSyn/Syn SNP Ratio All Strains": "snp_nonsyn_syn_ratio",
    "P.falciparum 3D7 piggyBac insertion mutagenesis - mutant fitness score": "piggybac_mfs",
    "P.falciparum 3D7 piggyBac insertion mutagenesis - mutagenesis index score": "piggybac_mis",
    "Interpro ID": "interpro_ids",
    "PFam ID": "pfam_ids",
    "SignalP Peptide": "signalp",
}

#: Columns that are numbers. Everything else stays a string, and the two identifier lists stay
#: semicolon-joined text because a gene has as many domains as it has.
NUMERIC = ("length", "transcript_length", "exon_count", "molecular_weight", "isoelectric_point",
           "n_tm", "ortholog_number", "paralog_number", "snp_total_all_strains",
           "snp_nonsynonymous", "snp_synonymous", "snp_noncoding", "snp_stop_codon",
           "snp_nonsyn_syn_ratio", "piggybac_mfs", "piggybac_mis")

#: PlasmoDB writes a missing value three ways.
BLANK = ("N/A", "null", "")


def _collapse(frame: pd.DataFrame) -> pd.DataFrame:
    """One row per gene. The report is per transcript and every quantity here is per gene."""
    if "gene_id" not in frame.columns:
        return pd.DataFrame()
    # Longest transcript wins where a gene has more than one, so protein length and exon count
    # describe the same transcript rather than being mixed between rows.
    order = frame["transcript_length"] if "transcript_length" in frame.columns else frame.index
    frame = frame.assign(_order=pd.to_numeric(order, errors="coerce").fillna(0))
    frame = frame.sort_values("_order", ascending=False).drop(columns=["_order"])
    return frame[~frame["gene_id"].duplicated()].sort_values("gene_id").reset_index(drop=True)


def build(report_path: str) -> pd.DataFrame:
    """The Plasmodium table from one PlasmoDB attributes report."""
    if not os.path.exists(report_path):
        return pd.DataFrame()
    d = pd.read_csv(report_path, sep="\t", dtype=str)
    if "Gene ID" not in d.columns:
        return pd.DataFrame()
    d = d.replace(dict.fromkeys(BLANK, None))
    keep = {header: name for header, name in COLUMNS.items() if header in d.columns}
    out = d[list(keep)].rename(columns=keep)
    for column in NUMERIC:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    out = _collapse(out)
    if out.empty:
        return out
    if "interpro_ids" in out.columns:
        # A count and a flag, so a slot can ask "how much domain content" without parsing a list.
        out["n_interpro"] = out["interpro_ids"].fillna("").apply(
            lambda s: len([p for p in str(s).split(";") if p]))
        out["has_domain"] = out["n_interpro"] > 0
    if "signalp" in out.columns:
        # A call or nothing. Nothing means not called, which is not the same as no peptide.
        out["has_signal_peptide"] = out["signalp"].notna()
        out = out.drop(columns=["signalp"])
    if "n_tm" in out.columns:
        out["is_tm"] = out["n_tm"] > 0
    if "paralog_number" in out.columns:
        out["has_paralog"] = out["paralog_number"] > 0
    return out


def load(base: str) -> pd.DataFrame:
    """The shipped Plasmodium table, or an empty frame if it has not been built."""
    path = os.path.join(base, "starplast", "data", TABLE)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
