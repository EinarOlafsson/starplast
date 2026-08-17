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


# --------------------------------------------------------------------------- expression
#: A second PlasmoDB report: transcript abundance per life stage, plus the polysomal fraction.
EXPRESSION_TABLE = "plasmodb_pf3d7_expression.tsv"

#: Header substring -> column. Matched on a substring because PlasmoDB's headers carry the study
#: name, the sample, the read-assignment rule and the dataset in one string, and the part that
#: identifies the sample is the middle.
EXPRESSION = (
    ("Su Seven Stages", "Ring", "expr_ring"),
    ("Su Seven Stages", "Early Trophozoite", "expr_early_trophozoite"),
    ("Su Seven Stages", "Late Trophozoite", "expr_late_trophozoite"),
    ("Su Seven Stages", "Schizont", "expr_schizont"),
    ("Su Seven Stages", "Gametocyte II", "expr_gametocyte_ii"),
    ("Su Seven Stages", "Gametocyte V", "expr_gametocyte_v"),
    ("Su Seven Stages", "Ookinete", "expr_ookinete"),
    # The polysomal fraction is what is ON ribosomes, which is a translation readout and not a
    # transcript level. Its steady-state partner from the same experiment is the transcript level,
    # and keeping both is the point: the pair is the only thing here that separates "more mRNA" from
    # "more translated".
    ("Polysomal and steady-state", "Polysomal ring", "polysomal_ring"),
    ("Polysomal and steady-state", "Polysomal troph", "polysomal_trophozoite"),
    ("Polysomal and steady-state", "Polysomal schiz", "polysomal_schizont"),
    ("Polysomal and steady-state", "Steady_state ring", "steady_state_ring"),
    ("Polysomal and steady-state", "Steady_state troph", "steady_state_trophozoite"),
    ("Polysomal and steady-state", "Steady_state schiz", "steady_state_schizont"),
    # Protein, and COMPOSITIONAL -- which is why these are not called `protein_ring`. The study is
    # TMT isobaric labelling of ring, trophozoite and schizont, and PlasmoDB serves the channels
    # row-normalised: every gene's three values sum to 12.07 +/- 0.20 and the three columns are
    # anti-correlated with each other (-0.46, -0.74, -0.14). So they say WHICH STAGE a protein sits
    # in, not how much of it there is, and they cannot answer a slot that asks for abundance. Named
    # `stage_share` so nobody has to rediscover that by wondering why protein disagrees with its own
    # transcript: the -0.25 correlation between ring protein and ring mRNA is the normalisation
    # showing through, not biology.
    ("Ring Ave", "", "protein_stage_share_ring"),
    ("Troph Ave", "", "protein_stage_share_trophozoite"),
    ("Schizont Ave", "", "protein_stage_share_schizont"),
    ("sense - asexual blood stages", "", "expr_asexual_blood"),
    ("sense - midgut oocysts", "", "expr_oocyst"),
    ("sense - salivary gland sporozoites", "", "expr_sporozoite"),
)


def expression(report_path: str) -> pd.DataFrame:
    """Per-stage transcript abundance and the polysomal fraction, keyed by gene."""
    if not os.path.exists(report_path):
        return pd.DataFrame()
    d = pd.read_csv(report_path, sep="\t", dtype=str)
    if "Gene ID" not in d.columns:
        return pd.DataFrame()
    d = d.replace(dict.fromkeys(BLANK, None))
    out = pd.DataFrame({"gene_id": d["Gene ID"].astype(str)})
    for study, sample, column in EXPRESSION:
        found = [c for c in d.columns if study in c and sample in c]
        if len(found) == 1:
            out[column] = pd.to_numeric(d[found[0]], errors="coerce")
    if len(out.columns) == 1:
        return pd.DataFrame()
    # Same transcript-versus-gene collapse as the attribute report, and the same reason.
    return out[~out["gene_id"].duplicated()].sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- export prediction
#: ExportPred calls a protein exported to the erythrocyte from its signal sequence and PEXEL motif.
#: PlasmoDB serves it as a SEARCH with a score threshold rather than as a per-gene attribute, so the
#: score is recovered by asking at several thresholds and recording the highest one a gene survives.
#: The scale saturates: asking for 20 returns nothing at all, so 10 -- the algorithm's own default --
#: is the top tier and not an arbitrary cut.
EXPORT_DIR = "exportpred"
EXPORT_TIERS = (1, 5, 10)


def export_prediction(folder: str) -> pd.DataFrame:
    """Predicted export to the host erythrocyte, as an ordinal tier.

    A tier rather than a boolean because the threshold matters and the default loses real biology:
    MESA and PfEMP3 are exported by any textbook and both fall below 10, while KAHRP and the FIKK
    kinases sit above it. Collapsing to the default would have called two of the best-known exported
    proteins in the organism not exported.

    Absence is a real negative here and not a gap. ExportPred is a sequence model, so it was
    evaluated on every protein; a gene missing from all three answers is predicted NOT exported,
    which is the opposite of the screen columns where absence means nobody looked.
    """
    if not os.path.isdir(folder):
        return pd.DataFrame()
    tiers = {}
    for rank, threshold in enumerate(EXPORT_TIERS, start=1):
        path = os.path.join(folder, f"exportpred_score_ge_{threshold}.tsv")
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t", dtype=str)
        if "Gene ID" not in d.columns:
            continue
        for gene in d["Gene ID"].astype(str):
            tiers[gene] = max(tiers.get(gene, 0), rank)
    if not tiers:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": sorted(tiers), "export_pred_tier": [tiers[g] for g in sorted(tiers)]})
    out["is_exported"] = out["export_pred_tier"] >= len(EXPORT_TIERS)
    return out


def build_all(dataset_root: str, log=print) -> pd.DataFrame:
    """The whole Plasmodium table from the three PlasmoDB reports, assembled in one place.

    Exists so the table is reproducible rather than the product of whatever was typed at a prompt,
    and so the one judgement call in the assembly is written down: a gene absent from ExportPred is
    recorded as tier 0 rather than as missing, because ExportPred is a sequence model evaluated on
    every protein and its silence is a prediction of "not exported". Every other absence in this
    table is ignorance and stays missing.
    """
    base = os.path.join(dataset_root, "reference", "plasmodb")
    nodes = build(os.path.join(base, "plasmodb_pf3d7_gene_attributes.tsv"))
    if nodes.empty:
        return nodes
    stages = expression(os.path.join(base, EXPRESSION_TABLE))
    if not stages.empty:
        nodes = nodes.merge(stages, on="gene_id", how="left")
    exported = export_prediction(os.path.join(base, EXPORT_DIR))
    if not exported.empty:
        nodes = nodes.merge(exported, on="gene_id", how="left")
        nodes["export_pred_tier"] = nodes["export_pred_tier"].fillna(0).astype(int)
        nodes["is_exported"] = nodes["export_pred_tier"] >= len(EXPORT_TIERS)
    log(f"Plasmodium table: {len(nodes):,} genes, {len(nodes.columns)} columns")
    return nodes
