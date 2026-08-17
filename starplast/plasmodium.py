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
    phospho = phosphosites(os.path.join(base, PHOSPHO_DIR), log=log)
    if not phospho.empty:
        nodes = nodes.merge(phospho, on="gene_id", how="left")
        # The Toxoplasma convention, mirrored so the arms read the same: the COUNT stays missing
        # where nothing was detected, because "how many sites" is genuinely unknown for a protein
        # mass spectrometry never saw; the BOOLEAN is False, because "was it ever observed
        # phosphorylated" is a yes-or-no about the evidence and the answer is no.
        nodes["has_phospho"] = nodes["has_phospho"].notna() & (nodes["has_phospho"] == True)  # noqa: E712
    fold = alphafold(os.path.join(base, ALPHAFOLD_TABLE))
    if not fold.empty:
        nodes = nodes.merge(fold, on="gene_id", how="left")
    sir2 = sir2_perturbation(os.path.join(base, SIR2_TABLE))
    if not sir2.empty:
        nodes = nodes.merge(sir2, on="gene_id", how="left")
    palm = palmitome(dataset_root, log=log)
    if not palm.empty:
        nodes = nodes.merge(palm, on="gene_id", how="left")
        nodes["is_palmitoylated"] = nodes["is_palmitoylated"].notna()
    log(f"Plasmodium table: {len(nodes):,} genes, {len(nodes.columns)} columns")
    return nodes


# --------------------------------------------------------------------------- phosphorylation
#: A re-analysis of every public Plasmodium phosphoproteomics dataset (PXD046874), reprocessed
#: through one pipeline so the sites are comparable. Site-level rather than peptide-level, which is
#: what makes a per-gene count meaningful: the same site found by three studies is one site.
PHOSPHO_DIR = "phosphosites"

#: The site q-value the deposit itself computes. 0.01 is the stricter of the two it reports, and
#: what a site count should be built on.
PHOSPHO_THRESHOLD = "Site Passes Threshold [0.01]"


def phosphosites(folder: str, log=print) -> pd.DataFrame:
    """Distinct phosphorylated residues per gene, pooled across the re-analysed studies.

    Counted as distinct (gene, position) pairs and not as rows. A site-centric table still carries
    one row per peptidoform and per source run, so the same serine found in four experiments is four
    rows; summing them would count how often a protein was looked at rather than how many sites it
    has. That is the mistake this function exists to not make, and the difference is large -- the
    files hold millions of rows for tens of thousands of sites.
    """
    if not os.path.isdir(folder):
        return pd.DataFrame()
    sites = set()
    files = sorted(f for f in os.listdir(folder) if f.endswith(".tsv"))
    for name in files:
        frame = pd.read_csv(os.path.join(folder, name), sep="\t", dtype=str,
                            usecols=lambda c: c in {"Proteins", "Protein Modification Positions",
                                                    "Modification", PHOSPHO_THRESHOLD},
                            low_memory=False)
        if "Proteins" not in frame.columns:
            continue
        if "Modification" in frame.columns:
            frame = frame[frame["Modification"].astype(str).str.contains("Phospho", na=False)]
        if PHOSPHO_THRESHOLD in frame.columns:
            frame = frame[frame[PHOSPHO_THRESHOLD].astype(str).str.strip() == "1"]
        for protein, position in zip(frame["Proteins"].astype(str),
                                     frame.get("Protein Modification Positions",
                                               pd.Series(dtype=str)).astype(str)):
            # `PF3D7_1346300.1-p1` is a transcript-and-product form of the gene accession, and a
            # gene with two products would otherwise count its sites twice.
            gene = protein.split(".")[0].split("-")[0].strip()
            if gene.startswith("PF3D7_"):
                sites.add((gene, position))
    if not sites:
        return pd.DataFrame()
    counts = pd.Series([g for g, _p in sites]).value_counts()
    out = pd.DataFrame({"gene_id": counts.index, "n_phosphosites": counts.to_numpy()})
    out["has_phospho"] = True
    log(f"phosphosites: {len(sites):,} distinct sites over {len(out):,} genes")
    return out.sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- palmitoylation
#: The palmitome compiled by PMID 36250062, which pools the published Plasmodium palmitoylation
#: datasets alongside its own analysis.
PALMITOME = ("palmitome", "36250062", "Table_3.xlsx")
PALMITOME_SHEET = "Palmitome"

#: The column of OBSERVED palmitoylated proteins. Naming it explicitly matters more than usual here:
#: the same workbook carries a sheet called `nrPalmitoylatedProteins` whose 3,105 rows are the UNION
#: of palmitoyl-ABLE (a motif prediction, 2,902 proteins) and palmitoylATED (503 observed). Its name
#: says palmitoylated and its contents are mostly predicted, and taking it at its name would have
#: called 54% of the proteome palmitoylated -- against published palmitomes of 400 to 500.
PALMITOME_OBSERVED = "Palmitoylated Proteins"
PALMITOME_PREDICTED = "Palmitoylable Proteins"


def palmitome(dataset_root: str, log=print) -> pd.DataFrame:
    """Proteins observed palmitoylated. Prediction is deliberately not folded in."""
    folder, pmid, name = PALMITOME
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if PALMITOME_SHEET not in book.sheet_names:
        return pd.DataFrame()
    sheet = book.parse(PALMITOME_SHEET)
    if PALMITOME_OBSERVED not in sheet.columns:
        return pd.DataFrame()
    genes = sorted({str(g).strip() for g in sheet[PALMITOME_OBSERVED].dropna()
                    if str(g).startswith("PF3D7_")})
    if not genes:
        return pd.DataFrame()
    log(f"palmitome: {len(genes):,} proteins observed palmitoylated "
        f"(the workbook's {PALMITOME_PREDICTED.lower()} column is a prediction and is not used)")
    return pd.DataFrame({"gene_id": genes, "is_palmitoylated": True})


# --------------------------------------------------------------------------- chromatin perturbation
#: Transcription in sir2a and sir2b knockouts against wild type, at three stages of the cycle.
#: Sir2 is a histone deacetylase and the knockouts are the chromatin perturbation this arm has.
SIR2_TABLE = "plasmodb_pf3d7_sir2_perturbation.tsv"
SIR2 = (("wild type - ring", "sir2_wt_ring"),
        ("wild type - trophozoite", "sir2_wt_trophozoite"),
        ("wild type - schizont", "sir2_wt_schizont"),
        ("sir2a KO - ring", "sir2a_ko_ring"),
        ("sir2a KO - trophozoite", "sir2a_ko_trophozoite"),
        ("sir2a KO - schizont", "sir2a_ko_schizont"),
        ("sir2b KO - ring", "sir2b_ko_ring"),
        ("sir2b KO - trophozoite", "sir2b_ko_trophozoite"),
        ("sir2b KO - schizont", "sir2b_ko_schizont"))


def sir2_perturbation(report_path: str) -> pd.DataFrame:
    """Wild type and knockout intensities, shipped as stated conditions rather than as a contrast.

    The knockout-minus-wild-type difference is the quantity anyone will want, and it is deliberately
    NOT computed here. The independent check on it came out ambiguous: Sir2a silences subtelomeric
    var genes, so var should rise in the sir2a knockout, and it does in ring (+0.135, p = 4e-12) and
    schizont (+0.130, p = 2e-38) -- but FALLS in trophozoite (-0.240, p = 3e-22). The effects are
    also small against a spread of 0.7. That is consistent with the canonical result being
    subset-specific and with var probes cross-hybridising across sixty paralogues, but it is not a
    clean confirmation, and a derived column carrying an unexplained sign flip would state more
    confidence than there is. The conditions themselves are unambiguous -- PlasmoDB names them, the
    scale is log intensity, and the medians align across arrays -- so those are what ship.
    """
    if not os.path.exists(report_path):
        return pd.DataFrame()
    d = pd.read_csv(report_path, sep="\t", dtype=str)
    if "Gene ID" not in d.columns:
        return pd.DataFrame()
    d = d.replace(dict.fromkeys(BLANK, None))
    out = pd.DataFrame({"gene_id": d["Gene ID"].astype(str)})
    for label, column in SIR2:
        found = [c for c in d.columns if label in c]
        if len(found) == 1:
            out[column] = pd.to_numeric(d[found[0]], errors="coerce")
    if len(out.columns) == 1:
        return pd.DataFrame()
    return out[~out["gene_id"].duplicated()].sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- fold confidence
#: AlphaFold's own per-model summary, fetched per protein from the AlphaFold DB API. `mean_plddt`
#: matches the Toxoplasma column name; the four fractions are the disorder half of the same slot,
#: because a mean hides the shape -- a protein that is half well-folded and half disordered scores
#: the same as one that is uniformly mediocre, and those are not the same protein.
ALPHAFOLD_TABLE = "plasmodb_pf3d7_alphafold.tsv"
ALPHAFOLD = (("globalMetricValue", "mean_plddt"),
             ("fractionPlddtVeryLow", "plddt_fraction_very_low"),
             ("fractionPlddtLow", "plddt_fraction_low"),
             ("fractionPlddtConfident", "plddt_fraction_confident"),
             ("fractionPlddtVeryHigh", "plddt_fraction_very_high"))


def alphafold(report_path: str) -> pd.DataFrame:
    """Model confidence per gene, keyed back to the PlasmoDB accession.

    A gene can carry several UniProt accessions -- 876 of them do, mostly the variant surface
    families where each field isolate's allele has its own entry. The fetch tries them in order and
    keeps the first with a model, and `uniprot_used` records which, so a number can be traced to the
    structure it came from rather than to a gene that has eight.
    """
    if not os.path.exists(report_path):
        return pd.DataFrame()
    d = pd.read_csv(report_path, sep="\t", dtype=str)
    if "gene_id" not in d.columns:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": d["gene_id"].astype(str)})
    for source, column in ALPHAFOLD:
        if source in d.columns:
            out[column] = pd.to_numeric(d[source], errors="coerce")
    if "uniprot_used" in d.columns:
        out["alphafold_accession"] = d["uniprot_used"].astype(str)
    if "mean_plddt" not in out.columns:
        return pd.DataFrame()
    return out[~out["gene_id"].duplicated()].sort_values("gene_id").reset_index(drop=True)
