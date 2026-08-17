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

import numpy as np
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
    iso = isoforms(dataset_root, log=log)
    if not iso.empty:
        # Left-joined without filling: a gene with no long-read model was not sequenced deeply
        # enough to say, which is not the same as having one transcript.
        nodes = nodes.merge(iso, on="gene_id", how="left")
    mapping = strain_map(os.path.join(base, NF54_TABLE), nodes)
    lac = lactylome(dataset_root, mapping, log=log) if mapping else pd.DataFrame()
    if not lac.empty:
        nodes = nodes.merge(lac, on="gene_id", how="left")
        nodes["has_lactyl"] = nodes["has_lactyl"].notna() & (nodes["has_lactyl"] == True)  # noqa: E712
    acetyl = acetylome(dataset_root, log=log)
    if not acetyl.empty:
        nodes = nodes.merge(acetyl, on="gene_id", how="left")
        # Same split as phosphorylation: the count is unknown where nothing was seen, the flag is no.
        nodes["has_acetyl"] = nodes["has_acetyl"].notna() & (nodes["has_acetyl"] == True)  # noqa: E712
    myr = myristoylome(dataset_root, log=log)
    if not myr.empty:
        # Left-joined and NOT filled: a gene outside the pulldown stays missing, because it was
        # never assayed. Only the 609 that were carry True or False.
        nodes = nodes.merge(myr, on="gene_id", how="left")
    palm = palmitome(dataset_root, log=log)
    if not palm.empty:
        nodes = nodes.merge(palm, on="gene_id", how="left")
        nodes["is_palmitoylated"] = nodes["is_palmitoylated"].notna()
    epitopes = bcell_epitopes(dataset_root, log=log)
    if not epitopes.empty:
        # Absent is absent: IEDB records what somebody tested.
        nodes = nodes.merge(epitopes, on="gene_id", how="left")
    derived = derived_labels(nodes, log=log)
    for column in derived.columns:
        nodes[column] = derived[column].to_numpy()
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


# --------------------------------------------------------------------------- N-myristoylation
#: Click-chemistry capture of myristoylated proteins, with and without an NMT inhibitor. The
#: quantity that identifies a substrate is not being pulled down -- background comes down too -- but
#: coming down LESS when the transferase is blocked.
MYRISTOYLOME = ("myristoylome", "34695132", "pbio.3001408.s011.xlsx")


def myristoylome(dataset_root: str, log=print) -> pd.DataFrame:
    """Proteins whose capture drops when N-myristoyltransferase is inhibited.

    Three states, not two, and the difference matters. A protein marked significant is a substrate;
    a protein in the table but not marked was assayed and is not one; a protein absent from the
    table was never in the pulldown at all and nothing is known about it. Collapsing the last two
    would claim the whole proteome had been tested for myristoylation by a single experiment that
    saw 609 proteins.
    """
    folder, pmid, name = MYRISTOYLOME
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    sheet = book.parse(book.sheet_names[0])
    if not {"Significant", "Difference", "Protein IDs"} <= set(sheet.columns):
        return pd.DataFrame()
    genes = sheet["Protein IDs"].astype(str).str.split(";").str[0].str.split(".").str[0]
    significant = sheet["Significant"].astype(str).str.strip() == "+"
    # Depletion on inhibition, not enrichment. A positive difference under a blocked transferase
    # would be a protein that came down MORE without it, which is not what a substrate does.
    depleted = pd.to_numeric(sheet["Difference"], errors="coerce") < 0
    out = pd.DataFrame({"gene_id": genes, "is_myristoylated": significant & depleted})
    out = out[out["gene_id"].str.startswith("PF3D7_")]
    if out.empty:
        return pd.DataFrame()
    out = out.groupby("gene_id", as_index=False)["is_myristoylated"].max()
    log(f"myristoylome: {int(out['is_myristoylated'].sum())} substrates of {len(out):,} assayed")
    return out.sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- acetylation
#: Proteome-wide lysine acetylation. The site list carries an Ascore -- the probability that the
#: acetyl group sits on the lysine named rather than a neighbouring one -- and it is not pre-filtered.
ACETYLOME = ("acetylome", "26813983", "srep19722-s2.xls")
ACETYLOME_SHEET = "1. Final Ac-K List w Motifs"

#: The site-count cut. Identifying an acetylated PEPTIDE is one claim and localising the acetyl group
#: to a particular lysine is a harder one, so the two columns are built to different standards: the
#: flag uses every identification, the count only sites localised at 0.75 or better.
ACETYL_LOCALISATION = 0.75


def acetylome(dataset_root: str, log=print) -> pd.DataFrame:
    """Acetylated lysines per gene, and whether the gene was seen acetylated at all."""
    folder, pmid, name = ACETYLOME
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if ACETYLOME_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(ACETYLOME_SHEET, header=1)
    if not {"Accession Number", "Site"} <= set(d.columns):
        return pd.DataFrame()
    d = d.assign(gene=d["Accession Number"].astype(str).str.strip())
    d = d[d["gene"].str.startswith("PF3D7_")]
    if d.empty:
        return pd.DataFrame()
    score = pd.to_numeric(d.get("Ascore Localization Probability"), errors="coerce")
    localised = d[score.fillna(0) >= ACETYL_LOCALISATION].drop_duplicates(["gene", "Site"])
    counts = localised.groupby("gene").size()
    out = pd.DataFrame({"gene_id": sorted(set(d["gene"]))})
    out["n_acetylsites"] = out["gene_id"].map(counts)
    out["has_acetyl"] = True
    log(f"acetylome: {len(out):,} genes seen acetylated, "
        f"{int(counts.sum()):,} sites localised at >= {ACETYL_LOCALISATION}")
    return out


# --------------------------------------------------------------------------- strain identity
#: NF54 is the line 3D7 was cloned from, and some studies report against its annotation instead.
#: Joining those to this table on the accession string would drop every one of them silently, which
#: is the failure the Toxoplasma identity layer exists to prevent.
NF54_TABLE = "plasmodb_nf54_orthogroups.tsv"

#: A pair whose two proteins differ in length by more than this is not the same gene, whatever the
#: orthogroup says, and is dropped rather than trusted.
STRAIN_LENGTH_TOLERANCE = 0.5


def strain_map(report_path: str, nodes: pd.DataFrame) -> dict:
    """NF54 accession -> 3D7 accession, for orthogroups holding exactly one gene on each side.

    Built from orthology and then CHECKED against protein length, because orthology is a claim about
    ancestry and this needs a claim about identity. The two coincide here: 96.8% of the pairs have
    exactly the same protein length and 99.2% are within 5%, which is what it should look like when
    one line was cloned from the other. Groups with more than one gene on either side are dropped
    rather than guessed at -- those are the paralogous surface families, where guessing would attach
    a measurement to the wrong member.
    """
    if not os.path.exists(report_path) or nodes.empty:
        return {}
    d = pd.read_csv(report_path, sep="\t", dtype=str).replace(dict.fromkeys(BLANK, None))
    if len(d.columns) < 3:
        return {}
    d.columns = ["nf54", "orthogroup", "length"][:len(d.columns)]
    d = d[d["orthogroup"].astype(str).str.startswith("OG6_")]
    d["length"] = pd.to_numeric(d["length"], errors="coerce")
    here = nodes[nodes["orthogroup"].astype(str).str.startswith("OG6_")]
    single_nf = set(d.groupby("orthogroup").size().pipe(lambda s: s[s == 1]).index)
    single_3d = set(here.groupby("orthogroup").size().pipe(lambda s: s[s == 1]).index)
    shared = single_nf & single_3d
    if not shared:
        return {}
    pairs = (d[d["orthogroup"].isin(shared)]
             .merge(here[here["orthogroup"].isin(shared)][["gene_id", "orthogroup", "length"]],
                    on="orthogroup", suffixes=("_nf", "_3d")))
    both = pairs.dropna(subset=["length_nf", "length_3d"])
    off = (both["length_nf"] - both["length_3d"]).abs() / both["length_3d"]
    keep = set(both.loc[off <= STRAIN_LENGTH_TOLERANCE, "nf54"]) | set(
        pairs.loc[pairs["length_nf"].isna() | pairs["length_3d"].isna(), "nf54"])
    return dict(zip(pairs.loc[pairs["nf54"].isin(keep), "nf54"],
                    pairs.loc[pairs["nf54"].isin(keep), "gene_id"]))


# --------------------------------------------------------------------------- lactylation
LACTYLOME = ("lactylome", "41417877", "pgen.1011991.s014.xlsx")
LACTYLOME_SHEET = "K(La)"
LACTYL_LOCALISATION = 0.75


def lactylome(dataset_root: str, mapping: dict, log=print) -> pd.DataFrame:
    """Lysine lactylation sites per gene, reported against NF54 and resolved to 3D7."""
    folder, pmid, name = LACTYLOME
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path) or not mapping:
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if LACTYLOME_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(LACTYLOME_SHEET)
    if "Gene ID" not in d.columns:
        return pd.DataFrame()
    genes = d["Gene ID"].astype(str).str.split(";").str[0].str.strip().map(mapping)
    score = pd.to_numeric(d.get("Localization prob"), errors="coerce")
    position = d.get("Position within protein ", d.get("Position within protein"))
    frame = pd.DataFrame({"gene_id": genes, "site": position, "score": score}).dropna(
        subset=["gene_id"])
    if frame.empty:
        return pd.DataFrame()
    localised = frame[frame["score"].fillna(0) >= LACTYL_LOCALISATION].drop_duplicates(
        ["gene_id", "site"])
    out = pd.DataFrame({"gene_id": sorted(set(frame["gene_id"]))})
    out["n_lactylsites"] = out["gene_id"].map(localised.groupby("gene_id").size())
    out["has_lactyl"] = True
    log(f"lactylome: {len(out)} genes resolved from NF54, "
        f"{int(out['n_lactylsites'].sum())} localised sites")
    return out


# --------------------------------------------------------------------------- isoforms
#: Long-read transcript models classified by SQANTI against the reference annotation.
ISOFORMS = ("isoforms", "40316999", "SupplementaryData2.xlsx")
ISOFORM_SHEET = "PF_sense"

#: SQANTI categories that mean the model is NOT the annotated transcript. `full-splice_match` is the
#: reference transcript recovered; everything here is a splice pattern the annotation does not have.
#: `novel_transcript_models` is named for the Toxoplasma column that answers the same question.
NOVEL_CATEGORIES = ("novel_in_catalog", "novel_not_in_catalog", "fusion")


def isoforms(dataset_root: str, log=print) -> pd.DataFrame:
    """Transcript models per gene, and how many of them the annotation does not contain."""
    folder, pmid, name = ISOFORMS
    path = os.path.join(dataset_root, "transcription", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if ISOFORM_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(ISOFORM_SHEET)
    if not {"associated_gene", "structural_category"} <= set(d.columns):
        return pd.DataFrame()
    # A fusion model names both genes it joins; the first is the one it is anchored on.
    gene = d["associated_gene"].astype(str).str.split("_novel").str[0].str.split("_").str[:2]
    gene = gene.str.join("_").str.strip()
    frame = pd.DataFrame({"gene_id": gene,
                          "category": d["structural_category"].astype(str).str.strip()})
    frame = frame[frame["gene_id"].str.match(r"^PF3D7_\w+$", na=False)]
    if frame.empty:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": sorted(set(frame["gene_id"]))})
    out["n_transcript_models"] = out["gene_id"].map(frame.groupby("gene_id").size())
    novel = frame[frame["category"].isin(NOVEL_CATEGORIES)]
    out["novel_transcript_models"] = out["gene_id"].map(
        novel.groupby("gene_id").size()).fillna(0).astype(int)
    log(f"isoforms: {len(out):,} genes, {int(out['n_transcript_models'].sum()):,} models, "
        f"{int(out['novel_transcript_models'].sum())} not in the annotation")
    return out


# --------------------------------------------------------------------------- derived stage labels
#: The life-cycle stages this arm can compare a gene across, and the columns that measure each.
#: Passed to `cellcycle.stage_enrichment` so both arms label stages by the same construction.
PF_STAGE_COLUMNS = {
    "ring": ("expr_ring",),
    "trophozoite": ("expr_early_trophozoite", "expr_late_trophozoite"),
    "schizont": ("expr_schizont",),
    "gametocyte": ("expr_gametocyte_ii", "expr_gametocyte_v"),
    "ookinete": ("expr_ookinete",),
    "oocyst": ("expr_oocyst",),
    "sporozoite": ("expr_sporozoite",),
}


def derived_labels(nodes: pd.DataFrame, log=print) -> pd.DataFrame:
    """The two DERIVED columns: peak expression across stages, and which stage a gene belongs to.

    Both are computed from the stage columns and declare it, so leakage closure can exclude them
    together with the measurements they come from. `expr_max` is the maximum on a log2 scale, matching
    the Toxoplasma column of the same name; the stage call and its margin come from the shared
    construction in `cellcycle`, which leaves a gene unlabelled unless one stage leads the next
    clearly -- a label that is really a coin toss looks like a measurement in every table it reaches.
    """
    from . import cellcycle
    present = [c for cols in PF_STAGE_COLUMNS.values() for c in cols if c in nodes.columns]
    if len(present) < 2:
        return pd.DataFrame(index=nodes.index)
    out = pd.DataFrame(index=nodes.index)
    out["expr_max"] = np.log2(nodes[present].max(axis=1) + 1)
    labels = cellcycle.stage_enrichment(nodes, log=log, stages=PF_STAGE_COLUMNS)
    for column in labels.columns:
        out[column] = labels[column]
    return out


# --------------------------------------------------------------------------- antibody epitopes
#: IEDB's antibody half, filtered to falciparum source antigens. `iedb` does the same for Toxoplasma
#: and reaches its genes through product descriptions, because IEDB's Toxoplasma antigen names are
#: verbatim ToxoDB descriptions. The falciparum names carry the UniProt accession instead, which is a
#: better key and is why this does not share that module's mapping.
UNIPROT_TABLE = "plasmodb_pf3d7_uniprot.tsv"

#: The two halves of IEDB, kept apart on purpose. `seroreactivity / antigenicity` is a question about
#: ANTIBODIES and `T-cell epitope content` is a question about T cells; answering either with the
#: other's number, or with the two pooled, would be answering a different question with a plausible
#: column. The numbers are not interchangeable either -- 434 antigens carry an antibody epitope and
#: only 44 carry a T-cell one.
IEDB_TABLES = {"n_bcell_epitopes": "iedb_pf_bcell_epitopes.tsv",
               "n_tcell_epitopes": "iedb_pf_tcell_epitopes.tsv"}
IEDB_TABLE = IEDB_TABLES["n_bcell_epitopes"]
IEDB_COLUMN = "n_bcell_epitopes"


def uniprot_index(report_path: str) -> dict:
    """UniProt accession -> gene, for accessions that name exactly one gene.

    209 accessions in PlasmoDB point at more than one gene. Those are dropped: an epitope belongs to
    a protein, and attaching it to whichever paralogue sorted first would be inventing the answer.
    """
    if not os.path.exists(report_path):
        return {}
    d = pd.read_csv(report_path, sep="\t", dtype=str).replace(dict.fromkeys(BLANK, None))
    if len(d.columns) < 2:
        return {}
    d.columns = ["gene_id", "uniprot"][:len(d.columns)]
    owners = {}
    for gene, cell in zip(d["gene_id"], d["uniprot"].fillna("")):
        for accession in str(cell).split(","):
            accession = accession.strip()
            if accession:
                owners.setdefault(accession, set()).add(gene)
    return {a: next(iter(g)) for a, g in owners.items() if len(g) == 1}


def bcell_epitopes(dataset_root: str, log=print) -> pd.DataFrame:
    """Distinct antibody epitope sequences per gene.

    DISTINCT sequences, not assay records: MSP1 alone carries 1,739 epitopes across 14,610 records,
    and counting records would rank antigens by how many groups have studied them rather than by how
    much of the protein antibodies recognise. Absent is absent and not zero -- IEDB records what
    somebody tested, and a protein nobody has raised an antibody against is not a protein without
    epitopes.
    """
    base = os.path.join(dataset_root, "reference", "plasmodb")
    owners = uniprot_index(os.path.join(base, UNIPROT_TABLE))
    if not owners:
        return pd.DataFrame()
    out = None
    for column, name in IEDB_TABLES.items():
        path = os.path.join(base, name)
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t", dtype=str)
        if not {"uniprot", "epitope"} <= set(d.columns):
            continue
        d = d.assign(gene_id=d["uniprot"].map(owners)).dropna(subset=["gene_id", "epitope"])
        if d.empty:
            continue
        counts = d.drop_duplicates(["gene_id", "epitope"]).groupby("gene_id").size()
        piece = pd.DataFrame({"gene_id": counts.index, column: counts.to_numpy()})
        out = piece if out is None else out.merge(piece, on="gene_id", how="outer")
        kind = "antibody" if "bcell" in name else "T-cell"
        log(f"{kind} epitopes (IEDB): {len(counts)} antigens, "
            f"{int(counts.sum()):,} distinct epitopes")
    return out if out is not None else pd.DataFrame()
