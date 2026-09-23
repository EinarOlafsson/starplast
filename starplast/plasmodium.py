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
import re

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
    # The minority strand from the same runs. Kept alongside its sense partner rather than as a ratio,
    # because the ratio's denominator is what makes it interpretable and a reader should see both.
    ("antisense - asexual blood stages", "", "antisense_asexual_blood"),
    ("antisense - midgut oocysts", "", "antisense_oocyst"),
    ("antisense - salivary gland sporozoites", "", "antisense_sporozoite"),
    ("sense - midgut oocysts", "", "expr_oocyst"),
    ("sense - salivary gland sporozoites", "", "expr_sporozoite"),
)


def _matches(column: str, label: str) -> bool:
    """Whether a PlasmoDB header names this sample, without matching a longer name that contains it.

    Plain containment is not enough and the failure is silent. `sense - asexual blood stages` is a
    SUBSTRING of `antisense - asexual blood stages`, so the moment the antisense columns were fetched
    each sense entry matched two headers, failed the one-match test, and three columns disappeared --
    adding antisense deleted sense. The label must start the header or be preceded by something that
    is not a letter.
    """
    at = column.find(label)
    if at < 0:
        return False
    return at == 0 or not column[at - 1].isalnum()


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
        found = [c for c in d.columns if _matches(c, study) and (not sample or sample in c)]
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


#: The Plasmodium CDS file, and the node table its CAI reference set is read from.
CDS_TABLE = "plasmodb_cds.tsv.gz"


def codon_usage(base: str, log=print) -> pd.DataFrame:
    """ENC, GC3 and CAI per gene, through the shared construction in `codons`.

    Not reimplemented. ENC and GC3 are definitions and the CAI reference set is "the ribosomal
    proteins" in both arms, so two implementations could only differ by being wrong in one of them.
    """
    from . import codons
    return codons.codon_usage(base, log=log, table=CDS_TABLE, nodes=TABLE)


def build_all(dataset_root: str, log=print, return_mentions: bool = False):
    """The whole Plasmodium table from the three PlasmoDB reports, assembled in one place.

    Exists so the table is reproducible rather than the product of whatever was typed at a prompt,
    and so the one judgement call in the assembly is written down: a gene absent from ExportPred is
    recorded as tier 0 rather than as missing, because ExportPred is a sequence model evaluated on
    every protein and its silence is a prediction of "not exported". Every other absence in this
    table is ignorance and stays missing.

    Assembly does not write package files. ``return_mentions=True`` returns
    ``(nodes, mentions)`` so a build command can publish both after its checks.
    """
    base = os.path.join(dataset_root, "reference", "plasmodb")
    nodes = build(os.path.join(base, "plasmodb_pf3d7_gene_attributes.tsv"))
    if nodes.empty:
        return (nodes,pd.DataFrame()) if return_mentions else nodes
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
    fever = febrile(os.path.join(base, FEBRILE_TABLE))
    if not fever.empty:
        nodes = nodes.merge(fever, on="gene_id", how="left")
    sir2 = sir2_perturbation(os.path.join(base, SIR2_TABLE))
    if not sir2.empty:
        nodes = nodes.merge(sir2, on="gene_id", how="left")
    iso = isoforms(dataset_root, log=log)
    if not iso.empty:
        # Left-joined without filling: a gene with no long-read model was not sequenced deeply
        # enough to say, which is not the same as having one transcript.
        nodes = nodes.merge(iso, on="gene_id", how="left")
    ribo = riboseq(dataset_root, log=log)
    if not ribo.empty:
        # Left-joined and never filled. Coverage is 61% and uneven across the cycle -- 2,182 genes
        # at the ring and 1,174 at the merozoite -- and every gap is a gene the deposit did not
        # report, which is silence about the measurement rather than an absence of ribosomes.
        nodes = nodes.merge(ribo, on="gene_id", how="left")
    partners = ip_ms(dataset_root, log=log)
    if not partners.empty:
        # Degree over the undirected pairs, so a bait counts its partners and a partner counts the
        # baits that pulled it down. Left missing for every gene outside the experiment: four
        # pulldowns are not a survey of the proteome, and a zero here would say "nothing binds this"
        # about a protein nobody tested.
        both = pd.concat([partners[["bait", "prey"]],
                          partners[["prey", "bait"]].rename(columns={"prey": "bait",
                                                                     "bait": "prey"})])
        degree = both.drop_duplicates().groupby("bait")["prey"].nunique()
        nodes["n_ip_ms_partners"] = nodes["gene_id"].map(degree)
    transferred = berghei_fitness(dataset_root, log=log)
    if not transferred.empty:
        # Left missing for a gene the berghei screen never carried: 2,448 of 5,720 have an ortholog
        # in it, and the rest are not dispensable, they are unscreened.
        nodes = nodes.merge(transferred, on="gene_id", how="left")
    similar = structure_similarity(dataset_root, log=log)
    if not similar.empty:
        # Degree over the undirected pairs. Zero where a gene has a model and no structural
        # neighbour -- which is a real answer, since the search looked -- and missing where the
        # proteome has no model for it at all.
        modelled = set(similar["gene_a"]) | set(similar["gene_b"])
        counted = pd.concat([similar[["gene_a", "gene_b"]],
                             similar[["gene_b", "gene_a"]].rename(
                                 columns={"gene_b": "gene_a", "gene_a": "gene_b"})])
        degree = counted.groupby("gene_a")["gene_b"].nunique()
        nodes["n_struct_similar"] = nodes["gene_id"].map(degree).where(
            nodes["gene_id"].isin(modelled) | nodes["alphafold_accession"].notna(), other=np.nan)
        nodes["n_struct_similar"] = nodes["n_struct_similar"].fillna(
            pd.Series(0, index=nodes.index).where(nodes["alphafold_accession"].notna()))
    liver = berghei_liver_fitness(dataset_root, log=log)
    if not liver.empty:
        nodes = nodes.merge(liver, on="gene_id", how="left")
    timing = idc_peak(dataset_root, log=log)
    if not timing.empty:
        # Left missing for a gene whose profile does not cycle: a flat series still has an angle,
        # and shipping it would put a number where there is no measurement.
        nodes = nodes.merge(timing, on="gene_id", how="left")
    kinase = kinase_substrates(dataset_root, log=log)
    if not kinase.empty:
        # Left-joined and not completed: a gene with no CDPK1-dependent site either has none or was
        # never quantified in those five replicates, and the file does not say which.
        nodes = nodes.merge(kinase, on="gene_id", how="left")
    vesicles = secretome(dataset_root, log=log)
    if not vesicles.empty:
        # Left-joined and deliberately NOT completed with a False flag, which is what the other
        # mass-spectrometry columns on this arm do. Those pool proteome-wide assays, where "never
        # observed" is an answer; this is one vesicle preparation from one isolate, and the genes
        # it did not report are six times less expressed than the ones it did. A False there would
        # be the absence-as-measurement this project exists to refuse -- the same fault as reading
        # hyperLOPIT `unassigned` as a compartment.
        nodes = nodes.merge(vesicles, on="gene_id", how="left")
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
    cplx = complexes(dataset_root, log=log)
    if not cplx.empty:
        nodes = nodes.merge(cplx, on="gene_id", how="left")
    from . import pf_graph
    degree = pf_graph.host_degree(nodes, dataset_root, log=log)
    if len(degree):
        nodes["n_host_targets"] = degree.to_numpy()
    ec = enzyme_classification(dataset_root, log=log)
    if not ec.empty:
        nodes = nodes.merge(ec, on="gene_id", how="left")
        nodes["has_ec"] = nodes["has_ec"].notna() & (nodes["has_ec"] == True)  # noqa: E712
    columns, _edges, mentions = literature_layer(nodes, os.path.dirname(os.path.abspath(dataset_root)),
                                                log=log)
    if not columns.empty:
        nodes = nodes.merge(columns, on="gene_id", how="left")
    codons_here = codon_usage(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), log=log)
    if not codons_here.empty:
        nodes = nodes.merge(codons_here, left_on="gene_id", right_index=True, how="left")
    derived = derived_labels(nodes, log=log)
    for column in derived.columns:
        nodes[column] = derived[column].to_numpy()
    log(f"Plasmodium table: {len(nodes):,} genes, {len(nodes.columns)} columns")
    return (nodes,mentions) if return_mentions else nodes


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


# --------------------------------------------------------------------------- gene identity
#: PlasmoDB's identity table: current accession, symbol, previous accessions, product. Committed
#: beside the ToxoDB ones and never merged with them -- see `fetch_names`.
IDENTITY_TABLE = "plasmodb_identity.tsv"


def previous_id_index(path: str) -> dict:
    """{any accession this genome has used, lowercased} -> the current `PF3D7_` id.

    The Plasmodium arm went without an identity layer for as long as every source happened to be
    keyed on current accessions, and the first one that was not joined nothing at all: a 2014
    ribosome-profiling deposit reports `PFE0630c` and `PF13_0222`, the chromosome-based ids that
    were current before the 2012 renaming, and the string join found zero of its 3,629 genes.

    Two rules, both taken from the Toxoplasma layer rather than reinvented:

    * **Ambiguity is withdrawn, never guessed.** 66 old ids are claimed by two current genes each,
      which is what a gene model being SPLIT looks like from the other side; a measurement made
      against the old model belongs to neither half in particular, so the string resolves to
      nothing and the cost is countable rather than invisible.
    * **A current accession outranks a previous one.** A file mixing both -- most do -- must not
      have a live id captured by some other gene's history.
    """
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path, sep="\t", dtype=str)
    if "gene_id" not in d.columns:
        return {}
    owners: dict = {}
    for gene, ids in zip(d["gene_id"], d.get("previous_ids", pd.Series([""] * len(d))).fillna("")):
        for token in re.split(r"[;,]", str(ids).replace("Previous IDs:", "")):
            token = token.strip().lower()
            if token and token != "n/a":
                owners.setdefault(token, set()).add(gene)
    index = {token: next(iter(genes)) for token, genes in owners.items() if len(genes) == 1}
    index.update({str(g).lower(): g for g in d["gene_id"]})
    return index


# --------------------------------------------------------------------------- ribosome profiling
#: Where the deposit lands under the dataset archive, and the five points of the cycle it covers.
#: The numbers are the deposit's own sample order; the names are the stages its sample records give.
RIBOSEQ = ("riboseq", "25493618")
RIBOSEQ_STAGES = {"1": "ring", "2": "early_trophozoite", "3": "late_trophozoite",
                  "4": "schizont", "5": "merozoite"}
RIBOSEQ_ARMS = {"mRNA": "mrna", "ribosome_footprints": "rpf"}
RIBOSEQ_FILE = re.compile(r"^GSM\d+_(mRNA|ribosome_footprints)_([1-5])_rpkm\.txt\.gz$")


def riboseq(dataset_root: str, log=print) -> pd.DataFrame:
    """Ribosome-footprint and mRNA density per gene at five points of the asexual cycle.

    The measurement `Pf_translation - per cell-cycle phase` asks for: how much ribosome is on a
    transcript, per gene, at each stage of the intraerythrocytic cycle. Both arms of the experiment
    are read and both are kept as CONDITIONS. The ratio between them -- translation efficiency -- is
    deliberately NOT computed here; `RIBOSEQ_TE_REFUSED` in the registry note says what happened
    when it was.

    Three things are refused rather than resolved, and together they cost about 100 of 3,629 ids:

    * an id claimed by two current genes (`previous_id_index` withdraws it);
    * the deposit's `-a` / `-b` split entries, which report two segments of one gene -- RPKM is
      already length-normalised, so neither summing nor averaging them means anything;
    * two source ids landing on one current gene inside one file, which is a MERGE seen from the
      other side and has the same problem.
    """
    folder, pmid = RIBOSEQ
    base = os.path.join(dataset_root, "translation", folder, pmid)
    if not os.path.isdir(base):
        return pd.DataFrame()
    from . import paths
    index = previous_id_index(paths.cache_file(IDENTITY_TABLE))
    if not index:
        log("riboseq: no PlasmoDB identity table -- run python -m starplast.fetch_names")
        return pd.DataFrame()
    columns, unresolved, collided = {}, set(), 0
    for name in sorted(os.listdir(base)):
        hit = RIBOSEQ_FILE.match(name)
        if not hit:
            continue
        arm, stage = RIBOSEQ_ARMS[hit.group(1)], RIBOSEQ_STAGES[hit.group(2)]
        d = pd.read_csv(os.path.join(base, name), sep="\t", header=None,
                        names=["source_id", "rpkm"], dtype={"source_id": str})
        d["source_id"] = d["source_id"].str.strip().str.lower()
        d["gene_id"] = d["source_id"].map(index)
        unresolved |= set(d.loc[d["gene_id"].isna(), "source_id"])
        d = d.dropna(subset=["gene_id"])
        duplicated = d["gene_id"].duplicated(keep=False)
        collided += int(duplicated.sum())
        d = d[~duplicated]
        columns[f"riboseq_{arm}_{stage}"] = pd.to_numeric(
            d.set_index("gene_id")["rpkm"], errors="coerce")
    if not columns:
        return pd.DataFrame()
    out = pd.DataFrame(columns).sort_index()
    out.index.name = "gene_id"
    log(f"riboseq: {len(out):,} genes over {len(columns)} arms x stages, "
        f"{len(unresolved)} ids unresolved, {collided} withdrawn for landing on one gene twice")
    return out.reset_index()


# --------------------------------------------------------------------------- secreted proteins
#: Extracellular vesicles purified from a Kenyan clinical isolate, and the same paper's compilation
#: of which proteins a second EV proteome also reported. The compilation is the sheet worth reading:
#: it is the union of two independent preparations with a membership column each, so "how many
#: studies saw this" is a fact in the file rather than a join someone has to get right.
SECRETOME = ("secretome", "28944300", "5c097a1c-efe5-4ed8-b97b-f9ba656268a6.xlsx")
SECRETOME_SHEET = "Combined list of PfEVs antigens"

#: The two EV proteomes the sheet crosses. Its other columns are seroreactivity and antibody-array
#: results from unrelated studies -- claims about immunity, not about vesicles -- and they belong to
#: other slots, so they are deliberately not read here.
SECRETOME_STUDIES = ("PfEVs_9605", "PfEVs_Mantel CHM 2013")


def secretome(dataset_root: str, log=print) -> pd.DataFrame:
    """Parasite proteins found in extracellular vesicles, and how many EV proteomes found them.

    The slot asks what this parasite puts outside itself. Vesicle proteomics is one instrument for
    that question and it is the one with data, so the column says `extracellular vesicle` rather
    than `secreted`: the second word would assert a route this measurement does not establish.

    ONE column, and no companion flag. A boolean completed with False for every other gene would
    have read as 5,720 genes tested and 5,536 negative, and the atlas would have graded the slot A
    at 100% coverage -- for an experiment that identified 184 proteins. The count being present IS
    the flag, and its absence is unknown.

    The sheet's trailing rows are its reference list, which is why resolution goes through the
    identity index rather than a `PF3D7_` regex -- the citations contain accessions.
    """
    folder, pmid, name = SECRETOME
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if SECRETOME_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(SECRETOME_SHEET)
    if "geneid" not in d.columns or not set(SECRETOME_STUDIES) <= set(d.columns):
        return pd.DataFrame()
    from . import paths
    index = previous_id_index(paths.cache_file(IDENTITY_TABLE))
    if not index:
        log("secretome: no PlasmoDB identity table -- run python -m starplast.fetch_names")
        return pd.DataFrame()
    d = d.assign(gene_id=d["geneid"].astype(str).str.strip().str.lower().map(index))
    d = d.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    if d.empty:
        return pd.DataFrame()
    seen = d[list(SECRETOME_STUDIES)].apply(pd.to_numeric, errors="coerce").fillna(0)
    out = pd.DataFrame({"gene_id": d["gene_id"].to_numpy(),
                        "ev_studies": seen.gt(0).sum(axis=1).astype(int).to_numpy()})
    out = out[out["ev_studies"] > 0]
    log(f"secretome: {len(out):,} proteins in extracellular vesicles, "
        f"{int((out['ev_studies'] > 1).sum())} of them in both preparations")
    return out.sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- kinase substrates
#: Sites that LOSE phosphorylation when PfCDPK1 is knocked down. The Toxoplasma arm answers the same
#: question from the other direction, with thiophosphorylation labelling a kinase's own substrates.
KINASE_SUBSTRATE = ("kinase_substrate", "28680058", "41467_2017_53_MOESM3_ESM.xls")
KINASE_SUBSTRATE_SHEET = "List of hypophosphyrlated sites"

#: The header sits on the third row, under the sheet's own title.
KINASE_SUBSTRATE_HEADER = 2

#: The window the deposit gives around each site, and the only usable key in the file.
KINASE_SUBSTRATE_WINDOW = "Phosphowindow"


def kinase_substrates(dataset_root: str, log=print) -> pd.DataFrame:
    """Phosphosites per gene that depend on PfCDPK1, keyed by SEQUENCE because nothing else works.

    The supplement reports sites as numeric ids from the 2017 annotation (`3885720(S422)`), which
    neither the current accessions nor PlasmoDB's previous-id list carry. What it does give is the
    15-residue window around each site, and a window is an identifier when it occurs in exactly one
    protein: 73 of 79 match one gene, 6 match none, and none matches two. A window matching more
    than one protein is dropped rather than assigned, the same rule the accession index follows.

    The column is named for DEPENDENCE, not for substrate. A site that loses phosphorylation when a
    kinase is knocked down may be phosphorylated by that kinase or by something downstream of it,
    and the file cannot tell the two apart.
    """
    folder, pmid, name = KINASE_SUBSTRATE
    path = os.path.join(dataset_root, "post_translation", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if KINASE_SUBSTRATE_SHEET not in book.sheet_names:
        return pd.DataFrame()
    try:
        d = book.parse(KINASE_SUBSTRATE_SHEET, header=KINASE_SUBSTRATE_HEADER)
    except (IndexError, ValueError, StopIteration):
        # The header row is at a fixed offset under the sheet's own title, so a re-issued supplement
        # with a different layout must report rather than take the build down with it.
        log("kinase substrates: the sheet is not laid out as this loader expects")
        return pd.DataFrame()
    if KINASE_SUBSTRATE_WINDOW not in d.columns:
        return pd.DataFrame()
    from . import codons, paths
    cds_path = paths.cache_file(CDS_TABLE)
    if not os.path.exists(cds_path):
        log("kinase substrates: no CDS table, so a window cannot be resolved to a gene")
        return pd.DataFrame()
    cds = pd.read_csv(cds_path, sep="\t")
    proteins = {gene: codons.translate(seq) for gene, seq in zip(cds["gene_id"], cds["cds"])}
    counts: dict = {}
    unmatched = ambiguous = 0
    for window in d[KINASE_SUBSTRATE_WINDOW].dropna():
        window = str(window).strip().upper()
        owners = [gene for gene, seq in proteins.items() if window in seq]
        if len(owners) != 1:
            unmatched += len(owners) == 0
            ambiguous += len(owners) > 1
            continue
        counts[owners[0]] = counts.get(owners[0], 0) + 1
    if not counts:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": sorted(counts)})
    out["cdpk1_dependent_sites"] = out["gene_id"].map(counts)
    log(f"kinase substrates: {len(out):,} genes carry {sum(counts.values())} CDPK1-dependent sites; "
        f"{unmatched} windows matched no protein and {ambiguous} matched more than one")
    return out


# --------------------------------------------------------------------------- measured binding
#: The EPIC interactome: four co-immunoprecipitations against their own controls, published as pages
#: of a supplementary PDF rather than as a table. `(folder, pmid, file)`.
IP_MS = ("ip_ms", "28691708", "ncomms16044-s1.pdf")

#: The accession of each bait that HAS one. `PfEMP1B` is a transgene built from a var gene and
#: appears in no row, so its pulldown cannot become a parasite-parasite pair -- see the loader.
IP_MS_BAITS = {"PV1-HA": "PF3D7_1129100", "PV2-HA": "PF3D7_1226900",
               "EXP3-HA": "PF3D7_1024800", "PfEMP1B": None}

#: A table starts on the page whose column header names its bait and its control; a page carrying
#: only the caption continues it. Both are read, and they have to agree.
IP_MS_HEADER = re.compile(r"(PfEMP1B|PV1-HA|PV2-HA|EXP3-HA)\s+(?:PfEMP1F|WT)")
IP_MS_CAPTION = re.compile(r"Supplementary Table \d\s*\|\s*Mass spectrometry analysis of "
                           r"(PfEMP1B|PV1-HA|PV2-HA|EXP3-HA)")
#: Anchored to the END of a record. Annotations carry numbers -- `exported protein 3`, `HSP70-2`,
#: `Pfj2` -- so matching the FIRST run of digits reads an annotation's own number as a count and
#: shifts every column by one, silently. Anchored, the eight counts are the last eight numbers in
#: the record, and a row carrying seven is refused rather than padded from its own name.
#:
#: A stricter version of this refused `heat shock protein DnaJ homologue, Pfj2`, whose annotation
#: legitimately ends in a digit -- one silent corruption traded for one honest loss, which is the
#: wrong trade when the anchor alone gets both right.
IP_MS_ROW = re.compile(r"^\s*(.*?)\s*((?:\d+\s+){7}\d+)\s*([<>]?[\d.]+)?\s*$")

#: The last row on a page is followed by the table's false-discovery-rate line, which is not part of
#: the record and would otherwise stop it anchoring to the end.
IP_MS_TAIL = re.compile(r"False discovery rate")


def ip_ms(dataset_root: str, log=print) -> pd.DataFrame:
    """Co-immunoprecipitation pairs: which parasite proteins came down with each tagged bait.

    The source is a PDF, and its CAPTIONS are wrong -- the table headed *parasite interacting
    proteins* holds human ones and the table headed *human* holds parasite ones. So nothing here
    reads a caption for what a table contains. The bait comes from the column header, which names
    the pulldown and its control, and the organism comes from the identifier space: `PF3D7_` rows are
    parasite, `*_HUMAN` rows are host and are not read here.

    **The check that licenses the whole assignment: a bait must be the first row of its own table.**
    PV1 tops the PV1 pulldown at 97 and 120 spectra, PV2 tops PV2's, EXP3 tops EXP3's. A page whose
    first row is not its bait is a CONTINUATION page, which is a different claim and is recorded as
    one; a page that starts a table and fails the check is refused, because the alternative is
    attributing one protein's partners to another.

    Two things are deliberately not represented. `PfEMP1B` is a transgene made from a var gene and
    has no accession in the table, so its 38 rows cannot become parasite-parasite pairs and are
    counted rather than guessed at. And a row carrying seven counts instead of eight (PIESP2, in the
    PV2 pulldown) is dropped: the missing number could be either arm.
    """
    folder, pmid, name = IP_MS
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        from pypdf import PdfReader
    except ImportError:                  # pragma: no cover - declared in pyproject, guarded anyway
        log("ip_ms: pypdf is not installed, so the interactome PDF cannot be read")
        return pd.DataFrame()
    reader = PdfReader(path)
    rows, bait, started, unnamed, malformed = [], None, set(), 0, 0
    for page in reader.pages:
        text = " ".join(page.extract_text().split())
        parts = re.split(r"(PF3D7_\w+)", text)
        if len(parts) < 3:
            continue
        head, cap = IP_MS_HEADER.search(text), IP_MS_CAPTION.search(text)
        if head:
            bait = head.group(1)
        elif cap:
            bait = cap.group(1)
        if bait is None:
            continue
        found = []
        for accession, body in zip(parts[1::2], parts[2::2]):
            hit = IP_MS_ROW.match(IP_MS_TAIL.split(body)[0])
            if not hit:
                malformed += 1
                continue
            counts = [int(n) for n in hit.group(2).split()]
            # The eight counts are TWO experiments of four: bait rep 1, bait rep 2, control rep 1,
            # control rep 2. Splitting them down the middle sums experiment 1 against experiment 2
            # instead of bait against control -- which reads as an "enriched" partner with 145
            # spectra in the untagged line, and is only obvious because that is impossible.
            bait_counts = counts[0] + counts[1] + counts[4] + counts[5]
            control_counts = counts[2] + counts[3] + counts[6] + counts[7]
            found.append((accession, hit.group(1).strip(), bait_counts, control_counts))
        if not found:
            continue
        own = IP_MS_BAITS[bait]
        if bait not in started:
            started.add(bait)
            if own is not None and found[0][0] != own:
                log(f"ip_ms: {bait}'s table does not start with {own}, so it is refused")
                continue
        if own is None:
            unnamed += len(found)
            continue
        for accession, product, in_bait, in_control in found:
            if accession != own:
                rows.append({"bait": own, "prey": accession, "product": product,
                             "spectra_bait": in_bait, "spectra_control": in_control})
    out = (pd.DataFrame(rows).drop_duplicates(["bait", "prey"]).reset_index(drop=True)
           if rows else pd.DataFrame())
    # Logged even when nothing was paired, because "no pairs" and "no pairs, and here is what was
    # read instead" are different reports and only the second one can be acted on.
    log(f"ip_ms: {len(out):,} pairs over {out.bait.nunique() if len(out) else 0} baits and "
        f"{out.prey.nunique() if len(out) else 0} partners; {unnamed} rows from the transgene bait "
        f"have no accession to pair, {malformed} rows carried the wrong number of counts")
    return out


# --------------------------------------------------------------------------- the cell cycle
#: A 48-hour intraerythrocytic time course sampled every three hours. `(folder, pmid)`; the deposit's
#: own metadata file sits beside the per-sample tables and is what says which sample is which.
IDC = ("idc_timecourse", "34668757")
IDC_METADATA = "GSE163144_family.soft"

#: The arm to read. The experiment crosses two haemoglobin genotypes with two parasite lines, and
#: only one cell of that square is an unperturbed reference: the laboratory strain in normal red
#: cells. Sickle-trait cells are the study's variable, and FUP is a different parasite.
IDC_GENOTYPE = "HbAA"
IDC_STRAIN = "3D7"


def idc_samples(path: str) -> dict:
    """{GSM: hours post-invasion} for the reference arm of the time course, from the deposit's own
    metadata rather than from the sample titles, which encode the genotype and not the hour."""
    if not os.path.exists(path):
        return {}
    out = {}
    for block in open(path, encoding="utf8", errors="replace").read().split("^SAMPLE = ")[1:]:
        gsm = block.split("\n", 1)[0].strip()
        traits = dict(re.findall(r"!Sample_characteristics_ch1 = ([^:]+): (.+)", block))
        hour = traits.get("hours post-invasion", "N/A").strip()
        if (traits.get("hemoglobin genotype") == IDC_GENOTYPE
                and traits.get("strain") == IDC_STRAIN and hour not in ("", "N/A")):
            out[gsm] = float(hour)
    return out


def idc_peak(dataset_root: str, log=print) -> pd.DataFrame:
    """When in the 48-hour cycle each gene's transcript peaks, from a three-hourly series.

    The slot asks for a cell-cycle timing label, and this is the measurement behind one: sixteen
    timepoints three hours apart, two replicates, one parasite line in unmodified red cells. The
    peak hour is the timepoint of the highest mean expression, which is a summary of a measurement
    rather than a model of one -- and it is left MISSING where a gene never rises, since the hour of
    a flat profile is noise wearing a number.
    """
    folder, pmid = IDC
    base = os.path.join(dataset_root, "transcription", folder, pmid)
    hours = idc_samples(os.path.join(base, IDC_METADATA))
    if not hours:
        return pd.DataFrame()
    columns = {}
    for name in sorted(os.listdir(base)):
        gsm = name.split("_", 1)[0]
        if gsm not in hours or not name.endswith(".csv.gz"):
            continue
        table = pd.read_csv(os.path.join(base, name))
        if table.shape[1] < 2:
            continue
        table.columns = ["gene_id"] + list(table.columns[1:])
        series = pd.to_numeric(table.iloc[:, 1], errors="coerce")
        columns.setdefault(hours[gsm], []).append(series.set_axis(table["gene_id"].astype(str)))
    if len(columns) < 8:
        log(f"idc: only {len(columns)} timepoints found, which is not a cycle")
        return pd.DataFrame()
    means = pd.DataFrame({hour: pd.concat(reps, axis=1).mean(axis=1)
                          for hour, reps in sorted(columns.items())})
    from . import cellcycle
    fitted = cellcycle.cyclic_phase(means, period=IDC_PERIOD)
    timed = fitted["amplitude"] >= IDC_MIN_AMPLITUDE
    out = pd.DataFrame({"gene_id": means.index[timed],
                        "idc_peak_hour": fitted.loc[timed, "phase"].to_numpy().round(1),
                        "idc_cycling_amplitude": fitted.loc[timed, "amplitude"].to_numpy()})
    log(f"idc: {len(out):,} genes are timed across {len(means.columns)} timepoints; "
        f"{int((~timed).sum())} do not cycle strongly enough to place")
    return out.reset_index(drop=True)


#: The cycle's length, which is what makes the axis wrap: hour 48 is hour 0 of the next round.
IDC_PERIOD = 48.0

#: How much of a gene's variation the first harmonic has to explain before its phase means anything.
#: A profile that does not cycle still has an angle, and shipping one would be a number where there
#: is no measurement.
IDC_MIN_AMPLITUDE = 0.4


# --------------------------------------------------------------------------- transferred fitness
#: The PlasmoGEM barcoded-knockout screen of *P. berghei*, whose own table carries the *falciparum*
#: ortholog per row. `(folder, pmid, file, sheet)`.
PB_TRANSFER = ("pb_transfer", "28708996", "mmc1.xlsx", "Table S1")

#: The screen's own verdict on a gene it could not call. Left missing rather than read as a middle
#: value: "insufficient data" is the absence of a measurement, not a slow-growth phenotype.
PB_UNCALLED = "Insufficient data"


def berghei_fitness(dataset_root: str, log=print) -> pd.DataFrame:
    """*P. berghei* knockout fitness, carried onto *falciparum* by the SOURCE's own ortholog column.

    This is a transfer, and instruction 39 requires transfers to be visible rather than folded into
    the measured slot -- so the columns say `pb_` and the slot says `transferred from Pb`. What makes
    this one safe is that the orthology is not mine: the screen's table names a *falciparum* gene per
    row, and no *falciparum* gene is named by two *berghei* ones, so the join is one-to-one and
    nothing has to be dropped for ambiguity.

    It is checked against the receiving arm's OWN screen, which is the check a transfer has to pass:
    genes the *berghei* screen calls essential have a median piggyBac mutagenesis index of 0.160,
    slow-growing ones 0.394 and dispensable ones 0.996 -- a monotonic ordering across two species and
    two unrelated methods (barcoded knockouts in mice against saturation mutagenesis in culture),
    with 65 of 71 ribosomal proteins essential. Backwards, it would be refused.
    """
    folder, pmid, name, sheet = PB_TRANSFER
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if sheet not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(sheet)
    needed = {"P. falciparum ID", "Phenotype", "Relative growth rate", "Confidence"}
    if not needed <= set(d.columns):
        return pd.DataFrame()
    # Forty rows name a TRANSCRIPT (`PF3D7_0108400.1`). The suffix is stripped, as the phosphosite
    # loader already does for the same source of ids -- kept whole, those forty genes would have been
    # dropped for not matching an accession pattern, which is a silent loss dressed as strictness.
    d = d.assign(gene_id=d["P. falciparum ID"].astype(str).str.strip()
                 .str.replace(r"\.\d+$", "", regex=True))
    d = d[d["gene_id"].str.match(r"^PF3D7_\w+$", na=False)]
    if d.empty:
        return pd.DataFrame()
    ambiguous = d["gene_id"].duplicated(keep=False)
    if ambiguous.any():
        log(f"pb transfer: {int(ambiguous.sum())} falciparum genes named by more than one berghei "
            f"gene, withdrawn")
        d = d[~ambiguous]
    uncalled = d["Phenotype"].astype(str).str.strip() == PB_UNCALLED
    out = pd.DataFrame({
        "gene_id": d["gene_id"].to_numpy(),
        "pb_transferred_phenotype": d["Phenotype"].where(~uncalled).to_numpy(),
        "pb_transferred_growth_rate": pd.to_numeric(
            d["Relative growth rate"], errors="coerce").where(~uncalled).to_numpy(),
        "pb_transfer_confidence": pd.to_numeric(d["Confidence"], errors="coerce").to_numpy()})
    log(f"pb transfer: {len(out):,} falciparum genes carry a berghei knockout phenotype; "
        f"{int(uncalled.sum())} the screen could not call are left missing")
    return out.reset_index(drop=True)


#: The liver-stage and transmission barcode screen, whose table is Pb-keyed and carries no
#: falciparum column -- so the mapping comes from the blood-stage screen above, which is the same
#: consortium's own pairing rather than an orthology this project derived.
PB_LIVER = ("pb_transfer", "31730853", "mmc2.xlsx", "Sheet1")

#: Column offsets in that sheet, which has two header rows and repeats `Log2-FC / SD / Power` per
#: transition. Read by POSITION because the second header row gives every group the same three
#: names; the first row is what says which transition a group is.
PB_LIVER_COLUMNS = {"pb_liver_log2fc": 21, "pb_liver_power": 23}

#: The transition that spans the liver: salivary gland to the second blood infection, corrected for
#: blood-stage fitness. Uncorrected, a gene needed in blood looks liver-essential because the
#: transition ENDS in blood.
PB_LIVER_TRANSITION = "SG-B2 data, normalized, BS fitness-corrected"

#: `no power` means the barcodes were too few to say anything, which is not a measurement of zero
#: effect. Those rows keep no value.
PB_NO_POWER = "no power"


def berghei_liver_fitness(dataset_root: str, log=print) -> pd.DataFrame:
    """Liver-stage fitness of *berghei* knockouts, carried onto *falciparum* orthologs.

    Two things make this safe, and both are borrowed rather than invented. The falciparum id comes
    from the blood-stage screen's own table -- the same consortium pairing the same mutants -- and
    the value is the SG-B2 transition corrected for blood-stage fitness by the authors, because the
    transition ends in blood and an uncorrected drop would call every blood-essential gene
    liver-essential.

    Validated on the genes the field would name: LISP1, the UIS/ETRAMP early transcribed membrane
    proteins, and perforin-like protein 1 all come out reduced across this transition, which is the
    textbook set for liver-stage development and hepatocyte egress.
    """
    folder, pmid, name, sheet = PB_LIVER
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    if raw.shape[1] < max(PB_LIVER_COLUMNS.values()) + 1 or len(raw) < 3:
        return pd.DataFrame()
    if str(raw.iloc[0, PB_LIVER_COLUMNS["pb_liver_log2fc"]]).strip() != PB_LIVER_TRANSITION:
        log("pb liver: the sheet's columns are not where this loader expects them")
        return pd.DataFrame()
    body = raw.iloc[2:].reset_index(drop=True)
    mapping = _berghei_to_falciparum(dataset_root)
    if not mapping:
        log("pb liver: no blood-stage table to map berghei ids onto falciparum ones")
        return pd.DataFrame()
    gene = body.iloc[:, 0].astype(str).str.strip().map(mapping)
    value = pd.to_numeric(body.iloc[:, PB_LIVER_COLUMNS["pb_liver_log2fc"]], errors="coerce")
    power = body.iloc[:, PB_LIVER_COLUMNS["pb_liver_power"]].astype(str).str.strip()
    keep = gene.notna() & value.notna() & (power != PB_NO_POWER)
    out = pd.DataFrame({"gene_id": gene[keep].to_numpy(),
                        "pb_transferred_liver_log2fc": value[keep].to_numpy(),
                        "pb_transferred_liver_reduced": (power[keep] == "reduced").to_numpy()})
    out = out.drop_duplicates("gene_id").reset_index(drop=True)
    log(f"pb liver: {len(out):,} falciparum genes carry a liver-stage phenotype; "
        f"{int(out['pb_transferred_liver_reduced'].sum())} are reduced, "
        f"{int((~keep & gene.notna()).sum())} dropped for having no power to say")
    return out


def _berghei_to_falciparum(dataset_root: str) -> dict:
    """{berghei id: falciparum id}, from the blood-stage screen's own ortholog column."""
    folder, pmid, name, sheet = PB_TRANSFER
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path):
        return {}
    d = pd.read_excel(path, sheet_name=sheet)
    if "P. berghei current ID" not in d.columns or "P. falciparum ID" not in d.columns:
        return {}
    pb = d["P. berghei current ID"].astype(str).str.strip()
    pf = d["P. falciparum ID"].astype(str).str.strip().str.replace(r"\.\d+$", "", regex=True)
    ok = pf.str.match(r"^PF3D7_\w+$", na=False)
    return dict(zip(pb[ok], pf[ok]))


# --------------------------------------------------------------------------- structural similarity
#: The Foldseek pair table, written by `scripts/run_foldseek.py` from the proteome's AlphaFold
#: models. `(folder, file)` under the PlasmoDB reference tree.
STRUCT_TABLE = ("structures", "pf_struct_pairs.tsv")

#: The TM-score two structures must reach to be one edge. 0.7 is the Toxoplasma layer's threshold and
#: is kept rather than tuned: the point of the second arm is that a number means the same thing in
#: both, and a layer built at a different cut-off would not be comparable to anything.
STRUCT_TM = 0.7


def structure_similarity(dataset_root: str, log=print) -> pd.DataFrame:
    """Gene pairs whose AlphaFold models superpose, from an all-against-all TM-align search.

    The layer that reaches what the others cannot: structural similarity needs no orthology, so it
    finds relatives among the lineage-specific proteins where a sequence search returns nothing --
    which in this genome is most of the exported families.

    The search itself is a build step (`scripts/run_foldseek.py`) and its output is archived, so the
    threshold lives HERE. Raising or lowering it is a re-read of a table rather than an hour of
    compute, and the table keeps every pair the search reported.
    """
    folder, name = STRUCT_TABLE
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    if not {"accession_a", "accession_b", "alntmscore"} <= set(d.columns):
        return pd.DataFrame()
    owners = uniprot_index(os.path.join(dataset_root, "reference", "plasmodb", UNIPROT_TABLE))
    if not owners:
        log("structures: no UniProt table, so model accessions cannot be resolved to genes")
        return pd.DataFrame()
    score = pd.to_numeric(d["alntmscore"], errors="coerce")
    keep = d.assign(tm=score)[score >= STRUCT_TM]
    keep = keep.assign(gene_a=keep["accession_a"].map(owners), gene_b=keep["accession_b"].map(owners))
    keep = keep.dropna(subset=["gene_a", "gene_b"])
    keep = keep[keep["gene_a"] != keep["gene_b"]]
    if keep.empty:
        return pd.DataFrame()
    # One row per unordered pair, keeping the better score: an all-against-all search reports A->B
    # and B->A, and TM-align is not symmetric, so the two differ slightly.
    pair = pd.DataFrame({"gene_a": keep[["gene_a", "gene_b"]].min(axis=1),
                         "gene_b": keep[["gene_a", "gene_b"]].max(axis=1),
                         "tm": keep["tm"].to_numpy()})
    out = pair.groupby(["gene_a", "gene_b"], as_index=False)["tm"].max()
    log(f"structures: {len(out):,} pairs at TM >= {STRUCT_TM} over "
        f"{len(set(out.gene_a) | set(out.gene_b)):,} genes")
    return out


# --------------------------------------------------------------------------- literature
#: The abstract corpus, outside the package for the same reason the Toxoplasma one is: it is 21 MB of
#: somebody else's text, it grows, and every number derived from it is a snapshot. Relative to the
#: dataset tree's parent, which is where the skills corpora live on this machine.
LITERATURE_CORPUS = os.path.join(".claude", "skills", "plasmodium-scientist", "corpus",
                                 "pubmed_plasmodium.jsonl")


def gene_index(nodes: pd.DataFrame, log=print):
    """The identity index for reading *falciparum* genes out of free text.

    The same builder the Toxoplasma arm uses, with this organism's accession shapes -- current
    `PF3D7_`, and the three older forms its literature still cites -- and this organism's identity
    table. Sharing the builder is the point: symbol matching in prose is where a literature layer
    goes wrong, and one implementation means one set of rules to get right, not two.
    """
    from . import identity, paths
    return identity.build_index(nodes["gene_id"], paths.cache_file(IDENTITY_TABLE), log=log,
                                accession_rx=identity.PF_ACC_RX, canonical_prefix="PF3D7",
                                symbol_prefix="Pf")


def literature_layer(nodes: pd.DataFrame, base: str, log=print) -> tuple:
    """Scan the abstract corpus for gene mentions; return node columns and co-mention edges.

    Composed entirely of the machinery the other arm already uses -- `identity` for who is named,
    `corpus` for what a document is, `literature` for the counting and the attention correction --
    so the two arms' attention numbers mean the same thing. What is Plasmodium-specific is the index
    and the corpus path.

    Returns `(columns, edges, mentions)`: a frame keyed by gene, the co-mention layer in the graph's
    array shape, and the auditable intermediate every figure can be recomputed from.
    """
    from collections import Counter
    from . import corpus, literature
    path = os.path.join(base, LITERATURE_CORPUS)
    if not os.path.exists(path):
        log(f"literature: no corpus at {path}; run scripts/fetch_pubmed_corpus.py")
        return pd.DataFrame(), {}, pd.DataFrame()
    index = gene_index(nodes, log=log)
    documents = list(corpus.iter_abstracts(path))
    log(f"literature: {len(documents):,} abstracts")
    mentions, co, meta = literature.scan(documents, index, log=log)
    if mentions.empty:
        return pd.DataFrame(), {}, mentions
    counts = literature.publication_counts(mentions, "abstract")
    depth = literature.attention_depth(mentions)
    columns = pd.DataFrame({"gene_id": nodes["gene_id"].to_numpy()})
    columns["n_publications"] = columns["gene_id"].map(counts).fillna(0).astype(int)
    for name in depth.columns:
        mapped = columns["gene_id"].map(depth[name])
        columns[name] = mapped.fillna(0).astype(int) if name.startswith("n_") else mapped.fillna("")
    # A tier this corpus cannot produce is not shipped. `incidental` means a mention in a body or a
    # caption, and this arm loads abstracts only -- so the column would be zero for all 5,720 genes,
    # which is a column that says the same thing about everything and reads as a measured absence.
    empty = [c for c in columns.columns
             if c.startswith("n_papers_") and not columns[c].any()]
    if empty:
        columns = columns.drop(columns=empty)
        log(f"literature: no {', '.join(t.replace('n_papers_', '') for t in empty)} tier from "
            f"abstracts alone, so those columns are not shipped")
    position = {gene: i for i, gene in enumerate(nodes["gene_id"])}
    edges = {}
    pairs = literature.comention_edges(co.get("abstract", Counter()), meta, "abstract")
    if pairs:
        edges["comention"] = (
            np.array([position[a] for a, _, _, _ in pairs], dtype=int),
            np.array([position[b] for _, b, _, _ in pairs], dtype=int),
            np.array([w for _, _, w, _ in pairs], dtype=float),
            np.array([r for _, _, _, r in pairs], dtype=float))
        log(f"comention: {len(pairs):,} edges (raw and attention-corrected)")
    named = int((columns["n_publications"] > 0).sum())
    log(f"literature: {named:,} of {len(columns):,} genes are named in at least one abstract")
    return columns, edges, mentions


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


# --------------------------------------------------------------------------- febrile stress
#: Transcription at 37 C and at the 41 C of a malarial fever, in wild type and in two mutants.
FEBRILE_TABLE = "plasmodb_pf3d7_febrile.tsv"
FEBRILE = (("WT 37C", "febrile_wt_37c"), ("WT 41C", "febrile_wt_41c"),
           ("delta-LRR5-37C", "febrile_lrr5ko_37c"), ("delta-LRR5-41C", "febrile_lrr5ko_41c"),
           ("delta-DHC-37C", "febrile_dhcko_37c"), ("delta-DHC-41C", "febrile_dhcko_41c"))


def febrile(report_path: str) -> pd.DataFrame:
    """The six conditions, as conditions. The 41-versus-37 contrast is left to the caller.

    Same restraint as `sir2_perturbation`, for a different reason. There the check on the contrast
    contradicted itself; here it came out NULL: heat shock proteins move by a median log2 of +0.08
    against -0.07 for everything else (p = 0.2), so a fever does not induce them measurably. That is
    consistent with what is known -- Plasmodium's chaperones are constitutively high rather than
    stress-induced -- and the genes that do rise, Maurer's cleft two-TM proteins and stevor at four
    to five log2, match published fever-driven surface remodelling. But a null result on the one
    prediction available is not a validation, and a derived column would imply it had passed one.
    The conditions themselves are unambiguous and are what ship.
    """
    if not os.path.exists(report_path):
        return pd.DataFrame()
    d = pd.read_csv(report_path, sep="\t", dtype=str).replace(dict.fromkeys(BLANK, None))
    if "Gene ID" not in d.columns:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": d["Gene ID"].astype(str)})
    for label, column in FEBRILE:
        found = [c for c in d.columns if label in c]
        if len(found) == 1:
            out[column] = pd.to_numeric(d[found[0]], errors="coerce")
    if len(out.columns) == 1:
        return pd.DataFrame()
    return out[~out["gene_id"].duplicated()].sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- complex membership
#: Complexes called from the crosslinking data. The clusters mix organisms -- seven of the 47 contain
#: both parasite and host proteins -- because the experiment crosslinked parasite inside erythrocyte.
COMPLEXES = ("complexes", "41966402", "mmc5.xlsx")
COMPLEXES_SHEET = "Clusters"


def complexes(dataset_root: str, log=print) -> pd.DataFrame:
    """Which crosslink-derived complex a gene belongs to, and how big it is.

    `complex_spans_host` is kept rather than dropped, and it is the informative column. Seven of the
    47 clusters contain human proteins as well as parasite ones, which is not contamination -- the
    experiment crosslinked parasite inside erythrocyte, so a complex reaching into the host is a
    finding. But a parasite gene in one of those has partners this table cannot name, and a reader
    counting `complex_size` without knowing that would over-count its parasite neighbours.
    """
    folder, pmid, name = COMPLEXES
    path = os.path.join(dataset_root, "reference", "plasmodb", folder, pmid, name)
    if not os.path.exists(path):
        return pd.DataFrame()
    book = pd.ExcelFile(path)
    if COMPLEXES_SHEET not in book.sheet_names:
        return pd.DataFrame()
    d = book.parse(COMPLEXES_SHEET)
    d.columns = [str(c).strip() for c in d.columns]
    number = next((c for c in d.columns if c.startswith("complex nr")), None)
    if not {"uniprot", "organism_name"} <= set(d.columns) or number is None:
        return pd.DataFrame()
    owners = uniprot_index(os.path.join(dataset_root, "reference", "plasmodb", UNIPROT_TABLE))
    if not owners:
        return pd.DataFrame()
    mixed = set(d.groupby(number)["organism_name"].nunique().pipe(lambda s: s[s > 1]).index)
    here = d[d["organism_name"].astype(str).str.contains("falciparum", na=False)].copy()
    here["gene_id"] = here["uniprot"].astype(str).str.strip().map(owners)
    here = here.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    if here.empty:
        return pd.DataFrame()
    out = pd.DataFrame({"gene_id": here["gene_id"].to_numpy(),
                        "complex_id": here[number].to_numpy(),
                        "complex_size": pd.to_numeric(here.get("csize"), errors="coerce").to_numpy(),
                        "complex_spans_host": here[number].isin(mixed).to_numpy()})
    log(f"complexes: {len(out)} genes in {out['complex_id'].nunique()} complexes, "
        f"{int(out['complex_spans_host'].sum())} of them in a complex that reaches the host")
    return out.sort_values("gene_id").reset_index(drop=True)


# --------------------------------------------------------------------------- enzyme classification
#: PlasmoDB serves two EC fields and they are different KINDS of evidence: one curated for this
#: organism, one inferred from its OrthoMCL group. They are kept in separate columns for that reason --
#: merging them would put inference where annotation is, and the merged column would be 1,584 genes
#: with no way to tell which 343 of them were never annotated here at all.
EC_TABLE = "plasmodb_pf3d7_ec.tsv"


def enzyme_classification(dataset_root: str, log=print) -> pd.DataFrame:
    """EC number per gene, curated and orthology-derived kept apart."""
    path = os.path.join(dataset_root, "reference", "plasmodb", EC_TABLE)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t", dtype=str).replace(dict.fromkeys(BLANK, None))
    if "Gene ID" not in d.columns or len(d.columns) < 2:
        return pd.DataFrame()
    d.columns = ["gene_id", "ec_number", "ec_number_orthology"][:len(d.columns)]
    out = d[~d["gene_id"].duplicated()].copy()
    out["has_ec"] = out["ec_number"].notna()
    curated = int(out["has_ec"].sum())
    derived = int(out["ec_number_orthology"].notna().sum()) if "ec_number_orthology" in out else 0
    log(f"enzyme classification: {curated} genes with a curated EC, {derived} with one derived from "
        f"orthology")
    return out.reset_index(drop=True)
