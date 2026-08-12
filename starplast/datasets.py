#!/usr/bin/env python3
"""The dataset registry — one machine-readable record of what starplast is built from.

This is the single source of truth for provenance. The download notebook fetches from it, the methods
documentation is generated from it, and the application shows it when asked where a number came from.
Previously that information lived in three prose documents that could drift apart; here a dataset is
described once.

    from starplast import datasets
    datasets.registry()                 # every entry
    datasets.registry(level="DNA")      # filter
    datasets.provenance("fit_ifng")     # which dataset produced a node column
    datasets.unresolved()               # entries whose citation still needs confirming

`level` follows the on-disk layout (`datasets/<level>/<type>/<PMID>/`):

    DNA              genetic perturbation or DNA-level readout   CRISPR screens, ChIP-seq
    transcription    RNA abundance                               RNAseq, scRNAseq
    translation      protein abundance                           LFQ / iBAQ / iTRAQ, Ribo-seq
    post_translation properties of the folded protein            IPMS, BioID, XLMS, phospho, LOPIT
    reference        not a study result                          orthology, domains, identity tables

`citation` is None where the originating publication is not recorded in the pipeline. That is a real gap,
not an oversight to paper over: the data is verified on disk, the citation is not, and `unresolved()`
lists exactly those so they cannot reach a manuscript unchecked.
"""
from __future__ import annotations

import json
import os

from dataclasses import dataclass, field, asdict


@dataclass(frozen=True)
class Dataset:
    """One dataset the map is built from, and everything known about where it came from.

    The unit of provenance for the whole project: the README table, the methods document and the
    per-dataset scripts are all generated from these, so a source is described once and cannot drift
    between three prose documents. `citation` is None where the originating publication has not been
    confirmed -- a real gap, listed by `unresolved()`, rather than a guess that could reach a paper.
    """
    key: str
    name: str
    level: str
    kind: str
    provides: str
    columns: tuple = ()          # node-table columns this produces
    coverage: str = ""           # genes covered, as measured at build time
    pmid: str | None = None
    accession: str | None = None
    citation: str | None = None  # None = confirm before citing
    url: str | None = None       # direct download, where one exists
    path: str | None = None      # where it lands under toxoplasma_projects/
    derived_from: tuple = ()     # node columns this was COMPUTED from, if it is a derivation
    note: str = ""


# Springer and PLOS serve supplementary directly; PMC's /bin/ path 404s.
SPRINGER = "https://static-content.springer.com/esm/art%3A10.1038%2F{doi}/MediaObjects/{f}"
PLOS = "https://journals.plos.org/plospathogens/article/file?id=10.1371/{doi}.{s}&type=supplementary"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
TOXODB = ("https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
          "/reports/attributesTabular")

REGISTRY = [
    # ------------------------------------------------------------------ transcription
    Dataset("xue_singlecell", "Single-parasite transcriptional atlas (cell cycle)", "transcription",
            "scRNAseq",
            "Measured cell-cycle phase per gene, and pseudotime cluster",
            ("cellcycle_phase", "cellcycle_pseudotime"),
            "873 genes phased, 7,499 clustered", pmid="32065584",
            citation="Xue Y et al. eLife 2020;9:e54129",
            url="https://cdn.elifesciences.org/articles/54129/elife-54129-supp3-v2.csv",
            path="datasets/transcription/scRNAseq/32065584/cellcycle_phase_RH.csv",
            note="Tab-separated despite the .csv extension. RH files use TGGT1_ accessions and the "
                 "Pru files in the same supplement use TGME49_; the prefix is the only thing that "
                 "distinguishes them. The only MEASURED discrete cell-cycle label in the project."),
    # kind names the ASSAY, not the provenance. This was "derived", which answers how the column was
    # produced rather than what kind of measurement is underneath it -- and left the one entry in the
    # registry whose type you could not read off its type field. It is bulk RNA-seq: the argmax of
    # three expression columns that come from GSE108740 and GSE206344. That it is a derivation is
    # said by the name, by the note, and machine-readably by derived_from being non-empty.
    Dataset("stage_enriched", "Life-cycle stage enrichment (DERIVED)", "transcription", "RNAseq",
            "Which stage a gene's own expression is highest in",
            ("stage_enriched_derived", "stage_margin_derived"),
            "1,911 of 8,140 genes called",
            derived_from=("expr_tachy", "expr_cyst", "expr_sporulated"),
            note="DERIVED, not measured: computed here from expr_tachy / expr_cyst / expr_sporulated "
                 "by z-scoring each and taking the argmax where it leads by 0.5 z. It is a "
                 "restatement of those columns, so holding it out against an embedding built on them "
                 "is circular by construction. Left unlabelled where no stage leads clearly."),

    # ------------------------------------------------------------------ reference
    Dataset("toxodb_identity", "ToxoDB gene identity", "reference", "identity",
            "Symbols, previous IDs, product descriptions", ("gene_id", "product"),
            "8,843 ME49 genes", accession="ToxoDB ME49", url=TOXODB,
            path="starplast/data/toxodb_identity.tsv",
            note="Retrieved 2026-08-11 via the REST API; strain tables for GT1 and VEG alongside."),
    Dataset("orthomcl", "OrthoMCL orthogroups", "reference", "orthology",
            "Orthogroup assignment and cross-species bridge", ("orthogroup",),
            "16,793 groups", accession="OrthoMCL release 6.21",
            url="https://orthomcl.org/common/downloads/release-6.21/groups_OrthoMCL-6.21.txt.gz",
            path="datasets/MASTER_parasite_wide_by_orthogroup.csv",
            note="The RELEASE is load-bearing, not decoration: all 16,793 groups reproduce exactly "
                 "from 6.21, while 6.20 differs in 180 cells and Current_Release (v7) renumbers "
                 "every group into an OG7_ namespace that matches nothing here. The shipped CSV is "
                 "derived from this file -- filter to tgon/pfal/cpar/tbrt, then pivot wide. Note the "
                 "T. brucei taxon code is tbrt, though the CSV column is tbru."),
    Dataset("interpro", "InterPro domains", "reference", "domains",
            "Domain identity and count", ("n_interpro", "interpro_id", "interpro_desc", "pfam_id"),
            "8,140",
            url=TOXODB + ("?organism=%5B%22Toxoplasma%20gondii%20ME49%22%5D&reportConfig=%7B%22"
                          "attributes%22%3A%5B%22gene_source_id%22%2C%22interpro_id%22%2C%22"
                          "interpro_description%22%2C%22pfam_id%22%2C%22pfam_description%22%5D%2C%22"
                          "includeHeader%22%3Atrue%2C%22attachmentType%22%3A%22plain%22%7D"),
            path="datasets/interpro_tgon.csv",
            note="Domain annotation partly records study effort, not conserved architecture."),
    Dataset("alphafold", "AlphaFold DB", "reference", "structure",
            "Per-gene mean pLDDT; coordinates fetched on demand", ("mean_plddt",),
            "6,480 (79.6%)", accession="UP000001529 (taxid 508771), AlphaFold DB",
            citation="Varadi et al. 2024 NAR (database); Jumper et al. 2021 Nature (method)",
            url="https://alphafold.ebi.ac.uk/api/prediction/{acc}",
            note="Missing for the largest proteins, which here are disproportionately secreted "
                 "effectors. NO ANONYMOUS BULK DOWNLOAD: the per-proteome tar exists at "
                 "gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-508771-0_v4"
                 ".tar but plain HTTPS returns 403, so it needs `gcloud storage cp` and a Google "
                 "account. The per-accession API above is the credential-free route and is what "
                 "structures.py uses. The database paper is Varadi et al., not Jumper et al. -- "
                 "Jumper is the method, and there is no paper by Jumper titled after the database."),

    # ------------------------------------------------------------------ localisation
    Dataset("lopit_tgon", "T. gondii hyperLOPIT", "post_translation", "LOPIT",
            "Subcellular compartment, MAP and MCMC, with posteriors",
            ("compartment", "lopit_map", "lopit_mcmc", "lopit_prob_map", "lopit_prob_mcmc"),
            "3,827 (47.0%)", pmid="33053376",
            citation="A Comprehensive Subcellular Atlas of the Toxoplasma Proteome via hyperLOPIT "
                     "(Barylyuk et al. 2020)",
            url="https://ars.els-cdn.com/content/image/1-s2.0-S193131282030514X-mmc5.xls",
            path="datasets/lopit_toxoplasma_gondii_ME49.csv",
            note="MAP and MCMC disagree for 980 of 3,827 (26%). Assignment tracks abundance, so the "
                 "unassigned half is biased toward low-abundance proteins."),
    Dataset("lopit_pfal", "P. falciparum LOPIT", "post_translation", "LOPIT",
            "Donor labels for orthoLOPIT transfer", (), "1,646 usable", pmid="42218142",
            citation="Chisholm SA et al., The spatial proteome of the Plasmodium falciparum "
                     "schizont. Nat Commun 2026;17:6192 -- CONFIRM against the file on disk",
            url="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2F"
                "s41467-026-73664-2/MediaObjects/41467_2026_73664_MOESM3_ESM.xlsx",
            path="datasets/lopit_plasmodium_falciparum_3D7.csv",
            note="The URL downloads and is the right kind of data, but this file predates the "
                 "registry entry, so that this paper is the source of THIS csv is inference, not "
                 "verification. Confirm against the file before citing."),
    Dataset("lopit_cpar", "C. parvum hyperLOPIT", "post_translation", "LOPIT",
            "Donor labels for orthoLOPIT transfer", (), "1,107 usable",
            citation="Guerin et al. 2023",
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1931312823001051-mmc4.xlsx",
            path="datasets/lopit_cryptosporidium_parvum_MEASURED_Guerin2023.csv",
            note="MASTER_parasite_wide_by_orthogroup.csv has an empty cpar_lopit_native column; this "
                 "data is joined from source instead."),

    # ------------------------------------------------------------------ transcription
    Dataset("gse108740", "Stage transcriptome", "transcription", "RNAseq",
            "Tachyzoite, day 3/5/7, in vivo tissue cyst (12 columns)",
            ("expr_tachy", "expr_cyst", "expr_max"), "7,739 (95.1%)", accession="GSE108740",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108740"),
    Dataset("gse206344", "Oocyst sporulation series", "transcription", "RNAseq",
            "Unsporulated / sporulating / sporulated, 2 replicates (6 columns)",
            ("expr_sporulated",), "7,974 (98.0%)", accession="GSE206344",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206344"),

    # ------------------------------------------------------------------ genetic screens
    Dataset("crispr_invitro", "In vitro CRISPR fitness (HFF)", "DNA", "CRISPR_screen",
            "Competitive growth in fibroblasts", ("fit_invitro_hff",), "7,325 (90.0%)",
            pmid="27594426",
            citation="A Genome-wide CRISPR Screen in Toxoplasma Identifies Essential Apicomplexan Genes "
                     "(Sidik et al. 2016)",
            url="https://ars.els-cdn.com/content/image/1-s2.0-S0092867416310704-mmc3.xlsx",
            note="Competitive growth, NOT essentiality. Predicted from protein features at R2 = 0.453, "
                 "while the other screens are predicted at -0.105 to +0.102."),
    Dataset("crispr_invivo_composite", "In vivo CRISPR composite scores", "DNA", "CRISPR_screen",
            "Peritoneum, lung, liver, spleen composite scores",
            ("fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver", "fit_invivo_spleen"),
            "7,395 (90.8%)", pmid="31481656",
            accession="ToxoDB tgonGt1CrisprFunc*",
            url=TOXODB + ("?organism=%5B%22Toxoplasma%20gondii%20GT1%22%5D&reportConfig=%7B%22"
                          "attributes%22%3A%5B%22primary_key%22%2C%22tgonGt1CrisprMeanPhenotype%22"
                          "%2C%22tgonGt1CrisprFuncPE%22%2C%22tgonGt1CrisprFuncLung%22%2C%22"
                          "tgonGt1CrisprFuncLiver%22%2C%22tgonGt1CrisprFuncSpleen%22%5D%2C%22"
                          "includeHeader%22%3Atrue%2C%22attachmentType%22%3A%22plain%22%7D"),
            note="Corresponds to the in vivo CRISPR platform paper; confirm before citing. The URL "
                 "pulls the tgonGt1CrisprFunc* tracks straight from ToxoDB, keyed on GT1."),
    Dataset("crispr_macrophage", "Macrophage CRISPR screens", "DNA", "CRISPR_screen",
            "Naive BMDM and IFN-gamma survival", ("fit_naive_bmdm", "fit_ifng"), "7,402 (90.9%)",
            pmid="33067458",
            citation="Wang Y et al., Genome-wide screens identify Toxoplasma gondii determinants of "
                     "parasite fitness in IFN-gamma-activated murine macrophages. Nat Commun "
                     "2020;11:5258",
            url=SPRINGER.format(doi="s41467-020-18991-8", f="41467_2020_18991_MOESM5_ESM.xlsx"),
            note="Sign convention is INVERTED relative to the other screens. THE PREVIOUS CITATION "
                 "WAS WRONG and its own 'CONFIRM this is the source' warning was justified: PMID "
                 "25867017 is a 2015 JoVE video protocol using CHEMICAL mutagenesis, verified "
                 "against PubMed -- not a CRISPR screen, and it predates the first one. The entry "
                 "had copied that protocol's title verbatim. A video protocol cannot be the source "
                 "of 7,402 per-gene fitness scores; Wang 2020 is genome-wide in IFN-gamma-activated "
                 "macrophages, which is exactly what these two columns are."),
    Dataset("crispr_young2019", "Young 2019 in vivo screen", "DNA", "CRISPR_screen",
            "In vivo fitness", ("fit_invivo_young2019",), "115", pmid="31481656",
            citation="Young J et al., A CRISPR platform for targeted in vivo screens identifies "
                     "Toxoplasma gondii virulence factors in mice. Nat Commun 2019;10:3963",
            url=SPRINGER.format(doi="s41467-019-11855-w", f="41467_2019_11855_MOESM6_ESM.xlsx"),
            note="Same paper as invivo_platform, a different supplementary table (MOESM6 vs MOESM5). "
                 "TARGETED, not genome-wide: the libraries are 200, 800 and 3200 gRNAs, which is why "
                 "this covers 115 genes rather than the proteome."),
    Dataset("gra17_synthlethal", "GRA17 synthetic-lethal screen", "DNA", "CRISPR_screen",
            "RH and RH-delta-gra17 phenotype by passage; MAGeCK p-values",
            ("crispr_gra17ko_phenotype", "crispr_gra17_synthlethal_delta", "crispr_gra17_candidate"),
            "7,553 (genome-wide)", pmid="37498952",
            citation="Genome-wide CRISPR screen identifies genes synthetically lethal with GRA17, "
                     "a nutrient channel encoding gene in Toxoplasma",
            url=PLOS.format(doi="journal.ppat.1011543", s="s001"),
            path="datasets/DNA/CRISPR_screen/37498952/"),
    Dataset("invivo_platform", "In vivo CRISPR platform", "DNA", "CRISPR_screen",
            "Mean log fold-change across replicates", ("crispr_invivo_platform_lfc",), "168",
            pmid="31481656",
            citation="A CRISPR platform for targeted in vivo screens identifies Toxoplasma gondii "
                     "virulence factors in mice",
            url=SPRINGER.format(doi="s41467-019-11855-w", f="41467_2019_11855_MOESM5_ESM.xlsx"),
            path="datasets/DNA/CRISPR_screen/31481656/",
            note="Cites pre-2012 TGME49_0xxxxx accessions for every gene; must go through identity.py."),
    Dataset("gra12", "GRA12 strains and mouse subspecies", "DNA", "CRISPR_screen",
            "Median L2FC in vitro and in vivo, DISCO score; two screens",
            ("crispr_gra12s1_l2fc_invivo", "crispr_gra12s2_l2fc_invivo"), "236 / 232",
            pmid="40240328",
            citation="GRA12 is a common virulence factor across Toxoplasma gondii strains and "
                     "mouse subspecies",
            url=SPRINGER.format(doi="s41467-025-58876-2", f="41467_2025_58876_MOESM5_ESM.xlsx"),
            path="datasets/DNA/CRISPR_screen/40240328/",
            note="The two screens are NOT replicates: in-vivo L2FC correlate at r = 0.41."),
    Dataset("hosttx_effectors", "Host-transcription effector screen", "DNA", "CRISPR_screen",
            "Hotelling T2 statistic per effector, adjusted p", ("hosttx_T2", "hosttx_padj"), "252",
            pmid="37827122",
            citation="High-throughput identification of Toxoplasma gondii effector proteins that "
                     "target host cell transcription",
            url=EPMC.format(pmcid="PMC12033024"),
            path="datasets/DNA/CRISPR_screen/37827122/"),

    # ------------------------------------------------------------------ protein level
    Dataset("proteome_pru", "Pru proteome and IP abundance", "translation", "proteomics",
            "Median log2 iBAQ across replicates", ("protein_ibaq_log2",), "748 (9.2%)",
            accession="PXD043808, PXD065585",
            url="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD065585",
            path="toxo_stage_atlas/data/proteomics/",
            note="Immunoprecipitation experiments of 424 and 594 proteins. Enrichment, NOT a deep "
                 "proteome; do not report as proteome-wide."),
    Dataset("phosphosites", "Phosphosite counts", "post_translation", "phosphoproteomics",
            "Count of phosphosites per protein, no positions", ("n_phosphosites",), "1,175 (14.4%)",
            citation="Treeck M et al. 2011 -- CONFIRM against the file on disk",
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1931312811002885-mmc2.xls",
            note="The URL downloads a real phosphoproteomics table, but that it is the source of "
                 "THIS column is inference rather than verification; confirm before citing. "
                 "Missing for 85.6% of genes; effectively an indicator of having been in a "
                 "phosphoproteomics experiment."),

    # ------------------------------------------------------------------ interactions
    Dataset("starpath_xlms", "StarPath crosslink MS", "post_translation", "XLMS",
            "Measured physical proximity; residue-level crosslinks and Chai-1 complexes",
            ("n_xlink_partners", "best_model_agreement"), "2,842 pairs / 1,630 genes", pmid="40874616",
            citation="Mapping a Toxoplasma gondii interactome by crosslinking mass spectrometry and "
                     "machine learning (2025)",
            url=EPMC.format(pmcid="PMC12505969"),
            path="starpath_crosslinks.json, starpath_dump/cifs/",
            note="RH88 accessions do NOT map to ME49 by suffix; use the alias column. 60% of predicted "
                 "complexes place no crosslink within reach; only 162 pairs are trustworthy."),
    Dataset("ipms_baits", "IP-MS of tagged baits", "post_translation", "IPMS",
            "Replicated pulldown vs untagged control", ("n_ipms_partners",), "64 pairs / 48 genes",
            accession="PXD043808, PXD065585",
            url="https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD043808/files",
            note="The PRIDE API lists the files for one accession; swap the accession in the path "
                 "for PXD065585. It returns JSON metadata, not the data -- follow the download "
                 "links it gives."),
    Dataset("foldseek_struct", "Foldseek structural similarity", "post_translation", "structure",
            "TM-align over Toxoplasma AlphaFold models, TM >= 0.7", ("n_struct_similar",),
            "11,684 pairs / 2,338 genes",
            derived_from=("mean_plddt",),
            note="COMPUTED HERE, so there is nothing to download: Foldseek all-vs-all over the "
                 "Toxoplasma AlphaFold models (see the alphafold entry for how to obtain those), "
                 "keeping pairs at TM >= 0.7. Needs no orthology, so it reaches lineage-specific effectors homology edges cannot."),
    Dataset("bioid_corpus", "Proximity-labelling corpus", "post_translation", "BioID",
            "42 BioID/TurboID/APEX studies with a tagged Toxoplasma protein", (),
            "28 studies with data, 127 files",
            url=EPMC.format(pmcid="{pmcid}"),
            path="datasets/post_translation/BioID/",
            note="Downloaded and indexed; NOT yet parsed into edges. HOW TO REFETCH: this is a "
                 "harvest of many papers, so there is no single URL -- but nothing is lost. Every "
                 "one of the 97 studies has its PMCID and every file its filename recorded in "
                 "interaction_studies.parquet and interaction_study_members.parquet, so the whole "
                 "corpus reconstructs by substituting each PMCID into the URL above."),
    Dataset("pulldown_corpus", "Pulldown corpus", "post_translation", "IPMS",
            "55 IP-MS / co-IP studies with a tagged Toxoplasma protein", (),
            "29 studies with data, 140 files",
            url=EPMC.format(pmcid="{pmcid}"),
            path="datasets/post_translation/IPMS/",
            note="Downloaded and indexed; NOT yet parsed into edges. Refetch as for bioid_corpus: "
                 "substitute each recorded PMCID into the URL above."),

    # ------------------------------------------------------------------ literature
    Dataset("pubmed", "PubMed abstracts", "reference", "literature",
            "Titles and abstracts for co-mention and attention",
            ("n_publications", "attention_depth"), "33,924 records",
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                "?db=pubmed&term=Toxoplasma&retmax=100000",
            path=".claude/skills/toxoplasma-scientist/corpus/pubmed_toxoplasma.jsonl",
            note="ASSEMBLED HERE from an E-utilities query rather than downloaded as a file: the "
                 "esearch above returns the PMID set, and efetch retrieves each record. The count "
                 "grows over time, so a rebuild will not reproduce 33,924 exactly -- record the "
                 "date. This corpus was built 2026."),
    Dataset("pmc_oa", "PubMed Central open-access full texts", "reference", "literature",
            "Sectioned JATS XML", ("n_fulltext",), "6,667 articles",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
            path="/mnt/wd4tb/skill_corpora/toxoplasma-scientist/",
            note="A biased subset: only what publishers deposited open access. ASSEMBLED HERE: for "
                 "each PMID in the pubmed corpus that has a PMCID, fetch the JATS from the URL "
                 "above. Machine-local by size, which is why the built cache is what ships."),
]

_BY_KEY = {d.key: d for d in REGISTRY}
_BY_COLUMN = {c: d for d in REGISTRY for c in d.columns}


def registry(level: str | None = None, kind: str | None = None) -> list:
    """Every dataset, optionally filtered by level or kind."""
    out = REGISTRY
    if level:
        out = [d for d in out if d.level == level]
    if kind:
        out = [d for d in out if d.kind == kind]
    return list(out)


def get(key: str) -> Dataset:
    """The dataset with this key. Raises KeyError naming the key if there is none."""
    return _BY_KEY[key]


def provenance(column: str) -> Dataset | None:
    """Which dataset produced a given node-table column."""
    return _BY_COLUMN.get(column)


def derived_sources(column: str) -> tuple:
    """The node columns a derived column was computed from, as DECLARED by whoever derived it.

    Declaration is needed in addition to measurement, not instead of it. Measured association is
    pairwise, and a label computed as the argmax of three columns is a JOINT function of them: each
    source individually explains only about 0.6 of it, under any sensible exclusion threshold, so a
    pairwise measure cannot see the derivation however good the statistic is. Measurement catches the
    undeclared leak -- a renamed copy -- and declaration catches the joint one. Neither alone is enough,
    and this project has already been burnt by trusting one of them.
    """
    for d in REGISTRY:
        if column in d.columns and d.derived_from:
            return tuple(d.derived_from)
    return ()


def unresolved() -> list:
    """Entries whose originating publication is not recorded. Confirm before citing."""
    return [d for d in REGISTRY if d.citation is None]


def downloadable() -> list:
    """Entries with a direct download URL, in the order the notebook should fetch them."""
    return [d for d in REGISTRY if d.url]


def as_table() -> "object":
    """The registry as a DataFrame, for the methods table and the application's provenance view."""
    import pandas as pd
    return pd.DataFrame([asdict(d) for d in REGISTRY])


# --------------------------------------------------------------------------- locating and fetching
def local_path(key: str) -> str | None:
    """Where this dataset's file is on THIS machine, or None if it is not here.

    The registry's `path` is recorded relative to the dataset root, which differs between a clone, the
    historical layout and a download cache -- so it is resolved through paths.find rather than opened
    directly. Opening it directly is what made the build work on exactly one machine.
    """
    from . import paths
    d = get(key)
    if not d.path:
        return None
    rel = d.path[len("datasets/"):] if d.path.startswith("datasets/") else d.path
    if d.path.startswith("starplast/data/"):        # part of the committed cache, not a raw input
        p = paths.cache_file(d.path[len("starplast/data/"):])
        return p if __import__("os").path.exists(p) else None
    return paths.find(rel)


# A URL that names a file can be fetched. A URL that names a landing page cannot -- GEO's acc.cgi and
# ProteomeXchange's GetDataset return HTML describing the data, and saving that HTML as though it were
# the dataset is the kind of failure that stays invisible until something tries to parse it.
_PAGE_HOSTS = ("ncbi.nlm.nih.gov/geo/query", "proteomecentral", "toxodb.org/toxo/service")
_FILE_SUFFIXES = (".xlsx", ".xls", ".csv", ".tsv", ".txt", ".zip", ".gz", ".docx", ".pdf")


def fetchable(key: str) -> tuple:
    """(can_fetch, how). Says which of the three situations a dataset is in, rather than guessing."""
    d = get(key)
    if not d.url:
        return False, "no download URL recorded"
    if d.url.startswith(TOXODB):
        # ToxoDB's tabular report is a POST with a JSON body, so the URL alone looks like a landing
        # page while the data is entirely fetchable -- fetch_names has done it all along. Reporting it
        # as unfetchable understated what a clean machine can rebuild by one dataset.
        return True, "toxodb"
    if d.accession and d.accession.startswith("GSE"):
        return True, "geo"                       # the FTP supplementary listing, not the landing page
    if any(h in d.url for h in _PAGE_HOSTS):
        return False, "URL is a landing page, not a file"
    if d.url.lower().endswith(_FILE_SUFFIXES) or "MediaObjects" in d.url or "type=supplementary" in d.url:
        return True, "direct"
    return False, "URL is not recognisably a file"


def ensure(key: str, log=print) -> str | None:
    """Return a local path for this dataset, downloading it if it is absent and can be fetched.

    Returns None rather than raising when the data genuinely cannot be obtained automatically. Several
    of these datasets exist only inside a paper's supplementary section and one exists only as raw
    instrument files; pretending otherwise would produce a path to something that is not the dataset.
    """
    from . import paths, sources

    p = local_path(key)
    if p:
        return p

    d = get(key)
    ok, how = fetchable(key)
    if not ok:
        log(f"{key}: cannot fetch automatically -- {how}")
        return None

    dest = os.path.join(paths.dataset_root(create=True), "_downloads", key)
    os.makedirs(dest, exist_ok=True)
    if how == "geo":
        got = sources.geo_supplementary(d.accession, dest, log=log)
        return got[0] if got else None

    if how == "toxodb":
        from . import fetch_names
        out = os.path.join(dest, os.path.basename(d.path or f"{key}.tsv"))
        try:
            fetch_names.write(fetch_names.fetch("Toxoplasma gondii ME49",
                                                ["primary_key", "gene_name", "gene_previous_ids",
                                                 "gene_product"]), out)
        except Exception as e:                    # noqa: BLE001 -- any transport failure is the same
            log(f"{key}: ToxoDB request failed ({e})")
            return None
        record_checksum(key, out)
        return out

    name = d.url.rsplit("/", 1)[-1].split("?")[0] or f"{key}.dat"
    out = os.path.join(dest, name)
    if not os.path.exists(out):
        try:
            data = sources._get(d.url)
        except Exception as e:                    # noqa: BLE001 -- any transport failure is the same answer
            log(f"{key}: download failed ({e})")
            return None
        # An HTML error page is a successful HTTP response, so size alone does not prove a file arrived.
        if data[:15].lstrip().lower().startswith(b"<!doctype html") or data[:6].lower() == b"<html>":
            log(f"{key}: server returned a web page, not a file")
            return None
        with open(out, "wb") as fh:
            fh.write(data)
        # Pinned only now, after the content has been accepted as a file rather than as an HTML error
        # page -- otherwise the first bad download becomes the truth every good one is measured against.
        record_checksum(key, out)
        log(f"{key}: fetched {len(data)/1e6:.1f} MB -> {out}")
    else:
        verify_checksum(key, out, log=log)
    return out


def missing() -> list:
    """Registry entries whose data is not on this machine. The honest first-run report."""
    return [d.key for d in REGISTRY if d.path and not local_path(d.key)]


# --------------------------------------------------------------------------- checksums
# Recorded on first fetch and checked on every one after. A publisher reissuing a supplement under the
# same URL is the failure this exists for: the file changes, the build re-runs, every number moves a
# little, and nothing anywhere says why.
#
# The checksum is recorded AFTER the content has been accepted as a file rather than as an HTML error
# page, because otherwise the first bad download becomes the pinned truth and every good one after it
# is reported as the corruption.
CHECKSUMS = "checksums.json"


def _checksum_path(paths_mod) -> str:
    return os.path.join(paths_mod.dataset_root(create=True), CHECKSUMS)


def digest(path: str) -> str:
    """SHA-256 of a file, streamed -- some of these are gigabytes."""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def recorded_checksums() -> dict:
    """The SHA-256 recorded for each downloaded file, so a reissued supplement is caught."""
    from . import paths
    p = _checksum_path(paths)
    if not os.path.exists(p):
        return {}
    try:
        with open(p) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return {}


def record_checksum(key: str, path: str) -> str:
    """Pin what was fetched. Returns the digest."""
    from . import paths
    d = digest(path)
    store = recorded_checksums()
    store[key] = {"sha256": d, "bytes": os.path.getsize(path),
                  "file": os.path.basename(path)}
    with open(_checksum_path(paths), "w") as fh:
        json.dump(store, fh, indent=1, sort_keys=True)
    return d


def verify_checksum(key: str, path: str, log=print) -> bool:
    """True when the file matches what was pinned, or when nothing was pinned yet.

    An unpinned file is not a failure -- most of this tree arrived before checksums existed -- so the
    honest answer for it is "no claim", and the digest is recorded so the NEXT fetch has something to
    check against.
    """
    store = recorded_checksums()
    if key not in store:
        record_checksum(key, path)
        return True
    want = store[key]["sha256"]
    got = digest(path)
    if got == want:
        return True
    log(f"{key}: CHECKSUM MISMATCH -- the file on disk is not the one that was pinned.\n"
        f"  pinned {want[:16]}...  now {got[:16]}...\n"
        f"  A publisher reissuing a supplement under the same URL looks exactly like this. Delete "
        f"{path} to re-fetch, or re-pin deliberately with record_checksum().")
    return False


# --------------------------------------------------------------------------- documentation
LEVEL_ORDER = ("DNA", "transcription", "translation", "post_translation", "reference")
LEVEL_TITLE = {
    "DNA": "DNA — genetic perturbation and DNA-level readouts",
    "transcription": "Transcription — RNA abundance",
    "translation": "Translation — protein abundance",
    "post_translation": "Post-translation — properties of the folded protein",
    "reference": "Reference — not a study result",
}


def _reference(d: "Dataset") -> str:
    """What a reader needs to judge the data: the publication, not the file path.

    A path is an implementation detail that changes with the layout; a PMID does not. Where neither a
    citation nor an accession is recorded, that is stated rather than left blank -- `unresolved()`
    exists so those cannot reach a manuscript unchecked, and a blank cell would hide them.
    """
    bits = []
    if d.citation:
        bits.append(d.citation)
    if d.pmid:
        bits.append(f"PMID [{d.pmid}](https://pubmed.ncbi.nlm.nih.gov/{d.pmid}/)")
    if d.accession:
        bits.append(f"`{d.accession}`")
    return "; ".join(bits) if bits else "*citation not yet confirmed*"


def readme_table() -> str:
    """The dataset table for the README, as markdown, generated from this registry.

    Generated rather than written by hand because a hand-written table drifts the moment a dataset is
    added, and a README that misstates which data is inside is worse than one that omits it. A test
    asserts the committed README matches what this produces.
    """
    lines = []
    for level in LEVEL_ORDER:
        entries = [d for d in REGISTRY if d.level == level]
        if not entries:
            continue
        lines.append(f"### {LEVEL_TITLE[level]}")
        lines.append("")
        lines.append("| Dataset | Type of data | Coverage | Reference |")
        lines.append("|---|---|---|---|")
        for d in sorted(entries, key=lambda x: x.name):
            lines.append(f"| {d.name} | {d.provides} | {d.coverage or '—'} | {_reference(d)} |")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
