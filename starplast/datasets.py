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

# Exact shipped columns for multi-column assays.  These are deliberately explicit rather than
# prefixes: provenance is a statement about which measurements a dataset produced, and a future
# column that happens to share a prefix must not silently inherit the wrong paper.  The same tuples
# are used by the slot/embedding audit, so adding a column to the cache without assigning its source
# becomes visible immediately.
GSE108740_COLUMNS = (
    "expr_tachy", "expr_cyst", "expr_max",
    "rna108740_Tachyzoites_T2_FPKM", "rna108740_Tachyzoites_T4_FPKM",
    "rna108740_Tissue_cysts_A_FPKM", "rna108740_Tissue_cysts_B_FPKM",
    "rna108740_Day3_CS4_FPKM", "rna108740_Day3_CS5_FPKM", "rna108740_Day3_CS6_FPKM",
    "rna108740_Day5_CS2_FPKM", "rna108740_Day5_K4_14_FPKM",
    "rna108740_Day7_CS7_FPKM", "rna108740_Day7_CS8_FPKM", "rna108740_Day7_CS9_FPKM",
)
GSE206344_COLUMNS = (
    "expr_sporulated", "rna206344_Unsporulated_R1", "rna206344_Unsporulated_R2",
    "rna206344_Sporulating_R1", "rna206344_Sporulating_R2",
    "rna206344_Sporulated_R1", "rna206344_Sporulated_R2",
)
GSE22258_COLUMNS = ("rna22258_tachyzoite", "rna22258_bradyzoite")
GSE168465_COLUMNS = tuple(
    f"brain168465_{day}_{measure}"
    for day in ("1d", "2d", "4d", "7d", "14d")
    for measure in ("base_mean", "lfc")
)
GSE99395_COLUMNS = tuple(
    f"{assay}99395_{context}_r{rep}"
    for assay in ("rpf", "rna", "te")
    for context in ("extracellular", "intracellular")
    for rep in (1, 2)
)
GSE129869_COLUMNS = tuple(
    f"{assay}129869_{context}_r{rep}"
    for assay in ("rpf", "rna", "te")
    for context in ("confluent", "subconfluent")
    for rep in (1, 2, 3)
)
GSE19092_COLUMNS = tuple(
    f"cellcycle19092_{condition}_r{rep}"
    for condition in ("async", "blocked", *(f"{hour}h" for hour in range(1, 13)))
    for rep in (1, 2)
)
GSE51780_COLUMNS = tuple(
    [f"rna51780_tachy_r{rep}" for rep in (1, 2)]
    + [f"rna51780_mero_r{rep}" for rep in (3, 4, 5, 6)]
)
GSE168155_COLUMNS = tuple(
    f"cpsf4rna168155_{condition}"
    for condition in ("ut_1", "ut_2", "iaa_7h_1", "iaa_7h_2",
                      "iaa_24h_1", "iaa_24h_2", "iaa_48h_1", "iaa_48h_2")
)
GSE200962_COLUMNS = tuple(f"restriction200962_{name}" for name in (
    "7th_gt_2_kk31_kim_s1_l001", "7th_gt_20_kk27_kim_s14_l001",
    "7th_gt_22_kk29_kim_s15_l001", "7th_gt_4_kk32_kim_s2_l001",
    "kk17_p2wo_ph8_2_s5", "kk18_p2wo_ph8_2_s1", "kk19_p2w_ph8_2_s22",
    "kk20_p2w_ph8_2_s16", "kk21_p5wo_ph8_2_s13", "kk22_p5wo_ph8_2_s15",
    "kk23_p5w_ph8_2_s12", "kk24_p5w_ph8_2_s11", "kk5_p2wo_ph7_4_s6",
    "kk6_p2wo_ph7_4_s19", "kk7_p2w_ph7_4_s14", "kk8_p2w_ph7_4_s23",
))
GSE245775_COLUMNS = tuple(
    f"{assay}245775_{context}_r{rep}"
    for assay in ("rpf", "rna", "te")
    for context in ("parent_tachy", "parent_prebrady", "eif12ko_tachy", "eif12ko_prebrady")
    for rep in (1, 2, 3)
)
GSE25388_COLUMNS = ("fit_hyperlopit_unassigned_invivo_lib1",
                    "fit_hyperlopit_unassigned_invivo_lib2")
INVIVO_BRAIN_COLUMNS = (
    "invivo_TZ_1", "invivo_TZ_2",
    "invivo_WholeBrain_Acute_1", "invivo_WholeBrain_Acute_2", "invivo_WholeBrain_Acute_3",
    "invivo_WholeBrain_Chronic_1", "invivo_WholeBrain_Chronic_2",
    "invivo_WholeBrain_Chronic_3", "invivo_BZ_28DPI_1", "invivo_BZ_28DPI_3",
    "invivo_BZ_90DPI_1", "invivo_BZ_90DPI_2", "invivo_BZ_120DPI_2",
    "invivo_BZ_120DPI_3",
)
HOST_SIGNATURE_COLUMNS = tuple(f"hosttx_signature_pc{i:02d}" for i in range(1, 21)) + (
    "hosttx_signature_norm", "hosttx_signature_n_de")
STRESS_COLUMNS = (
    "stress_5-1_unstress_red_24", "stress_5-2_unstress_red_24",
    "stress_5-24_stress_green_48", "stress_5-27_stress_green_48",
    "stress_6-1_unstress_red_24", "stress_6-4_unstress_red_24",
    "stress_6-16_stress_green_48", "stress_6-18_stress_green_48",
    "stress_7-1_unstress_red_24", "stress_7-18_stress_green_48",
)
MORC_COLUMNS = (
    "morc_MORC_UT1", "morc_MORC_UT2", "morc_MORC_UT3",
    "morc_MORC_IAA1", "morc_MORC_IAA2", "morc_MORC_IAA3",
    "morc_MORC-KD-BFD1-KO_UT1", "morc_MORC-KD-BFD1-KO_UT2",
    "morc_MORC-KD-BFD1-KO_UT3", "morc_MORC-KD-BFD1-KO_IAA1",
    "morc_MORC-KD-BFD1-KO_IAA2", "morc_MORC-KD-BFD1-KO_IAA3",
    "morc_DD-BFD1-Ty_UT1", "morc_DD-BFD1-Ty_UT2", "morc_DD-BFD1-Ty_UT3",
    "morc_DD-BFD1-Ty_Shield1", "morc_DD-BFD1-Ty_Shield2", "morc_DD-BFD1-Ty_Shield3",
)
TOTAL_PROTEOME_COLUMNS = (
    "proteome_log2_normalized_and_imputated_abundances_UT_R1",
    "proteome_log2_normalized_and_imputated_abundances_UT_R2",
    "proteome_log2_normalized_and_imputated_abundances_UT_R3",
    "proteome_log2_normalized_and_imputated_abundances_T_24h_R1",
    "proteome_log2_normalized_and_imputated_abundances_T_24h_R2",
    "proteome_log2_normalized_and_imputated_abundances_T_24h_R3",
    "proteome_log2_normalized_and_imputated_abundances_T_32h_R1",
    "proteome_log2_normalized_and_imputated_abundances_T_32h_R2",
    "proteome_log2_normalized_and_imputated_abundances_T_32h_R3",
    "proteome_log2_normalized_and_imputated_abundances_T_48h_R1",
    "proteome_log2_normalized_and_imputated_abundances_T_48h_R2",
    "proteome_log2_normalized_and_imputated_abundances_T_48h_R3",
    "proteome_lfc_UT_Vs_T_24h_log2_fold_change",
    "proteome_lfc_UT_Vs_T_32h_log2_fold_change_1",
    "proteome_lfc_UT_Vs_T_48h_log2_fold_change_2",
)
QUANTITATIVE_PHOSPHO_COLUMNS = (
    "phospho_up_sites", "phospho_up_ratio", "phospho_down_sites",
    "phospho_down_ratio", "phospho_sites_measured",
)
OOCYST_ITRAQ_COLUMNS = (
    "oocyst_itraq_115_113", "oocyst_itraq_116_113", "oocyst_itraq_115_114",
    "oocyst_itraq_116_114", "oocyst_itraq_115_113_1", "oocyst_itraq_116_113_1",
    "oocyst_itraq_115_114_1", "oocyst_itraq_116_114_1",
)

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
                 "is circular by construction. Left unlabeled where no stage leads clearly."),

    # ------------------------------------------------------------------ reference
    Dataset("toxodb_identity", "ToxoDB gene identity", "reference", "identity",
            "Symbols, previous IDs, product descriptions", ("gene_id", "product"),
            "8,843 ME49 genes", accession="ToxoDB ME49", url=TOXODB,
            path="starplast/data/toxodb_identity.tsv",
            note="Retrieved 2026-08-11 via the REST API; strain tables for GT1 and VEG alongside."),
    Dataset("toxodb_strain_snps", "Strain variation (ToxoDB HTS SNPs)", "reference", "variation",
            "SNPs per gene across every sequenced strain, split by effect",
            ("snp_total_all_strains", "snp_nonsynonymous", "snp_synonymous", "snp_noncoding",
             "snp_stop_codon"), "8,140 genes (100%)", accession="ToxoDB ME49",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_strain_snps.tsv",
            note="Retrieved 2026-08-16 through the same REST report as the identity table, asking "
                 "for the five gene_hts_*_snps attributes. Verified against known biology rather "
                 "than against a metadata field: nonsynonymous SNPs per kb come out at 112 for the "
                 "SRS surface antigens, 55 for the ROP5/ROP18/GRA15 virulence loci, 30 across all "
                 "genes and 2.9 for ribosomal proteins. That ordering -- what immunity sees, then "
                 "the strain-typing markers, then the conserved core -- is the check. Zero is a "
                 "measurement here, not a gap: 690 genes carry no SNP in any sequenced strain."),
    Dataset("toxodb_codon_usage", "Codon usage bias (COMPUTED)", "reference", "sequence",
            "Effective number of codons, GC3, and codon adaptation index",
            ("codon_enc", "codon_gc3", "codon_cai_ribosomal"), "8,140 genes (100%)",
            accession="ToxoDB ME49",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_cds.tsv.gz",
            note="COMPUTED here from the coding sequences fetched 2026-08-16. ENC is Wright's "
                 "effective number of codons and GC3 the synonymous third-position GC, both "
                 "reference-free. CAI's reference set is the 158 ribosomal proteins, chosen by "
                 "product annotation and NOT by this map's expression columns -- the usual choice, "
                 "'the most highly expressed genes', would have built a sequence column out of an "
                 "expression column and then found them correlated. Verified by the signs "
                 "translational selection predicts: ribosomal proteins are more biased than the rest "
                 "(ENC 46.5 against 54.0), and CAI rises with transcription (rho +0.37) and with "
                 "protein abundance (rho +0.24) while ENC falls with both. Those correlations are a "
                 "finding here rather than a construction."),
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

    # ------------------------------------------------------------------ localization
    Dataset("lopit_tgon", "T. gondii hyperLOPIT", "post_translation", "LOPIT",
            "Subcellular compartment, MAP and MCMC, with posteriors",
            # Every column this experiment produces, including the ones `localization.py` computes
            # from it. Listing them all is not bookkeeping: `search.excluded_for` reads provenance to
            # decide what a map may SEE when the target came from here, and the five that were
            # missing are the five that leaked -- `lopit_prob_map`, `lopit_prob_mcmc` and
            # `lopit_methods_agree` together recover `compartment` better than any measurement block
            # in the cache.
            ("compartment", "compartment_best", "compartment_source", "lopit_map", "lopit_mcmc",
             "lopit_prob_map", "lopit_prob_mcmc", "lopit_methods_agree", "lopit_confident",
             "lopit_unified"),
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
            GSE108740_COLUMNS, "7,739 (95.1%)", accession="GSE108740",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108740"),
    Dataset("gse206344", "Oocyst sporulation series", "transcription", "RNAseq",
            "Unsporulated / sporulating / sporulated, 2 replicates (6 columns)",
            GSE206344_COLUMNS, "7,974 (98.0%)", accession="GSE206344",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206344"),
    Dataset("gse22258", "Pru tachyzoite / 72-hour bradyzoite stage array", "transcription",
            "microarray", "Matched tachyzoite and alkaline-induced bradyzoite expression",
            GSE22258_COLUMNS, "7,253 genes", accession="GSE22258",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22258",
            path="datasets/stagetranscriptome_GSE22258_series_matrix.txt.gz",
            note="Already keyed by TGME49 accessions. Kept separate from the newer RNA-seq stage "
                 "series; it is not averaged as though microarray intensity were FPKM."),
    Dataset("gse168465", "Primary brain-cell parasite differentiation time course",
            "transcription", "RNAseq",
            "Parasite base mean and log2 fold-change at days 1, 2, 4, 7 and 14",
            GSE168465_COLUMNS, "measured at build time", pmid="34610266", accession="GSE168465",
            citation="Mouveaux T et al., Primary brain cell infection by Toxoplasma gondii reveals "
                     "spontaneous bradyzoite differentiation and modification of neuron biology",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168465",
            path="datasets/stagetranscriptome_GSE168465_DESeq2-Toxo-all-time-points.xlsx",
            note="Dual host-parasite RNA-seq; only the workbook explicitly containing Toxoplasma "
                 "gene results enters this map. p-values remain evidence metadata, not features."),
    Dataset("gse99395", "Intracellular/extracellular ribosome profiling", "translation",
            "Ribo-seq", "Ribosome footprints, matched RNA and relative translation efficiency",
            GSE99395_COLUMNS, "measured at build time", pmid="29228904", accession="GSE99395",
            citation="Hassan MA et al., Comparative ribosome profiling uncovers a dominant role "
                     "for translational control in Toxoplasma gondii. BMC Genomics 2017;18:961",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99395/suppl/"
                "GSE99395_Raw_counts.txt.gz",
            path="datasets/toxoplasma_acquisition_2026_08_14/GSE99395_Raw_counts.txt.gz"),
    Dataset("gse129869", "Host-context parasite ribosome profiling", "translation", "Ribo-seq",
            "Parasite ribosome footprints, RNA and translation efficiency in two HFF states",
            GSE129869_COLUMNS, "measured at build time", pmid="31167946", accession="GSE129869",
            citation="Holmes MJ et al., Simultaneous Ribosome Profiling of Human Host Cells "
                     "Infected with Toxoplasma gondii. mSphere 2019;4:e00292-19",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129869/suppl/"
                "GSE129869_RAW.tar",
            path="datasets/toxoplasma_acquisition_2026_08_14/GSE129869_RAW.tar"),
    Dataset("gse19092", "Synchronized tachyzoite cell-cycle transcriptome", "transcription",
            "microarray", "Two replicates across blocked, asynchronous and hourly release states",
            GSE19092_COLUMNS, "measured at build time", pmid="20865045", accession="GSE19092",
            citation="Behnke MS et al., Coordinated progression through two subtranscriptomes "
                     "underlies the tachyzoite cycle of Toxoplasma gondii. PLoS ONE 2010;5:e12354",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19092/matrix/"
                "GSE19092_series_matrix.txt.gz",
            path="datasets/toxoplasma_acquisition_2026_08_14/GSE19092_series_matrix.txt.gz",
            note="Legacy GPL7186 probes are mapped through the platform ToxoDB field and the "
                 "project's previous-ID resolver."),
    Dataset("gse51780", "Feline merozoite transcriptome", "transcription", "microarray",
            "Merozoite expression with matched tachyzoite comparators", GSE51780_COLUMNS,
            "measured at build time", pmid="24885521", accession="GSE51780",
            citation="Behnke MS et al., Toxoplasma gondii merozoite gene expression analysis with "
                     "comparison to the life cycle. BMC Genomics 2014;15:350",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51780/matrix/"
                "GSE51780_series_matrix.txt.gz",
            path="datasets/toxoplasma_acquisition_2026_08_14/GSE51780_series_matrix.txt.gz"),
    Dataset("gse168155", "CPSF4 RNA-processing perturbation transcriptome", "transcription",
            "RNAseq", "RNA response at 7, 24 and 48 hours after CPSF4 depletion",
            GSE168155_COLUMNS, "measured at build time", pmid="34263725", accession="GSE168155",
            citation="Farhat DC et al., A plant-like mechanism coupling m6A reading to "
                     "polyadenylation safeguards transcriptome integrity. eLife 2021;10:e68312",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168155/suppl/"
                "GSE168155_Matrix_table_processed_data.xlsx",
            path="datasets/toxoplasma_acquisition_2026_08_14/"
                 "GSE168155_Matrix_table_processed_data.xlsx",
            note="This is a perturbation-response transcriptome, not a direct gene-wise m6A map."),
    Dataset("gse200962", "Bradyzoite restriction-checkpoint transcriptome", "transcription",
            "RNAseq", "Cyclin perturbations in tachyzoite and bradyzoite conditions",
            GSE200962_COLUMNS, "measured at build time", accession="GSE200962",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200962/suppl/"
                "GSE200962_gene_count_matrix_geo.csv.gz",
            path="datasets/toxoplasma_acquisition_2026_08_14/"
                 "GSE200962_gene_count_matrix_geo.csv.gz",
            note="No publication is linked from GEO; raw sample identifiers are retained verbatim "
                 "in column names rather than assigned conditions by guesswork."),
    Dataset("gse253884_5", "In-vivo fitness of hyperLOPIT-unassigned proteins", "DNA",
            "CRISPR_screen", "Two targeted libraries tested during mouse infection",
            GSE25388_COLUMNS, "measured at build time", pmid="39082802", accession="GSE253884;GSE253885",
            citation="Tachibana Y et al., CRISPR screens identify genes essential for in vivo "
                     "virulence among proteins of hyperLOPIT-unassigned localization. mBio 2024",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253884",
            path="datasets/toxoplasma_acquisition_2026_08_14/",
            note="Only the newly measured in-vivo fitness values enter; copied comparator columns "
                 "in the summary workbook are not duplicated."),
    Dataset("invivo_brain_transcriptome", "In vivo brain-stage transcriptome", "transcription",
            "RNAseq", "Tachyzoites, acute/chronic whole brain, and purified bradyzoites",
            INVIVO_BRAIN_COLUMNS, "7,663 (94.1%)", pmid="31726967",
            citation="Garfoot AL et al., Proteomic and transcriptomic analyses of early and "
                     "late-chronic Toxoplasma gondii infection shows novel and stage specific "
                     "transcripts. BMC Genomics 2019;20:859",
            url=SPRINGER.format(doi="s12864-019-6213-0",
                                f="12864_2019_6213_MOESM4_ESM.csv"),
            path="datasets/translation/proteomics/31726967/12864_2019_6213_MOESM4_ESM.csv",
            note="These columns are FPKM-derived transcript abundance, despite the mixed "
                 "transcriptome/proteome paper and the legacy proteomics directory. They belong "
                 "to transcription slots, never fitness."),
    Dataset("gse132248_stress", "Alkaline-stress differentiation transcriptome", "transcription",
            "RNAseq", "Unstressed tachyzoites and alkaline-stressed bradyzoites",
            STRESS_COLUMNS, "7,880 (96.8%)", pmid="31955846", accession="GSE132248",
            citation="Waldman BS et al., Identification of a Master Regulator of Differentiation "
                     "in Toxoplasma. Cell 2020;180:359-372.e16",
            url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132248",
            path="toxo_stage_atlas/data/transcriptomics/GSE132248_STAR_counts_matrix.tsv"),
    Dataset("morc_depletion", "MORC depletion and BFD1 perturbation transcriptome",
            "transcription", "RNAseq",
            "MORC knockdown, BFD1 knockout and BFD1 stabilization series", MORC_COLUMNS,
            "7,841 (96.3%)", accession="PXD058095",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD058095",
            path="toxo_stage_atlas/data/proteomics/"
                 "PXD058095_supp_DatasetEV1_MORC_RNAseq_counts_TPM.xlsx",
            note="The accession and file location sit in a proteomics collection, but the shipped "
                 "workbook is explicitly RNA-seq counts/TPM and is normalized as transcription. "
                 "The originating publication still needs confirmation before citation."),

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
            ("crispr_gra12s1_l2fc_invitro", "crispr_gra12s1_l2fc_invivo",
             "crispr_gra12s1_disco", "crispr_gra12s2_l2fc_invitro",
             "crispr_gra12s2_l2fc_invivo", "crispr_gra12s2_disco"), "236 / 232",
            pmid="40240328",
            citation="GRA12 is a common virulence factor across Toxoplasma gondii strains and "
                     "mouse subspecies",
            url=SPRINGER.format(doi="s41467-025-58876-2", f="41467_2025_58876_MOESM5_ESM.xlsx"),
            path="datasets/DNA/CRISPR_screen/40240328/",
            note="The two screens are NOT replicates: in-vivo L2FC correlate at r = 0.41."),
    Dataset("hosttx_effectors", "Host-transcription effector screen", "DNA", "CRISPR_screen",
            "Hotelling T2 plus full per-effector host-response signature",
            ("hosttx_T2", "hosttx_padj") + HOST_SIGNATURE_COLUMNS, "252 screened / 22 full signatures",
            pmid="37827122",
            citation="High-throughput identification of Toxoplasma gondii effector proteins that "
                     "target host cell transcription",
            url=EPMC.format(pmcid="PMC12033024"),
            path="datasets/DNA/CRISPR_screen/37827122/",
            note="The 737,726-row host differential-expression table is represented by 20 PCA "
                 "coordinates, its L2 norm and substantial-DE count; PCA is a dimensional summary, "
                 "not a host-gene measurement."),

    Dataset("bioid_corpus_membership", "BioID/TurboID supplement membership corpus",
            "post_translation", "proximity_labelling",
            "Number of downloaded proximity-labeling studies whose supplement names each gene",
            ("n_bioid_studies",), "measured at build time",
            url=EPMC.format(pmcid="{pmcid}"),
            path="datasets/post_translation/BioID/",
            note="Membership is not enrichment and is never converted to an interaction edge."),
    Dataset("ipms_corpus_membership", "IP-MS supplement membership corpus",
            "post_translation", "IP-MS",
            "Number of downloaded pulldown studies whose supplement names each gene",
            ("n_ipms_studies",), "measured at build time",
            url=EPMC.format(pmcid="{pmcid}"),
            path="datasets/post_translation/IPMS/",
            note="Membership is not enrichment and is never converted to an interaction edge."),

    # ------------------------------------------------------------------ protein level
    Dataset("proteome_pru", "Pru proteome and IP abundance", "translation", "proteomics",
            "Median log2 iBAQ across replicates", ("protein_ibaq_log2",), "748 (9.2%)",
            accession="PXD043808, PXD065585",
            url="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD065585",
            path="toxo_stage_atlas/data/proteomics/",
            note="Immunoprecipitation experiments of 424 and 594 proteins. Enrichment, NOT a deep "
                 "proteome; do not report as proteome-wide."),
    Dataset("proteome_total", "AP2XII-1/AP2XI-2 perturbation total proteome", "translation",
            "proteomics", "Replicate abundance and log2 fold-change during pre-sexual conversion",
            TOTAL_PROTEOME_COLUMNS, "3,005 (36.9%)", pmid="38093015",
            accession="PXD039400, PXD042658",
            citation="Antunes AV et al., In vitro production of cat-restricted Toxoplasma "
                     "pre-sexual stages. Nature 2024;625:366-376",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD039400",
            path="toxo_stage_atlas/data/proteomics/"
                 "PXD039400_PXD042658_supp_SupplTable3_total_proteome.xlsx",
            note="An actual total proteome, not IP enrichment. The abundance and fold-change "
                 "columns are distinct quantification types and expression.total_proteome "
                 "normalizes them separately before they enter the cache."),
    Dataset("phosphosites", "Phosphosite counts", "post_translation", "phosphoproteomics",
            "Count of phosphosites per protein, no positions", ("n_phosphosites", "has_phospho"),
            "1,175 (14.4%)",
            citation="Treeck M et al. 2011 -- CONFIRM against the file on disk",
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1931312811002885-mmc2.xls",
            note="The URL downloads a real phosphoproteomics table, but that it is the source of "
                 "THIS column is inference rather than verification; confirm before citing. "
                 "Missing for 85.6% of genes; effectively an indicator of having been in a "
                 "phosphoproteomics experiment."),
    Dataset("phospho_quantitative", "Oocyst-versus-tachyzoite phosphoproteome",
            "post_translation", "phosphoproteomics",
            "Measured-site counts and strongest up/down phosphosite ratios",
            QUANTITATIVE_PHOSPHO_COLUMNS, "1,603 (19.7%)", pmid="35164288",
            accession="PXD017032",
            citation="Wang Z-X et al., Comparative Phosphoproteomic Analysis of Sporulated "
                     "Oocysts and Tachyzoites of Toxoplasma gondii Reveals Stage-Specific "
                     "Patterns. Molecules 2022;27:1109",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD017032",
            path="toxo_stage_atlas/data/proteomics/"
                 "PXD017032_supp_TableS1_upregulated_phosphosites_oocyst_vs_tachy.xlsx",
            note="A per-gene cache cannot retain residue positions. It carries how many sites were "
                 "measured, how many moved each way and the median ratios; the source workbooks "
                 "remain the residue-level record."),
    Dataset("oocyst_itraq", "Oocyst developmental-stage iTRAQ proteome", "translation",
            "proteomics", "iTRAQ abundance ratios across oocyst developmental stages",
            OOCYST_ITRAQ_COLUMNS, "2,079 (25.5%)", pmid="28626452", accession="PXD003765",
            citation="Possenti A et al., Proteomic Differences between Developmental Stages of "
                     "Toxoplasma gondii Revealed by iTRAQ-Based Quantitative Proteomics. Front "
                     "Microbiol 2017;8:1732",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD003765",
            path="toxo_stage_atlas/data/proteomics/"
                 "PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls",
            note="The legacy XLS needs conversion before pandas can read it. Ratios are retained as "
                 "ratios rather than logged or centered, because moving their reference changes "
                 "the measurement."),

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
    Dataset("bioid_corpus", "Proximity-labeling corpus", "post_translation", "BioID",
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
            # The per-tier paper counts belong here too, and `attention_depth` is declared as what it
            # is: np.select over exactly those three columns (`literature.attention_depth`). Without
            # the declaration the NEGATIVE control was embedding its own inputs -- the counts score
            # 0.25-0.43 against the tiering, under any workable threshold, because a count is not a
            # restatement of a tier while determining it completely.
            ("n_publications", "n_papers_focal", "n_papers_substantive", "n_papers_incidental",
             "attention_depth"), "33,924 records",
            derived_from=("n_papers_focal", "n_papers_substantive", "n_papers_incidental"),
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

    Dataset("gse223620_bfd2_rip", "BFD2-bound transcriptome (RIP-seq, COMPUTED)", "transcription",
            "RIPseq", "Enrichment of each transcript in the BFD2 immunoprecipitation",
            ("bfd2_rip_log2_ip_over_input",), "7,463 genes (92%)", accession="GSE223620",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223620/suppl/"
                "GSE223620_ProcessedDataFile_BFD2.RIPseq.xls.gz",
            path="datasets/quarantine/2026_08_16_unverified/Tg/rna_binding_protein_targets/",
            note="COMPUTED: the deposit publishes read counts for the IP and the input, not the "
                 "ratio. Both are scaled to a common library size and the log2 taken, at a floor of "
                 "20 reads across the pair -- deliberately low, because the point of a RIP is the "
                 "enriched tail and a stricter floor would drop the genes the slot asks about. "
                 "Verified against the published mechanism: BFD2 binds and stabilises the BFD1 "
                 "transcript, and BFD1 comes out at +3.52, rank 32 of 8,090, the top 0.4%. BFD2's "
                 "own transcript is unremarkable at +0.31, which is what says the enrichment is not "
                 "an artefact of the tagged locus. One column and not a set: the slot holds one "
                 "protein's targets, and a second RIP would sit beside this rather than be averaged "
                 "into it."),
    Dataset("gse245775", "Differentiation ribosome profiling (eIF1.2)", "translation", "RiboSeq",
            "RPF and RNA counts, and their ratio, in tachyzoites and pre-bradyzoites",
            GSE245775_COLUMNS, "7,880 genes (97%)", pmid="38782906", accession="GSE245775",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245775/suppl/"
                "GSE245775_RAW.tar",
            path="datasets/quarantine/2026_08_16_unverified/Tg/stage_conversion_phenotype/",
            note="The `prebrady` arms are 48 hours in RPMI pH 8.3 at ambient CO2, which the "
                 "submitters label `cell type: pre-bradyzoites`. Named for that and not for "
                 "`bradyzoite`: these are not tissue cysts. Verified by reproducing the paper's own "
                 "result from the counts -- BFD1 rises 3.6 log2 on conversion in the parental line "
                 "and 2.6 in the knockout, and BFD2 rises 1.4 and 0.1, so the knockout's failure to "
                 "induce BFD2 is visible in the column itself. LDH2, BAG1 and SRS also rise on "
                 "conversion, which is what confirms the arms are not swapped. Reached the map "
                 "proposed for `stage-conversion phenotype`, which it is not."),

    Dataset("gse313048_atac", "Promoter accessibility (ATAC-seq, COMPUTED)", "DNA", "ATACseq",
            "Mean ATAC coverage over the promoter, relative to the genome mean",
            ("atac_promoter_ut",), "7,988 genes (98%)", accession="GSE313048",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313048/suppl/"
                "GSE313048_ATACseq_GCN5b-KD_UT.bw",
            path="datasets/quarantine/2026_08_16_unverified/Tg/acetylation/",
            note="COMPUTED here: GEO serves this as a bigWig with no peak calls and no per-gene "
                 "table, so the summary is the mean coverage 1 kb either side of the transcription "
                 "start, over the genome mean. No peak calling and no thresholds -- those would be "
                 "modelling choices invented here rather than taken from the authors. Only the "
                 "UNTREATED arm is read; the deposit is a GCN5b knockdown and the depleted arm "
                 "answers a different question. Verified against expression: rho +0.50, and the top "
                 "decile of expressed genes carries 2.0 log2 more promoter signal than the bottom "
                 "decile. Correlates with fitness at rho +0.02, so it is not merely tracking "
                 "essentiality."),
    Dataset("gse277553_cuttag", "HDAC3 occupancy (CUT&TAG, COMPUTED)", "DNA", "CUTandTAG",
            "Mean HDAC3 CUT&TAG coverage over the promoter, relative to the genome mean",
            ("cuttag_hdac3_promoter_ut",), "8,140 genes (100%)", accession="GSE277553",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277553/suppl/"
                "GSE277553_RAW.tar",
            path="datasets/quarantine/2026_08_16_unverified/Tg/chromatin_accessibility/",
            note="COMPUTED the same way as the ATAC column, from the untreated arm; the deposit's "
                 "other arm is an AP2XII-5 knockout. Two replicates and not three: "
                 "GSM8524430_UT_2.bw begins with eight 0xFF bytes and is not a bigWig, identically "
                 "whether taken from the series tar or fetched from GEO as a sample file, so the "
                 "corruption is in the deposit. It is skipped with a message rather than silently, "
                 "because a replicate dropped without saying so makes the mean smaller than the "
                 "note beside the column claims."),

    Dataset("toxodb_palmitome", "S-palmitoylome (Foe 2015, via ToxoDB)", "post_translation",
            "proteomics", "17-ODYA enrichment per gene, against hydroxylamine and against palmitate",
            ("palmitome_odya_vs_hydroxylamine_log2", "palmitome_odya_vs_palmitate_log2"),
            "470 and 488 genes", pmid="26468752", accession="ToxoDB Foe palmitome",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByProteomicsDirecttgonGT1_quantitativeMassSpec_Foe_Lipidome_Palmitoylome_RSRC"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_palmitome_hydroxylamine.tsv",
            note="The paper is not open access and PMC serves its supplementary spreadsheets only "
                 "through a download interstitial, so the numbers come from ToxoDB's own query "
                 "service for the same dataset -- the authors' fold differences, not a re-analysis. "
                 "Two comparisons and only one is palmitoylation: hydroxylamine cleaves thioester "
                 "bonds, which is the bond an S-palmitoyl group makes, so that column is "
                 "thioester-specific; the palmitate competition shows only that the label is "
                 "fatty-acid-dependent and includes N-myristoylated proteins. Verified against known "
                 "substrates: ROP5 +2.20, GAP45 +1.34, AMA1 +0.93, MLC1 +0.71, IMC proteins +0.34, "
                 "against a measured-gene median of -0.17. ToxoDB reports a SIGNED fold difference "
                 "and not a ratio -- -3.12 means three-fold down -- so reading it as a ratio would "
                 "have made every depleted protein NaN and dropped half the table."),

    Dataset("iedb_bcell", "Antibody epitopes (IEDB)", "reference", "immunity",
            "Distinct antibody epitope sequences per gene", ("n_bcell_epitopes",),
            "34 genes, 222 distinct epitopes", accession="IEDB bcell_search",
            url="https://query-api.iedb.org/bcell_search"
                "?parent_source_antigen_source_org_name=ilike.*Toxoplasma*",
            path="starplast/data/iedb_bcell_epitopes.tsv",
            note="From IEDB directly, NOT through ToxoDB, because ToxoDB's epitope integration is "
                 "not split by type and antigenicity is a question about antibodies. 189 of the 221 "
                 "genes in the ToxoDB epitope column have no antibody record at all, so the two are "
                 "genuinely different measurements filling different slots. Distinct epitope "
                 "SEQUENCES and not assay records: a protein studied by twenty groups accumulates "
                 "twenty records for one peptide, and counting records would rank antigens by "
                 "fashion. Verified by what tops it -- SRS29B (SAG1) 47, GRA6 29, GRA1 26, GRA4 20, "
                 "GRA7 16, MIC3 11, which is the panel commercial Toxoplasma serodiagnostic kits "
                 "are built from. Antigens are mapped by their product description, because IEDB "
                 "names them verbatim from ToxoDB; the trailing-symbol route resolves ten fewer and "
                 "loses SRS29B, the most studied antigen in the organism."),
    Dataset("toxodb_epitopes", "IEDB epitopes mapped to genes (via ToxoDB)", "reference",
            "immunity", "How many IEDB epitopes ToxoDB maps to this gene",
            ("iedb_epitope_count",), "221 genes", accession="ToxoDB / IEDB",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesWithEpitopes/reports/attributesTabular",
            path="starplast/data/toxodb_epitopes.tsv",
            note="ToxoDB's own join of IEDB against the ME49 proteome, at all three confidence "
                 "levels. NOT split by epitope type -- the integration does not expose that -- so "
                 "the count is T-cell and B-cell epitopes alike, which is why the column is named "
                 "iedb_ and not t_cell_. Verified by what comes out on top: SRS29B (SAG1) with 45, "
                 "then GRA6, GRA7, GRA2 and ROP18. Those are the canonical Toxoplasma serology "
                 "antigens, in the order a serologist would put them."),
    Dataset("toxodb_h4_acetylation", "Histone H4 acetylation (ChIP-chip, via ToxoDB)", "DNA",
            "ChIPchip", "Genome-wide H4 K5/K8/K12/K16 acetylation score within 1 kb of the gene",
            ("h4_acetylation_chip_score",), "7,515 genes (92%)",
            accession="ToxoDB Hakimi/Ali H4 acetylation",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByChIPchiptgonME49_chipChipExper_Hakimi_Ali_RSRC/reports/attributesTabular",
            path="starplast/data/toxodb_h4_acetylation.tsv",
            note="Verified as an ACTIVE mark must behave: rho +0.36 with transcription, +0.43 with "
                 "promoter ATAC, -0.02 with fitness, and genes in its top decile are expressed four "
                 "times as highly as those in its bottom. The Einstein H3K4me1 report from the same "
                 "site, the same assay type and the same query shape does the opposite -- its marked "
                 "genes have LESS accessible promoters -- and is in quarantine. This entry is the "
                 "counter-example that says that refusal is about the data and not about the reader."),
    Dataset("toxodb_macrophage", "Expression in infected macrophages (via ToxoDB)",
            "transcription", "RNAseq",
            "Expression percentile in ME49-infected murine macrophages",
            ("macrophage_expression_percentile",), "8,140 genes (100%)",
            accession="ToxoDB Saeij 29 strains",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByRNASeqtgonME49_Saeij_Jeroen_strains_rnaSeq_RSRCPercentile"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_macrophage.tsv",
            note="The ME49 arm of a 29-strain panel, so the strain matches the rest of the map. "
                 "Correlates at rho +0.89 with the fibroblast transcriptome (`expr_tachy`). That is "
                 "a FINDING -- the parasite's transcriptional programme is largely independent of "
                 "which host cell it is in -- and not a construction: it is an independent "
                 "measurement in a different host context, with nothing to declare in derived_from. "
                 "It is written down here so that nobody counts the two as independent evidence when "
                 "they agree, which they mostly will."),
    Dataset("toxodb_ec_numbers", "Enzyme classification (ToxoDB)", "reference", "annotation",
            "EC number per gene, and whether it has one", ("ec_number", "has_ec"),
            "1,313 enzymes of 8,140 genes", accession="ToxoDB ME49",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_ec_numbers.tsv",
            note="An annotation rather than a measurement, with the same standing as the InterPro "
                 "domains that fill `domain content`. It matters because it is the ONLY gene-indexed "
                 "metabolic datum there is -- every other metabolism question in the catalog is "
                 "about metabolites, and a metabolite is not a gene. Verified against conservation: "
                 "enzymes have a Plasmodium ortholog 59.6% of the time against 30.4% for other "
                 "genes (odds 3.39, p = 8e-88), are lineage-specific a third as often, and are more "
                 "costly to lose in vitro. `has_ec` is 0 and not missing where ToxoDB reports no EC: "
                 "the whole proteome was asked, so no assignment is an answer about the gene."),
    Dataset("crosslink_interactome", "Crosslinking MS interactome", "post_translation", "XLMS",
            "How many proteins this one crosslinks to", ("n_crosslink_partners",),
            "494 proteins", pmid="40874616", accession="mBio 02159-25 supplementary file s0004",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12505969/supplementaryFiles",
            path="starplast/data/crosslink_partners.tsv",
            note="395 high-confidence protein pairs, counted per protein. Verified against the two "
                 "largest obligate complexes any cell has: 21 of 32 proteasome subunits are in the "
                 "interactome (odds 30.8, p = 1e-18) and 57 of 158 ribosomal proteins (odds 9.7, "
                 "p = 4e-30). Crosslinking finds stable abundant complexes, and if it did not find "
                 "those two it would not be finding complexes. Absent is absent: a protein with no "
                 "partner here may be in no complex or may simply not have crosslinked."),
    Dataset("pvm_proximity", "PVM proximity labelling", "post_translation", "proteomics",
            "Whether the study placed this protein at the parasitophorous vacuole membrane",
            ("pvm_proximity_positive",), "1,274 genes (73 positive)", pmid="34749525",
            accession="mBio 00260-21 Data Set S1",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8576527/supplementaryFiles",
            path="starplast/data/pvm_proximity.tsv",
            note="One of the few binaries in the map with REAL zeros: the study publishes a "
                 "likely-negative list of 1,201 genes beside its 73 positives, so a zero here is a "
                 "measurement and not a gap -- for those 1,201 and only for them. Everything else "
                 "is NaN. Verified by what the positives are: 53 of 73 are dense granule proteins "
                 "against 0 of 1,201 negatives (Fisher p = 2e-77), and dense granule proteins are "
                 "exactly what Toxoplasma secretes into the vacuole and inserts into the membrane "
                 "it shares with the host cytosol."),
    Dataset("cdpk1_substrates", "CDPK1 substrates (thiophosphate labelling)", "post_translation",
            "proteomics", "Thiophosphorylated peptides per gene from analog-sensitive CDPK1",
            ("cdpk1_thiophospho_peptides",), "361 genes", pmid="37933960",
            accession="eLife 85654 supplementary file 6",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10629828/supplementaryFiles",
            path="starplast/data/cdpk1_substrates.tsv",
            note="The enriched fraction only, one gene per peptide from the master accession the "
                 "search engine assigned -- counting every protein a shared peptide maps to would "
                 "credit ambiguous peptides several times. The count tracks abundance, as every "
                 "phosphoproteomic count does, and the top of it is HSP70, HSP90 and BiP. What says "
                 "it is nonetheless CDPK1's substrate set is the enrichment: microneme proteins are "
                 "11-fold over-represented (Fisher p = 2e-05) and CDPK1 is the kinase that governs "
                 "microneme secretion, myosin A is in it, and so is the HOOK protein that the paper "
                 "exists to report."),
    Dataset("mrna_stability", "mRNA stability after actinomycin D", "transcription", "RNAseq",
            "Proportion of transcript remaining after five hours of transcription block",
            ("mrna_remaining_5h_actinomycin",), "412 genes", pmid="39899594",
            accession="PLoS Pathogens 1012857 Table S12",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11801735/supplementaryFiles",
            path="starplast/data/mrna_stability.tsv",
            note="A direct measurement: block transcription, wait, see what is left. Untreated "
                 "parasites at five hours, so the column is stability and not the iron response the "
                 "paper is about. BIASED BY CONSTRUCTION and the bias is worth stating -- the table "
                 "is the 426 transcripts that fell below 75% remaining, so it describes the unstable "
                 "tail and a gene absent from it is stable OR was not measured, which the column "
                 "cannot distinguish. Consistent with that: ribosomal-protein transcripts, which are "
                 "classically stable, are under-represented among the responders at odds 0.37."),
    Dataset("myristoylome", "N-myristoylated proteome", "post_translation", "proteomics",
            "The authors' confidence that this protein is myristoylated, 3 high to 1 low",
            ("myristoylation_confidence",), "65 substrates", pmid="32618271",
            accession="eLife 57861 supplementary file 4",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7373427/supplementaryFiles",
            path="starplast/data/myristoylome.tsv",
            note="A category the catalog did not have. N-myristoylation is co-translational and "
                 "irreversible, has its own enzyme in NMT and its own drug programme, and its "
                 "substrates are published. Verified against chemistry rather than against "
                 "annotation: myristoylation happens on an N-terminal glycine, and all 65 of 65 "
                 "substrates have glycine at position 2 against 5.8% of every other gene "
                 "(Fisher p = 4e-79). No other column in the map can be checked that cleanly. "
                 "Taken from the paper rather than from PXD019677, its PRIDE deposit, which ships "
                 "MaxQuant archives of 250-340 MB apiece; the answer is a 65-row table in "
                 "supplementary file 4."),
    Dataset("toxodb_arginine_methylation", "Monomethylarginine proteome (via ToxoDB)",
            "post_translation", "proteomics",
            "Monomethylarginine sites reported per gene", ("n_arginine_methylation_sites",),
            "368 genes", accession="ToxoDB Yakubu monomethylarginine",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByPTM"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_arginine_methylation.tsv",
            note="Yakubu et al.'s RH proteomics, through ToxoDB's PTM search. Verified by substrate "
                 "class rather than by a metadata field: RNA-binding, RRM and helicase proteins are "
                 "enriched 2.9-fold among the methylated (Fisher p = 0.002) and transporters and "
                 "membrane proteins are depleted at odds 0.43. That is the PRMT substrate profile -- "
                 "RG and RGG motifs sit in RNA-binding proteins. The slot it fills did not exist "
                 "before: arginine methylation has its own writers, its own substrate class and its "
                 "own deposit, and its absence from the catalog was a gap rather than a lack of "
                 "data."),
    Dataset("toxodb_nanopore_isoforms", "Novel transcript models (Nanopore, via ToxoDB)",
            "transcription", "LongRead",
            "How many novel TALON transcript models long reads support for this gene",
            ("novel_transcript_models",), "798 genes (10%)",
            accession="ToxoDB Stuart/Ralph nanopore",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByLongReadEvidence_tgonME49_Stuart_Ralph_nanopore_rnaSeqNextflow_RSRC"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_nanopore_isoforms.tsv",
            note="Incomplete-splice-match, novel-in-collection and novel-not-in-collection models, "
                 "at five supporting reads or more. `Known` is excluded because the annotated model "
                 "says nothing about isoform use, and `Genomic` because it is unspliced. Its power "
                 "is in PRESENCE rather than magnitude -- most genes that have a novel model have "
                 "one. Verified on that basis: genes with a novel model have a median of 6 exons "
                 "against 4 for genes without, Mann-Whitney p = 9e-43, which is the relationship "
                 "alternative splicing has to produce. Absent is NOT zero: a gene with no novel "
                 "model here may simply not have been sequenced deeply enough."),
    Dataset("toxodb_enteroepithelial", "Enteroepithelial stage transcriptome (via ToxoDB)",
            "transcription", "RNAseq",
            "Expression in the feline enteroepithelial stages against tachyzoites",
            ("ees_vs_tachyzoite_log2",), "7,739 genes (95%)",
            accession="ToxoDB Ramakrishnan enteroepithelial",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByRNASeqtgonME49_Ramakrishnan_enteroepithelial_stages_ebi_rnaSeq_RSRC"
                "/reports/attributesTabular",
            path="starplast/data/toxodb_enteroepithelial.tsv",
            note="EES1-5 averaged against tachyzoites, sense strand, no fold-change floor. The only "
                 "life stage in the map that happens inside a cat. Verified against stage markers: "
                 "GRA11B, which is merozoite-specific, comes out at +8.99 log2 and the family A/B/C "
                 "merozoite antigens at +2.66, against a genome median of -0.14. ToxoDB reports a "
                 "SIGNED fold difference, so -1257.4 means 1257-fold down and is converted rather "
                 "than logged."),

    # ------------------------------------------------------------------ PRIDE deposits
    # Counted from the submitters' own search output by `proteomics.deposit_counts`, never from the
    # raw spectra. Each column is "how many sites did THIS STUDY report on this gene", which is a
    # smaller claim than "how many sites does this gene have" and the only one the deposit supports.
    # Genes the study never saw are NaN rather than zero.
    Dataset("pride_acetylation", "Lysine acetylome (GCN5b)", "post_translation", "proteomics",
            "Acetylation sites reported per gene", ("n_acetylation_sites",),
            "3,921 genes measured", accession="PXD079431",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD079431",
            path="datasets/quarantine/2026_08_16_pride/Tg/acetylation/",
            note="Read from the deposit's own MaxQuant Sites tables. Held in quarantine rather than "
                 "the dataset archive: a deposit is promoted by being checked, not by being "
                 "downloaded, and the check here was reading the files and finding both Toxoplasma "
                 "genes and the modification named."),
    Dataset("pride_proximity", "Proximity labelling", "post_translation", "proteomics",
            "Proximity partners reported per gene", ("n_proximity_partners",),
            "1,734 genes measured", accession="PXD059579",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD059579",
            path="datasets/quarantine/2026_08_16_pride/Tg/interaction_proximity_labelling/",
            note="mzIdentML rather than MaxQuant, and shipped as a lone .gz -- which is one "
                 "compressed file and not an archive, a distinction that read as an empty deposit "
                 "until it was handled."),
    Dataset("pride_glycosylation", "O-fucosylated glycoproteins (AAL pulldown)",
            "post_translation", "proteomics",
            "Peptide identifications in the AAL lectin pulldown, per gene",
            ("n_o_fucosyl_peptides",), "394 genes", pmid=None, accession="PXD004426",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD004426",
            path="datasets/quarantine/2026_08_16_pride/Tg/glycosylation/",
            note="Aleuria aurantia lectin affinity purification, so the count is peptide "
                 "IDENTIFICATIONS and not sites -- the same standing as the proximity-labelling "
                 "column, which is also a claim about what came down rather than about a residue. "
                 "Verified against compartment: the pulldown is enriched for nucleus-chromatin at "
                 "odds 3.71 (p = 1.5e-22) and cytosol at 2.61, and not enriched for mitochondrion. "
                 "O-fucosylation through SPY is a nucleocytoplasmic modification and the paper "
                 "describes punctiform signal beside the nuclei, so that is the right answer. Keyed "
                 "on TGGT1_ accessions throughout."),
    Dataset("pride_nitrosylation", "S-nitrosylation (iodoTMT)", "post_translation", "proteomics",
            "S-nitrosylation sites reported per gene", ("n_nitrosylation_sites",),
            "660 genes measured", accession="PXD046083",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD046083",
            path="datasets/quarantine/2026_08_16_pride/Tg/S_nitrosylation/",
            note="Counted only from the iodoTMT tables. Counting the whole txt folder put 90% of the "
                 "proteome in this slot, which is what a modification measured on nearly every gene "
                 "should always look like: a bug."),
    Dataset("pride_lactylation", "Lysine lactylome", "post_translation", "proteomics",
            "Lactylation sites reported per gene", ("n_lactylation_sites",),
            "515 genes measured", accession="PXD031526",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD031526",
            path="datasets/quarantine/2026_08_16_pride/Tg/lactylation/",
            note="The second lactylation deposit tried. PXD022700 ships a RAR that bsdtar cannot "
                 "open, and the 0 genes that produced was a fact about the reader, not about a study "
                 "that names 537 proteins. This one reads, and its `La (K)Sites` table is MaxQuant's "
                 "lactylation search. Keyed entirely on TGGT1_ accessions -- 515 of its 524 genes "
                 "reach the map through the identity layer and would reach none without it."),
    Dataset("pride_ubiquitination", "Ubiquitination / SUMOylation (GlyGly)", "post_translation", "proteomics",
            "GlyGly sites reported per gene", ("n_ubiquitination_sites",),
            "128 genes measured", accession="PXD042937",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD042937",
            path="datasets/quarantine/2026_08_16_pride/Tg/ubiquitination_SUMOylation/"),

    Dataset("oxidative_stress_screen", "Oxidative-stress CRISPR screen", "DNA", "CRISPR_screen",
            "Screening score per gene under oxidative challenge",
            ("oxidative_stress_screen_score",), "7,384 genes (91%)", pmid="34163449",
            accession="PMC8216390 Data Sheet 1",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8216390/supplementaryFiles",
            path="datasets/quarantine/2026_08_16_unverified/Tg/fitness_oxidative_stress/",
            note="The authors' own `Screening score` sheet, not a recomputation from the guide "
                 "counts beside it in the same workbook -- they published the score, so it is "
                 "theirs to define. Negative is required. Verified by the sign and the extreme: "
                 "catalase comes out at -6.15, essentially the bottom of the whole screen, and it "
                 "is the enzyme that disposes of hydrogen peroxide. Peroxiredoxin (-1.57), "
                 "superoxide dismutase (-1.01), glutaredoxin (-0.71) and thioredoxin (-0.63) all "
                 "sit below the genome median of -0.38."),
    # UNPUBLISHED. Manuscript under submission; this ships inside the data cache, so it is the first
    # thing to remove before any package release. Flagged here rather than only in a note because
    # `datasets.registry()` is what a release check would read.
    Dataset("spacr_escrt_screen", "Host ESCRT recruitment screen (UNPUBLISHED)", "DNA",
            "CRISPR_screen",
            "Per-gene effect on host TSG101 recruitment to the vacuole, by two models",
            ("escrt_recruitment_xgboost", "escrt_recruitment_maxvit"),
            "13 and 8 genes", accession="spaCR screen, bioRxiv 10.64898/2026.07.08.737057",
            citation="Olafsson EB et al., A pooled image-based CRISPR screen identifies EAF1 as a "
                     "T. gondii modulator of ESCRT subversion. bioRxiv 2026 (under submission)",
            url="https://doi.org/10.64898/2026.07.08.737057",
            path="datasets/quarantine/2026_08_16_unverified/Tg/escrt_screen/",
            note="A pooled image-based screen of secretory proteins, deconvolved to gene effects by "
                 "regression. The phenotype is host ESCRT recruitment and NOT invasion or egress, "
                 "so it fills a slot of its own rather than the invasion slot whose "
                 "high-content-imaging context it matches -- a slot names a question, not a method. "
                 "Both deconvolution models are shipped because their agreement is the "
                 "verification, and the agreement is close: TGGT1_244480 is rank 1 in each, "
                 "TGGT1_409250 and GRA14 fill out the top three of each, and EAF1 -- which the "
                 "PRIDE deposit PXD080696 identifies as TGGT1_225160 -- is rank 7 in BOTH. Two "
                 "independent deconvolutions landing the same gene at the same rank is a stronger "
                 "statement than any single ordering. An earlier version of this note called "
                 "TGGT1_244480 EAF1 and said EAF1 was rank 1; both were wrong, and the deposit's "
                 "own abstract is what settled it. "
                 "The MYR1 host bridge added the same day is verified by ESCRT machinery topping "
                 "it, so two unrelated datasets in this map now point at the same biology."),
    Dataset("myr1_host_ip", "MYR1 host interactome (bridge)", "post_translation", "IPMS",
            "Host proteins co-immunoprecipitating with the parasite protein MYR1",
            ("bridge:host",), "219 host proteins, 1 parasite gene", pmid="32075880",
            accession="PXD016383",
            url="https://www.ebi.ac.uk/pride/archive/projects/PXD016383",
            path="starplast/data/host_bridges.parquet",
            note="The first BRIDGE in the project: pairs whose two ends are in different tables, "
                 "which `graph.npz` cannot hold because its edges are index pairs into the parasite "
                 "table. Verified twice over. The parasite side: MYR1 is rank 1 of 325 in its own "
                 "pulldown at +34.4, with MYR3 at 47 and GRA44, GRA7 and GRA9 in the top ten. The "
                 "host side against independent literature: PDCD6, which is ALG-2, is rank 1 of "
                 "219, and a 2026 paper reports Toxoplasma GRA8 engaging host ALG-2 at the vacuole; "
                 "its partner ALIX is rank 12 and VPS28 rank 62, so ESCRT as a class sits at "
                 "p = 0.006 against the rest of the host proteins. Two filters do the work and both "
                 "were got wrong first: a group is a contaminant group if ANY entry in it is one -- "
                 "MaxQuant prefixes only on the leading entry, so keratin hides mid-group and is "
                 "otherwise the four most enriched host proteins -- and two unique peptides are "
                 "required in BOTH bait replicates, which takes 674 host groups to 219."),
    Dataset("pv_host_uptake", "Host proteins at the vacuole", "post_translation", "proteomics",
            "How enriched a host protein is at the parasitophorous vacuole",
            ("pv_enrichment_log2",), "12 host proteins", pmid="34898650",
            accession="PLoS Pathogens 1010138 supplementary table",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8700025/supplementaryFiles",
            path="starplast/data/host_proteins.parquet",
            note="A property OF a host protein rather than a bridge, because the bait is the "
                 "compartment and not a named parasite gene -- a bridge needs a parasite gene at one "
                 "end. Averaged over three infection contexts: tachyzoite-infected fibroblast, "
                 "bradyzoite-infected fibroblast and neuron. The sheet lists parasite and host "
                 "proteins together, which is how the authors show the experiment worked -- the "
                 "dense granule proteins top it -- and only the host rows are kept. Top of those: "
                 "PDCD6/ALG-2 at +5.90, VPS37C at +4.47, then CHMP4B, PEF1 and VPS28. That is the "
                 "FOURTH independent dataset in this map to put ALG-2 at the host-parasite "
                 "interface, after the MYR1, EAF1 and GRA35 pulldowns."),
    Dataset("metabolome_iron", "Metabolome and isotope labelling under iron deprivation",
            "reference", "metabolomics",
            "Steady-state metabolite levels and the fraction labelled from glucose or glutamine",
            ("metabolite_level_log2fc_iron_depleted", "metabolite_level_padj",
             "labelled_fraction_glucose", "labelled_fraction_glutamine"),
            "1,102 metabolites", pmid="41925342",
            accession="mBio 03788-25 Tables S3 and S5",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13170339/supplementaryFiles",
            path="starplast/data/metabolites.parquet",
            note="Rows are COMPOUNDS, not genes -- the first table in the project that is not the "
                 "node table, and the shape instruction 39 describes for host tables. The flux "
                 "column is one minus the unlabelled isotopologue, which is the only labelling "
                 "readout comparable between molecules of different carbon number. The archive's "
                 "two labelling sheets differ in whether they carry a title row, and assuming they "
                 "did not silently dropped the glucose arm: the sheet read fine and had no column "
                 "called `Metabolite`. Joining a second study means matching compound NAMES, which "
                 "is lossy; that cost is unpaid with one study and is the first thing to fix when a "
                 "second arrives."),
    Dataset("lipidome_vesicles", "Membrane lipid composition of parasite vesicles",
            "reference", "lipidomics",
            "Lipid species abundance, and its proportion against the host cell",
            ("lipid_ev_level_log2", "lipid_ev_vs_host_clr"),
            "194 lipid species", pmid="41716462",
            accession="Front Cell Infect Microbiol 1745625 Tables 1-3",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12913473/supplementaryFiles",
            path="starplast/data/metabolites.parquet",
            note="The SECOND study in the metabolite table, and the name-matching cost the entry "
                 "above warns about is not paid here: lipid shorthand and polar compound names are "
                 "different naming systems and none of the 194 species collides with the 1,102 "
                 "compounds. Rows are appended, not joined. These are read as PARASITE lipids "
                 "because the vesicles came from post-egress tachyzoites in host-cell-free medium, "
                 "and because the composition does not move when the host does: across four host "
                 "backgrounds the host cells differ in 1,018-1,362 species and the vesicles the "
                 "same parasite released in them differ in 0-4. A lipidome of an infected culture "
                 "would not have supported this slot at all -- most of that lipid is host. The "
                 "`vs_host` column is COMPUTED here, sample-centred so it compares proportion "
                 "rather than amount; the archive's own EV-minus-cell column is not used because "
                 "its transform could not be reproduced to better than 3 log units."),
    Dataset("splitcas9_imaging_screen", "Arrayed splitCas9 imaging screen", "DNA",
            "imaging_screen",
            "What a parasite looks like when a gene is off: egress, actin, apicoplast, replication",
            ("screen_egress_phenotype", "screen_actin_phenotype", "screen_apicoplast_phenotype",
             "screen_replication_phenotype", "screen_any_phenotype", "screen_scorers_agree"),
            "319 genes screened, 99 with a phenotype, 35 at egress", pmid="35538310",
            accession="Nat Microbiol 41564-2022-01114 Supplementary Tables 2 and 3B",
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-022-01114-y/"
                "MediaObjects/41564_2022_1114_MOESM4_ESM.xlsx",
            path="datasets/DNA/imaging_screen/35538310/41564_2022_1114_MOESM4_ESM.xlsx",
            note="The only per-gene invasion-or-egress phenotype table published for Toxoplasma, "
                 "and it covers the EGRESS half: the screen's own figure legend calls it a screen "
                 "for actin dynamics, apicoplast segregation and egress, and invasion is a property "
                 "its hits were shown to have afterwards rather than a category anything was scored "
                 "into. Two earlier passes read this slot as blocked because the other candidates "
                 "promise invasion and egress in their titles and deliver it by characterising one "
                 "gene. The category codes had to be EARNED: the workbook ships no legend, it lives "
                 "in a figure that is an image, so `E` is read as egress because the paper names "
                 "exactly four categories in three places and because the two genes it names as its "
                 "egress mutants, CGP TGGT1_240380 and SLF TGGT1_208420, both carry an E -- which "
                 "the test asserts, so the mapping can fail. The SUBSCRIPT is deliberately not "
                 "read: E3 and E4 differ in something no accessible text defines, and a severity "
                 "invented from a digit is a number with no measurement behind it. Missingness "
                 "carries the other half of the meaning -- a screened gene with no egress call was "
                 "looked at and was normal, and the 7,800 unscreened genes stay missing, because "
                 "collapsing those would tell the map that nearly every gene has been checked."),
    Dataset("plasmodb_pf3d7_attributes", "Plasmodium falciparum 3D7 gene attributes",
            "reference", "annotation",
            "The second species: sequence, orthology, domains, strain SNPs and piggyBac fitness",
            ("length", "molecular_weight", "isoelectric_point", "transcript_length", "exon_count",
             "n_tm", "is_tm", "has_signal_peptide", "ortholog_number", "orthogroup",
             "paralog_number", "has_paralog", "n_interpro", "has_domain", "interpro_ids",
             "pfam_ids", "snp_total_all_strains", "snp_nonsynonymous", "snp_synonymous",
             "snp_noncoding", "snp_stop_codon", "snp_nonsyn_syn_ratio",
             "piggybac_mis", "piggybac_mfs"),
            "5,720 P. falciparum genes",
            accession="PlasmoDB GenesByTaxon attributesTabular",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByTaxon/reports/attributesTabular",
            path="starplast/data/pf_nodes.parquet",
            note="The FIRST parasite table beyond Toxoplasma, and it is its own table -- nothing is "
                 "merged, because PF3D7 and TGME49 identifiers do not align and neither does the "
                 "data behind them, so a merged frame would encode which species was convenient to "
                 "work on as though it were biology. Column names deliberately match the "
                 "Toxoplasma table where the quantity is the same, so a slot pattern reads the "
                 "same on both arms, but the Pf slots name their own patterns rather than "
                 "inheriting Toxoplasma's -- inheriting them would have claimed falciparum slots "
                 "with gondii numbers. The report is served by the TRANSCRIPT record type, so "
                 "5,791 rows describe 5,720 genes and the rows are collapsed on the longest "
                 "transcript; taking the row count at face value double-weights 71 genes. The "
                 "piggyBac direction is the other trap: a LOW mutagenesis index means the gene "
                 "resists disruption and is therefore essential, and inverting it would swap the "
                 "essential and dispensable genomes without crashing, so the test checks it "
                 "against biology -- ribosomal proteins come out at median MIS 0.15 and the var, "
                 "rifin and stevor families at 0.94."),
    Dataset("cotranslation_edges", "Co-translation layer (COMPUTED)", "translation", "RiboSeq",
            "Gene pairs whose ribosome footprints covary", ("edge:cotranslation",),
            "6,231 edges over 7,437 genes",
            derived_from=("rpf99395_intracellular_r1", "rpf129869_confluent_r1",
                          "rpf245775_parent_tachy_r1"),
            note="COMPUTED here: the same construction co-expression uses over transcripts, applied "
                 "to the 22 ribosome-footprint columns, at r >= 0.95 and the top 25 neighbours. "
                 "DERIVED, so it is declared -- an embedding built on the RPF columns must not then "
                 "be validated against this layer. It is not a copy of co-expression: 97% of its "
                 "edges are not co-expression edges, Jaccard 0.003. What says it is co-TRANSLATION "
                 "is that ribosomal proteins pair with each other 464 times where chance gives 3; "
                 "they are made together stoichiometrically, which is the textbook case of "
                 "co-translational regulation. 85 of its edges are also measured crosslink "
                 "contacts."),
    Dataset("m6a_peaks", "m6A methylome (MeRIP peaks)", "transcription", "MeRIP",
            "How many m6A peaks the authors called on this gene in tachyzoites",
            ("n_m6a_peaks",), "837 genes (10%)", pmid="34324585",
            accession="PLoS Pathogens 1009335 Table S3A",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8354455/supplementaryFiles",
            path="starplast/data/m6a_peaks.tsv",
            note="GEO carries no MeRIP for Toxoplasma -- its deposit for this study is depletion "
                 "RNA-seq, which says which transcripts DEPEND on m6A and not which CARRY it. The "
                 "peaks are in the paper. This slot had already been written up as unservable in "
                 "instruction 41 when the supplement turned up, and that entry is now struck "
                 "through rather than deleted. 866 of 8,922 genes carry a peak, which is the right "
                 "order for m6A. Verified against the paper's own second dataset: marked genes are "
                 "enriched among those responding to METTL3 depletion, odds 1.35, p = 2e-03 -- "
                 "modest because removing a writer has broad indirect effects, but the direction a "
                 "writer's own substrates have to take."),
    Dataset("sexual_stages", "Sexual development in the cat (single-cell atlas)", "transcription",
            "scRNAseq", "Enrichment at 8 days post-infection, when gametogony happens",
            ("sexual_stage_8dpi_log2fc",), "4,463 genes (55%)", pmid="41929010",
            accession="PMC13042011 supplementary media-2",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13042011/supplementaryFiles",
            path="starplast/data/sexual_stages.tsv",
            note="The only stage in the map that is both inside a cat and sexual; the "
                 "enteroepithelial column holds the asexual stages that precede it. Verified by "
                 "what rises: oocyst wall protein sits at the 97th percentile, and the oocyst wall "
                 "is built at the end of the sexual cycle, while ribosomal housekeeping genes sit "
                 "at the 23rd. Accessions arrive as `DEAD/DEAHboxhelicase-TGME49-220860` -- product "
                 "description glued to the accession with hyphens for underscores -- so they are "
                 "extracted and normalised rather than matched."),
    Dataset("secretome_partition", "Secreted-fraction partition", "post_translation", "proteomics",
            "How a secreted protein splits between the soluble and vesicular fractions",
            ("secretome_soluble_over_vesicle_log2",), "165 proteins", pmid="40874616",
            accession="ToxoDB Ramirez-Flores vesicles",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByProteomicstgonGT1_quantitativeMassSpec_Ramirez_Flores_Vesicles_RSRC"
                "/reports/attributesTabular",
            path="starplast/data/secretome_partition.tsv",
            note="REFUSED earlier in the same session and then accepted, which is worth stating "
                 "rather than reversing quietly. It was first asked which proteins are enriched in "
                 "secreted VESICLES, and the answer was incoherent -- dense granule proteins "
                 "depleted, ribosomal proteins the most enriched thing in it. The two fractions it "
                 "compares are BOTH secreted material, so the question it can answer is how a "
                 "secreted protein partitions between them, and asked that way it behaves: "
                 "micronemes, which dominate classical excretory-secretory antigen preparations, at "
                 "+3.72, and the GPI-anchored surface antigens at -0.68. Its limit is that there is "
                 "no negative list -- 171 proteins were seen in secreted material and nothing says "
                 "what was looked for and missed, so absence is not evidence."),
    Dataset("antisense_level", "Antisense transcription (via ToxoDB)", "transcription", "RNAseq",
            "Percentile of antisense signal at this gene, across the life cycle",
            ("antisense_expression_percentile",), "8,140 genes (100%)",
            accession="ToxoDB full life-cycle transcriptome, Antisense",
            url="https://toxodb.org/toxo/service/record-types/transcript/searches/"
                "GenesByRNASeqtgonME49_tgme49_spor_ocyst_rnaseq_ebi_rnaSeq_RSRCPercentile"
                "/reports/attributesTabular",
            path="starplast/data/antisense_level.tsv",
            note="The LEVEL of antisense transcription and not its change. An earlier attempt at "
                 "this slot used the sense/antisense CHANGE across a tachyzoite time course and was "
                 "refused because two strain time courses of it shared 12 of their top 200 genes "
                 "where chance gives 28. The level reproduces: against the independent Gregory ME49 "
                 "series it shares 197 of its top 500 where chance gives 31, a six-fold enrichment. "
                 "How much antisense a gene has is a property of the gene; how much it changed was a "
                 "property of the run. Correlates with sense transcription at only rho = +0.23, so "
                 "it is not a restatement of expression."),
    Dataset("melting_temperature", "Protein melting temperature (mineCETSA)", "post_translation",
            "proteomics", "Where this protein's melting curve sits, in degrees",
            ("melting_temperature_tm",), "3,120 proteins (38%)", pmid="35976251",
            accession="eLife 80336 supplementary file 3",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9436416/supplementaryFiles",
            path="starplast/data/melting_temperature.tsv",
            note="The same file as the thermal-shift column and a different question: that one asks "
                 "how far the curve MOVES when calcium is added, this asks where it SITS. Median "
                 "55.1 C, which is where protein melting temperatures live. Curves with R2 below "
                 "0.8 or a Tm outside 30-80 C are dropped -- the fit reports values up to 8,563, "
                 "which is a failed fit and not a thermophile. Verified by reproducing across "
                 "independent replicates at rho = +0.78 over 1,623 proteins; a Tm that did not "
                 "reproduce would be describing the run."),
    Dataset("thermal_shift_cetsa", "Calcium thermal-shift proteome (mineCETSA)", "post_translation",
            "proteomics", "How far a protein's melting curve moves when calcium is added",
            ("cetsa_calcium_ed_score",), "2,348 proteins", pmid="35976251",
            accession="PMC9436416 Supplementary file 3",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9436416/supplementaryFiles",
            path="datasets/quarantine/2026_08_16_unverified/Tg/thermal_shift/",
            note="The authors' own Euclidean-distance score, not recomputed from the ten "
                 "temperature points published beside it. Verified against the calcium sensors: "
                 "CAM1 and CAM2 sit at the 98th percentile and CAM3 at the 83rd, and a protein "
                 "whose melting curve does not move when calcium is added is not calcium-binding. "
                 "One caveat for anyone comparing against the paper -- its headline conclusion is "
                 "about PP1, and PP1 is unremarkable in THIS column. That claim comes from the "
                 "zaprinast time course in the same paper, a different experiment; this is the "
                 "calcium mineCETSA sheet. PXD033642, the deposit for the same study, publishes "
                 "only identifications and could not have filled this slot."),
    Dataset("cyst_wall_interactome", "Cyst wall interactome", "post_translation", "IPMS",
            "Strongest bait signal and how many baits saw the protein",
            ("cyst_wall_max_spectral", "cyst_wall_n_baits"), "56 proteins", pmid="32019789",
            accession="PMC7002340 Data Set S1",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7002340/supplementaryFiles",
            path="datasets/quarantine/2026_08_16_unverified/Tg/cyst_wall/",
            note="Two numbers because the bait count is the more honest one: a protein found by "
                 "thirteen independent pulldowns is in the cyst wall in a way a single strong hit "
                 "is not. Verified by what comes out on top -- MAG1 and MAG2, the canonical cyst "
                 "matrix proteins, with MCP3, MCP4 and SRS44 beside them. 57 of the table's 265 "
                 "rows name a Toxoplasma accession; the rest are human, because the pulldowns were "
                 "done on infected cultures and the table lists everything identified."),

    # ------------------------------------------------------------------ differentiation
    Dataset("second_background_fitness", "Fitness in the reporter strain (COMPUTED)", "DNA",
            "CRISPR_screen", "Guide depletion over eight passages of ordinary growth",
            ("crispr_reporter_strain_p8_log2",), "262 genes", accession="GSE132237",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132237/suppl/"
                "GSE132237_RAW.tar",
            path="datasets/quarantine/2026_08_16_pride/Tg/"
                 "essentiality_in_a_second_background/",
            note="The SAME archive as the differentiation screen, answering its other question. Its "
                 "passage arms -- p8 against the input library -- are ordinary tachyzoite growth, "
                 "and the reporter line is not the type I RH the genome-wide screens use, which is "
                 "what makes this a second background rather than a repeat. A targeted library of "
                 "nucleic-acid binding genes, so 262 genes and not the genome. Verified by agreeing "
                 "with the RH screen where it should: rho = +0.62 against fit_invitro_hff over 130 "
                 "shared genes, close enough that the direction and the join are right and far "
                 "enough that it is not a copy. The deposit sat in a folder named for this slot all "
                 "day while only its differentiation arms were read."),
    Dataset("differentiation_screen", "Differentiation reporter CRISPR screen (COMPUTED)",
            "DNA", "CRISPR_screen",
            "Guide enrichment in reporter-positive parasites against the bulk population",
            ("diff_reporter_log2_mNG_over_bulk",), "235 genes",
            accession="GSE132237",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132237/suppl/"
                "GSE132237_RAW.tar",
            path="datasets/quarantine/2026_08_16_pride/Tg/"
                 "essentiality_in_a_second_background/GSE132237_RAW.tar",
            note="COMPUTED here: the deposit publishes guide counts, not the ratio. Guides are summed "
                 "per gene, scaled to a common library size, and log2((mNG+ + 1)/(bulk + 1)) is "
                 "averaged over the two reporter lines at 10 days. It is the comparison the authors' "
                 "design names, and the column says ratio rather than phenotype so nobody mistakes "
                 "it for a number they reported. A targeted screen against nucleic-acid binding "
                 "proteins, so 240 genes is its full extent and not a coverage failure. Which member "
                 "is which sample comes from the series matrix, never from the file name."),
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


def derived_dependents(columns) -> tuple:
    """Transitive derived outputs that depend on any supplied column.

    Leakage travels both directions in the dependency graph: holding out a derived target removes its
    raw inputs, while holding out a raw measurement must also remove every summary or label computed
    from it.  The registry is small, so an explicit fixed-point walk is clearer than a cached graph.
    """
    seen = set(columns if not isinstance(columns, str) else (columns,))
    changed = True
    while changed:
        changed = False
        for dataset in REGISTRY:
            if dataset.derived_from and set(dataset.derived_from) & seen:
                new = set(dataset.columns) - seen
                if new:
                    seen.update(new)
                    changed = True
    return tuple(sorted(seen - set(columns if not isinstance(columns, str) else (columns,))))


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
    return False, "URL is not recognizably a file"


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
