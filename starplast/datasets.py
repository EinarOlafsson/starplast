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
    # GEO deposits whose processed tables sit under each SAMPLE rather than under the series: the
    # suffix that names them. Declared rather than inferred, because "fetch every per-sample file"
    # would pull the coverage tracks too -- 139 MB of wiggle for 250 kB of per-gene numbers.
    geo_file_suffix: str | None = None


# Springer and PLOS serve supplementary directly; PMC's /bin/ path 404s.
SPRINGER = "https://static-content.springer.com/esm/art%3A10.1038%2F{doi}/MediaObjects/{f}"
PLOS = "https://journals.plos.org/plospathogens/article/file?id=10.1371/{doi}.{s}&type=supplementary"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
TOXODB = ("https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
          "/reports/attributesTabular")
PLASMODB = ("https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon"
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
    Dataset("curated_enteric_fitness", "Enteric / sexual-cycle fitness per gene (CURATED)",
            "reference", "literature",
            "Gene disruptions carried through the feline stage with oocyst output measured",
            ("enteric_oocyst_yield", "enteric_sporulation", "enteric_measurements"),
            "8 genes, 4 studies", pmid="28288194",
            accession="PMIDs 28288194, 30728393, 36809045 and PMC12942651",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5363998/fullTextXML",
            note="THIRD curated source, and it exists because the sweeps that closed this slot "
                 "looked for a POOLED SCREEN through the enteroepithelial stages and correctly "
                 "found none -- nobody has put a barcoded library through a cat. But the slot asks "
                 "whether disrupting a gene costs the parasite oocysts, and feeding one knockout to "
                 "a cat answers that one gene at a time. Two papers were found and REFUSED for "
                 "failing the bar: one says the cat experiment 'should be carried out', one says "
                 "oocysts were seen but the numbers were 'not quantified'. Four of the ten rows are "
                 "unchanged, and they are the strongest rows here -- deleting all four LEA genes at "
                 "once left oocyst yield alone (30 against 34 million from paired kittens), which a "
                 "single knockout could not have established because redundancy could have hidden "
                 "it. Yield and sporulation are separate columns because HAP2 sheds a few oocysts "
                 "that never sporulate while Grx5 sheds fewer that sporulate poorly, and those are "
                 "different events. Magnitudes stay in the evidence text: they are not comparable "
                 "across cats, strains and inocula, and one numeric column would invent a precision "
                 "the experiments do not have."),
    Dataset("curated_drug_sensitivity", "Drug sensitivity per gene (CURATED)", "reference",
            "literature",
            "Knockouts with a measured shift in sensitivity to a named compound",
            ("drug_compounds_tested", "drug_sensitivity_shifts", "drug_sensitivity_directions"),
            "3 genes, 2 compounds", pmid="41025776", accession="mBio, PMID 41025776",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12607627/fullTextXML",
            note="SECOND curated source, and it exists because six catalogue sweeps looked for a "
                 "genome-wide chemogenomic screen and correctly found none -- while the sentence they "
                 "tested, 'no screen exists', is not the same statement as 'the question cannot be "
                 "answered'. That is the same error the lipid slot took three passes to notice. The "
                 "bar is a MEASURED shift under a NAMED compound, not 'an inhibitor of this protein "
                 "kills the parasite', which is target engagement and has its own slot; a test "
                 "asserts no row's evidence reads that way. `unchanged` rows are KEPT -- a "
                 "transporter deleted with no effect on analog sensitivity is a result, and dropping "
                 "those would leave the column looking like a list of hits. TgENT3's two rows were "
                 "measured in a ΔTgAT1 background and say so, because reading a double mutant's "
                 "phenotype off one of its genes is its own error. The per-row product check earned "
                 "its keep immediately: the annotation calls TGME49_244440 'adenosine transporter "
                 "AT1', which is independent confirmation the accession is TgAT1."),
    Dataset("curated_resistance_alleles", "Validated resistance-conferring mutations (CURATED)",
            "reference", "literature",
            "Mutations shown to CAUSE drug resistance by putting them back into a clean background",
            ("resistance_allele_count", "resistance_compound_count", "resistance_substitutions",
             "resistance_compounds"),
            "1 gene, 3 substitutions, 3 compounds", pmid="24533298",
            accession="Int J Parasitol Drugs Drug Resist, PMIDs 24533298 and 25941623",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3862444/fullTextXML",
            note="The ONLY curated source in the map: Toxoplasma has no resistome to parse, so this "
                 "is read out of papers one allele at a time and every row carries its paper, "
                 "substitution, compound and how causality was shown. The bar is that the mutation "
                 "was put BACK -- introduced into a clean background and shown to produce the "
                 "resistance -- and that bar is why this holds one gene rather than ten. Nine genes "
                 "were curated from the artemisinin and auranofin in-vitro-evolution studies and "
                 "then DROPPED: the artemisinin table is titled 'Mutations found in candidate "
                 "genes', and the auranofin paper names SOD2 as its likeliest locus and then "
                 "reports that SOD2 L201P was not sufficient to confer resistance when introduced "
                 "into wild-type parasites. Three well-known alleles are deliberately absent and "
                 "the module says why: DHFR-TS pyrimethamine alleles (primary text pre-PMC and "
                 "unreachable), DHODH N302S (primary not open access), and cytochrome b atovaquone "
                 "alleles (mitochondrially encoded, so no row exists in a table of nuclear genes -- "
                 "the eighteen nuclear cytochrome b hits are b-c1 subunits and would be the wrong "
                 "gene). Absence here is ignorance rather than a negative result: nobody selected "
                 "resistance in most genes, so only curated genes carry a value."),
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
    Dataset("plasmodb_pf3d7_expression", "Plasmodium falciparum life-stage and polysomal RNA",
            "transcription", "RNAseq",
            "Transcript abundance across seven life stages, and what is on ribosomes",
            ("expr_ring", "expr_early_trophozoite", "expr_late_trophozoite", "expr_schizont",
             "expr_gametocyte_ii", "expr_gametocyte_v", "expr_ookinete", "expr_asexual_blood",
             "expr_oocyst", "expr_sporozoite", "polysomal_ring", "polysomal_trophozoite",
             "polysomal_schizont", "steady_state_ring", "steady_state_trophozoite",
             "steady_state_schizont", "protein_stage_share_ring",
             "protein_stage_share_trophozoite", "protein_stage_share_schizont",
             "antisense_asexual_blood", "antisense_oocyst", "antisense_sporozoite"),
            "5,720 P. falciparum genes",
            accession="PlasmoDB: Su seven stages, Bunnik polysomal IDC, Gomez-Diaz mosquito stages",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByTaxon/reports/attributesTabular",
            path="starplast/data/pf_nodes.parquet",
            note="Three studies in one report, split across slots so that no column answers two "
                 "questions: the seven-stage study answers the individual stages, and the polysomal "
                 "study is an IDC time course at 0h, 18h and 36h whose steady-state arm answers "
                 "cell-cycle phase while its polysomal arm answers translation. Polysome-associated "
                 "RNA is what is ON ribosomes rather than a transcript level, and keeping its "
                 "steady-state partner is what makes that distinction measurable instead of "
                 "assumed. Stage labels are checked against marker genes rather than trusted, "
                 "because the columns are matched by substring out of one wide report and a "
                 "mislabelling would silently shift a stage: CSP peaks in sporozoite, MSP1 in "
                 "schizont, Pfs16 in gametocyte II. Pfs25 is the informative one -- its TRANSCRIPT "
                 "peaks in gametocyte V rather than in the ookinete where the protein acts, which "
                 "is the textbook translational-repression stockpile and would look like an "
                 "off-by-one error to anyone who checked the protein instead. The antisense columns "
                 "sit beside their sense partners rather than being reduced to a ratio, because the "
                 "denominator is what makes the ratio interpretable. Adding them exposed a silent "
                 "matcher fault worth recording: `sense - asexual blood stages` is a SUBSTRING of "
                 "`antisense - asexual blood stages`, so plain containment made each sense entry "
                 "match two headers, fail its one-match test and vanish -- fetching antisense DELETED "
                 "sense and the slot count fell by two while a slot was being added. The matcher now "
                 "requires the label to start the header or follow a non-alphanumeric character, and "
                 "a test asserts no declared sample goes missing. The TMT proteome in "
                 "the same report is COMPOSITIONAL and is named for it: PlasmoDB serves the "
                 "channels row-normalised, so a gene's three values sum to a constant and the "
                 "columns anti-correlate by construction. They say which stage a protein sits in, "
                 "not how much there is, so `protein abundance · asexual blood stage` is left EMPTY "
                 "rather than filled with a share. The tell would have been well hidden: ring "
                 "protein correlates -0.25 with ring mRNA, which reads as a biological puzzle and "
                 "is only the normalisation showing through."),
    Dataset("pf_alphafold_confidence", "Plasmodium model confidence and disorder (AlphaFold DB)",
            "reference", "structure",
            "Mean pLDDT per protein, and the fraction of it at each confidence band",
            ("mean_plddt", "plddt_fraction_very_low", "plddt_fraction_low",
             "plddt_fraction_confident", "plddt_fraction_very_high", "alphafold_accession"),
            "5,098 of 5,720 genes",
            accession="AlphaFold DB API, per UniProt accession",
            url="https://alphafold.ebi.ac.uk/api/prediction/",
            path="datasets/reference/plasmodb/plasmodb_pf3d7_alphafold.tsv",
            note="Fetched per protein because the bulk proteome archive for this organism is not "
                 "where the documented path says it is. The FRACTIONS matter as much as the mean: a "
                 "protein half well-folded and half disordered has the same mean as one uniformly "
                 "mediocre, and the slot asks about disorder as well as confidence. A gene can "
                 "carry several UniProt accessions -- 876 do, mostly the variant surface families "
                 "where each field isolate's allele has its own entry -- so the fetch tries them in "
                 "order and `alphafold_accession` records which one supplied the model, or a number "
                 "could not be traced back to a structure. Validated on an ordering rather than a "
                 "total: proteins carrying a recognised InterPro domain model at median pLDDT 73.5 "
                 "against 56.0 for those without (p = 2e-159), because a domain is a thing that "
                 "folds. The correlation with protein length is NEGATIVE at -0.555, which is not a "
                 "fault -- it is this proteome's low-complexity asparagine insertions, which are "
                 "long and disordered."),
    Dataset("pf_febrile_stress", "Plasmodium transcription at febrile temperature",
            "transcription", "RNAseq",
            "Wild type and two mutants at 37 C and at the 41 C of a malarial fever",
            ("febrile_wt_37c", "febrile_wt_41c", "febrile_lrr5ko_37c", "febrile_lrr5ko_41c",
             "febrile_dhcko_37c", "febrile_dhcko_41c"), "5,791 genes",
            accession="PlasmoDB Pfal3D7 Febrile temps RNA-Seq",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByTaxon/reports/attributesTabular",
            path="datasets/reference/plasmodb/plasmodb_pf3d7_febrile.tsv",
            note="Shipped as CONDITIONS rather than as a 41-versus-37 contrast, the same restraint as "
                 "the Sir2 entry and for a different reason. There the check on the contrast "
                 "contradicted itself; here it came out NULL. Heat shock proteins move by a median "
                 "log2 of +0.08 against -0.07 for everything else (p = 0.2), so a fever does not "
                 "measurably induce them -- which is consistent with what is known, since this "
                 "organism's chaperones are constitutively high rather than stress-induced, and the "
                 "genes that do rise are Maurer's cleft two-TM proteins and stevor at four to five "
                 "log2, matching published fever-driven surface remodelling. But a null result on "
                 "the one available prediction is not a validation, and a derived column would imply "
                 "it had passed one. A test asserts the arms are on a comparable scale, which is the "
                 "precondition for the caller making the contrast themselves."),
    Dataset("pf_sir2_perturbation", "Plasmodium transcription under Sir2 knockout",
            "transcription", "microarray",
            "Wild type and sir2a / sir2b knockout at ring, trophozoite and schizont",
            ("sir2_wt_ring", "sir2_wt_trophozoite", "sir2_wt_schizont", "sir2a_ko_ring",
             "sir2a_ko_trophozoite", "sir2a_ko_schizont", "sir2b_ko_ring",
             "sir2b_ko_trophozoite", "sir2b_ko_schizont"), "5,615 genes",
            accession="PlasmoDB Sir2 KO Marray",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByTaxon/reports/attributesTabular",
            path="datasets/reference/plasmodb/plasmodb_pf3d7_sir2_perturbation.tsv",
            note="Shipped as stated CONDITIONS and not as a knockout-minus-wild-type contrast, "
                 "which is the quantity anyone will want. The independent check on that contrast "
                 "came out ambiguous and the column would have stated more confidence than there "
                 "is: Sir2a silences subtelomeric var genes, so var should rise in the sir2a "
                 "knockout, and it does in ring (+0.135, p = 4e-12) and schizont (+0.130, "
                 "p = 2e-38) but FALLS in trophozoite (-0.240, p = 3e-22), with effects small "
                 "against a spread of 0.7. That fits the canonical result being subset-specific "
                 "and var probes cross-hybridising across sixty paralogues, but it is not a clean "
                 "confirmation. The conditions themselves are unambiguous -- PlasmoDB names them, "
                 "the values are log intensities, and the medians align across arrays within 0.1, "
                 "which a test asserts because it is the precondition that makes differencing them "
                 "meaningful at all."),
    # `kind` names the measurement underneath, not the fact of derivation -- the same rule the
    # Toxoplasma `stage_enriched` entry records. It is RNA-seq: the argmax of nine expression columns.
    Dataset("pf_iedb_bcell", "Plasmodium antibody epitopes (IEDB)", "reference", "immunity",
            "Distinct antibody epitope sequences per gene", ("n_bcell_epitopes",),
            "434 antigens, 7,366 distinct epitopes", accession="IEDB bcell_search",
            url="https://query-api.iedb.org/bcell_search"
                "?parent_source_antigen_source_org_name=ilike.*Plasmodium%20falciparum*",
            path="datasets/reference/plasmodb/iedb_pf_bcell_epitopes.tsv",
            note="DISTINCT sequences, not assay records. MSP1 alone carries 1,739 epitopes out of "
                 "14,610 records, so counting records would rank antigens by how many groups have "
                 "studied them rather than by how much of the protein antibodies recognise. Reached "
                 "through UniProt rather than through product descriptions -- the Toxoplasma arm has "
                 "to match descriptions because IEDB's Toxoplasma antigen names are verbatim ToxoDB "
                 "text, while the falciparum names carry the accession, which is a better key. The "
                 "209 PlasmoDB accessions naming more than one gene are dropped: an epitope belongs "
                 "to a protein, and attaching it to whichever paralogue sorted first would be "
                 "inventing the answer. 434 of 444 antigens resolve. Absent is absent and not zero, "
                 "because IEDB records what somebody tested. Validated on the history of the field: "
                 "MSP1 is the top antigen and CSP, the RTS,S vaccine antigen, is present."),
    Dataset("pf_iedb_tcell", "Plasmodium T-cell epitopes (IEDB)", "reference", "immunity",
            "Distinct T-cell epitope sequences per gene", ("n_tcell_epitopes",),
            "44 antigens, 1,542 distinct epitopes", accession="IEDB tcell_search",
            url="https://query-api.iedb.org/tcell_search"
                "?parent_source_antigen_source_org_name=ilike.*Plasmodium%20falciparum*",
            path="datasets/reference/plasmodb/iedb_pf_tcell_epitopes.tsv",
            note="Kept apart from the antibody half on purpose, and the numbers show why they are "
                 "not interchangeable: 434 antigens carry an antibody epitope and only 44 carry a "
                 "T-cell one. Pooling them, or filling either slot with the other, would answer one "
                 "question with the other's number. Same UniProt keying and same distinct-sequence "
                 "counting as the antibody table, and the loader reads whichever halves are present "
                 "so one fetch failing does not cost the other column."),
    Dataset("pf_derived_stage_labels", "Plasmodium peak expression and stage label (DERIVED)",
            "transcription", "RNAseq",
            "Maximum expression across stages, and which stage a gene belongs to",
            ("expr_max", "stage_enriched_derived", "stage_margin_derived"),
            "5,720 genes for the maximum, 310 labelled",
            derived_from=("expr_ring", "expr_early_trophozoite", "expr_late_trophozoite",
                          "expr_schizont", "expr_gametocyte_ii", "expr_gametocyte_v",
                          "expr_ookinete", "expr_oocyst", "expr_sporozoite"),
            path="starplast/data/pf_nodes.parquet",
            note="COMPUTED from the stage columns and declaring it, so leakage closure excludes them "
                 "together. The stage call reuses `cellcycle.stage_enrichment` rather than "
                 "reimplementing it, which is deliberate: if the two arms' stage labels are ever "
                 "compared, a difference should mean the biology differs and not that one z-scored "
                 "and the other did not. Only 310 of 5,720 genes are labelled, because a gene is "
                 "left unlabelled unless one stage leads the next by half a z-unit -- a label that "
                 "is really a coin toss looks like a measurement in every table it reaches. Read the "
                 "class counts with the same caveat the Toxoplasma arm carries: ookinete takes 181 "
                 "of the 310 not because it uses more genes but because ring, trophozoite and "
                 "schizont are highly correlated with one another and rarely win by a margin, while "
                 "the mosquito stages are separable. The margin rule is working; the interpretation "
                 "is what needs care."),
    Dataset("plasmodb_identity", "PlasmoDB gene identity", "reference", "identity",
            "Symbols, previous IDs, product descriptions for the Plasmodium arm",
            ("gene_id", "product"), "5,791 P. falciparum 3D7 transcripts", accession="PlasmoDB 3D7",
            url=PLASMODB, path="starplast/data/plasmodb_identity.tsv",
            note="The Plasmodium arm's identity layer, fetched 2026-08-18 through the same WDK "
                 "report as the ToxoDB tables and kept in its own file -- two identifier spaces in "
                 "one index is the merge this project refuses everywhere else. It exists because "
                 "the first Pf source keyed on anything but current accessions joined ZERO rows: "
                 "the 2014 ribosome-profiling deposit reports the pre-2012 chromosome-based ids "
                 "(`PFE0630c`, `PF13_0222`), and a string join found none of its 3,629 genes. 9,106 "
                 "previous ids resolve; 66 are claimed by two current genes each -- a gene model "
                 "SPLIT, seen from the other side -- and those are withdrawn rather than assigned "
                 "to whichever came first, the same rule the Toxoplasma layer applies to 153 "
                 "strings."),
    Dataset("pf_riboseq", "Plasmodium ribosome profiling across the asexual cycle",
            "translation", "riboseq",
            "Ribosome-footprint and mRNA density per gene at five points of the blood-stage cycle",
            ("riboseq_rpf_ring", "riboseq_rpf_early_trophozoite", "riboseq_rpf_late_trophozoite",
             "riboseq_rpf_schizont", "riboseq_rpf_merozoite",
             "riboseq_mrna_ring", "riboseq_mrna_early_trophozoite",
             "riboseq_mrna_late_trophozoite", "riboseq_mrna_schizont", "riboseq_mrna_merozoite"),
            "3,501 genes (61%), 2,182 at the ring and 1,174 at the merozoite", pmid="25493618",
            accession="GSE58402", geo_file_suffix="_rpkm.txt.gz",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58402/suppl/",
            path="datasets/translation/riboseq/25493618/",
            note="The first MEASURED translation on this arm: the slot was previously answerable "
                 "only through polysome-associated RNA, which is what is on ribosomes rather than "
                 "how much ribosome is on it. Both arms of the experiment ship as conditions and "
                 "the ratio between them does NOT, which is the second time this project has "
                 "computed a Plasmodium translation efficiency and refused it -- and this time on a "
                 "different instrument, a different strain and a different decade. Ribosomal "
                 "proteins carry far more footprint than the rest at every stage (median log1p 5.6 "
                 "to 7.9 against 3.5 to 4.1, p <= 2e-18), which is the check that says the "
                 "measurement behaves; but the ratio puts them BELOW the rest at the schizont, and "
                 "correlates negatively with codon adaptation at every stage (rho -0.01 to -0.12), "
                 "where the textbook expectation is positive. Two checks disagreeing is the "
                 "contradictory case, whose answer is to ship the conditions. What the replication "
                 "adds is where to look: `codon_cai_ribosomal` does not separate the very "
                 "ribosomal proteins it is built from in this genome (0.710 against 0.717, "
                 "p = 0.37) while it does in Toxoplasma (0.771 against 0.714, p = 4e-22), so the "
                 "quantity that fails to behave is the codon index, not the footprints. Stage "
                 "labels are the deposit's own and are CHECKED against the independent PlasmoDB "
                 "stage series: ring, early trophozoite, late trophozoite and schizont each "
                 "correlate highest with their own stage (rho 0.65 to 0.77). The merozoite arm has "
                 "no counterpart there and lands on the ring, which is the neighbouring point of "
                 "the cycle rather than a contradiction. Strain W2, not 3D7, so the surface-antigen "
                 "families are the place to distrust it. Keyed on pre-2012 accessions and resolved "
                 "through `plasmodb_identity`; the deposit's `-a`/`-b` split entries are dropped "
                 "rather than summed, since RPKM is already length-normalised."),
    Dataset("pf_berghei_liver_transfer", "P. berghei liver-stage fitness, transferred to falciparum",
            "DNA", "CRISPR_screen",
            "How a berghei knockout fares through the liver, carried onto its falciparum ortholog",
            ("pb_transferred_liver_log2fc", "pb_transferred_liver_reduced"),
            "754 falciparum genes; 180 reduced", pmid="31730853",
            accession="Cell 2019 Table S2 (PlasmoGEM liver stage)",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6904910/supplementaryFiles",
            path="datasets/reference/plasmodb/pb_transfer/31730853/mmc2.xlsx",
            note="Two borrowings, both deliberate. The falciparum id comes from the BLOOD-STAGE "
                 "screen's table -- the same consortium pairing the same mutants -- rather than "
                 "from an orthology this project derived. And the value is the authors' own "
                 "blood-stage-CORRECTED figure for the salivary-gland-to-blood transition, because "
                 "that transition ends in blood and the uncorrected column would call every "
                 "blood-essential gene liver-essential. The sheet has two header rows and repeats "
                 "`Log2-FC / SD / Power` per transition, so columns are read by position and the "
                 "loader refuses the file if the transition is not where it expects it. 507 genes "
                 "are dropped for `no power`: too few barcodes to say anything, which is not a "
                 "measurement of no effect. Validated on the genes the field would name -- LISP1, "
                 "the UIS/ETRAMP early transcribed membrane proteins and perforin-like protein 1 "
                 "all come out reduced, which is the textbook set for liver development and "
                 "hepatocyte egress. The same file's two MOSQUITO transitions are NOT shipped: the "
                 "markers available to check them (P25, P28, SOAP, chitinase) are the redundant "
                 "ones, so nothing in the data confirms the direction, and a transmission slot "
                 "filled on an unchecked axis is what this campaign refuses."),
    Dataset("host_erythrocyte_proteome", "Human red blood cell proteome, by fraction",
            "reference", "proteomics",
            "Which human proteins are present in the cell the blood stage lives in",
            ("rbc_membrane_psms", "rbc_cytoplasm_psms"),
            "5,264 human proteins: 4,777 membrane, 2,350 cytoplasmic", pmid="41654503",
            accession="Sci Data 06792 Supplementary Table S1",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12992542/supplementaryFiles",
            path="datasets/host/erythrocyte/41654503/41597_2026_6792_MOESM2_ESM.xlsx",
            note="The first HOST TISSUE reference in this project, and the thing instruction 39's "
                 "host slots have been waiting for: a proteome of the cell the parasite lives in "
                 "answers a question no pulldown can, since a pulldown says what a bait touched and "
                 "this says what was there to touch. The two fractions are kept apart because they "
                 "are different measurements -- a protein in the membrane extract is at the surface "
                 "the merozoite invades through, one in the cytoplasm is in the haemoglobin around "
                 "it. Self-validating as a fractionation should be: spectrin beta heads the "
                 "membrane list and haemoglobin alpha the cytoplasmic one. Rows are HUMAN proteins "
                 "keyed by UniProt accession and live in the host table, never in a parasite one; a "
                 "row can name several genes (`HBA1; HBA2`) and the string is kept as given rather "
                 "than one of them chosen."),
    Dataset("pf_literature", "Plasmodium falciparum abstract corpus (COMPUTED layer)",
            "reference", "literature",
            "Who is named in the malaria literature, how deeply, and which genes appear together",
            ("n_publications", "n_papers_focal", "n_papers_substantive", "attention_depth",
             "edge:comention"),
            "43,482 abstracts; 732 genes named in the first 10,000", accession="PubMed",
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            path=".claude/skills/plasmodium-scientist/corpus/pubmed_plasmodium.jsonl",
            note="Fetched by `scripts/fetch_pubmed_corpus.py`, which slices the query by YEAR "
                 "because NCBI stops at ten thousand twice over: esearch will not page past it and "
                 "efetch answers 400 for a retstart beyond it. Both limits are silent -- the first "
                 "version of that script fetched 9,989 abstracts of 43,482 and reported success. "
                 "The scan itself is the Toxoplasma arm's, unchanged: `identity` for who is named, "
                 "`corpus` for what a document is, `literature` for the counting and the attention "
                 "correction, so the two arms' attention numbers mean the same thing. What is "
                 "organism-specific is the accession shapes -- this literature cites `PF3D7_`, "
                 "`PFA0110w`, `PF13_0222` and `MAL1P4.01` in the same paragraph -- and the `Pf` "
                 "symbol prefix, both of which are now arguments to `identity.build_index` rather "
                 "than constants in it. Hard-coded to Toxoplasma they registered 9 of 9,106 "
                 "previous accessions and the corpus read as one that never mentions a gene. "
                 "Abstracts only: there is no `incidental` tier, since that means a mention in a "
                 "body or a caption, and no full-text corpus is loaded for this arm."),
    Dataset("pf_berghei_transfer", "P. berghei knockout fitness, transferred to falciparum",
            "DNA", "CRISPR_screen",
            "Relative growth of berghei knockouts, carried onto their falciparum orthologs",
            ("pb_transferred_phenotype", "pb_transferred_growth_rate", "pb_transfer_confidence"),
            "2,448 falciparum genes of 2,578 berghei mutants", pmid="28708996",
            accession="Cell 2017 Table S1 (PlasmoGEM)",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5509546/supplementaryFiles",
            path="datasets/reference/plasmodb/pb_transfer/28708996/mmc1.xlsx",
            note="A TRANSFER, and it says so in every column name, because instruction 39 requires "
                 "one to be visible rather than folded into the measured slot -- transfer berghei "
                 "fitness onto falciparum, hold out falciparum fitness and recover it, and you have "
                 "measured orthology. What makes this one safe is that the orthology is not ours: "
                 "the screen's own table names a falciparum gene per row, and no falciparum gene is "
                 "named by two berghei ones, so nothing is dropped for ambiguity and nothing is "
                 "derived. Forty rows name a transcript rather than a gene and are stripped, the "
                 "same suffix the phosphosite loader handles; kept whole they would have vanished "
                 "for not matching an accession. CHECKED AGAINST THE RECEIVING ARM'S OWN SCREEN, "
                 "which is the check a transfer has to pass: berghei-essential genes have a median "
                 "piggyBac mutagenesis index of 0.160, slow ones 0.394 and dispensable ones 0.996 "
                 "-- monotonic across two species and two unrelated methods, barcoded knockouts in "
                 "mice against saturation mutagenesis in culture, p = 7e-107 -- and 65 of 71 "
                 "ribosomal proteins come out essential. The 12 mutants the screen calls "
                 "`Insufficient data` keep their confidence and lose their phenotype and growth "
                 "rate: that phrase is the absence of a measurement, not a middle value. The other "
                 "3,272 falciparum genes are UNSCREENED, not dispensable, and stay missing."),
    Dataset("pf_idc_timing", "Plasmodium intraerythrocytic cycle timing", "transcription", "RNAseq",
            "When in the 48-hour cycle each transcript peaks, and how strongly it cycles",
            ("idc_peak_hour", "idc_cycling_amplitude"),
            "5,038 genes timed of 5,499; 461 do not cycle strongly enough to place",
            pmid="34668757", accession="GSE163144", geo_file_suffix=".csv.gz",
            url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163144/suppl/",
            # NOT a `derived_from`, though a computation sits inside the loader. That field names
            # shipped COLUMNS this was computed from, so the leakage closure can exclude an input
            # together with its output; these columns are computed from the deposit's own sample
            # files, which are not columns of anything. Declaring the deposit there would point the
            # guard at a string no dataset provides.
            path="datasets/transcription/idc_timecourse/34668757/",
            note="Sixteen timepoints three hours apart across 48 hours, two replicates, read from "
                 "ONE cell of the deposit's square: the 3D7 line in normal (HbAA) red cells. The "
                 "study's variable is sickle-trait haemoglobin and its second line is FUP, so the "
                 "other three cells are perturbations rather than a reference -- and which sample is "
                 "which comes from the deposit's own metadata file rather than from the sample "
                 "titles, which encode the genotype and not the hour. Samples are TPM and sum to "
                 "1e6 exactly, checked rather than assumed. The peak is a PHASE, not the largest "
                 "column: the axis wraps, since hour 48 is hour 0 of the next cycle, and an argmax "
                 "splits the invasion peak between the last timepoint and the first -- it put AMA1 "
                 "and PTRAMP, textbook invasion transcripts, at hour 3, where a culture "
                 "synchronised at invasion is still carrying the merozoite's mRNA. Fitted through "
                 "the first Fourier harmonic (`cellcycle.cyclic_phase`, shared rather than written "
                 "here) the markers land where the biology says: MSP1 45.7 h, SERA5 42.9 h, RhopH2 "
                 "36.2 h, KAHRP 25.1 h, SBP1 15.1 h -- and AMA1 at 2.5 h, which is five hours from "
                 "MSP1 across the wrap and not the other way round. READ IT AS A CIRCLE. Genes "
                 "whose first harmonic explains less than 40% of their variation are left missing, "
                 "because a flat profile still has an angle."),
    Dataset("pf_ip_ms", "Plasmodium co-immunoprecipitation interactome (EPIC)",
            "post_translation", "IP-MS",
            "Which parasite proteins came down with each tagged bait, against its own control",
            ("n_ip_ms_partners",), "98 edges over 65 partners and 3 baits", pmid="28691708",
            accession="Nat Commun 16044 Supplementary Tables 1-5",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5508133/supplementaryFiles",
            path="datasets/reference/plasmodb/ip_ms/28691708/ncomms16044-s1.pdf",
            note="Published as PAGES OF A PDF rather than as a table, and its captions are wrong: "
                 "the table headed `parasite interacting proteins` holds human ones and the table "
                 "headed `human` holds parasite ones. So no caption is read for what a table "
                 "contains -- the bait comes from the column header, which names the pulldown and "
                 "its control, and the organism from the identifier space. The check that licenses "
                 "the assignment is that a BAIT MUST TOP ITS OWN TABLE: PV1 leads the PV1 pulldown "
                 "at 97 and 120 spectra, PV2 leads PV2's, EXP3 leads EXP3's, and a page whose first "
                 "row is not its bait is a continuation rather than a new table. Read two ways "
                 "before it was trusted -- pypdf and poppler independently give 38, 44, 13, 11 and "
                 "37 rows on the five table pages. The count layout is the trap: eight numbers are "
                 "TWO experiments of four (bait, bait, control, control), so splitting them down "
                 "the middle compares experiment 1 with experiment 2 and reports an enriched "
                 "partner with 145 spectra in the untagged line. Caught because that is impossible; "
                 "read correctly, no pair has more spectra in its control than in its bait and the "
                 "median control is 0. Two things are deliberately not carried: the 38 rows of the "
                 "PfEMP1B pulldown, whose bait is a var-gene transgene with no accession in the "
                 "table, and one PV2 row (PIESP2) that carries seven counts instead of eight, since "
                 "the missing number could be either arm. Degree is missing outside the experiment "
                 "-- four pulldowns are not a survey, and a zero would say `nothing binds this` "
                 "about a protein nobody tested."),
    Dataset("pf_cdpk1_dependent_sites", "Plasmodium CDPK1-dependent phosphosites",
            "post_translation", "phosphoproteomics",
            "Phosphosites per gene that are lost when PfCDPK1 is knocked down",
            ("cdpk1_dependent_sites",), "62 genes, 73 sites", pmid="28680058",
            accession="Nat Commun 00053 Supplementary Data 2a",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5498596/supplementaryFiles",
            path="datasets/post_translation/kinase_substrate/28680058/"
                 "41467_2017_53_MOESM3_ESM.xls",
            note="Keyed by SEQUENCE, because nothing else in the file can be resolved: the sites "
                 "are numbered in a 2017 annotation (`3885720(S422)`) that neither current "
                 "accessions nor PlasmoDB's previous-id list carries. The 15-residue window around "
                 "each site is an identifier when it occurs in exactly one protein -- 73 of 79 "
                 "match one gene, 6 match none, none matches two -- and the translation comes from "
                 "`codons.translate` over the CDS table this project already ships, so the "
                 "genetic code is not written down twice. The mapping is then CHECKED against a "
                 "field it did not use: 69 of 73 rows agree with the current product description, "
                 "and the four that do not are re-annotations rather than wrong genes (a "
                 "`conserved membrane protein` now named basal complex protein bleb, a `formin 2` "
                 "now an Eps15-like protein, and two that my word matcher split on a digit). The "
                 "set then reproduces the paper's own conclusion from the other side: 14.5% of the "
                 "62 genes are motor, IMC or invasion machinery against 1.6% of the proteome, and "
                 "GAP45, myosin A, actin I and IMC1c/1g are all in it. Named for DEPENDENCE and "
                 "not for substrate -- a site lost under knockdown may be phosphorylated by this "
                 "kinase or by something downstream of it, and the file cannot tell them apart. "
                 "The paper's PfPKA-R result is not in this sheet and is not claimed here."),
    Dataset("pf_secretome", "Plasmodium extracellular vesicle proteome",
            "post_translation", "proteomics",
            "Parasite proteins found in extracellular vesicles, and how many preparations found them",
            ("ev_studies",), "184 proteins, 53 of them in both preparations", pmid="28944300",
            accession="Wellcome Open Res 11910 S2 (PRIDE PXD006925)",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5583745/supplementaryFiles",
            path="datasets/post_translation/secretome/28944300/"
                 "5c097a1c-efe5-4ed8-b97b-f9ba656268a6.xlsx",
            note="The columns say `extracellular vesicle` and not `secreted`, because the second "
                 "word would assert a route this measurement does not establish. The sheet read is "
                 "the paper's own compilation: the union of two independent EV preparations with a "
                 "membership column each, so 'how many studies saw this' is a fact in the file "
                 "rather than a join. Its other columns are seroreactivity and antibody-array "
                 "results from unrelated studies -- claims about immunity, not about vesicles -- "
                 "and are deliberately not read. The deposit itself is raw-only (24 RAW files, no "
                 "RESULT), which is the usual shape here: per-protein numbers come from the "
                 "supplement. Checked in the direction a secretome should go -- 20.1% carry a "
                 "signal peptide against 10.2% of the proteome (p = 8e-05), exported proteins run "
                 "6.0% against 3.3% (p = 0.06, same direction and not significant at 184 genes) -- "
                 "and RESA, KAHRP, MSP1 and Ag332 are all present. The confound that must travel "
                 "with the column is abundance: EV genes have a median blood-stage expression of "
                 "71.3 against 12.2 for the rest (p = 8e-31), so this is what mass spectrometry "
                 "found in a vesicle preparation and not a list of what the parasite exports. Same "
                 "caveat as hyperLOPIT assignment on the other arm, and for the same reason. "
                 "Absence is unknown and stays missing. A companion boolean completed with "
                 "False would have read as 5,720 genes tested and 5,536 negative, and graded the "
                 "slot A at 100% for an experiment that identified 184 proteins -- so this one "
                 "ships a single column whose presence is the evidence."),
    Dataset("pf_isoforms", "Plasmodium long-read transcript models", "transcription", "nanopore",
            "Transcript models per gene, and how many the annotation does not contain",
            ("n_transcript_models", "novel_transcript_models"),
            "1,857 genes, 2,498 models, 238 novel", pmid="40316999",
            accession="Malar J 05376 Supplementary Data 2",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12046715/supplementaryFiles",
            path="datasets/transcription/isoforms/40316999/SupplementaryData2.xlsx",
            note="SQANTI classifications of long-read models. `full-splice_match` is the reference "
                 "transcript recovered and is NOT counted as novel -- only novel_in_catalog, "
                 "novel_not_in_catalog and fusion are, which is what the Toxoplasma column of the "
                 "same name counts. Absence is sequencing depth rather than a statement that a gene "
                 "has one transcript, so unseen genes stay missing instead of reading as 1. Found "
                 "while looking for something else: this paper was opened for its m6A data, whose "
                 "Pf arm turned out to be a 43-gene intersection with P. vivax rather than a "
                 "methylome, and was refused for that -- the isoform table beside it is the usable "
                 "one. The obvious correlation holds: more expressed genes yield more models "
                 "(rho +0.36), which is detection depth and is why the count is not read as "
                 "isoform diversity."),
    Dataset("pf_lactylome", "Plasmodium lysine lactylome (resolved from NF54)",
            "post_translation", "lactylation",
            "Lactylated lysines per gene, reported against NF54 and resolved to 3D7",
            ("n_lactylsites", "has_lactyl"), "144 genes", pmid="41417877",
            accession="PLoS Genet 1011991 S1",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12742760/supplementaryFiles",
            path="datasets/post_translation/lactylome/41417877/pgen.1011991.s014.xlsx",
            note="Reported against the NF54 annotation, not 3D7, so joining on the accession string "
                 "would have dropped all 186 genes without a word -- the exact failure the "
                 "Toxoplasma identity layer exists to prevent, met here for the first time on this "
                 "arm. `plasmodium.strain_map` resolves it through orthogroups holding exactly one "
                 "gene on each side, which drops the paralogous surface families rather than "
                 "guessing which member a measurement belongs to. Orthology is a claim about "
                 "ancestry, so it is CHECKED against one about identity: 96.8% of the 4,310 pairs "
                 "have exactly the same protein length and 99.2% are within 5%, which is what it "
                 "should look like when one line was cloned from the other, and a test fails if a "
                 "future release breaks it. 144 of 186 genes resolve; the rest are in multi-gene "
                 "groups. Site counts use a 0.75 localisation cut and the flag does not, the same "
                 "split as acetylation."),
    Dataset("pf_acetylome", "Plasmodium lysine acetylome", "post_translation", "acetylation",
            "Acetylated lysines per gene, and whether the gene was seen acetylated at all",
            ("n_acetylsites", "has_acetyl"), "1,145 genes, 2,163 localised sites", pmid="26813983",
            accession="Sci Rep 19722 S2",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4728587/supplementaryFiles",
            path="datasets/post_translation/acetylome/26813983/srep19722-s2.xls",
            note="Two columns built to two standards, because identifying an acetylated PEPTIDE and "
                 "localising the acetyl group to a particular lysine are different claims. The FLAG "
                 "uses every identification; the COUNT uses only sites with an Ascore of 0.75 or "
                 "better, since a site count is meaningless if you do not know which lysine. The "
                 "list is titled Final and is not pre-filtered on localisation -- Ascores run down "
                 "to 0 -- so taking its length as a site count would have been wrong by about a "
                 "quarter. Self-validating: fourteen histones appear, and the most heavily "
                 "acetylated proteins are the PHD finger proteins, the MYST acetyltransferase and "
                 "the coactivator ADA2, which is to say the acetylation machinery itself."),
    Dataset("pf_myristoylome", "Plasmodium N-myristoylome (NMT-inhibitor sensitive)",
            "post_translation", "myristoylation",
            "Proteins whose click-chemistry capture drops when N-myristoyltransferase is blocked",
            ("is_myristoylated",), "16 substrates of 609 assayed", pmid="34695132",
            accession="PLoS Biol 3001408 S11",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8544853/supplementaryFiles",
            path="datasets/post_translation/myristoylome/34695132/pbio.3001408.s011.xlsx",
            note="The evidence for a substrate is not being pulled down -- background comes down "
                 "too -- but coming down LESS when the transferase is inhibited, so the loader "
                 "requires significance AND a negative difference. A positive difference under a "
                 "blocked transferase would be a protein that came down MORE without it, which is "
                 "not what a substrate does. THREE states, not two: 16 substrates, 593 assayed and "
                 "not substrates, and 5,111 genes never in the pulldown, which stay missing -- "
                 "collapsing the last two would claim the whole proteome had been tested for "
                 "myristoylation by one experiment that saw 609 proteins. Sparse because the "
                 "biology is: Plasmodium has roughly thirty predicted NMT substrates. The list "
                 "validates itself -- GAP45, ARO, CDPK1, Rab-5B, ARF1 and ISP3 are the canonical "
                 "N-myristoylated families in apicomplexans and all are present."),
    Dataset("pf_palmitome", "Plasmodium palmitome (observed only)", "post_translation",
            "palmitoylation",
            "Proteins observed S-palmitoylated, with the motif prediction deliberately excluded",
            ("is_palmitoylated",), "503 proteins", pmid="36250062",
            accession="Front Cell Infect Microbiol Table 3",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9556994/supplementaryFiles",
            path="datasets/post_translation/palmitome/36250062/Table_3.xlsx",
            note="Read from the OBSERVED column and not from the sheet named for it. The same "
                 "workbook carries `nrPalmitoylatedProteins`, whose name says palmitoylated and "
                 "whose 3,105 rows are the union of palmitoyl-ABLE -- a motif prediction over 2,902 "
                 "proteins -- and the 503 actually observed. Taking that sheet at its name would "
                 "have called 54% of the proteome palmitoylated, against published palmitomes of "
                 "400 to 500, and the first rows being PfEMP1 and rifin is what gave it away. "
                 "Validated on substrates and on mechanism rather than on a total: GAP45 and CDPK1, "
                 "the canonical Plasmodium substrates, are both present, and membrane proteins are "
                 "enriched 2.1-fold among the palmitoylated (44% against 27%, p = 7e-15), which is "
                 "what a membrane-anchoring modification has to do. ARO is a known miss -- no "
                 "palmitome is complete, and absence here means not observed."),
    Dataset("pf_phosphoproteome_meta", "Plasmodium phosphosites (re-analysis of all public data)",
            "post_translation", "phosphoproteomics",
            "Distinct phosphorylated residues per gene, pooled across every public study",
            ("n_phosphosites", "has_phospho"), "16,318 sites over 2,503 genes",
            accession="PXD046874",
            url="https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/11/PXD046874/",
            path="datasets/reference/plasmodb/phosphosites/",
            note="A re-analysis of every public Plasmodium phosphoproteomics dataset through one "
                 "pipeline, which is what makes a per-gene count meaningful: the same serine found "
                 "by three groups is one site rather than three. Counted as distinct (gene, "
                 "position) pairs and NOT as rows -- a site-centric table still carries one row per "
                 "peptidoform and per source run, so summing rows would count how often a protein "
                 "was looked at instead of how many sites it has, and the files hold millions of "
                 "rows for 16,318 sites. Accessions are stripped of their transcript and product "
                 "suffix (`PF3D7_1346300.1-p1`), or a gene with two products counts its sites "
                 "twice. Where a source has a merged table the merged one is used and its per-run "
                 "siblings are skipped. Missingness mirrors the Toxoplasma arm: the COUNT stays "
                 "missing where nothing was detected, because how many sites a protein has is "
                 "genuinely unknown if mass spectrometry never saw it, while the FLAG is False, "
                 "because whether it was ever observed phosphorylated is a question about the "
                 "evidence and the answer is no. Validated on orderings rather than totals: 67% of "
                 "kinases carry a site against 44% of genes at large, and site count rises with "
                 "protein length at rho +0.46."),
    Dataset("plasmodb_pf3d7_exportpred", "Plasmodium export prediction (ExportPred)",
            "reference", "annotation",
            "Predicted export to the erythrocyte, as an ordinal confidence tier",
            ("export_pred_tier", "is_exported"), "440 genes called at some threshold, 191 at the default",
            accession="PlasmoDB GenesByExportPrediction",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByExportPrediction/reports/attributesTabular",
            path="datasets/reference/plasmodb/exportpred/",
            note="PREDICTED, not measured, which is why it answers `export / PEXEL trafficking` and "
                 "not `exposure to host cytosol` -- that slot wants a measured exportome and a "
                 "sequence model filling it would be a model answering for an experiment. PlasmoDB "
                 "serves ExportPred as a search with a score threshold rather than as a per-gene "
                 "attribute, so the score is recovered by asking at several thresholds and keeping "
                 "the highest a gene survives; the scale saturates, since asking for 20 returns "
                 "nothing, so 10 is the algorithm's own default and the top tier rather than an "
                 "arbitrary cut. It is a TIER and not a boolean because the default loses real "
                 "biology: MESA and PfEMP3 are exported by any textbook and both fall below 10, "
                 "while KAHRP and the FIKK kinases sit above it. Absence is a real negative here "
                 "and not a gap -- a sequence model was evaluated on every protein, so its silence "
                 "is a prediction of not-exported, which is the opposite of the screen columns."),
    Dataset("pf_host_degree", "Plasmodium host interaction degree (COMPUTED)", "reference",
            "crosslink_MS",
            "How many host proteins a gene was crosslinked to, where it was looked at",
            ("n_host_targets",), "117 genes seen, 10 with a host partner",
            derived_from=("bridge:host",),
            path="starplast/data/pf_nodes.parquet",
            note="THREE states, and the middle one is why this is not a fillna(0). The Toxoplasma "
                 "column of the same name is 0 everywhere without a curated host target, which is "
                 "right there because its source is a curated table covering the literature. This "
                 "source is ONE experiment: a gene it never detected has not been shown to lack host "
                 "partners. So the 117 genes seen in the crosslink data carry a count -- zero "
                 "included, because being crosslinked only to parasite proteins is a real observation "
                 "-- and the other 5,603 stay missing. Classification is three-way for a reason a "
                 "test found: a PF3D7 accession the node table does not carry is a PARASITE protein "
                 "with no row, and reading it as host inflated this count. The shipped numbers were "
                 "unaffected, because every accession in this file is in the table, but the fix is "
                 "what stops the next file from being wrong."),
    Dataset("pf_enzyme_classification", "Plasmodium enzyme classification (PlasmoDB)",
            "reference", "annotation",
            "EC number per gene, curated and orthology-derived kept apart",
            ("ec_number", "has_ec", "ec_number_orthology"),
            "1,220 curated, 335 more from orthology",
            accession="PlasmoDB ec_numbers and ec_numbers_derived",
            url="https://plasmodb.org/plasmo/service/record-types/transcript/searches/"
                "GenesByTaxon/reports/attributesTabular",
            path="datasets/reference/plasmodb/plasmodb_pf3d7_ec.tsv",
            note="PlasmoDB serves two EC fields and they are different KINDS of evidence -- one "
                 "curated for this organism, one inferred from the gene's OrthoMCL group -- so they "
                 "are separate columns and `has_ec` counts only the curated one. Merged they would be "
                 "1,584 genes with no way to tell which 335 were never annotated here at all, which "
                 "is inference standing where annotation should. The slot lists the curated column "
                 "first and its policy is `one`, so the leading candidate wins and the derived field "
                 "is there to be chosen deliberately rather than by default."),
    Dataset("pf_codon_usage", "Plasmodium codon usage (COMPUTED)", "reference", "annotation",
            "Effective number of codons, GC3, and CAI against the ribosomal proteins",
            ("codon_enc", "codon_gc3", "codon_cai_ribosomal"), "5,318 genes",
            accession="PlasmoDB-68 Pfalciparum3D7 AnnotatedCDSs",
            url="https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/"
                "PlasmoDB-68_Pfalciparum3D7_AnnotatedCDSs.fasta",
            path="starplast/data/plasmodb_cds.tsv.gz",
            derived_from=("length",),
            note="COMPUTED through the SAME code as the Toxoplasma arm rather than reimplemented -- "
                 "`codons.codon_usage` gained a table and node-table parameter for it. ENC and GC3 "
                 "are definitions and the CAI reference set is 'the ribosomal proteins' in both arms, "
                 "so two implementations could only differ by being wrong in one of them. Sharing it "
                 "buys the first measurement the two arms can be COMPARED on, and the comparison is "
                 "the validation: GC3 median 0.150 here against 0.583 in Toxoplasma, and ENC 37.6 "
                 "against 53.9 -- P. falciparum has the most AT-rich genome of any eukaryote, so "
                 "extreme codon bias is what has to appear, and a test fails if the two arms ever "
                 "converge. The sequence report API returned 422, 400 and 500 to three different "
                 "request shapes; the static release FASTA is what works."),
    Dataset("pf_host_bridge_xlms", "Plasmodium to human contacts (crosslinking MS)", "reference",
            "crosslink_MS",
            "Parasite protein to erythrocyte protein, measured as a crosslink",
            ("bridge:host",), "10 pairs, 10 parasite genes, 7 human proteins", pmid="41966402",
            accession="Cell Rep mmc1 sheet D",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles",
            path="starplast/data/pf_host_bridges.parquet",
            note="The half of the crosslink file the edge layer throws away, and it is a BRIDGE "
                 "rather than an edge for the reason instruction 39 gives: the pair's two ends live "
                 "in different tables and a human protein has no index in this one. Second bridge in "
                 "the project after the Toxoplasma host IP-MS one, and the first thing it needed was "
                 "a species-aware bridge lookup -- both arms key their bridge `host`, because both "
                 "cross to a human protein, so the name cannot say whose contacts these are and only "
                 "the parasite end can. Validated on an interaction that is in the textbooks: MESA "
                 "(PF3D7_0500800) crosslinks to erythrocyte ankyrin, and the rest of the human side "
                 "is stomatin, calpain, actin and spectrin beta -- the membrane skeleton, which is "
                 "what an exported parasite protein should be touching."),
    Dataset("pf_complexes", "Plasmodium complexes from crosslinking MS", "reference",
            "crosslink_MS",
            "Which crosslink-derived complex a gene belongs to, and whether it reaches the host",
            ("complex_id", "complex_size", "complex_spans_host"),
            "128 genes in 42 complexes", pmid="41966402", accession="Cell Rep mmc5 Clusters",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles",
            path="datasets/reference/plasmodb/complexes/41966402/mmc5.xlsx",
            note="`complex_spans_host` is the informative column and is kept rather than dropped. "
                 "Seven of the 47 clusters contain human proteins as well as parasite ones, which is "
                 "not contamination -- the experiment crosslinked parasite inside erythrocyte, so a "
                 "complex reaching into the host is a finding. But a parasite gene in one of those "
                 "has partners this table cannot name, and a reader taking `complex_size` at face "
                 "value would over-count its parasite neighbours. Only parasite members get a row; "
                 "the host members belong to a bridge table. Same study as the crosslink edge layer "
                 "and a different question: that one is which pairs touch, this one is which "
                 "assembly a protein sits in."),
    Dataset("pf_crosslink_ms", "Plasmodium crosslinking MS contacts", "reference", "crosslink_MS",
            "Protein pairs joined by a measured crosslink", ("edge:xlms",),
            "79 parasite-parasite pairs", pmid="41966402",
            accession="Cell Rep mmc1 sheet D",
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles",
            path="datasets/reference/plasmodb/crosslink/41966402/mmc1.xlsx",
            note="Two filters, neither optional. The experiment crosslinked PARASITE INSIDE "
                 "ERYTHROCYTE, so a third of the 106 protein pairs have a human protein at one or "
                 "both ends -- spectrin, band 3, protein 4.2. Those are real contacts and they are a "
                 "HOST BRIDGE rather than a parasite-parasite edge, so they are dropped from this "
                 "layer instead of being indexed against a table with no row for them. And a pair "
                 "whose ends resolve to one gene is a homomeric crosslink: evidence the protein "
                 "self-associates, not an edge, and drawing it would put a zero-length line in the "
                 "graph. Validated on complexes that have to be there: EXP2, PTEX150 and HSP101 -- "
                 "three subunits of the PTEX translocon -- crosslink to one another, prohibitin 1 to "
                 "prohibitin 2, and RAP1 to RAP2. A contact map that missed those would not be "
                 "measuring contacts. Which end is host and which is parasite goes through the "
                 "ACCESSION INDEX and not through a pattern: the first version matched PF3D7_ in the "
                 "mapping field and shipped 73 edges, but the source writes some rows with a UniProt "
                 "symbol instead -- `sp|Q6ZMA7|Pfs16` is PF3D7_0406200, a parasite gene -- so six "
                 "real contacts were dropped as host-at-one-end and a parasite protein was on its "
                 "way into a host bridge. Resolving recovers all 79, and a test fails if the count "
                 "ever drops back."),
    Dataset("pf_relation_layers", "Plasmodium relation layers (COMPUTED)", "reference", "graph",
            "Gene pairs sharing an orthogroup or a domain, and pairs whose stages covary",
            ("edge:orthogroup", "edge:domain", "edge:coexpression"),
            "1,741 + 24,123 + 63,158 pairs",
            derived_from=("orthogroup", "interpro_ids", "expr_ring", "expr_schizont",
                          "expr_sporozoite"),
            path="starplast/data/pf_graph.npz",
            note="A SECOND graph file, because an edge is a pair of indices into a table and a "
                 "falciparum gene has no index in the Toxoplasma one. Constructions are copied from "
                 "the Toxoplasma arm rather than re-invented, so that a difference between the arms "
                 "means the biology differs and not that the edges were drawn by different rules. "
                 "The one deliberate change is the domain weight. Two proteins can share more than "
                 "one InterPro domain, and emitting the pair once per shared domain draws the same "
                 "edge repeatedly, which reads as repeated evidence; the naive construction gave "
                 "10,764 duplicate emissions among 34,887 here, nearly all inside the var, rifin "
                 "and stevor families that share whole multi-domain architectures. The count is now "
                 "the weight, so a pair sharing eleven domains says so. The Toxoplasma arm emits no "
                 "duplicates at all today -- checked rather than assumed -- and would acquire the "
                 "same fault the moment its annotation gained a pair sharing two domains."),
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
    if d.url.startswith((TOXODB, PLASMODB)):
        # A POST with a JSON body, so the URL alone looks like a landing page while the data is
        # entirely fetchable -- fetch_names has done it all along. It now needs an API key: VEuPathDB
        # closed these reports to anonymous use in August 2026, and a guest session is not enough.
        # Saying "yes, with a key you do not have" would be the same overstatement the old
        # "unfetchable" was, one direction over.
        from . import fetch_names
        site = "toxodb" if d.url.startswith(TOXODB) else "plasmodb"
        if not fetch_names.api_key():
            return False, f"needs a VEuPathDB API key ({fetch_names.API_KEY_ENV})"
        return True, site
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
        got = (sources.geo_sample_files(d.accession, dest, d.geo_file_suffix, log=log)
               if d.geo_file_suffix else sources.geo_supplementary(d.accession, dest, log=log))
        return got[0] if got else None

    if how in ("toxodb", "plasmodb"):
        from . import fetch_names
        organism, url = (("Toxoplasma gondii ME49", fetch_names.URL) if how == "toxodb"
                         else ("Plasmodium falciparum 3D7", fetch_names.PLASMODB_URL))
        out = os.path.join(dest, os.path.basename(d.path or f"{key}.tsv"))
        try:
            fetch_names.write(fetch_names.fetch(organism,
                                                ["primary_key", "gene_name", "gene_previous_ids",
                                                 "gene_product"], url), out)
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
