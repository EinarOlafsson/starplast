#!/usr/bin/env python3
"""Small per-gene tables: one file, two columns, one column of the map.

The shape is the point. A ToxoDB search report, a paper's supplementary table and a list extracted
from a spreadsheet all reduce to gene-and-number, and once they do there is no reason for each to
have its own loader. `SOURCES` is the list, and the fourth field of each entry records what was
asked for -- which query, which sheet -- because a search that can be asked five ways produces five
different columns and the note beside the data has to say which one this is.

Named `toxodb_evidence` until it held eLife supplements too, at which point the name was a claim
about provenance that half its rows did not meet.

Every column is the value the source gives, unmodified except for the sign convention below.

## Which way round a fold change points

`fold_change_avg` is **comp/ref**, not ref/comp, and nothing in the API says so. It was established
against a case with only one possible answer: asking the enteroepithelial dataset for tachyzoites as
reference and tissue cysts as comparison puts BAG1 at +21.6 and LDH2 at +16.1 -- both
bradyzoite-specific -- and SRS29B, which is SAG1, at -14.3. Positive is therefore higher in the
comparison group.

That check is worth repeating whenever a new fold-change search is added here. Getting it backwards
does not raise anything; it silently inverts a column.

## The sign convention, again

Like the palmitome, the RNA-seq searches report a signed fold DIFFERENCE: -1257.4 means
1257-fold down, not a negative expression. It is converted with the same `signed_log2`, for the same
reason -- passing it to a logarithm would turn every down-regulated gene into NaN.
"""
from __future__ import annotations

import os

import pandas as pd

from .palmitome import signed_log2

#: Each downloaded report, the column it becomes, and whether its value is a signed fold difference.
#: `query` records what was asked for, because a search that can be asked five ways produces five
#: different columns and the note beside the data must say which one this is.
#: The H4 acetylation ChIP is the counter-example that makes the refusal below firm. Same site, same
#: assay type, same pipeline, same query shape -- and it behaves the way an active mark must:
#: rho +0.36 with expression, +0.43 with promoter ATAC, -0.02 with fitness, and genes in its top
#: decile are expressed four times as highly as those in its bottom. So the H3K4me1 result is not an
#: artefact of how these reports are read here.
#:
#: Also NOT here: the Gregory sense/antisense CHANGE analysis, fetched for the same slot the
#: antisense LEVEL now fills. The distinction is the whole of it -- how much antisense a gene has is
#: a property of the gene and reproduces across independent datasets at a top-500 overlap of 197
#: where chance gives 31; how much its antisense CHANGED over a time course did not reproduce at all.
#: The rejected analysis is described below.
#:
#: The Gregory sense/antisense analysis, fetched for `noncoding and antisense
#: transcription`. Its `max_FC_product` -- the strongest sense-down/antisense-up coupling across a
#: tachyzoite time course -- is independent of expression (rho -0.13), which was encouraging, and
#: does not reproduce. The same analysis on the ME49 and the GT1 time course agrees at rho +0.155
#: over 1,413 shared genes, and the two top-200 lists share 12 genes where chance alone would give
#: 28. A measurement that anti-correlates with its own replicate is measuring the run, not the gene.
#: The report is in `datasets/quarantine/2026_08_16_toxodb/`.
#:
#: Also NOT here: the Ramirez-Flores self-assembled vesicle proteome, fetched for `secretome /
#: excreted`. Taking exosomes and ectosomes against the remaining supernatant, the dense granule
#: proteins come out at -0.85 and the microneme proteins at -3.72 -- DEPLETED from the vesicle
#: fraction -- while ribosomal proteins, which nothing secretes, are the most enriched thing in it at
#: +0.49. The direction was checked both ways round and the two are exact negations, so this is not
#: an orientation mistake. The likeliest reading is that GRAs and MICs are secreted as soluble
#: protein and stay in the supernatant, which would make the dataset a correct measurement of vesicle
#: partitioning and still not an answer to what the parasite secretes. Correct-measurement-of-a-
#: different-thing is not something this can distinguish from wrong, so the slot stays empty.
#:
#: NOT here, and deliberately: the Einstein H3K4me1 ChIP-on-chip. Its 231 "marked" genes have LESS
#: accessible promoters than unmarked genes (-0.51 against +0.39, Mann-Whitney p = 3e-50) and its
#: score correlates with promoter ATAC at rho = -0.37. H3K4me1 marks active and poised chromatin, so
#: that is backwards, and no reading of the assay makes it right. Whether the fault is in the
#: peak-to-gene assignment, in the array's coverage, or in what the dataset actually contains is not
#: something this project can settle -- and an unexplained backwards correlation is exactly the shape
#: of a column that would look like data and rank genes wrongly. The report is in
#: `datasets/quarantine/2026_08_16_toxodb/`, and `chromatin state - histone marks` stays empty.
SOURCES = (
    ("toxodb_epitopes.tsv", "iedb_epitope_count", False,
     "GenesWithEpitopes, organism=Toxoplasma gondii ME49, confidence High+Medium+Low"),
    ("toxodb_h4_acetylation.tsv", "h4_acetylation_chip_score", False,
     "GenesByChIPchip Hakimi/Ali genome-wide H4 K5-K8-K12-K16 acetylation, within 1 kb, no floor"),
    ("toxodb_macrophage.tsv", "macrophage_expression_percentile", False,
     "GenesByRNASeq Saeij 29 strains, ME49-infected murine macrophages, percentile, channel 1"),
    ("antisense_level.tsv", "antisense_expression_percentile", False,
     "GenesByRNASeq full life-cycle, Antisense profileset, percentile across tachyzoite, cyst and "
     "sporulated"),
    ("melting_temperature.tsv", "melting_temperature_tm", False,
     "eLife 80336 supplementary file 3, sheet 4.5_mineCETSA_curve_fits, median Tm where R2 > 0.8"),
    ("crosslink_partners.tsv", "n_crosslink_partners", False,
     "mBio 02159-25 supplementary file s0004, the high-confidence PPI table, partners per protein"),
    ("pvm_proximity.tsv", "pvm_proximity_positive", False,
     "mBio 00260-21 Data Set S1, the authors' likely-positive and likely-negative PVM lists"),
    ("cdpk1_substrates.tsv", "cdpk1_thiophospho_peptides", False,
     "eLife 85654 supplementary file 6, sheet 6.2_ThioP_enriched, master accession per peptide"),
    ("mrna_stability.tsv", "mrna_remaining_5h_actinomycin", False,
     "PLoS Pathogens 1012857 supplementary Table S12, untreated parasites at 5 hours"),
    ("myristoylome.tsv", "myristoylation_confidence", False,
     "eLife 57861 supplementary file 4, Substrate List, confidence High=3 Medium=2 Low=1"),
    ("toxodb_arginine_methylation.tsv", "n_arginine_methylation_sites", False,
     "GenesByPTM monomethylarginine, Yakubu et al. RH proteomics, at least one site"),
    ("toxodb_nanopore_isoforms.tsv", "novel_transcript_models", False,
     "GenesByLongReadEvidence Stuart/Ralph nanopore, ISM + NIC + NNC novelty, >=5 reads"),
    ("toxodb_enteroepithelial.tsv", "ees_vs_tachyzoite_log2", True,
     "GenesByRNASeq Ramakrishnan enteroepithelial, EES1-5 against tachyzoites, sense strand"),
)


def enzyme_classification(base: str, resolve=None, log=print) -> pd.DataFrame:
    """EC number per gene, and whether it has one at all.

    An annotation rather than a measurement, and treated as one -- the same standing as InterPro
    domains, which fill `domain content` on the sequence axis. What makes it worth a metabolism slot
    is that it is the ONLY gene-indexed metabolic datum there is: every other metabolism question in
    the catalog is about metabolites, and a metabolite is not a gene. See `slots.py` on why those
    three slots are `unit=metabolite`.

    `has_ec` is 0 for genes ToxoDB reports without an EC number, not missing: the whole proteome was
    asked, so "no EC assigned" is an answer about the gene.
    """
    path = os.path.join(base, "starplast", "data", "toxodb_ec_numbers.tsv")
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    if d.shape[1] < 2:
        return pd.DataFrame()
    genes = d[d.columns[0]].astype(str)
    if resolve is not None:
        genes = genes.map(lambda g: resolve(g) or g)
    raw = d[d.columns[1]].astype(str).str.strip()
    text = raw.where(~raw.isin(("N/A", "nan", "")), None)
    out = pd.DataFrame({"ec_number": text.to_numpy(),
                        "has_ec": text.notna().astype(float).to_numpy()},
                       index=pd.Index(genes.to_numpy(), name="gene_id"))
    # A gene listed twice keeps the row that HAS an annotation.
    out = out.sort_values("has_ec", ascending=False)
    out = out[~out.index.duplicated()]
    log(f"enzyme classification: {int(out['has_ec'].sum()):,} of {len(out):,} genes carry an EC "
        f"number")
    return out


def evidence(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Every downloaded ToxoDB search report as one column each.

    The value column is taken by POSITION -- the second -- rather than by name. ToxoDB's tabular
    report ships display names as the header (`Epitope Count`, `Score`, `Fold Change (Avg)`), and
    those change with the site's release while the report's shape does not.
    """
    out = pd.DataFrame()
    for filename, column, is_fold, _query in SOURCES:
        path = os.path.join(base, "starplast", "data", filename)
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t")
        if d.shape[1] < 2:
            continue
        genes = d[d.columns[0]].astype(str)
        if resolve is not None:
            genes = genes.map(lambda g: resolve(g) or g)
        raw = d[d.columns[1]]
        values = pd.Series([signed_log2(v) for v in raw] if is_fold
                           else pd.to_numeric(raw, errors="coerce").to_numpy(),
                           index=genes.to_numpy(), dtype=float)
        # Strongest rather than mean, for the same reason the palmitome takes the maximum: two rows
        # for one gene are one gene measured twice, and averaging a hit with a non-hit erases it.
        series = values.groupby(level=0).max()
        out = out.join(series.to_frame(column), how="outer") if len(out) else series.to_frame(column)
        log(f"small table: {column}, {int(series.notna().sum()):,} genes")
    return out
