"""Per-gene tables derived from the deposits added in the 2026-09 data audit.

Each source here is a published deposit -- a supplementary workbook, a GEO processed table -- that
needs a computation between the file and the column: a mean of two replicate screens, a log, a
moderated contrast. The computation lives here and runs against the raw deposit under the dataset
root; `derive_all` writes the result as a small keyed table into the cache
(`starplast/data/deposit_<key>.tsv`), and the merge reads only that. So the map rebuilds without
the raw deposits, and the raw deposits re-derive every shipped number.

`notebooks/derive_deposits_2026_09.ipynb` runs every derivation below with the sanity checks that
admitted it (built by `scripts/derive_deposits.py`).

What each one is, and the check it had to pass:

* **In vivo heart and brain fitness** (Giuliano et al. 2024, PMID 38977907). The four tissues
  already shipped were cited to the 2019 platform paper; they are this study's composite scores,
  identical to the supplement to 5e-8. The same sheet carries heart and brain, the only in vivo
  heart and brain CRISPR fitness that exists.
* **Serum restriction** (Bitew et al. 2025, PMID 41407671): two independent genome-wide screens in
  10% and 1% serum. Fitness in either serum agrees with the fibroblast screen at rho 0.79-0.85 --
  it is fibroblast fitness again -- while the 10%-minus-1% differential is orthogonal to it (rho
  0.06): a new axis, lipid dependence. GRA38, the paper's gene, is -6.66 in 10% and +0.33 in 1%.
* **Translation efficiency, GSE302107**: the most reproducible TE in the organism (replicates rho
  0.973, against 0.87 and 0.67 for the two already shipped); agrees with them at 0.78-0.79.
  Shipped as log2 of the authors' linear footprint/RNA ratio, the scale the other TE columns use.
* **mRNA decay, GSE329845**: wild-type log2(actinomycin D / vehicle) after 4 h from a moderated
  linear model with replicate blocking. No spike-ins, so it is decay RELATIVE to the median
  transcript, not a half-life. Replicates agree at rho 0.96, ribosomal-protein mRNAs are stable.
  Genome-wide (6,406 genes) where the shipped column is the 412-gene unstable tail of another
  study; the two do not correlate (rho -0.03), which the tail's selection explains.
The second pass, 2026-09-26 (instruction 52), added six Plasmodium deposits and the checks that
admitted them:

* **The spatial proteome** (PMID 42218142): the first subcellular localization in the Plasmodium
  table, 1,646 of 3,000 proteins into 24 niches, reproduced exactly, with RAP1, MAHRP1, ACP, EXP2
  and GAPDH each where they belong. The same paper's field pN/pS and between-species dN/dS come with
  it.
* **mRNA synthesis and decay rates** (PMID 29985403): 4-thiouracil labelling every hour of the
  blood-stage cycle. Fluxes in transcripts per minute, NOT half-lives, so not comparable with the
  Toxoplasma actinomycin columns; the paper's claim of transcription in every stage reproduces, and
  the timing agrees with an independent expression series.
* **The blood-stage proteome and Hsp90 dependence** (bioRxiv 10.64898/2026.08.28.747854): abundance
  in a DMSO control and what two chemically distinct Hsp90 inhibitors do to it. The paper's 131
  chaperone-dependent proteins reproduce exactly from its stated rule.
* **RNA dependence** (PMID 38355719): 898 of 3,671 proteins shift towards the light fractions when
  RNA is removed -- the paper's number, from its own q-values.
* **The sexually committed proteome** (PMID 41482054): commitment, a cycle before the gametocyte
  proteome already in the map.
* **The resistome** (PMID 39607932): 724 clones evolved against 118 compounds. Per gene, how many
  selections and compounds mutated it -- the ranking is PfMDR1, PfATP4, PfCRT, PI4K -- plus field
  variation over 5,970 Pf6 isolates.

* **Host**: the human fibroblast response to Toxoplasma (GSE335016, moderated contrast of infected
  against uninfected on log FPKM; CXCL8, IL6 and the interferon genes up, GAPDH flat), baseline mouse
  bone-marrow macrophage expression (GSE267544, M0 arm only -- LPS+IFN-gamma is not infection), and
  the human hepatocyte response to P. falciparum (GSE263643, donor and treatment blocked; weak: 28
  genes at padj < 0.05, and the uninfected wells' mock material is not described).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- statistics
def bh(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    if n:
        order = np.argsort(q)
        ranked = q[order] * n / np.arange(1, n + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        adj = np.empty(n)
        adj[order] = np.minimum(ranked, 1.0)
        out[ok] = adj
    return out


def _trigamma_inverse(x: float) -> float:
    from scipy import special
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        step = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + step
        if -step / y < 1e-8:
            break
    return float(y)


def moderated_t(Y, X, contrast) -> tuple:
    """limma's lmFit + eBayes (Smyth 2004, no intensity trend) for one contrast.

    Y genes x samples on a log scale, X samples x parameters. Returns (coefficient, t, p, info).
    Squeezing each gene's variance toward the common prior is what makes three-replicate contrasts
    usable: a gene whose three values happen to agree is not declared certain.
    """
    from scipy import special, stats
    Y = np.asarray(Y, dtype=float)
    X = np.asarray(X, dtype=float)
    c = np.asarray(contrast, dtype=float)
    XtXi = np.linalg.pinv(X.T @ X)
    B = Y @ X @ XtXi
    resid = Y - B @ X.T
    d = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (resid ** 2).sum(1) / d
    coef = B @ c
    v = float(c @ XtXi @ c)
    ok = s2 > 0
    z = np.log(s2[ok])
    e = z - special.digamma(d / 2) + np.log(d / 2)
    mean_e = e.mean()
    var_e = ((e - mean_e) ** 2).sum() / (len(e) - 1) - special.polygamma(1, d / 2)
    if var_e > 0:
        d0 = 2 * _trigamma_inverse(var_e)
        s02 = float(np.exp(mean_e + special.digamma(d0 / 2) - np.log(d0 / 2)))
        post = (d0 * s02 + d * s2) / (d0 + d)
        df = d + d0
    else:
        d0, s02 = np.inf, float(np.exp(mean_e))
        post = np.full_like(s2, s02)
        df = 1e6
    t = coef / np.sqrt(post * v)
    p = 2 * stats.t.sf(np.abs(t), df)
    return coef, t, p, {"df_residual": int(d), "df_prior": float(d0), "s2_prior": s02}


def _file(root: str, *parts) -> str | None:
    path = os.path.join(root, *parts)
    return path if os.path.exists(path) else None


# --------------------------------------------------------------------------- Toxoplasma
GIULIANO = ("DNA", "CRISPR_screen", "38977907",
            "Supplementary Data 5 (Genome-wide mouse screen scores).xlsx")
GIULIANO_SHEET = "Genome-Wide Differential"
GIULIANO_COLUMNS = {"PE differential composite score": "fit_invivo_PE",
                    "Liver differential composite score": "fit_invivo_liver",
                    "Spleen differential composite score": "fit_invivo_spleen",
                    "Lung differential composite score": "fit_invivo_lung",
                    "Heart differential composite score": "fit_invivo_heart",
                    "Brain differential composite score": "fit_invivo_brain"}
#: The two tissues the ToxoDB route never carried. All six are shipped from the supplement: the
#: four the table already had are identical there (to 5e-8), ToxoDB dropped the tracks in release
#: 71, and resolving the supplement's A/B-split GT1 loci adds 65 genes the ToxoDB route lost.
GIULIANO_NEW = ("fit_invivo_heart", "fit_invivo_brain")


def giuliano_invivo(root: str, all_tissues: bool = True) -> pd.DataFrame:
    """Composite in vivo fitness per tissue, GT1 accessions. Empty if the deposit is absent."""
    path = _file(root, *GIULIANO)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name=GIULIANO_SHEET, header=1)
    d = d.rename(columns={"gene": "gene_id", **GIULIANO_COLUMNS})
    d["gene_id"] = d["gene_id"].astype(str).str.strip()
    # TGGT1_000000 is the library's non-targeting control, not a gene.
    d = d[d["gene_id"].str.match(r"^TGGT1_\d{6}[A-Z]?$") & (d["gene_id"] != "TGGT1_000000")]
    keep = list(GIULIANO_COLUMNS.values()) if all_tissues else list(GIULIANO_NEW)
    out = d[["gene_id"] + keep].copy()
    for c in keep:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.reset_index(drop=True)


SERUM = ("DNA", "CRISPR_screen", "41407671", "41467_2025_66137_MOESM4_ESM.xlsx")
SERUM_SHEET = "10%vs1%FBS_CRISPR-screen_DATA"


def serum_restriction(root: str) -> pd.DataFrame:
    """Fitness in lipid-rich (10% serum) and lipid-limited (1%) medium, and the difference.

    Two independent screens (Exp1 at passage 8; Exp2 at passages 4, 5 and 8), averaged where they
    measured the same passage, because replicate screens of one condition are replicates. Sign as
    published: negative fitness = depleted = the gene is needed; a negative differential means the
    gene is needed in lipid-RICH medium.
    """
    path = _file(root, *SERUM)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name=SERUM_SHEET)
    d.columns = [" ".join(str(c).split()) for c in d.columns]
    num = lambda c: pd.to_numeric(d[c], errors="coerce")          # noqa: E731
    out = pd.DataFrame({"gene_id": d["ToxoDB_ID"].astype(str).str.strip()})
    out["fit_lipid_rich_p8"] = pd.concat([num("10%FBS_P8 Exp1 Fitness score"),
                                          num("10%FBS_P8 Exp2")], axis=1).mean(axis=1)
    out["fit_lipid_limited_p8"] = pd.concat([num("1%FBS_P8 Exp1 Fitness score"),
                                             num("1%FBS_P8 Exp2")], axis=1).mean(axis=1)
    out["fit_lipid_rich_p4p5"] = pd.concat([num("10%FBS_P4 Exp2"), num("10%FBS_P5 Exp2")],
                                           axis=1).mean(axis=1)
    out["fit_lipid_limited_p4p5"] = pd.concat([num("1%FBS_P4 Exp2"), num("1%FBS_P5 Exp2")],
                                              axis=1).mean(axis=1)
    out["fit_serum_differential_p8"] = num("P8 Mean (Exp1, Exp 2) Phenotype (10%-1%)")
    out["fit_serum_differential_p4p5"] = num("P4/P5 Mean Phenotype (10%-1%)")
    out = out[out["gene_id"].str.match(r"^TGGT1_\d{6}[A-Z]?$")]
    return out.reset_index(drop=True)


GSE302107 = ("translation", "riboseq", "GSE302107", "GSE302107_RPKM_and_TE.xlsx")


def riboseq_302107(root: str) -> pd.DataFrame:
    """log2 translation efficiency per replicate (footprint RPKM over total-RNA RPKM).

    A zero ratio is a gene with no footprints counted, which is 'not measured' rather than an
    infinitely low efficiency, so it becomes NaN rather than -inf.
    """
    path = _file(root, *GSE302107)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name="toxo_TE")
    out = pd.DataFrame({"gene_id": d["Gene_ID"].astype(str).str.strip()})
    for rep in ("1", "2"):
        te = pd.to_numeric(d[f"TE_rep{rep}"], errors="coerce")
        out[f"te302107_tachy_r{rep}"] = np.log2(te.where(te > 0))
    return out


GSE329845 = ("transcription", "RNAseq", "GSE329845", "GSE329845_Processed_Data.csv.gz")


def mrna_decay_329845(root: str, return_fit: bool = False):
    """Wild-type mRNA remaining after 4 h of actinomycin D, log2 relative to vehicle.

    Median-of-ratios size factors over all twelve libraries (DESeq2's), genes kept at a mean of 10
    normalized counts in wild-type vehicle, then a moderated linear model on log2(count + 1) with
    replicate as a block. Relative to the median transcript: with no spike-in, a global drop in RNA
    is normalized away, so 0 is 'as stable as the typical mRNA', not 'fully stable'.
    """
    path = _file(root, *GSE329845)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path, encoding="utf-8-sig").set_index("Gene_ID")
    positive = d[(d > 0).all(axis=1)]
    logs = np.log(positive)
    size = np.exp(logs.sub(logs.mean(axis=1), axis=0).median())
    norm = d / size
    vehicle = [f"cWT_vehicle_REP{i}" for i in (1, 2, 3)]
    treated = [f"cWT_ActD_REP{i}" for i in (1, 2, 3)]
    keep = norm[vehicle].mean(axis=1) >= 10
    L = np.log2(norm.loc[keep, vehicle + treated] + 1)
    X = np.column_stack([np.ones(6), [0, 0, 0, 1, 1, 1], [0, 1, 0, 0, 1, 0],
                         [0, 0, 1, 0, 0, 1]]).astype(float)
    coef, t, p, info = moderated_t(L.to_numpy(), X, [0, 1, 0, 0])
    out = pd.DataFrame({"gene_id": L.index.astype(str),
                        "mrna_log2_remaining_4h_actinomycin": coef})
    if return_fit:
        return out, {"size_factors": size.to_dict(), "t": t, "p": p, **info,
                     "normalized": norm}
    return out


# --------------------------------------------------------------------------- host
def _index(root: str, species: str, symbols: bool = False) -> dict:
    from . import host
    ix = host.uniprot_index(root, species, log=lambda *a: None)
    return ix if symbols else ix["by_ensembl"]


def _keyed(frame: pd.DataFrame, ensembl: pd.Series, names: pd.Series, index: dict,
           order: pd.Series) -> pd.DataFrame:
    """Rows keyed by reviewed UniProt accession; of two rows on one accession the better
    expressed is kept (a second Ensembl id on one protein is almost always a minor locus)."""
    acc = ensembl.astype(str).str.split(".").str[0].map(index)
    out = frame.assign(host_id=acc.to_numpy(), host_name=names.astype(str).to_numpy(),
                       _order=order.to_numpy())
    out = out[out["host_id"].notna()].sort_values("_order", ascending=False)
    return out.drop_duplicates("host_id").drop(columns="_order").reset_index(drop=True)


HFF = ("host", "fibroblast_infection", "GSE335016", "GSE335016_245740-HOMO-gene_fpkm.txt.gz")


def hff_tg_infection(root: str) -> pd.DataFrame:
    """Human foreskin fibroblasts infected with wild-type Toxoplasma against uninfected, log2.

    Moderated t on log2(FPKM + 1) over the three arms (uninfected, wild-type, TGGT1_245740
    knockout; three replicates each), genes with a mean FPKM of 1. Only the wild-type contrast is
    shipped: the knockout answers a question about one parasite gene, not about infection.
    """
    path = _file(root, *HFF)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    ctrl = [c for c in d.columns if "Control" in c]
    wt = [c for c in d.columns if "WTinf" in c]
    ko = [c for c in d.columns if "KO245740" in c]
    if not (len(ctrl) == len(wt) == len(ko) == 3):
        return pd.DataFrame()
    expressed = d[ctrl + wt + ko].mean(axis=1) >= 1
    d = d[expressed].reset_index(drop=True)
    L = np.log2(d[ctrl + wt + ko] + 1)
    X = np.column_stack([np.ones(9), [0] * 3 + [1] * 3 + [0] * 3,
                         [0] * 6 + [1] * 3]).astype(float)
    coef, _t, p, _info = moderated_t(L.to_numpy(), X, [0, 1, 0])
    frame = pd.DataFrame({"hff_tg_infection_log2fc": coef, "hff_tg_infection_padj": bh(p)})
    index = _index(root, "human")
    if not index:
        return pd.DataFrame()
    return _keyed(frame, d["gene_id"], d["gene_name"], index, d[ctrl + wt].mean(axis=1))


BMDM = ("host", "bmdm", "GSE267544", "GSE267544_fpkm_BMDM.txt.gz")


def bmdm_baseline(root: str) -> pd.DataFrame:
    """Unstimulated (M0) mouse bone-marrow macrophage expression, TPM from FPKM, 3 replicates.

    Keyed by Ensembl, never by symbol: the deposit's symbols are old (Emr1 for Adgre1, Irg1 for
    Acod1) and a symbol join would silently lose exactly the macrophage markers.
    """
    path = _file(root, *BMDM)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    m0 = [c for c in d.columns if c.startswith("M0_")]
    if not m0:
        return pd.DataFrame()
    tpm = d[m0] / d[m0].sum() * 1e6
    frame = pd.DataFrame({"bmdm_tpm": tpm.mean(axis=1)})
    index = _index(root, "mouse")
    if not index:
        return pd.DataFrame()
    return _keyed(frame, d["geneID"], d["Gene_name"], index, frame["bmdm_tpm"])


HEPATOCYTE = ("host", "hepatocyte_infection", "GSE263643",
              "GSE263643_deseq2-normalized_counts.txt.gz")


def hepatocyte_pf_infection(root: str) -> pd.DataFrame:
    """Primary human hepatocytes infected with P. falciparum against uninfected, log2.

    The paper never makes this comparison; it is made here from its DESeq2-normalized counts with
    treatment and donor as blocks (freshly isolated cells excluded), genes with a mean of 10. The
    infected wells hold few infected cells, so the contrast is diluted and weak by construction.
    """
    path = _file(root, *HEPATOCYTE)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t").set_index("Gene")
    cols = [c for c in d.columns if not c.startswith("fresh")]
    infected = np.array([("infected" in c) and ("uninfected" not in c) for c in cols], float)
    treatment = [c.split("_")[0] for c in cols]
    replicate = [c.split("_")[-1] for c in cols]
    X = np.column_stack([np.ones(len(cols)), [t == "iwp2" for t in treatment],
                         [t == "5c" for t in treatment], [r == "rep2" for r in replicate],
                         infected]).astype(float)
    expressed = d[cols].mean(axis=1) >= 10
    L = np.log2(d.loc[expressed, cols] + 1)
    coef, _t, p, _info = moderated_t(L.to_numpy(), X, [0, 0, 0, 0, 1])
    frame = pd.DataFrame({"hepatocyte_pf_infection_log2fc": coef,
                          "hepatocyte_pf_infection_padj": bh(p)})
    full = _index(root, "human", symbols=True)
    index = full["by_ensembl"]
    if not index:
        return pd.DataFrame()
    ids = pd.Series(L.index)
    # The deposit carries Ensembl ids only; the name shown is the reviewed entry's gene symbol.
    symbol_of = {acc: sym for sym, acc in full["by_symbol"].items()}
    names = ids.str.split(".").str[0].map(index).map(symbol_of).fillna(ids)
    return _keyed(frame, ids, names, index, pd.Series(d.loc[expressed, cols].mean(axis=1).values))


# --------------------------------------------------------------------------- the second wave
GLUCOSE = ("DNA", "CRISPR_screen", "glucose_limitation_2025", "TableS1.xlsx")


def glucose_limitation(root: str) -> pd.DataFrame:
    """Fitness with glucose or with glutamine withdrawn (Uboldi et al., bioRxiv 2025, Table S1).

    After three passages in complete medium the library was split into glucose-only and
    glutamine-only medium for two more, three replicates each. Each arm is measured against the
    post-selection library, so negative means depleted there. The dependence column is the authors'
    edgeR contrast, glutamine-only minus glucose-only: negative means the gene is needed when
    glucose is absent. GDH1 ranks first of 8,155 and PEPCK second, which is the paper's result, and
    ribosomal proteins are needed in complete medium (median -2.21 against -0.02).

    Noisy: the three replicate differentials agree at rho 0.10-0.17 only, so the FDR ships beside
    the value and neither should be read gene by gene without it.
    """
    path = _file(root, *GLUCOSE)
    if path is None:
        return pd.DataFrame()
    x = pd.read_excel(path, sheet_name="Fig 1 Gln vs Glc")
    post = [f"LogCPM(Post Selection (Sample 2) Rep {i})" for i in (1, 2, 3)]
    glc = [f"LogCPM(Glucose (Sample 4) Rep {i})" for i in (1, 2, 3)]
    gln = [f"LogCPM(Glutamine (Sample 3) Rep {i})" for i in (1, 2, 3)]
    out = pd.DataFrame({"gene_id": x["Accession Number"].astype(str).str.strip()})
    out["fit_complete_medium_2025"] = pd.to_numeric(
        x["logFC (Post Selection (2) vs Input (2))"], errors="coerce")
    out["fit_no_glutamine"] = x[glc].mean(axis=1) - x[post].mean(axis=1)
    out["fit_no_glucose"] = x[gln].mean(axis=1) - x[post].mean(axis=1)
    out["fit_glucose_dependence"] = pd.to_numeric(
        x["logFC (Glutamine (3) vs Glucose (4))"], errors="coerce")
    out["fit_glucose_dependence_fdr"] = pd.to_numeric(
        x["FDR (Glutamine (3) vs Glucose (4))"], errors="coerce")
    return out[out["gene_id"].str.match(r"^TG[A-Z0-9]+_\d{5,6}[A-Z]?$")].reset_index(drop=True)


UTR5 = ("translation", "riboseq", "GSE302108",
        "SupplementaryData4_TE_info_and_UTRfeatures.xlsx")


def utr5_architecture(root: str) -> pd.DataFrame:
    """Each gene's 5' untranslated region, reannotated from long reads (Peters et al., bioRxiv
    2025, Supplementary Data 4): length, upstream AUGs and ORFs, overlapping ORFs, in-frame
    extensions, and how closely the start context matches a Kozak sequence.

    Sequence properties rather than a measurement of the gene, and they behave as the biology says
    they should: upstream AUGs go with LOW translation efficiency (rho -0.47 against the mean of
    the TE columns) and Kozak strength with high (+0.22). The study's reporter assay is not a
    per-gene column -- its 30,235 scored sequences are variants of twelve endogenous UTRs.
    """
    path = _file(root, *UTR5)
    if path is None:
        return pd.DataFrame()
    f = pd.read_excel(path, sheet_name="Toxoplasma")
    out = pd.DataFrame({"gene_id": f["gene_ID"].astype(str).str.strip(),
                        "utr5_length": f["fiveprimeUTR_length"],
                        "utr5_n_uaugs": f["n_uAUGs"], "utr5_n_uorfs": f["uORFs"],
                        "utr5_n_oorfs": f["oORF"],
                        "utr5_n_inframe_ext": f["in_frame_extension"],
                        "utr5_kozak_score": f["Kozak_similarity_score"]})
    for c in out.columns[1:]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out[out["gene_id"].str.startswith("TGME49_")].reset_index(drop=True)


BRADY = ("transcription", "scRNAseq", "41580398", "41467_2026_68489_MOESM3_ESM.xlsx")


def bradyzoite_subtypes(root: str) -> pd.DataFrame:
    """Average expression in each of five in vivo bradyzoite subtypes (Ulu et al., Nat Commun
    2026, Supplementary Data 1): single cells from brain cysts of mice infected for 30 days.

    The subtypes are the point -- a tissue cyst is not one transcriptional state. Every group has
    BAG1, LDH2, ENO1, SRS9 and CST1 above its 92nd percentile and SAG1 below the 41st, so all five
    are bradyzoites; SRS22A is 3.27 in Group B against 0.34-0.66 elsewhere, which is the paper's
    Group B signature.
    """
    path = _file(root, *BRADY)
    if path is None:
        return pd.DataFrame()
    raw = pd.read_excel(path, header=None, skiprows=1)
    out = pd.DataFrame({"gene_id": raw[0].astype(str).str.strip()})
    for i, group in enumerate("ABCDE"):
        base = 2 + 4 * i
        if not (raw[base].astype(str) == f"Group {group}").all():
            raise ValueError(f"Supplementary Data 1 layout changed at Group {group}")
        out[f"bzsub_{group}_expr"] = pd.to_numeric(raw[base + 1], errors="coerce")
    return out[out["gene_id"].str.startswith("TGME49_")].reset_index(drop=True)


IRON = ("translation", "proteomics", "41925342")


def iron_depletion(root: str) -> pd.DataFrame:
    """What 24 hours without iron does to each protein, and to each transcript (Hanna et al., mBio
    2026, Tables S1 and S2). Positive means higher when iron is withheld.

    The proteome is genome-scale: 5,052 protein groups, three replicates each side. The RNA
    contrast covers only the 3,113 genes of the paper's joint analysis, which is a
    significance-filtered subset, so a missing RNA value means 'not reported' rather than
    'unchanged'. The paper's 61 iron-sulfur proteins shift down on balance (median -0.04 against
    +0.02; 12 of the 16 that change significantly fall), a modest effect. Despite the paper's title there
    is no ribosome profiling in it -- translation is measured by microscopy -- so this fills protein
    abundance under stress, not translation.
    """
    folder = _file(root, *IRON)
    if folder is None:
        return pd.DataFrame()
    p = pd.read_excel(os.path.join(folder, "mbio.03788-25-s0002.xlsx"), sheet_name=0, header=1)
    d = pd.read_excel(os.path.join(folder, "mbio.03788-25-s0003.xlsx"), sheet_name=0, header=1)
    prot = pd.DataFrame({"gene_id": p["Protein_Accessions"].astype(str).str.strip(),
                         "iron_depletion_protein_log2fc": pd.to_numeric(p["logFC"],
                                                                        errors="coerce"),
                         "iron_depletion_protein_padj": pd.to_numeric(p["adjusted_p"],
                                                                      errors="coerce")})
    prot = prot[prot["gene_id"].str.match(r"^TG[A-Z0-9]+_\d{5,6}[A-Z]?$")]
    prot = prot.groupby("gene_id").agg({"iron_depletion_protein_log2fc": "mean",
                                        "iron_depletion_protein_padj": "min"})
    # The two halves are keyed differently on purpose: the proteome gives only GT1 accessions, and
    # the transcripts are keyed on the ME49 accessions the paper pairs them with, because keyed on
    # GT1 ten of them no longer resolve. They meet after the identity layer resolves both, so no row
    # of THIS table holds a protein beside its own transcript; the notebook pairs them to check it.
    rna = pd.DataFrame({"gene_id": d["ME49"].astype(str).str.strip(),
                        "iron_depletion_rna_log2fc": pd.to_numeric(d["RNA_log2FoldChange"],
                                                                   errors="coerce")})
    rna = rna[rna["gene_id"].str.startswith("TGME49_")].groupby("gene_id").mean()
    return prot.join(rna, how="outer").reset_index()


TURBOID = ("post_translation", "BioID", "organelle_surface_turboid_2026", "media-1.xlsx")
#: bait -> (sheet of everything it detected, sheet of the hits that also beat the cytosolic bait)
TURBOID_BAITS = {"apicoplast": ("Apicoplast - and + biotin", "Apicoplast vs spatial reference"),
                 "mitochondrion": ("Mitochodrion - and + biotin", "Mitoch vs spatial reference"),
                 "er": ("ER - and + biotin", "ER vs spatial reference")}


def _product_index() -> dict:
    """product text -> gene ids, from the shipped table, to repair a broken accession."""
    path = os.path.join(os.path.dirname(__file__), "data", "nodes.parquet")
    try:
        n = pd.read_parquet(path, columns=["gene_id", "product"])
    except (OSError, ValueError, KeyError):
        return {}
    return n.groupby("product")["gene_id"].apply(list).to_dict()


def _accession(value, product, products: dict):
    """An accession, or the gene its product text names when the accession is malformed.

    One row of this deposit carries a backtick where its accession should be. Repaired only when
    the product matches exactly one gene, because a guess here attaches a measurement to the wrong
    protein and nothing downstream could tell.
    """
    import re
    acc = str(value).strip()
    if re.fullmatch(r"TG[A-Z0-9]+_\d{5,6}[A-Z]?", acc):
        return acc
    hits = products.get(str(product).replace("product=", "").strip(), [])
    return hits[0] if len(hits) == 1 else None


def organelle_surface(root: str) -> pd.DataFrame:
    """What sits against the cytosolic face of the apicoplast, mitochondrion and ER (Parker & Huet,
    bioRxiv 2026, Table S1): TurboID on baits anchored in each outer membrane, tail out.

    Per bait, the enrichment of each protein with biotin over without, and whether it is one of the
    paper's stringent hits -- also enriched over a cytosolic TurboID, so not merely abundant near
    everything. Against hyperLOPIT the stringent mitochondrial hits are mitochondrial at odds 13.8
    and the ER hits ER at 7.6; the apicoplast bait shows no enrichment for apicoplast proteins
    (odds 1.2, p = 0.8), so that bait reports proximity, not location.

    A protein a bait never detected keeps a missing stringent flag rather than a zero: that bait
    did not test it.
    """
    path = _file(root, *TURBOID)
    if path is None:
        return pd.DataFrame()
    sheets = {k.strip(): v.rename(columns=lambda c: str(c).strip())
              for k, v in pd.read_excel(path, sheet_name=None).items()}
    products = _product_index()

    def ids(frame):
        blank = [""] * len(frame)
        return [_accession(a, p, products) for a, p in
                zip(frame["Accession.Number"], frame.get("Product_Description", blank))]

    frames = []
    for bait, (detected, stringent) in TURBOID_BAITS.items():
        a = sheets[detected]
        t = pd.DataFrame({"gene_id": ids(a),
                          f"surface_{bait}_log2fc": pd.to_numeric(a["log2foldchange"],
                                                                  errors="coerce")})
        t = t.dropna(subset=["gene_id"]).groupby("gene_id").mean()
        hits = {g for g in ids(sheets[stringent]) if g}
        t[f"surface_{bait}_stringent"] = t.index.isin(hits).astype(float)
        frames.append(t)
    out = pd.concat(frames, axis=1)
    for bait in TURBOID_BAITS:
        missing = out[f"surface_{bait}_log2fc"].isna()
        out.loc[missing, f"surface_{bait}_stringent"] = np.nan
    return out.reset_index()


K562 = ("host", "k562_rhoptry_screen", "v2_media-1.xlsx")


def k562_rhoptry_screen(root: str) -> pd.DataFrame:
    """Which HUMAN genes a Toxoplasma rhoptry needs in order to discharge into the cell it invades
    (Valleau et al., bioRxiv 2025, Table S1): a genome-wide CRISPR screen in K562 read out by a
    fluorescent reporter of discharge.

    Cells were sorted into those still susceptible to discharge and those that had become
    resistant; each was compared with the unsorted population. A POSITIVE combined score means
    knocking the gene out makes discharge fail, so the host gene is required for it. The screen
    recovers its own biology: SLC35A2, the transporter the paper is about, ranks first of 20,010,
    and a panel of twenty N-glycan genes has a median rank of 135 with 85% in the top tenth, while
    B4GALT1, FUT8 and ST6GAL1 -- glycosylation the paper argues is NOT involved -- are not hits.

    A host table, so rows are keyed by reviewed UniProt accession. The deposit gives gene symbols;
    a symbol that resolves to no reviewed human protein (155 of them are microRNA loci) is dropped.
    """
    path = _file(root, *K562)
    if path is None:
        return pd.DataFrame()
    t = pd.read_excel(path, sheet_name="Rhoptry discharge MAgECK-MLE")
    symbols = t["Gene"].astype(str).str.strip()
    frame = pd.DataFrame({
        "rhoptry_discharge_score": pd.to_numeric(t["Combined Rhoptry Score"], errors="coerce"),
        "rhoptry_discharge_beta": pd.to_numeric(t["greenN|beta"], errors="coerce"),
        # The Wald FDR, not the permutation one: the permutation FDR of this deposit is quantised
        # into a few values, which would read as ties between genes that are not tied.
        "rhoptry_discharge_fdr": pd.to_numeric(t["greenN|wald-fdr"], errors="coerce")})
    index = _index(root, "human", symbols=True)
    by_symbol = index.get("by_symbol") or {}
    if not by_symbol:
        return pd.DataFrame()
    accession = symbols.map(by_symbol)
    out = frame.assign(host_id=accession.to_numpy(), host_name=symbols.to_numpy())
    out = out[out["host_id"].notna()]
    # One protein reached by two symbols keeps the stronger score rather than their mean: the two
    # rows are two knockouts of one protein, and averaging a hit with a non-hit erases it.
    return out.reindex(out["rhoptry_discharge_score"].abs().sort_values(
        ascending=False).index).drop_duplicates("host_id").reset_index(drop=True)


MELTOME = ("translation", "proteomics", "41315737", "pxd056075_tm_gene.tsv")


def pf_melting_temperature(root: str) -> pd.DataFrame:
    """The temperature at which each P. falciparum protein comes out of solution (Pazicky et al.,
    Nat Microbiol 2025, Supplementary Table 1; PRIDE PXD056075).

    Intact-cell CETSA over ten temperatures, seven points of the blood-stage cycle, three
    replicates. No melting temperature is published anywhere -- the deposit is raw mass
    spectrometry and the paper tabulates abundances -- so it is FITTED, by
    `scripts/fit_meltome.py`, with the acceptance the shipped Toxoplasma column already uses
    (R-squared at least 0.8, melting point between 30 and 80 degrees). That script writes the table
    this reads; the fit is not repeated inside the build.

    It behaves as a melting point must: replicates agree at rho 0.84 with a spread of 1.15 degrees,
    subunits of one complex melt together (0.19 degrees across the T-complex, 0.42 across the
    proteasome's alpha ring, against 5.6 for size-matched random sets, p = 9e-5), glycolytic
    enzymes and the proteasome are thermostable at about 62.5 degrees while HSP70 and HSP90 are
    labile at 50.4. It also agrees with the independent Toxoplasma measurement through orthology
    (+0.28). Median 53.0 degrees against 55.2 for Toxoplasma.

    One column, not seven: a protein's melting point barely moves across the cycle (3.7 degrees of
    range against 1.15 of replicate noise), so per-stage columns would claim a resolution the data
    does not have.
    """
    path = _file(root, *MELTOME)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    return pd.DataFrame({"gene_id": d["gene_id"].astype(str).str.strip(),
                         "melting_temperature_tm": pd.to_numeric(d["tm_median"], errors="coerce"),
                         "melting_temperature_sd": pd.to_numeric(d["tm_sd"], errors="coerce")})


FERTILITY = ("DNA", "transmission_screen", "39541984", "mmc2.xlsx")


def pb_fertility(root: str) -> pd.DataFrame:
    """Whether a gene is needed for male or for female fertility, carried from Plasmodium berghei
    (Russell et al., Cell Syst 2024, Table S1 sheet A): barcoded knockouts crossed and scored
    through mosquito transmission.

    Measured in the rodent parasite and transferred by orthology -- the table names the falciparum
    ortholog itself, so nothing is guessed here. Negative means the mutant lost fertility on that
    side. The result is sex-specific and the controls prove it: HAP2, P48/45, P230, CDPK4 and MAPK2
    fail in males only, P47, NEK4 and DMC1 in females only, and the redundant P25/P28 pair is
    correctly not called. The union of reduced mutants is 348, the paper's own number.

    A transfer between species, so it belongs to the orthology axis for leakage: an ortholog's
    measurement is not an independent measurement of this gene.
    """
    path = _file(root, *FERTILITY)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name="(A) Fertility screen results")
    columns = {c: str(c) for c in d.columns}

    def column(prefix):
        for c, name in columns.items():
            if name.startswith(prefix):
                return d[c]
        raise KeyError(prefix)
    ortholog = column("P. falciparum 3D7 orthologue").astype(str).str.strip()
    out = pd.DataFrame({"gene_id": ortholog,
                        "fertility_female": pd.to_numeric(column("Female fertility [log2FC] ("),
                                                          errors="coerce"),
                        "fertility_male": pd.to_numeric(column("Male fertility [log2FC] ("),
                                                        errors="coerce")})
    out = out[out["gene_id"].str.match(r"^PF3D7_\d{7}$")]
    # A falciparum gene named by two berghei mutants, or a mutant covering two genes, is dropped
    # rather than averaged: which mutant measured which gene is then unknown.
    counts = out["gene_id"].value_counts()
    return out[out["gene_id"].isin(counts[counts == 1].index)].reset_index(drop=True)


# --------------------------------------------------------------------------- Plasmodium, third wave
PF_LATENCY = ("transcription", "scRNAseq", "pf_latency_2026")


def pf_latency(root: str) -> pd.DataFrame:
    """What a parasite entering a drug-tolerant latent state transcribes (bioRxiv 2026,
    Supplementary Tables 7 and 24): single cells after artemisinin or nutrient deprivation.

    Positive means higher in the latent cells. The paper's own 200-gene latency classifier ships
    beside it as a membership flag, so a user can ask 'is my gene one of the ones that defines this
    state' without rebuilding the classifier.

    Its FDR is deliberately NOT shipped: cells were treated as replicates, so 87% of genes reach
    padj 0.05 and the column would say 'significant' about almost everything. Absence is not
    evidence either -- ribosomal-protein genes and most var and rifin genes are missing from the
    table altogether, so a gene without a value was not necessarily unchanged.
    """
    folder = _file(root, *PF_LATENCY)
    if folder is None:
        return pd.DataFrame()
    de = pd.read_excel(os.path.join(folder, "ST07_media-4.xlsx"), header=1)
    out = pd.DataFrame({"gene_id": de["Gene"].astype(str).str.strip(),
                        "latency_log2fc": pd.to_numeric(de["log2FC"], errors="coerce")})
    classifier = pd.read_excel(os.path.join(folder, "ST24_media-20.xlsx"), header=4)
    members = classifier["PF3D7 accession"].astype(str).str.strip()
    members = set(members[members.str.startswith("PF3D7")])
    out["latency_classifier_member"] = out["gene_id"].isin(members).astype(float)
    out = out[out["gene_id"].str.match(r"^PF3D7_\w+$")]
    return out.groupby("gene_id", as_index=False).mean()


PF_M6A = ("transcription", "m6a", "pf_m6a_2026", "SupplementaryTables_media-2.xlsx")
#: Sheet per time point of the intra-erythrocytic cycle, and the read depth a transcript needs
#: before its stoichiometry is believed.
PF_M6A_SHEETS = {28: "Supplementary Table 6", 32: "Supplementary Table 7",
                 36: "Supplementary Table 8"}
PF_M6A_MIN_DEPTH = 20


def pf_m6a(root: str) -> pd.DataFrame:
    """How many methylated adenosines a transcript carries, and how fully they are methylated
    (bioRxiv 2026, Supplementary Tables 6-8): nanopore direct RNA sequencing at three time points.

    Sites are the union over the time points, because a site detected at one hour is a site. The
    stoichiometry is the weighted average methylation of the canonical sites, averaged over the time
    points where the transcript was read at least 20 times -- below that the fraction is a
    coin-flip. Where a gene has several transcripts, the deepest one represents it.

    Only the control condition is used. The deposit's two groups are unlabelled and the second has
    about half the methylation of the first, which is consistent with the paper's knock-sideways of
    the methyltransferase; that inference is good enough to pick the control arm but not good enough
    to ship a difference column, so no difference is shipped.
    """
    import ast
    path = _file(root, *PF_M6A)
    if path is None:
        return pd.DataFrame()
    frames = []
    for hour, sheet in PF_M6A_SHEETS.items():
        d = pd.read_excel(path, sheet_name=sheet, header=2)
        d["gene_id"] = d["ID"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        d["sites"] = d["canonical_mods"].astype(str).map(
            lambda v: set(ast.literal_eval(v)) if v.startswith("[") else set())
        d["hour"] = hour
        frames.append(d)
    a = pd.concat(frames, ignore_index=True)
    sites = a.groupby("gene_id")["sites"].agg(lambda s: set().union(*s))
    out = pd.DataFrame({"m6a_n_canonical_sites": sites.map(len).astype(float)})
    deepest = a.sort_values("average_depth_g1", ascending=False).drop_duplicates(
        ["gene_id", "hour"])
    read = deepest[deepest["average_depth_g1"] >= PF_M6A_MIN_DEPTH]
    methylated = read[read["sites"].map(len) > 0]
    out["m6a_canonical_stoichiometry"] = methylated.groupby("gene_id")["canonical_wam_g1"].mean()
    out.index.name = "gene_id"
    return out.reset_index()


PF_GAMETOCYTE = ("translation", "proteomics", "pf_gametocyte_2026", "Extended Datasets.xlsx")


def pf_gametocyte_proteome(root: str) -> pd.DataFrame:
    """What a mature (stage V) gametocyte's proteome contains, and what it is still making
    (bioRxiv 2026, Extended Data Tables 1 and 3): label-free quantification, and click chemistry on
    newly made protein.

    The abundance column is clean -- replicates agree at about 0.99. Newly made means found with the
    label and in no control: the deposit's group 4, which is the paper's 705 proteins. The
    translatome is shipped only as membership, not as an enrichment: there is one pooled sample per condition and most of the
    proteins it names are absent from the controls, so a ratio would be a ratio to nothing. A gene
    the proteome quantified but the translatome did not name is a zero rather than a blank, because
    that experiment did look for it.
    """
    import re
    path = _file(root, *PF_GAMETOCYTE)
    if path is None:
        return pd.DataFrame()
    valid = re.compile(r"^PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})$")
    total = pd.read_excel(path, sheet_name="Extended Data Table 1", header=2)
    total = total[total["Gene ID"].astype(str).str.match(valid)]
    replicates = [c for c in total.columns if str(c).strip().startswith("Replicate")]
    if len(replicates) != 3:
        raise ValueError(f"expected three replicate columns, found {replicates}")
    mean = total[replicates].mean(axis=1)
    out = pd.DataFrame({"gene_id": total["Gene ID"].astype(str).to_numpy(),
                        "gametocyte_proteome_log2": np.log2(mean.where(mean > 0)).to_numpy()})
    out = out.groupby("gene_id", as_index=False).mean()
    made = pd.read_excel(path, sheet_name="Extended Data Table 3", header=2)
    made = made[made["Gene ID"].astype(str).str.match(valid)]
    # Group 4 is the set found ONLY with the label -- newly made protein. Every other group is
    # found in a control too: 7 and 6 also appear without the label, 5 also appears when synthesis
    # is blocked, and 1 and 2 are control-only. Taking the whole sheet would call 1,179 proteins
    # newly made where the paper calls 705.
    made = set(made.loc[made["Group"] == 4, "Gene ID"].astype(str))
    # The union of the two tables, not the abundance table alone: a protein the label found but the
    # abundance run did not quantify is still newly made, and dropping it lost 16 of the paper's 705.
    out = out.set_index("gene_id").reindex(sorted(set(out["gene_id"]) | made))
    out.index.name = "gene_id"
    out["gametocyte_newly_made"] = out.index.isin(made).astype(float)
    return out.reset_index()


PF_TPP = ("post_translation", "thermal", "pf_tpp_2026", "SupplTable3_media-6.csv")


def pf_target_engagement(root: str) -> pd.DataFrame:
    """How many antimalarials measurably engage each protein (bioRxiv 2026, Supplementary Table 3):
    thermal proteome profiling against 25 compounds, as a dose series at fixed temperatures.

    Counts, not a continuous response, and that is the finding rather than a simplification: the
    continuous response agrees between replicates at r 0.0-0.55, so a per-compound number would be
    mostly noise, while the call of hit or not reproduces the paper (99 stabilised hits, and
    cladosporine-KRS, MMV665915-ACS10 and KAF156-prohibitin all in the top thirteen). The
    denominator ships with the count: a protein seen in five assays and hit twice is not the same
    evidence as one seen in 25 and hit twice.

    NOT a melting temperature. This design has no temperature series, so it cannot give one; the
    melting points come from the separate MAP-X meltome.
    """
    path = _file(root, *PF_TPP)
    if path is None:
        return pd.DataFrame()
    d = pd.read_csv(path)
    d["gene_id"] = d["id"].astype(str).str.strip()
    d = d[d["gene_id"].str.match(r"^PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})$")]
    grouped = d.groupby("gene_id")
    out = pd.DataFrame({"engaged_n_compounds_tested": grouped.size().astype(float),
                        "engaged_n_compounds_hit": grouped["hit"].apply(
                            lambda s: float((s == "hit").sum()))})
    out.index.name = "gene_id"
    return out.reset_index()


PF_FEBRILE = ("post_translation", "phosphosites", "pf_febrile_2026", "elife-107860-supp1-v1.xlsx")
PF_FEBRILE_FC = "Log2FC [WT (Post Heat Stress) /WT (No Heat Stress)]"


def pf_febrile_phospho(root: str) -> pd.DataFrame:
    """How each protein's phosphorylation changes in a fever (eLife 2026, Supplementary File 1):
    infected red cells held at 39 degrees, then returned to 37 and harvested.

    Per protein, the strongest change at any of its sites and how many sites moved in each
    direction. The paper's regulated set reproduces exactly at 143 up and 53 down -- but those are
    counts of ROWS, and a row is a site at one phosphorylation multiplicity, so a site seen doubly
    and singly phosphorylated is counted twice. Collapsed to unique sites, each represented by its
    strongest change, it is 128 up and 51 down, which is what the columns here count. The response
    is specific rather than global: 60% of the proteins with a rising site are exported to the host
    cell, against 7% of everything quantified.

    Thirty-nine degrees, not forty, whatever the preprint's abstract says -- the methods give the
    temperature the cultures were held at. A phosphosite whose peptide could belong to more than one
    gene is dropped rather than assigned to the first.
    """
    import re
    path = _file(root, *PF_FEBRILE)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name="1 Pf Phosphoproteome")
    valid = re.compile(r"PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})")
    d["gene_id"] = d["Leading Protein ID"].astype(str).str.extract(f"({valid.pattern})")[0]
    matches = [c for c in d.columns if str(c).strip() == "Possible Protein ID Matches"]
    if matches:
        genes = d[matches[0]].astype(str).map(lambda s: len(set(valid.findall(s))))
        d = d[genes <= 1]
    change = pd.to_numeric(d[PF_FEBRILE_FC], errors="coerce")
    d = d.assign(change=change).dropna(subset=["gene_id", "change"])
    position = [c for c in d.columns if str(c).strip().startswith("Positions within protein")]
    d["site"] = d["gene_id"] + "_" + d[position[0]].astype(str) if position else d["gene_id"]
    # One site measured at several multiplicities is one site: its strongest change represents it.
    per_site = d.groupby(["gene_id", "site"])["change"].agg(
        lambda x: x.loc[x.abs().idxmax()])
    by_gene = per_site.groupby(level=0)
    out = pd.DataFrame({"febrile_phospho_max_abs_log2fc": by_gene.apply(lambda x: x.abs().max()),
                        "febrile_phospho_n_sites_up": by_gene.apply(
                            lambda x: float((x > 1).sum())),
                        "febrile_phospho_n_sites_down": by_gene.apply(
                            lambda x: float((x < -1).sum()))})
    out.index.name = "gene_id"
    return out.reset_index()


PF_CHROMATIN = ("post_translation", "BioID", "pf_chromatin_proxiome_2026")
#: bait -> (file, sheet, enrichment column, -log10 p column, unique-peptide column, the paper's own
#: enrichment threshold for a high-confidence hit). Four of the seven baits ship: the two marks that
#: define euchromatin, the heterochromatin protein, and the centromere. The thresholds are the
#: paper's Methods verbatim, and they differ per bait because the baits differ in signal -- one
#: threshold for all of them would be tidier and wrong.
PF_CHROMATIN_BAITS = {
    "hp1": ("TableS1_HP1_proximity_proteome_media-3.xlsx", "Sheet C. HP1 protA-turboID",
            "t-test Difference", "-Log t-test p value", "Unique peptides", 1.5),
    "h3k27ac": ("TableS5_euchromatin_proximity_proteome_media-7.xlsx", "Sheet B. H3K27ac",
                "t-test Difference", "-Log t-test p value", "Unique peptides", 1.5),
    "h3k4me3": ("TableS5_euchromatin_proximity_proteome_media-7.xlsx", "Sheet C. H3K4me3",
                "t-test Difference", "-Log t-test p value", "Unique peptides", 0.8),
    "centromere": ("TableS7_centromeric_proximity_proteome_media-9.xlsx",
                   "Sheet A. CENH3 (mT)-DiQ-BioIDs", "T-test Difference", "-Log T-test p-value",
                   "Unique peptides", 1.0),
}
PF_CHROMATIN_MIN_P = 1.3            # -log10 p, the paper's cut
PF_CHROMATIN_MIN_PEPTIDES = 4


def pf_chromatin_proxiome(root: str) -> pd.DataFrame:
    """Which proteins sit in each chromatin environment (bioRxiv 2026, Tables S1, S5 and S7):
    proximity labelling from baits that mark heterochromatin, active chromatin and the centromere.

    Per bait, how enriched a protein is over the bait's own control and whether it clears the
    paper's threshold for that bait. The high-confidence counts reproduce the paper exactly
    (H3K27ac 99, H3K4me3 48, heterochromatin 61), and the positive controls land where they must:
    HP1, HDA2, GDV1 and AP2-HC in heterochromatin, CENH3 at the centromere, BDP1 and GCN5 in active
    chromatin. SIR2A is in none of them -- it was not detected by any bait, which is absence of
    evidence, not evidence of absence.

    A protein group spanning more than one gene is dropped rather than assigned to the first of
    them; where one gene appears twice in a bait, the row with more unique peptides represents it.
    A protein a bait never detected keeps a missing flag, not a zero.
    """
    import re
    folder = _file(root, *PF_CHROMATIN)
    if folder is None:
        return pd.DataFrame()
    valid = re.compile(r"PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})")

    def gene_of(value):
        genes = sorted(set(valid.findall(str(value))))
        return genes[0] if len(genes) == 1 else None

    frames = []
    for bait, (name, sheet, diff, p, peptides, threshold) in PF_CHROMATIN_BAITS.items():
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            continue
        d = pd.read_excel(path, sheet_name=sheet)
        ids = [c for c in d.columns if "ajority" in str(c) and "rotein" in str(c)][0]
        d = d[d[ids].notna()].copy()
        d["gene_id"] = d[ids].map(gene_of)
        d = d[d["gene_id"].notna()]
        d = d.sort_values(peptides, ascending=False).drop_duplicates("gene_id")
        enrichment = pd.to_numeric(d[diff], errors="coerce")
        significance = pd.to_numeric(d[p], errors="coerce")
        unique = pd.to_numeric(d[peptides], errors="coerce")
        t = pd.DataFrame({f"chromprox_{bait}_log2fc": enrichment.to_numpy(),
                          f"chromprox_{bait}_hit": (
                              (enrichment >= threshold) & (significance > PF_CHROMATIN_MIN_P)
                              & (unique >= PF_CHROMATIN_MIN_PEPTIDES)).astype(float).to_numpy()},
                         index=pd.Index(d["gene_id"].to_numpy(), name="gene_id"))
        frames.append(t)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, axis=1)
    for bait in PF_CHROMATIN_BAITS:
        if f"chromprox_{bait}_log2fc" in out.columns:
            out.loc[out[f"chromprox_{bait}_log2fc"].isna(), f"chromprox_{bait}_hit"] = np.nan
    return out.reset_index()


# ------------------------------------------------------------------- Plasmodium, fourth wave
#: A nuclear, apicoplast or mitochondrial gene of P. falciparum 3D7. Written once here because five
#: of the six deposits below filter on it.
PF_ACCESSION = r"^PF3D7_(?:\d{7}|API\d{5}|MIT\d{5})$"

PF_LOPIT = ("post_translation", "LOPIT", "42218142")
PF_LOPIT_SUMMARY = "41467_2026_73664_MOESM3_ESM.xlsx"
PF_LOPIT_VARIATION = "41467_2026_73664_MOESM5_ESM.xlsx"
#: The two-experiment map, not the three-experiment one. S1-S2 is the paper's headline classifier --
#: 1,646 of 3,000 proteins into 24 niches -- while S1-S2-S3 covers 2,496 proteins and classifies
#: fewer of them (1,323), because adding a third replicate to the concatenation costs proteins that
#: were not quantified in it. Both are in the deposit; shipping both would be one measurement twice.
PF_LOPIT_LOCATION = "Final location (S1-S2)"
PF_LOPIT_SCORE = "svm.score (S1-S2)"


def pf_spatial_proteome(root: str) -> pd.DataFrame:
    """Where each P. falciparum protein sits in a late schizont (Chisholm et al., Nat Commun 2026,
    Supplementary Data 1 and 3): hyperLOPIT over density gradients, classified by a support-vector
    machine into 24 cellular niches.

    The first subcellular localization of any kind in the Plasmodium table, which until now had no
    compartment at all. `lopit_pf_location` is the label and `lopit_pf_svm_score` the classifier's
    confidence in it; `unknown` becomes a MISSING label rather than a 25th niche, because the map's
    first rule is that absence must not read as a measurement.

    Reproduces the paper exactly: 3,000 proteins, 1,646 classified, 24 niches. The markers land
    where they must -- RAP1 in the rhoptries, MAHRP1 in the Maurer's clefts, ACP in the apicoplast,
    EXP2 and HSP101 at the parasitophorous vacuole membrane, GAPDH in the cytosol.

    The same paper's variation table is a separate deposit (`pf_field_variation`) rather than more
    columns here, because holding out a localization must not also withdraw a dN/dS: the closure
    excludes columns that share provenance, and a genomic analysis printed in the same paper is not
    an output of this gradient experiment.
    """
    folder = _file(root, *PF_LOPIT)
    if folder is None:
        return pd.DataFrame()
    summary = pd.read_excel(os.path.join(folder, PF_LOPIT_SUMMARY),
                            sheet_name="PlasmoLOPIT_data_summary")
    label = summary[PF_LOPIT_LOCATION].astype(str).str.strip()
    out = pd.DataFrame({"gene_id": summary["Gene accession"].astype(str).str.strip(),
                        "lopit_pf_location": label.where(label.str.lower() != "unknown"),
                        "lopit_pf_svm_score": pd.to_numeric(summary[PF_LOPIT_SCORE],
                                                            errors="coerce")})
    return out[out["gene_id"].str.match(PF_ACCESSION)].reset_index(drop=True)


def pf_field_variation(root: str) -> pd.DataFrame:
    """How variable each gene is in the field, and how fast it has evolved between species
    (Chisholm et al., Nat Commun 2026, Supplementary Data 3).

    `field_pnps_adj` and `field_variant_fraction` are computed over FIELD isolates, which is the
    question the map had no data for -- the SNP columns already shipped are PlasmoDB's counts over
    laboratory strains, and the two agree at rho 0.29, a related quantity on a different population
    rather than a copy of it. `dnds_laverania` and `dnds_plasmodium` are ratios against orthologs in
    other species: selection since the species split, not variation inside this one.
    """
    folder = _file(root, *PF_LOPIT)
    if folder is None:
        return pd.DataFrame()
    variation = pd.read_excel(os.path.join(folder, PF_LOPIT_VARIATION),
                              sheet_name="Plasmodium falciparum variation")
    var = pd.DataFrame({"gene_id": variation["Gene accession"].astype(str).str.strip(),
                        "field_pnps_adj": pd.to_numeric(variation["Adj. pN/pS"], errors="coerce"),
                        "field_variant_fraction": pd.to_numeric(variation["Variant fraction"],
                                                                errors="coerce"),
                        "dnds_laverania": pd.to_numeric(variation["Avg. dN/dS (L)"],
                                                        errors="coerce"),
                        "dnds_plasmodium": pd.to_numeric(variation["Avg. dN/dS (P)"],
                                                         errors="coerce")})
    return var[var["gene_id"].str.match(PF_ACCESSION)].reset_index(drop=True)


PF_4TU = ("transcription", "RNAseq", "29985403", "41467_2018_4966_MOESM5_ESM.xlsx")


def pf_mrna_dynamics(root: str) -> pd.DataFrame:
    """How fast each transcript is made and how fast it disappears (Painter et al., Nat Commun 2018,
    Supplementary Data 2): 4-thiouracil labelling every hour through the 48-hour blood-stage cycle,
    with the rates read off a model of total and nascent RNA.

    Transcripts per minute at the gene's own peak, so these are FLUXES and not half-lives: an
    abundant transcript decays faster in transcripts per minute than a scarce stable one, which is
    why the decay column correlates with abundance (rho 0.42) and why it is not comparable with the
    Toxoplasma actinomycin columns, which are fractions remaining. Signs as published: transcription
    positive (4,414 of 4,416), decay negative (4,443 of 4,462).

    The paper's claim is active transcription in every stage of the cycle, and the deposit carries
    it: 616 to 962 genes take their peak transcription in each of the six windows. The timing checks
    out against a series this study had nothing to do with -- of the genes whose transcription peaks
    in a ring window, 87% to 92% also peak in the shipped ring expression column.

    The peak-stage labels are in the deposit and are deliberately NOT shipped: which stage a gene's
    transcription peaks in is what the stage expression series in the table already says, and a
    second copy of it would be a column the leakage closure has to guess at.
    """
    path = _file(root, *PF_4TU)
    if path is None:
        return pd.DataFrame()
    rate = pd.read_excel(path, sheet_name="Transcription Rate")
    decay = pd.read_excel(path, sheet_name="Decay Rate")
    rate.columns = ["gene_id", "peak_stage", "transcription_rate_4tu"]
    decay.columns = ["gene_id", "peak_stage", "mrna_decay_rate_4tu"]
    frames = []
    for frame, column in ((rate, "transcription_rate_4tu"), (decay, "mrna_decay_rate_4tu")):
        t = pd.DataFrame({"gene_id": frame["gene_id"].astype(str).str.strip(),
                          column: pd.to_numeric(frame[column], errors="coerce")})
        frames.append(t[t["gene_id"].str.match(r"^PF3D7_\w+$")].set_index("gene_id"))
    out = frames[0].join(frames[1], how="outer")
    out.index.name = "gene_id"
    return out.reset_index()


PF_HSP90 = ("translation", "proteomics", "pf_hsp90_2026")
PF_HSP90_TABLE = "TableS1_media-2.xlsx"
PF_HSP90_HITS = "TableS2_media-3.xlsx"
PF_HSP90_REPORT = "20240914_094033_10643_hu_pf_loc_091424_Report_protein_abundances.tsv"
#: The paper's own rule for a chaperone-dependent protein, from its Table S2A: down by at least half
#: a log2 and p < 0.05 under BOTH inhibitors. It reproduces their list exactly -- 131 proteins, all
#: 131 of them -- which is why it is written here as a rule rather than the list being copied.
PF_HSP90_MIN_DROP = -0.5
PF_HSP90_MAX_P = 0.05


def pf_hsp90_chemoproteome(root: str) -> pd.DataFrame:
    """The asexual blood-stage proteome, and what falls when Hsp90 is inhibited (Ibrasheva et al.,
    bioRxiv 2026, Tables S1 and S2; PRIDE PXD079493).

    Three columns for three questions. `proteome_blood_log2` is how much of each protein a
    DMSO-treated culture contains -- log2 of the mean of the four control replicates' normalized
    abundances -- which is the first protein abundance in the map for the stage the parasite spends
    its life in. The two fold changes are what geldanamycin and XL888 do to it, two inhibitors with
    different scaffolds; their continuous responses agree only at rho 0.17, which is why the hit
    call demands BOTH -- what survives that conjunction is the chaperone's, while either column
    alone is largely one compound's. `hsp90_dependent` is the paper's own hit call, reproduced from
    the fold changes and p-values by the rule its methods state.

    The deposit gives UniProt accessions and the map is keyed by gene, so the mapping comes from the
    deposit's OWN Spectronaut report (`PG.UniProtIds` beside `PG.Genes`), never from a lookup that
    could have drifted: 3,050 of 3,265 rows resolve, and a protein group naming more than one gene
    is dropped rather than assigned to the first.

    The PRIDE record says 133 hits and the supplementary table lists 131. The 131 are what is
    reproducible, so they are what ships, and the difference is recorded rather than reconciled.
    """
    folder = _file(root, *PF_HSP90)
    if folder is None:
        return pd.DataFrame()
    valid = re.compile(PF_ACCESSION)
    report = pd.read_csv(os.path.join(folder, PF_HSP90_REPORT), sep="\t",
                         usecols=["PG.Genes", "PG.UniProtIds"])
    mapping = {}
    for genes, accessions in zip(report["PG.Genes"].fillna("").astype(str),
                                 report["PG.UniProtIds"].fillna("").astype(str)):
        named = {g.split(".")[0] for g in genes.split(";") if valid.match(g.split(".")[0])}
        if len(named) != 1:
            continue                      # a group spanning two genes names neither
        gene = named.pop()
        for accession in accessions.split(";"):
            mapping[accession.strip()] = gene
    table = pd.read_excel(os.path.join(folder, PF_HSP90_TABLE), sheet_name=0, header=1)
    control = [c for c in table.columns if "DMSO" in str(c) and "N_abund" in str(c)]
    ga = pd.to_numeric(table["Average lg2fc (GA)"], errors="coerce")
    xl = pd.to_numeric(table["Average lg2fc (XL)"], errors="coerce")
    p_ga = pd.to_numeric(table["p-value (GA)"], errors="coerce")
    p_xl = pd.to_numeric(table["p-value (XL)"], errors="coerce")
    mean = table[control].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    out = pd.DataFrame({"gene_id": table["UniProt ID"].astype(str).str.strip().map(mapping),
                        "proteome_blood_log2": np.log2(mean.where(mean > 0)),
                        "hsp90_inhibition_ga_log2fc": ga,
                        "hsp90_inhibition_xl_log2fc": xl,
                        "hsp90_dependent": ((ga < PF_HSP90_MIN_DROP) & (p_ga < PF_HSP90_MAX_P)
                                            & (xl < PF_HSP90_MIN_DROP)
                                            & (p_xl < PF_HSP90_MAX_P)).astype(float)})
    out = out.dropna(subset=["gene_id"])
    # The hit call takes the MAXIMUM where two protein groups resolve to one gene: averaging a yes
    # with a no would write 0.5, which the column's own definition has no meaning for.
    return out.groupby("gene_id", as_index=False).agg(
        {"proteome_blood_log2": "mean", "hsp90_inhibition_ga_log2fc": "mean",
         "hsp90_inhibition_xl_log2fc": "mean", "hsp90_dependent": "max"})


PF_RDEEP = ("post_translation", "RDeeP", "38355719", "41467_2024_45519_MOESM4_ESM.xlsx")
PF_RDEEP_SHEET = "Left-shifting"
PF_RDEEP_HEADER = 5
PF_RDEEP_MAX_Q = 0.05


def pf_rna_dependence(root: str) -> pd.DataFrame:
    """Whether a protein needs RNA to stay in the complex it is in (Hollin et al., Nat Commun 2024,
    Supplementary Data 2): R-DeeP -- the same lysate run down a sucrose gradient with and without
    RNase, 25 fractions each, and a protein whose distribution shifts towards the light fractions
    when the RNA is gone was held there by RNA.

    A question nothing in the map asked. It is not RNA binding: a protein can shift because its
    PARTNER binds RNA, which is why the study calls them RNA-DEPENDENT rather than RNA-binding.

    `rna_dependent` is 1 for the 898 significantly left-shifted proteins and 0 for the rest of the
    3,671 the run quantified -- a zero here means tested and not shifted, which is why the flag
    ships beside the q-value rather than as a list. The 898 at q < 0.05 is the paper's own number,
    reproduced from the deposit.

    The classes come out the right way round: RNA helicases are enriched among the dependent (28 of
    61, odds 2.7, p 2e-4) and the proteasome, a large complex that needs no RNA, has 0 of 14.
    Ribosomal proteins are DEPLETED rather than enriched (23 of 137, odds 0.6) -- a ribosome is
    mostly RNA and does not come apart in this assay the way an mRNA-bridged complex does -- so this
    column must not be read as 'binds RNA'.
    """
    path = _file(root, *PF_RDEEP)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name=PF_RDEEP_SHEET, header=PF_RDEEP_HEADER)
    q = pd.to_numeric(d["q-value"], errors="coerce")
    out = pd.DataFrame({"gene_id": d["Protein"].astype(str).str.strip(),
                        "rna_dependence_qvalue": q,
                        "rna_dependent": (q < PF_RDEEP_MAX_Q).astype(float)})
    out = out[out["gene_id"].str.match(PF_ACCESSION)]
    # A protein quantified twice keeps its strongest evidence rather than the mean of two q-values.
    return out.sort_values("rna_dependence_qvalue").drop_duplicates(
        "gene_id").reset_index(drop=True)


PF_COMMITTED = ("translation", "proteomics", "41482054", "mmc4.xlsx")


def pf_committed_proteome(root: str) -> pd.DataFrame:
    """What changes in a parasite that has committed to becoming a gametocyte but still looks
    asexual (Venugopal et al., Mol Cell Proteomics 2026, Table S3): parasites sorted on MSRP1, a
    marker this study establishes, and compared with their uncommitted siblings.

    Commitment, not the gametocyte -- the map already has a stage V gametocyte proteome, and this is
    the decision that precedes it by a cycle. Positive means more abundant in the committed cells.
    The columns are the study's combined contrast over both reporter lines and both sorting
    directions (`logfc.all`, `fdr.all`), because the four arms are the same comparison done four
    ways.

    MSRP1 itself is +1.68 at FDR 0, which is a positive control and nothing more: it is the
    protein the sort was done on. The check that means something is the paper's finding that
    merozoite surface proteins separate the two populations: MSP1 (+0.34, FDR 0) and MSP2 (+0.67)
    are up in the committed cells.
    """
    path = _file(root, *PF_COMMITTED)
    if path is None:
        return pd.DataFrame()
    d = pd.read_excel(path, sheet_name="stats")
    out = pd.DataFrame({"gene_id": d["protein"].astype(str).str.strip(),
                        "committed_vs_asexual_log2fc": pd.to_numeric(d["logfc.all"],
                                                                     errors="coerce"),
                        "committed_vs_asexual_fdr": pd.to_numeric(d["fdr.all"], errors="coerce")})
    out = out[out["gene_id"].str.match(PF_ACCESSION)]
    return out.groupby("gene_id", as_index=False).mean()


PF_RESISTOME = ("DNA", "in_vitro_evolution", "39607932")
PF_RESISTOME_SNVS = "SupplementaryData3_SNVs-INDELs.xlsx"
PF_RESISTOME_GENES = "SupplementaryData6_Pf6Genes.xlsx"
PF_RESISTOME_TARGETS = "SupplementaryData5_CompoundTarget.xlsx"
#: The two classes in which the paper's hypergeometric test supports the call. `Other Evidence`,
#: `Unclear`, `Ambiguous CNV` and the copy-number-only class are left out: a count that mixed a
#: statistically supported target with a compound whose own authors wrote `Unclear` would be a
#: curated column that curates nothing.
PF_RESISTOME_SUPPORTED = ("Overrepresented SNVs", "Overrepresented SNVs and CNVs")
#: The deposit names the three mitochondrial-genome genes by their pre-2010 identifiers, which match
#: no accession pattern and would be dropped -- and one of them is apocytochrome b, the atovaquone
#: resistance gene, so dropping it would lose the best-known entry in the table. Each is mapped on
#: the deposit's OWN gene description matching exactly one product in the shipped table:
#: `apocytochrome B` -> the single gene whose product is `cytochrome b`, `cytochrome c oxidase
#: subunit I (cox1)` and `... subunit 3 (COX3)` likewise. Written out rather than pattern-matched,
#: because a legacy name resolved by guesswork attaches a measurement to the wrong protein.
PF_RESISTOME_LEGACY = {"mal_mito_1": "PF3D7_MIT01400",     # cytochrome c oxidase subunit 3
                       "mal_mito_2": "PF3D7_MIT02100",     # cytochrome c oxidase subunit 1
                       "mal_mito_3": "PF3D7_MIT02300"}     # apocytochrome b
#: Variant classes that change the protein. An intergenic or intron variant is recorded against the
#: nearest gene by the deposit's annotation and is not evidence that this gene's product changed, so
#: the counts are of coding variants and the rest are left out rather than silently included.
PF_RESISTOME_CODING = ("missense_variant", "frameshift_variant", "stop_gained", "inframe_deletion",
                       "disruptive_inframe_deletion", "inframe_insertion",
                       "disruptive_inframe_insertion", "stop_lost", "start_lost",
                       "splice_acceptor_variant", "splice_donor_variant")


def pf_resistome(root: str) -> pd.DataFrame:
    """Which genes mutate when the parasite is selected for drug resistance, and how variable they
    are in the field (Luth et al., Science 2024, Supplementary Data 3 and 6): 724 clones, each
    evolved in vitro against one of 118 compounds and sequenced whole-genome.

    Four columns. Three are counts over CODING variants only: how many independently selected clones
    carry one, how many distinct compounds those selections used, and how many distinct variants
    there are. The fourth, `resistance_target_compounds`, is the paper's own call -- for how many
    compounds is this gene the classified target or resistance gene, counting only the classes its
    hypergeometric test supports -- and it is the column that answers the slot, because a COUNT of
    mutations is not a claim that the gene confers resistance.

    The difference is visible at the top of each column and is the reason both ship. By raw
    selections the leaders are AP2-G (13 compounds) and PfEMP1, which mutate under prolonged culture
    whatever the drug: losing gametocyte production is cheap in a flask. By the paper's
    classification the leaders are PfATP4 and PfMDR1 (5 compounds each), the prodrug-activating
    esterase (4), then cytochrome b, CARL, PI4K beta and the Niemann-Pick C1-related protein at 3,
    with DHODH, PfCRT, the tRNA ligases and DHFR-TS behind them -- the canonical antimalarial
    resistance set, arrived at without being told what it should contain.

    A gene with no selected mutation gets NO value rather than a zero. That is deliberate and it is
    the opposite of the choice made for the R-DeeP flag: a screen that tested every gene can say
    'not a hit', whereas 118 compounds are not a test of the other 4,600 genes.

    Beside them, the same paper's field variation from 5,970 Pf6 isolates: `pf6_field_dnds` and the
    non-synonymous SNV count. The deposit writes -1 where a ratio could not be computed, which
    becomes missing here -- carried through as a number it would be a strong negative selection
    signal that nobody measured.
    """
    folder = _file(root, *PF_RESISTOME)
    if folder is None:
        return pd.DataFrame()
    valid = re.compile(PF_ACCESSION)
    snv = pd.read_excel(os.path.join(folder, PF_RESISTOME_SNVS), header=1)
    snv["gene_id"] = snv["Gene Name"].astype(str).str.strip().replace(PF_RESISTOME_LEGACY)
    snv = snv[snv["gene_id"].str.match(valid)
              & snv["Effect"].astype(str).str.strip().isin(PF_RESISTOME_CODING)]
    counts = snv.groupby("gene_id").agg(
        resistance_selection_clones=("Clone Name", "nunique"),
        resistance_selection_compounds=("Compound", "nunique"),
        resistance_selection_variants=("Variant ID", "nunique")).astype(float)
    genes = pd.read_excel(os.path.join(folder, PF_RESISTOME_GENES),
                          sheet_name="evolved_vs_pf6_dNdS", header=1)
    field = pd.DataFrame({"gene_id": genes["Gene ID"].astype(str).str.strip(),
                          "pf6_field_dnds": pd.to_numeric(genes["Pf6 dN/dS"], errors="coerce"),
                          "pf6_field_nonsyn_snvs": pd.to_numeric(
                              genes["Numer of Pf6 Non-synonymous unique SNVs"], errors="coerce")})
    field = field[field["gene_id"].str.match(valid)].set_index("gene_id")
    field["pf6_field_dnds"] = field["pf6_field_dnds"].where(field["pf6_field_dnds"] >= 0)
    targets = pd.read_excel(os.path.join(folder, PF_RESISTOME_TARGETS), header=1)
    supported = targets[targets["Evidence classification"].astype(str).str.strip()
                        .isin(PF_RESISTOME_SUPPORTED)].copy()
    names = supported["Target or Resistance Gene ID"].astype(str).str.strip()
    supported["gene_id"] = names.replace(PF_RESISTOME_LEGACY)
    supported = supported[supported["gene_id"].str.match(valid)]
    called = supported.groupby("gene_id")["Compound"].nunique().astype(float).rename(
        "resistance_target_compounds")
    out = counts.join(field, how="outer").join(called, how="outer")
    out.index.name = "gene_id"
    return out.reset_index()


# --------------------------------------------------------------------------- the list
@dataclass(frozen=True)
class Deposit:
    """One deposit: which registry entry it belongs to, whose table it lands in, and how it is read.

    A dataclass rather than a dict so that adding a deposit is adding a line, and so the organism
    cannot be forgotten -- a Plasmodium table merged into the Toxoplasma one would join zero rows
    and look exactly like a dataset with no coverage.
    """
    key: str                 # the registry entry it belongs to
    organism: str            # "Tg", "Pf" or "host"
    derive: Callable         # dataset_root -> frame; first column gene_id, or host_id + host_name
    combine: str = "mean"    # rows resolving to one gene: "mean" for replicate-like rows


DEPOSITS = (
    Deposit("crispr_invivo_composite", "Tg", giuliano_invivo),
    Deposit("crispr_serum_restriction", "Tg", serum_restriction),
    Deposit("crispr_glucose_limitation", "Tg", glucose_limitation),
    Deposit("gse302107_riboseq", "Tg", riboseq_302107),
    Deposit("gse302108_utr5", "Tg", utr5_architecture),
    Deposit("mrna_decay_gse329845", "Tg", mrna_decay_329845),
    Deposit("brady_subtypes", "Tg", bradyzoite_subtypes),
    Deposit("iron_depletion_proteome", "Tg", iron_depletion),
    Deposit("organelle_surface_turboid", "Tg", organelle_surface),
    Deposit("pf_meltome", "Pf", pf_melting_temperature),
    Deposit("pf_pb_fertility_transfer", "Pf", pb_fertility),
    Deposit("pf_latency_transcriptome", "Pf", pf_latency),
    Deposit("pf_m6a_nanopore", "Pf", pf_m6a),
    Deposit("pf_gametocyte_proteome", "Pf", pf_gametocyte_proteome),
    Deposit("pf_target_engagement", "Pf", pf_target_engagement),
    Deposit("pf_febrile_phospho", "Pf", pf_febrile_phospho),
    Deposit("pf_chromatin_proxiome", "Pf", pf_chromatin_proxiome),
    Deposit("pf_spatial_proteome", "Pf", pf_spatial_proteome),
    Deposit("pf_field_variation", "Pf", pf_field_variation),
    Deposit("pf_mrna_dynamics", "Pf", pf_mrna_dynamics),
    Deposit("pf_hsp90_chemoproteome", "Pf", pf_hsp90_chemoproteome),
    Deposit("pf_rna_dependence", "Pf", pf_rna_dependence),
    Deposit("pf_committed_proteome", "Pf", pf_committed_proteome),
    Deposit("pf_resistome", "Pf", pf_resistome),
    Deposit("host_hff_tg_infection", "host", hff_tg_infection),
    Deposit("host_bmdm_baseline", "host", bmdm_baseline),
    Deposit("host_hepatocyte_pf_infection", "host", hepatocyte_pf_infection),
    Deposit("host_k562_rhoptry_screen", "host", k562_rhoptry_screen),
)


def table_path(base: str, key: str) -> str:
    """Where a derived deposit table is cached, so the merge can run without the raw files."""
    return os.path.join(base, "starplast", "data", f"deposit_{key}.tsv")


def derive_all(dataset_root: str, base: str, keys=None, log=print) -> dict:
    """Derive every deposit present under `dataset_root` and write its keyed table into the cache.

    A deposit absent from the root is skipped and its shipped table, if any, left alone -- a
    machine without the raw files must not erase what one with them derived.
    """
    out = {}
    for dep in DEPOSITS:
        if keys and dep.key not in keys:
            continue
        frame = dep.derive(dataset_root)
        if frame is None or frame.empty:
            log(f"deposit {dep.key}: raw files not under {dataset_root}; kept what is shipped")
            continue
        frame.to_csv(table_path(base, dep.key), sep="\t", index=False, float_format="%.10g")
        out[dep.key] = frame
        log(f"deposit {dep.key}: {len(frame):,} rows x {frame.shape[1] - 1} columns")
    return out


def _read(base: str, key: str) -> pd.DataFrame:
    path = table_path(base, key)
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else pd.DataFrame()


def parasite_columns(base: str, organism: str, ids, resolve=None, log=print) -> pd.DataFrame:
    """Every shipped deposit column for one organism, aligned to `ids`.

    Accessions go through `resolve` (the identity layer: GT1 to ME49, old to current). Where rows
    land on one gene: an accession that is the gene itself (TGGT1_224540 for TGME49_224540) wins
    over the A/B halves of a locus GT1 splits (TGGT1_224540A/B); only a gene reached by halves
    alone gets their mean. Averaging all three once moved TGME49_224540 from -3.99 to +0.20 --
    the halves are different gene models, not replicates of the whole.

    A LABEL column (a compartment name, a stage) has no mean, so where two rows land on one gene it
    is carried only if they agree; where they disagree the gene gets no label rather than the first
    row's. Added for the Plasmodium spatial proteome, whose measurement IS a label -- before that
    every deposit column was numeric and `mean(numeric_only=True)` dropped anything else in silence.
    """
    ids = pd.Index([str(i) for i in ids])
    out = pd.DataFrame(index=ids)
    for dep in DEPOSITS:
        if dep.organism != organism:
            continue
        frame = _read(base, dep.key)
        if frame.empty:
            continue
        source = frame["gene_id"].astype(str)
        genes = source.map(lambda g: resolve(g) or g) if resolve is not None else source
        half = source.str.contains(r"\d[A-Z]$").to_numpy()
        whole_genes = set(genes[~half])
        keep = ~half | ~genes.isin(whole_genes).to_numpy()
        values = frame.drop(columns="gene_id")[keep].set_index(genes[keep].to_numpy())
        labels = [c for c in values.columns if not pd.api.types.is_numeric_dtype(values[c])]
        numeric = values.groupby(level=0).mean(numeric_only=True)
        if labels:
            agreed = values[labels].groupby(level=0).agg(
                lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else None)
            numeric = numeric.join(agreed)[[c for c in values.columns]]
        aligned = numeric.reindex(ids)
        for c in aligned.columns:
            out[c] = aligned[c].to_numpy()
        log(f"deposit {dep.key}: " + ", ".join(
            f"{c} {int(aligned[c].notna().sum()):,}" for c in aligned.columns))
    return out


def host_columns(base: str) -> pd.DataFrame:
    """The host deposits merged on the host key, for `host.merge_tissue`."""
    from . import host
    out = pd.DataFrame()
    for dep in DEPOSITS:
        if dep.organism == "host":
            out = host.merge_tissue(out, _read(base, dep.key))
    return out
