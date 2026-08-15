#!/usr/bin/env python3
"""Published CRISPR screens and mass-spectrometry abundance, normalized to one gene per row.

Everything here is read from `datasets/crispr_screens/` and `toxo_stage_atlas/data/proteomics/`; see
`datasets/crispr_screens/SOURCES.md` for provenance. Four screens and two proteomes, each with a very
different scope, and the scope is the thing to keep in view:

* **GRA17 synthetic-lethal** (PMC10409377) -- the only genome-wide one. Phenotype scores for RH and
  RHdelta-gra17 by passage, so the *difference* is the synthetic-lethality signal.
* **GRA12 strains/subspecies** (PMC12003902) -- two targeted pooled screens of ~235 genes. They are not
  replicates: their in-vivo L2FCs correlate at r = 0.41, so they are kept separate.
* **in vivo CRISPR platform** (PMC6722137) -- targeted, 35-147 genes. Its peritoneum/lung/liver/spleen
  composite scores are almost certainly already in the tree as `fit_invivo_*`, so this is corroboration
  rather than new signal, and is loaded under its own names so the two can be compared instead of merged.
* **host-transcription effectors** (PMC12033024) -- a pooled single-cell screen scoring each effector by
  how much it perturbs host transcription (Hotelling T2).

**The screens do not share a sign convention or a scale.** That is established for the seven screens
already in the node table, and nothing here changes it: always rank-normalize and check direction before
pooling anything with anything.

**Mass spectrometry is thin.** The only abundance data in this tree is two Pru immunoprecipitation
experiments covering 424 and 594 proteins -- under 8% of the proteome, and enrichment rather than a deep
whole-cell proteome. It is loaded because it is real, and flagged because it is not coverage.
"""
from __future__ import annotations

import os
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

ACC = re.compile(r"(TGME49_\d{5,6})")


_RESOLVE = None      # set by crispr_screens/proteomics; see identity.GeneIndex


def _acc(series: pd.Series) -> pd.Series:
    """Pull a TGME49_ accession out of the supplement and resolve it to a current gene id.

    Supplements cite whatever accession was current when the paper was written. The 2019 in vivo screen
    uses pre-2012 `TGME49_0xxxxx` ids throughout, which match the regex, exist in no current node table,
    and silently drop every one of its 181 genes. Routing through the identity layer is the whole reason
    that layer exists.
    """
    raw = series.astype(str).str.extract(ACC, expand=False)
    if _RESOLVE is None:
        return raw
    return raw.map(lambda s: _RESOLVE(s) if isinstance(s, str) else s)


def _find(root: str, filename: str):
    """Locate a supplement whether it sits flat in `root` or one PMID directory down."""
    direct = os.path.join(root, filename)
    if os.path.exists(direct):
        return direct
    for sub in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        cand = os.path.join(root, sub, filename)
        if os.path.exists(cand):
            return cand
    return direct


def _read(path, **kw):
    """Read a sheet, or None if the file or that sheet is not there.

    A workbook that exists but lacks the expected sheet used to raise ValueError out of here and take
    the whole build down. Publishers reissue supplements with sheets renamed or replicates dropped, and
    one missing worksheet should cost that one screen, not every screen after it.
    """
    if not os.path.exists(path):
        return None
    try:
        return pd.read_excel(path, **kw)
    except ValueError:
        return None


def crispr_screens(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Return one row per gene, one column per published screen measurement.

    `resolve` maps a raw accession to a current gene id (see identity.GeneIndex); without it, papers
    that cite superseded accessions contribute nothing.
    """
    global _RESOLVE
    _RESOLVE = resolve
    # The screens were migrated into datasets/<level>/<type>/<PMID>/. Both layouts are searched, since
    # a flat crispr_screens/ is what an older checkout has and silently finding nothing is the failure
    # this exact move already caused once.
    roots = [os.path.join(base, "datasets", "DNA", "CRISPR_screen"),
             os.path.join(base, "datasets", "crispr_screens")]
    D = next((r for r in roots if os.path.isdir(r)), roots[0])
    out = []

    # --- GRA17 synthetic lethality: genome-wide, header on the second row
    f = _find(D, "gra17_synthlethal_PMC10409377_S1_phenotypes.xlsx")
    d = _read(f, sheet_name="TableS1_ALL_DATA", header=1)
    if d is not None:
        d = d.rename(columns=lambda c: str(c).strip())
        g = _acc(d["ME49ID"])
        cols = {}
        for src, dst in (("RHΔgra17 AVG_P3/P4", "crispr_gra17ko_phenotype"),
                         ("Δgra17-RH AVG_P3/P4", "crispr_gra17_synthlethal_delta")):
            if src in d.columns:
                cols[dst] = pd.to_numeric(d[src], errors="coerce")
        if "Synthetic lethality/viability candidate" in d.columns:
            cols["crispr_gra17_candidate"] = (
                d["Synthetic lethality/viability candidate"].notna()
                & d["Synthetic lethality/viability candidate"].astype(str).str.strip().ne("")).astype(int)
        t = pd.DataFrame(cols).assign(gene_id=g).dropna(subset=["gene_id"])
        out.append(t.groupby("gene_id").mean())
        log(f"screens: GRA17 synthetic-lethal (genome-wide) -> {len(out[-1]):,} genes")

    # --- GRA12: two targeted screens, deliberately not merged
    for tag, fn, sheet in (("s1", "gra12_PMC12003902_D3_gene_L2FC_screen1.xlsx", "2D.Gene L2FCs"),
                           ("s2", "gra12_PMC12003902_D4_gene_L2FC_screen2.xlsx", "3D.Gene L2FCs")):
        d = _read(_find(D, fn), sheet_name=sheet)
        if d is None:
            continue
        t = pd.DataFrame({
            f"crispr_gra12{tag}_l2fc_invitro": pd.to_numeric(d.MEDIAN_L2FC_IN_VITRO, errors="coerce"),
            f"crispr_gra12{tag}_l2fc_invivo": pd.to_numeric(d.MEDIAN_L2FC_IN_VIVO, errors="coerce"),
            f"crispr_gra12{tag}_disco": pd.to_numeric(d.DISCO_SCORE, errors="coerce"),
        }).assign(gene_id=_acc(d.GENE)).dropna(subset=["gene_id"])
        out.append(t.groupby("gene_id").mean())
        log(f"screens: GRA12 {tag} (targeted) -> {len(out[-1]):,} genes")

    # --- in vivo platform: several library sizes, pooled into one column per gene
    rows = []
    for fn, sheets in (("invivo_platform_PMC6722137_D2_phenotype_scores.xlsx",
                        ["200 - Phenotype scores", "800 - Phenotype scores"]),
                       ("invivo_platform_PMC6722137_D3_phenotype_scores.xlsx",
                        ["Phenotype scores"])):
        for sh in sheets:
            for hdr in (0, 1):
                d = _read(_find(D, fn), sheet_name=sh, header=hdr)
                if d is None or "Gene" not in [str(c).strip() for c in d.columns]:
                    continue
                d = d.rename(columns=lambda c: str(c).strip())
                lfc = [i for i, c in enumerate(d.columns) if "mean lfc" in str(c).lower()]
                if not lfc:
                    continue
                # These sheets repeat "mean lfc across replicates" once per condition, so selecting by
                # name returns a frame rather than a column. Take them positionally and average.
                v = d.iloc[:, lfc].apply(pd.to_numeric, errors="coerce").mean(axis=1)
                gene = d.iloc[:, [i for i, c in enumerate(d.columns)
                                  if str(c).strip() == "Gene"][0]]
                rows.append(pd.DataFrame({"gene_id": _acc(gene),
                                          "crispr_invivo_platform_lfc": v}))
                break
    if rows:
        t = pd.concat(rows).dropna(subset=["gene_id"])
        out.append(t.groupby("gene_id").mean())
        log(f"screens: in vivo platform (targeted) -> {len(out[-1]):,} genes")

    # --- hyperLOPIT-unassigned library: targeted in-vivo screen, GSE253884
    hyperlopit = []
    for accession, library in (("GSE253884", 1), ("GSE253885", 2)):
        acquired = os.path.join(base, "datasets", "toxoplasma_acquisition_2026_08_14",
                                f"{accession}_unassigned_{library}_summary.xlsx")
        d = _read(acquired, sheet_name="Summary")
        if d is None or not {"ME49_ID", "In vivo fitness score"} <= set(d.columns):
            continue
        library_frame = pd.DataFrame({
            "gene_id": _acc(d["ME49_ID"]),
            f"fit_hyperlopit_unassigned_invivo_lib{library}": pd.to_numeric(
                d["In vivo fitness score"], errors="coerce"),
        }).dropna(subset=["gene_id"]).groupby("gene_id").mean()
        hyperlopit.append(library_frame)
    if hyperlopit:
        t = pd.concat(hyperlopit, axis=1)
        out.append(t.groupby(level=0).mean())
        log(f"screens: hyperLOPIT-unassigned in vivo (GSE253884/5) -> "
            f"{len(out[-1]):,} genes")

    # --- host-transcription effectors: Target -> gene via the sgRNA map the paper ships
    gmap = _read(_find(D, "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx"))
    if gmap is not None:
        gmap = gmap.rename(columns=lambda c: str(c).strip())
        idcol = next((c for c in ("Gene_ID_Updated", "Gene_ID_Old") if c in gmap.columns), None)
        namecol = "Gene_Name" if "Gene_Name" in gmap.columns else None
        lookup = {}
        if idcol:
            for nm, gid in zip(gmap[namecol] if namecol else gmap[idcol], _acc(gmap[idcol])):
                if isinstance(gid, str):
                    lookup[str(nm).strip().upper()] = gid
        recs = []
        for fn in ("hosttx_effectors_PMC12033024_D4A_T2_statistic.xlsx",
                   "hosttx_effectors_PMC12033024_D5_T2_statistic.xlsx"):
            d = _read(_find(D, fn))
            if d is None:
                continue
            d = d.rename(columns=lambda c: str(c).strip())
            if "Target" not in d.columns:
                continue
            gid = _acc(d["Target"]).fillna(
                d["Target"].astype(str).str.strip().str.upper().map(lookup))
            recs.append(pd.DataFrame({
                "gene_id": gid,
                "hosttx_T2": pd.to_numeric(d.get("T2_Statistic"), errors="coerce"),
                "hosttx_padj": pd.to_numeric(d.get("P_Value_Adjusted"), errors="coerce"),
            }))
        if recs:
            t = pd.concat(recs).dropna(subset=["gene_id"])
            if len(t):
                out.append(t.groupby("gene_id").mean())
                log(f"screens: host-transcription effectors -> {len(out[-1]):,} genes")

    if not out:
        log("screens: none found")
        return pd.DataFrame()
    return pd.concat(out, axis=1)


def host_transcription_signatures(base: str, resolve=None, components: int = 20,
                                  log=print) -> pd.DataFrame:
    """Per-effector host-response signatures from the dual perturb-seq differential expression.

    The published table contains 33,531 host genes for each significant parasite effector. Putting
    thirty-three thousand mostly redundant columns into the node table would make the assay dominate
    by width, so the centered log-fold-change matrix is represented by its reproducible principal
    components plus the number and norm of substantial host effects. The parasite target still goes
    through the identity layer; symbols and the historical numeric target are not guessed as ME49 ids.
    """
    root = os.path.join(base, "datasets", "DNA", "CRISPR_screen", "37827122")
    path = os.path.join(root, "hosttx_effectors_DE_host_genes.csv")
    mapping_path = os.path.join(root, "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx")
    if not os.path.exists(path) or not os.path.exists(mapping_path):
        log("screens: full dual perturb-seq host signatures absent")
        return pd.DataFrame()
    data = pd.read_csv(path)
    required = {"target", "gene", "avg_log2FC", "p_val_bh"}
    if not required <= set(data.columns):
        log("screens: dual perturb-seq signature table has unexpected columns")
        return pd.DataFrame()
    mapping = pd.read_excel(mapping_path).rename(columns=lambda column: str(column).strip())
    lookup = {}
    for row in mapping.itertuples(index=False):
        name = str(getattr(row, "Gene_Name", "")).strip().upper()
        accession = str(getattr(row, "Gene_ID_Updated", "")).strip()
        if name and accession and accession.lower() != "nan":
            lookup[name] = accession

    def target_gene(value):
        text = str(value).strip()
        accession = lookup.get(text.upper(), text if text.upper().startswith("TG") else "")
        return resolve(accession) if resolve and accession else accession or None

    data["gene_id"] = data.target.map(target_gene)
    data["avg_log2FC"] = pd.to_numeric(data.avg_log2FC, errors="coerce")
    data["p_val_bh"] = pd.to_numeric(data.p_val_bh, errors="coerce")
    data = data.dropna(subset=["gene_id", "gene", "avg_log2FC"])
    if data.empty:
        return pd.DataFrame()
    matrix = data.pivot_table(index="gene_id", columns="gene", values="avg_log2FC",
                              aggfunc="mean", fill_value=0.0)
    from sklearn.decomposition import PCA
    n = int(min(max(components, 1), max(min(matrix.shape) - 1, 1)))
    scores = PCA(n_components=n, random_state=42).fit_transform(matrix.to_numpy(dtype=float))
    out = pd.DataFrame(scores, index=matrix.index,
                       columns=[f"hosttx_signature_pc{i + 1:02d}" for i in range(n)])
    out["hosttx_signature_norm"] = np.sqrt(np.square(matrix.to_numpy(dtype=float)).sum(axis=1))
    substantial = data[(data.p_val_bh <= 0.05) & (data.avg_log2FC.abs() >= 0.5)]
    out["hosttx_signature_n_de"] = substantial.groupby("gene_id").gene.nunique().reindex(
        out.index, fill_value=0).astype(float)
    log(f"screens: full host-response signatures -> {len(out)} effectors x {matrix.shape[1]:,} "
        f"host genes, represented by {n} PCs")
    return out


def proteomics(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Median log2 iBAQ across replicates, from the Pru IP experiments.

    Coverage is small and it is enrichment, not a deep proteome -- see the module docstring. Reported as
    a measured abundance where present and left missing elsewhere, never imputed to zero: absent from a
    594-protein IP says almost nothing about a gene's true abundance.
    """
    global _RESOLVE
    _RESOLVE = resolve
    P = os.path.join(base, "toxo_stage_atlas", "data", "proteomics")
    frames = []
    for fn, sheets in (("PXD065585_supp_S1_Pru_proteome_and_IP_iBAQ.xlsx",
                        ["Pru_Rep1", "Pru_Rep2", "Pru_Rep3"]),
                       ("PXD043808_supp_MOESM5_Pru_proteome_and_AP2XI2_IP_abundance.xlsx",
                        ["Pru_1", "Pru_2", "Pru_3"])):
        for sh in sheets:
            d = _read(os.path.join(P, fn), sheet_name=sh)
            if d is None or "iBAQ" not in d.columns:
                continue
            v = pd.to_numeric(d["iBAQ"], errors="coerce")
            frames.append(pd.DataFrame({"gene_id": _acc(d["Protein_ID"]),
                                        "ibaq": np.log2(v.where(v > 0))}).dropna())
    if not frames:
        log("proteomics: no iBAQ tables found")
        return pd.DataFrame()
    t = pd.concat(frames).groupby("gene_id").ibaq.median().to_frame("protein_ibaq_log2")
    log(f"proteomics: median log2 iBAQ for {len(t):,} proteins "
        f"(Pru IP experiments -- enrichment, not a deep proteome)")
    return t
