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
#: Either strain namespace, for tables that mix them or use type I throughout.
ACC_ANY = re.compile(r"(TGME49_\d{6}|TGGT1_\d{6})")


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


#: The differentiation reporter screen, GSE132237. Two libraries against putative nucleic-acid
#: binding proteins in a strain carrying a bradyzoite reporter; guides enriched in parasites that
#: FAILED to switch the reporter on mark genes the switch needs.
DIFFERENTIATION_SCREEN = "GSE132237_RAW.tar"

#: The comparison the deposit exists to support: reporter-positive parasites against the bulk
#: population at the same timepoint. Everything else in the archive -- the input library and the
#: passages -- measures ordinary growth, which this map already has from the fitness screens.
DIFFERENTIATION_ARMS = (("L1 mNG+ 10d", "L1 bulk brady 10d"), ("L2 mNG+ 10d", "L2 bulk brady 10d"))


#: The mineCETSA Euclidean-distance score: how far a protein's melting curve moves when calcium is
#: added. One number per protein, from the paper's own fit -- not recomputed from the ten temperature
#: points beside it, which are also published.
THERMAL_SHIFT = "PMC9436416_mineCETSA_ED_score.tsv"


def thermal_shift(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Calcium-induced thermal-shift score per protein.

    Verified against the calcium sensors: CAM1 and CAM2 sit at the 98th percentile and CAM3 at the
    83rd. A protein whose melting curve does not move when calcium is added is not calcium-binding,
    so those three had to be near the top or the column would be measuring something else.

    Note for anyone comparing this against the paper: its headline conclusion is about PP1, and PP1
    is unremarkable HERE. That claim comes from the zaprinast time course in the same paper, which is
    a different experiment; this is the calcium mineCETSA sheet.
    """
    path = _find(base, THERMAL_SHIFT)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    if d.shape[1] < 2:
        return pd.DataFrame()
    gid = d[d.columns[0]].astype(str)
    gid = gid.map(lambda g: resolve(g) or g) if resolve is not None else _acc(gid)
    out = pd.DataFrame({"cetsa_calcium_ed_score":
                        pd.to_numeric(d[d.columns[1]], errors="coerce").to_numpy()})
    out["gene_id"] = gid.to_numpy()
    out = out.dropna(subset=["gene_id"]).groupby("gene_id").max()
    log(f"thermal shift (PMC9436416): {len(out):,} proteins")
    return out


#: The cyst wall interactome table. Its bait columns are named for the protein pulled down, and
#: everything that is not one of these nine bookkeeping columns is a bait.
CYST_WALL = "PMC7002340_cyst_wall_interactome.tsv"
CYST_WALL_META = ("#", "Visible?", "Starred?", "Identified Proteins (248/260)",
                  "Accession Number", "Alternate ID", "Molecular Weight",
                  "Protein Grouping Ambiguity", "Taxonomy")


def cyst_wall_interactome(base: str, log=print, resolve=None) -> pd.DataFrame:
    """What co-purifies with the cyst wall, across every bait in the study.

    Two numbers per protein: the strongest normalised spectral count it reached with any bait, and
    how many baits saw it at all. The second is the more honest of the two -- a protein found by
    thirteen independent pulldowns is in the cyst wall in a way that a single strong hit is not.

    Most rows are human: the pulldowns were done on infected cultures and the table lists everything
    identified. Only rows naming a Toxoplasma accession are kept, which is 57 of 265.
    """
    path = _find(base, CYST_WALL)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t", low_memory=False)
    if "Accession Number" not in d.columns:
        return pd.DataFrame()
    baits = [c for c in d.columns if c not in CYST_WALL_META]
    if not baits:
        return pd.DataFrame()
    gene = d["Accession Number"].astype(str).str.extract(ACC_ANY, expand=False)
    if resolve is not None:
        gene = gene.map(lambda g: (resolve(g) or g) if isinstance(g, str) else g)
    values = d[baits].apply(pd.to_numeric, errors="coerce")
    out = pd.DataFrame({"cyst_wall_max_spectral": values.max(axis=1),
                        "cyst_wall_n_baits": (values > 0).sum(axis=1).astype(float)})
    out["gene_id"] = gene.to_numpy()
    out = out.dropna(subset=["gene_id"]).groupby("gene_id").max()
    log(f"cyst wall interactome (PMC7002340): {len(out):,} Toxoplasma proteins")
    return out


#: The screen's own score column, one row per gene, from the paper's Data Sheet 1.
OXIDATIVE_SCREEN = "PMC8216390_screening_score.tsv"


def oxidative_stress_screen(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Per-gene score from the genome-wide CRISPR screen for defence against oxidative stress.

    Read from the authors' own `Screening score` sheet rather than recomputed from the guide counts
    beside it: they published the score, so it is theirs to define.

    Negative is required -- a gene whose disruption is depleted under oxidative challenge. The check
    that says the sign is right is catalase at -6.15, which is the extreme of the whole screen and is
    the enzyme that disposes of hydrogen peroxide.
    """
    path = _find(base, OXIDATIVE_SCREEN)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    if d.shape[1] < 2:
        return pd.DataFrame()
    gid = _acc(d[d.columns[0]]) if resolve is None else pd.Series(
        [resolve(g) or g for g in d[d.columns[0]].astype(str)])
    out = pd.DataFrame({"oxidative_stress_screen_score":
                        pd.to_numeric(d[d.columns[1]], errors="coerce").to_numpy()})
    out["gene_id"] = gid.to_numpy()
    out = out.dropna(subset=["gene_id"]).groupby("gene_id").mean()
    log(f"oxidative stress screen (PMC8216390): {len(out):,} genes")
    return out

def _gse132237_counts(base: str) -> dict:
    """Guide counts per gene for every sample in the differentiation-screen archive.

    Shared by the two questions this deposit answers -- what happens on differentiation, and what
    happens over eight passages of ordinary growth. Which member is which sample comes from the
    series matrix; the file names carry the GSM accession and the submitter's own suffix, so
    matching a title against a name is a guess about punctuation.
    """
    import gzip
    import tarfile

    path = _find(base, DIFFERENTIATION_SCREEN)
    matrix = _find(base, DIFFERENTIATION_SCREEN.replace("_RAW.tar", "_series_matrix.txt.gz"))
    if not os.path.exists(path) or not os.path.exists(matrix):
        return {}
    with gzip.open(matrix, "rt", errors="replace") as fh:
        text = fh.read(400_000)
    names = re.search(r"!Sample_title\t(.*)", text)
    accessions = re.search(r"!Sample_geo_accession\t(.*)", text)
    if not (names and accessions):
        return {}
    by_gsm = dict(zip([v.strip().strip('"') for v in accessions.group(1).split("\t")],
                      [v.strip().strip('"') for v in names.group(1).split("\t")]))
    out = {}
    with tarfile.open(path) as archive:
        for member in archive.getnames():
            gsm = re.match(r"(GSM\d+)_", os.path.basename(member))
            title = by_gsm.get(gsm.group(1)) if gsm else None
            if not title:
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            blob = handle.read()
            body = gzip.decompress(blob).decode("utf8", "replace") if member.endswith(".gz") \
                else blob.decode("utf8", "replace")
            per_gene = {}
            for line in body.splitlines():
                parts = line.split("\t")
                if len(parts) < 2 or "_" not in parts[0]:
                    continue
                gene = parts[0].rsplit("-", 1)[0]
                try:
                    per_gene[gene] = per_gene.get(gene, 0.0) + float(parts[1])
                except ValueError:
                    continue
            out[title] = pd.Series(per_gene, dtype=float)
    return out


def _cpm(series: pd.Series) -> pd.Series:
    """Counts as a share of their own library. Guards the empty library rather than dividing by it."""
    total = series.sum()
    return series / total * 1e6 if total else series


def _ratio_arms(counts: dict, arms) -> list:
    """log2 of one arm over another, per gene, for each pair that is present."""
    out = []
    for numerator, denominator in arms:
        a, b = counts.get(numerator), counts.get(denominator)
        if a is None or b is None:
            continue
        shared = a.index.intersection(b.index)
        if not len(shared):
            continue
        out.append(np.log2((_cpm(a)[shared] + 1.0) / (_cpm(b)[shared] + 1.0)))
    return out


def _finish(ratios: list, resolve, column: str, log, what: str) -> pd.DataFrame:
    """Average the arms, resolve accessions, drop duplicates, name the column."""
    out = pd.concat(ratios, axis=1).mean(axis=1).to_frame(column)
    out.index = _acc(pd.Series(out.index)) if resolve is None else \
        pd.Index([resolve(i) or i for i in out.index])
    out = out[~out.index.duplicated()]
    log(f"{what} (GSE132237): {len(out):,} genes")
    return out


#: Eight passages of ordinary tachyzoite growth against the input library. The reporter strain is not
#: the type I RH line the genome-wide screens use, which is what makes this a SECOND background
#: rather than a repeat of them.
PASSAGE_ARMS = (("L1 p8", "L1 input library"), ("L2 p8", "L2 input library"))


def second_background_fitness(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Fitness cost of disrupting each gene over eight passages in the reporter strain.

    A targeted library of nucleic-acid binding genes, so it covers 262 of them and not the genome.
    That is its honest extent, and for those genes it answers something the genome-wide RH screens
    cannot: whether essentiality holds in another background.

    Verified by agreeing with them where it should. Against `fit_invitro_hff` it comes out at
    rho = +0.62 over 130 shared genes -- close enough that the direction and the join are right, far
    enough that the column is not a copy of one already here.
    """
    counts = _gse132237_counts(base)
    ratios = _ratio_arms(counts, PASSAGE_ARMS) if counts else []
    if not ratios:
        log("second-background fitness: no matching passage arms in the archive")
        return pd.DataFrame()
    return _finish(ratios, resolve, "crispr_reporter_strain_p8_log2", log,
                   "second-background fitness")


def differentiation_screen(base: str, log=print, resolve=None) -> pd.DataFrame:
    """Per-gene differentiation phenotype from the reporter screen, as a log2 ratio.

    Guide counts are summed per gene before the ratio is taken, which is the standard readout and
    also the only honest one at this depth: a single guide's count is noisy enough that a per-guide
    ratio averaged afterwards is dominated by whichever guide happened to be sampled least.

    Counts are scaled to a common library size first. Two arms are averaged. Genes seen in neither
    arm are absent rather than zero -- a gene no guide covered is not a gene with no phenotype.

    This is COMPUTED here rather than read from the authors' table, because the deposit publishes
    counts and not the ratio. It is the comparison their design names, and the column is called a
    ratio rather than a phenotype so nobody mistakes it for something they reported.
    """
    counts = _gse132237_counts(base)
    ratios = _ratio_arms(counts, DIFFERENTIATION_ARMS) if counts else []
    if not ratios:
        log("differentiation screen: no matching arms in the archive")
        return pd.DataFrame()
    return _finish(ratios, resolve, "diff_reporter_log2_mNG_over_bulk", log,
                   "differentiation screen")
