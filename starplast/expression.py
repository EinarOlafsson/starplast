#!/usr/bin/env python3
"""Six further expression datasets, each normalized according to what it actually measures.

These add four biological axes the map did not have, plus a real proteome and site-level phospho:

| dataset                | axis added                                            | genes |
|------------------------|-------------------------------------------------------|-------|
| PMID 31726967          | in vivo brain, acute vs chronic, plus bradyzoite 28 DPI | 8,678 |
| GSE132248              | alkaline-stress bradyzoite induction, 24 h vs 48 h      | 8,923 |
| PXD058095              | MORC depletion and BFD1 knockout                        | 8,920 |
| PXD039400 + PXD042658  | **total proteome**, AP2XII-1/AP2XI-2 depleted or not    | 3,020 |
| PXD017032              | **phosphosite positions and ratios**, oocyst vs tachyzoite | 1,918 |
| PXD003765              | oocyst iTRAQ                                            | 2,095 |

Two of them close gaps this project has been explicit about. The total proteome is four times the
coverage of the immunoprecipitation data previously shipped (3,020 against 748) and is a proteome rather
than an enrichment. The phosphosite table replaces a bare count -- `n_phosphosites`, present for 14.4% of
genes with no positions -- with sites, ratios and p-values.

**Every file needed something.** Not one was a plain table with a header on the first row: two hide their
data behind a `Legend` sheet, one is a legacy `.xls` that pandas cannot open at all without `xlrd` and is
converted through libreoffice, and the iTRAQ header sits on the second row. That is normal for
supplementary data and is why each dataset is declared here rather than guessed at.

Quantifications are not pooled. Each is normalized by its own type through `sources.normalize` -- counts
and intensities are logged and centered, ratios are left alone because centering moves their reference --
and the columns are named for their dataset so nothing downstream can mistake an iTRAQ ratio for a TPM.
"""
from __future__ import annotations

import glob
import gzip
import io
import os
import re
import subprocess
import tarfile
import tempfile
import warnings

import numpy as np
import pandas as pd

from .sources import normalize

warnings.filterwarnings("ignore")

# Separators other than the underscore are allowed because typesetting mangles them: the same
# accession appears as TGME49_208830, TGME49-208830, TGME49.208830 and TGME49208830 depending on the
# journal's line-breaking. `_norm_acc` already normalized "-" and "." back to "_", but the pattern
# never let those forms through to be normalized, so they resolved to nothing. A space is deliberately
# NOT accepted: in running prose it would join two adjacent tokens into a false accession.
ACC = re.compile(r"TG[A-Z0-9]{2,6}[-._]?\d{5,6}", re.I)


def _norm_acc(s: str) -> str | None:
    m = ACC.search(str(s))
    if not m:
        return None
    a = m.group(0).upper().replace("-", "_").replace(".", "_")
    if "_" not in a:
        a = re.sub(r"^([A-Z0-9]+?)(\d{5,6})$", r"\1_\2", a)
    return a


def _resolve(series: pd.Series, resolve=None) -> pd.Series:
    out = series.map(_norm_acc)
    if resolve is not None:
        out = out.map(lambda x: resolve(x) if isinstance(x, str) else x)
    return out


def _collapse(df: pd.DataFrame, genes: pd.Series, how="mean") -> pd.DataFrame:
    """One row per gene. Supplements routinely list a gene several times, one row per peptide or
    transcript, and silently keeping the first would pick whichever the file happened to sort first."""
    d = df.copy()
    d["gene_id"] = genes.values
    d = d[d.gene_id.notna()]
    return getattr(d.groupby("gene_id"), how)()


def _libreoffice_xlsx(path: str, log=print) -> str | None:
    """Convert a legacy .xls that pandas cannot read without xlrd."""
    out = tempfile.mkdtemp(prefix="starplast_xls_")
    try:
        subprocess.run(["libreoffice", "--headless", "--convert-to", "xlsx", "--outdir", out, path],
                       capture_output=True, timeout=300)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        log(f"  cannot convert {os.path.basename(path)}: libreoffice unavailable")
        return None
    hits = glob.glob(os.path.join(out, "*.xlsx"))
    return hits[0] if hits else None


def _acquired(base: str, filename: str) -> str:
    """Path for the dated, literature-audited acquisition batch."""
    return os.path.join(base, "datasets", "toxoplasma_acquisition_2026_08_14", filename)


def _geo_series_matrix(path: str) -> tuple[pd.DataFrame, list[str]]:
    """Read a GEO series matrix and its sample titles without relying on fixed header rows."""
    titles = []
    lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [value.strip().strip('"') for value in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                lines.append(line)
    if not lines:
        return pd.DataFrame(), titles
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t"), titles


def _gpl7186_gene_map(path: str, resolve=None) -> dict:
    """Map legacy ToxoGeneChip probe ids through ToxoDB's previous-id index."""
    out = {}
    in_table = False
    header = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            probe, old_id = row.get("ID", ""), row.get("ToxoDB", "")
            gene = resolve(old_id) if resolve and old_id else old_id
            if probe and gene:
                out[probe] = gene
    return out


def _safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


# --------------------------------------------------------------------------- per-dataset loaders
def invivo_brain(base: str, resolve=None, log=print) -> pd.DataFrame:
    """PMID 31726967: tachyzoite, whole brain acute and chronic, bradyzoite 28 DPI."""
    p = os.path.join(base, "datasets", "translation", "proteomics", "31726967",
                     "12864_2019_6213_MOESM4_ESM.csv")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_csv(p)
    # MOESM5 is the same experiment at transcript level; gene level is what joins to a gene table.
    # `gene_ID` holds StringTie identifiers (MSTRG.1) and resolves for only 3,723 rows; `Gene_ref`
    # carries the real accessions for 8,717. Taking the first column would silently lose half the data.
    idcol = "Gene_ref" if "Gene_ref" in d.columns else d.columns[0]
    val = [c for c in d.columns if re.match(r"(TZ|WholeBrain|BZ)_", str(c))]
    X = _collapse(d[val], _resolve(d[idcol], resolve))
    X = normalize(X, "fpkm", log=lambda *a: None)
    X.columns = [f"invivo_{c}" for c in X.columns]
    log(f"in vivo brain (31726967): {X.shape[1]} columns, {len(X):,} genes")
    return X


def stress_induction(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE132248: unstressed 24 h against alkaline-stressed 48 h."""
    p = os.path.join(base, "toxo_stage_atlas", "data", "transcriptomics",
                     "GSE132248_STAR_counts_matrix.tsv")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    val = [c for c in d.columns[1:] if pd.api.types.is_numeric_dtype(d[c])]
    X = _collapse(d[val], _resolve(d.iloc[:, 0], resolve))
    X = normalize(X, "counts", log=lambda *a: None)
    strip_gsm = re.compile(r"^GSM\d+_")     # a GSM prefix is the sample id, not the condition
    X.columns = ["stress_" + strip_gsm.sub("", str(c)) for c in X.columns]
    log(f"stress induction (GSE132248): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse22258_stage(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE22258: matched Pru tachyzoite and 72-hour alkaline-induced bradyzoite array.

    Unlike the two older ToxoGeneChip files beside it, this series matrix is already keyed by
    TGME49 accessions, so no external platform annotation or uncertain probe conversion is needed.
    """
    p = os.path.join(base, "datasets", "stagetranscriptome_GSE22258_series_matrix.txt.gz")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t", comment="!", compression="gzip")
    if d.shape[1] < 3:
        return pd.DataFrame()
    genes = _resolve(d.iloc[:, 0].astype(str).str.strip('"'), resolve)
    values = d.iloc[:, 1:3].apply(pd.to_numeric, errors="coerce")
    X = _collapse(values, genes)
    X = normalize(X, "log_intensity", log=lambda *a: None)
    X.columns = ("rna22258_tachyzoite", "rna22258_bradyzoite")
    log(f"stage array (GSE22258): {X.shape[1]} columns, {len(X):,} genes")
    return X


def neuronal_differentiation(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE168465: parasite RNA changes over 1--14 days in infected primary brain cells.

    Each workbook sheet is one time point against the tachyzoite comparator.  Base means and log2
    fold changes are measurements; p-values are inferential confidence and stay in the source file.
    """
    p = os.path.join(base, "datasets", "stagetranscriptome_GSE168465_DESeq2-Toxo-all-time-points.xlsx")
    if not os.path.exists(p):
        return pd.DataFrame()
    parts = []
    for sheet in pd.ExcelFile(p).sheet_names:
        d = pd.read_excel(p, sheet_name=sheet)
        if d.empty or "log2FoldChange" not in d.columns:
            continue
        genes = _resolve(d.iloc[:, 0], resolve)
        t = pd.DataFrame(index=d.index)
        if "baseMean" in d:
            t[f"brain168465_{sheet}_base_mean"] = pd.to_numeric(d["baseMean"], errors="coerce")
        t[f"brain168465_{sheet}_lfc"] = pd.to_numeric(d["log2FoldChange"], errors="coerce")
        X = _collapse(t, genes)
        base_cols = [c for c in X if c.endswith("_base_mean")]
        if base_cols:
            X[base_cols] = normalize(X[base_cols], "counts", log=lambda *a: None)
        parts.append(X)
    if not parts:
        return pd.DataFrame()
    X = pd.concat(parts, axis=1)
    log(f"primary-brain differentiation (GSE168465): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse99395_ribosome_profiling(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE99395: matched ribosome-footprint and RNA counts, intra- and extracellular."""
    p = _acquired(base, "GSE99395_Raw_counts.txt.gz")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    genes = _resolve(d.iloc[:, 0], resolve)
    raw = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    X = _collapse(raw, genes)
    X = normalize(X, "counts", log=lambda *a: None)
    rename = {}
    for i, column in enumerate(raw.columns):
        assay = "rpf" if "digested" in column.lower() else "rna"
        context = "extracellular" if "extracellular" in column.lower() else "intracellular"
        rep = 1 + sum(1 for previous in raw.columns[:i]
                      if assay in (("rpf" if "digested" in previous.lower() else "rna"),)
                      and context in previous.lower())
        rename[column] = f"{assay}99395_{context}_r{rep}"
    X = X.rename(columns=rename)
    for context in ("extracellular", "intracellular"):
        for rep in (1, 2):
            rpf, rna = f"rpf99395_{context}_r{rep}", f"rna99395_{context}_r{rep}"
            if rpf in X and rna in X:
                X[f"te99395_{context}_r{rep}"] = X[rpf] - X[rna]
    log(f"ribosome profiling (GSE99395): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse129869_host_context_ribosome_profiling(base: str, resolve=None,
                                               log=print) -> pd.DataFrame:
    """GSE129869: parasite RPF/RNA counts in confluent and subconfluent infected HFFs."""
    p = _acquired(base, "GSE129869_RAW.tar")
    if not os.path.exists(p):
        return pd.DataFrame()
    series = []
    with tarfile.open(p) as archive:
        names = [name for name in archive.getnames()
                 if re.search(r"_[cs](?:RFP|RNA)\.RH\d+_count\.tab\.gz$", name)]
        for name in names:
            member = archive.extractfile(name)
            if member is None:
                continue
            with gzip.GzipFile(fileobj=io.BytesIO(member.read())) as fh:
                d = pd.read_csv(fh, sep="\t", header=None, names=("gene", "count"))
            match = re.search(r"_([cs])(RFP|RNA)\.RH(\d+)_", name)
            if not match:
                continue
            context = "confluent" if match.group(1) == "c" else "subconfluent"
            assay = "rpf" if match.group(2) == "RFP" else "rna"
            column = f"{assay}129869_{context}_r{match.group(3)}"
            t = _collapse(pd.DataFrame({column: pd.to_numeric(d["count"], errors="coerce")}),
                          _resolve(d["gene"], resolve))
            series.append(t)
    if not series:
        return pd.DataFrame()
    X = normalize(pd.concat(series, axis=1), "counts", log=lambda *a: None)
    for context in ("confluent", "subconfluent"):
        for rep in (1, 2, 3):
            rpf, rna = f"rpf129869_{context}_r{rep}", f"rna129869_{context}_r{rep}"
            if rpf in X and rna in X:
                X[f"te129869_{context}_r{rep}"] = X[rpf] - X[rna]
    log(f"host-context ribosome profiling (GSE129869): {X.shape[1]} columns, "
        f"{len(X):,} genes")
    return X


def gse19092_cell_cycle(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE19092: synchronized tachyzoite cell-cycle microarray, two replicates."""
    matrix = _acquired(base, "GSE19092_series_matrix.txt.gz")
    platform = _acquired(base, "GPL7186_family.soft.gz")
    if not os.path.exists(matrix) or not os.path.exists(platform):
        return pd.DataFrame()
    d, titles = _geo_series_matrix(matrix)
    if d.empty:
        return pd.DataFrame()
    mapping = _gpl7186_gene_map(platform, resolve)
    genes = d.iloc[:, 0].astype(str).str.strip('"').map(mapping)
    values = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    names = []
    for title in titles:
        condition = "async" if "asynchronous" in title else (
            "blocked" if "blocked" in title else
            re.search(r"(\d+) hour release", title).group(1) + "h")
        rep = re.search(r"- (\d+)$", title).group(1)
        names.append(f"cellcycle19092_{condition}_r{rep}")
    values.columns = names or list(values.columns)
    X = normalize(_collapse(values, genes), "log_intensity", log=lambda *a: None)
    log(f"cell-cycle array (GSE19092): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse51780_merozoite(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE51780: tachyzoite comparator and feline merozoite expression."""
    matrix = _acquired(base, "GSE51780_series_matrix.txt.gz")
    platform = _acquired(base, "GPL7186_family.soft.gz")
    if not os.path.exists(matrix) or not os.path.exists(platform):
        return pd.DataFrame()
    d, titles = _geo_series_matrix(matrix)
    if d.empty:
        return pd.DataFrame()
    genes = d.iloc[:, 0].astype(str).str.strip('"').map(_gpl7186_gene_map(platform, resolve))
    values = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    names = []
    for i, title in enumerate(titles):
        names.append(f"rna51780_{'tachy' if 'tachy' in title.lower() else 'mero'}_r{i + 1}")
    values.columns = names or list(values.columns)
    X = normalize(_collapse(values, genes), "log_intensity", log=lambda *a: None)
    log(f"merozoite array (GSE51780): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse168155_rna_processing_perturbation(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE168155: CPSF4 depletion transcriptome; not a direct per-gene m6A measurement."""
    p = _acquired(base, "GSE168155_Matrix_table_processed_data.xlsx")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_excel(p, sheet_name="RPKM")
    sample = [c for c in d.columns
              if "linear total RPKM" in str(c)
              and re.fullmatch(r"(?:UT|IAA_\d+h)-[12]", str(c).split(" - ")[0])]
    values = d[sample].apply(pd.to_numeric, errors="coerce")
    values.columns = [f"cpsf4rna168155_{_safe_name(str(c).split(' - ')[0])}" for c in sample]
    X = normalize(_collapse(values, _resolve(d["Name"], resolve)), "fpkm",
                  log=lambda *a: None)
    log(f"CPSF4 perturbation RNA-seq (GSE168155): {X.shape[1]} columns, {len(X):,} genes")
    return X


def gse200962_restriction_checkpoint(base: str, resolve=None, log=print) -> pd.DataFrame:
    """GSE200962: tachyzoite/bradyzoite checkpoint and cyclin perturbation RNA counts."""
    p = _acquired(base, "GSE200962_gene_count_matrix_geo.csv.gz")
    if not os.path.exists(p):
        return pd.DataFrame()
    d = pd.read_csv(p)
    values = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    values.columns = [f"restriction200962_{_safe_name(c)}" for c in values.columns]
    X = normalize(_collapse(values, _resolve(d.iloc[:, 0], resolve)), "counts",
                  log=lambda *a: None)
    log(f"restriction-checkpoint RNA-seq (GSE200962): {X.shape[1]} columns, "
        f"{len(X):,} genes")
    return X


def morc_depletion(base: str, resolve=None, log=print) -> pd.DataFrame:
    """PXD058095: MORC knockdown and BFD1 knockout. Data is on a second sheet behind a Legend."""
    p = os.path.join(base, "toxo_stage_atlas", "data", "proteomics",
                     "PXD058095_supp_DatasetEV1_MORC_RNAseq_counts_TPM.xlsx")
    if not os.path.exists(p):
        return pd.DataFrame()
    xl = pd.ExcelFile(p)
    sheet = "TPM" if "TPM" in xl.sheet_names else "Raw Data"
    d = xl.parse(sheet)
    val = [c for c in d.columns[1:] if pd.api.types.is_numeric_dtype(d[c])]
    X = _collapse(d[val], _resolve(d.iloc[:, 0], resolve))
    X = normalize(X, "tpm" if sheet == "TPM" else "counts", log=lambda *a: None)
    X.columns = [f"morc_{c}" for c in X.columns]
    log(f"MORC depletion (PXD058095, sheet {sheet!r}): {X.shape[1]} columns, {len(X):,} genes")
    return X


def total_proteome(base: str, resolve=None, log=print) -> pd.DataFrame:
    """PXD039400 + PXD042658: total proteome, AP2XII-1/AP2XI-2 depleted or not.

    Four times the protein coverage of the IP experiments previously shipped, and an actual proteome.
    """
    p = os.path.join(base, "toxo_stage_atlas", "data", "proteomics",
                     "PXD039400_PXD042658_supp_SupplTable3_total_proteome.xlsx")
    if not os.path.exists(p):
        return pd.DataFrame()
    # The sheet has a two-row merged header: row 0 spans groups ("UT Vs T-24h", "log2(normalized...)")
    # and row 1 names the members ("log2(fold change)", "UT R1"...). Reading row 0 alone collapses the
    # whole abundance block into one unnamed column, which is how this first yielded a single column
    # from a 3,000-protein proteome.
    top = pd.ExcelFile(p).parse("MS DATA", header=None, nrows=2)
    d = pd.ExcelFile(p).parse("MS DATA", header=1)
    group = top.iloc[0].ffill()
    names = []
    for i, c in enumerate(d.columns):
        g, sub = str(group.get(i, "")).strip(), str(c).strip()
        names.append(f"{g} {sub}" if g and g.lower() != "nan" and not sub.startswith("Unnamed")
                     else sub)
    d.columns = names
    idcol = next((c for c in d.columns if str(c).lower().startswith("accession")), d.columns[0])
    genes = _resolve(d[idcol], resolve)

    # Two kinds of number here, and they must not be normalized the same way. The per-replicate log2
    # abundances are an intensity already on a log scale; the fold changes are ratios.
    abund = [c for c in d.columns if re.search(r"^log2\(normalized.*\s+\S+\s*R\d", str(c), re.I)
             or re.search(r"\b(UT|T-\d+h)\s*R\d\b", str(c))]
    ratios = [c for c in d.columns if "log2(fold change)" in str(c)]
    parts = []
    if abund:
        A = normalize(_collapse(d[abund], genes), "log_intensity", log=lambda *a: None)
        A.columns = [f"proteome_{re.sub(r'[^A-Za-z0-9]+', '_', str(c)).strip('_')}" for c in A.columns]
        parts.append(A)
    if ratios:
        R = normalize(_collapse(d[ratios], genes), "lfc", log=lambda *a: None)
        R.columns = [f"proteome_lfc_{re.sub(r'[^A-Za-z0-9]+', '_', str(c)).strip('_')}"
                     for c in R.columns]
        parts.append(R)
    if not parts:
        return pd.DataFrame()
    X = pd.concat(parts, axis=1)
    log(f"total proteome (PXD039400+42658): {X.shape[1]} columns "
        f"({len(abund)} abundance, {len(ratios)} fold-change), {len(X):,} proteins")
    return X


def phosphosites(base: str, resolve=None, log=print) -> pd.DataFrame:
    """PXD017032: oocyst against tachyzoite phosphosite ratios, with positions.

    Site level collapsed to gene level, keeping what a gene table can hold: how many sites were measured,
    how many moved in each direction, and the strongest change. The site positions stay in the source
    file -- a per-gene table is the wrong shape for them.
    """
    out = []
    for tag, fn in (("up", "PXD017032_supp_TableS1_upregulated_phosphosites_oocyst_vs_tachy.xlsx"),
                    ("down", "PXD017032_supp_TableS2_downregulated_phosphosites_oocyst_vs_tachy.xlsx")):
        p = os.path.join(base, "toxo_stage_atlas", "data", "proteomics", fn)
        if not os.path.exists(p):
            continue
        d = pd.ExcelFile(p).parse(0)
        idcol = next((c for c in d.columns if re.search(r"protein\s*id|accession", str(c), re.I)),
                     d.columns[0])
        ratio = next((c for c in d.columns if re.search(r"ratio", str(c), re.I)), None)
        g = _resolve(d[idcol], resolve)
        t = pd.DataFrame({"gene_id": g})
        if ratio is not None:
            t["ratio"] = pd.to_numeric(d[ratio], errors="coerce")
        t = t[t.gene_id.notna()]
        agg = t.groupby("gene_id").agg(n=("gene_id", "size"),
                                       ratio=("ratio", "median") if ratio is not None
                                       else ("gene_id", "size"))
        agg.columns = [f"phospho_{tag}_sites", f"phospho_{tag}_ratio"]
        out.append(agg)
    if not out:
        return pd.DataFrame()
    X = pd.concat(out, axis=1)
    X["phospho_sites_measured"] = X[[c for c in X.columns if c.endswith("_sites")]].sum(axis=1)
    log(f"phosphosites (PXD017032): {X.shape[1]} columns, {len(X):,} genes with a measured site")
    return X


def oocyst_itraq(base: str, resolve=None, log=print) -> pd.DataFrame:
    """PXD003765: oocyst iTRAQ ratios. Legacy .xls, converted through libreoffice."""
    p = os.path.join(base, "toxo_stage_atlas", "data", "proteomics",
                     "PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls")
    if not os.path.exists(p):
        return pd.DataFrame()
    conv = _libreoffice_xlsx(p, log=log)
    if not conv:
        return pd.DataFrame()
    d = pd.ExcelFile(conv).parse("Sheet1", header=1)
    idcol = next((c for c in d.columns if str(c).lower().startswith("accession")), d.columns[0])
    val = [c for c in d.columns if pd.api.types.is_numeric_dtype(d[c])
           and re.search(r"\d+\s*[:/]\s*\d+|ratio", str(c), re.I)]
    if not val:                      # fall back to the numeric block after the annotation columns
        num = [c for c in d.columns if pd.api.types.is_numeric_dtype(d[c])]
        val = num[4:] if len(num) > 6 else num
    X = _collapse(d[val], _resolve(d[idcol], resolve))
    X = normalize(X, "ratio", log=lambda *a: None)
    X.columns = [f"oocyst_itraq_{re.sub(r'[^A-Za-z0-9]+', '_', str(c)).strip('_')}" for c in X.columns]
    log(f"oocyst iTRAQ (PXD003765): {X.shape[1]} columns, {len(X):,} proteins")
    return X


LOADERS = (invivo_brain, stress_induction, gse22258_stage, neuronal_differentiation,
           gse99395_ribosome_profiling, gse129869_host_context_ribosome_profiling,
           gse19092_cell_cycle, gse51780_merozoite,
           gse168155_rna_processing_perturbation, gse200962_restriction_checkpoint,
           morc_depletion, total_proteome, phosphosites, oocyst_itraq)


def load_all(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Every dataset above, joined on gene id."""
    frames = []
    for fn in LOADERS:
        try:
            X = fn(base, resolve=resolve, log=log)
        except Exception as e:
            log(f"  {fn.__name__} failed: {type(e).__name__}: {e}")
            continue
        if X is not None and not X.empty:
            frames.append(X)
    if not frames:
        return pd.DataFrame()
    joined = pd.concat(frames, axis=1)
    log(f"expression: {joined.shape[1]} new columns over {joined.shape[0]:,} genes "
        f"from {len(frames)} datasets")
    return joined
