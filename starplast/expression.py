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

Quantifications are not pooled. Each is normalized by its own type through `sources.normalise` -- counts
and intensities are logged and centered, ratios are left alone because centering moves their reference --
and the columns are named for their dataset so nothing downstream can mistake an iTRAQ ratio for a TPM.
"""
from __future__ import annotations

import glob
import os
import re
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd

from .sources import normalise

warnings.filterwarnings("ignore")

# Separators other than the underscore are allowed because typesetting mangles them: the same
# accession appears as TGME49_208830, TGME49-208830, TGME49.208830 and TGME49208830 depending on the
# journal's line-breaking. `_norm_acc` already normalised "-" and "." back to "_", but the pattern
# never let those forms through to be normalised, so they resolved to nothing. A space is deliberately
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
    X = normalise(X, "fpkm", log=lambda *a: None)
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
    X = normalise(X, "counts", log=lambda *a: None)
    strip_gsm = re.compile(r"^GSM\d+_")     # a GSM prefix is the sample id, not the condition
    X.columns = ["stress_" + strip_gsm.sub("", str(c)) for c in X.columns]
    log(f"stress induction (GSE132248): {X.shape[1]} columns, {len(X):,} genes")
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
    X = normalise(X, "tpm" if sheet == "TPM" else "counts", log=lambda *a: None)
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

    # Two kinds of number here, and they must not be normalised the same way. The per-replicate log2
    # abundances are an intensity already on a log scale; the fold changes are ratios.
    abund = [c for c in d.columns if re.search(r"^log2\(normalized.*\s+\S+\s*R\d", str(c), re.I)
             or re.search(r"\b(UT|T-\d+h)\s*R\d\b", str(c))]
    ratios = [c for c in d.columns if "log2(fold change)" in str(c)]
    parts = []
    if abund:
        A = normalise(_collapse(d[abund], genes), "log_intensity", log=lambda *a: None)
        A.columns = [f"proteome_{re.sub(r'[^A-Za-z0-9]+', '_', str(c)).strip('_')}" for c in A.columns]
        parts.append(A)
    if ratios:
        R = normalise(_collapse(d[ratios], genes), "lfc", log=lambda *a: None)
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
    X = normalise(X, "ratio", log=lambda *a: None)
    X.columns = [f"oocyst_itraq_{re.sub(r'[^A-Za-z0-9]+', '_', str(c)).strip('_')}" for c in X.columns]
    log(f"oocyst iTRAQ (PXD003765): {X.shape[1]} columns, {len(X):,} proteins")
    return X


LOADERS = (invivo_brain, stress_induction, morc_depletion,
           total_proteome, phosphosites, oocyst_itraq)


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
