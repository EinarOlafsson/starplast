#!/usr/bin/env python3
"""Downloading expression data from its source repositories, and normalizing it comparably.

Two repositories, two very different situations, and the difference decides what can be automated.

**GEO deposits processed matrices.** `GSE108740_FPKM.xlsx`, `GSE206344_Normalized_data...xlsx` — the
per-gene quantification is right there, so transcriptomics can be fetched and normalized end to end.

**PRIDE deposits raw instrument files.** PXD019729 holds 24 `.raw` files, two search outputs and no
per-protein table at all. Turning those into quantification means running MaxQuant or Proteome Discoverer,
which is a different project. So for proteomics this module fetches what PRIDE labels `RESULT` where it
exists, and otherwise says plainly that only raw was deposited and the numbers must come from the paper's
supplementary tables. Pretending otherwise would be the useful-looking lie here.

## Why normalization needs the quantification type recorded

Different quantifications are not interchangeable, and the failure is silent in both directions. The
same series demonstrates each. GEO's `GSE108740_FPKM.xlsx` reaches **16,520** and is genuine FPKM: skip
the log and one gene dominates every distance in the matrix. The columns derived from it in this
project's node table, still named `rna108740_*_FPKM`, top out at **9.7** because they were logged
upstream and the name was kept: log them again and real variation compresses into nothing.

Neither error is visible in the output. So the range decides, not the filename, and `infer_quant`
reports what it concluded.

So every dataset carries a `quant` label, and the transform follows from it:

    counts        log2(x + 1), then per-sample median centering
    fpkm / tpm    log2(x + 1), then per-sample median centering
    log_intensity already logged: center only, never log twice
    ratio / lfc   already relative: leave alone, centering would destroy the reference
    ibaq / lfq    log2(x + 1), then per-sample median centering

Cross-dataset comparison uses ranks, because that is the only defensible way to put an FPKM column and an
iBAQ column on one axis.
"""
from __future__ import annotations

import os
import re
import urllib.error
import urllib.request

import numpy as np
import pandas as pd

UA = {"User-Agent": "starplast (research; einar.olafsson@gmail.com)"}

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{gse}/suppl/"
PRIDE_FILES = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}/files"

QUANT_TYPES = ("counts", "fpkm", "tpm", "log_intensity", "ratio", "lfc", "ibaq", "lfq")
# Which transforms apply. `log` = take log2(x+1) first; `centre` = subtract each sample's median.
TRANSFORM = {
    "counts": ("log", "centre"), "fpkm": ("log", "centre"), "tpm": ("log", "centre"),
    "ibaq": ("log", "centre"), "lfq": ("log", "centre"),
    "log_intensity": (None, "centre"),
    "ratio": (None, None), "lfc": (None, None),
}


# --------------------------------------------------------------------------- fetching
def _get(url: str, timeout=300) -> bytes:
    """Fetch one URL, and record what came back.

    Logged at both ends rather than only on failure: "which URL did this table actually come from"
    is the first question asked of any number derived from a download, and the answer has to survive
    the session that made it.
    """
    from .logging_util import get_logger
    log = get_logger(__name__)
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
            data = r.read()
    except Exception as exc:
        # WARNING, not DEBUG. A fetch that fails and is handled quietly is how this project shipped
        # a version-pinned URL that returned 404 for months while reporting "no model available".
        log.warning("fetch failed: %s -- %s: %s", url, type(exc).__name__, exc)
        raise
    log.info("fetched %s (%.1f kB)", url, len(data) / 1024)
    return data


def geo_supplementary(gse: str, out_dir: str, log=print) -> list:
    """Download a GEO series' processed supplementary matrices."""
    stub = gse[:-3]
    try:
        page = _get(GEO_SUPPL.format(stub=stub, gse=gse)).decode("utf8", "replace")
    except (urllib.error.URLError, OSError) as e:
        log(f"{gse}: listing failed ({e})")
        return []
    names = [n for n in re.findall(r'href="([^"]+)"', page)
             if not n.startswith(("/", "http")) and n not in ("../",)]
    os.makedirs(out_dir, exist_ok=True)
    got = []
    for n in names:
        dest = os.path.join(out_dir, n)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            got.append(dest)
            continue
        try:
            data = _get(GEO_SUPPL.format(stub=stub, gse=gse) + n)
        except (urllib.error.URLError, OSError) as e:
            log(f"  {n}: {e}")
            continue
        open(dest, "wb").write(data)
        got.append(dest)
        log(f"  {gse}/{n}  {len(data)/1e6:.1f} MB")
    return got


def pride_files(accession: str, log=print) -> pd.DataFrame:
    """What a PRIDE project actually contains, by category.

    The point is usually to discover that it contains no processed table. Raw files are not silently
    downloaded: they are tens of gigabytes and useless without a search engine.
    """
    try:
        import json
        d = json.loads(_get(PRIDE_FILES.format(acc=accession)))
    except Exception as e:
        log(f"{accession}: {e}")
        return pd.DataFrame()
    files = d if isinstance(d, list) else d.get("_embedded", {}).get("files", [])
    rows = [{"accession": accession,
             "category": (f.get("fileCategory") or {}).get("value", "?"),
             "name": f.get("fileName", ""),
             "size_mb": round((f.get("fileSizeBytes") or 0) / 1e6, 1),
             "url": next((p.get("value") for p in (f.get("publicFileLocations") or [])
                          if "FTP" in str(p.get("name", "")).upper()), None)}
            for f in files]
    t = pd.DataFrame(rows)
    if not t.empty:
        counts = t.category.value_counts().to_dict()
        log(f"{accession}: {len(t)} files {counts}")
        if "RESULT" not in counts:
            log("  no RESULT files -- only raw was deposited; per-protein numbers must come from "
                "the paper's supplementary tables")
    return t


# --------------------------------------------------------------------------- normalizing
def infer_quant(values, name: str = "") -> str:
    """Guess a quantification type when it was not recorded, and say so where it is a guess.

    Deliberately conservative. The one inference worth making is spotting data that is *already* logged
    despite a name claiming otherwise -- the `rna108740_*_FPKM` columns top out at 9.13, and logging them
    a second time would flatten real differences into noise.
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if s.empty:
        return "unknown"
    lo, hi = float(s.min()), float(s.max())
    n = name.lower()
    if "lfc" in n or "log2fc" in n or "l2fc" in n or (lo < -1 and hi < 30):
        return "lfc"
    if "ratio" in n:
        return "ratio"
    if hi <= 32 and lo >= 0:
        # Bounded like a log even where the name claims otherwise. This is the case worth catching: the
        # node table's `rna108740_*_FPKM` columns top out at 9.13 and were logged upstream, while the
        # GEO source of the same series reaches 16,520 and is genuine FPKM. Trust the range, not the name.
        return "log_intensity"
    if "ibaq" in n:
        return "ibaq"
    if "tpm" in n:
        return "tpm"
    if float((s % 1 == 0).mean()) == 1.0 and hi > 32:
        return "counts"
    # Wide, positive, non-integer: a linear intensity that still needs logging. Returning
    # "log_intensity" here -- as an earlier version did for anything non-integer -- would skip the log
    # on real FPKM and leave one gene at 16,520 dominating every distance in the matrix.
    return "fpkm"


def normalize(df: pd.DataFrame, quant: str, log=print) -> pd.DataFrame:
    """Standardise one dataset's numeric matrix according to its quantification type."""
    if quant not in TRANSFORM:
        log(f"unknown quantification {quant!r}; leaving values untouched")
        return df
    take_log, centre = TRANSFORM[quant]
    X = df.apply(pd.to_numeric, errors="coerce")
    if take_log == "log":
        X = np.log2(X.clip(lower=0) + 1)
    if centre == "centre":
        # Per-sample median, not per-gene: it corrects for how much material was loaded in each run,
        # which is the systematic difference between columns. Centring per gene would erase the
        # between-gene differences that are the signal.
        X = X.sub(X.median(axis=0), axis=1)
    return X


def rank_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Per-column ranks scaled to [-0.5, 0.5].

    The only defensible way to put an FPKM column and an iBAQ column on one axis: it keeps the ordering,
    which both quantifications support, and discards the units, which are not comparable.
    """
    return df.rank(axis=0, pct=True, na_option="keep") - 0.5


def harmonise(datasets: dict, log=print) -> pd.DataFrame:
    """Normalize several datasets and join them on gene id, ranked onto a common scale.

    `datasets` maps name -> (DataFrame indexed by gene_id, quant type).
    """
    out = []
    for name, (df, quant) in datasets.items():
        q = quant or infer_quant(df.stack(), name)
        norm = rank_normalize(normalize(df, q, log=log))
        norm.columns = [f"{name}__{c}" for c in norm.columns]
        log(f"{name}: {df.shape[1]} columns as {q}, rank-normalized")
        out.append(norm)
    if not out:
        return pd.DataFrame()
    joined = pd.concat(out, axis=1)
    log(f"harmonised: {joined.shape[0]:,} genes x {joined.shape[1]} columns "
        f"from {len(datasets)} datasets")
    return joined
