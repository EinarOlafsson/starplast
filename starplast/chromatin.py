#!/usr/bin/env python3
"""Chromatin signal per gene, summarised from coverage tracks over promoters.

GEO serves these two experiments as bigWig only -- no peak calls, no per-gene table -- so the number
has to be computed here. What is computed is deliberately the dullest possible summary: the mean
coverage over a fixed window around each gene's transcription start, normalised to the genome mean.
No peak calling, no background model, no thresholds. Every one of those is a modelling choice that
would be invented here rather than taken from the authors, and the project's rule is that a number
which looks like a measurement must be one.

## Why the promoter and not the gene body

Accessibility and factor occupancy are promoter phenomena; over a gene body they mostly report how
long the gene is. The window is `PROMOTER` bases either side of the start, taken from the strand, so
a reverse-strand gene's promoter is at its higher coordinate. Getting that backwards is silent -- the
column still looks like data -- which is why the loader is tested on a reverse-strand gene.

## What only the unperturbed arm is used

Both deposits are knockdowns, and both include an untreated arm. Only the untreated arm is read: it
is the measurement of what chromatin looks like in a normal parasite, which is what the slot asks.
The perturbed arms measure what a specific depletion does, which is a different question and would
need its own slot rather than being averaged into this one.
"""
from __future__ import annotations

import io
import os
import re
import tarfile

import numpy as np
import pandas as pd

#: Bases either side of the transcription start. 1 kb is the conventional promoter window and it is
#: wide enough to survive the annotation being a few hundred bases out, which for a parasite genome
#: is a real risk.
PROMOTER = 1000

LOCATION_TABLE = "toxodb_gene_location.tsv"

#: `TGME49_chrXII:2,245,476..2,248,187(-)` -- ToxoDB's display form, commas and all.
LOCATION = re.compile(r"^(?P<chrom>[^:]+):(?P<start>[\d,]+)\.\.(?P<end>[\d,]+)\((?P<strand>[+-])\)")


def gene_windows(base: str, resolve=None) -> pd.DataFrame:
    """Promoter windows per gene: sequence, start, end, taken from the strand."""
    path = os.path.join(base, "starplast", "data", LOCATION_TABLE)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    column = next((c for c in d.columns if "Location" in c or c == "location"), None)
    if column is None:
        return pd.DataFrame()
    rows = []
    for gene, text in zip(d[d.columns[0]].astype(str), d[column].astype(str)):
        match = LOCATION.match(text)
        if not match:
            continue
        start = int(match.group("start").replace(",", ""))
        end = int(match.group("end").replace(",", ""))
        # The transcription start is the LOW coordinate on the forward strand and the HIGH one on
        # the reverse. Taking `start` for both would put half of every promoter set at the far end
        # of the gene, and nothing downstream would say so.
        tss = start if match.group("strand") == "+" else end
        rows.append((resolve(gene) or gene if resolve is not None else gene,
                     match.group("chrom"), max(tss - PROMOTER, 0), tss + PROMOTER))
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows, columns=("gene_id", "chrom", "start", "end"))
    return out.drop_duplicates("gene_id").set_index("gene_id")


def track_means(path: str, windows: pd.DataFrame) -> pd.Series:
    """Mean coverage in each window, as a fraction of the track's genome-wide mean.

    Normalised to the track rather than left raw, because two tracks differ in sequencing depth by
    whatever the submitters happened to load. A ratio to the track's own mean is comparable between
    them and is what makes averaging replicates meaningful.
    """
    import pyBigWig
    bw = pyBigWig.open(path)
    try:
        header = bw.header()
        covered = header.get("nBasesCovered") or 0
        genome_mean = (header.get("sumData") or 0) / covered if covered else 0.0
        chroms = bw.chroms()
        values = []
        for chrom, start, end in zip(windows["chrom"], windows["start"], windows["end"]):
            limit = chroms.get(chrom)
            if not limit or start >= limit:
                values.append(np.nan)
                continue
            got = bw.stats(chrom, int(start), int(min(end, limit)), type="mean")
            values.append(got[0] if got and got[0] is not None else np.nan)
    finally:
        bw.close()
    series = pd.Series(values, index=windows.index, dtype=float)
    return series / genome_mean if genome_mean else series


def readable(path: str) -> bool:
    """Whether this file is a bigWig at all, by its magic number.

    `GSM8524430_UT_2.bw` in GSE277553 begins with eight 0xFF bytes and is not a bigWig. It is the
    same 6,943,536 bytes whether taken from the series tar or fetched from GEO as a sample file, so
    the corruption is in the deposit and not in the download -- which is the only reason it is
    correct to skip it rather than to re-fetch. It is skipped LOUDLY: a replicate silently dropped
    is a mean over fewer samples than the note beside the column claims.
    """
    with open(path, "rb") as fh:
        return fh.read(4).hex() in ("26fc8f88", "888ffc26")


def _from_tar(path: str, pattern: str, windows: pd.DataFrame, log=print) -> list:
    """Every matching bigWig inside an archive, summarised.

    Written out first because bigWig is a random-access format: the reader seeks to an index at the
    end of the file and back, which a tar member stream cannot do. A temporary directory rather than
    somewhere under the dataset tree, so a run that dies halfway leaves nothing behind that a later
    run would mistake for a download.
    """
    import tempfile
    out = []
    with tarfile.open(path) as archive, tempfile.TemporaryDirectory() as scratch:
        for name in sorted(archive.getnames()):
            if not re.search(pattern, os.path.basename(name)):
                continue
            member = archive.extractfile(name)
            if member is None:
                continue
            where = os.path.join(scratch, os.path.basename(name))
            with open(where, "wb") as fh:
                fh.write(member.read())
            if not readable(where):
                log(f"chromatin: {os.path.basename(name)} is not a bigWig, skipped")
                continue
            out.append(track_means(where, windows))
    return out


def chromatin_signals(base: str, resolve=None, log=print) -> pd.DataFrame:
    """ATAC accessibility and HDAC3 occupancy over promoters, from the untreated arms."""
    windows = gene_windows(base, resolve=resolve)
    if windows.empty:
        return pd.DataFrame()
    root = os.path.join(base, "datasets", "quarantine", "2026_08_16_unverified", "Tg")
    out = pd.DataFrame(index=windows.index)

    atac = os.path.join(root, "acetylation", "GSE313048_ATACseq_GCN5b-KD_UT.bw")
    if os.path.exists(atac) and readable(atac):
        out["atac_promoter_ut"] = np.log2(track_means(atac, windows) + 0.01)
        log(f"chromatin: ATAC promoter signal (GSE313048), "
            f"{int(out['atac_promoter_ut'].notna().sum()):,} genes")

    cuttag = os.path.join(root, "chromatin_accessibility", "GSE277553_RAW.tar")
    if os.path.exists(cuttag):
        parts = _from_tar(cuttag, r"_UT_\d+\.bw$", windows, log=log)
        if parts:
            out["cuttag_hdac3_promoter_ut"] = np.log2(
                pd.concat(parts, axis=1).mean(axis=1) + 0.01)
            log(f"chromatin: HDAC3 CUT&TAG promoter signal (GSE277553), {len(parts)} replicates, "
                f"{int(out['cuttag_hdac3_promoter_ut'].notna().sum()):,} genes")
    return out.dropna(axis=1, how="all")
