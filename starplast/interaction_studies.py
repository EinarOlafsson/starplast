#!/usr/bin/env python3
"""Parsing published proximity-labelling and pulldown supplements into gene lists.

97 studies with a tagged *Toxoplasma* protein were catalogued from the abstract corpus and their
supplementary files downloaded. This module turns the machine-readable ones into a table of
`(study, gene)` membership plus a guess at which protein was tagged.

What it can and cannot do, measured rather than assumed:

* **33 of 97 studies carry accessions in a spreadsheet.** The rest publish their hit lists in PDF or DOCX,
  or by gene name only. Those are recorded with `parsed = False` rather than silently omitted.
* **`TGGT1_` accessions are more common than `TGME49_` here** (40 files against 39), so everything goes
  through the identity layer. Parsing on `TGME49_` alone would lose more than half of what is available.
* **A supplement is not a hit list.** These files mix hits with controls, background, primers and
  full-proteome backgrounds, and the sheet that means "enriched" is named differently in every paper.
  This module therefore reports *membership* — this gene appears in this study's supplement — and does
  not claim enrichment. Turning membership into interaction edges needs per-paper curation, and
  pretending otherwise would manufacture thousands of false interactions.

The bait guess comes from the paper title matched against ToxoDB symbols, which works for the many
papers titled after the protein they tagged ("...IMC29...", "...PDI8...") and fails quietly otherwise.
"""
from __future__ import annotations

import glob
import json
import os
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ACC = re.compile(r"TG[A-Z0-9]{2,6}_\d{5,6}", re.I)
MAX_SHEETS = 8
MAX_ROWS = 20000


def _read_any(path: str) -> str:
    """Flatten a spreadsheet or delimited file to text. Returns '' for formats we cannot read."""
    low = path.lower()
    try:
        if low.endswith((".xlsx", ".xls")):
            xl = pd.ExcelFile(path)
            return "".join(
                xl.parse(sh, header=None, nrows=MAX_ROWS).astype(str).to_string()
                for sh in xl.sheet_names[:MAX_SHEETS])
        if low.endswith((".csv", ".tsv", ".txt")):
            return pd.read_csv(path, sep=None, engine="python", header=None,
                               nrows=MAX_ROWS, dtype=str).astype(str).to_string()
    except Exception:
        return ""
    return ""


def parse_studies(root: str, resolve=None, log=print) -> pd.DataFrame:
    """One row per (study, gene) found in a study's supplementary tables."""
    rows, seen_studies = [], []
    for meta_path in sorted(glob.glob(os.path.join(root, "*", "*", "META.json"))):
        d = os.path.dirname(meta_path)
        meta = json.load(open(meta_path))
        pmid = str(meta.get("pmid") or os.path.basename(d))
        method = meta.get("method", "")
        files = [f for f in sorted(os.listdir(d)) if f != "META.json"]
        tables = [f for f in files if f.lower().endswith((".xlsx", ".xls", ".csv", ".tsv", ".txt"))]

        found = {}
        for f in tables:
            txt = _read_any(os.path.join(d, f))
            if not txt:
                continue
            for a in set(ACC.findall(txt)):
                g = resolve(a) if resolve else (a.upper() if a.upper().startswith("TGME49") else None)
                if g:
                    found.setdefault(g, set()).add(f)

        seen_studies.append({
            "pmid": pmid, "pmcid": meta.get("pmcid", ""), "year": meta.get("year", ""),
            "method": "BioID" if method == "proximity" else "IPMS",
            "title": meta.get("title", ""), "n_files": len(files), "n_tables": len(tables),
            "n_genes": len(found), "parsed": bool(found),
            "unreadable_formats": ";".join(sorted({os.path.splitext(f)[1].lower()
                                                   for f in files if f not in tables})),
        })
        for g, src in found.items():
            rows.append({"pmid": pmid, "gene_id": g,
                         "method": "BioID" if method == "proximity" else "IPMS",
                         "n_files": len(src), "files": ";".join(sorted(src))})

    members = pd.DataFrame(rows)
    studies = pd.DataFrame(seen_studies)
    if not studies.empty:
        ok = int(studies.parsed.sum())
        log(f"interaction studies: {len(studies)} catalogued, {ok} parsed "
            f"({len(studies) - ok} publish hit lists only as PDF/DOCX or by gene name)")
        if not members.empty:
            log(f"  {len(members):,} (study, gene) memberships over "
                f"{members.gene_id.nunique():,} distinct genes")
    return members, studies


def guess_baits(studies: pd.DataFrame, symbol_to_gene: dict, log=print) -> pd.DataFrame:
    """Guess the tagged protein from the paper title, via ToxoDB symbols.

    Titles like "IMC29 Plays an Important Role..." or "...disulfide isomerase PDI8..." name the bait
    directly. This is a heuristic and is labelled as one: `bait_confidence` is `title` when a symbol was
    found and empty otherwise, so nothing downstream can treat a guess as an annotation.
    """
    if studies.empty:
        return studies
    out = studies.copy()
    baits, conf = [], []
    for t in out.title.astype(str):
        hit = None
        for tok in re.findall(r"[A-Za-z][A-Za-z0-9-]{2,}", t):
            g = symbol_to_gene.get(re.sub(r"[-\s]", "", tok).upper())
            if g:
                hit = g
                break
        baits.append(hit)
        conf.append("title" if hit else "")
    out["bait_gene"] = baits
    out["bait_confidence"] = conf
    log(f"  bait identified from the title for {int(out.bait_gene.notna().sum())} of {len(out)} studies")
    return out


def study_gene_counts(members: pd.DataFrame, nodes_index) -> pd.Series:
    """How many interaction studies name each gene -- a per-gene node column."""
    if members.empty:
        return pd.Series(0, index=nodes_index, dtype=int)
    c = members.groupby("gene_id").pmid.nunique()
    return pd.Series(nodes_index).map(c).fillna(0).astype(int).values
