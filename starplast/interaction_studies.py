#!/usr/bin/env python3
"""Parsing published proximity-labelling and pulldown supplements into gene lists.

97 studies with a tagged *Toxoplasma* protein were catalogued from the abstract corpus and their
supplementary files downloaded. This module turns the machine-readable ones into a table of
`(study, gene)` membership plus a guess at which protein was tagged.

What it can and cannot do, measured rather than assumed:

* **38 of 97 studies yield identifiers.** 32 carry accessions in a spreadsheet; reading the PDFs and
  DOCX adds six more, and recovers 1,124 accessions plus 465 symbols from documents that a spreadsheet
  parser saw as empty. The remaining 59 are recorded with `parsed = False` rather than silently omitted —
  40 of them because Europe PMC held no supplementary files at all, not because parsing failed.
* **Some papers name proteins only by symbol.** One yielded 68 symbols and zero accessions. Symbol
  matching carries the usual guard: a digit-free symbol must appear in upper case, since HOOK, CLAMP,
  CLIP, SPARK and REMIND are all real symbols and all real English words.
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


def _pdf_text(path: str, timeout=180) -> str:
    """PDF to text with columns preserved.

    `-layout` is not optional: without it pdftotext reflows a table into prose and the columns
    interleave, so every extracted row is a blend of two different rows. The output looks fine and the
    data is wrong, which is the worst combination.
    """
    import subprocess
    for args in (["pdftotext", "-layout", "-nopgbrk", path, "-"],
                 ["pdftotext", "-nopgbrk", path, "-"]):
        try:
            r = subprocess.run(args, capture_output=True, timeout=timeout)
            if r.returncode == 0 and r.stdout.strip():
                return r.stdout.decode("utf8", "replace")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return ""


def _docx_text(path: str) -> str:
    """DOCX to text: the format is a zip of XML, so no third-party reader is needed.

    Cell and row boundaries are converted before tags are stripped; otherwise every cell in a row runs
    together and a gene/description table becomes one unusable line.
    """
    import zipfile
    try:
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf8", "replace")
    except (zipfile.BadZipFile, KeyError, OSError):
        return ""
    # A cell's closing paragraph ends the CELL, not the line. Replacing </w:p> with a newline
    # unconditionally fires inside every table cell before the cell boundary is applied, so a
    # gene/description row came out split across two lines -- the same corruption this function exists
    # to prevent, in the other direction. Accession extraction was unaffected, because that scans the
    # whole text, but anything reading the table as rows saw nonsense.
    xml = re.sub(r"</w:p>\s*</w:tc>", "\t", xml)
    xml = xml.replace("</w:tc>", "\t").replace("</w:tr>", "\n").replace("</w:p>", "\n")
    return re.sub(r"<[^>]+>", "", xml)


def _read_any(path: str) -> str:
    """Flatten a document or spreadsheet to text. Returns '' for formats we cannot read."""
    low = path.lower()
    try:
        if low.endswith(".pdf"):
            return _pdf_text(path)
        if low.endswith(".docx"):
            return _docx_text(path)
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
        # PDF and DOCX included: 65 of 97 studies publish their hit list only in those, and
        # 20 of them yield identifiers once the documents are actually read.
        tables = [f for f in files
                  if f.lower().endswith((".xlsx", ".xls", ".csv", ".tsv", ".txt",
                                         ".pdf", ".docx"))]

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


# --------------------------------------------------------------------------- host targets
def host_interactions(base: str, resolve=None, log=print) -> pd.DataFrame:
    """Curated Toxoplasma-protein to host-protein interactions.

    Deliberately the curated table rather than anything mined from the 97 supplements. Those publish
    complete quantification tables -- median 754 genes, largest 7,866 -- so extracting host targets from
    them automatically would invent thousands of interactions. 31 curated pairs with a stated mode,
    confidence and reference are worth more than that, and can be defended one by one.

    The parasite side is given by symbol (ROP16, GRA24), so it goes through the identity layer.
    """
    path = os.path.join(base, "known_host_parasite_interactions.csv")
    if not os.path.exists(path):
        log("no curated host-parasite table found")
        return pd.DataFrame()
    d = pd.read_csv(path)
    d = d[d.organism.astype(str).str.contains("gondii", case=False, na=False)]
    if resolve is not None:
        d["gene_id"] = d.parasite_protein.astype(str).map(
            lambda s: resolve(re.sub(r"[-\s]", "", s.strip()).upper()))
    else:
        d["gene_id"] = None
    ok = d.gene_id.notna()
    log(f"host interactions: {len(d)} curated pairs, {int(ok.sum())} resolved to a gene "
        f"({d.loc[ok, 'gene_id'].nunique()} distinct parasite proteins, "
        f"{d.loc[ok, 'host_target'].nunique()} host targets)")
    if (~ok).any():
        log(f"  unresolved parasite symbols: {', '.join(sorted(set(d.loc[~ok, 'parasite_protein'].astype(str)))[:8])}")
    return d[ok].reset_index(drop=True)
