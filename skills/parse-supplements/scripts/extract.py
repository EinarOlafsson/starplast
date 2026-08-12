#!/usr/bin/env python3
"""Extract text and gene identifiers from supplementary documents.

Publishers put hit lists in whatever they like. This reads PDF, DOCX, spreadsheets and plain text
uniformly, using only tools that are already on a normal Linux box -- `pdftotext`, `libreoffice`, and
Python's own zipfile for DOCX -- so nothing has to be pip-installed to read a paper's tables.

    python extract.py <dir-or-file> [--organism toxo] [--out found.tsv]

Two things it does that a naive text dump does not:

* **`pdftotext -layout`**, not the default flow. Without `-layout` a two-column table becomes interleaved
  prose and every row is corrupted; with it, columns stay in columns.
* **Symbols as well as accessions.** The documents that fail accession extraction usually fail because
  the paper lists proteins by name (GRA16, MIC2) rather than by ID. Matching symbols needs the same care
  as anywhere else: digit-free symbols must appear in upper case, because HOOK, CLAMP, CLIP, SPARK and
  REMIND are all real gene symbols and all real English words.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import zipfile

# Accession shapes, per organism. Toxoplasma is the default; the others are here because the same
# supplementary-parsing job recurs across the apicomplexan projects in this tree.
ACCESSIONS = {
    "toxo": re.compile(r"\bTG[A-Z0-9]{2,6}[_\-. ]?\d{5,6}\b", re.I),
    "plasmodium": re.compile(r"\bPF3D7[_\-. ]?\d{6,7}\b", re.I),
    "crypto": re.compile(r"\bcgd\d[_\-. ]?\d{3,5}\b", re.I),
}
# A token that could be a gene symbol. Unicode-aware on purpose: an ASCII class truncates accented
# words mid-token and invents matches out of them.
TOKEN = re.compile(r"[^\W\d_][\w-]{2,}", re.UNICODE)


def pdf_text(path: str, timeout=120) -> str:
    """PDF to text, preserving column layout.

    `-layout` matters more than it sounds: without it pdftotext reflows a table into prose and the
    columns interleave, so every extracted row is a mixture of two.
    """
    for args in (["pdftotext", "-layout", "-nopgbrk", path, "-"],
                 ["pdftotext", "-nopgbrk", path, "-"]):
        try:
            r = subprocess.run(args, capture_output=True, timeout=timeout)
            if r.returncode == 0 and r.stdout.strip():
                return r.stdout.decode("utf8", "replace")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return ""


def docx_text(path: str) -> str:
    """DOCX to text without python-docx: the format is a zip of XML.

    Tables are the point, so `</w:tc>` becomes a tab and `</w:tr>` a newline; otherwise every cell in a
    row runs together and a two-column table of gene/description is unreadable.
    """
    try:
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist()
                     if n.startswith("word/") and n.endswith(".xml")
                     and ("document" in n or "table" in n)]
            out = []
            for n in names or ["word/document.xml"]:
                try:
                    xml = z.read(n).decode("utf8", "replace")
                except KeyError:
                    continue
                xml = xml.replace("</w:tc>", "\t").replace("</w:tr>", "\n").replace("</w:p>", "\n")
                out.append(re.sub(r"<[^>]+>", "", xml))
            return "\n".join(out)
    except (zipfile.BadZipFile, OSError):
        return ""


def office_text(path: str, timeout=180) -> str:
    """Legacy .doc/.xls via libreoffice, which is slow but already installed."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        try:
            subprocess.run(["libreoffice", "--headless", "--convert-to", "txt:Text",
                            "--outdir", td, path], capture_output=True, timeout=timeout)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return ""
        for f in os.listdir(td):
            if f.endswith(".txt"):
                return open(os.path.join(td, f), encoding="utf8", errors="replace").read()
    return ""


def sheet_text(path: str) -> str:
    try:
        import pandas as pd
        xl = pd.ExcelFile(path)
        return "\n".join(xl.parse(sh, header=None, nrows=20000).astype(str).to_string()
                         for sh in xl.sheet_names[:10])
    except Exception:
        return office_text(path)


def text_of(path: str) -> str:
    """Whatever this file is, its text."""
    low = path.lower()
    if low.endswith(".pdf"):
        return pdf_text(path)
    if low.endswith(".docx"):
        return docx_text(path)
    if low.endswith((".xlsx", ".xls")):
        return sheet_text(path)
    if low.endswith((".txt", ".csv", ".tsv")):
        return open(path, encoding="utf8", errors="replace").read()
    if low.endswith(".doc"):
        return office_text(path)
    return ""


def find_accessions(text: str, organism="toxo") -> set:
    """Accessions, normalised to the underscore form publishers keep mangling."""
    rx = ACCESSIONS.get(organism, ACCESSIONS["toxo"])
    out = set()
    for m in rx.finditer(text or ""):
        s = re.sub(r"[\-. ]", "_", m.group(0).upper())
        s = re.sub(r"__+", "_", s)
        if "_" not in s:                       # TGME49208830 -> TGME49_208830
            s = re.sub(r"^([A-Z0-9]+?)(\d{5,6})$", r"\1_\2", s)
        out.add(s)
    return out


def find_symbols(text: str, symbols: set) -> set:
    """Gene symbols, with the guard that stops English words becoming genes.

    A symbol containing no digit must appear in upper case. HOOK, CLAMP, CLIP, SPARK and REMIND are all
    valid Toxoplasma symbols and all ordinary words; without this the extraction fills up with prose.
    """
    if not symbols:
        return set()
    upper_only = {s for s in symbols if not any(c.isdigit() for c in s)}
    out = set()
    for m in TOKEN.finditer(text or ""):
        tok = m.group(0)
        key = re.sub(r"[-\s]", "", tok).upper()
        if key in symbols and (key not in upper_only or tok.isupper()):
            out.add(key)
    return out


def walk(target: str):
    if os.path.isfile(target):
        yield target
    else:
        for root, _, files in os.walk(target):
            for f in sorted(files):
                if f != "META.json":
                    yield os.path.join(root, f)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="a file, or a directory to walk")
    ap.add_argument("--organism", default="toxo", choices=sorted(ACCESSIONS))
    ap.add_argument("--symbols", help="optional TSV of gene_id<TAB>symbol to also match by name")
    ap.add_argument("--out", help="write a TSV of (file, kind, identifier)")
    ap.add_argument("--text-dir", help="also dump the extracted text per file, for inspection")
    a = ap.parse_args(argv)

    symbols = {}
    if a.symbols and os.path.exists(a.symbols):
        for line in open(a.symbols, encoding="utf8", errors="replace"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].strip() and parts[1].strip().upper() != "N/A":
                symbols[re.sub(r"[-\s]", "", parts[1].strip()).upper()] = parts[0].strip()

    rows, n_files, n_empty = [], 0, 0
    for path in walk(a.target):
        txt = text_of(path)
        n_files += 1
        if not txt.strip():
            n_empty += 1
            continue
        if a.text_dir:
            os.makedirs(a.text_dir, exist_ok=True)
            flat = path.replace(os.sep, "__").lstrip("_.")
            open(os.path.join(a.text_dir, flat + ".txt"), "w").write(txt)
        for acc in sorted(find_accessions(txt, a.organism)):
            rows.append((path, "accession", acc))
        for sym in sorted(find_symbols(txt, set(symbols))):
            rows.append((path, "symbol", f"{sym}\t{symbols[sym]}"))

    print(f"{n_files} files read, {n_empty} yielded no text", file=sys.stderr)
    accs = {r[2] for r in rows if r[1] == "accession"}
    syms = {r[2] for r in rows if r[1] == "symbol"}
    print(f"{len(accs)} distinct accessions, {len(syms)} distinct symbols", file=sys.stderr)

    out = open(a.out, "w") if a.out else sys.stdout
    out.write("file\tkind\tidentifier\n")
    for p, k, v in rows:
        out.write(f"{p}\t{k}\t{v}\n")
    if a.out:
        out.close()
        print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
