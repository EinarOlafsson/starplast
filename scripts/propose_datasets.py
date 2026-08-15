#!/usr/bin/env python3
"""Propose accession-bearing papers for empty biological slots via NCBI E-utilities.

Queries are derived from the slot axis, context and organism; they are stored with every proposal so
the result is auditable. This writes candidates, never measurements, and the slot generator merges
the JSON into its rendered tables on the next run.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_slot_table as slots  # noqa: E402

ACCESSION = re.compile(r"\b(?:GSE\d+|PXD\d+|PRJ(?:NA|EB)\d+)\b", re.I)
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def query_for(slot: dict) -> str:
    """The reproducible PubMed query for one slot."""
    organism = slots.ORGANISMS[slot["organism"]]["query"]
    assay = slots.ASSAY_TERMS.get(slot["axis"], "")
    context = re.sub(r"[^A-Za-z0-9 -]", " ", slot["context"]).strip()
    words = [word for word in context.split() if len(word) >= 4 and word.lower() not in {"gene"}]
    context_term = ("(" + " OR ".join(f"{word}[Title/Abstract]" for word in words) + ")"
                    if words else "")
    terms = [organism, assay, context_term]
    return " AND ".join(term for term in terms if term)


def _get(endpoint: str, params: dict) -> bytes:
    params = {**params, "tool": "starplast", "email": "einar.olafsson@gmail.com"}
    url = f"{BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as response:
        return response.read()


def _papers(query: str, retmax: int = 8) -> list:
    search = json.loads(_get("esearch.fcgi", {"db": "pubmed", "term": query,
                                               "retmode": "json", "retmax": retmax}))
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    root = ET.fromstring(_get("efetch.fcgi", {"db": "pubmed", "id": ",".join(ids),
                                               "retmode": "xml"}))
    out = []
    for article in root.findall(".//PubmedArticle"):
        pmid = "".join(article.findtext(".//PMID", default=""))
        title = "".join(article.find(".//ArticleTitle").itertext()) \
            if article.find(".//ArticleTitle") is not None else ""
        abstract = " ".join("".join(node.itertext()) for node in article.findall(".//AbstractText"))
        accessions = sorted(set(ACCESSION.findall(f"{title} {abstract}")), key=str.upper)
        year = (article.findtext(".//PubDate/Year") or article.findtext(".//ArticleDate/Year") or "")
        journal = article.findtext(".//Journal/Title", default="")
        out.append({"pmid": pmid, "title": title, "year": year, "journal": journal,
                    "accession": ", ".join(a.upper() for a in accessions),
                    "note": "automatically proposed; verify assay and parasite-gene shape"})
    return sorted(out, key=lambda paper: (not bool(paper["accession"]), paper["pmid"]))


def main(argv=None) -> int:
    """Search empty slots and write the mergeable candidate catalog."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=slots.CANDIDATES_JSON)
    parser.add_argument("--retmax", type=int, default=8)
    parser.add_argument("--max-slots", type=int, default=0,
                        help="limit work for a smoke test; zero searches every empty slot")
    parser.add_argument("--delay", type=float, default=0.34,
                        help="seconds between NCBI calls (3 requests/s without an API key)")
    args = parser.parse_args(argv)
    definitions = [row for row in slots.all_slots() if not row["patterns"]
                   and row["axis"] != "NEVER a feature"]
    if args.max_slots:
        definitions = definitions[:args.max_slots]
    proposals = {}
    for index, row in enumerate(definitions, 1):
        query = query_for(row)
        key = f"{row['organism']}::{row['name']}"
        try:
            found = _papers(query, args.retmax)
        except (OSError, ValueError, ET.ParseError) as exc:
            print(f"{index}/{len(definitions)} {key}: {type(exc).__name__}")
            found = []
        if found:
            proposals[key] = [{**paper, "query": query} for paper in found]
        with_accession = sum(bool(paper["accession"]) for paper in found)
        print(f"{index}/{len(definitions)} {key}: {len(found)} papers, "
              f"{with_accession} with an abstract accession")
        time.sleep(max(args.delay, 0.0))
    with open(args.out, "w", encoding="utf8") as fh:
        json.dump(proposals, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(f"{sum(map(len, proposals.values()))} proposals -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
