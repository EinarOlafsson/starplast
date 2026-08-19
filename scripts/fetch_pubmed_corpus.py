#!/usr/bin/env python3
"""Fetch a PubMed abstract corpus for one organism, in the shape the literature layer reads.

The Toxoplasma arm's attention layer -- who has been studied, how deeply, and which genes are named
together -- is built from 33,924 abstracts sitting outside the repository. The Plasmodium arm has no
such corpus, which is why five of its slots are empty for a reason no acquisition can fix: there is
nothing to count. This fetches the same thing for any organism, so the two arms are built from one
construction rather than two.

    python scripts/fetch_pubmed_corpus.py --query '"Plasmodium falciparum"[Title/Abstract]' \
        --out .claude/skills/plasmodium-scientist/corpus/pubmed_plasmodium.jsonl

Written one record per line, `{pmid, year, journal, title, abstract, mesh, pubtypes}` -- the format
`corpus.py` already parses, chosen by reading a record of the existing file rather than by inventing
one. Resumable: an existing output is read first and its PMIDs are skipped, because 43,000 abstracts
is twenty minutes and a network hiccup should not cost all of it.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: NCBI asks for a tool and an address on every call, and gives 3 requests a second without a key.
CONTACT = {"tool": "starplast", "email": "einar.olafsson@gmail.com"}
DELAY = 0.34

#: How many PMIDs to ask for in one efetch. 200 is what NCBI's own documentation suggests for
#: efetch by POST; larger batches time out rather than fail cleanly, which is the worse failure.
BATCH = 200


def _get(endpoint: str, params: dict, tries: int = 3) -> bytes:
    """One E-utilities call, retried, because a 500 in the middle of 200 batches is normal."""
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode({**params, **CONTACT})}"
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except (urllib.error.URLError, OSError, TimeoutError):
            if attempt == tries - 1:
                raise
            time.sleep(2 ** attempt)
    return b""


def search(query: str, log=print) -> list:
    """Every PMID matching the query, gathered one YEAR at a time.

    NCBI stops at ten thousand, twice over: esearch will not page past it, and efetch answers
    `400 Bad Request` for a `retstart` beyond it. Both limits are silent in their own way -- the
    first returns an empty page, which reads as the end of the list, and the first version of this
    script fetched 9,989 abstracts of 43,482 and reported success.

    A year is the natural slice: no year of this literature comes near the cap, and the partition is
    reproducible, so the same query re-run tomorrow asks the same questions. A year that DOES exceed
    the cap is reported rather than truncated silently.
    """
    ids = []
    for year in range(1960, 2031):
        found = json.loads(_get("esearch.fcgi", {
            "db": "pubmed", "term": f'({query}) AND ("{year}"[PDAT] : "{year}"[PDAT])',
            "retmode": "json", "retmax": 9999}), strict=False)["esearchresult"]
        chunk = found.get("idlist", [])
        if int(found.get("count", 0)) > 9999:
            log(f"  {year}: {found['count']} papers, more than one request returns -- TRUNCATED")
        if chunk:
            ids.extend(chunk)
            log(f"  {year}: {len(chunk):,} ({len(ids):,} so far)")
        time.sleep(DELAY)
    unique = list(dict.fromkeys(ids))
    log(f"{len(unique):,} distinct papers")
    return unique


def _text(node) -> str:
    """All text under a node, including the tails of any markup inside it."""
    return "".join(node.itertext()).strip() if node is not None else ""


def records(pmids: list, log=print):
    """Yield one dict per paper, in the shape the corpus file already uses."""
    for start in range(0, len(pmids), BATCH):
        batch = pmids[start:start + BATCH]
        try:
            root = ET.fromstring(_get("efetch.fcgi", {"db": "pubmed", "retmode": "xml",
                                                      "id": ",".join(batch)}))
        except (ET.ParseError, urllib.error.URLError, OSError) as exc:
            log(f"  batch at {start}: {type(exc).__name__}, skipped")
            continue
        for article in root.findall(".//PubmedArticle"):
            pmid = _text(article.find(".//PMID"))
            if not pmid:
                continue
            abstract = " ".join(_text(part) for part in article.findall(".//AbstractText"))
            yield {"pmid": pmid,
                   "year": _text(article.find(".//PubDate/Year")) or _text(
                       article.find(".//PubDate/MedlineDate"))[:4],
                   "journal": _text(article.find(".//Journal/ISOAbbreviation")),
                   "title": _text(article.find(".//ArticleTitle")),
                   "abstract": abstract,
                   "mesh": [_text(m.find("DescriptorName"))
                            for m in article.findall(".//MeshHeading")],
                   "pubtypes": [_text(p) for p in article.findall(".//PublicationType")]}
        log(f"  {min(start + BATCH, len(pmids)):,}/{len(pmids):,} fetched")
        time.sleep(DELAY)


def main(argv=None, log=print) -> int:
    """Fetch, skipping whatever the output file already holds."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="a PubMed query")
    parser.add_argument("--out", required=True, help="JSONL to write, appended if it exists")
    parser.add_argument("--limit", type=int, default=0, help="stop after this many, for a smoke test")
    args = parser.parse_args(argv)

    have = set()
    if os.path.exists(args.out):
        with open(args.out, encoding="utf8") as fh:
            for line in fh:
                try:
                    have.add(json.loads(line)["pmid"])
                except (ValueError, KeyError):
                    continue
        log(f"{len(have):,} already in {args.out}")
    pmids = [p for p in search(args.query, log=log) if p not in have]
    if args.limit:
        pmids = pmids[:args.limit]
    log(f"{len(pmids):,} to fetch")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    written = skipped = 0
    with open(args.out, "a", encoding="utf8") as fh:
        for record in records(pmids, log=log):
            if record["pmid"] in have:
                skipped += 1
                continue
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")
            written += 1
    if skipped:
        log(f"{skipped:,} already had")
    log(f"wrote {written:,} records to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
