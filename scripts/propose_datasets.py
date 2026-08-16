#!/usr/bin/env python3
"""Propose accession-bearing papers for empty biological slots via NCBI E-utilities.

Queries are derived from the slot axis, context and organism; they are stored with every proposal so
the result is auditable. This writes candidates, never measurements, and the slot generator merges
the JSON into its rendered tables on the next run.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import sys
import time
import urllib.error
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


#: Organism names as GEO indexes them, for the `gds` search below.
GEO_ORGANISM = {"Tg": "Toxoplasma gondii", "Pf": "Plasmodium"}


def _series(slot: dict, retmax: int = 8) -> list:
    """Candidates from GEO itself rather than from papers about it.

    Worth having in addition to the PubMed pass, and on this data worth more. An accession reaches
    a PubMed abstract only when the authors put it there; most put it in data availability, which
    the abstract does not include -- so an abstract search finds the datasets that happened to be
    advertised rather than the datasets that exist. GEO's own index is keyed BY accession, so every
    hit is by construction something that can be fetched.
    """
    organism = GEO_ORGANISM.get(slot["organism"], "")
    assay = slots.ASSAY_TERMS.get(slot["axis"], "")
    context = re.sub(r"[^A-Za-z0-9 -]", " ", slot["context"]).strip()
    words = [w for w in context.split() if len(w) >= 4 and w.lower() != "gene"][:2]
    terms = [f'"{organism}"[Organism]', assay]
    if words:
        terms.append("(" + " OR ".join(f"{w}[All Fields]" for w in words) + ")")
    query = " AND ".join(term for term in terms if term)
    try:
        found = json.loads(_get("esearch.fcgi", {"db": "gds", "term": query, "retmode": "json",
                                                 "retmax": retmax}))
        uids = found.get("esearchresult", {}).get("idlist", [])
        if not uids:
            return []
        summary = json.loads(_get("esummary.fcgi", {"db": "gds", "id": ",".join(uids),
                                                    "retmode": "json"}))
    except (ValueError, urllib.error.URLError, urllib.error.HTTPError):
        return []
    out = []
    for uid in uids:
        record = summary.get("result", {}).get(uid) or {}
        accession = str(record.get("accession", ""))
        if not accession.startswith("GSE"):
            continue
        out.append({"pmid": str((record.get("pubmedids") or [""])[0]),
                    "title": str(record.get("title", "")),
                    "year": str(record.get("pdat", ""))[:4],
                    "journal": "GEO", "accession": accession,
                    "samples": int(record.get("n_samples", 0) or 0),
                    "query": query,
                    "note": "from the GEO index; verify assay and parasite-gene shape"})
    return out


#: Which repository actually holds each kind of measurement. GEO indexes sequencing; it does not
#: index mass spectrometry or metabolomics, and asking it for a palmitoylome returns the nearest
#: sequencing study instead of nothing. That is not a hypothetical: six different Plasmodium
#: modification slots and seven Toxoplasma ones all came back pointing at one lactylation paper and
#: one chromatin paper respectively, and 90 of 102 downloads failed verification outright.
#:
#: The keyword comes from the slot's NAME rather than its axis, because the name is where the
#: measurement is written -- "acetylation" and "palmitoylation" share an axis and are not the same
#: experiment.
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects"
PRIDE_SPECIES = {"Tg": "Toxoplasma gondii", "Pf": "Plasmodium falciparum"}

#: Slot-name fragments that mean "this is mass spectrometry", and the PRIDE keyword for each.
MASS_SPEC = {
    "phosphoryl": "phosphoproteome", "acetyl": "acetylation", "lactyl": "lactylation",
    "nitrosyl": "nitrosylation", "ubiquitin": "ubiquitination", "sumoyl": "SUMOylation",
    "glycosyl": "glycosylation", "palmitoyl": "palmitoylation",
    "protein abundance": "proteome", "turnover": "protein turnover",
    "crosslink ms": "crosslinking mass spectrometry", "ip-ms": "interactome",
    "proximity labelling": "proximity labeling", "secretome": "secretome",
    "thermal shift": "thermal proteome profiling", "with host proteins": "host interactome",
}


def _pride(slot: dict, retmax: int = 8) -> list:
    """Candidates from PRIDE, for the axes GEO cannot serve.

    Returned only for slots that are ABOUT mass spectrometry. A transcription slot has no business
    here, and offering it a proteomics dataset would be the same mistake in the other direction.
    """
    name = str(slot["name"]).lower()
    keyword = next((kw for fragment, kw in MASS_SPEC.items() if fragment in name), "")
    species = PRIDE_SPECIES.get(slot["organism"], "")
    if not keyword or not species:
        return []
    url = (f"{PRIDE}?keyword={urllib.parse.quote(keyword)}"
           f"&filter=organisms=={urllib.parse.quote(species)}&pageSize={retmax}")
    try:
        with urllib.request.urlopen(url, timeout=60) as fh:
            found = json.loads(fh.read().decode("utf8", "replace"))
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError):
        return []
    items = found if isinstance(found, list) else found.get("_embedded", {}).get("compactprojects", [])
    # PRIDE's keyword search is full-text and loose: `acetylation` against a Toxoplasma filter
    # returns the same top projects as a bare organism search. So the modification has to appear in
    # the RECORD, not merely in the query -- the identical mistake to the one that put six PTM slots
    # on a single lactylation paper, caught here instead of three steps downstream.
    stem = keyword.split()[0].rstrip("e")
    out = []
    for project in items:
        accession = str(project.get("accession", ""))
        if not accession.startswith("PXD"):
            continue
        blurb = f"{project.get('title', '')} {project.get('projectDescription', '')}".lower()
        if stem.lower() not in blurb:
            continue
        out.append({"pmid": "", "title": str(project.get("title", "")),
                    "year": str(project.get("publicationDate", ""))[:4], "journal": "PRIDE",
                    "accession": accession, "query": url,
                    "note": "from the PRIDE index; verify the modification and the organism"})
    return out


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
    # Accession-bearing only, not merely accession-first. A citation nobody can download is not a
    # filled slot, and offering one is how an empty slot comes to look filled -- which is the one
    # reading this table has to be trusted for. A slot with no downloadable candidate stays empty,
    # and empty is a true statement about what has been measured.
    return sorted((paper for paper in out if paper["accession"]), key=lambda paper: paper["pmid"])


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
            # Both sources: GEO's index first, because every hit there is fetchable by
            # construction, then the literature for anything GEO does not carry -- proteomics in
            # PRIDE, and the studies whose data went somewhere GEO does not index.
            found = (_pride(row, args.retmax) + _series(row, args.retmax)
                     + _papers(query, args.retmax))
        except (OSError, ValueError, ET.ParseError) as exc:
            print(f"{index}/{len(definitions)} {key}: {type(exc).__name__}")
            found = []
        seen, unique = set(), []
        for paper in found:
            if paper["accession"] in seen:
                continue
            seen.add(paper["accession"])
            unique.append(paper)
        found = unique
        if found:
            proposals[key] = [{**paper, "query": paper.get("query", query)} for paper in found]
        where = collections.Counter(paper.get("journal") or "literature" for paper in found)
        print(f"{index}/{len(definitions)} {key}: {len(found)} candidates "
              + ", ".join(f"{n} from {src}" for src, n in where.most_common()))
        time.sleep(max(args.delay, 0.0))
    with open(args.out, "w", encoding="utf8") as fh:
        json.dump(proposals, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(f"{sum(map(len, proposals.values()))} proposals -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
