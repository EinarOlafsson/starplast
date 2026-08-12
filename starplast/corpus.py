#!/usr/bin/env python3
"""Documents: one uniform stream over the abstract corpus and the open-access full texts.

Two very different sources have to be read the same way but must never be silently merged, because they
support different claims:

* **abstracts** -- 33,924 PubMed records. Complete coverage of the field, but a gene named only in a
  paper's body is invisible. This is the source behind the "at most 7.4% of the proteome is named"
  figure, and it is a lower bound on attention.
* **fulltext** -- the open-access XML on local disk. Far deeper per paper, but a biased subset: only
  papers whose publisher deposited them. Coverage computed over full texts answers a different question
  and cannot be quoted as if it covered the field.

Every Document therefore carries its `source`, and every Section carries its `kind`, so downstream code
chooses granularity explicitly instead of inheriting whatever the parser happened to concatenate.

Section kinds: ``title``, ``abstract``, ``body``, ``caption``.

Reference lists are excluded. A reference list reproduces the titles of cited papers, so including it
would credit a paper with mentioning every gene named in its bibliography -- inflating both coverage and
co-mention with the citation graph rather than the paper's own content.
"""
from __future__ import annotations

import glob
import json
import os
import re
from dataclasses import dataclass
from xml.etree import ElementTree as ET


@dataclass(frozen=True)
class Section:
    """One titled section of a full text, with its body."""
    kind: str
    text: str


@dataclass(frozen=True)
class Document:
    """One article: its identifiers, its sections, and the text they contain."""
    doc_id: str            # "pmid:12345678" or "pmc:PMC10000077"
    source: str            # "abstract" | "fulltext"
    pmid: str | None       # lets the two sources be linked and de-duplicated
    year: str | None
    sections: tuple


# --------------------------------------------------------------------------- abstracts
def iter_abstracts(path: str):
    """Yield one Document per PubMed record in the JSONL corpus."""
    if not os.path.exists(path):
        return
    with open(path, encoding="utf8", errors="replace") as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except Exception:
                continue
            pmid = str(r.get("pmid") or "").strip() or None
            secs = []
            if r.get("title"):
                secs.append(Section("title", r["title"]))
            if r.get("abstract"):
                secs.append(Section("abstract", r["abstract"]))
            if not secs:
                continue
            yield Document(f"pmid:{pmid}", "abstract", pmid, str(r.get("year") or "") or None,
                           tuple(secs))


# --------------------------------------------------------------------------- full texts
def _text(el) -> str:
    """Flatten an element's descendant text, which JATS scatters across <italic>, <sup> and friends."""
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def _strip(root, tags=("ref-list",)) -> None:
    """Drop reference lists (see module docstring) in place, anywhere in the tree."""
    # Collect first, then remove: mutating while ElementTree.iter() walks the tree skips siblings.
    doomed = [(parent, child) for parent in root.iter()
              for child in list(parent) if child.tag in tags]
    for parent, child in doomed:
        parent.remove(child)


def parse_jats(path: str) -> Document | None:
    """Parse one PMC JATS file into a sectioned Document, or None if it is unreadable."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    # Strip references from the whole tree before reading anything: <ref-list> lives in <back>, and the
    # title/abstract scans below walk the full tree.
    _strip(root)

    pmid = pmcid = year = None
    for aid in root.iter("article-id"):
        t = aid.get("pub-id-type")
        if t == "pmid":
            pmid = (aid.text or "").strip() or None
        elif t == "pmcid":
            pmcid = (aid.text or "").strip() or None
    for d in root.iter("pub-date"):
        y = d.find("year")
        if y is not None and (y.text or "").strip().isdigit():
            year = y.text.strip()
            break

    secs = []
    for t in root.iter("article-title"):
        if _text(t):
            secs.append(Section("title", _text(t)))
        break
    for ab in root.iter("abstract"):
        if _text(ab):
            secs.append(Section("abstract", _text(ab)))

    body = root.find("body")
    if body is not None:
        # Figure and table captions wrap their text in <p>, so an unguarded body scan emits caption text
        # twice -- once as body, once as caption -- double-counting every gene named in a caption and the
        # co-mention unit it sits in.
        in_caption = {id(p) for cap in body.iter("caption") for p in cap.iter("p")}
        # Paragraph-level sections: co-mention inside one paragraph is a claim about a relation, whereas
        # co-occurrence anywhere in a 10,000-word paper mostly is not.
        for p in body.iter("p"):
            s = _text(p)
            if id(p) not in in_caption and len(s) > 40:
                secs.append(Section("body", s))
        for cap in body.iter("caption"):
            s = _text(cap)
            if len(s) > 20:
                secs.append(Section("caption", s))

    if not secs:
        return None
    doc_id = f"pmc:{pmcid}" if pmcid else f"file:{os.path.basename(path)}"
    return Document(doc_id, "fulltext", pmid, year, tuple(secs))


def iter_fulltexts(directory: str, limit: int | None = None):
    """Yield one Document per PMC XML file in `directory`."""
    files = sorted(glob.glob(os.path.join(directory, "*.xml")))
    if limit:
        files = files[:limit]
    for f in files:
        doc = parse_jats(f)
        if doc is not None:
            yield doc


def count_fulltexts(directory: str) -> int:
    """How many open-access full texts are present locally. A biased subset, by construction."""
    return len(glob.glob(os.path.join(directory, "*.xml")))
