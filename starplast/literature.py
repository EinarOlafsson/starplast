#!/usr/bin/env python3
"""Mentions: the one auditable intermediate between documents and every literature-derived number.

Before this module, per-gene publication counts and co-mention pairs were accumulated in a single pass
inside build_graph, with no intermediate. Nothing could be checked, and no figure could be recomputed at a
different confidence threshold without re-scanning the whole corpus. Everything downstream -- coverage,
attention correction, co-mention edges -- now derives from one tidy table:

    gene_id | doc_id | source | match_kind | section | n_hits

`source` keeps abstracts and full texts separate (they answer different questions -- see corpus.py) and
`match_kind` keeps confidence tiers separate (see identity.py), so a coverage claim can always name the
evidence it rests on.

Co-mention is counted per *unit*, not per document. For an abstract the unit is the record; for a full
text it is the paragraph, because two genes named in the same paragraph plausibly stand in a relation
while two genes named anywhere in a 10,000-word paper mostly do not.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

# Confidence tiers over identity.KINDS. Accession matches are near-certain; a symbol matched in a
# full-text body is the form that can collide with a host-gene name.
TIER = {"accession": "accession", "accession_prev": "accession", "accession_strain": "accession",
        "symbol": "symbol", "alias": "symbol"}

# A unit naming this many genes is a list, a table dump or a screen summary. It evidences membership of a
# set, not a pairwise relation, so it is excluded from co-mention (but still counts for coverage).
MAX_GENES_PER_UNIT = 12

# Section kinds that count as the paper's own claims.
SECTION_ORDER = ("title", "abstract", "body", "caption")


def scan(docs, index, max_genes_per_unit: int = MAX_GENES_PER_UNIT, log=print):
    """Scan Documents against a GeneIndex.

    Returns (mentions, co, meta) where `mentions` is the tidy DataFrame described above, `co` maps
    source -> Counter of (gene_a, gene_b) -> number of units naming both, and `meta` carries the unit
    counts that the attention correction needs as its denominator.
    """
    rows = defaultdict(int)                      # (gene, doc, source, kind, section) -> hits
    co = defaultdict(Counter)                    # source -> Counter[(g1,g2)]
    unit_hits = defaultdict(Counter)             # source -> Counter[gene] -> units naming it
    n_docs = Counter()
    n_units = Counter()
    n_units_skipped = Counter()

    for doc in docs:
        n_docs[doc.source] += 1
        # An abstract record is one unit; a full text is one unit per section, so co-mention keeps its
        # meaning as the corpus gets deeper.
        if doc.source == "abstract":
            units = [doc.sections]
        else:
            units = [(s,) for s in doc.sections]

        for unit in units:
            in_unit = set()
            for sec in unit:
                for gene, kind in index.find(sec.text):
                    rows[(gene, doc.doc_id, doc.pmid or "", doc.source, kind, sec.kind)] += 1
                    in_unit.add(gene)
            if len(in_unit) > max_genes_per_unit:
                n_units_skipped[doc.source] += 1
                continue
            # Counted after the list exclusion, so the attention correction's observed and expected
            # counts come from one population. Counting n_units over every unit while counting
            # co-mention over only the eligible ones inflates every expectation for genes that appear
            # mostly in screen tables, and quietly drags their residuals below zero.
            n_units[doc.source] += 1
            for gene in in_unit:
                unit_hits[doc.source][gene] += 1
            if len(in_unit) < 2:
                continue
            g = sorted(in_unit)
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    co[doc.source][(g[i], g[j])] += 1

    m = pd.DataFrame(
        [(g, did, pmid, src, kind, sec, n)
         for (g, did, pmid, src, kind, sec), n in rows.items()],
        columns=["gene_id", "doc_id", "pmid", "source", "match_kind", "section", "n_hits"])
    m["tier"] = m.match_kind.map(TIER)

    for src, n in sorted(n_docs.items()):
        got = m[m.source == src]
        log(f"{src}: {n:,} documents scanned ({n_units[src]:,} co-mention units); "
            f"{got.gene_id.nunique():,} genes mentioned; "
            f"{got.doc_id.nunique():,} documents with >=1 gene; "
            f"{len(co[src]):,} co-mention pairs ({n_units_skipped[src]:,} units skipped as lists)")
    meta = {"n_docs": dict(n_docs), "n_units": dict(n_units),
            "n_units_skipped": dict(n_units_skipped),
            "unit_hits": {k: dict(v) for k, v in unit_hits.items()}}
    return m, co, meta


def comention_edges(co: Counter, meta: dict, source: str, min_count: int = 2):
    """Raw co-mention counts and their attention-corrected residuals for one source.

    Attention correction (design decision 3, and a correctness feature rather than decoration): two genes
    named in n1 and n2 of N units are expected under independence to share n1*n2/N units. What the app
    shows by default is the residual log2((observed + 0.5) / (expected + 0.5)), because raw co-mention
    reproduces the field's popularity contest -- GRA16 and ROP18 become bright hubs on attention alone.

    The denominator is the number of *units*, which is what independence is defined over. An earlier
    version divided by the sum of per-gene counts, a larger and differently-scaled quantity, which
    inflated every expectation and shrank every residual toward zero.
    """
    hits = meta["unit_hits"].get(source, {})
    n = meta["n_units"].get(source, 0) or 1
    out = []
    for (g1, g2), c in co.items():
        if c < min_count:                       # a single shared unit is not a relation
            continue
        exp = hits.get(g1, 0) * hits.get(g2, 0) / n
        out.append((g1, g2, float(c), float(math.log2((c + 0.5) / (exp + 0.5)))))
    return out


# --------------------------------------------------------------------------- derived figures
def coverage(mentions: pd.DataFrame) -> pd.DataFrame:
    """Genes named, broken out by source and confidence tier -- never as one undifferentiated number."""
    out = []
    for src in sorted(mentions.source.unique()):
        s = mentions[mentions.source == src]
        for tier in ("accession", "symbol"):
            t = s[s.tier == tier]
            out.append({"source": src, "tier": tier, "genes": t.gene_id.nunique(),
                        "documents": t.doc_id.nunique()})
        out.append({"source": src, "tier": "any", "genes": s.gene_id.nunique(),
                    "documents": s.doc_id.nunique()})
    # Documents are de-duplicated by PMID for the union row: an open-access paper is usually also a
    # PubMed record, so summing the two document counts would count it twice.
    key = mentions.pmid.where(mentions.pmid.astype(str) != "", mentions.doc_id)
    out.append({"source": "union", "tier": "any", "genes": mentions.gene_id.nunique(),
                "documents": key.nunique()})
    return pd.DataFrame(out)


def publication_counts(mentions: pd.DataFrame, source: str | None = None) -> pd.Series:
    """Distinct documents naming each gene. This is the denominator of the attention correction."""
    m = mentions if source is None else mentions[mentions.source == source]
    return m.groupby("gene_id").doc_id.nunique()


# Where in a paper a gene is named says what the paper did with it. These tiers are read off document
# structure rather than assigned as weights, because an invented weighting is exactly the kind of made-up
# tier this project refuses elsewhere.
#
#   focal        named in a title            -> the paper is about this gene
#   substantive  named in an abstract        -> it is a stated part of the paper's claims
#   incidental   named only in body/caption  -> mentioned in passing, or listed in a screen's hit table
#
# The distinction matters because coverage counts every tier alike. A gene reached only through a
# supplementary table of 2,000 hits is *named*, not *studied*, and reporting it as attention would repeat
# the error the attention correction exists to prevent.
DEPTH_OF = {"title": "focal", "abstract": "substantive", "body": "incidental",
            "caption": "incidental"}
DEPTH_ORDER = ("focal", "substantive", "incidental")


def attention_depth(mentions: pd.DataFrame) -> pd.DataFrame:
    """Per gene: how many papers name it where, and the strongest tier any paper affords it.

    Papers are de-duplicated by PMID across the two sources. Note that the abstract corpus carries only
    title and abstract text, so it can never contribute an `incidental` paper -- incidental evidence is
    visible only where a full text is available. That asymmetry limits the counts, not the tiering: a
    title mention means the same thing whichever source it came from.
    """
    if mentions.empty:
        return pd.DataFrame(columns=["n_papers_focal", "n_papers_substantive",
                                     "n_papers_incidental", "attention_depth"])
    m = mentions.copy()
    m["paper"] = m.pmid.where(m.pmid.astype(str) != "", m.doc_id)
    m["depth"] = m.section.map(DEPTH_OF).fillna("incidental")

    # One tier per (gene, paper): the strongest place that paper names the gene.
    rank = {d: i for i, d in enumerate(DEPTH_ORDER)}
    best = (m.assign(r=m.depth.map(rank)).sort_values("r")
            .drop_duplicates(["gene_id", "paper"]))
    wide = (best.pivot_table(index="gene_id", columns="depth", values="paper",
                             aggfunc="nunique", fill_value=0))
    for d in DEPTH_ORDER:
        if d not in wide.columns:
            wide[d] = 0
    out = pd.DataFrame({f"n_papers_{d}": wide[d] for d in DEPTH_ORDER})
    out["attention_depth"] = np.select(
        [out.n_papers_focal > 0, out.n_papers_substantive > 0, out.n_papers_incidental > 0],
        list(DEPTH_ORDER), default="")
    return out
