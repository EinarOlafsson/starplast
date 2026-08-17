#!/usr/bin/env python3
"""Gene identity: resolve every string the literature uses for a gene to one canonical gene_id.

This layer exists because the literature does not use one identifier. The same gene appears as
``TGME49_208830``, as the pre-2012 accession ``TGME49_008830``, as the GT1-strain accession
``TGGT1_208830``, as the symbol ``GRA16``, and as ``TgGRA16``. Before this module, build_graph matched
only current ME49 accessions and ToxoDB symbols, so every other form was silently discarded -- it did not
even fail loudly, because unmatched accessions were dropped by an intersection with the node table.

Design decisions:

* **Ambiguity is recorded, never guessed.** A string resolving to more than one gene is excluded from
  lookup AND kept in ``ambiguous`` with its collision set, so the cost of the exclusion is reportable
  instead of invisible. This follows the project's existing rule for ambiguous symbols.
* **Every match carries a kind.** Accession matches are near-certain; symbol matches in a full-text body
  are the ones that can collide with host-gene names. Downstream code can therefore report coverage by
  confidence tier rather than as one undifferentiated number.
* **Cross-strain accessions map by numeric suffix.** ``TGGT1_208830`` -> ``TGME49_208830``. This is the
  VEuPathDB naming convention for syntenic orthologs, and it is verified, not assumed: of the ME49/GT1
  pairs sharing a suffix and an OrthoMCL release, 99.53% share an orthogroup (VEG: 99.46%). The mapping is
  accepted only when the suffix exists in the node table.
* **Symbols need a shape that survives free text.** A symbol is matchable only if it contains a digit or
  is four or more capitals -- "AAP" would otherwise match prose, while "GRA16" and "MORC" would not.

Build the index from ``data/toxodb_identity.tsv`` (see fetch_names.py) plus the node table.
"""
from __future__ import annotations

import os
import re
from collections import defaultdict
from dataclasses import dataclass, field

import pandas as pd

# Match kinds, most trustworthy first. Downstream tiering relies on this order.
KINDS = ("accession", "accession_prev", "accession_strain", "symbol", "alias")

# Any Toxoplasma-style accession: TGME49_208830, TGGT1_208830, TGVEG_208830, old TGME49_008830.
ACC_RX = re.compile(r"\bTG[A-Z0-9]{2,6}_(\d{5,6})[A-Za-z]?\b", re.I)

# Word-ish tokens for symbol lookup. Hyphens are kept so "GRA-16" survives tokenisation and is then
# normalized to "GRA16"; a giant regex alternation over ~3,000 symbols is far slower than set lookup.
#
# The character class is Unicode-aware on purpose. An ASCII-only class ([A-Za-z]) truncates accented words
# mid-token, and this corpus carries French, Portuguese and Spanish abstracts: "Santé" tokenised as "Sant"
# and handed 23 abstracts to TGME49_321440, whose symbol is SANT.
TOKEN_RX = re.compile(r"[^\W\d_][\w-]{2,}", re.UNICODE)


def norm(s: str) -> str:
    """Normalize a candidate string to a lookup key: upper case, no hyphens/spaces/underscores."""
    return re.sub(r"[-\s_.]", "", (s or "").strip().upper())


# Strings that are ToxoDB symbols but read as something else in this literature. Toxoplasma strain and
# lineage designations are the whole problem: "the GT1 strain" appears in thousands of papers, and
# TGME49_214320 is symbolled GT1 (glucose transporter 1), so an unfiltered match hands that one gene a
# large and entirely fictional share of the field's attention. The list is deliberately restricted to
# closed, well-known naming systems -- strains, and the ME49/GT1/VEG reference isolates -- rather than
# growing into a general stopword list, which would quietly suppress real genes.
BLOCKED_SYMBOLS = frozenset(norm(s) for s in (
    "ME49", "GT1", "VEG", "RH", "PRU", "PLK", "CTG", "COUG", "VAND", "MAS", "ARI", "FOU", "GAB", "RUB",
    "CEP", "BOF", "P89", "B41", "B73", "WH3", "WH6", "QHO", "TgCkAu", "TgCatBr5", "TgToucan",
))


def matchable_symbol(sym: str) -> bool:
    """True if a symbol is specific enough to match in running prose without catching English."""
    sym = (sym or "").strip()
    if len(sym) < 3 or sym.upper() in ("N/A", "NA", "NONE"):
        return False
    if norm(sym) in BLOCKED_SYMBOLS:
        return False
    return any(c.isdigit() for c in sym) or (sym.isupper() and len(sym) >= 4)


@dataclass
class GeneIndex:
    """Canonical gene ids plus every string that resolves to exactly one of them."""

    canonical: set = field(default_factory=set)
    lookup: dict = field(default_factory=dict)          # key -> (gene_id, kind)
    ambiguous: dict = field(default_factory=dict)       # key -> {gene_id, ...}
    suffix: dict = field(default_factory=dict)          # accession numeric suffix -> gene_id
    strict: set = field(default_factory=set)            # keys that only match an upper-case surface form

    # ------------------------------------------------------------------ construction
    def _add(self, key: str, gene_id: str, kind: str, strict: bool = False) -> None:
        key = norm(key)
        if not key or gene_id not in self.canonical:
            return
        if strict:
            self.strict.add(key)
        if key in self.ambiguous:
            self.ambiguous[key].add(gene_id)
            return
        prev = self.lookup.get(key)
        if prev is None:
            self.lookup[key] = (gene_id, kind)
            return
        if prev[0] == gene_id:
            # same gene, keep the more trustworthy kind
            if KINDS.index(kind) < KINDS.index(prev[1]):
                self.lookup[key] = (gene_id, kind)
            return
        # Two genes claim one string. Ambiguity is declared only WITHIN a kind tier: a gene's own current
        # accession must not be withdrawn because some other gene once carried that string as a previous
        # id. Across tiers the more specific identifier wins; within a tier we refuse to guess.
        pk, nk = KINDS.index(prev[1]), KINDS.index(kind)
        if nk < pk:
            self.lookup[key] = (gene_id, kind)
            return
        if nk > pk:
            return
        self.ambiguous[key] = {prev[0], gene_id}
        del self.lookup[key]

    def stats(self) -> dict:
        """Counts of how many accessions resolved, and how many were ambiguous."""
        by_kind = defaultdict(int)
        for _, kind in self.lookup.values():
            by_kind[kind] += 1
        return {"genes": len(self.canonical), "keys": len(self.lookup),
                "ambiguous_keys": len(self.ambiguous), **{f"k_{k}": by_kind[k] for k in KINDS}}

    # ------------------------------------------------------------------ matching
    def find(self, text: str):
        """Yield (gene_id, kind) for every gene mentioned in `text`. Duplicates are yielded per hit."""
        if not text:
            return
        spans = []
        for m in ACC_RX.finditer(text):
            spans.append(m.span())
            # Registered accessions (current, previous, and strain forms) resolve through `lookup`, which
            # carries the right kind. The suffix map is only a fallback for a strain accession that the
            # strain tables did not cover; an unregistered ME49-shaped id resolves to nothing on purpose.
            hit = self.lookup.get(norm(m.group(0)))
            if hit is None:
                gid = self.suffix.get(m.group(1))
                if gid is None:
                    continue
                hit = (gid, "accession" if norm(m.group(0)).startswith("TGME49")
                       else "accession_strain")
            yield hit
        for m in TOKEN_RX.finditer(text):
            # An accession's own prefix must not be re-read as a symbol. Without this, every "TGGT1_209030"
            # in the corpus also scored a mention of TGME49_214320, whose symbol is "GT1".
            if any(s <= m.start() < e for s, e in spans):
                continue
            key = norm(m.group(0))
            hit = self.lookup.get(key)
            if hit is None:
                continue
            # Digit-free symbols are ordinary English words as often as they are genes -- HOOK, CLAMP,
            # CLIP, SPARK, REMIND are all real Toxoplasma symbols and all real words. A paper writing the
            # gene capitalises it, so those keys require an upper-case surface form; symbols carrying a
            # digit ("GRA16") cannot collide with prose and stay case-insensitive.
            if key in self.strict and not m.group(0).isupper():
                continue
            yield hit


def build_index(node_ids, identity_tsv: str, log=print) -> GeneIndex:
    """Build a GeneIndex over `node_ids` from the ToxoDB identity table."""
    ix = GeneIndex(canonical=set(node_ids))

    # Accession suffix map: current ME49 accessions define the suffix space that every other accession
    # form (old ids, GT1/VEG) resolves into.
    for gid in ix.canonical:
        part = gid.split("_")[-1]
        if part.isdigit():
            ix.suffix.setdefault(part, gid)
        ix._add(gid, gid, "accession")

    if not os.path.exists(identity_tsv):
        log(f"identity table absent ({identity_tsv}) -- run python -m starplast.fetch_names; "
            "only current accessions will resolve")
        return ix

    t = pd.read_csv(identity_tsv, sep="\t", dtype=str, keep_default_na=False)
    need = {"gene_id", "gene_name", "previous_ids"}
    if not need.issubset(t.columns):
        raise ValueError(f"{identity_tsv} missing columns: {need - set(t.columns)}")

    n_prev = n_sym = 0
    for row in t.itertuples(index=False):
        gid = row.gene_id
        if gid not in ix.canonical:
            continue
        for p in str(row.previous_ids).split(";"):
            p = p.replace("Previous IDs:", "").strip()
            # Current/previous accessions serve literature matching.  The old ``12.m00123``
            # ToxoGeneChip models are also registered for structured dataset joins (GPL7186); the
            # free-text matcher cannot mistake them for prose because its token regex does not span
            # the leading number and dot. Other historical gene-model systems remain excluded.
            m = re.fullmatch(r"TG[A-Z0-9]{2,6}_(\d{5,6})[A-Za-z]?", p, re.I)
            old_chip = re.fullmatch(r"\d+\.m\d+", p, re.I)
            if (m or old_chip) and norm(p) != norm(gid):
                ix._add(p, gid, "accession_prev")
                n_prev += 1
        sym = str(row.gene_name).strip()
        if matchable_symbol(sym):
            ix._add(sym, gid, "symbol", strict=not any(c.isdigit() for c in sym))
            n_sym += 1
            # The Tg- prefix convention: papers write TgGRA16 as often as GRA16. The prefixed form is
            # never an English word, so it does not need the upper-case restriction.
            if not norm(sym).startswith("TG"):
                ix._add("TG" + sym, gid, "alias")

    log(f"identity: {len(ix.canonical):,} genes; {len(ix.lookup):,} resolvable strings "
        f"({n_sym:,} symbols, {n_prev:,} previous accessions offered); "
        f"{len(ix.ambiguous):,} strings dropped as ambiguous")
    return ix


#: The accession prefix each strain writes. Used to complete the suffix mapping below.
STRAIN_PREFIXES = {"GT1": "TGGT1", "VEG": "TGVEG"}


def add_strain_accessions(ix: GeneIndex, strain_tsvs: dict, log=print) -> GeneIndex:
    """Register GT1/VEG accessions, which map to ME49 by numeric suffix (see module docstring).

    Two passes, and the second one matters. The first walks the strain TSVs, which is how the
    accessions carrying a trailing letter get registered. But the rule this module documents is that a
    cross-strain accession is accepted *when its suffix exists in the node table* -- it says nothing
    about the accession having been enumerated in a file. Relying on the file alone made resolution
    depend on that file's completeness, and it is not complete: `TGGT1_212960`, `_251570`, `_297960`
    and `_310430` are absent from the 8,637-row GT1 list while all four ME49 genes are in the node
    table, so four genes of the splitCas9 imaging screen silently failed to join. The second pass
    closes that by registering the implied accession for every suffix the node table knows.

    It can only add. Each suffix maps to exactly one ME49 gene, so no existing resolution changes.
    """
    for strain, path in strain_tsvs.items():
        if not os.path.exists(path):
            continue
        t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        n = 0
        for gid in t.gene_id:
            m = re.fullmatch(r"TG[A-Z0-9]{2,6}_(\d{5,6})[A-Za-z]?", str(gid), re.I)
            if not m:
                continue
            target = ix.suffix.get(m.group(1))
            if target is not None:
                ix._add(gid, target, "accession_strain")
                n += 1
        log(f"identity: {n:,} {strain} accessions mapped to ME49 by suffix")
    # Second pass: complete the convention for every suffix in the node table, whether or not the
    # strain list happened to enumerate it.
    for strain, prefix in STRAIN_PREFIXES.items():
        added = 0
        for suffix, target in ix.suffix.items():
            acc = f"{prefix}_{suffix}"
            if norm(acc) not in ix.lookup:
                ix._add(acc, target, "accession_strain")
                added += 1
        log(f"identity: {added:,} further {strain} accessions implied by suffix")
    return ix
