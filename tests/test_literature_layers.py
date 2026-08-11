#!/usr/bin/env python3
"""Regression tests for the identity / corpus / literature layers.

Every test here pins a bug that was real. The literature layer is the part of starplast most able to be
confidently wrong -- a false match does not crash anything, it just quietly hands a gene attention it
never had -- so the guards against that are tested rather than trusted.

Run: pytest tests/ -q   (no display, no network, no machine-local corpora required)
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import corpus, identity, literature  # noqa: E402

GENES = ["TGME49_208830", "TGME49_205250", "TGME49_214320", "TGME49_259720", "TGME49_321440"]


@pytest.fixture
def idtable(tmp_path):
    """A miniature ToxoDB identity table covering the interesting cases."""
    p = tmp_path / "identity.tsv"
    pd.DataFrame([
        # GRA16: an ordinary symbol carrying a digit, plus a real previous accession
        ("TGME49_208830", "GRA16", "Previous IDs: TGME49_008830;TgTwinScan_1234", "dense granule"),
        ("TGME49_205250", "ROP18", "Previous IDs: TGME49_005250", "rhoptry kinase"),
        # GT1: a symbol that is also a strain name and an accession prefix
        ("TGME49_214320", "GT1", "N/A", "glucose transporter"),
        # REMIND: a symbol that is also an English word
        ("TGME49_259720", "REMIND", "N/A", "regulator of membrane-interacting domain"),
        # SANT: a symbol that ASCII tokenisation can carve out of "Santé"
        ("TGME49_321440", "SANT", "N/A", "SWI2/SNF2 ISWI-like"),
    ], columns=["gene_id", "gene_name", "previous_ids", "product"]).to_csv(p, sep="\t", index=False)
    return str(p)


@pytest.fixture
def ix(idtable, tmp_path):
    strain = tmp_path / "gt1.tsv"
    pd.DataFrame({"gene_id": ["TGGT1_208830"], "gene_name": ["N/A"]}).to_csv(
        strain, sep="\t", index=False)
    i = identity.build_index(GENES, idtable, log=lambda *a: None)
    identity.add_strain_accessions(i, {"GT1": str(strain)}, log=lambda *a: None)
    return i


def genes(ix, text):
    return {g for g, _ in ix.find(text)}


def kinds(ix, text):
    return {(g, k) for g, k in ix.find(text)}


# --------------------------------------------------------------------------- identity: resolution
def test_current_accession_resolves(ix):
    assert kinds(ix, "we deleted TGME49_208830 here") == {("TGME49_208830", "accession")}


def test_previous_accession_resolves_to_current_gene(ix):
    """Old ToxoDB ids are still cited. They used to match the regex and then be dropped as unknown."""
    assert kinds(ix, "the locus TGME49_008830") == {("TGME49_208830", "accession_prev")}


def test_strain_accession_maps_by_suffix(ix):
    assert kinds(ix, "in the GT1 strain, TGGT1_208830 is intact") >= \
        {("TGME49_208830", "accession_strain")}


def test_symbol_and_tg_prefixed_alias(ix):
    assert genes(ix, "GRA16 binds") == {"TGME49_208830"}
    assert genes(ix, "TgGRA16 binds") == {"TGME49_208830"}


def test_hyphenated_symbol_normalises(ix):
    assert genes(ix, "GRA-16 was detected") == {"TGME49_208830"}


# --------------------------------------------------------------------------- identity: precision
def test_accession_prefix_is_not_read_as_a_symbol(ix):
    """TGGT1_208830 contains the token 'TGGT1'; TGME49_214320's symbol is 'GT1'."""
    assert "TGME49_214320" not in genes(ix, "the allele TGGT1_208830 was sequenced")


def test_strain_name_is_not_a_gene(ix):
    """'the GT1 strain' is the commonest phrase in this literature; it is not a glucose transporter."""
    assert "TGME49_214320" not in genes(ix, "we compared the GT1 strain with ME49")


def test_lowercase_english_word_is_not_a_gene(ix):
    """REMIND, HOOK, CLAMP, CLIP and SPARK are all real symbols and all real words."""
    assert genes(ix, "these results remind us of earlier work") == set()
    assert genes(ix, "the REMIND domain negatively regulates binding") == {"TGME49_259720"}


def test_accented_word_is_not_truncated_into_a_symbol(ix):
    """An ASCII-only token class carved 'Sant' out of French 'Santé' and gave SANT 23 abstracts."""
    assert genes(ix, "Institut de Santé publique, Paris") == set()
    assert genes(ix, "the SANT domain was resolved") == {"TGME49_321440"}


def test_ambiguity_is_recorded_not_guessed(idtable, tmp_path):
    """A string claimed by two genes at the same tier is withdrawn, and kept for reporting."""
    p = tmp_path / "dup.tsv"
    pd.DataFrame([("TGME49_208830", "DUPE1", "N/A", "x"),
                  ("TGME49_205250", "DUPE1", "N/A", "y")],
                 columns=["gene_id", "gene_name", "previous_ids", "product"]).to_csv(
        p, sep="\t", index=False)
    i = identity.build_index(GENES, str(p), log=lambda *a: None)
    assert identity.norm("DUPE1") not in i.lookup
    assert i.ambiguous[identity.norm("DUPE1")] == {"TGME49_208830", "TGME49_205250"}


def test_current_accession_survives_another_genes_previous_id(idtable, tmp_path):
    """Cross-tier collisions resolve to the more specific identifier instead of withdrawing both."""
    p = tmp_path / "clash.tsv"
    pd.DataFrame([("TGME49_208830", "N/A", "N/A", "x"),
                  ("TGME49_205250", "N/A", "Previous IDs: TGME49_208830", "y")],
                 columns=["gene_id", "gene_name", "previous_ids", "product"]).to_csv(
        p, sep="\t", index=False)
    i = identity.build_index(GENES, str(p), log=lambda *a: None)
    assert i.lookup[identity.norm("TGME49_208830")] == ("TGME49_208830", "accession")


# --------------------------------------------------------------------------- corpus
JATS = """<?xml version="1.0"?><article><front><article-meta>
<article-id pub-id-type="pmid">12345678</article-id>
<article-id pub-id-type="pmcid">PMC999</article-id>
<title-group><article-title>GRA16 and the host</article-title></title-group>
<abstract><p>We study TgGRA16.</p></abstract>
<pub-date><year>2021</year></pub-date>
</article-meta></front><body>
<sec><p>The protein ROP18 was measured in a long enough paragraph to be kept by the parser.</p>
<fig><caption><p>Figure 1. Localisation of GRA16 in infected cells.</p></caption></fig></sec>
</body><back><ref-list><ref><mixed-citation>Smith et al. TGME49_205250 is a kinase.
</mixed-citation></ref></ref-list></back></article>"""


@pytest.fixture
def jats_doc(tmp_path):
    p = tmp_path / "PMC999.xml"
    p.write_text(JATS)
    return corpus.parse_jats(str(p))


def test_jats_parses_ids_and_sections(jats_doc):
    assert jats_doc.doc_id == "pmc:PMC999"
    assert jats_doc.pmid == "12345678" and jats_doc.year == "2021"
    assert {s.kind for s in jats_doc.sections} == {"title", "abstract", "body", "caption"}


def test_reference_list_is_excluded(jats_doc):
    """A bibliography names every gene its cited papers named; crediting those inflates coverage."""
    assert all("Smith" not in s.text for s in jats_doc.sections)


def test_body_paragraphs_are_separate_sections(jats_doc):
    """Co-mention units are paragraphs, so paragraphs must not be concatenated."""
    assert sum(1 for s in jats_doc.sections if s.kind == "body") == 1


# --------------------------------------------------------------------------- literature
def _doc(doc_id, source, texts, kind="abstract", pmid=None):
    return corpus.Document(doc_id, source, pmid, None,
                           tuple(corpus.Section(kind, t) for t in texts))


def test_union_document_count_deduplicates_by_pmid():
    """An open-access paper is usually also a PubMed record; the union must not count it twice."""
    m = pd.DataFrame({"gene_id": ["g1", "g1"], "doc_id": ["pmid:7", "pmc:PMC7"],
                      "pmid": ["7", "7"], "source": ["abstract", "fulltext"],
                      "match_kind": ["symbol"] * 2, "section": ["abstract", "body"],
                      "n_hits": [1, 1], "tier": ["symbol"] * 2})
    cov = literature.coverage(m).set_index(["source", "tier"])
    assert cov.loc[("union", "any"), "documents"] == 1


def test_scan_builds_a_tidy_auditable_table(ix):
    docs = [_doc("pmid:1", "abstract", ["GRA16 interacts with ROP18"]),
            _doc("pmid:2", "abstract", ["GRA16 alone"])]
    m, co, meta = literature.scan(docs, ix, log=lambda *a: None)
    assert set(m.columns) == {"gene_id", "doc_id", "pmid", "source", "match_kind", "section",
                              "n_hits", "tier"}
    assert m[m.gene_id == "TGME49_208830"].doc_id.nunique() == 2
    assert co["abstract"][("TGME49_205250", "TGME49_208830")] == 1


def test_sources_are_kept_separate(ix):
    docs = [_doc("pmid:1", "abstract", ["GRA16 here"]),
            _doc("pmc:1", "fulltext", ["ROP18 here"], kind="body")]
    m, _, _ = literature.scan(docs, ix, log=lambda *a: None)
    cov = literature.coverage(m).set_index(["source", "tier"]).genes
    assert cov[("abstract", "any")] == 1 and cov[("fulltext", "any")] == 1
    assert cov[("union", "any")] == 2


def test_list_units_are_excluded_from_comention(ix):
    """A unit naming many genes evidences set membership, not a pairwise relation."""
    many = " ".join(f"TGME49_{i}" for i in [208830, 205250, 214320, 259720, 321440])
    m, co, meta = literature.scan([_doc("pmid:1", "abstract", [many])], ix,
                                  max_genes_per_unit=3, log=lambda *a: None)
    assert co["abstract"] == {}
    assert meta["n_units_skipped"]["abstract"] == 1
    assert m.gene_id.nunique() == 5          # still counts for coverage


def test_attention_correction_denominator_uses_units(ix):
    """Two genes co-mentioned exactly as often as independence predicts get a residual near zero."""
    docs = ([_doc(f"pmid:a{i}", "abstract", ["GRA16 alone"]) for i in range(8)] +
            [_doc(f"pmid:b{i}", "abstract", ["ROP18 alone"]) for i in range(8)] +
            [_doc("pmid:c1", "abstract", ["GRA16 and ROP18"]),
             _doc("pmid:c2", "abstract", ["GRA16 and ROP18"])])
    _, co, meta = literature.scan(docs, ix, log=lambda *a: None)
    edges = literature.comention_edges(co["abstract"], meta, "abstract", min_count=2)
    assert len(edges) == 1
    g1, g2, raw, resid = edges[0]
    assert raw == 2.0
    # n1 = n2 = 10 over N = 18 units -> expected 5.56, observed 2 -> clearly negative
    assert resid < 0


def _mention(gene, paper, section, source="fulltext"):
    return {"gene_id": gene, "doc_id": f"pmc:{paper}", "pmid": paper, "source": source,
            "match_kind": "symbol", "section": section, "n_hits": 1, "tier": "symbol"}


def test_attention_depth_reads_tiers_off_document_structure():
    m = pd.DataFrame([_mention("g_title", "1", "title"),
                      _mention("g_abs", "2", "abstract"),
                      _mention("g_body", "3", "body"),
                      _mention("g_cap", "4", "caption")])
    d = literature.attention_depth(m)
    assert d.loc["g_title", "attention_depth"] == "focal"
    assert d.loc["g_abs", "attention_depth"] == "substantive"
    assert d.loc["g_body", "attention_depth"] == "incidental"
    assert d.loc["g_cap", "attention_depth"] == "incidental"


def test_attention_depth_takes_the_strongest_tier_per_paper():
    """A gene in a paper's title and its body is focal for that paper, counted once."""
    m = pd.DataFrame([_mention("g", "1", "title"), _mention("g", "1", "body"),
                      _mention("g", "1", "caption")])
    d = literature.attention_depth(m)
    assert d.loc["g", "n_papers_focal"] == 1
    assert d.loc["g", "n_papers_incidental"] == 0


def test_attention_depth_separates_listed_from_studied():
    """The distinction the tiering exists for: a gene in 40 hit tables is named, not studied."""
    listed = [_mention("g_listed", str(i), "caption") for i in range(40)]
    studied = [_mention("g_studied", "99", "title")]
    d = literature.attention_depth(pd.DataFrame(listed + studied))
    assert d.loc["g_listed", "attention_depth"] == "incidental"
    assert d.loc["g_listed", "n_papers_incidental"] == 40
    assert d.loc["g_studied", "attention_depth"] == "focal"


def test_attention_depth_deduplicates_papers_across_sources():
    """The same paper as a PubMed record and an open-access full text is one paper."""
    m = pd.DataFrame([_mention("g", "7", "abstract", source="abstract"),
                      _mention("g", "7", "body", source="fulltext")])
    d = literature.attention_depth(m)
    assert d.loc["g", "n_papers_substantive"] == 1
    assert d.loc["g", "n_papers_incidental"] == 0


def _edges(**kw):
    """Build an edges dict of the shape build_graph passes around: label -> (a, b, w, r)."""
    import numpy as np
    out = {}
    for k, ps in kw.items():
        if not ps:
            continue
        a = np.array([p[0] for p in ps]); b = np.array([p[1] for p in ps])
        w = np.ones(len(ps))
        out[k] = (a, b, w, w)
    return out


def _nodes(n=6):
    return pd.DataFrame({"gene_id": [f"g{i}" for i in range(n)],
                         "attention_depth": ["focal"] * n})


def test_structural_hole_needs_two_independent_families():
    """orthogroup + domain is one fact, not two: paralogs almost always share domains.

    Before this collapse, 53 of the 66 sharpest candidates were nothing but paralogy.
    """
    from starplast import build_graph as B
    e = _edges(orthogroup=[(0, 1)], domain=[(0, 1)])
    assert B.structural_holes(e, _nodes()) is None

    e = _edges(coexpression=[(0, 1)], cofitness=[(0, 1)])
    hole = B.structural_holes(e, _nodes())
    assert hole is not None and len(hole[0]) == 1


def test_homology_cannot_be_one_of_the_two_legs():
    """Paralogs co-express *because* they are paralogs, so homology corroborates but never qualifies.

    Allowing it admitted 291 expression+homology pairs, 76% of them same-orthogroup, and put one
    protein family at the top of every ranking.
    """
    from starplast import build_graph as B
    assert B.structural_holes(_edges(coexpression=[(0, 1)], orthogroup=[(0, 1)]), _nodes()) is None
    assert B.structural_holes(_edges(cofitness=[(0, 1)], domain=[(0, 1)]), _nodes()) is None
    # both phenotypes present: homology may ride along and raises the weight
    hole = B.structural_holes(
        _edges(coexpression=[(0, 1)], cofitness=[(0, 1)], domain=[(0, 1)]), _nodes())
    assert hole is not None and hole[2][0] == 3.0


def test_structural_hole_requires_the_literature_to_be_silent():
    """A pair discussed anywhere -- any abstract, any open-access paragraph -- is not a hole."""
    from starplast import build_graph as B
    base = dict(coexpression=[(0, 1)], cofitness=[(0, 1)])
    assert B.structural_holes(_edges(**base), _nodes()) is not None
    assert B.structural_holes(_edges(comention=[(0, 1)], **base), _nodes()) is None
    assert B.structural_holes(_edges(comention_ft=[(0, 1)], **base), _nodes()) is None


def test_structural_hole_counts_per_gene():
    from starplast import build_graph as B
    n = _nodes()
    e = _edges(coexpression=[(0, 1), (0, 2)], cofitness=[(0, 1), (0, 2)])
    B.structural_holes(e, n)
    assert n.n_holes.tolist() == [2, 1, 1, 0, 0, 0]


def test_unwritten_interaction_is_measured_binding_the_literature_missed():
    """A stronger claim than a hole: the interaction was observed, not predicted."""
    from starplast import build_graph as B
    n = _nodes()
    e = _edges(xlms=[(0, 1), (2, 3)], comention=[(2, 3)])
    a, b, w, _ = B.unwritten_interactions(e, n)
    assert list(zip(a, b)) == [(0, 1)], "a co-mentioned pair is not unwritten"


def test_unwritten_interaction_merges_xlms_and_ipms_explicitly():
    """The merge is labelled, and both source types stay separately toggleable."""
    from starplast import build_graph as B
    e = _edges(xlms=[(0, 1)], ip_ms=[(2, 3)])
    a, b, _, _ = B.unwritten_interactions(e, _nodes())
    assert sorted(zip(a, b)) == [(0, 1), (2, 3)]
    assert "xlms" in e and "ip_ms" in e


def test_unwritten_interaction_needs_measured_evidence():
    """Correlational edge types alone never qualify -- that is what structural_hole is for."""
    from starplast import build_graph as B
    assert B.unwritten_interactions(_edges(coexpression=[(0, 1)], cofitness=[(0, 1)]),
                                    _nodes()) is None


def test_compartment_alone_never_makes_a_hole():
    """Sharing one of 27 hyperLOPIT classes is too unspecific, and tracks abundance."""
    from starplast import build_graph as B
    assert B.structural_holes(_edges(compartment=[(0, 1)], coexpression=[(0, 1)]), _nodes()) is None


def test_attention_correction_rewards_unexpected_pairs(ix):
    """Two rarely-studied genes that always appear together score positive despite a small raw count."""
    docs = ([_doc(f"pmid:x{i}", "abstract", ["nothing here"]) for i in range(200)] +
            [_doc("pmid:y1", "abstract", ["GRA16 and ROP18"]),
             _doc("pmid:y2", "abstract", ["GRA16 and ROP18"])])
    _, co, meta = literature.scan(docs, ix, log=lambda *a: None)
    edges = literature.comention_edges(co["abstract"], meta, "abstract", min_count=2)
    assert edges and edges[0][3] > 0
