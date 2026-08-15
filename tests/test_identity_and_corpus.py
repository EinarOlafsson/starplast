#!/usr/bin/env python3
"""Gene identity resolution, and the document stream it runs over.

Identity is the layer every other join depends on, and its failures are silent by construction: an
accession that resolves to nothing produces a dataset with no coverage, which looks exactly like a
dataset nobody measured. A 2019 in vivo screen using pre-2012 ids contributed 0 of 8,140 rows until it
was routed through here, then 168.

The symbol matching is where it gets dangerous. HOOK, CLAMP, CLIP, SPARK and REMIND are all real
*Toxoplasma* gene symbols and all ordinary English words, and GT1 is both a strain designation and the
symbol of a glucose transporter. Every guard below exists because removing it fills the corpus with
prose.
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import corpus as CO  # noqa: E402
from starplast import identity as ID  # noqa: E402


def _identity_tsv(tmp_path, rows):
    p = tmp_path / "toxodb_identity.tsv"
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return str(p)


# --------------------------------------------------------------------------- normalization
def test_normalization_is_case_and_punctuation_insensitive():
    assert ID.norm("GRA-16") == ID.norm("gra16") == ID.norm(" GRA 16 ")


# --------------------------------------------------------------------------- what may be a symbol
def test_a_digit_free_symbol_must_be_upper_case():
    """HOOK, CLAMP, CLIP, SPARK and REMIND are real symbols and ordinary English words. Matched
    case-insensitively they fill the corpus with prose: SANT went from 24 mentions to 1, REMIND from
    12 to 3."""
    assert ID.matchable_symbol("HOOK")
    assert not ID.matchable_symbol("hook")


def test_a_symbol_containing_a_digit_is_safe_either_way():
    assert ID.matchable_symbol("GRA16")
    assert ID.matchable_symbol("gra16")


def test_a_symbol_that_is_too_short_is_refused():
    assert not ID.matchable_symbol("AB")


def test_strain_designations_are_blocked_outright():
    """GT1 is both a strain name and the symbol of TGME49_214320, a glucose transporter. Every
    occurrence in the literature is the strain."""
    for s in ID.BLOCKED_SYMBOLS:
        assert not ID.matchable_symbol(s)


# --------------------------------------------------------------------------- building the index
def test_a_missing_identity_table_is_reported_with_what_to_run(tmp_path):
    msgs = []
    ix = ID.build_index(["TGME49_200010"], str(tmp_path / "absent.tsv"), log=msgs.append)
    assert any("fetch_names" in m for m in msgs)
    assert ix.lookup


def test_a_table_missing_required_columns_is_an_explicit_error(tmp_path):
    p = _identity_tsv(tmp_path, [{"something": "x"}])
    with pytest.raises(ValueError, match="missing columns"):
        ID.build_index(["TGME49_200010"], p, log=lambda *_: None)


def test_an_accession_resolves_to_itself(tmp_path):
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "",
                                  "previous_ids": "", "product": "hypothetical"}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    assert ix.lookup[ID.norm("TGME49_200010")][0] == "TGME49_200010"


def test_a_previous_accession_resolves_forward(tmp_path):
    """Papers cite whatever was current when they were written."""
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "",
                                  "previous_ids": "TGME49_008830", "product": ""}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    hit = ix.lookup[ID.norm("TGME49_008830")]
    assert hit[0] == "TGME49_200010" and hit[1] == "accession_prev"


def test_a_legacy_toxogenechip_model_resolves_for_structured_dataset_joins(tmp_path):
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "",
                                  "previous_ids": "1.m00014", "product": ""}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    assert ix.lookup[ID.norm("1.m00014")] == ("TGME49_200010", "accession_prev")


def test_a_symbol_resolves_and_is_tagged_as_a_symbol(tmp_path):
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "GRA16",
                                  "previous_ids": "", "product": ""}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    assert ix.lookup[ID.norm("GRA16")] == ("TGME49_200010", "symbol")


def test_a_stronger_kind_of_evidence_wins_over_a_weaker_one(tmp_path):
    """An accession beats a symbol for the same string, because it is a stronger claim."""
    p = _identity_tsv(tmp_path, [
        {"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200020", "gene_name": "TGME49_200010", "previous_ids": "", "product": ""}])
    ix = ID.build_index(["TGME49_200010", "TGME49_200020"], p, log=lambda *_: None)
    assert ix.lookup[ID.norm("TGME49_200010")][0] == "TGME49_200010"


def test_a_string_claimed_by_two_genes_at_the_same_strength_is_dropped(tmp_path):
    """Ambiguity is recorded, not guessed. Picking one would invent an assignment."""
    p = _identity_tsv(tmp_path, [
        {"gene_id": "TGME49_200010", "gene_name": "SAME1", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200020", "gene_name": "SAME1", "previous_ids": "", "product": ""}])
    ix = ID.build_index(["TGME49_200010", "TGME49_200020"], p, log=lambda *_: None)
    assert ID.norm("SAME1") not in ix.lookup
    assert ID.norm("SAME1") in ix.ambiguous


def test_genes_absent_from_the_node_table_are_not_registered(tmp_path):
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_999999", "gene_name": "GHOST1",
                                  "previous_ids": "", "product": ""}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    assert ID.norm("GHOST1") not in ix.lookup


def test_the_index_reports_what_it_holds(tmp_path):
    p = _identity_tsv(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "GRA16",
                                  "previous_ids": "TGME49_008830", "product": ""}])
    ix = ID.build_index(["TGME49_200010"], p, log=lambda *_: None)
    stats = ix.stats()
    assert stats["k_symbol"] >= 1 and stats["k_accession_prev"] >= 1
    assert stats["genes"] == 1


# --------------------------------------------------------------------------- finding genes in text
def _index(tmp_path, rows, node_ids=None):
    p = _identity_tsv(tmp_path, rows)
    return ID.build_index(node_ids or [r["gene_id"] for r in rows], p, log=lambda *_: None)


def test_an_accession_is_found_in_running_text(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    assert [h[0] for h in ix.find("we deleted TGME49_200010 in this study")] == ["TGME49_200010"]


def test_a_symbol_inside_an_accession_is_not_matched_again(tmp_path):
    """Without excluding the accession's own span, the digits and letters inside it re-match as a
    symbol and one mention is counted twice."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "ME49", "previous_ids": "",
                            "product": ""}])
    assert len(list(ix.find("TGME49_200010 was tagged"))) == 1


def test_a_lower_case_english_word_that_is_also_a_symbol_is_not_matched(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "HOOK", "previous_ids": "",
                            "product": ""}])
    assert list(ix.find("the parasite uses a hook to attach")) == []
    assert len(list(ix.find("HOOK was tagged"))) == 1


def test_tokenisation_is_unicode_aware(tmp_path):
    """An ASCII-only class truncates French Sante to Sant, which is a real symbol, and manufactures a
    match from nothing."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "SANT", "previous_ids": "",
                            "product": ""}])
    assert list(ix.find("publié dans Santé publique")) == []


def test_finding_in_empty_text_returns_nothing(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "GRA16", "previous_ids": "",
                            "product": ""}])
    assert list(ix.find("")) == []
    assert list(ix.find(None)) == []


# --------------------------------------------------------------------------- strain accessions
def test_strain_accessions_map_to_me49_by_numeric_suffix(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    strain = tmp_path / "gt1.tsv"
    pd.DataFrame({"gene_id": ["TGGT1_200010"], "gene_name": [""]}).to_csv(strain, sep="\t",
                                                                         index=False)
    ID.add_strain_accessions(ix, {"GT1": str(strain)}, log=lambda *_: None)
    assert ix.lookup[ID.norm("TGGT1_200010")] == ("TGME49_200010", "accession_strain")


def test_a_strain_accession_with_no_me49_counterpart_is_not_invented(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    strain = tmp_path / "gt1.tsv"
    pd.DataFrame({"gene_id": ["TGGT1_999999"], "gene_name": [""]}).to_csv(strain, sep="\t",
                                                                         index=False)
    ID.add_strain_accessions(ix, {"GT1": str(strain)}, log=lambda *_: None)
    assert ID.norm("TGGT1_999999") not in ix.lookup


def test_a_missing_strain_table_is_skipped(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    ID.add_strain_accessions(ix, {"GT1": str(tmp_path / "absent.tsv")}, log=lambda *_: None)
    assert ix.lookup


def test_a_malformed_strain_identifier_is_skipped(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    strain = tmp_path / "gt1.tsv"
    pd.DataFrame({"gene_id": ["not an accession"], "gene_name": [""]}).to_csv(strain, sep="\t",
                                                                             index=False)
    ID.add_strain_accessions(ix, {"GT1": str(strain)}, log=lambda *_: None)


# --------------------------------------------------------------------------- the abstract stream
def test_abstracts_become_documents_with_title_and_abstract_sections(tmp_path):
    p = tmp_path / "abs.jsonl"
    p.write_text(json.dumps({"pmid": "123", "year": "2020", "title": "GRA16 does something",
                             "abstract": "We show that GRA16 ..."}) + "\n")
    docs = list(CO.iter_abstracts(str(p)))
    assert len(docs) == 1
    assert [s.kind for s in docs[0].sections] == ["title", "abstract"]
    assert docs[0].pmid == "123" and docs[0].source == "abstract"


def test_a_record_with_neither_title_nor_abstract_is_skipped(tmp_path):
    p = tmp_path / "abs.jsonl"
    p.write_text(json.dumps({"pmid": "123"}) + "\n")
    assert list(CO.iter_abstracts(str(p))) == []


def test_a_malformed_line_does_not_stop_the_stream(tmp_path):
    """One bad line in a 47 MB corpus should cost that record, not the corpus."""
    p = tmp_path / "abs.jsonl"
    p.write_text("{not json\n" + json.dumps({"pmid": "1", "title": "ok"}) + "\n")
    assert len(list(CO.iter_abstracts(str(p)))) == 1


def test_a_missing_corpus_yields_nothing_rather_than_raising(tmp_path):
    assert list(CO.iter_abstracts(str(tmp_path / "absent.jsonl"))) == []


# --------------------------------------------------------------------------- more identity edges
def test_a_symbol_that_is_not_matchable_is_never_registered(tmp_path):
    """A lower-case digit-free name would otherwise enter the index and match prose everywhere."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "hook", "previous_ids": "",
                            "product": ""}])
    assert ID.norm("hook") not in ix.lookup


def test_a_string_already_ambiguous_stays_ambiguous(tmp_path):
    """Once two genes have claimed a string at equal strength, a third claim must not resolve it."""
    ix = _index(tmp_path, [
        {"gene_id": "TGME49_200010", "gene_name": "SAME1", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200020", "gene_name": "SAME1", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200030", "gene_name": "SAME1", "previous_ids": "", "product": ""}])
    assert ID.norm("SAME1") not in ix.lookup
    assert len(ix.ambiguous[ID.norm("SAME1")]) >= 2


def test_a_weaker_claim_does_not_displace_a_stronger_one(tmp_path):
    """An alias must never overwrite an accession, whichever order the table happens to list them in."""
    ix = _index(tmp_path, [
        {"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200020", "gene_name": "", "previous_ids": "TGME49_200010",
         "product": ""}])
    assert ix.lookup[ID.norm("TGME49_200010")] == ("TGME49_200010", "accession")


def test_the_same_gene_claiming_a_string_twice_is_not_an_ambiguity(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "GRA16",
                            "previous_ids": "GRA16", "product": ""}])
    assert ix.lookup[ID.norm("GRA16")][0] == "TGME49_200010"
    assert ID.norm("GRA16") not in ix.ambiguous


def test_an_unregistered_accession_pattern_in_text_resolves_by_suffix(tmp_path):
    """A strain accession the index has never seen still maps to ME49 by its numeric suffix."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    hits = list(ix.find("we used TGVEG_200010 for this"))
    assert [h[0] for h in hits] == ["TGME49_200010"]


# --------------------------------------------------------------------------- JATS full texts
def _jats(tmp_path, body, name="a.xml"):
    p = tmp_path / name
    p.write_text(body, encoding="utf8")
    return str(p)


def test_a_jats_article_becomes_a_sectioned_document(tmp_path):
    p = _jats(tmp_path, """<article>
      <front><article-meta>
        <article-id pub-id-type="pmid">123</article-id>
        <article-id pub-id-type="pmc">PMC9</article-id>
        <title-group><article-title>GRA16 in <italic>Toxoplasma</italic></article-title></title-group>
        <abstract><p>We show GRA16 does something.</p></abstract>
        <pub-date><year>2020</year></pub-date>
      </article-meta></front>
      <body><sec><title>Results</title><p>GRA16 was tagged and localized to the dense granules in every replicate we examined.</p></sec></body>
    </article>""")
    doc = CO.parse_jats(p)
    assert doc is not None
    assert doc.pmid == "123" and doc.source == "fulltext"
    kinds = [s.kind for s in doc.sections]
    assert "title" in kinds and "abstract" in kinds and "body" in kinds


def test_text_scattered_across_inline_tags_is_flattened(tmp_path):
    """JATS puts italics and superscripts mid-sentence, and a gene symbol can straddle them."""
    p = _jats(tmp_path, """<article><front><article-meta>
        <title-group><article-title>The <italic>GRA</italic>16 protein</article-title></title-group>
      </article-meta></front></article>""")
    doc = CO.parse_jats(p)
    assert "GRA16" in doc.sections[0].text.replace(" ", "")


def test_reference_lists_are_removed_before_anything_is_read(tmp_path):
    """Otherwise every gene named in a cited paper's title counts as a mention of this paper's."""
    p = _jats(tmp_path, """<article>
      <front><article-meta>
        <title-group><article-title>A study</article-title></title-group>
      </article-meta></front>
      <body><sec><p>Nothing here.</p></sec></body>
      <back><ref-list><ref><element-citation><article-title>GRA16 elsewhere</article-title>
      </element-citation></ref></ref-list></back>
    </article>""")
    doc = CO.parse_jats(p)
    assert all("GRA16" not in s.text for s in doc.sections)


def test_an_unparseable_file_returns_none(tmp_path):
    assert CO.parse_jats(_jats(tmp_path, "<article><unclosed>")) is None


def test_a_document_with_no_readable_text_is_none_or_empty(tmp_path):
    doc = CO.parse_jats(_jats(tmp_path, "<article></article>"))
    assert doc is None or not doc.sections


def test_caption_text_is_not_counted_twice(tmp_path):
    """Figure and table captions wrap their text in <p>, so an unguarded body scan emits caption text
    once as body and once as caption -- double-counting every gene named in a caption, and the
    co-mention unit it sits in."""
    p = _jats(tmp_path, """<article>
      <front><article-meta><title-group><article-title>A study</article-title></title-group>
      </article-meta></front>
      <body><fig><caption><p>Figure 1. GRA16 localizes to the nucleus in all conditions tested.</p>
      </caption></fig></body></article>""")
    doc = CO.parse_jats(p)
    texts = [s.text for s in doc.sections if "GRA16" in s.text]
    assert len(texts) == 1
    assert [s.kind for s in doc.sections if "GRA16" in s.text] == ["caption"]


def test_very_short_paragraphs_are_skipped(tmp_path):
    """A three-word paragraph is a heading fragment or a table cell, not a claim about a relation."""
    p = _jats(tmp_path, """<article>
      <front><article-meta><title-group><article-title>A study</article-title></title-group>
      </article-meta></front>
      <body><sec><p>See below.</p></sec></body></article>""")
    doc = CO.parse_jats(p)
    assert all(s.kind != "body" for s in doc.sections)


def test_the_publication_year_is_taken_from_the_first_usable_pub_date(tmp_path):
    p = _jats(tmp_path, """<article>
      <front><article-meta>
        <pub-date><year>not-a-year</year></pub-date>
        <pub-date><year>2019</year></pub-date>
        <title-group><article-title>A study</article-title></title-group>
      </article-meta></front></article>""")
    assert CO.parse_jats(p).year == "2019"


def test_a_document_without_a_pmcid_is_identified_by_its_filename(tmp_path):
    """So two files without accessions never collapse into one document."""
    p = _jats(tmp_path, """<article><front><article-meta>
      <title-group><article-title>A study</article-title></title-group>
      </article-meta></front></article>""", name="unique_name.xml")
    assert CO.parse_jats(p).doc_id == "file:unique_name.xml"


def test_full_texts_are_streamed_from_a_directory(tmp_path):
    for i in range(3):
        _jats(tmp_path, f"""<article><front><article-meta>
          <title-group><article-title>Study {i}</article-title></title-group>
          </article-meta></front></article>""", name=f"{i}.xml")
    assert len(list(CO.iter_fulltexts(str(tmp_path)))) == 3
    assert len(list(CO.iter_fulltexts(str(tmp_path), limit=2))) == 2


def test_an_unreadable_file_does_not_stop_the_stream(tmp_path):
    """6,667 files on disk; one truncated download should cost that paper, not the corpus."""
    _jats(tmp_path, "<article><unclosed>", name="bad.xml")
    _jats(tmp_path, """<article><front><article-meta>
      <title-group><article-title>Good one</article-title></title-group>
      </article-meta></front></article>""", name="good.xml")
    assert len(list(CO.iter_fulltexts(str(tmp_path)))) == 1


def test_an_empty_key_or_unknown_gene_registers_nothing(tmp_path):
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    before = len(ix.lookup)
    ix._add("", "TGME49_200010", "symbol")
    ix._add("SOMETHING", "TGME49_NOT_A_GENE", "symbol")
    assert len(ix.lookup) == before


def test_the_same_gene_upgrades_its_own_entry_to_a_stronger_kind(tmp_path):
    """Registered first as an alias and later as an accession, the stronger claim should win --
    whichever order the table happened to list them in."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    ix._add("SOMEKEY", "TGME49_200010", "alias")
    assert ix.lookup[ID.norm("SOMEKEY")][1] == "alias"
    ix._add("SOMEKEY", "TGME49_200010", "accession")
    assert ix.lookup[ID.norm("SOMEKEY")][1] == "accession"
    ix._add("SOMEKEY", "TGME49_200010", "alias")
    assert ix.lookup[ID.norm("SOMEKEY")][1] == "accession", "a weaker claim must not downgrade it"


def test_a_more_specific_identifier_wins_across_tiers(tmp_path):
    """A gene's own current accession must not be withdrawn because another gene once carried that
    string as a previous id."""
    ix = _index(tmp_path, [
        {"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "", "product": ""},
        {"gene_id": "TGME49_200020", "gene_name": "", "previous_ids": "", "product": ""}])
    ix._add("SHARED", "TGME49_200010", "accession_prev")
    ix._add("SHARED", "TGME49_200020", "accession")
    assert ix.lookup[ID.norm("SHARED")] == ("TGME49_200020", "accession")


def test_a_pmcid_in_the_article_metadata_is_used_as_the_document_id(tmp_path):
    p = _jats(tmp_path, """<article><front><article-meta>
      <article-id pub-id-type="pmcid">PMC12345</article-id>
      <title-group><article-title>A study</article-title></title-group>
      </article-meta></front></article>""")
    assert CO.parse_jats(p).doc_id == "pmc:PMC12345"


def test_the_corpus_can_be_counted_without_being_parsed(tmp_path):
    """Used to report corpus size before a build that takes minutes."""
    for i in range(3):
        _jats(tmp_path, "<article/>", name=f"{i}.xml")
    assert CO.count_fulltexts(str(tmp_path)) == 3


def test_an_unregistered_me49_shaped_accession_resolves_to_nothing(tmp_path):
    """The suffix map is a fallback for strain accessions the strain tables missed. An ME49-shaped id
    that is in no table is not a gene, and inventing one would be worse than missing it."""
    ix = _index(tmp_path, [{"gene_id": "TGME49_200010", "gene_name": "", "previous_ids": "",
                            "product": ""}])
    assert list(ix.find("we also looked at TGME49_777777 here")) == []
