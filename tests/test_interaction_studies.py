#!/usr/bin/env python3
"""The 97 tagged-protein studies, and the reading of the documents they publish their data in.

The single most important thing this module does is refuse to overclaim. A supplement is usually the
paper's COMPLETE quantification table rather than its hit list -- median 754 genes per parsed study, the
largest 7,866, essentially the whole proteome -- so what comes out is *membership*, and calling it
interaction would manufacture roughly fifty thousand false edges.

The document readers are tested for the failure that produces wrong data rather than an error: without
`-layout`, pdftotext reflows a table into prose and columns interleave, so every extracted row is a
blend of two different rows, and the output looks perfectly fine.
"""
from __future__ import annotations

import json
import os
import sys
import zipfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import interaction_studies as IS  # noqa: E402


# --------------------------------------------------------------------------- reading documents
def test_pdftotext_is_asked_for_layout_first(monkeypatch, tmp_path):
    """Without -layout a two-column table reflows into prose and the columns interleave. The corruption
    is invisible unless you read the output."""
    import subprocess
    seen = []

    class R:
        returncode = 0
        stdout = b"TGME49_200010  some protein"

    def fake_run(args, **kw):
        seen.append(args)
        return R()

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert "TGME49_200010" in IS._pdf_text(str(tmp_path / "x.pdf"))
    assert "-layout" in seen[0]


def test_pdf_extraction_falls_back_without_layout_if_the_first_attempt_is_empty(monkeypatch, tmp_path):
    import subprocess
    calls = []

    class Empty:
        returncode = 0
        stdout = b"   "

    class Good:
        returncode = 0
        stdout = b"TGME49_200010"

    def fake_run(args, **kw):
        calls.append(args)
        return Empty() if "-layout" in args else Good()

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert "TGME49_200010" in IS._pdf_text(str(tmp_path / "x.pdf"))
    assert len(calls) == 2


def test_a_missing_pdftotext_yields_empty_text_rather_than_raising(monkeypatch, tmp_path):
    import subprocess

    def boom(*a, **k):
        raise FileNotFoundError("pdftotext")

    monkeypatch.setattr(subprocess, "run", boom)
    assert IS._pdf_text(str(tmp_path / "x.pdf")) == ""


def test_docx_cell_and_row_boundaries_survive_tag_stripping(tmp_path):
    """Strip the tags without converting boundaries first and every cell in a row runs together, which
    turns a gene/description table into one unusable line."""
    p = tmp_path / "x.docx"
    xml = ("<w:document><w:body>"
           "<w:tr><w:tc><w:p><w:t>TGME49_200010</w:t></w:p></w:tc>"
           "<w:tc><w:p><w:t>hypothetical</w:t></w:p></w:tc></w:tr>"
           "<w:tr><w:tc><w:p><w:t>TGME49_200020</w:t></w:p></w:tc></w:tr>"
           "</w:body></w:document>")
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("word/document.xml", xml)
    text = IS._docx_text(str(p))
    lines = [l for l in text.splitlines() if l.strip()]
    assert any("TGME49_200010" in l and "hypothetical" in l for l in lines)
    assert not any("TGME49_200010" in l and "TGME49_200020" in l for l in lines)


def test_a_corrupt_docx_yields_empty_text(tmp_path):
    p = tmp_path / "x.docx"
    p.write_text("this is not a zip")
    assert IS._docx_text(str(p)) == ""


def test_a_docx_without_a_document_part_yields_empty_text(tmp_path):
    p = tmp_path / "x.docx"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("other.xml", "<x/>")
    assert IS._docx_text(str(p)) == ""


def test_spreadsheets_and_delimited_text_are_flattened(tmp_path):
    xlsx = tmp_path / "a.xlsx"
    pd.DataFrame({"id": ["TGME49_200010"]}).to_excel(xlsx, index=False)
    assert "TGME49_200010" in IS._read_any(str(xlsx))

    csv = tmp_path / "b.csv"
    csv.write_text("id\nTGME49_200020\n")
    assert "TGME49_200020" in IS._read_any(str(csv))


def test_an_unreadable_format_yields_empty_text_rather_than_raising(tmp_path):
    p = tmp_path / "x.zip"
    p.write_text("nope")
    assert IS._read_any(str(p)) == ""
    broken = tmp_path / "b.xlsx"
    broken.write_text("not a workbook")
    assert IS._read_any(str(broken)) == ""


# --------------------------------------------------------------------------- parsing studies
def _study(root, pmid, method="proximity", files=None, **meta):
    d = root / method / str(pmid)
    d.mkdir(parents=True)
    json.dump({"pmid": pmid, "method": method, **meta}, open(d / "META.json", "w"))
    for name, content in (files or {}).items():
        (d / name).write_text(content)
    return d


def test_membership_is_extracted_from_a_supplement(tmp_path):
    _study(tmp_path, 111, files={"S1.csv": "id\nTGME49_200010\nTGME49_200020\n"})
    members, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert set(members.gene_id) == {"TGME49_200010", "TGME49_200020"}
    assert studies.loc[0, "n_genes"] == 2
    assert bool(studies.loc[0, "parsed"])


def test_the_method_is_recorded_as_bioid_or_ipms(tmp_path):
    _study(tmp_path, 111, method="proximity", files={"S1.csv": "id\nTGME49_200010\n"})
    _study(tmp_path, 222, method="pulldown", files={"S1.csv": "id\nTGME49_200010\n"})
    _, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert set(studies.method) == {"BioID", "IPMS"}


def test_a_study_with_no_files_is_catalogued_as_unparsed(tmp_path):
    """40 of 97 are in exactly this state, and they must appear in the catalogue rather than vanish
    from it -- otherwise the coverage figure silently improves."""
    _study(tmp_path, 111)
    members, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert members.empty
    assert len(studies) == 1
    assert not bool(studies.loc[0, "parsed"])
    assert studies.loc[0, "n_files"] == 0


def test_formats_that_could_not_be_read_are_named(tmp_path):
    """So a future session knows whether "unparsed" means "no data" or "a format we skipped"."""
    _study(tmp_path, 111, files={"S1.zip": "x", "S2.csv": "id\nTGME49_200010\n"})
    _, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert studies.loc[0, "unreadable_formats"] == ".zip"


def test_which_file_each_gene_came_from_is_kept(tmp_path):
    """Curating membership into interactions later needs to know which sheet a gene was found in."""
    _study(tmp_path, 111, files={"S1.csv": "id\nTGME49_200010\n",
                                 "S2.csv": "id\nTGME49_200010\n"})
    members, _ = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert members.loc[0, "n_files"] == 2
    assert members.loc[0, "files"] == "S1.csv;S2.csv"


def test_accessions_are_routed_through_the_identity_layer(tmp_path):
    """TGGT1_ accessions turned out MORE common in these files than TGME49_."""
    _study(tmp_path, 111, files={"S1.csv": "id\nTGGT1_200010\n"})
    members, _ = IS.parse_studies(str(tmp_path), resolve=lambda a: a.replace("TGGT1", "TGME49"),
                                  log=lambda *_: None)
    assert list(members.gene_id) == ["TGME49_200010"]


def test_without_a_resolver_only_me49_accessions_are_kept(tmp_path):
    """Conservative by design: an unresolved strain accession is not silently treated as an ME49 gene."""
    _study(tmp_path, 111, files={"S1.csv": "id\nTGGT1_200010\nTGME49_200020\n"})
    members, _ = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert list(members.gene_id) == ["TGME49_200020"]


def test_an_empty_root_yields_empty_frames(tmp_path):
    members, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert members.empty and studies.empty


def test_a_file_that_reads_as_nothing_is_skipped(tmp_path):
    _study(tmp_path, 111, files={"S1.txt": ""})
    members, studies = IS.parse_studies(str(tmp_path), log=lambda *_: None)
    assert members.empty
    assert studies.loc[0, "n_tables"] == 1


# --------------------------------------------------------------------------- baits
def test_a_bait_guessed_from_the_title_is_labelled_a_guess():
    """`bait_confidence` exists so nothing downstream can treat a heuristic as an annotation."""
    studies = pd.DataFrame({"title": ["IMC29 Plays an Important Role in Daughter Budding"]})
    out = IS.guess_baits(studies, {"IMC29": "TGME49_243200"}, log=lambda *_: None)
    assert out.loc[0, "bait_gene"] == "TGME49_243200"
    assert out.loc[0, "bait_confidence"] == "title"


def test_a_title_naming_no_known_symbol_leaves_the_bait_empty():
    studies = pd.DataFrame({"title": ["A study of something else entirely"]})
    out = IS.guess_baits(studies, {"IMC29": "TGME49_243200"}, log=lambda *_: None)
    assert out.loc[0, "bait_gene"] is None
    assert out.loc[0, "bait_confidence"] == ""


def test_hyphenated_symbols_in_a_title_are_matched():
    studies = pd.DataFrame({"title": ["the disulfide isomerase PDI-8 interactome"]})
    out = IS.guess_baits(studies, {"PDI8": "TGME49_211680"}, log=lambda *_: None)
    assert out.loc[0, "bait_gene"] == "TGME49_211680"


def test_guessing_baits_on_nothing_returns_nothing():
    assert IS.guess_baits(pd.DataFrame(), {}, log=lambda *_: None).empty


# --------------------------------------------------------------------------- per-gene counts
def test_studies_are_counted_per_gene_and_deduplicated_by_pmid():
    members = pd.DataFrame({"pmid": ["1", "1", "2"], "gene_id": ["g1", "g1", "g1"]})
    out = IS.study_gene_counts(members, ["g1", "g2"])
    assert list(out) == [2, 0]


def test_a_gene_named_by_no_study_counts_zero_not_missing():
    """Zero studies is a measurement about attention; missing would read as "not looked at"."""
    out = IS.study_gene_counts(pd.DataFrame(), ["g1", "g2"])
    assert list(out) == [0, 0]


# --------------------------------------------------------------------------- curated host targets
def _host_table(tmp_path, rows):
    pd.DataFrame(rows).to_csv(tmp_path / "known_host_parasite_interactions.csv", index=False)


def test_a_missing_curated_table_is_reported_and_yields_nothing(tmp_path):
    msgs = []
    assert IS.host_interactions(str(tmp_path), log=msgs.append).empty
    assert any("no curated host-parasite table" in m for m in msgs)


def test_only_toxoplasma_rows_are_kept(tmp_path):
    """The curated table covers several parasites; the others belong to other projects."""
    _host_table(tmp_path, [
        {"organism": "Toxoplasma gondii", "parasite_protein": "ROP16", "host_target": "STAT3"},
        {"organism": "Plasmodium falciparum", "parasite_protein": "PfEMP1", "host_target": "CD36"}])
    out = IS.host_interactions(str(tmp_path), resolve=lambda s: f"TGME49_{s}",
                               log=lambda *_: None)
    assert list(out.parasite_protein) == ["ROP16"]


def test_parasite_symbols_go_through_the_identity_layer(tmp_path):
    """The parasite side is given by symbol (ROP16, GRA24), not by accession."""
    _host_table(tmp_path, [{"organism": "T. gondii", "parasite_protein": "GRA-24",
                            "host_target": "p38"}])
    seen = []

    def resolve(s):
        seen.append(s)
        return "TGME49_230180"

    out = IS.host_interactions(str(tmp_path), resolve=resolve, log=lambda *_: None)
    assert seen == ["GRA24"], "punctuation and spacing are normalized before lookup"
    assert out.loc[0, "gene_id"] == "TGME49_230180"


def test_unresolved_symbols_are_dropped_and_named(tmp_path):
    """Named rather than only counted, so a missing symbol can be added to the identity layer instead
    of being rediscovered."""
    _host_table(tmp_path, [{"organism": "T. gondii", "parasite_protein": "ROP16",
                            "host_target": "STAT3"},
                           {"organism": "T. gondii", "parasite_protein": "MYSTERY1",
                            "host_target": "X"}])
    msgs = []
    out = IS.host_interactions(str(tmp_path),
                               resolve=lambda s: None if s == "MYSTERY1" else "TGME49_1",
                               log=msgs.append)
    assert list(out.parasite_protein) == ["ROP16"]
    assert any("MYSTERY1" in m for m in msgs)


def test_without_a_resolver_nothing_is_claimed(tmp_path):
    _host_table(tmp_path, [{"organism": "T. gondii", "parasite_protein": "ROP16",
                            "host_target": "STAT3"}])
    assert IS.host_interactions(str(tmp_path), log=lambda *_: None).empty


def test_read_any_dispatches_pdf_and_docx_by_extension(tmp_path, monkeypatch):
    """65 of 97 studies publish their hit list only as PDF or DOCX, so these two branches carry most of
    what is recoverable at all."""
    monkeypatch.setattr(IS, "_pdf_text", lambda p, **k: "TGME49_200010 from pdf")
    monkeypatch.setattr(IS, "_docx_text", lambda p: "TGME49_200020 from docx")
    (tmp_path / "a.pdf").write_text("")
    (tmp_path / "b.DOCX").write_text("")
    assert "from pdf" in IS._read_any(str(tmp_path / "a.pdf"))
    assert "from docx" in IS._read_any(str(tmp_path / "b.DOCX"))
