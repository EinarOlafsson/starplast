"""A broad audit must preserve missing metadata and use the service's actual query semantics."""
import io
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import audit_dataset_selection as A
from starplast import datasets as D, organisms as O, slots as S
from starplast.dataset_selection import DatasetCandidate


def test_recorded_doi_extraction_does_not_include_supplement_file_paths():
    source = D.Dataset("synthetic", "synthetic", "DNA", "synthetic", "synthetic",
                       organism=O.TOXOPLASMA, url="https://example.org/art%3A10.1234%2Fsynthetic/MediaObjects/file.xlsx")
    assert A._doi(source) == "10.1234/synthetic"
    source = D.Dataset("synthetic", "synthetic", "DNA", "synthetic", "synthetic",
                       organism=O.TOXOPLASMA, citation="No identifier recorded")
    assert A._doi(source) == ""


def test_query_order_is_not_silently_treated_as_free_text(monkeypatch, tmp_path):
    captured = []
    def fetch(request, timeout):
        captured.append(request.full_url)
        return io.BytesIO(json.dumps({"hitCount": 0, "resultList": {"result": []}}).encode())
    monkeypatch.setattr(A.urllib.request, "urlopen", fetch)
    monkeypatch.setattr(A.time, "sleep", lambda *_: None)
    client = A.PublicationClient(tmp_path)
    client.search('TITLE_ABS:"synthetic"', sort="CITED desc")
    client.search('TITLE_ABS:"synthetic"', sort="FIRST_PDATE desc")
    queries = [A.urllib.parse.parse_qs(A.urllib.parse.urlparse(url).query)["query"][0] for url in captured]
    assert queries == ['TITLE_ABS:"synthetic" sort_cited:y', 'TITLE_ABS:"synthetic" sort_date:y']
    client.search('TITLE_ABS:"synthetic"', sort="CITED desc")
    assert len(captured) == 2  # exact primary response reused
    with pytest.raises(ValueError):
        client.search("synthetic", sort="unknown")


def test_failed_queries_remain_failed_instead_of_inventing_zero_citations(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise TimeoutError("synthetic service unavailable")
    monkeypatch.setattr(A.urllib.request, "urlopen", fail)
    monkeypatch.setattr(A.time, "sleep", lambda *_: None)
    response = A.PublicationClient(tmp_path).search("synthetic")
    assert response["status"] == "error" and "TimeoutError" in response["error"]
    assert A._results(response) == []
    saved = json.loads(next((tmp_path / "requests").glob("*.json")).read_text())
    assert saved["response"] == {} and saved["error"]


def test_publication_identity_and_retraction_flags_are_provider_records():
    record = {"source": "MED", "id": "SYNTHETIC", "title": "Synthetic paper", "citedByCount": 0,
              "firstPublicationDate": "2020-01-01", "pubTypeList": {"pubType": ["Retracted Publication"]}}
    publication = A._publication(record)
    assert publication["publication_id"] == "MED:SYNTHETIC" and publication["retracted"]
    assert publication["citations"] == 0
    assert A._publication({"source": "PPR", "id": "synthetic"})["citations"] is None


def test_every_slot_has_an_explicit_organism_scoped_query():
    for slot in S.all_slots():
        query = A.query_for(slot)
        assert "TITLE_ABS:" in query and "[Title/Abstract]" not in query
        assert any(name in query for name in ("Toxoplasma gondii", "Plasmodium falciparum", "Homo sapiens", "Mus musculus", "Anopheles"))


def test_host_baselines_are_distinct_from_infection_questions():
    baseline = next(s for s in S.all_slots(O.TOXOPLASMA) if s.name == "host transcriptome · human fibroblast")
    infection = next(s for s in S.all_slots(O.TOXOPLASMA) if s.name == "host response to infection · human fibroblast")
    assert O.get(O.TOXOPLASMA).species not in A.query_for(baseline)
    assert O.get(O.TOXOPLASMA).species in A.query_for(infection)
    assert A.query_for(baseline) != A.query_for(infection)


def test_different_measurements_do_not_share_only_a_broad_assay_query():
    a = S.Slot(O.TOXOPLASMA, "palmitoylation", "PTM", "tachyzoite", "gene", (), "one")
    b = S.Slot(O.TOXOPLASMA, "phosphorylation", "PTM", "tachyzoite", "gene", (), "one")
    assert A.query_for(a) != A.query_for(b)
    assert '"palmitoylation"' in A.query_for(a) and '"phosphorylation"' in A.query_for(b)


def test_storage_fullness_never_admits_an_output_or_proposes_a_winner():
    candidate = DatasetCandidate("synthetic", "scope", "synthetic", "2020-01-01", 100, 1,
                                 "synthetic", "2026-10-07", "accepted", "Metadata diagnostic only")
    result = A.comparison_factors(candidate, "2026-10-07")
    assert result["diagnostic_score"] > 0
    assert result["score"] is None and result["admission"] == "pending"


def test_vector_queries_do_not_substitute_the_parasite_for_the_host():
    slot = next(s for s in S.all_slots(O.FALCIPARUM) if s.name == "host transcriptome · Anopheles midgut")
    query = A.query_for(slot)
    assert 'TITLE_ABS:"Anopheles"' in query
    assert O.get(O.FALCIPARUM).species not in query


def test_multiword_assays_are_quoted_and_host_tissue_synonyms_are_recorded():
    slot = next(s for s in S.all_slots(O.TOXOPLASMA) if s.name == "host proteome · mouse bone-marrow macrophage")
    query = A.query_for(slot)
    assert '"mass spectrometry"' in query
    assert '"bone marrow-derived macrophage"' in query and "BMDM" in query


def test_discovered_publications_report_rates_without_inventing_comprehensiveness():
    publication = A._publication({"source": "MED", "id": "synthetic", "citedByCount": 12,
                                  "firstPublicationDate": "2025-10-07"})
    result = A.citation_factors(publication, "2026-10-07")
    assert result["citations_per_year"] == pytest.approx(12, rel=.001)
    assert result["comprehensiveness"] is None and result["selection_score"] is None
    publication["citations"] = None
    assert A.citation_factors(publication, "2026-10-07")["citations_per_year"] is None
