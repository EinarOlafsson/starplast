#!/usr/bin/env python3
"""Structure resolution: the one thing the app fetches on demand.

Metadata ships and coordinates do not, so this is the only code path where "not here" is normal rather
than broken. Every test below is about that distinction: a missing structure must come back as None
with a reason, never as an exception, because AlphaFold DB skips the largest proteins and in this
proteome those are disproportionately the secreted effectors people actually click on.

No test touches the network. The download path is exercised against a stubbed opener, because asserting
that urllib works is not this project's job.
"""
from __future__ import annotations

import os
import sys
import urllib.error

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import structures as ST  # noqa: E402


# --------------------------------------------------------------------------- local models
def test_a_model_named_by_gene_id_is_found(tmp_path, monkeypatch):
    (tmp_path / "TGME49_200010.pdb").write_text("ATOM")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010") == str(tmp_path / "TGME49_200010.pdb")


def test_a_model_named_by_uniprot_is_preferred_when_one_is_known(tmp_path, monkeypatch):
    (tmp_path / "AF-Q9XYZ1-F1.pdb").write_text("ATOM")
    (tmp_path / "TGME49_200010.pdb").write_text("ATOM")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010", "Q9XYZ1").endswith("AF-Q9XYZ1-F1.pdb")


def test_the_cif_form_is_accepted_as_well_as_pdb(tmp_path, monkeypatch):
    (tmp_path / "AF-Q9XYZ1-F1-model_v4.cif").write_text("data_")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010", "Q9XYZ1").endswith(".cif")


def test_an_af3_job_directory_is_searched(tmp_path, monkeypatch):
    """AlphaFold 3 writes a directory per job, named in lower case, rather than a flat file."""
    job = tmp_path / "tgme49_200010"
    job.mkdir()
    (job / "ranked_0.cif").write_text("data_")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010").endswith("ranked_0.cif")


def test_a_job_directory_with_no_coordinates_is_not_a_hit(tmp_path, monkeypatch):
    job = tmp_path / "tgme49_200010"
    job.mkdir()
    (job / "notes.txt").write_text("nothing here")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010") is None


def test_local_directories_that_do_not_exist_are_skipped(tmp_path, monkeypatch):
    """These mirrors are machine-specific by nature; absence is not an error."""
    (tmp_path / "TGME49_200010.pdb").write_text("ATOM")
    monkeypatch.setattr(ST, "LOCAL_DIRS", ["/no/such/place", str(tmp_path)])
    assert ST.local_structure("TGME49_200010") is not None


def test_nothing_anywhere_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    assert ST.local_structure("TGME49_200010") is None


# --------------------------------------------------------------------------- AlphaFold DB
def test_a_download_is_written_to_the_cache_and_reused(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "CACHE", str(tmp_path))
    calls = {"n": 0}

    class Resp:
        def read(self):
            calls["n"] += 1
            return b"data_AF"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(ST.urllib.request, "urlopen", lambda *a, **k: Resp())
    p = ST.fetch_alphafold("Q9XYZ1")
    assert p and os.path.exists(p)
    assert ST.fetch_alphafold("Q9XYZ1") == p
    assert calls["n"] == 1, "a cached model must not be downloaded again"


def test_being_offline_is_reported_as_absence_not_as_a_crash(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "CACHE", str(tmp_path))

    def boom(*a, **k):
        raise urllib.error.URLError("no route to host")

    monkeypatch.setattr(ST.urllib.request, "urlopen", boom)
    assert ST.fetch_alphafold("Q9XYZ1") is None


def test_an_empty_response_is_not_written_as_a_structure(tmp_path, monkeypatch):
    """A zero-byte file in the cache would then be returned forever as though it were a model."""
    monkeypatch.setattr(ST, "CACHE", str(tmp_path))

    class Empty:
        def read(self):
            return b""

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(ST.urllib.request, "urlopen", lambda *a, **k: Empty())
    assert ST.fetch_alphafold("Q9XYZ1") is None
    assert not os.listdir(tmp_path)


def test_no_uniprot_means_no_download_attempt(monkeypatch):
    monkeypatch.setattr(ST.urllib.request, "urlopen",
                        lambda *a, **k: pytest.fail("must not call the network without an accession"))
    assert ST.fetch_alphafold("") is None
    assert ST.fetch_alphafold(None) is None


def test_a_zero_byte_cached_file_is_refetched(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "CACHE", str(tmp_path))
    open(os.path.join(tmp_path, "AF-Q9XYZ1-F1-model_v4.cif"), "wb").close()

    class Resp:
        def read(self):
            return b"data_AF"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(ST.urllib.request, "urlopen", lambda *a, **k: Resp())
    p = ST.fetch_alphafold("Q9XYZ1")
    assert p and os.path.getsize(p) > 0


# --------------------------------------------------------------------------- the resolver
def test_a_local_model_wins_and_says_so(tmp_path, monkeypatch):
    (tmp_path / "TGME49_200010.pdb").write_text("ATOM")
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    monkeypatch.setattr(ST.urllib.request, "urlopen",
                        lambda *a, **k: pytest.fail("must not fetch what is already here"))
    p, origin = ST.find_structure("TGME49_200010", "Q9XYZ1")
    assert p and origin == "local"


def test_with_the_network_disabled_the_reason_says_so(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    p, why = ST.find_structure("TGME49_200010", "Q9XYZ1", allow_network=False)
    assert p is None and "network disabled" in why


def test_a_gene_with_no_uniprot_accession_says_that_rather_than_failing(tmp_path, monkeypatch):
    """Distinct from "AlphaFold has no model": one is a gap in our mapping, the other in theirs."""
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    p, why = ST.find_structure("TGME49_200010", None)
    assert p is None and "UniProt" in why


def test_no_alphafold_model_is_a_distinct_reason(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    monkeypatch.setattr(ST, "fetch_alphafold", lambda *a, **k: None)
    p, why = ST.find_structure("TGME49_200010", "Q9XYZ1")
    assert p is None and "no AlphaFold model" in why


def test_a_fetched_model_is_labelled_with_its_origin(tmp_path, monkeypatch):
    monkeypatch.setattr(ST, "LOCAL_DIRS", [str(tmp_path)])
    monkeypatch.setattr(ST, "fetch_alphafold", lambda *a, **k: "/tmp/AF.cif")
    p, origin = ST.find_structure("TGME49_200010", "Q9XYZ1")
    assert p == "/tmp/AF.cif" and origin == "AlphaFold DB"


# --------------------------------------------------------------------------- crosslink complexes
class _Row:
    def __init__(self, model_dir=None, model_files=None):
        self.model_dir = model_dir
        self.model_files = model_files


def test_only_complexes_that_exist_locally_are_returned(tmp_path):
    d = tmp_path / "models"
    d.mkdir()
    (d / "a.cif").write_text("data_")
    out = ST.crosslink_model_paths(_Row("models", "a.cif;b.cif"), str(tmp_path))
    assert out == [str(d / "a.cif")]


def test_a_row_with_no_models_returns_an_empty_list(tmp_path):
    assert ST.crosslink_model_paths(_Row(None, None), str(tmp_path)) == []
    assert ST.crosslink_model_paths(_Row("models", None), str(tmp_path)) == []


def test_a_row_missing_the_attributes_entirely_does_not_raise(tmp_path):
    assert ST.crosslink_model_paths(object(), str(tmp_path)) == []
