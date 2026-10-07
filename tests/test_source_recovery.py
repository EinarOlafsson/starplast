"""Downloads and relocated inputs require explicit source identity and valid content."""
import hashlib
import io
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import recover_source_files as R
from starplast import datasets as D, organisms as O


def test_relocated_file_needs_url_species_name_and_matching_receipt_hash(tmp_path):
    folder = tmp_path / "relocated"
    folder.mkdir()
    file = folder / "processed.tsv"
    file.write_text("id\tvalue\na\t1\n")
    source = D.Dataset("synthetic", "synthetic", "reference", "RNAseq", "test", organism=O.HUMAN,
                       path="datasets/original/processed.tsv", url="https://example.org/processed.tsv")
    receipt = {"space": O.HUMAN, "urls": [source.url], "files": [{"name": file.name,
               "path": "relocated/processed.tsv", "sha256": hashlib.sha256(file.read_bytes()).hexdigest()}]}
    binding, status = R.receipt_binding(source, tmp_path, [(tmp_path / "receipt.json", receipt)])
    assert status == "receipt_verified" and binding.path == str(file)
    receipt["space"] = O.MOUSE
    assert R.receipt_binding(source, tmp_path, [(tmp_path / "receipt.json", receipt)]) == (None, "not_located")
    receipt["space"] = O.HUMAN
    receipt["files"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="hash mismatch"):
        R.receipt_binding(source, tmp_path, [(tmp_path / "receipt.json", receipt)])


@pytest.mark.parametrize("payload", [b"<html>error</html>", b"<!doctype html>error", b"<Error>denied</Error>", b""])
def test_successful_http_error_pages_never_land_as_datasets(tmp_path, monkeypatch, payload):
    monkeypatch.setattr(R.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(payload))
    target = tmp_path / "dataset.tsv"
    with pytest.raises(ValueError):
        R.download("https://example.org/dataset.tsv", target)
    assert not target.exists() and not list(tmp_path.glob(".partial-*"))


def test_downloads_are_bounded_and_do_not_overwrite_an_input(tmp_path, monkeypatch):
    monkeypatch.setattr(R.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(b"a\t1\n"))
    target = tmp_path / "dataset.tsv"
    with pytest.raises(ValueError, match="scope"):
        R.download("https://example.org/dataset.tsv", target, limit=1)
    assert not target.exists()
    R.download("https://example.org/dataset.tsv", target)
    assert target.read_text() == "a\t1\n"
    with pytest.raises(ValueError, match="overwrite"):
        R.download("https://example.org/dataset.tsv", target)


def test_spreadsheet_downloads_cannot_be_plain_text_renamed_as_xlsx(tmp_path, monkeypatch):
    monkeypatch.setattr(R.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(b"this is not Excel"))
    target = tmp_path / "dataset.xlsx"
    with pytest.raises(Exception):
        R.download("https://example.org/dataset.xlsx", target)
    assert not target.exists()
