"""Publisher redirects cannot justify guessing or silently relabelling a legacy input."""
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from retrieve_primary_table import retrieve
from starplast import datasets as D


@pytest.mark.parametrize("change", [{"source_id": "another_source"}, {"recorded_registry_url": "https://invalid.example/table"},
                                    {"status": 404}, {"final_url": "https://invalid.example/unnamed"},
                                    {"final_url": "https://invalid.example/raw.tar"}])
def test_mismatched_or_unnamed_publisher_resources_never_download(tmp_path, monkeypatch, change):
    source = D.get("host_macrophage_surfaceome")
    metadata = {"source_id": source.key, "recorded_registry_url": source.url, "status": 200,
                "final_url": "https://publisher.example/table.xlsx", **change}
    path = tmp_path / "headers.json"
    path.write_text(json.dumps(metadata))
    def unexpected(*args, **kwargs):
        raise AssertionError("Invalid source metadata triggered a download")
    monkeypatch.setattr("retrieve_primary_table.download", unexpected)
    with pytest.raises(ValueError):
        retrieve(source.key, tmp_path, path, tmp_path / "output")
    assert not (tmp_path / "output").exists()
