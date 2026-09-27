"""Publishing one sweep must not erase measurements or provenance from another."""
import json

import pytest

from starplast import calibration as C


def test_partial_publication_preserves_other_species_and_strategies(tmp_path):
    path = str(tmp_path / "cal.json")
    old = {"Hs": {"one": {"runs": 5}, "two": {"runs": 10}},
           "Mm": {"one": {"runs": 20}}}
    C.write(old, path, {"run_dir": "first"})
    C.load(path)  # publishing invalidates this cached result
    C.write({"Hs": {"one": {"runs": 7}}}, path, {"run_dir": "second"})
    result = C.load(path)
    assert result["organisms"] == {**old, "Hs": {"one": {"runs": 7}, "two": {"runs": 10}}}
    assert result["provenance"]["Mm"]["one"] == {"run_dir": "first"}
    assert result["provenance"]["Hs"]["two"] == {"run_dir": "first"}
    assert result["provenance"]["Hs"]["one"] == {"run_dir": "second"}
    assert old["Hs"]["one"]["runs"] == 5


def test_legacy_provenance_is_retained_and_corruption_is_refused(tmp_path):
    path = tmp_path / "cal.json"
    path.write_text(json.dumps({"meta": {"run_dir": "legacy"},
                               "organisms": {"Hs": {"one": {"runs": 5}}}}))
    C.write({"Mm": {"two": {"runs": 6}}}, str(path), {"run_dir": "new"})
    assert C.load(str(path))["provenance"]["Hs"]["one"] == {"run_dir": "legacy"}
    path.write_text("{broken")
    with pytest.raises(ValueError):
        C.write({"Mm": {"two": {}}}, str(path))
    assert path.read_text() == "{broken"


def test_failed_atomic_publication_leaves_existing_file_intact(tmp_path, monkeypatch):
    path = str(tmp_path / "cal.json")
    C.write({"Hs": {"one": {"runs": 1}}}, path)
    before = (tmp_path / "cal.json").read_bytes()

    def fail(*args):
        raise OSError("disk failure")

    monkeypatch.setattr(C.os, "replace", fail)
    with pytest.raises(OSError, match="disk failure"):
        C.write({"Mm": {"one": {"runs": 2}}}, path)
    assert (tmp_path / "cal.json").read_bytes() == before
    assert list(tmp_path.iterdir()) == [tmp_path / "cal.json"]


def test_publication_pages_include_spaces_outside_the_original_pair():
    from scripts import calibrate_strategies as script
    assert script._published_spaces({"Hs": {}, "Mm": {}}) == [("Hs", "Hs"), ("Mm", "Mm")]
    text = script.calibration_doc({"Hs": {"unknown_strategy": {}}},
                                  {"date": "today", "runs": 1, "run_dir": "here"})
    assert "Hs" in text and "latest publication" in text
