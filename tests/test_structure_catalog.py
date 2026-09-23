"""AF3 provenance, residue mapping and missing-data guarantees."""
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from starplast import structure_catalog as C


def test_exact_sequence_collision_abstains_even_with_alias():
    sequences = {"g1": "ACDE", "g2": "ACDE"}
    assert C.match_sequence("ACDE", sequences, {"ACDE": ["g1", "g2"]}, "g1")[0] is None


def test_fragment_requires_unique_exact_subsequence_and_explicit_alias():
    sequences = {"g1": "MACDEF"}
    assert C.match_sequence("ACD", sequences, {}, "g1") == ("g1", 1, 4, "exact_fragment_sequence")
    assert C.match_sequence("ACD", sequences, {}, None)[0] is None
    assert C.match_sequence("ACD", {"g1": "ACDACD"}, {}, "g1")[0] is None
    assert C.match_sequence("ACE", sequences, {}, "g1")[0] is None


def test_seed_replicates_are_not_independent_models(tmp_path):
    job = tmp_path / "af3_Q1"
    job.mkdir()
    root = job / "af3_Q1_model.cif"
    root.touch()
    seed = job / "seed-1_sample-1"
    seed.mkdir()
    (seed / "model.cif").touch()
    assert list(C.model_files([tmp_path, tmp_path])) == [root]


def test_overlap_uses_residue_confidence_once_and_fragments_have_no_global_shape(tmp_path, monkeypatch):
    for name in ["af3_Q1_f1", "af3_Q1_f2", "bait_g1__host_X"]:
        directory = tmp_path / name
        directory.mkdir()
        (directory / (name + "_model.cif")).write_text(name)
    def chain(path):
        if "f1" in str(path):
            return "ACDE", np.array([10, 20, 30, 40]), np.ones((4, 3))
        return "DEFG", np.array([80, 90, 100, 60]), np.ones((4, 3))
    monkeypatch.setattr(C, "read_chain", chain)
    models, features = C.inventory([tmp_path], {"g1": "ACDEFG"}, {"Q1": "g1"}, log=lambda *a, **k: None)
    row = features.iloc[0]
    assert row.af3_sequence_coverage == 1
    assert row.af3_mean_plddt == pytest.approx(np.mean([10, 20, 80, 90, 100, 60]))
    assert row.af3_confident_sequence_fraction == .5
    assert "af3_confident_rg_angstrom" not in row
    assert models.model_type.value_counts()["fragment"] == 2
    assert models.model_type.value_counts()["other_job"] == 1


def test_missing_residues_are_not_low_confidence_measurements(tmp_path, monkeypatch):
    p = tmp_path / "AF3-Q1.cif"
    p.write_text("fixture")
    monkeypatch.setattr(C, "read_chain", lambda p: ("ACD", np.full(3, 80), np.ones((3, 3))))
    _, features = C.inventory([tmp_path], {"g1": "MACDEF"}, {"Q1": "g1"}, log=lambda *a, **k: None)
    assert features.iloc[0].af3_sequence_coverage == .5
    assert features.iloc[0].af3_mean_plddt == 80
    assert features.iloc[0].af3_low_confidence_modelled_fraction == 0


def test_feature_join_preserves_gene_order_and_missingness(tmp_path):
    p = tmp_path / "features.parquet"
    pd.DataFrame({"gene_id": ["g2", "g1"], "af3_mean_plddt": [80, 90]}).to_parquet(p)
    nodes = pd.DataFrame({"gene_id": ["g1", "g3", "g2"]}, index=[7, 3, 2])
    result = C.attach_features(nodes, p)
    assert result.index.tolist() == [7, 3, 2]
    assert result.af3_mean_plddt.iloc[0] == 90
    assert pd.isna(result.af3_mean_plddt.iloc[1])
    assert "af3_mean_plddt" not in nodes


def test_indexed_viewer_prefers_complete_existing_file(tmp_path):
    paths = [tmp_path / name for name in ["partial.cif", "full.cif", "missing.cif"]]
    paths[0].touch(); paths[1].touch()
    data = pd.DataFrame({"gene_id": ["g1"] * 3, "model_type": ["fragment", "full_protein", "full_protein"],
                         "sequence_coverage": [.5, 1, 1], "mean_plddt": [90, 60, 99],
                         "path": [str(p) for p in paths]})
    index = tmp_path / "index.parquet"
    data.to_parquet(index)
    assert C.local_model("g1", index) == str(paths[1])
    assert C.local_model("g2", index) is None


def test_confident_geometry_excludes_uncertain_residues():
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [100, 0, 0]])
    values = C.geometry(coords, np.array([90, 90, 90, 10]))
    assert values['confident_rg_angstrom'] == pytest.approx(np.sqrt(2 / 3))
    assert values['confident_contacts_per_residue'] == 0


def test_unreadable_model_is_auditable_not_silently_dropped(tmp_path, monkeypatch):
    (tmp_path / "AF3-Q1.cif").write_text("bad")
    def fail(p):
        raise ValueError("invalid coordinates")
    monkeypatch.setattr(C, "read_chain", fail)
    models, features = C.inventory([tmp_path], {}, log=lambda *a, **k: None)
    assert features.empty
    assert models.iloc[0].mapping_status == "parse_error"
    assert "invalid coordinates" in models.iloc[0].error
