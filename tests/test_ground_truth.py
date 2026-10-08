"""Truth contracts must preserve unknowns, organism identity and cohort integrity."""
from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, scorecard as C
from starplast.ground_truth import (BenchmarkEntry, StrategyBenchmark, cohort_digest,
                                   eligibility_mask, read_registry, write_registry)
from starplast.provenance import SourceFile


def _entry(tmp_path, **overrides):
    ids = ["protein_a", "protein_b", "protein_c"]
    universe = tmp_path / "universe.json"
    universe.write_text(json.dumps(ids))
    truth = tmp_path / "truth.tsv"
    truth.write_text("id\tvalue\nprotein_a\t0\nprotein_b\t1\nprotein_c\t\n")
    mask = [True, True, False]
    entry = BenchmarkEntry("test:v1", O.HUMAN, "value", C.T_VALUES, "protein", "direct_experiment",
                           "synthetic fixture testing software only", ("test_source",),
                           SourceFile.inspect(truth, "installed_cache"), SourceFile.inspect(universe, "mapping_reference"),
                           np.packbits(mask).tobytes().hex(), cohort_digest(ids, mask), 3, 2, None,
                           ("fixture",), "fixture_units", ("missing withheld",), "No certified negatives", gaps=("assay review",))
    return replace(entry, **overrides)


def test_eligibility_preserves_zero_false_and_unknowns_without_filling():
    values = pd.Series([0, False, None, np.nan, "unknown", "UnAssigned", "", "present"], index=list("abcdefgh"))
    assert eligibility_mask(values, "categorical").tolist() == [True, True, False, False, False, False, False, True]
    assert eligibility_mask(pd.Series([0, -1, np.inf, -np.inf, None, "bad"], index=list("abcdef")), "numeric").tolist() == [True, True, False, False, False, False]
    assert eligibility_mask(values, "categorical", excluded_ids=["a", "h"]).sum() == 1


def test_eligibility_refuses_ambiguous_entity_universe():
    with pytest.raises(ValueError, match="unique"):
        eligibility_mask(pd.Series([1, 2], index=["same", "same"]), "numeric")
    with pytest.raises(ValueError, match="nonmissing"):
        eligibility_mask(pd.Series([1], index=[None]), "numeric")
    with pytest.raises(ValueError, match="Declare"):
        eligibility_mask(pd.Series([1]), "guess")


def test_cohort_identity_prevents_reordering_or_replacing_entities(tmp_path):
    entry = _entry(tmp_path)
    assert entry.mask(["protein_a", "protein_b", "protein_c"]).tolist() == [True, True, False]
    for ids in (["protein_b", "protein_a", "protein_c"], ["protein_a", "protein_b", "other"]):
        with pytest.raises(ValueError, match="universe"):
            entry.mask(ids)


@pytest.mark.parametrize("change", [{"eligible_population": 3}, {"eligibility_hex": "e0"},
                                    {"eligibility_hex": "c1"}, {"eligibility_hex": "not_hex"},
                                    {"stored_population": -1}, {"cohort_sha256": "x" * 64}])
def test_corrupted_masks_counts_and_hashes_are_rejected(tmp_path, change):
    with pytest.raises(ValueError):
        _entry(tmp_path, **change)


def test_predictions_transfers_and_synthetic_truth_cannot_be_admitted_as_biology(tmp_path):
    for grade in ("prediction", "orthology_transfer", "derived_quantity", "synthetic_control", "unresolved"):
        with pytest.raises(ValueError, match="admission"):
            _entry(tmp_path, evidence_grade=grade, measured_population=3, status="admitted", gaps=())
    with pytest.raises(ValueError, match="admission"):
        _entry(tmp_path, status="admitted")


def test_registry_references_are_organism_and_task_qualified(tmp_path):
    entry = _entry(tmp_path)
    for organism, task in ((O.MOUSE, C.T_VALUES), (O.HUMAN, C.T_LABEL)):
        reference = StrategyBenchmark(organism, "fixture", task, (entry.benchmark_id,), ())
        with pytest.raises(ValueError, match="reference"):
            write_registry(tmp_path / "invalid.json", [entry], [reference])
    with pytest.raises(ValueError, match="references or explicit gaps"):
        StrategyBenchmark(O.HUMAN, "fixture", C.T_VALUES, (), ())


def test_registry_roundtrip_and_input_mutation_detection(tmp_path):
    entry = _entry(tmp_path)
    ref = StrategyBenchmark(O.HUMAN, "fixture", C.T_VALUES, (entry.benchmark_id,), ("not biological truth",))
    path = tmp_path / "registry.json"
    write_registry(path, [entry], [ref])
    assert read_registry(path) == ((entry,), (ref,))
    with pytest.raises(FileExistsError):
        write_registry(path, [entry], [ref])
    with pytest.raises(ValueError, match="Duplicate"):
        write_registry(tmp_path / "duplicate.json", [entry, entry], [ref])
    with pytest.raises(ValueError, match="Duplicate"):
        write_registry(tmp_path / "duplicate_ref.json", [entry], [ref, ref])
    from pathlib import Path
    Path(entry.truth_file.path).write_text("changed")
    with pytest.raises(ValueError, match="content"):
        read_registry(path)


def test_frozen_census_keeps_real_truth_gaps_and_covers_every_strategy():
    from pathlib import Path
    from starplast import strategies as S

    path = Path(__file__).resolve().parents[1] / "results/ground_truth_registry_2026_10_07_v2/registry.json"
    # The snapshot is a local review artifact; installed wheels may omit results.
    if not path.exists():
        pytest.skip("Review artifact is not available")
    entries, refs = read_registry(path, relocate_snapshots_to=path.parent)
    assert len(refs) == 4 * len(S.catalog())
    for organism in (O.TOXOPLASMA, O.FALCIPARUM, O.HUMAN, O.MOUSE):
        assert {ref.strategy for ref in refs if ref.organism == organism} == {s.key for s in S.catalog()}
    assert all(entry.status == "candidate" and entry.measured_population is None for entry in entries)
    assert all(ref.gaps for ref in refs)
    assert all(not ref.benchmark_ids for ref in refs if ref.strategy == "link_prediction" or ref.task == C.T_REPL)
    grades = {(entry.organism, entry.target): entry.evidence_grade for entry in entries}
    assert grades[O.FALCIPARUM, "is_exported"] == "prediction"
    assert grades[O.FALCIPARUM, "pb_transferred_phenotype"] == "orthology_transfer"
    assert grades[O.FALCIPARUM, "pb_transferred_liver_reduced"] == "orthology_transfer"
    assert grades[O.FALCIPARUM, "mean_plddt"] == "prediction"
    assert grades[O.TOXOPLASMA, "fit_invivo_PE"] == "derived_quantity"
    assert grades[O.TOXOPLASMA, "compartment_best"] == "orthology_transfer"
    assert grades[O.TOXOPLASMA, "stage_enriched_derived"] == "derived_quantity"


def test_explicit_snapshot_relocation_requires_original_content(tmp_path):
    import shutil
    from pathlib import Path

    original = tmp_path / "original"
    original.mkdir()
    entry = _entry(original)
    path = tmp_path / "registry.json"
    write_registry(path, [entry], [])
    moved = tmp_path / "moved"
    moved.mkdir()
    for source in (entry.truth_file, entry.universe_file):
        shutil.copyfile(source.path, moved / Path(source.path).name)
    relocated, _ = read_registry(path, relocate_snapshots_to=moved)
    assert relocated[0].truth_file.sha256 == entry.truth_file.sha256
    assert Path(relocated[0].truth_file.path).parent == moved
    assert relocated[0].cohort_sha256 == entry.cohort_sha256
    with pytest.raises(ValueError, match="requires content verification"):
        read_registry(path, verify_files=False, relocate_snapshots_to=moved)
    (moved / Path(entry.truth_file.path).name).write_text("wrong source")
    with pytest.raises(ValueError, match="content"):
        read_registry(path, relocate_snapshots_to=moved)
