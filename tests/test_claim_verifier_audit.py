"""The verifier audit cannot swap evidence, mix genes or promote outside-range coverage."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

import audit_claim_verifiers as A  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402


@pytest.fixture
def ctx():
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(63)],
                          "label": ["a"] * 30 + ["b"] * 30 + [None] * 3})
    # The first unknown has a real residual-supported link; the other two have only
    # popular, below-expected links. Raw counts would make every unknown seem covered.
    graph = {"comention__a": np.array([*range(29), *range(30, 59), 0, 1, 2]),
             "comention__b": np.array([*range(1, 30), *range(31, 60), 60, 61, 62]),
             "comention__w": np.full(61, 100.0),
             "comention__r": np.array([1.0] * 59 + [-1.0, 0.0])}
    return S.Context(nodes, graph=graph, organism="example")


def test_corrected_graph_uses_positive_residuals_without_mutating_the_source(ctx):
    before = ctx.operator("comention").copy()
    candidate, reason = A.literature_context(ctx, "label", "comention")
    assert reason == "measured"
    assert candidate.graph["comention__w"][-3:].tolist() == [1.0, 0.0, 0.0]
    assert (ctx.graph["comention__w"] == 100.0).all()
    assert (ctx.operator("comention") != before).nnz == 0
    pred, _ = S.propagate(candidate.operator("comention"), ctx.truth("label"), [60, 61, 62])
    assert pred.iloc[60] == "a" and pred.iloc[[61, 62]].isna().all()
    raw, _ = A.literature_context(ctx, "label", "comention", corrected=False)
    raw_pred, _ = S.propagate(raw.operator("comention"), ctx.truth("label"), [60, 61, 62])
    assert raw_pred.iloc[[60, 61, 62]].notna().all()


def test_banned_layer_is_refused_even_when_another_graph_is_available(ctx, monkeypatch):
    monkeypatch.setattr(ctx, "banned_layers", lambda _target: {"comention"})
    candidate, reason = A.literature_context(ctx, "label", "comention")
    assert candidate is None and reason == "banned layer"


@pytest.mark.parametrize("case, reason", [
    ("missing_layer", "layer unavailable"),
    ("missing_residual", "attention residuals unavailable"),
    ("nonfinite", "nonfinite edge weights"),
    ("nonpositive", "no positive edges"),
])
def test_unusable_evidence_is_a_refusal_not_a_raw_count_fallback(ctx, case, reason):
    if case == "missing_layer":
        ctx.graph.pop("comention__a")
    elif case == "missing_residual":
        ctx.graph.pop("comention__r")
    elif case == "nonfinite":
        ctx.graph["comention__r"][0] = np.nan
    else:
        ctx.graph["comention__r"][:] = -1
    candidate, why = A.literature_context(ctx, "label", "comention")
    assert candidate is None and why == reason


def test_only_literature_graphs_can_be_audited_here(ctx):
    with pytest.raises(ValueError, match="not a literature layer"):
        A.literature_context(ctx, "label", "compartment")


@pytest.fixture
def reference(ctx):
    return T.evaluate(ctx, "layer_propagation", "label", settings={"layer": "comention"})


@pytest.mark.parametrize("column,value", [("gene_id", "foreign"), ("fold", 99),
                                         ("truth", "foreign"), ("seed", 99),
                                         ("organism", "foreign"), ("target", "foreign")])
def test_stale_or_misaligned_records_cannot_be_compared(reference, column, value):
    measured = reference.copy()
    measured.loc[0, column] = value
    with pytest.raises(ValueError, match="differ in identities"):
        A.validate_reference(reference, measured, "layer_propagation")


def test_duplicate_missing_and_reordered_rows_are_handled_explicitly(reference):
    A.validate_reference(reference, reference.iloc[::-1], "layer_propagation")
    with pytest.raises(ValueError, match="one held-out row"):
        A.validate_reference(reference, pd.concat([reference, reference.iloc[:1]]), "layer_propagation")
    with pytest.raises(ValueError, match="differ in identities"):
        A.validate_reference(reference, reference.iloc[1:], "layer_propagation")


def test_new_coverage_excludes_existing_checks_and_outside_range():
    claims = pd.DataFrame({"gene_id": ["a", "b", "c", "d"],
                          "status": ["tested", "untested", "outside tested range", "untested"]})
    calls = pd.DataFrame({"gene_id": ["a", "b", "c", "d", "e"],
                         "prediction": ["x", "y", "z", None, "x"]})
    coverage = A.coverage(claims, calls)
    assert coverage["unknown_reached"] == 4 and coverage["claims_reached"] == 3
    assert coverage["new_in_range"] == 1 and coverage["outside_range_reached"] == 1
    assert coverage["tested_with_candidate"] == 2


def test_audit_preserves_folds_and_reports_real_unknown_gene_reach(ctx, reference):
    recipe = pd.Series({"generator": "layer_propagation", "verifiers": ""})
    claims = pd.DataFrame({"gene_id": ["g60", "g61", "g62"], "status": ["untested"] * 3})
    row, measured, calls = A.audit_candidate(ctx, "label", "comention", True,
                                            reference, recipe, claims)
    A.validate_reference(reference, measured, "layer_propagation")
    assert row["status"] == "measured" and len(measured) == 60
    assert row["unknown_reached"] == 1 and row["new_in_range"] == 1
    assert calls.loc[calls.prediction.notna(), "gene_id"].tolist() == ["g60"]
    assert not row["eligible"], "a flawless generator supplies no mistakes to test independence"


def test_raw_control_is_never_eligible_even_if_its_errors_look_independent(ctx, reference, monkeypatch):
    monkeypatch.setattr(A.C, "shared_mistakes", lambda *_: {"ratio_high": 1.0})
    recipe = pd.Series({"generator": "layer_propagation", "verifiers": ""})
    claims = pd.DataFrame({"gene_id": [], "status": []})
    row, _, _ = A.audit_candidate(ctx, "label", "comention", False, reference, recipe, claims)
    assert row["independent"] and not row["eligible"]
    assert "combined_calibration_error" not in row


@pytest.mark.parametrize("new_reach, confident, precision, calibrated, eligible", [
    (True, 100, 0.9, True, True), (False, 100, 0.9, True, False),
    (True, 20, 0.9, True, False), (True, 100, 0.7, True, False),
    (True, 100, 0.9, False, False),
])
def test_eligibility_needs_calibration_confident_support_and_new_reach(
        ctx, reference, monkeypatch, new_reach, confident, precision, calibrated, eligible):
    from types import SimpleNamespace
    monkeypatch.setattr(A.C, "shared_mistakes", lambda *_: {"ratio_high": 1.0})
    model = SimpleNamespace(calibration_error=0.02 if calibrated else 0.08,
                            calibrated=calibrated, certainty=SimpleNamespace(calibrated=True),
                            _oof=np.full(confident, 0.9),
                            _y=np.array([1] * int(confident * precision) +
                                        [0] * (confident - int(confident * precision))))
    monkeypatch.setattr(A.C, "verified_model", lambda *_: model)
    recipe = pd.Series({"generator": "layer_propagation", "verifiers": ""})
    claims = pd.DataFrame({"gene_id": ["g60"], "status": ["untested" if new_reach else "tested"]})
    row, _, _ = A.audit_candidate(ctx, "label", "comention", True, reference, recipe, claims)
    assert row["eligible"] == eligible


def test_cli_finishes_and_saves_a_reexecutable_notebook(tmp_path, monkeypatch, capsys):
    import json
    from notebook_runner import ExecutedNotebook
    summary = pd.DataFrame([{
        "organism": "example", "target": "label", "verifier": "literature", "status": "measured",
        "corrected": True, "independent": False, "eligible": False, "ratio": 2.0,
        "ratio_low": 1.8, "ratio_high": 2.2, "agrees_n": 40, "agrees_rate": 0.8,
        "disagrees_rate": 0.3, "unknown_genes": 10, "unknown_reached": 5,
        "new_in_range": 2, "outside_range_reached": 1,
    }])
    monkeypatch.setattr(A, "audit", lambda: (summary, pd.DataFrame({"gene": [0]}),
                                            pd.DataFrame({"gene_id": ["g60"]})))
    monkeypatch.setattr(A, "fingerprints", lambda: {"input": "hash"})
    out = tmp_path / "audit"
    assert A.main(["--out", str(out)]) == 0
    assert "Wrote" in capsys.readouterr().out
    notebook = json.loads((out / "audit.ipynb").read_text())
    assert (out / "manifest.json").exists() and (out / "summary.csv").exists()
    assert pd.read_parquet(out / "unknown_calls.parquet").gene_id.tolist() == ["g60"]
    # Notebook inputs must be defined in its own cells, not only in the CLI's injected namespace.
    replay = ExecutedNotebook("Replay")
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            replay.code("".join(cell["source"]))
    assert replay.ns["summary"].equals(summary)
    with pytest.raises(FileExistsError):
        A.main(["--out", str(out)])
