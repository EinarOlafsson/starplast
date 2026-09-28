"""Strategy cards: the same four bars for every strategy, its test explained, and real examples.

The card is what a person reads first, so it is held to three rules here: the headline bars are the
same four, in the same places, for every strategy of a task (and computed from the shipped
calibration, never invented); every strategy says what it does, how it is evaluated, how it fails
and how it succeeds; and the worked examples are real runs from the calibration sweep, re-checked
against the record they came from.
"""
from __future__ import annotations

import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import calibration as CAL  # noqa: E402
from starplast import scorecard as SC  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import strategy_explainers as EX  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "results", "calibration_2026-09-26b", "runs.jsonl")


# --------------------------------------------------------------------------- the headline mapping
PLAIN = {
    "label calls": ("Right calls", "Fair across classes", "accuracy", "macro_f1"),
    "ranking": ("True ones ranked first", "Clean top of the list", "auroc", "auprc_lift"),
    "values": ("Order predicted", "Variance explained", "spearman", "r2"),
    "cluster recovery": ("Label falls out as a cluster", "Partition agreement",
                         "weighted_f1_clusters", "ari"),
    "set retrieval": ("Returned genes that are real", "Members found", "precision", "recall"),
    "replication": ("Findings that hold", "Beyond chance", "replication_rate", "replication_lift"),
}


def test_every_task_has_the_same_four_positions():
    assert set(SC.HEADLINE) == set(SC.TASKS) == set(SC.REACH)
    for task, (l3, l4, k3, k4) in PLAIN.items():
        bars = SC.headline(task)
        assert len(bars) == 4
        assert bars[0].label == "Better than chance" and bars[0].key == "skill"
        assert bars[0].chance == 0.0
        assert bars[1].label == "Reach"
        assert (bars[2].label, bars[3].label) == (l3, l4)
        assert (bars[2].key, bars[3].key) == (k3, k4)
        for h in bars[1:]:
            assert h.key in SC.TASKS[task].metrics, (task, h.key)
            assert h.technical and h.reading and h.scale in SC.SCALES
        if task not in ("label calls", "values"):
            # A task without coverage says which analogue stands in for it.
            assert bars[1].note, task


def test_bar_positions_and_values_read_plainly():
    assert SC.position("unit", 0.5) == 0.5
    assert SC.position("unit", -0.3) == 0.0 and SC.position("unit", 2) == 1.0
    assert SC.position("lift", 1.0) == pytest.approx(1 / 6)
    assert SC.position("lift", 32) == 1.0 and SC.position("lift", 0.1) == 0.0
    assert SC.position("count", 1) == 0.0 and SC.position("count", 1000) == 1.0
    assert math.isnan(SC.position("unit", None))
    assert SC.fmt_value("unit", 0.1234) == "0.12"
    assert SC.fmt_value("lift", 7.73) == "x7.7" and SC.fmt_value("lift", 15.2) == "x15"
    assert SC.fmt_value("count", 1234.0) == "1,234"
    assert SC.fmt_value("unit", float("nan")) == "--"


def test_chance_levels_are_fixed_measured_or_derived():
    # The verdict's own measured chance when a bar shows the very metric the verdict rests on.
    bars = SC.headline_bars("label calls", {"accuracy": 0.6, "coverage": 0.9, "macro_f1": 0.4},
                            {"skill": 0.3, "observed": 0.6, "chance": 0.45})
    assert [b["label"] for b in bars] == ["Better than chance", "Reach", "Right calls",
                                          "Fair across classes"]
    assert bars[2]["chance"] == 0.45 and math.isnan(bars[3]["chance"])
    assert math.isnan(bars[1]["chance"])                  # coverage has no chance level
    # Replication's chance is the scrambled-evidence rate on the same card.
    rep = SC.headline_bars("replication", {"replication_rate": 0.9, "findings": 60,
                                           "null_rate": 0.05, "replication_lift": 18},
                           {"skill": 0.9})
    assert rep[2]["chance"] == 0.05 and rep[3]["chance"] == 1.0 and rep[1]["value"] == 60
    # A set's chance levels follow from fold enrichment: precision / fold is the members' share.
    st = SC.headline_bars("set retrieval", {"precision": 0.2, "recall": 0.5,
                                            "fold_enrichment": 10.0, "returned": 100}, {})
    assert st[2]["chance"] == pytest.approx(0.02) and st[3]["chance"] == pytest.approx(0.05)
    # Clusters show the clustered share, the unclustered share turned round, interval included.
    cl = SC.headline_bars("cluster recovery", {"noise_share": {"mean": 0.3, "low": 0.2,
                                                               "high": 0.4}}, {})
    assert cl[1]["value"] == pytest.approx(0.7)
    assert (cl[1]["low"], cl[1]["high"]) == (pytest.approx(0.6), pytest.approx(0.8))


def test_every_calibrated_strategy_fills_its_bars_from_the_shipped_calibration():
    for org in ("Tg", "Pf"):
        for s in S.catalog():
            e = CAL.entry(s.key, org) or {}
            d = e.get("default") or {}
            bars = SC.headline_bars(s.task, d.get("scorecard") or {}, d)
            assert len(bars) == 4
            if d.get("skill") is not None and math.isfinite(d["skill"]):
                assert bars[0]["value"] == pytest.approx(d["skill"])
                assert bars[0]["low"] == pytest.approx(d["skill_low"])
            card = d.get("scorecard") or {}
            for b in bars[2:]:
                if b["key"] in card:
                    assert b["value"] == pytest.approx(card[b["key"]]["mean"])
                    assert b["low"] == pytest.approx(card[b["key"]]["low"])
            for b in bars:
                assert b["reading"].startswith(b["label"])
                p = b["position"]
                assert math.isnan(p) or 0.0 <= p <= 1.0


# --------------------------------------------------------------------------- about this test
def test_every_strategy_explains_its_test_in_four_fields():
    keys = {s.key for s in S.catalog()}
    assert set(EX.EXPLAINERS) == keys
    assert [f for f, _t in EX.FIELDS] == ["does", "evaluated", "failure", "success"]
    seen = {}
    for key in keys:
        e = EX.explainer(key)
        for field, _title in EX.FIELDS:
            text = e[field].strip()
            assert text, f"{key}: empty {field}"
            assert 80 <= len(text) <= 700, (key, field, len(text))
            assert EX.first_sentence(text) and len(EX.first_sentence(text)) <= len(text)
            assert text not in seen, f"{key} {field} repeats {seen.get(text)}"
            seen[text] = key
        assert EX.markdown(key).count("**") == 8


def test_first_sentence_stops_at_the_first_full_stop():
    assert EX.first_sentence("It does one thing well. Then more.") == "It does one thing well."
    assert EX.first_sentence("No stop at all") == "No stop at all"


# --------------------------------------------------------------------------- worked examples
@pytest.fixture(scope="module")
def shipped_examples():
    with open(EX.EXAMPLES) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def runs():
    if not os.path.exists(RUNS):
        pytest.skip("the calibration runs are not in this checkout")
    with open(RUNS) as fh:
        return {r["id"]: r for r in (json.loads(line) for line in fh if line.strip())}


def test_every_strategy_ships_a_failure_and_a_success_per_organism(shipped_examples):
    for org in ("Tg", "Pf"):
        assert set(shipped_examples[org]) == {s.key for s in S.catalog()}, org
        for key, e in shipped_examples[org].items():
            f, ok = e["failure"], e["success"]
            assert e["task"] == S.get(key).task
            assert f["source"] in ("calibration", "noise table") and f["explanation"]
            if f["source"] == "calibration":
                assert f["verdict"] == "FAIL", (org, key)
            else:
                assert f["why_noise"].startswith("none of its"), (org, key)
            assert ok["explanation"]
            if ok["source"] == "calibration":
                assert ok["verdict"] == "PASS", (org, key)
                new = ok.get("new") or {}
                assert new.get("rows") or new.get("skipped") or new.get("what"), (org, key)
                for r in new.get("rows") or []:
                    assert {"gene_id", "product", "call", "support"} <= set(r)
                assert len(new.get("rows") or []) <= 5
            else:
                assert ok["source"] == "none" and e["passes"] == 0, (org, key)


def test_examples_are_the_calibration_runs_they_claim_to_be(shipped_examples, runs):
    """Every real example is re-read from runs.jsonl: same run, same verdict, same numbers."""
    checked = 0
    for org in ("Tg", "Pf"):
        for key, e in shipped_examples[org].items():
            for ex in (e["failure"], e["success"]):
                if ex.get("source") != "calibration":
                    continue
                r = runs[ex["run_id"]]
                assert (r["organism"], r["strategy"], r["verdict"]) == (org, key, ex["verdict"])
                assert r["seed"] == ex["seed"] and r["settings"] == ex["settings"]
                assert ex["observed"] == pytest.approx(r["observed"], abs=1e-3)
                assert ex["chance"] == pytest.approx(r["null_mean"], abs=1e-3)
                for m, v in ex["card"].items():
                    if v is not None and r["scorecard"].get(m) is not None:
                        assert v == pytest.approx(r["scorecard"][m], abs=1e-3)
                checked += 1
    assert checked > 100


def test_a_strategy_that_never_failed_on_real_data_fails_on_noise(shipped_examples, runs):
    """The noise-table fallback is used exactly when no real run of that strategy failed."""
    def scored(r):          # a skill can be computed: the choice rule only takes those
        o, c = r.get("observed"), r.get("null_mean")
        return o is not None and c is not None and abs(1 - c) > 1e-9

    failed = {(r["organism"], r["strategy"]) for r in runs.values()
              if r.get("verdict") == "FAIL" and scored(r)}
    for org in ("Tg", "Pf"):
        for key, e in shipped_examples[org].items():
            real = e["failure"]["source"] == "calibration"
            assert real == ((org, key) in failed), (org, key)


def test_the_docs_page_has_every_card():
    path = os.path.join(ROOT, "docs", "strategy_cards.md")
    text = open(path, encoding="utf8").read()
    for s in S.catalog():
        assert f"## {s.number:02d} · {s.title}" in text
    for task in SC.TASKS:
        for h in SC.headline(task):
            assert h.label in text
    assert "What failure looks like" in text and "Works when" in text
