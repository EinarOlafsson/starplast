"""Calibration: the summary must not flatter the sweep it summarises.

Planted runs throughout, because the point of these tests is the arithmetic and the discipline --
that choosing a setting and reporting it use different seeds, that an interval widens when the runs
within a target are correlated, and that a grade follows the numbers rather than the other way
round. The shipped file is then checked for the properties a reader relies on.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import calibration as C  # noqa: E402


def _runs(strategy="planted", organism="Tg", settings=None, targets=("a", "b"), seeds=(1, 2, 3, 4, 5),
          observed=0.6, chance=0.2, verdict="PASS"):
    rows = []
    for target in targets:
        for seed in seeds:
            rows.append({"organism": organism, "strategy": strategy,
                         "settings": {"target": target, **(settings or {})}, "seed": seed,
                         "verdict": verdict, "metric": "correct calls", "observed": observed,
                         "null_mean": chance, "null_sd": 0.05, "null_high": chance + 0.1,
                         "p_value": 0.01, "effect": observed - chance, "min_effect": 0.05,
                         "n_hidden": 100, "hidden": "25%", "null_kind": "shuffled",
                         "note": "", "wall_seconds": 1.0})
    return rows


# --------------------------------------------------------------------------- the arithmetic
def test_skill_is_zero_at_chance_and_one_at_perfect():
    df = C.runs_frame(_runs(observed=0.2, chance=0.2))
    assert abs(df["skill"].mean()) < 1e-9
    df = C.runs_frame(_runs(observed=1.0, chance=0.2))
    assert abs(df["skill"].mean() - 1.0) < 1e-9
    # Worse than the null is negative, not clipped: a strategy can be worse than shuffling.
    df = C.runs_frame(_runs(observed=0.1, chance=0.5))
    assert df["skill"].mean() < 0


def test_the_target_is_split_off_from_the_setting():
    """Which label was held out is not a setting of the strategy: pooling them would report one
    number for a strategy that works on localization and fails on the cell cycle."""
    df = C.runs_frame(_runs(settings={"k": 15}))
    assert set(df["target"]) == {"a", "b"}
    assert set(df["setting"]) == {json.dumps({"k": 15})}


def test_a_run_that_could_not_run_is_counted_but_not_scored():
    rows = _runs() + [{"organism": "Tg", "strategy": "planted", "settings": {"target": "a"},
                       "seed": 1, "verdict": "NOT RUN", "metric": None, "observed": None,
                       "null_mean": None, "note": "the layer is built from the target"}]
    df = C.runs_frame(rows)
    summary = C.summarise(df)["Tg"]["planted"]
    assert summary["runs"] == len(rows)
    assert summary["default"]["conclusive"] == len(rows) - 1
    assert summary["default"]["not_run"] == 1


# --------------------------------------------------------------------------- the intervals
def test_the_interval_widens_when_runs_within_a_target_agree_with_each_other():
    """Runs on one held-out label share its labels and are not independent observations.

    A bootstrap over runs alone would treat twenty correlated runs as twenty, and report an
    interval several times too narrow. Resampling targets first is what stops that, and this is the
    test that it does: the same spread, arranged as two targets that differ, must give a wider
    interval than as twenty runs that differ individually.
    """
    rng = np.random.default_rng(0)
    values = np.concatenate([np.full(10, 0.2), np.full(10, 0.8)])         # all the variance is
    clustered = ["a"] * 10 + ["b"] * 10                                   # between the targets
    scattered = [f"t{i}" for i in range(20)]
    _m, low_c, high_c = C.cluster_bootstrap(values, clustered, n_boot=1500, seed=1)
    _m, low_s, high_s = C.cluster_bootstrap(rng.permutation(values), scattered, n_boot=1500, seed=1)
    assert (high_c - low_c) > 2 * (high_s - low_s)


def test_one_run_gives_a_point_and_not_an_interval():
    mean, low, high = C.cluster_bootstrap([0.5], ["a"])
    assert mean == low == high == 0.5


def test_a_wilson_interval_is_inside_the_unit_range_at_the_extremes():
    for k, n in ((0, 5), (5, 5), (1, 200)):
        low, high = C.wilson(k, n)
        assert 0.0 <= low <= high <= 1.0
    assert all(math.isnan(x) for x in C.wilson(0, 0))


# --------------------------------------------------------------------------- choosing honestly
def test_the_tuned_setting_is_chosen_on_some_seeds_and_reported_on_others():
    """A setting picked as the best of many is optimistic on the runs that picked it.

    Here one setting is lucky only on the choosing seeds and ordinary on the rest. The tuned figure
    must come out ordinary, because it is reported on seeds that had no part in the choice.
    """
    rows = []
    for setting, on_choice, on_report in (({"k": 5}, 0.9, 0.3), ({"k": 15}, 0.5, 0.5)):
        for seed in C.CHOOSE_SEEDS:
            rows += _runs(settings=setting, seeds=(seed,), observed=on_choice, chance=0.0)
        for seed in C.REPORT_SEEDS:
            rows += _runs(settings=setting, seeds=(seed,), observed=on_report, chance=0.0)
    entry = C.summarise(C.runs_frame(rows))["Tg"]["planted"]
    assert entry["tuned"]["setting"] == {"k": 5}                      # the lucky one was chosen
    assert abs(entry["tuned"]["skill"] - 0.3) < 1e-6                  # and reported honestly
    assert entry["tuned"]["chosen_on_seeds"] == list(C.CHOOSE_SEEDS)


def test_the_default_row_is_the_default_setting_when_one_is_named():
    rows = _runs(settings={"k": 5}, observed=0.9) + _runs(settings={"k": 15}, observed=0.3)
    df = C.runs_frame(rows)
    defaults = {("Tg", "planted"): json.dumps({"k": 15})}
    entry = C.summarise(df, defaults)["Tg"]["planted"]
    assert entry["default"]["setting"] == {"k": 15}
    assert abs(entry["default"]["observed"] - 0.3) < 1e-6


# --------------------------------------------------------------------------- the grades
@pytest.mark.parametrize("observed, chance, verdict, expected", [
    (0.9, 0.2, "PASS", "reliable"),
    (0.21, 0.2, "FAIL", "no skill"),
])
def test_a_grade_follows_the_numbers(observed, chance, verdict, expected):
    rows = _runs(observed=observed, chance=chance, verdict=verdict)
    entry = C.summarise(C.runs_frame(rows))["Tg"]["planted"]
    assert entry["grade"] == expected


def test_too_few_conclusive_runs_is_untestable_rather_than_a_verdict():
    rows = _runs(seeds=(1,), targets=("a",)) + _runs(seeds=(2,), targets=("a",),
                                                    verdict="INCONCLUSIVE")
    entry = C.summarise(C.runs_frame(rows))["Tg"]["planted"]
    assert entry["grade"] == "untestable"


def test_a_strategy_that_only_works_tuned_is_graded_as_such():
    rows = []
    for setting, value in (({"k": 5}, 0.85), ({"k": 15}, 0.2)):
        rows += _runs(settings=setting, observed=value, chance=0.2,
                      verdict="PASS" if value > 0.5 else "FAIL")
    df = C.runs_frame(rows)
    entry = C.summarise(df, {("Tg", "planted"): json.dumps({"k": 15})})["Tg"]["planted"]
    assert entry["grade"] == "works when tuned"


# --------------------------------------------------------------------------- writing and reading
def test_what_is_written_can_be_read_back_and_carries_its_provenance(tmp_path):
    rows = _runs()
    summary = C.summarise(C.runs_frame(rows))
    path = C.write(summary, str(tmp_path / "cal.json"), meta={"date": "2026-09-26", "runs": 10})
    loaded = C.load(path)
    assert loaded["meta"]["date"] == "2026-09-26"
    assert C.entry("planted", "Tg", path)["grade"] == summary["Tg"]["planted"]["grade"]
    assert C.entry("nothing_like_this", "Tg", path) is None


def test_a_not_a_number_is_written_as_null_rather_than_as_nan(tmp_path):
    """`NaN` is not JSON, and a file with a bare NaN in it cannot be read by anything else."""
    summary = C.summarise(C.runs_frame(_runs(seeds=(1,), targets=("a",))))
    path = C.write(summary, str(tmp_path / "cal.json"))
    text = open(path).read()
    assert "NaN" not in text and "Infinity" not in text
    json.loads(text)


def test_a_missing_calibration_reads_as_empty_rather_than_raising(tmp_path):
    missing = str(tmp_path / "not-here.json")
    assert C.load(missing)["organisms"] == {}
    assert C.entry("anything", "Tg", missing) is None
    assert C.tuned_settings("anything", "Tg", missing) == {}
    assert C.sentence("anything", "Tg", missing) == ""


def test_the_tuned_settings_offered_to_the_form_leave_the_target_to_the_user(tmp_path):
    rows = _runs(settings={"k": 5})
    path = C.write(C.summarise(C.runs_frame(rows)), str(tmp_path / "cal.json"))
    tuned = C.tuned_settings("planted", "Tg", path)
    assert tuned == {"k": 5}
    assert "target" not in tuned and "exclude" not in tuned


def test_the_sentence_names_the_grade_and_the_scale(tmp_path):
    path = C.write(C.summarise(C.runs_frame(_runs())), str(tmp_path / "cal.json"))
    text = C.sentence("planted", "Tg", path)
    assert "RELIABLE" in text and "skill" in text and "0 = the shuffled-data null" in text


def test_the_markdown_table_has_one_row_per_calibrated_strategy(tmp_path):
    from starplast import strategies as S
    key = S.catalog()[0].key
    rows = _runs(strategy=key)
    path = C.write(C.summarise(C.runs_frame(rows)), str(tmp_path / "cal.json"))
    table = C.markdown_table("Tg", path)
    assert table.count("\n") == 3                       # header, rule, one strategy
    assert S.get(key).title in table
    assert C.markdown_table("Pf", path).startswith("_No calibration")


# --------------------------------------------------------------------------- the shipped file
@pytest.fixture(scope="module")
def shipped():
    return C.load()


def test_the_shipped_calibration_says_when_and_how_it_was_measured(shipped):
    if not shipped["organisms"]:
        pytest.skip("no calibration shipped on this machine")
    meta = shipped["meta"]
    assert meta.get("date") and meta.get("runs", 0) > 100
    assert "skill" in meta.get("method", "") and meta.get("run_dir")


def test_every_shipped_strategy_carries_a_grade_the_module_knows(shipped):
    if not shipped["organisms"]:
        pytest.skip("no calibration shipped on this machine")
    for organism, entries in shipped["organisms"].items():
        for key, entry in entries.items():
            assert entry["grade"] in C.GRADES, (organism, key, entry["grade"])
            assert entry["runs"] >= 1
            low, high = entry["default"].get("skill_low"), entry["default"].get("skill_high")
            if low is not None and high is not None:
                assert low <= high


def test_the_shipped_calibration_covers_the_strategies_it_claims(shipped):
    from starplast import strategies as S
    if not shipped["organisms"]:
        pytest.skip("no calibration shipped on this machine")
    known = {s.key for s in S.catalog()}
    for organism, entries in shipped["organisms"].items():
        assert set(entries) <= known, set(entries) - known
