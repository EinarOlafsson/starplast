"""The track record: every labelled gene held out exactly once, and the arithmetic of the levels.

The properties worth holding are the ones that make the record trustworthy rather than merely
present: complete coverage, folds that never split an orthogroup, abstentions kept apart from
errors, and no percentage printed from too few genes to mean one.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture(scope="module")
def ledger(planted):
    return T.evaluate(planted, "feature_knn", "compartment")


def test_every_labelled_gene_is_held_out_exactly_once(planted, ledger):
    truth = planted.truth("compartment")
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    expected = set(np.flatnonzero(kept.notna().to_numpy()).tolist())
    assert set(ledger["gene"]) == expected, "a labelled gene was missed or scored twice"
    assert not ledger["gene"].duplicated().any()


def test_a_fold_never_splits_an_orthogroup(planted):
    """A gene recovered through a paralog that stayed visible is not the inference being tested."""
    truth = planted.truth("compartment")
    labelled = np.flatnonzero(truth.notna().to_numpy())
    folds = T.folds_by_group(planted, labelled)
    groups = np.asarray(planted.groups())[labelled]
    per_group = pd.DataFrame({"group": groups, "fold": folds}).groupby("group")["fold"].nunique()
    assert (per_group == 1).all(), "an orthogroup was split across folds"


def test_the_folds_are_balanced_enough_to_be_useful(planted):
    truth = planted.truth("compartment")
    labelled = np.flatnonzero(truth.notna().to_numpy())
    counts = pd.Series(T.folds_by_group(planted, labelled)).value_counts()
    assert len(counts) == T.FOLDS
    assert counts.max() <= 3 * counts.min(), "one fold swallowed the genes"


def test_an_abstention_is_neither_right_nor_wrong(ledger):
    for _, row in ledger.iterrows():
        if row["abstained"]:
            assert row["prediction"] is None and row["correct"] is None
        else:
            assert isinstance(row["prediction"], str) and row["correct"] in (True, False)
    pooled = T.summary(ledger, "strategy").iloc[0]
    assert pooled["answered"] == pooled["genes"] - pooled["abstained"]
    assert pooled["right"] + pooled["wrong"] == pooled["answered"]


def test_the_rate_is_over_answered_genes_and_carries_its_interval(ledger):
    row = T.summary(ledger, "strategy").iloc[0]
    assert row["rate"] == pytest.approx(row["right"] / row["answered"])
    assert row["rate_low"] <= row["rate"] <= row["rate_high"]
    assert bool(row["enough"]) == bool(row["answered"] >= T.MIN_FOR_RATE)


def test_too_few_answered_genes_is_said_rather_than_shown_as_a_percentage():
    rows = [{"organism": "Tg", "strategy": "s", "target": "t", "setting_key": "", "seed": 0,
             "fold": 0, "mode": "together", "gene": i, "gene_id": f"g{i}", "truth": "a",
             "prediction": "a" if i < 2 else None, "correct": True if i < 2 else None,
             "abstained": i >= 2, "support": 1.0} for i in range(6)]
    row = T.summary(pd.DataFrame(rows), "strategy").iloc[0]
    assert row["answered"] == 2 and not bool(row["enough"])


def test_every_level_pools_the_same_genes(ledger):
    for level in ("gene", "class", "target", "strategy"):
        out = T.summary(ledger, level)
        assert len(out), level
        assert out["genes"].sum() == len(ledger), level
        assert {"rate", "rate_low", "rate_high", "enough"} <= set(out.columns), level


def test_the_class_level_says_what_a_label_is_confused_with(planted):
    led = T.evaluate(planted, "feature_knn", "compartment")
    wrong = led[led["correct"] == False]                        # noqa: E712
    classes = T.summary(led, "class")
    if len(wrong):
        row = classes[classes["truth"] == wrong["truth"].iloc[0]].iloc[0]
        assert row["confused_with"], "a class with errors must name what it was called instead"


def test_a_gene_can_be_asked_about_by_name(planted, ledger):
    gene = ledger["gene_id"].iloc[0]
    one = T.for_gene(ledger, gene)
    assert len(one) == 1 and set(one["strategy"]) == {"feature_knn"}
    line = T.sentence(ledger, gene)
    assert line.startswith("Known ") and "recovered by" in line or "abstained" in line


def test_a_strategy_that_cannot_speak_about_a_target_returns_nothing_rather_than_guessing(planted):
    """A refusal is a fact about the pairing; it must not become an empty-looking success."""
    out = T.evaluate(planted, "holdout_search", "compartment")      # a map walk, not a gene caller
    assert out.empty and list(out.columns) == list(T.COLUMNS)


def test_the_supported_strategies_are_real_and_call_labels():
    keys = T.supported()
    assert len(keys) >= 10
    known = {s.key: s for s in S.catalog()}
    for key in keys:
        assert key in known, key
        assert known[key].task == "label calls", key
