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


# --------------------------------------------------------------------------- sets held out together
def test_a_set_is_hidden_all_at_once(planted):
    """Every member of the set must be invisible while any of them is predicted."""
    sets = T.class_sets(planted, "compartment")
    name, members = next(iter(sets.items()))
    led = T.evaluate_sets(planted, "feature_knn", "compartment", {name: members})
    assert len(led) == len(members)
    assert set(led["mode"]) == {"set"} and set(led["set_name"]) == {name}
    assert set(led["set_size"]) == {len(members)}, "the rows disagree about how many were hidden"


def test_hiding_a_whole_class_leaves_nothing_to_copy(planted):
    """A vote among visible labels cannot name a class when every example of it is gone -- so the
    rate is zero, and what the genes are called INSTEAD is the finding."""
    sets = T.class_sets(planted, "compartment")
    led = T.evaluate_sets(planted, "feature_knn", "compartment", sets)
    summary = T.set_summary(led)
    assert (summary["right"] == 0).all()
    assert summary["placed_at"].str.len().gt(0).any()
    assert ((summary["together"] >= 0) & (summary["together"] <= 1)).all()


def test_the_degradation_curve_covers_the_sizes_asked_for(planted):
    sets = T.random_sets(planted, "compartment", sizes=(1, 5, 20), repeats=3, seed=2)
    assert len(sets) == 9
    assert sorted({len(v) for v in sets.values()}) == [1, 5, 20]
    led = T.evaluate_sets(planted, "feature_knn", "compartment", sets)
    assert len(led) == 3 * (1 + 5 + 20)


def test_set_rows_and_fold_rows_live_in_one_table(planted, ledger):
    sets = T.random_sets(planted, "compartment", sizes=(5,), repeats=2, seed=3)
    both = pd.concat([ledger, T.evaluate_sets(planted, "feature_knn", "compartment", sets)],
                     ignore_index=True)
    assert set(both["mode"]) == {"together", "set"}
    assert len(T.summary(both, "strategy")) == 1, "the levels pool across both kinds of hold-out"
    assert len(T.set_summary(both)) == 2, "and the set level reports only the sets"


# --------------------------------------------------------------------------- the shipped record
@pytest.fixture(scope="module")
def built():
    return T.shipped()


def test_the_shipped_record_covers_both_organisms_and_every_supported_strategy(built):
    if not len(built):
        pytest.skip("the track record has not been built on this machine")
    folds = built[built["mode"] == "together"]
    for organism in built["organism"].unique():
        here = folds[folds["organism"] == organism]
        assert here["strategy"].nunique() >= 10, organism
        # Every labelled gene of the target appears exactly once per strategy.
        per = here.groupby("strategy")["gene"].agg(["count", "nunique"])
        assert (per["count"] == per["nunique"]).all(), organism


def test_the_shipped_record_kept_the_set_holdouts(built):
    if not len(built):
        pytest.skip("not built here")
    sets = built[built["mode"] == "set"]
    assert len(sets), "no set hold-outs were kept"
    names = set(sets["set_name"].astype(str))
    assert any(n.startswith("random ") for n in names), "the degradation curve is missing"
    assert T.set_summary(built)["together"].notna().any()


def test_a_gene_card_section_is_html_or_nothing(built):
    if not len(built):
        pytest.skip("not built here")
    organism = str(built["organism"].iloc[0])
    gene = str(built[built["mode"] == "together"]["gene_id"].iloc[0])
    html = T.gene_html(gene, organism)
    assert "If this gene were unknown" in html and "<table" in html
    assert T.gene_html("no such gene", organism) == "", "an unknown gene adds no section"


def test_the_class_level_names_the_best_strategy_and_the_confusions(built):
    if not len(built):
        pytest.skip("not built here")
    folds = built[built["mode"] == "together"]
    organism = str(folds["organism"].iloc[0])
    target = str(folds["target"].iloc[0])
    label = str(folds["truth"].value_counts().index[0])
    html = T.class_html(target, label, organism)
    assert label in html and "<table" in html and "rate [95%]" in html
    assert T.class_html(target, "no such class", organism) == ""


def test_a_strategy_line_names_where_it_is_weakest_but_only_where_judgeable(built):
    if not len(built):
        pytest.skip("not built here")
    organism = str(built["organism"].iloc[0])
    strategy = str(built[built["mode"] == "together"]["strategy"].iloc[0])
    line = T.weakest(strategy, organism)
    assert "held-out genes" in line
    # Nothing is called weakest on a count too small to judge.
    rows = T.summary(built[(built["mode"] == "together")
                           & (built["strategy"] == strategy)], "class")
    for row in rows[~rows["enough"].astype(bool)].itertuples():
        assert f"{row.truth} (" not in line, row.truth
