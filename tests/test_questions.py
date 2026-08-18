#!/usr/bin/env python3
"""The question catalogue: verdicts that are measured rather than asserted."""
import json
import os

import pandas as pd
import pytest

from starplast import questions as Q

RNA = ["biology", ["gene expression", "RNA abundance"]]
FITNESS = ["biology", ["parasite phenotype", "fitness and essentiality"]]
LOCALIZATION = ["biology", ["cell organization", "localization and topology"]]


@pytest.fixture(scope="module")
def nodes():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "starplast", "data", "nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    return pd.read_parquet(path).sample(n=400, random_state=0).reset_index(drop=True)


def test_a_missing_catalogue_is_empty_rather_than_fatal(tmp_path):
    assert Q.load(str(tmp_path / "nothing.json")) == []


def test_the_catalogue_round_trips(tmp_path):
    path = tmp_path / "q.json"
    path.write_text(json.dumps([{"question": "q", "verdict": "keep"}]))
    assert Q.load(str(path))[0]["question"] == "q"


def test_a_catalogue_entry_becomes_a_recipe_without_carrying_its_verdict():
    """A projection rather than a cast. `Recipe.from_dict` ignores what it does not know, and letting
    it do so here too would let a typo'd key ride along unnoticed in both directions."""
    r = Q.as_recipe({"question": "q", "inputs": [RNA], "holdout": "compartment",
                     "verdict": "keep", "reason": "because", "holdout_bins": 3})
    assert r.question == "q" and r.holdout_bins == 3
    assert not hasattr(r, "verdict") and not hasattr(r, "reason")


def test_a_question_its_author_kept_but_closure_refuses_is_reported_as_both(nodes):
    """The disagreement this table exists to surface: a hand-written verdict is a claim, and the
    closure's refusal is a measurement."""
    table = Q.verify(nodes, [
        {"question": "clean", "axis": "fitness and essentiality", "inputs": [FITNESS],
         "holdout": "compartment", "verdict": "keep"},
        {"question": "circular", "axis": "localisation and export", "inputs": [LOCALIZATION],
         "holdout": "compartment", "verdict": "keep", "reason": "looked fine to its author"},
    ])
    assert list(table.declared) == ["keep", "keep"]
    assert list(table.closure_ok) == [True, False]
    assert "restate a holdout" in table.loc[1, "refusal"]
    assert table.loc[0, "blocks"] > 0 and table.loc[1, "blocks"] == 0


def test_verify_counts_the_labelled_genes_a_holdout_actually_has(nodes):
    """The number that decides whether a question is scoreable at all, so it is counted rather than
    quoted from whoever wrote the question."""
    table = Q.verify(nodes, [{"question": "q", "inputs": [FITNESS], "holdout": "compartment"},
                             {"question": "binned", "inputs": [FITNESS], "holdout": "length",
                              "holdout_bins": 3}])
    assert 0 < table.loc[0, "holdout_labelled"] <= len(nodes)
    assert table.loc[1, "holdout_labelled"] == len(nodes)


def test_a_holdout_that_is_not_here_counts_nothing_rather_than_raising(nodes):
    table = Q.verify(nodes, [{"question": "q", "inputs": [FITNESS], "holdout": "no_such_column"}])
    assert table.loc[0, "holdout_labelled"] == 0 and not table.loc[0, "closure_ok"]


def test_only_questions_that_pass_closure_may_be_offered(nodes):
    catalogue = [
        {"question": "clean", "inputs": [FITNESS], "holdout": "compartment", "verdict": "keep",
         "shipped": True},
        {"question": "circular", "inputs": [LOCALIZATION], "holdout": "compartment",
         "verdict": "keep", "shipped": True},
        {"question": "sound but covered", "inputs": [FITNESS], "holdout": "compartment",
         "verdict": "keep"},
        {"question": "dropped", "inputs": [FITNESS], "holdout": "compartment", "verdict": "drop"},
    ]
    assert [q["question"] for q in Q.shipped(nodes, catalogue)] == ["clean"]
    # Without a table there is nothing to measure against, so the author's verdict is all there is.
    assert len(Q.shipped(None, catalogue)) == 2


def test_a_recipe_expected_to_fail_is_shipped_deliberately(nodes):
    """Instruction 45 asks for at least one recipe expected to FAIL, kept and labelled as such. A
    library where everything works is a library that has been fitted to its answers."""
    catalogue = [{"question": "known circular, kept as a demonstration", "inputs": [LOCALIZATION],
                  "holdout": "compartment", "verdict": "drop", "expect_refusal": True}]
    assert len(Q.shipped(nodes, catalogue)) == 1


def test_the_shipped_catalogue_is_consistent_with_itself():
    """Guards the real file: every entry names an axis the project declares, nothing is both kept
    and expected to fail, and no two shipped questions hold out the same label -- which is the rule
    that stops twenty recipes being one recipe run twenty times."""
    catalogue = Q.load()
    assert len(catalogue) >= 100, "the hundred questions instruction 45 asked for"
    for q in catalogue:
        assert q.get("axis", "") in Q.AXES, f"{q['question']}: unknown axis {q.get('axis')!r}"
        assert not (q.get("verdict") == "keep" and q.get("expect_refusal")), q["question"]
    shipped = [q for q in catalogue if q.get("shipped")]
    holdouts = [q["holdout"] for q in shipped]
    assert len(holdouts) == len(set(holdouts)), "two shipped questions hold out the same label"
    assert len(shipped) == 20, f"{len(shipped)} shipped, expected 20"
    assert sum(1 for q in catalogue if q.get("expect_refusal")) >= 1, \
        "no recipe is expected to fail; a library where everything works is fitted to its answers"
