"""Claims: certainty that means what it says, verification that is measured to be independent, and
claims that say plainly when they were not tested or lie outside what was measured.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import claims as C  # noqa: E402


def _ledger(truth, calls: dict, support=None, folds=5):
    """A ledger in the track record's shape: one row per (strategy, gene)."""
    rows = []
    n = len(truth)
    for key, pred in calls.items():
        for g in range(n):
            rows.append({"organism": "Tg", "strategy": key, "target": "t", "mode": "together",
                         "fold": g % folds, "gene": g, "gene_id": f"g{g}", "truth": truth[g],
                         "prediction": pred[g], "abstained": pred[g] is None,
                         "correct": None if pred[g] is None else pred[g] == truth[g],
                         "support": (support[key][g] if support else 0.5)})
    return pd.DataFrame(rows)


def _noisy(truth, rng, right_rate, classes):
    return [t if rng.random() < right_rate else rng.choice([c for c in classes if c != t])
            for t in truth]


def test_a_copy_of_the_generator_is_not_an_independent_check():
    rng = np.random.default_rng(0)
    classes = list("abcde")
    truth = list(rng.choice(classes, 2000))
    gen = _noisy(truth, rng, 0.5, classes)
    fresh = _noisy(truth, rng, 0.5, classes)
    led = _ledger(truth, {"gen": gen, "copy": list(gen), "fresh": fresh})
    same = C.shared_mistakes(led, "gen", "copy")
    other = C.shared_mistakes(led, "gen", "fresh")
    assert same["ratio"] > 3, "a copy repeats every mistake"
    assert other["ratio_low"] < 1 < other["ratio_high"], other
    assert not C.independent(led, "gen", "copy") and C.independent(led, "gen", "fresh")


def test_certainty_is_calibrated_when_support_is_the_true_rate():
    rng = np.random.default_rng(1)
    n = 3000
    support = rng.uniform(0.1, 0.95, n)
    truth = ["a"] * n
    pred = ["a" if rng.random() < s else "b" for s in support]
    led = _ledger(truth, {"gen": pred}, support={"gen": support})
    m = C.certainty_model(led, "gen")
    assert m.calibration_error < 0.05 and m.calibrated
    assert abs(float(m([0.8])[0]) - 0.8) < 0.1


def test_agreement_raises_and_disagreement_lowers_the_verified_certainty():
    rng = np.random.default_rng(2)
    classes = list("abcd")
    truth = list(rng.choice(classes, 3000))
    gen = _noisy(truth, rng, 0.5, classes)
    ver = _noisy(truth, rng, 0.7, classes)
    led = _ledger(truth, {"gen": gen, "ver": ver},
                  support={"gen": rng.uniform(0.3, 0.7, 3000), "ver": np.full(3000, 0.5)})
    m = C.verified_model(led, "gen", ("ver",))
    base = float(m.certainty([0.5])[0])
    up = float(m([base], {"ver": np.array(["agrees"])})[0])
    down = float(m([base], {"ver": np.array(["disagrees"])})[0])
    assert down < base < up
    v = C.verification(led, "gen", "ver")
    assert v["agrees"]["rate"] > v["base"]["rate"] > v["disagrees"]["rate"]


@pytest.fixture(scope="module")
def shipped():
    frame = C.shipped()
    if not len(frame):
        pytest.skip("claims have not been built here")
    return frame


def test_every_claim_states_whether_it_was_tested(shipped):
    assert set(shipped["status"].astype(str)) <= {"tested", "untested", "outside tested range"}
    outside = shipped[shipped["status"].astype(str) == "outside tested range"]
    assert outside["confidence"].isna().all(), "a claim outside the measured range carries a number"
    tested = shipped[shipped["status"].astype(str) == "tested"]
    verdicts = [c for c in shipped.columns if str(c).endswith(" verdict")]
    assert (tested[verdicts].astype(str).isin(["agrees", "disagrees"]).any(axis=1)).all()


def test_only_proven_recipes_made_claims():
    rec = C.recipes()
    if not len(rec):
        pytest.skip("not built here")
    assert (rec[~rec["proven"].astype(bool)]["claims"] == 0).all()


def test_the_lift_is_the_confidence_over_the_prior(shipped):
    ok = shipped["confidence"].notna()
    lift = shipped.loc[ok, "confidence"] / shipped.loc[ok, "prior"]
    assert np.allclose(lift, shipped.loc[ok, "lift"])


def test_discoveries_respect_their_thresholds(shipped):
    for organism in shipped["organism"].astype(str).unique():
        d = C.discoveries(organism, min_confidence=0.8, min_lift=2.0)
        assert (d["confidence"] >= 0.8).all() and (d["lift"] >= 2.0).all()
        assert set(d["status"].astype(str)) <= {"tested"}


def test_a_gene_card_shows_its_claims_and_opens_one(shipped):
    from PyQt6 import QtCore, QtWidgets
    _app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast import app as A
    w = A.Window()
    try:
        d = C.discoveries("Tg")
        if not len(d):
            pytest.skip("no Toxoplasma discoveries built")
        ids = w.nodes.gene_id.astype(str)
        row = int(ids[ids == str(d["gene_id"].iloc[0])].index[0])
        w.show_detail(row)
        html = w.detail.toHtml()
        assert "What Starplast claims" in html and "starplast://claim/" in html
        link = "starplast://claim/" + html.split("starplast://claim/", 1)[1].split('"', 1)[0]
        w._detail_link(QtCore.QUrl(link))
        shown = w.detail.toHtml()
        assert "Independent checks" in shown and "back to the gene" in shown
    finally:
        w.close()
