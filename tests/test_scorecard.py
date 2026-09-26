"""Scorecards: the standard metrics must be the textbook ones, and must say when they are missing.

The arithmetic is checked against scikit-learn on planted numbers, so a metric in the application
means what the same name means everywhere else. The glossary is checked for completeness, because
a number the reader cannot interpret is not a result.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pandas as pd
import pytest
from sklearn import metrics as skm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import scorecard as SC  # noqa: E402
from starplast import strategies as S  # noqa: E402
from starplast import techniques as TQ  # noqa: E402


@pytest.fixture(scope="module")
def calls():
    rng = np.random.default_rng(0)
    t = pd.Series(rng.choice(list("abc"), 400, p=[0.6, 0.3, 0.1]))
    p = t.copy()
    flip = rng.random(400) < 0.3
    p[flip] = rng.choice(list("abc"), int(flip.sum()))
    p[rng.random(400) < 0.1] = np.nan
    scores = pd.DataFrame(rng.random((400, 3)), columns=list("abc"))
    scores["a"] += (t == "a") * 0.5
    return t, p, scores


def test_label_calls_match_scikit_learn(calls):
    t, p, scores = calls
    card = SC.label_calls(p, t, np.arange(len(t)), scores)
    pc = p.fillna("none")
    kw = dict(labels=list("abc"), zero_division=0)
    assert card["accuracy"] == pytest.approx((p == t).mean())
    assert card["kappa"] == pytest.approx(skm.cohen_kappa_score(t, pc))
    assert card["mcc"] == pytest.approx(skm.matthews_corrcoef(t, pc))
    assert card["macro_recall"] == pytest.approx(skm.recall_score(t, pc, average="macro", **kw))
    assert card["macro_precision"] == pytest.approx(
        skm.precision_score(t, pc, average="macro", **kw))
    assert card["macro_f1"] == pytest.approx(skm.f1_score(t, pc, average="macro", **kw))
    assert card["weighted_f1"] == pytest.approx(skm.f1_score(t, pc, average="weighted", **kw))
    assert card["macro_auroc"] == pytest.approx(
        np.mean([skm.roc_auc_score(t == c, scores[c]) for c in "abc"]))
    assert card["macro_auprc"] == pytest.approx(
        np.mean([skm.average_precision_score(t == c, scores[c]) for c in "abc"]))


def test_an_abstention_is_wrong_for_accuracy_but_not_for_precision_of_calls(calls):
    t, p, _ = calls
    card = SC.label_calls(p, t, np.arange(len(t)))
    called = p.notna()
    assert card["coverage"] == pytest.approx(called.mean())
    assert card["precision_of_calls"] == pytest.approx((p[called] == t[called]).mean())
    assert card["accuracy"] == pytest.approx(card["coverage"] * card["precision_of_calls"])


def test_macro_auroc_is_missing_rather_than_guessed_without_class_scores(calls):
    t, p, _ = calls
    card = SC.label_calls(p, t, np.arange(len(t)))
    assert math.isnan(card["macro_auroc"]) and math.isnan(card["macro_auprc"])
    assert list(card) == list(SC.TASKS[SC.T_LABEL].metrics)


def test_ranking_matches_scikit_learn_and_reads_against_prevalence():
    rng = np.random.default_rng(1)
    y = rng.random(2000) < 0.05
    s = rng.random(2000) + y * 0.8
    card = SC.ranking(s, y)
    assert card["auroc"] == pytest.approx(skm.roc_auc_score(y, s))
    assert card["auprc"] == pytest.approx(skm.average_precision_score(y, s))
    assert card["ndcg"] == pytest.approx(skm.ndcg_score([y], [s]))
    assert card["pauroc_10"] == pytest.approx(skm.roc_auc_score(y, s, max_fpr=0.1))
    assert card["prevalence"] == pytest.approx(y.mean())
    assert card["auprc_lift"] == pytest.approx(card["auprc"] / y.mean())


def test_a_random_ranking_sits_at_its_chance_levels():
    rng = np.random.default_rng(2)
    y = rng.random(20000) < 0.1
    card = SC.ranking(rng.random(20000), y)
    assert card["auroc"] == pytest.approx(0.5, abs=0.02)
    assert card["auprc_lift"] == pytest.approx(1.0, abs=0.1)
    assert card["recall_at_10pct"] == pytest.approx(0.1, abs=0.02)


def test_an_unscored_item_is_ranked_last_not_dropped():
    y = np.array([True, True, False, False])
    card = SC.ranking(np.array([np.nan, 0.9, 0.1, 0.2]), y)
    assert card["auroc"] == pytest.approx(0.5)


def test_values_match_scipy_and_scikit_learn():
    from scipy import stats
    rng = np.random.default_rng(3)
    y = rng.normal(size=300)
    p = y + rng.normal(size=300)
    card = SC.values(p, y)
    assert card["spearman"] == pytest.approx(stats.spearmanr(p, y).statistic)
    assert card["kendall"] == pytest.approx(stats.kendalltau(p, y).statistic)
    assert card["r2"] == pytest.approx(skm.r2_score(y, p))
    assert card["mae"] == pytest.approx(skm.mean_absolute_error(y, p))
    assert card["nrmse"] == pytest.approx(math.sqrt(skm.mean_squared_error(y, p)) / np.std(y))
    assert card["value_coverage"] == 1.0


def test_a_good_ranking_on_the_wrong_scale_has_high_rho_and_negative_r2():
    y = np.linspace(0, 1, 100)
    card = SC.values(10 * y + 5, y)
    assert card["spearman"] == pytest.approx(1.0) and card["r2"] < 0


def test_cluster_recovery_matches_scikit_learn_and_counts_noise_apart():
    rng = np.random.default_rng(4)
    cl = rng.integers(0, 4, 300)
    t = np.array(list("abcd"))[cl]
    t[rng.random(300) < 0.2] = "a"
    card = SC.cluster_recovery(cl, t, {"a": 0, "b": 1, "c": 2, "d": 3})
    assert card["ari"] == pytest.approx(skm.adjusted_rand_score(t, cl))
    assert card["nmi"] == pytest.approx(skm.normalized_mutual_info_score(t, cl))
    assert card["noise_share"] == 0.0
    noisy = SC.cluster_recovery(np.full(300, -1), t, {})
    assert noisy["noise_share"] == 1.0 and noisy["weighted_f1_clusters"] == 0.0


def test_set_retrieval_counts_what_it_says():
    card = SC.set_retrieval([1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8], np.arange(100))
    assert card["precision"] == pytest.approx(0.5) and card["recall"] == pytest.approx(0.6)
    assert card["jaccard"] == pytest.approx(3 / 8)
    assert card["fold_enrichment"] == pytest.approx(0.5 / 0.05)
    y = np.isin(np.arange(100), [4, 5, 6, 7, 8])
    r = np.isin(np.arange(100), [1, 2, 3, 4, 5, 6])
    assert card["set_mcc"] == pytest.approx(skm.matthews_corrcoef(y, r))


def test_replication_lift_is_rate_over_chance():
    card = SC.replication(6, 10, [0.2, 0.3])
    assert card["replication_rate"] == 0.6 and card["replication_lift"] == pytest.approx(2.4)


# --------------------------------------------------------------------------- the glossaries
def test_every_metric_says_what_it_is_what_chance_gives_and_how_to_read_it():
    for key, m in SC.METRICS.items():
        assert m.task in SC.TASKS, key
        assert len(m.definition.split()) >= 6 and m.range and m.chance, key
        assert len(m.reading.split()) >= 5, key
        assert SC.explain(key).startswith(m.label)


def test_every_task_lists_its_metrics_once_and_in_a_fixed_order():
    seen = []
    for task in SC.TASKS.values():
        assert task.metrics and len(set(task.metrics)) == len(task.metrics)
        seen += list(task.metrics)
    assert sorted(seen) == sorted(SC.METRICS)
    table = SC.glossary()
    assert {"definition", "chance", "reading"} <= set(table.columns)
    assert len(table) == len(SC.METRICS) + len(SC.VERDICT)


def test_every_strategy_names_its_task_and_techniques_from_the_glossaries():
    used = set()
    for s in S.catalog():
        assert s.task in SC.TASKS, s.key
        assert s.techniques and all(t in TQ.TECHNIQUES for t in s.techniques), s.key
        used |= set(s.techniques)
        assert s.metrics == SC.TASKS[s.task].metrics
        assert set(s.scorecard_table()["key"]) == set(s.metrics)
        assert list(s.techniques_table()["key"]) == list(s.techniques)
    assert used == set(TQ.TECHNIQUES), set(TQ.TECHNIQUES) - used


def test_every_technique_says_what_it_does_and_why():
    for t in TQ.TECHNIQUES.values():
        assert len(t.what.split()) >= 10 and len(t.why.split()) >= 6, t.key


def test_a_strategy_without_a_task_or_with_an_unknown_technique_is_refused():
    s = S.catalog()[0]
    with pytest.raises(ValueError, match="task"):
        S.register(S.Strategy(**{**s.__dict__, "key": "no_task", "task": "guessing"}))
    with pytest.raises(ValueError, match="techniques"):
        S.register(S.Strategy(**{**s.__dict__, "key": "no_tech", "techniques": ("magic",)}))


def test_the_api_exposes_both_glossaries():
    assert len(S.metrics()) == len(SC.METRICS) + len(SC.VERDICT)
    assert list(S.metrics("ranking")["key"]) == list(SC.TASKS["ranking"].metrics)
    assert len(S.techniques()) == len(TQ.TECHNIQUES)
    assert "task" in S.overview("Tg").columns


def test_a_test_result_carries_its_card_through_json():
    r = S.judge("x", "m", 0.8, [0.2, 0.3], min_effect=0.05, n_hidden=50, hidden="h",
                null_kind="n", t0=0.0, task=SC.T_VALUES,
                scorecard=SC.values(np.arange(50.0), np.arange(50.0)))
    d = r.to_dict()
    assert d["task"] == "values" and d["scorecard"]["spearman"] == pytest.approx(1.0)
    assert d["skill"] == pytest.approx((0.8 - 0.25) / 0.75)
    card = r.card()
    assert list(card["section"].unique()) == ["verdict", "values"]
