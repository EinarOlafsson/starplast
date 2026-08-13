"""What "good structure" means, and the degenerate solution that games each answer.

Every test here is about a way of being wrong. An objective that cannot be gamed by a two-cluster
solution or by singletons is worth having; one that can, and is used without a guard, will find the
degenerate answer, because a hyperparameter walk is an optimiser and optimisers find whatever the
objective actually rewards.
"""
import numpy as np
import pandas as pd
import pytest

from starplast import objectives as O

TRUTH = pd.Series(["a"] * 60 + ["b"] * 60 + ["c"] * 60)
PERFECT = np.repeat([0, 1, 2], 60)
ONE_BIG = np.zeros(180, int)
SINGLETONS = np.arange(180)


def s(labels, objective, **kw):
    return O.score(labels, TRUTH, objective=objective, min_cluster=1, **kw)["score"]


def test_a_perfect_clustering_scores_one_under_every_objective():
    for ob in ("mean_precision", "mean_recall", "mean_f1", "best_precision", "best_f1",
               "v_measure"):
        assert s(PERFECT, ob) == pytest.approx(1.0), ob


def test_one_giant_cluster_games_recall_and_nothing_else():
    """The failure that actually happened: on the full proteome the positive control scored 0.675,
    the highest of four targets, from a two-cluster solution whose per-label recalls were all 1.000.
    """
    assert s(ONE_BIG, "mean_recall") == pytest.approx(1.0)
    assert s(ONE_BIG, "mean_precision") < 0.5
    assert s(ONE_BIG, "mean_f1") < 0.75
    assert s(ONE_BIG, "v_measure") == pytest.approx(0.0, abs=1e-6)


def test_singletons_game_precision_and_nothing_else():
    """A cluster of one is perfectly pure, so any precision objective is maximised by shattering."""
    assert s(SINGLETONS, "mean_precision") == pytest.approx(1.0)
    assert s(SINGLETONS, "best_precision") == pytest.approx(1.0)
    assert s(SINGLETONS, "mean_recall") < 0.1
    assert s(SINGLETONS, "mean_f1") < 0.1


def test_the_minimum_cluster_size_is_what_stops_the_singleton_exploit():
    """Which is why the floor is enforced here rather than left to the caller."""
    assert O.score(SINGLETONS, TRUTH, objective="best_precision", min_cluster=10)["score"] == 0.0


def test_chance_corrected_agreement_resists_both_degenerate_solutions():
    """The one family immune to both, and therefore the check on any other objective."""
    assert O.agreement(PERFECT, TRUTH)["adjusted_rand"] == pytest.approx(1.0)
    assert O.agreement(ONE_BIG, TRUTH)["adjusted_rand"] == pytest.approx(0.0, abs=1e-6)
    assert abs(O.agreement(SINGLETONS, TRUTH)["adjusted_rand"]) < 0.05


def test_precision_over_clusters_and_recall_over_labels_are_different_questions():
    """Precision asks of a CLUSTER whether its members share a label; recall asks of a LABEL whether
    its genes share a cluster. A split label shows the difference."""
    split = np.array([0] * 30 + [1] * 30 + [2] * 60 + [3] * 60)   # "a" split across two pure clusters
    r = O.score(split, TRUTH, objective="mean_precision", min_cluster=1)["score"]
    q = O.score(split, TRUTH, objective="mean_recall", min_cluster=1)["score"]
    assert r == pytest.approx(1.0), "every cluster is pure, so precision is perfect"
    assert q < 1.0, "but 'a' is in two clusters, so recall is not"


def test_best_precision_finds_one_pure_cluster_among_a_mess():
    """The GRA-hunting objective: is there a pure cluster ANYWHERE, not is the map good."""
    rng = np.random.default_rng(0)
    labels = rng.integers(2, 8, size=180)
    labels[:40] = 0                                  # one clean cluster of "a"
    r = O.score(labels, TRUTH, objective="best_precision", min_cluster=10)
    assert r["score"] == pytest.approx(1.0)
    assert "is 100% a" in r["detail"]


def test_precision_at_recall_refuses_a_pure_but_useless_cluster():
    """Three co-located genes are precision 1.0 and nothing to annotate from."""
    labels = np.full(180, 9)
    labels[:12] = 0                                  # tiny, pure, holds only 20% of "a"
    loose = O.score(labels, TRUTH, objective="best_precision", min_cluster=10)["score"]
    strict = O.score(labels, TRUTH, objective="precision_at_recall", min_cluster=10,
                     min_recall=0.5)
    assert loose == pytest.approx(1.0), "the tiny pure cluster wins on precision alone"
    # The recall floor rejects it and falls back to a cluster that actually holds a label, which is
    # far less pure -- a worse-looking number describing something you could annotate from.
    assert strict["score"] < loose
    assert "cluster 9" in strict["detail"]


def test_precision_at_recall_reports_nothing_when_no_cluster_holds_enough():
    labels = np.arange(180) // 3                     # 60 tiny clusters, none holding much
    r = O.score(labels, TRUTH, objective="precision_at_recall", min_cluster=1, min_recall=0.5)
    assert r["score"] == 0.0
    assert "at least 50%" in r["detail"]


def test_macro_weighting_does_not_bury_a_rare_class():
    """search.score_recovery weights by label size, so on this proteome nucleus-chromatin at 769
    genes dominates and dense granules at 167 barely registers -- and someone hunting dense granules
    is optimising against themselves."""
    truth = pd.Series(["big"] * 300 + ["rare"] * 20)
    labels = np.array([0] * 300 + [1] * 10 + [0] * 10)     # "rare" only half recovered
    macro = O.score(labels, truth, objective="mean_f1", weighting="macro", min_cluster=1,
                    min_label=5)["score"]
    sized = O.score(labels, truth, objective="mean_f1", weighting="size", min_cluster=1,
                    min_label=5)["score"]
    assert sized > macro, "size weighting should hide the rare class's poor score"


def test_a_single_category_can_be_optimised_on_its_own():
    r = O.score(PERFECT, TRUTH, objective="best_f1", category="b", min_cluster=1)
    assert r["score"] == pytest.approx(1.0)
    assert "b:" in r["detail"]


def test_an_unknown_objective_is_refused_rather_than_silently_defaulted():
    with pytest.raises(ValueError, match="unknown objective"):
        O.score(PERFECT, TRUTH, objective="vibes")


def test_every_score_reports_coverage_and_cluster_count():
    """A score without them is not interpretable: the attention control's winner had 52% noise."""
    noisy = PERFECT.copy()
    noisy[:90] = -1
    r = O.score(noisy, TRUTH, objective="mean_f1", min_cluster=1)
    assert r["coverage"] == pytest.approx(0.5)
    assert r["n_clusters"] == 2


def test_absence_labels_are_never_a_target():
    truth = pd.Series(["unassigned"] * 60 + ["a"] * 60)
    P = O.pairs(np.repeat([0, 1], 60), truth, min_label=5, min_cluster=5)
    assert set(P.label) == {"a"}


def test_nothing_scorable_says_so_rather_than_returning_zero_silently():
    r = O.score(np.full(20, -1), TRUTH.head(20), objective="mean_f1")
    assert r["score"] == 0.0
    assert r["detail"]


def test_every_named_objective_is_documented_and_runs():
    for ob in O.OBJECTIVES:
        assert O.OBJECTIVES[ob].strip()
        O.score(PERFECT, TRUTH, objective=ob, min_cluster=1)


def test_asking_for_a_category_no_cluster_holds_says_so(): 
    """Naming a category that never clears the floors must explain itself, not score zero silently."""
    r = O.score(PERFECT, TRUTH, objective="best_f1", category="not-a-compartment", min_cluster=1)
    assert r["score"] == 0.0
    assert "no cluster holds enough" in r["detail"]


def test_agreement_and_v_measure_are_undefined_rather_than_zero_on_nothing():
    """Zero would read as "no agreement"; there is no agreement to measure at all."""
    empty = pd.Series(["unassigned"] * 10)
    assert np.isnan(O.v_measure(np.zeros(10, int), empty))
    a = O.agreement(np.zeros(10, int), empty)
    assert np.isnan(a["adjusted_rand"]) and np.isnan(a["adjusted_mutual_info"])


# --------------------------------------------------------------------------- the permutation null
def test_the_null_is_what_this_objective_gives_for_nothing_on_this_data():
    """Three balanced classes: shuffled labels recover about a third of each, so mean F1 lands near
    0.33 rather than at 0. Quoting a raw 0.33 as a result is the mistake this exists to prevent."""
    null = O.null_score(PERFECT, TRUTH, objective="mean_f1", min_cluster=1, n_permutations=10)
    assert 0.2 < null < 0.45, null


def test_the_null_rises_with_fewer_classes_which_is_the_whole_point():
    """A score is not comparable across targets without it. Two classes are easier to hit by chance
    than twenty-four, so the same 0.30 means different things -- and a search that ranks
    configurations across targets is otherwise comparing against different nulls in silence."""
    two = pd.Series(["a"] * 90 + ["b"] * 90)
    many = pd.Series([f"c{i % 24}" for i in range(180)])
    labels = np.repeat(np.arange(6), 30)
    n2 = O.null_score(labels, two, objective="mean_f1", min_cluster=1, n_permutations=10)
    n24 = O.null_score(labels, many, objective="mean_f1", min_cluster=1, n_permutations=10)
    assert n2 > n24


def test_the_null_is_reproducible_from_its_seed():
    """A number nobody can rebuild is not a result -- including this one."""
    kw = dict(objective="mean_f1", min_cluster=1, n_permutations=5, seed=7)
    assert O.null_score(PERFECT, TRUTH, **kw) == O.null_score(PERFECT, TRUTH, **kw)


def test_at_least_one_permutation_is_run_however_few_are_asked_for():
    """Zero permutations would make the null 0.0 -- an absent correction wearing the same name as a
    measured one, which is worse than not having the option."""
    assert O.null_score(PERFECT, TRUTH, objective="mean_f1", min_cluster=1, n_permutations=0) > 0.0


def test_adjusted_puts_chance_at_zero_and_perfect_at_one():
    r = O.adjusted(PERFECT, TRUTH, objective="mean_f1", min_cluster=1, n_permutations=10)
    assert r["adjusted"] == pytest.approx(1.0)
    assert r["n_labels"] == 3 and 0.0 < r["null"] < r["score"]


def test_structure_worse_than_shuffling_is_reported_negative_not_clipped():
    """Worse than chance is a real result about the map, and clipping it to zero would hide the one
    case where the answer is "this configuration is actively misleading"."""
    r = O.adjusted(np.tile([0, 1, 2], 60), TRUTH, objective="mean_f1", min_cluster=1,
                   n_permutations=10)
    assert r["adjusted"] < 0.0


def test_a_null_of_one_leaves_the_correction_undefined_rather_than_infinite():
    """`(score - null) / (1 - null)` divides by zero when chance already scores perfectly -- one
    class, where every shuffle is the same shuffle. NaN says "no correction is defined here"; a
    number would claim one."""
    one_class = pd.Series(["a"] * 180)
    r = O.adjusted(ONE_BIG, one_class, objective="mean_f1", min_cluster=1, n_permutations=3)
    assert r["null"] == pytest.approx(1.0) and np.isnan(r["adjusted"])


def test_unlabelled_genes_are_not_counted_as_a_class():
    """`n_labels` is what the null is there to make sense of, so it has to be the number of classes
    the score was computed over. Counting the unlabelled as a class says "3 classes" for a score
    computed over 2 -- and unassigned is the largest group in the real column."""
    truth = pd.Series(["a"] * 60 + ["b"] * 60 + [None] * 60)
    r = O.adjusted(PERFECT, truth, objective="mean_f1", min_cluster=1, n_permutations=5)
    assert r["n_labels"] == 2, "unlabelled is not a compartment"
    assert 0.0 < r["null"] < 1.0


def test_a_class_too_small_to_score_is_not_counted_either():
    """`pairs` skips a label below `min_label`, so counting it would describe the score as spread
    over more classes than it could possibly have used."""
    truth = pd.Series(["a"] * 88 + ["b"] * 88 + ["rare"] * 4)
    r = O.adjusted(PERFECT, truth, objective="mean_f1", min_cluster=1, min_label=15,
                   n_permutations=5)
    assert r["n_labels"] == 2
