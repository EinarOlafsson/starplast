"""Validating an annotation: hide labels we already have, and see whether they come back.

The point of this module is to attach an error rate to a candidate list. These tests are mostly
about the ways it could produce a flattering number: scoring against genes that were visible when
the cluster was chosen, scoring a label the embedding was built from, or reporting one average
across categories that behave completely differently.
"""
import numpy as np
import pandas as pd
import pytest

from starplast import validate as V


def _labels_and_truth(n_per=40):
    """Three clean clusters, each one category, plus unlabelled genes mixed in."""
    labels = np.repeat([0, 1, 2], n_per)
    truth = pd.Series(["a"] * n_per + ["b"] * n_per + ["c"] * n_per)
    return labels, truth


def test_a_perfect_clustering_recovers_what_was_hidden():
    labels, truth = _labels_and_truth()
    r = V.masked_recovery(labels, truth, "a", folds=3, seed=0)
    assert r.recall == pytest.approx(1.0), "hidden genes were in the right cluster and not found"
    assert r.precision == pytest.approx(1.0)


def test_a_clustering_that_ignores_the_category_scores_badly():
    """The check that the metric can fail: labels assigned at random must not look like recovery."""
    rng = np.random.default_rng(0)
    truth = pd.Series(["a"] * 60 + ["b"] * 60)
    labels = rng.integers(0, 3, size=120)
    r = V.masked_recovery(labels, truth, "a", folds=5, seed=1)
    assert r.precision < 0.6, f"random clusters scored precision {r.precision}"


def test_hidden_genes_do_not_influence_which_cluster_is_chosen():
    """Otherwise the test marks its own homework: the cluster would be picked using the answer."""
    labels = np.array([0] * 10 + [1] * 10)
    truth = pd.Series(["a"] * 10 + ["a"] * 10)
    r = V.masked_recovery(labels, truth, "a", folds=4, seed=0, hold_frac=0.5)
    # Both clusters are all "a", so whichever is chosen from the visible half, the hidden half is
    # split across both -- recall must be well under 1, not 1.
    assert r.recall < 1.0


def test_scoring_a_label_the_embedding_used_is_refused():
    """A cluster matching a feature the map was built from is circular and measures nothing.

    The test is on the label COLUMN. It used to be on the category VALUE -- "a" against a list of
    column names -- which is never true of real inputs, and the application passed the embedding's
    BLOCK names, so the guard compared a compartment against "expression_summary" and never fired
    once in the tab that exists to prevent exactly this."""
    labels, truth = _labels_and_truth()
    truth = truth.rename("compartment")
    with pytest.raises(ValueError, match="circular"):
        V.masked_recovery(labels, truth, "a", used_columns=["compartment"])
    with pytest.raises(ValueError, match="circular"):
        V.validate_all(labels, truth, used_columns=["compartment"], min_size=2)


def test_scoring_a_label_against_its_own_experiment_s_other_outputs_is_refused():
    """Naming the label column is not enough. A map built on `lopit_prob_map` is a map built on
    hyperLOPIT's own output, and scoring hyperLOPIT's compartment against it is the same circularity
    by a longer route -- the one the search guard caught and this one did not, in the tab whose job
    is to say how much to believe a cluster."""
    labels, truth = _labels_and_truth()
    truth = truth.rename("compartment")
    with pytest.raises(ValueError, match="same experiment"):
        V.masked_recovery(labels, truth, "a", used_columns=["lopit_prob_map", "expr_tachy"])


def test_scoring_a_label_against_another_estimate_of_the_same_thing_is_refused():
    """`ortholopit_label` is localization transferred from another species: a different experiment,
    the same quantity, and no more valid as evidence about the map."""
    labels, truth = _labels_and_truth()
    truth = truth.rename("compartment")
    with pytest.raises(ValueError, match="same thing"):
        V.masked_recovery(labels, truth, "a", used_columns=["ortholopit_label"])


def test_an_unrelated_column_does_not_trip_the_wider_guard():
    """It has to refuse the right maps and only those: a guard that refused everything would make
    the tab useless while looking careful."""
    labels, truth = _labels_and_truth()
    truth = truth.rename("compartment")
    r = V.masked_recovery(labels, truth, "a", used_columns=["expr_tachy", "fit_invitro_hff"])
    assert r.precision >= 0.0


def test_the_column_can_be_named_when_the_series_does_not_carry_it():
    """A Series sliced out of a frame keeps its name; one built by hand may not, and the caller
    knows what it is."""
    labels, truth = _labels_and_truth()
    with pytest.raises(ValueError, match="circular"):
        V.masked_recovery(labels, truth, "a", used_columns=["compartment"],
                          target_column="compartment")


def test_a_label_the_embedding_never_saw_is_scored_normally():
    """The guard must not refuse everything: it fires on the column the map used, and nothing else."""
    labels, truth = _labels_and_truth()
    truth = truth.rename("cellcycle_phase")
    r = V.masked_recovery(labels, truth, "a", used_columns=["compartment", "expression_summary"],
                          folds=2)
    assert r.folds, "a held-out label was refused"


def test_a_category_with_too_few_genes_says_so_rather_than_scoring():
    labels = np.array([0, 0, 1])
    truth = pd.Series(["rare", "rare", "other"])
    r = V.masked_recovery(labels, truth, "rare", folds=3)
    assert r.folds == []
    assert "too few" in r.note
    assert np.isnan(r.f1)


def test_noise_is_never_chosen_as_the_annotating_cluster():
    """"It is in the noise" is not an annotation."""
    labels = np.array([-1] * 30 + [0] * 10)
    truth = pd.Series(["a"] * 40)
    r = V.masked_recovery(labels, truth, "a", folds=3, seed=0)
    assert all(f.cluster >= 0 for f in r.folds)


def test_absence_labels_are_never_treated_as_a_category():
    labels = np.repeat([0, 1], 40)
    truth = pd.Series(["unassigned"] * 40 + ["a"] * 40)
    t = V.validate_all(labels, truth, folds=3, min_size=5)
    assert "unassigned" not in set(t.category)


def test_every_category_is_reported_separately():
    """One global number would hide that a method works for one compartment and not another."""
    labels, truth = _labels_and_truth()
    t = V.validate_all(labels, truth, folds=3, min_size=5)
    assert set(t.category) == {"a", "b", "c"}
    assert list(t.f1) == sorted(t.f1, reverse=True), "not ranked"


def test_a_table_with_nothing_scorable_has_the_right_columns():
    """An empty frame with no columns breaks every caller that reads one, and the empty one must
    have the same columns as a full one -- including `refit`, without which a reader cannot tell a
    number that re-embedded per fold from one that did not."""
    labels = np.array([0, 1])
    truth = pd.Series(["a", "b"])
    t = V.validate_all(labels, truth, min_size=100)
    assert list(t.columns) == ["category", "n_labelled", "n_folds",
                               "precision", "recall", "f1", "refit", "note"]
    full = V.validate_all(np.repeat([0, 1], 30), pd.Series(["a"] * 30 + ["b"] * 30),
                          min_size=5, folds=2)
    assert list(full.columns) == list(t.columns)


def test_candidates_carry_the_numbers_that_say_how_much_to_believe_them():
    """A candidate must never travel without its cluster's composition."""
    labels = np.array([0] * 10)
    truth = pd.Series(["a"] * 6 + ["b"] * 2 + ["unassigned"] * 2)
    genes = [f"g{i}" for i in range(10)]
    c = V.candidates(labels, truth, "a", 0, genes)
    assert list(c.gene_id) == ["g8", "g9"], "only the unlabelled members are candidates"
    assert c.cluster_frac_category.iloc[0] == pytest.approx(0.6)
    assert c.cluster_frac_contradicting.iloc[0] == pytest.approx(0.2)
    # And what that fraction has to be read against. A cluster that is 10% dense granules sounds
    # like something until you notice that 9% of every labelled gene is.
    # Over the same denominator as the cluster fraction -- every gene, not only the labelled ones --
    # because the two are meant to be divided by each other.
    assert c.overall_frac_category.iloc[0] == pytest.approx(0.6)
    assert c.enrichment.iloc[0] == pytest.approx(1.0)


def test_the_summary_says_whether_it_re_embedded():
    """A number that did not re-fit is a weaker claim than one that did, and must say so."""
    labels, truth = _labels_and_truth()
    r = V.masked_recovery(labels, truth, "a", folds=2)
    assert "hidden only from scoring" in r.summary()
    r.refit = True
    assert "re-embedded per fold" in r.summary()


def test_missing_values_are_absence_not_a_category():
    labels = np.repeat([0, 1], 30)
    truth = pd.Series([None] * 30 + ["a"] * 30, dtype="object")
    t = V.validate_all(labels, truth, folds=2, min_size=5)
    assert set(t.category) <= {"a"}


# --------------------------------------------------------------------------- re-fitting
def test_asking_to_re_embed_without_a_way_to_do_it_is_an_error():
    """`refit=True` used to change only the sentence the result prints, so a run that re-fit
    nothing described itself as "re-embedded per fold". A number that overstates how it was obtained
    is worse than a weaker number that describes itself accurately."""
    labels, truth = _labels_and_truth()
    with pytest.raises(ValueError, match="rebuild"):
        V.masked_recovery(labels, truth, "a", refit=True)


def test_re_fitting_scores_against_the_fold_s_own_clustering():
    """Otherwise re-fitting is a label on a result that was computed the old way."""
    labels, truth = _labels_and_truth()
    seen = []

    def rebuild(fold):
        seen.append(fold)
        # A clustering that puts every gene in one cluster: recall stays high, precision collapses,
        # and both differ from what the real clustering would have scored.
        return np.zeros(len(labels), dtype=int)

    r = V.masked_recovery(labels, truth, "a", folds=3, refit=True, rebuild=rebuild)
    assert seen == [0, 1, 2], "the fold's own map was not built"
    assert r.refit and "re-embedded per fold" in r.summary()
    assert r.precision < V.masked_recovery(labels, truth, "a", folds=3).precision


def test_a_rebuild_that_returns_the_wrong_number_of_labels_is_an_error():
    """Silently misaligned, every gene would be scored against another gene's cluster."""
    labels, truth = _labels_and_truth()
    with pytest.raises(ValueError, match="labels for"):
        V.masked_recovery(labels, truth, "a", folds=2, refit=True,
                          rebuild=lambda k: np.zeros(3, dtype=int))


def test_whether_it_re_fit_is_recorded_on_every_row():
    """The difference between a fresh map per fold and one fixed map is invisible in the numbers."""
    labels, truth = _labels_and_truth()
    t = V.validate_all(labels, truth, folds=2, min_size=5)
    assert set(t.refit) == {False}


# --------------------------------------------------------------------------- the annotating cluster
def test_the_candidate_cluster_is_the_one_the_score_was_measured_on():
    """A candidate list from a different cluster than the one that was validated carries a number
    that is about something else."""
    labels, truth = _labels_and_truth()
    assert V.best_cluster(labels, truth, "a") == 0
    assert V.best_cluster(labels, truth, "c") == 2


def test_a_category_in_no_cluster_has_no_annotating_cluster():
    labels = np.array([-1] * 10)
    truth = pd.Series(["a"] * 10)
    assert V.best_cluster(labels, truth, "a") is None


# --------------------------------------------------------------------------- orthogonal support
def _nodes_with_groups():
    return pd.DataFrame({
        "gene_id": ["g0", "g1", "g2", "g3"],
        "orthogroup": ["OG1", "OG1", "OG2", ""],
        "pfam_id": ["PF1;PF2", "PF2", "PF9", None],
        "compartment": ["a", "", "", ""],
    })


def test_support_counts_the_labelled_genes_a_candidate_shares_a_group_with():
    """Not proof -- a paralog of a dense-granule protein need not be one -- but independent of the
    map, and countable."""
    nodes = _nodes_with_groups()
    truth = pd.Series(["a", "", "", ""], name="compartment")
    s = V.orthogonal_support(nodes, ["g1", "g2"], truth, "a", used_columns=(),
                             log=lambda *_: None)
    row = s.set_index("gene_id")
    assert row.loc["g1", "shares_orthogroup"] == 1     # same OG as the labelled g0
    assert row.loc["g1", "shares_pfam_id"] == 1        # shares PF2 with g0
    assert row.loc["g2", "shares_orthogroup"] == 0
    assert row.loc["g2", "shares_pfam_id"] == 0


def test_a_labelled_gene_sharing_two_domains_is_counted_once():
    nodes = pd.DataFrame({"gene_id": ["g0", "g1"], "pfam_id": ["PF1;PF2", "PF1;PF2"]})
    truth = pd.Series(["a", ""], name="compartment")
    s = V.orthogonal_support(nodes, ["g1"], truth, "a", columns=("pfam_id",), log=lambda *_: None)
    assert s.shares_pfam_id.iloc[0] == 1


def test_a_column_that_fed_the_embedding_is_refused_rather_than_counted():
    """Its agreement would be the map agreeing with itself, which is what this tab exists to
    measure rather than to repeat."""
    said = []
    s = V.orthogonal_support(_nodes_with_groups(), ["g1"], pd.Series(["a", "", "", ""]), "a",
                             columns=("orthogroup", "pfam_id"), used_columns=("orthogroup",),
                             log=said.append)
    assert "shares_orthogroup" not in s.columns and "shares_pfam_id" in s.columns
    assert any("not independent evidence" in m for m in said)


def test_a_column_the_table_does_not_have_is_skipped_quietly():
    s = V.orthogonal_support(_nodes_with_groups(), ["g1"], pd.Series(["a", "", "", ""]), "a",
                             columns=("orthogroup", "not_a_column"), log=lambda *_: None)
    assert list(s.columns) == ["gene_id", "shares_orthogroup"]


def test_a_candidate_that_is_not_in_the_table_gets_zero_rather_than_a_crash():
    """Gene ids come from whatever rows the embedding covered, and a mismatch must not take the
    panel down mid-click."""
    s = V.orthogonal_support(_nodes_with_groups(), ["nonesuch"], pd.Series(["a", "", "", ""]), "a",
                             columns=("orthogroup",), log=lambda *_: None)
    assert s.shares_orthogroup.iloc[0] == 0


def test_zero_support_is_reported_as_zero_never_as_a_blank():
    """Most of this proteome has no labelled orthogroup neighbour, and a candidate with no
    independent support is exactly the one to be careful about."""
    s = V.orthogonal_support(_nodes_with_groups(), ["g3"], pd.Series(["a", "", "", ""]), "a",
                             log=lambda *_: None)
    assert (s.drop(columns="gene_id").fillna(-1) >= 0).all().all()


# --------------------------------------------------------------------------- nothing to score
def test_hiding_every_gene_of_a_category_scores_nothing_rather_than_guessing():
    """With none left visible there is no cluster a user could have picked, so there is nothing to
    score against -- and scoring it anyway would be scoring the hidden genes against themselves."""
    labels, truth = _labels_and_truth(n_per=10)
    r = V.masked_recovery(labels, truth, "a", folds=3, hold_frac=1.0)
    assert r.folds == []
    assert "not enough labeled genes" in r.summary()


def test_a_clustering_that_is_all_noise_annotates_nothing():
    """"It is in the noise" is not an annotation, so no fold is scored rather than one being
    invented from the noise bucket."""
    labels = np.full(40, -1)
    truth = pd.Series(["a"] * 40)
    assert V.masked_recovery(labels, truth, "a", folds=2).folds == []


def test_a_run_reports_each_category_as_it_goes():
    """It is the only sign of life during a long run -- and with `refit` on, a full embedding per
    fold per category, it is also where a stop takes effect."""
    labels, truth = _labels_and_truth()
    said = []
    V.validate_all(labels, truth, folds=2, min_size=5, log=said.append)
    assert any("validating 3 categories" in m for m in said)
    assert sum(m.strip().startswith(("1/3", "2/3", "3/3")) for m in said) == 3


def test_a_cluster_no_richer_in_the_category_than_the_proteome_shows_an_enrichment_of_one():
    """The number that says a candidate list is not worth having."""
    labels = np.zeros(20, int)
    truth = pd.Series(["a"] * 10 + ["b"] * 5 + ["unassigned"] * 5)
    c = V.candidates(labels, truth, "a", 0, [f"g{i}" for i in range(20)])
    assert c.enrichment.iloc[0] == pytest.approx(1.0)


def test_a_category_nothing_carries_does_not_divide_by_zero():
    labels = np.zeros(6, int)
    truth = pd.Series(["unassigned"] * 6)
    c = V.candidates(labels, truth, "a", 0, [f"g{i}" for i in range(6)])
    assert len(c) == 6 and np.isnan(c.enrichment.iloc[0])
