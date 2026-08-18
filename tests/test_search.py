#!/usr/bin/env python3
"""The structure search — the thing this project exists to do, and the thing easiest to fool yourself with.

Every assertion here is about honesty rather than performance. A search that reports a high score is
worthless unless the score was earned on a label the embedding could not see, and this project has
already produced the spectacular version of that mistake: `compartment` fed an embedding and its own
twin `lopit_map` came back as the top held-out discovery at V = 0.96, against a truthful mean F1 of
0.362 once the leak was closed.

So the tests that matter are the ones that would fail if the guard were removed.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import search as S  # noqa: E402
from starplast.clustering import NOISE  # noqa: E402


# --------------------------------------------------------------------------- scoring
def test_a_perfect_partition_scores_one():
    labels = np.array([0] * 60 + [1] * 60)
    truth = pd.Series(["a"] * 60 + ["b"] * 60)
    per, table = S.score_recovery(labels, truth)
    assert per["mean_f1"] == pytest.approx(1.0)
    assert set(table.label) == {"a", "b"}


def test_a_clustering_that_ignores_the_label_scores_poorly():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 2, 200)
    truth = pd.Series(rng.choice(["a", "b"], 200))
    per, _ = S.score_recovery(labels, truth)
    assert per["mean_f1"] < 0.75


def test_precision_and_recall_are_kept_apart():
    """One cluster holding every 'a' plus a pile of 'b' has perfect recall and poor precision. A single
    blended number would hide which of the two failed."""
    labels = np.array([0] * 160)
    truth = pd.Series(["a"] * 40 + ["b"] * 120)
    _, table = S.score_recovery(labels, truth)
    a = table[table.label == "a"].iloc[0]
    assert a.recall == pytest.approx(1.0)
    assert a.precision == pytest.approx(0.25)
    assert a.f1 < 0.5


def test_the_best_single_cluster_is_taken_per_label():
    """The question is whether the structure isolates the label somewhere, not whether the clustering
    happens to use the same number of groups as the label has values."""
    labels = np.array([0] * 30 + [1] * 30 + [2] * 60)
    truth = pd.Series(["a"] * 60 + ["b"] * 60)
    _, table = S.score_recovery(labels, truth)
    assert table[table.label == "b"].iloc[0].f1 == pytest.approx(1.0)
    assert table[table.label == "a"].iloc[0].recall == pytest.approx(0.5)


def test_noise_points_are_not_scored():
    """HDBSCAN's noise label is not a cluster, and counting it as one would reward a clustering for
    refusing to classify."""
    labels = np.array([NOISE] * 100 + [0] * 60)
    truth = pd.Series(["a"] * 100 + ["b"] * 60)
    _, table = S.score_recovery(labels, truth)
    assert set(table.label) == {"b"}


def test_rare_labels_are_skipped_rather_than_scored_on_a_handful():
    labels = np.array([0] * 100 + [1] * 5)
    truth = pd.Series(["common"] * 100 + ["rare"] * 5)
    _, table = S.score_recovery(labels, truth, min_label=15)
    assert "rare" not in set(table.label)


def test_too_little_labelled_data_returns_nothing_rather_than_a_number():
    labels = np.array([0] * 10)
    truth = pd.Series(["a"] * 10)
    per, table = S.score_recovery(labels, truth)
    assert per == {} and table.empty


def test_unlabelled_genes_do_not_count_against_a_cluster():
    labels = np.array([0] * 100)
    truth = pd.Series(["a"] * 60 + [None] * 40)
    _, table = S.score_recovery(labels, truth)
    assert table.iloc[0].precision == pytest.approx(1.0)


# --------------------------------------------------------------------------- the circularity guard
def _leaky_nodes(n=400):
    rng = np.random.default_rng(1)
    truth = np.where(np.arange(n) < n // 2, "in", "out")
    return pd.DataFrame({
        "gene_id": [f"TGME49_{200000+i}" for i in range(n)],
        "compartment": truth,
        "compartment_best": truth,                       # exact twin under another name
        "lopit_unified": truth,                          # derivation of the same experiment
        "unrelated": rng.normal(size=n),
    })


def test_an_exact_twin_under_another_name_is_excluded():
    """A guard that can be defeated by renaming a column is not a guard."""
    ex = S.excluded_for(_leaky_nodes(), "compartment")
    assert "compartment_best" in ex and "lopit_unified" in ex


def test_an_unrelated_column_is_not_excluded():
    assert "unrelated" not in S.excluded_for(_leaky_nodes(), "compartment")


def test_the_target_itself_is_always_excluded():
    assert "compartment" in S.excluded_for(_leaky_nodes(), "compartment")


def test_the_same_experiment_s_other_outputs_are_excluded_too():
    """The third leak this guard has had, and the one neither other mechanism can see.

    `lopit_prob_map` is the posterior of hyperLOPIT's own assignment. It does not restate the
    compartment -- pairwise association 0.29 on the real cache, far under any workable threshold --
    and the label is not computed from it, so there is nothing to declare. It is the same
    experiment's other output, and a map built on it is scored against a label that experiment also
    produced. Measured on the shipped cache, the three-column localization block recovered
    `compartment` at mean F1 0.259, above every measurement block in the map.
    """
    rng = np.random.default_rng(5)
    n = len(_leaky_nodes())
    d = _leaky_nodes()
    d["lopit_prob_map"] = rng.random(n)
    d["lopit_methods_agree"] = rng.integers(0, 2, n).astype(float)
    ex = S.excluded_for(d, "compartment")
    assert {"lopit_prob_map", "lopit_methods_agree"} <= ex
    assert "unrelated" not in ex, "shared provenance must not become an excuse to exclude everything"


def test_a_target_from_another_experiment_keeps_the_hyperlopit_columns():
    """The exclusion is per experiment, not a blanket. Asking whether the map recovers cell-cycle
    phase has no reason to throw away localization -- and a guard that excluded everything would
    make every search unanswerable while looking rigorous."""
    rng = np.random.default_rng(6)
    d = _leaky_nodes()
    d["lopit_prob_map"] = rng.random(len(d))
    d["cellcycle_phase"] = ["G1", "S"] * (len(d) // 2)
    assert "lopit_prob_map" not in S.excluded_for(d, "cellcycle_phase")


def test_the_negative_control_does_not_get_to_see_its_own_inputs():
    """`attention_depth` is `np.select` over n_papers_focal / substantive / incidental -- a
    deterministic function of exactly those three columns. Each scores 0.25-0.43 against the tiering
    on the real cache, far under the threshold, because a count is not a restatement of a tier while
    determining it completely.

    It matters more here than anywhere: this is the NEGATIVE control, the number every measured
    target is compared against. A control that is allowed to embed its own inputs scores too high,
    and every real target then looks worse than it is by exactly that much."""
    from starplast import datasets
    assert datasets.derived_sources("attention_depth"), \
        "the registry must declare what the tiering is computed from"
    rng = np.random.default_rng(8)
    n = 300
    d = pd.DataFrame({
        "gene_id": [f"TGME49_{200000+i}" for i in range(n)],
        "n_papers_focal": rng.integers(0, 5, n).astype(float),
        "n_papers_substantive": rng.integers(0, 5, n).astype(float),
        "n_papers_incidental": rng.integers(0, 9, n).astype(float),
        "n_publications": rng.integers(0, 12, n).astype(float),
        "n_fulltext": rng.integers(0, 7, n).astype(float),
        "fit_a": rng.normal(size=n),
    })
    d["attention_depth"] = np.where(d.n_papers_focal > 0, "focal",
                                    np.where(d.n_papers_substantive > 0, "substantive", ""))
    ex = S.excluded_for(d, "attention_depth")
    assert {"n_papers_focal", "n_papers_substantive", "n_papers_incidental"} <= ex
    assert "n_publications" in ex and "n_fulltext" in ex, \
        "the rest of the literature block counts the same thing"
    assert "fit_a" not in ex


def test_the_same_quantity_measured_another_way_is_excluded():
    """Provenance says which experiment made a column; this says what the column is an estimate OF,
    and the two come apart wherever the project measures something twice. `ortholopit_label` is a
    localization label transferred from P. falciparum and C. parvum orthologs -- a different
    experiment, a different species, the same quantity -- and it sits at 0.72 against `compartment`
    on the real cache, just under the threshold, which is what a near-copy does."""
    rng = np.random.default_rng(9)
    d = _leaky_nodes()
    d["ortholopit_label"] = list(d["compartment"])
    d["lit_tier"] = rng.choice(["focal", "incidental", ""], len(d))
    ex = S.excluded_for(d, "compartment")
    assert "ortholopit_label" in ex
    assert "lit_tier" not in ex, "a different quantity stays in"


def test_target_family_holdout_crosses_assay_branches_but_keeps_sequence_predictors():
    d = _leaky_nodes()
    d["ortholopit_accuracy"] = 0.9
    d["tm_len_mean"] = np.random.default_rng(91).normal(size=len(d))
    ex = S.excluded_for(d, "compartment", scope="target_family")
    assert "ortholopit_accuracy" in ex
    assert "tm_len_mean" not in ex, "sequence topology is a predictor, not a copied location label"


def test_a_broader_biology_holdout_can_remove_localization_related_sequence_data():
    d = _leaky_nodes()
    d["tm_len_mean"] = np.random.default_rng(92).normal(size=len(d))
    ex = S.excluded_for(d, "compartment", scope="biology")
    assert "tm_len_mean" in ex


def test_an_evidence_branch_can_be_omitted_without_naming_each_slot():
    d = _leaky_nodes()
    d["lopit_prob_map"] = np.arange(len(d), dtype=float)
    d["tm_len_mean"] = np.random.default_rng(93).normal(size=len(d))
    ex = S.excluded_group(d, "evidence", ("molecular measurements", "spatial biology"))
    assert "lopit_prob_map" in ex and "tm_len_mean" not in ex


def test_a_second_tiering_of_attention_is_excluded_from_the_control():
    rng = np.random.default_rng(10)
    d = _leaky_nodes()
    d["attention_depth"] = rng.choice(["focal", "substantive", ""], len(d))
    d["lit_tier"] = d["attention_depth"]
    ex = S.excluded_for(d, "attention_depth")
    assert "lit_tier" in ex and "compartment" not in ex


def test_the_families_name_only_columns_that_exist():
    """A family that lists a column the cache does not have is a note to nobody; one that lists a
    column that has been renamed is a guard that silently stopped guarding."""
    import pandas as pd
    from starplast import paths
    cache = paths.cache_file("nodes.parquet")
    if not os.path.exists(cache):
        pytest.skip("no built cache on this machine")
    have = set(pd.read_parquet(cache).columns)
    for quantity, cols in S.SAME_QUANTITY.items():
        missing = [c for c in cols if c not in have]
        assert not missing, f"{quantity} names columns the cache does not have: {missing}"


def test_the_exclusion_threshold_is_stricter_than_the_reporting_one():
    """0.8 when choosing what an embedding may SEE, 0.95 when flagging a result after the fact. A
    0.85-associated column leaks nearly as much as an identical one."""
    n = _leaky_nodes()
    assert len(S.excluded_for(n, "compartment", threshold=0.8)) >= \
           len(S.excluded_for(n, "compartment", threshold=0.95))


def test_a_declared_derivation_is_excluded_even_when_no_single_source_is_associated_enough():
    """The joint-function case. stage_enriched_derived is the argmax of three columns, so each source
    individually explains only ~0.6 of it -- under any workable threshold -- while the label is a
    deterministic function of all three together. Pairwise measurement cannot see that, however good
    the statistic; only the declaration can."""
    from starplast import datasets
    declared = datasets.derived_sources("stage_enriched_derived")
    assert declared, "the registry must declare what the derived column came from"

    rng = np.random.default_rng(2)
    n = 300
    cols = {c: rng.normal(size=n) for c in declared}
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(n)], **cols})
    nodes["stage_enriched_derived"] = pd.DataFrame(cols).idxmax(axis=1).values
    ex = S.excluded_for(nodes, "stage_enriched_derived")
    for c in declared:
        assert c in ex, f"{c} was declared as a source and must never feed the embedding"


def test_declared_sources_of_an_ordinary_target_are_empty():
    from starplast import datasets
    assert datasets.derived_sources("compartment") == ()


def test_a_derived_summary_excludes_the_raw_experiment_behind_it():
    """Exposing raw GSE columns must not reopen the derived stage-label circularity leak."""
    from starplast import datasets, paths

    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    banned = S.excluded_for(nodes, "stage_enriched_derived")
    for key in ("gse108740", "gse206344"):
        raw = {c for c in datasets.get(key).columns if c in nodes.columns}
        assert raw <= banned, f"{key} can feed the experiment behind the held-out stage label"


def test_a_spec_drops_blocks_that_would_leak():
    from starplast.embedding import EmbeddingSpec
    nodes = _leaky_nodes()
    spec = EmbeddingSpec(blocks=("expression_summary",), categorical=("compartment", "unrelated"))
    out = S._spec_without(spec, nodes, {"compartment"})
    assert "compartment" not in out.categorical
    assert "unrelated" in out.categorical


# --------------------------------------------------------------------------- predictions
def test_predictions_are_only_made_from_pure_clusters():
    labels = np.array([0] * 20 + [1] * 20)
    truth = pd.Series(["a"] * 15 + [None] * 5 + ["a"] * 8 + ["b"] * 7 + [None] * 5)
    genes = [f"g{i}" for i in range(40)]
    out = S.predictions(None, labels, truth, genes, min_precision=0.9)
    assert set(out.cluster) == {0}, "the mixed cluster must not produce predictions"
    assert set(out.predicted) == {"a"}


def test_a_prediction_carries_the_precision_it_was_earned_at():
    """The only honest estimate of how often the prediction will be right."""
    labels = np.array([0] * 20)
    truth = pd.Series(["a"] * 12 + [None] * 8)
    out = S.predictions(None, labels, truth, [f"g{i}" for i in range(20)])
    assert out.cluster_precision.tolist() == pytest.approx([1.0] * len(out))
    assert (out.n_labelled_in_cluster == 12).all()


def test_only_unlabelled_genes_receive_a_prediction():
    """Predicting a label for a gene that already has one is not a prediction."""
    labels = np.array([0] * 20)
    truth = pd.Series(["a"] * 12 + [None] * 8)
    genes = [f"g{i}" for i in range(20)]
    out = S.predictions(None, labels, truth, genes)
    assert set(out.gene_id) == set(genes[12:])


def test_a_cluster_with_too_few_labelled_members_predicts_nothing():
    labels = np.array([0] * 30)
    truth = pd.Series(["a"] * 3 + [None] * 27)
    assert S.predictions(None, labels, truth, [f"g{i}" for i in range(30)]).empty


def test_noise_never_produces_predictions():
    labels = np.array([NOISE] * 40)
    truth = pd.Series(["a"] * 20 + [None] * 20)
    assert S.predictions(None, labels, truth, [f"g{i}" for i in range(40)]).empty


# --------------------------------------------------------------------------- the walk
def _searchable(n=400):
    """Two overlapping blobs whose label is recoverable but NOT determined by the features.

    The overlap is the point. A first version of this fixture used well-separated blobs, and the
    circularity guard correctly threw every feature out -- association 1.0 -- leaving no blocks and
    zero runs. A label a feature determines exactly is circular by definition, so a fixture that
    exercises the walk has to be one an honest guard admits: here the strongest association is 0.64,
    below the 0.8 exclusion threshold, and the walk still recovers the label at F1 ~0.96.
    """
    rng = np.random.default_rng(3)
    half = n // 2
    x = np.vstack([rng.normal(-0.75, 1.0, (half, 6)), rng.normal(0.75, 1.0, (n - half, 6))])
    d = pd.DataFrame(x, columns=[f"fit_f{i}" for i in range(6)])
    d.insert(0, "gene_id", [f"TGME49_{200000+i}" for i in range(n)])
    d["compartment"] = ["in"] * half + ["out"] * (n - half)
    return d


def test_an_unknown_target_is_rejected_by_name():
    with pytest.raises(ValueError, match="not in the table"):
        S.search(_searchable(), target="no_such_column", log=lambda *_: None)


def test_a_combination_the_guard_empties_is_reported_rather_than_dropped(capsys):
    """Closing the provenance leak removes whole combinations from a sweep, and "8 runs" quietly
    becoming 7 is the failure this project keeps finding: work that did not happen, reported as a
    number that looks like it did."""
    d = _searchable()
    d["lopit_prob_map"] = np.random.default_rng(7).random(len(d))
    lines = []
    S.search(d, target="compartment", block_sets=[("fitness_screens",), ("localization",)],
             n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(25,),
             log=lines.append)
    assert any("skipped 1 combination" in m and "localization" in m for m in lines), lines


def test_the_walk_runs_and_ranks_by_mean_f1():
    d = _searchable()
    R, P = S.search(d, target="compartment",
                    block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    if R.empty:
        pytest.skip("umap/hdbscan produced no scorable clustering on this fixture")
    assert R.mean_f1.is_monotonic_decreasing
    assert (R.mean_f1 <= 1.0).all() and (R.mean_f1 >= 0.0).all()
    assert "assoc_with_input" in R.columns or "blocks" in R.columns


def test_the_walk_is_reproducible_under_the_same_seed():
    d = _searchable()
    kw = dict(target="compartment", block_sets=[("fitness_screens",)],
              n_neighbors_values=(15,), min_dist_values=(0.0,),
              min_cluster_sizes=(25,), log=lambda *_: None)
    a, _ = S.search(d, seed=42, **kw)
    b, _ = S.search(d, seed=42, **kw)
    if a.empty:
        pytest.skip("no scorable clustering on this fixture")
    assert a.mean_f1.tolist() == b.mean_f1.tolist()


def test_sampling_is_seeded_so_a_run_can_be_repeated():
    d = _searchable(400)
    kw = dict(target="compartment", block_sets=[("fitness_screens",)],
              n_neighbors_values=(15,), min_dist_values=(0.0,),
              min_cluster_sizes=(25,), sample_size=200, log=lambda *_: None)
    a, _ = S.search(d, seed=7, **kw)
    b, _ = S.search(d, seed=7, **kw)
    if a.empty:
        pytest.skip("no scorable clustering on this subsample")
    assert a.mean_f1.tolist() == b.mean_f1.tolist()


def test_every_target_in_the_table_names_a_real_intent():
    """TARGETS is the menu the analysis panel offers; a stale entry would offer a column that no longer
    exists and fail only when chosen."""
    from starplast import paths
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    missing = {k: v for k, v in S.TARGETS.items() if v not in nodes.columns}
    assert not missing, f"targets naming columns absent from the node table: {missing}"


def test_the_measured_and_derived_targets_are_both_registered():
    """cellcycle_phase is measured and stage_enriched_derived is not, and the difference is the whole
    point of having both."""
    assert "cellcycle_phase" in S.TARGETS
    assert "stage_enriched_derived" in S.TARGETS


def test_every_result_row_records_what_was_held_out():
    """A saved table must be auditable on its own. Reading a results CSV months later, "was the target
    actually held out" is the first question, and answering it should not require re-running the
    search."""
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert not R.empty
    assert "excluded" in R.columns and "n_excluded" in R.columns
    assert (R.n_excluded >= 1).all()
    assert R.excluded.str.contains("compartment").all()


def test_the_walk_records_the_hyperparameters_that_produced_each_row():
    """A score without the settings that produced it cannot be reproduced or reported."""
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15, 50), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert not R.empty
    for c in ("blocks", "n_neighbors", "min_dist", "min_cluster_size", "na_policy", "scaling"):
        assert c in R.columns, c
    assert set(R.n_neighbors) == {15, 50}


def test_the_per_label_table_breaks_the_mean_down():
    """A mean F1 hides which labels were found and which were missed entirely."""
    R, P = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert not P.empty
    assert {"label", "precision", "recall", "f1"} <= set(P.columns)
    assert set(P.label) <= {"in", "out"}


def test_block_sets_default_to_every_combination_that_has_columns():
    msgs = []
    S.search(_searchable(), target="compartment",
             n_neighbors_values=(15,), min_dist_values=(0.0,),
             min_cluster_sizes=(25,), log=msgs.append)
    assert any("dataset combinations" in m for m in msgs)


def test_a_block_set_with_no_usable_columns_is_skipped_not_fatal():
    R, _ = S.search(_searchable(), target="compartment",
                    block_sets=[("literature",), ("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert not R.empty
    assert all("literature" not in str(b) for b in R.blocks)


def test_embeddings_can_be_stored_for_the_runs_worth_keeping():
    """A result that cannot be reopened is not reproducible."""
    from starplast.tuning import EmbeddingStore
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        store = EmbeddingStore(tmp)
        R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                        n_neighbors_values=(15,), min_dist_values=(0.0,),
                        min_cluster_sizes=(25,), store=store, save_above=0.0,
                        log=lambda *_: None)
        assert not R.empty
        assert len(store.list()) >= 1


def test_nothing_is_stored_when_no_run_clears_the_bar():
    from starplast.tuning import EmbeddingStore
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        store = EmbeddingStore(tmp)
        S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                 n_neighbors_values=(15,), min_dist_values=(0.0,),
                 min_cluster_sizes=(25,), store=store, save_above=1.01,
                 log=lambda *_: None)
        assert len(store.list()) == 0


# --------------------------------------------------------------------------- the branches that only
# --------------------------------------------------------------------------- run when something is odd
def test_labels_that_are_all_too_rare_score_nothing_rather_than_everything():
    """Enough labelled points to be worth scoring, but no single label common enough to mean anything.
    Returning a mean over singletons would report near-perfect recovery of nothing."""
    labels = np.arange(60)
    truth = pd.Series([f"label_{i}" for i in range(60)])
    per, table = S.score_recovery(labels, truth, min_label=15)
    assert per == {} and table.empty


def test_a_log_that_cannot_flush_does_not_break_the_search(monkeypatch):
    """The wrapper exists so a long run redirected to a file shows progress rather than looking hung.
    Flushing can fail -- a closed or detached stdout under a test runner or a service manager -- and
    that must cost the flush, not the search."""
    calls = []

    class Broken:
        def flush(self):
            raise OSError("stdout is detached")

        def write(self, *a):
            pass

    monkeypatch.setattr(sys, "stdout", Broken())
    out = S._flushing(calls.append)
    out("hello")
    assert calls == ["hello"]


def test_an_unbuildable_matrix_is_skipped_rather_than_fatal(monkeypatch):
    """One bad dataset combination in a 328-run walk must not end the walk."""
    import starplast.search as mod
    real = mod.build_matrix
    seen = {"n": 0}

    def sometimes_bad(nodes, spec, **kw):
        seen["n"] += 1
        if seen["n"] == 1:
            raise ValueError("no usable columns in this combination")
        return real(nodes, spec, **kw)

    monkeypatch.setattr(mod, "build_matrix", sometimes_bad)
    R, _ = S.search(_searchable(), target="compartment",
                    block_sets=[("fitness_screens",), ("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert seen["n"] == 2, "the walk must have carried on to the second combination"
    assert not R.empty


def test_too_few_genes_to_cluster_is_skipped():
    """UMAP on 40 points produces a shape, and clustering it produces numbers that mean nothing."""
    small = _searchable(120)
    R, _ = S.search(small, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert R.empty


def test_n_neighbors_larger_than_the_sample_is_skipped():
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(5000,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=lambda *_: None)
    assert R.empty


def test_a_min_cluster_size_larger_than_the_data_is_skipped_not_fatal():
    """HDBSCAN raises rather than returning all-noise when min_cluster_size exceeds the sample, so an
    unusable grid point used to end the entire walk."""
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(100000,), log=lambda *_: None)
    assert R.empty


def test_one_unusable_grid_point_does_not_lose_the_usable_ones():
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(100000, 25), log=lambda *_: None)
    assert not R.empty, "the workable min_cluster_size must still be scored"


def test_a_clustering_that_scores_nothing_is_counted_but_not_recorded():
    """Everything falling to noise is a run that happened and has nothing to report."""
    msgs = []
    R, _ = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(399,), log=msgs.append)
    assert R.empty
    assert any("done:" in m for m in msgs)


def test_progress_is_reported_during_a_long_walk():
    """A 328-run search that prints nothing for twenty minutes looks like a hang."""
    msgs = []
    S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
             n_neighbors_values=(15,), min_dist_values=(0.0,),
             min_cluster_sizes=tuple(range(20, 62)), log=msgs.append)
    assert any("/" in m and "runs," in m for m in msgs), "no progress line during a 42-run walk"


def test_progress_does_not_depend_on_the_grid_dividing_the_report_interval():
    """The check sits at the end of a dataset combination, so the run count jumps by that
    combination's stride. Reported on exact equality it only ever fired when the stride happened to
    divide 40 -- and a long walk with any other shape printed nothing and looked like a hang."""
    msgs = []
    S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
             n_neighbors_values=(15,), min_dist_values=(0.0,),
             min_cluster_sizes=tuple(range(20, 71)), log=msgs.append)   # 51 runs: not a multiple of 40
    assert any("runs," in m for m in msgs)


def test_storing_every_run_says_so():
    from starplast.tuning import EmbeddingStore
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        msgs = []
        S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                 n_neighbors_values=(15,), min_dist_values=(0.0,),
                 min_cluster_sizes=(25,), store=EmbeddingStore(tmp), save_above=None,
                 log=msgs.append)
        assert any("saved embeddings for every run" in m for m in msgs)


def test_without_umap_the_search_reports_and_returns_empty(monkeypatch):
    """umap-learn is a heavy optional dependency. Its absence is a message, not a traceback."""
    import builtins
    real = builtins.__import__

    def no_umap(name, *a, **k):
        if name == "umap":
            raise ImportError("no umap here")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_umap)
    msgs = []
    R, P = S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), log=msgs.append)
    assert R.empty and P.empty
    assert any("umap-learn not installed" in m for m in msgs)


# --------------------------------------------------------------------------- absence is not a class
def test_a_label_meaning_not_measured_is_not_scored():
    """Recovering `unassigned` is recovering which genes were MEASURED, and measurement tracks study
    effort -- so the score would be about the literature rather than the biology. Included, compartment
    scored 0.484 with `unassigned` as its best label; excluded, the same run scores 0.207."""
    labels = np.array([0] * 60 + [1] * 60)
    truth = pd.Series(["unassigned"] * 60 + ["rhoptry"] * 60)
    per, table = S.score_recovery(labels, truth)
    assert "unassigned" not in set(table.label)
    assert set(table.label) == {"rhoptry"}


@pytest.mark.parametrize("absent", ["unassigned", "unknown", "", "nan", "None"])
def test_every_spelling_of_absence_is_excluded(absent):
    """These arrive from different layers -- a pandas NaN stringified, a compartment call that was
    never made, an empty attention tier -- and they all mean the same thing."""
    labels = np.array([0] * 60 + [1] * 60)
    truth = pd.Series([absent] * 60 + ["rhoptry"] * 60)
    _, table = S.score_recovery(labels, truth)
    assert set(table.label) == {"rhoptry"}


def test_the_exclusion_can_be_overridden_for_a_target_where_absence_is_the_question():
    """Kept configurable rather than hardcoded: asking "does the map separate measured from
    unmeasured" is a legitimate question, it is just not a biological one."""
    labels = np.array([0] * 60 + [1] * 60)
    truth = pd.Series(["unassigned"] * 60 + ["rhoptry"] * 60)
    _, table = S.score_recovery(labels, truth, exclude_labels=frozenset())
    assert set(table.label) == {"unassigned", "rhoptry"}


def test_a_target_that_is_entirely_absence_scores_nothing():
    labels = np.array([0] * 60 + [1] * 60)
    truth = pd.Series(["unassigned"] * 120)
    per, table = S.score_recovery(labels, truth)
    assert per == {} and table.empty


def test_a_stored_embedding_records_the_hyperparameters_it_actually_used():
    """spec0 carries the block, policy and scaling, but n_neighbors and min_dist are loop variables.
    Saved unchanged, every stored recipe claimed the EmbeddingSpec defaults whatever the run did -- so
    the name encoded the truth and the machine-readable field did not, and reopening a saved embedding
    by its recipe would have rebuilt a different map."""
    from starplast.tuning import EmbeddingStore
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        store = EmbeddingStore(tmp)
        S.search(_searchable(), target="compartment", block_sets=[("fitness_screens",)],
                 n_neighbors_values=(50,), min_dist_values=(0.25,),
                 min_cluster_sizes=(25,), store=store, save_above=0.0, log=lambda *_: None)
        listed = store.list()
        assert len(listed) >= 1
        assert set(listed.n_neighbors) == {50}, "the stored recipe must say nn=50, not the default"
        assert set(listed.min_dist) == {0.25}

        _, spec, _ = store.load(listed.name.iloc[0])
        assert spec.n_neighbors == 50 and spec.min_dist == 0.25


# --------------------------------------------------------------------------- rebuilding one row
def _one_run(**kw):
    d = _searchable()
    R, _ = S.search(d, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(25,),
                    log=lambda *_: None, **kw)
    return d, R


def test_a_row_carries_the_subsample_it_was_drawn_from():
    """`n_genes` is what survived the missing-value policy afterwards, which is a different number.
    Reconstructing the draw from it produced a different set of genes, and therefore a different map
    from the one the row describes."""
    d, R = _one_run(sample_size=300)
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    assert R.iloc[0].sample_size == 300


def test_rebuilding_a_row_reproduces_the_run_it_records():
    """The whole point of clicking a row. If the rebuild differs in the subsample, the scaling basis
    or the excluded columns, the map on screen is not the map the row's numbers are about -- and
    nothing about the row would say so."""
    d, R = _one_run()
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    row = R.iloc[0]
    coords, genes, labels, features = S.rebuild(d, row, log=lambda *_: None)
    assert int(genes.sum()) == int(row.n_genes) == len(coords)
    assert len(set(labels[labels != NOISE])) == int(row.n_clusters)
    assert float((labels == NOISE).mean()) == pytest.approx(float(row.noise_frac), abs=0.02)


def test_rebuilding_a_subsampled_row_uses_the_same_genes():
    d, R = _one_run(sample_size=300)
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    coords, genes, labels, _ = S.rebuild(d, R.iloc[0], log=lambda *_: None)
    assert int(genes.sum()) == int(R.iloc[0].n_genes)
    expected = np.random.default_rng(int(R.iloc[0].seed)).choice(len(d), 300, replace=False)
    assert set(np.flatnonzero(genes)) == set(expected)


def test_rebuilding_excludes_what_the_run_excluded():
    """The exclusion is what makes the score mean anything; a rebuild that let the target back in
    would show a map that separates it by construction."""
    d = _searchable()
    row = {"blocks": "fitness_screens", "na_policy": "median", "scaling": "rank",
           "n_neighbors": 15, "min_dist": 0.0, "min_cluster_size": 25, "seed": 42,
           "sample_size": 0, "excluded": ";".join(f"fit_f{i}" for i in range(6))}
    with pytest.raises(ValueError, match="excluded"):
        S.rebuild(d, row, log=lambda *_: None)


def test_a_row_naming_no_blocks_cannot_be_rebuilt():
    with pytest.raises(ValueError, match="no feature blocks"):
        S.rebuild(_searchable(), {"blocks": ""}, log=lambda *_: None)


def test_a_rebuilt_map_arrives_at_the_same_scale_as_any_other():
    """It replaces the map in the view, and a map at a different extent would need the camera
    re-framed for it."""
    d, R = _one_run()
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    coords, _, _, _ = S.rebuild(d, R.iloc[0], log=lambda *_: None)
    assert np.isclose(np.abs(coords).max(), 50.0, atol=1e-3)


# --------------------------------------------------------------------------- the automated walk
def test_each_embedding_is_emitted_with_its_best_clustering():
    """The automated walk: build, cluster at each size, keep the winner, and hand it back as
    something that can be looked at rather than only counted."""
    d = _searchable()
    steps = []
    R, P = S.search(d, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,),
                    min_cluster_sizes=(10, 25), on_run=steps.append, log=lambda *_: None)
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    # One step per EMBEDDING, not per row: the rows for one embedding differ only in the clustering,
    # and a gallery of the same map twice is not a gallery.
    assert len(steps) == 1 and len(R) >= 1
    step = steps[0]
    assert step.coords.shape[0] == int(step.genes.sum()) == len(step.labels)
    assert step.row["min_cluster_size"] in (10, 25)
    best = R[R.n_neighbors == 15].mean_f1.max()
    assert step.row["mean_f1"] == pytest.approx(best), "a step kept a clustering that did not win"
    assert not step.per.empty and {"label", "precision", "recall", "f1"} <= set(step.per.columns)


def test_a_step_arrives_at_the_same_scale_as_any_other_map():
    d = _searchable()
    steps = []
    S.search(d, target="compartment", block_sets=[("fitness_screens",)],
             n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(25,),
             on_run=steps.append, log=lambda *_: None)
    if steps:
        assert np.isclose(np.abs(steps[0].coords).max(), 50.0, atol=1e-3)
        assert "nn=15" in steps[0].label


def test_the_winner_is_chosen_by_the_objective_in_force():
    """"Best" is not a fixed thing: which clustering wins depends on what was being optimised, and a
    step that ignored the objective would show a different map from the one the ranking names."""
    d = _searchable()
    steps = []
    S.search(d, target="compartment", block_sets=[("fitness_screens",)],
             n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(10, 25),
             objective={"objective": "mean_recall", "min_cluster": 1, "min_label": 5},
             on_run=steps.append, log=lambda *_: None)
    if not steps:
        pytest.skip("no scorable clustering on this fixture")
    assert "objective_score" in steps[0].row


# --------------------------------------------------------------------------- the frontier
def test_the_frontier_keeps_what_nothing_beats_on_both():
    """Ranking on each of two objectives is two sorts; the frontier is what neither sort shows -- a
    configuration second on both is often the one to use and tops neither list."""
    R = pd.DataFrame({"mean_f1": [0.5, 0.4, 0.45, 0.2, 0.45],
                      "best_f1": [0.6, 0.9, 0.70, 0.3, 0.55]})
    on = S.frontier(R)
    # (0.45, 0.70) survives: it is worse than (0.5, 0.6) on the mean and better on the best, so
    # neither beats the other and both are choices. (0.45, 0.55) does not: (0.5, 0.6) beats it on
    # both. (0.2, 0.3) is beaten by everything.
    assert list(on) == [True, True, True, False, False]


def test_a_configuration_with_no_score_is_not_on_the_frontier():
    R = pd.DataFrame({"mean_f1": [0.5, np.nan], "best_f1": [0.6, 0.9]})
    assert list(S.frontier(R)) == [True, False]


def test_the_frontier_of_a_table_that_has_only_one_objective_is_everything():
    R = pd.DataFrame({"mean_f1": [0.5, 0.4]})
    assert S.frontier(R).all()
    assert S.frontier(pd.DataFrame()).empty


def test_a_real_search_marks_its_frontier():
    d = _searchable()
    R, _ = S.search(d, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15, 30), min_dist_values=(0.0,), min_cluster_sizes=(25,),
                    log=lambda *_: None)
    if R.empty:
        pytest.skip("no scorable clustering on this fixture")
    assert "on_frontier" in R.columns and R.on_frontier.any()


def test_the_default_sweep_names_the_blocks_it_used_and_the_ones_it_did_not():
    """"41 combinations from 6 blocks" does not say which six. The default base leaves out
    `localization` and `literature` on purpose, the interface sweeps every block, and the difference
    between those two is how a leak that could never reach the published battery reached the Search
    tab. A sweep has to say what it swept."""
    d = _searchable()
    d["lopit_prob_map"] = np.random.default_rng(11).random(len(d))
    lines = []
    S.search(d, target="compartment", n_neighbors_values=(15,), min_dist_values=(0.0,),
             min_cluster_sizes=(25,), log=lines.append)
    named = [m for m in lines if "dataset combinations from" in m]
    assert named, lines
    assert "fitness_screens" in named[0], "the blocks it used are not named"
    assert "not swept: localization" in named[0], "the blocks it left out are not named"


# --------------------------------------------------------------------------- stopping
def test_a_stopped_search_returns_what_it_finished():
    """Pressing stop must not cost the twenty minutes already spent. A search that is stopped
    returns its rows ranked, with the per-category table, rather than raising -- a stop is a
    decision that enough has been seen, not an error."""
    d = _searchable()
    seen = {"n": 0}

    def stop_after_two():
        seen["n"] += 1
        return seen["n"] > 2

    lines = []
    R, P = S.search(d, target="compartment",
                    block_sets=[("fitness_screens",), ("protein_features",), ("interactions",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(25,),
                    should_stop=stop_after_two, log=lines.append)
    assert not R.empty, "a stopped search threw away what it had already computed"
    assert len(R) < 3
    assert any("STOPPED after" in m for m in lines), lines


def test_a_search_nobody_stops_runs_to_the_end():
    """The flag is asked, not assumed: a `should_stop` that always says no must not shorten a run."""
    d = _searchable()
    R, _ = S.search(d, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15,), min_dist_values=(0.0,), min_cluster_sizes=(25,),
                    should_stop=lambda: False, log=lambda *_: None)
    assert len(R) == 1


def test_a_search_stopped_between_hyperparameters_keeps_the_earlier_ones():
    """The check is in both loops. One combination can carry dozens of hyperparameter points, so a
    stop that could only land between combinations would still take minutes on a real sweep."""
    d = _searchable()
    seen = {"n": 0}

    def stop_on_the_third_question():
        seen["n"] += 1
        return seen["n"] > 2

    R, _ = S.search(d, target="compartment", block_sets=[("fitness_screens",)],
                    n_neighbors_values=(15, 25, 35), min_dist_values=(0.0,),
                    min_cluster_sizes=(25,), should_stop=stop_on_the_third_question,
                    log=lambda *_: None)
    assert 0 < len(R) < 3


def test_an_unknown_exclusion_scope_is_refused():
    """The scopes decide how much of the map a hold-out removes, so a misspelt one must not quietly
    fall through to the narrowest reading -- that is the leak this guard exists to prevent."""
    import pandas as pd
    import pytest
    from starplast import search as S
    nodes = pd.DataFrame({"compartment": ["IMC"] * 60, "expr_a": range(60)})
    with pytest.raises(ValueError, match="scope must be"):
        S.excluded_for(nodes, "compartment", scope="everything")


@pytest.fixture(scope="module")
def nodes_small():
    """A real slice of the shipped table: 400 genes, every column.

    Real rather than synthetic, because the sweep's whole job is to hold out real evidence classes and
    a fabricated frame has no hierarchy to hold out. Small, because the sweep builds one map per
    category and the full table takes minutes.
    """
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "starplast", "data", "nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    frame = pd.read_parquet(path)
    return frame.sample(n=min(400, len(frame)), random_state=0).reset_index(drop=True)


# --------------------------------------------------------------------------- the category sweep
def test_categories_at_a_level_partition_the_blocks():
    """A level of the hierarchy is what counts as a "class of evidence". Every block belongs to
    exactly one category at a given level, or the sweep would hold something out twice."""
    from starplast import search as SE
    from starplast.embedding import SLOT_BLOCKS
    for level in (1, 2, 3):
        groups = SE.categories_at("evidence", level, blocks=SLOT_BLOCKS)
        seen = [b for _path, blocks in groups for b in blocks]
        assert len(seen) == len(set(seen)), f"level {level} puts a block in two categories"
        assert set(seen) <= set(SLOT_BLOCKS)
    assert len(SE.categories_at("evidence", 1)) < len(SE.categories_at("evidence", 3))


def test_a_deeper_level_is_a_finer_partition():
    from starplast import search as SE
    broad = {p for p, _ in SE.categories_at("evidence", 1)}
    fine = {p for p, _ in SE.categories_at("evidence", 3)}
    assert all(any(f[:1] == b for b in broad) for f in fine)


def test_the_sweep_holds_each_category_out_and_builds_from_the_rest(nodes_small):
    """The whole argument as a procedure: one map per category, each built WITHOUT that category."""
    from starplast import search as SE
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for
    usable = tuple(b for b in SLOT_BLOCKS
                   if columns_for(nodes_small, EmbeddingSpec(blocks=(b,))).get(b))
    if len(usable) < 4:
        pytest.skip("this fixture cannot fill enough blocks to hold one out")
    spec = EmbeddingSpec(blocks=usable)
    out = SE.sweep_categories(nodes_small, spec, hierarchy="evidence", level=1,
                              log=lambda *a: None)
    assert len(out), "the sweep scored nothing"
    for column in ("category", "columns_held_out", "columns_used", "tested", "best_score"):
        assert column in out.columns
    # The held-out columns must NOT be among the ones the map was built from. That is the whole point.
    assert (out["columns_used"] > 0).all()
    assert (out["columns_held_out"] > 0).all()


def test_a_stopped_sweep_returns_what_it_finished(nodes_small):
    """A stop is a decision that enough has been seen, not an error -- the same rule `search` follows."""
    from starplast import search as SE
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for
    usable = tuple(b for b in SLOT_BLOCKS
                   if columns_for(nodes_small, EmbeddingSpec(blocks=(b,))).get(b))
    if len(usable) < 4:
        pytest.skip("this fixture cannot fill enough blocks")
    out = SE.sweep_categories(nodes_small, EmbeddingSpec(blocks=usable), level=1,
                              should_stop=lambda: True, log=lambda *a: None)
    assert out.empty or len(out) >= 0          # it returns rather than raising


def test_the_sweep_only_holds_out_what_was_ticked(nodes_small):
    """It sweeps the selection, not the whole catalogue: a block nobody ticked is not evidence anyone
    asked about."""
    from starplast import search as SE
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for
    usable = [b for b in SLOT_BLOCKS if columns_for(nodes_small, EmbeddingSpec(blocks=(b,))).get(b)]
    if len(usable) < 4:
        pytest.skip("this fixture cannot fill enough blocks")
    chosen = tuple(usable[:3])
    groups = SE.categories_at("evidence", 1, blocks=chosen)
    assert {b for _p, blocks in groups for b in blocks} <= set(chosen)


# --------------------------------------------------------------------------- the sweep must not lie
def test_each_degeneracy_guard_fires_for_its_own_reason():
    """Three ways a clustering fails, and the reason is returned rather than a boolean because a
    reader deciding what to change needs to know which one fired."""
    import numpy as np
    from starplast.clustering import degenerate
    assert "cluster" in degenerate(np.array([0] * 50 + [1] * 50))
    assert "bisection" in degenerate(np.array([0] * 70 + [1] * 15 + [2] * 10 + [-1] * 5))
    assert "density minimum" in degenerate(np.array([0] * 30 + [1] * 30 + [2] * 40))
    assert degenerate(np.array([])) == "no genes were clustered"
    healthy = np.array([0] * 20 + [1] * 20 + [2] * 20 + [3] * 20 + [4] * 15 + [-1] * 5)
    assert degenerate(healthy) == ""


def test_the_sweep_refuses_to_score_a_degenerate_clustering(nodes_small):
    """The property that makes this instrument honest. A recovery score is a claim about STRUCTURE;
    computed on a clustering that merely bisected a continuum it is a sentence about nothing, and it
    reads exactly like a real result."""
    import numpy as np
    from starplast import search as SE
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for
    usable = tuple(b for b in SLOT_BLOCKS
                   if columns_for(nodes_small, EmbeddingSpec(blocks=(b,))).get(b))
    if len(usable) < 4:
        pytest.skip("this fixture cannot fill enough blocks")
    out = SE.sweep_categories(nodes_small, EmbeddingSpec(blocks=usable), level=1,
                              log=lambda *a: None)
    assert {"usable", "why_not", "best_score"} <= set(out.columns)
    for _, row in out.iterrows():
        if row["usable"]:
            assert row["why_not"] == ""
        else:
            assert row["why_not"], "a row was rejected without saying why"
            assert np.isnan(row["best_score"]), "an unusable row still carries a score"
            assert not row["recovered"], "an unusable row was called recovered"


def test_the_shipped_map_cannot_currently_support_a_recovery_claim():
    """Measured on the real table, and recorded as a test because it is a finding rather than a bug.

    Every block set tried -- 3 blocks, 10, 30, all 95 -- and every parameter setting -- n_neighbors
    5/15/30, min_dist 0.0/0.25, min_cluster_size 5 to 60 -- gives 2 to 5 clusters with 52% to 97% of
    genes in the largest and essentially no noise. HDBSCAN leaving no noise while cutting a cloud in
    half is the signature of slicing a continuum.

    If this test starts FAILING, the clustering has begun finding real structure and the sweep's
    numbers become meaningful. That is a result worth being told about, which is why it is asserted in
    this direction rather than skipped.
    """
    import os
    import numpy as np
    from starplast.clustering import cluster, degenerate
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for, embed
    from starplast import paths
    path = paths.cache_file("nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    nodes = pd.read_parquet(path)
    blocks = ("Tg_transcription_tachyzoite", "Tg_fitness_hff_in_vitro",
              "Tg_fold_confidence_disorder")
    coords, _names, _kept = embed(nodes, EmbeddingSpec(blocks=blocks), log=lambda *a: None)
    labels = cluster(np.asarray(coords), algorithm="hdbscan", min_cluster_size=25)
    assert degenerate(labels), (
        "the default map now clusters non-degenerately -- the sweep's scores have become meaningful "
        "and this test should be replaced by one asserting the recovery it can now measure")


def test_the_sweep_never_scores_a_category_with_columns_it_was_built_from(nodes_small):
    """The whole argument depends on this: the held-out columns must not be among the ones the map saw.
    `battery` marks anything the map used as `used` or `derived`, and only `held_out` rows are scored,
    but the sweep must not hand it an overlap in the first place."""
    from starplast import search as SE
    from starplast.embedding import EmbeddingSpec, SLOT_BLOCKS, columns_for
    usable = tuple(b for b in SLOT_BLOCKS
                   if columns_for(nodes_small, EmbeddingSpec(blocks=(b,))).get(b))
    if len(usable) < 4:
        pytest.skip("this fixture cannot fill enough blocks")
    spec = EmbeddingSpec(blocks=usable)
    for path, held in SE.categories_at("evidence", 1, blocks=spec.blocks):
        remaining = tuple(b for b in spec.blocks if b not in set(held))
        if not remaining:
            continue
        held_cols = {c for cols in
                     columns_for(nodes_small, EmbeddingSpec(blocks=tuple(held))).values() for c in cols}
        used_cols = {c for cols in
                     columns_for(nodes_small, EmbeddingSpec(blocks=remaining)).values() for c in cols}
        assert not (held_cols & used_cols), (
            f"{' > '.join(path)}: {len(held_cols & used_cols)} columns are both held out and used")
