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
