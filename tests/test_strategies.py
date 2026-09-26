"""Strategies: the catalogue is complete, the machinery is right, and every self-test means something.

Three layers, because each catches what the others cannot:

* **the catalogue** -- thirty or more strategies, each with the prose a person needs (a tooltip, an
  explanation, a walkthrough, a description of its test) and parameters that explain themselves;
* **the machinery** -- hiding labels, the metrics, the nulls, the context's leakage guard;
* **the self-tests themselves**, on a planted organism where the answer is known. Every strategy
  must PASS there -- the mechanism can find what was planted -- and must NOT pass on null tables
  built with every column and edge dealt out at random. A self-test that passes on noise is a
  self-test of nothing, and this is the check that there is no such test in the catalogue.

The null check allows one chance pass in three tables: a 95th-percentile bar lets 5% of true nulls
through by design, and requiring zero would make the suite fail on a fair coin.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import strategies as S  # noqa: E402
from starplast.jobs import Stopped  # noqa: E402

CATALOG = S.catalog()
KEYS = [s.key for s in CATALOG]
#: Strategies whose planted structure needs a bigger table than the default, and the settings
#: that point them at it. Pattern-5 strategies split the genes in half and need enough of each
#: half for a finding to be significant; paralog divergence needs enough labelled paralog pairs.
BIG = {"split_clusters": {"b": "cellcycle_phase"},
       "conjunctions": {"b": "stage_enriched_derived", "min_cluster_size": 8},
       "paralog_divergence": {}}


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture(scope="module")
def planted_big():
    return S.planted_context(n=1500)


@pytest.fixture(scope="module")
def nulls():
    return [(S.planted_context(null=True, seed=s), S.planted_context(n=1500, null=True, seed=s))
            for s in (7, 8, 9)]


def _ctx(key, small, big):
    return big if key in BIG else small


def _genes_for(ctx, strategy):
    """A gene list for strategies that take one: a planted category's members."""
    if any(p.kind == "genes" for p in strategy.params):
        cat, members = S.example_set(ctx, "compartment")
        return {"genes": "\n".join(ctx.gene_ids[members]), "exclude": "compartment"}
    return {}


# =========================================================================== the catalogue
def test_there_are_at_least_thirty_strategies_with_unique_keys_and_numbers():
    assert len(CATALOG) >= 30
    assert len(set(KEYS)) == len(KEYS)
    assert [s.number for s in CATALOG] == list(range(1, len(CATALOG) + 1))


def test_the_two_founding_questions_are_the_first_two_strategies():
    """Hold a category out and search; hand over a list and hunt for its cluster."""
    assert KEYS[:2] == ["holdout_search", "geneset_hunt"]
    assert any(p.kind == "category" for p in CATALOG[0].params)
    assert any(p.kind == "genes" for p in CATALOG[1].params)


def test_the_families_are_grouped_and_each_holds_several_strategies():
    fams = S.families()
    assert len(fams) == 8
    for f in fams:
        members = [s for s in CATALOG if s.family == f]
        assert len(members) >= 2, f
        numbers = [s.number for s in members]
        assert numbers == list(range(numbers[0], numbers[0] + len(numbers))), f


@pytest.mark.parametrize("key", KEYS)
def test_every_strategy_explains_itself(key):
    s = S.get(key)
    assert len(s.tooltip.split()) >= 25, "a tooltip says what it infers and why"
    assert s.question.endswith("?")
    paras = [p for p in s.explanation.split("\n\n") if p.strip()]
    assert len(paras) >= 2 and len(s.explanation.split()) >= 120
    assert len(s.walkthrough) >= 4
    assert "null" in s.test_description.lower() and "pass" in s.test_description.lower()
    assert s.cost in ("seconds", "a minute", "minutes")
    assert callable(s.runner) and callable(s.tester)
    text = s.help_text()
    assert s.title in text and "Walkthrough" in text and "How it is tested" in text


@pytest.mark.parametrize("key", KEYS)
def test_every_parameter_says_why_it_exists(key):
    """The same bar the application sets for every control: a reason, not a name."""
    for p in S.get(key).params:
        assert p.kind in S.PARAM_KINDS
        assert len(p.tip.split()) >= 15, f"{key}.{p.name}"
        assert p.label


def test_a_second_registration_under_one_key_is_refused():
    s = CATALOG[0]
    with pytest.raises(ValueError):
        S.register(s)


def test_an_unknown_parameter_kind_is_refused():
    s = CATALOG[0]
    bad = S.Strategy(**{**s.__dict__, "key": "not_registered", "params": (
        S.Param("x", "colour", "x", "a tip long enough to satisfy nobody in particular today"),)})
    with pytest.raises(ValueError, match="unknown parameter kind"):
        S.register(bad)


def test_an_unknown_strategy_names_the_ones_that_exist():
    with pytest.raises(KeyError, match="holdout_search"):
        S.get("no_such_strategy")


def test_settings_refuse_a_parameter_the_strategy_does_not_take(planted):
    with pytest.raises(ValueError, match="has no parameter"):
        S.get("feature_knn").settings(planted, colour="blue")


def test_defaults_are_resolved_for_the_context(planted):
    d = S.get("feature_knn").defaults(planted)
    assert d["target"] == "compartment" and d["k"] == 15


# =========================================================================== small numerics
def test_auroc_ranks_ties_and_missing_scores():
    assert S.auroc([1, 2, 3, 4], [0, 0, 1, 1]) == 1.0
    assert S.auroc([4, 3, 2, 1], [0, 0, 1, 1]) == 0.0
    assert S.auroc([1, 1, 1, 1], [0, 1, 0, 1]) == 0.5
    assert S.auroc([np.nan, 5, 1, 2], [1, 1, 0, 0]) == 0.5      # a missing score ranks last
    assert np.isnan(S.auroc([1, 2], [1, 1]))
    assert S.auroc([np.nan, np.nan], [1, 0]) == 0.5


def test_spearman_needs_five_varying_pairs():
    assert np.isnan(S.spearman([1, 2, 3], [1, 2, 3]))
    assert np.isnan(S.spearman([1, 1, 1, 1, 1], [1, 2, 3, 4, 5]))
    assert S.spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]) == pytest.approx(1.0)


def test_bh_is_monotone_and_bounded():
    q = S.bh([0.01, 0.04, 0.03, 0.5])
    assert list(q) == pytest.approx([0.04, 0.04 * 4 / 3 * 3 / 3, 0.04, 0.5], abs=0.02)
    assert (q <= 1).all() and len(S.bh([])) == 0


def test_hypergeometric_edges():
    assert S.hypergeom_sf(0, 10, 5, 100) == 1.0
    assert S.hypergeom_sf(3, 0, 5, 100) == 1.0
    assert S.hypergeom_sf(5, 5, 5, 100) < 1e-6


def test_grids_parse_text_and_lists():
    assert S.parse_grid("15, 50,") == (15.0, 50.0)
    assert S.parse_grid("eom, leaf", str) == ("eom", "leaf")
    assert S.parse_grid([1, 2], int) == (1, 2)


def test_expand_writes_values_into_every_gene():
    out = S.expand(np.array([3, 4]), [1, 3], 5)
    assert list(out) == [-1, 3, -1, 4, -1]
    assert list(S.expand(np.array([], dtype=int), [], 3)) == [-1, -1, -1]


# =========================================================================== hiding and scoring
def _labels(n=200, seed=0):
    rng = np.random.default_rng(seed)
    v = rng.choice(["a", "b", "c", "rare"], size=n, p=[0.4, 0.3, 0.27, 0.03])
    s = pd.Series(v, dtype=object)
    s.iloc[:20] = np.nan
    return s


def test_hide_is_stratified_and_drops_classes_too_small_to_split():
    t = _labels()
    vis, hidden = S.hide(t, 0.25, seed=1)
    assert "rare" not in set(vis.dropna()) and "rare" not in set(t.iloc[hidden])
    for c in ("a", "b", "c"):
        n = int((t == c).sum())
        assert int((t.iloc[hidden] == c).sum()) == round(0.25 * n)
    assert vis.iloc[hidden].isna().all()


def test_hide_keeps_groups_whole():
    t = _labels()
    groups = np.array([f"g{i // 2}" for i in range(len(t))])       # pairs, like paralogs
    vis, hidden = S.hide(t, 0.3, seed=2, groups=groups)
    hid = set(hidden)
    for i in hidden:
        mate = i + 1 if i % 2 == 0 else i - 1
        if isinstance(t.iloc[mate], str) and t.iloc[mate] != "rare":
            assert mate in hid, "a paralog was left visible beside its hidden copy"


def test_shuffling_keeps_the_labels_and_where_they_are_missing():
    t = _labels()
    out = S.shuffled(t, np.random.default_rng(0))
    assert out.isna().equals(t.isna())
    assert sorted(out.dropna()) == sorted(t.dropna())


def test_an_abstention_is_a_miss_but_not_a_wrong_call():
    truth = pd.Series(["a", "b", "a", "b"], dtype=object)
    pred = pd.Series(["a", np.nan, "b", "b"], dtype=object)
    assert S.correct_rate(pred, truth, [0, 1, 2, 3]) == 0.5
    assert S.call_precision(pred, truth, [0, 1, 2, 3]) == (pytest.approx(2 / 3), 3)
    assert np.isnan(S.correct_rate(pred, truth, []))
    assert np.isnan(S.call_precision(pd.Series([np.nan] * 4, dtype=object), truth, [0, 1])[0])
    table = S.per_class(pred, truth, [0, 1, 2, 3])
    assert set(table["class"]) == {"a", "b"}
    assert table.set_index("class").loc["b", "called"] == 2


def test_the_analytic_null_is_the_product_of_the_class_shares():
    truth = pd.Series(["a"] * 50 + ["b"] * 50, dtype=object)
    pred = pd.Series(["a"] * 100, dtype=object)
    mean, sd = S.analytic_null(pred, truth, np.arange(100))
    assert mean == pytest.approx(0.5) and sd > 0
    assert np.isnan(S.analytic_null(pred, truth, [])[0])


def test_neighbour_votes_never_count_a_gene_as_its_own_neighbour():
    X = np.array([[0.0], [0.01], [5.0], [5.01]])
    vis = pd.Series(["a", "b", "b", "b"], dtype=object)
    pred, share = S.knn_vote(X, vis, k=1)
    assert pred.iloc[0] == "b"            # its nearest OTHER labelled gene
    empty, _ = S.knn_vote(X, pd.Series([np.nan] * 4, dtype=object))
    assert empty.isna().all()


def test_a_gene_with_no_labelled_partner_is_not_called():
    import scipy.sparse as sp
    A = sp.csr_matrix(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float))
    pred, support = S.graph_vote(A, pd.Series(["x", np.nan, np.nan], dtype=object))
    assert pred.iloc[1] == "x" and pd.isna(pred.iloc[2])
    none, _ = S.graph_vote(A, pd.Series([np.nan] * 3, dtype=object))
    assert none.isna().all()


def test_propagation_reaches_along_edges_and_nowhere_else():
    import scipy.sparse as sp
    A = sp.csr_matrix(np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 0]], float))
    d = np.asarray(A.sum(axis=1)).ravel()
    inv = np.divide(1.0, np.sqrt(d), out=np.zeros_like(d), where=d > 0)
    op = sp.diags(inv) @ A @ sp.diags(inv)
    pred, strength = S.propagate(op, pd.Series(["x", np.nan, np.nan, np.nan], dtype=object))
    assert pred.iloc[2] == "x" and pd.isna(pred.iloc[3])
    assert S.propagate(op, pd.Series([np.nan] * 4, dtype=object))[0].isna().all()
    field = S.diffuse(op, [0])
    assert field[1] > field[3] == 0


def test_a_set_is_scored_against_its_best_single_cluster():
    labels = np.array([0, 0, 0, 1, 1, 1, -1])
    members = np.array([1, 1, 0, 0, 0, 0, 1], bool)
    k, prec, rec, f1 = S.set_f1(labels, members)
    assert k == 0 and prec == pytest.approx(2 / 3) and rec == pytest.approx(2 / 3)
    assert S.set_f1(labels, np.zeros(7, bool))[0] == S.NOISE
    assert S.set_f1(np.full(3, -1), np.ones(3, bool))[3] == 0.0


def test_cluster_majorities_need_enough_labelled_members():
    labels = np.array([0, 0, 0, 1, 1, -1])
    vis = pd.Series(["a", "a", "b", "c", np.nan, "a"], dtype=object)
    assert S.cluster_majority(labels, vis, 0.5, 2) == {0: ("a", pytest.approx(2 / 3), 3)}
    pred, table = S.predict_from_clusters(labels, np.arange(6), vis, 6, 0.5, 2)
    assert list(pred.iloc[:3]) == ["a", "a", "a"] and pd.isna(pred.iloc[3])
    assert len(table) == 1


def test_a_knn_graph_is_symmetric_and_normalised():
    X = np.random.default_rng(0).normal(size=(30, 3))
    op = S.knn_operator(X, k=4)
    assert abs(op - op.T).max() < 1e-12 and op.diagonal().sum() == 0


# =========================================================================== verdicts
def _result(observed, null, **kw):
    import time
    return S.judge("k", "metric", observed, null, min_effect=kw.pop("min_effect", 0.05),
                   n_hidden=kw.pop("n_hidden", 100), hidden="h", null_kind="n",
                   t0=time.monotonic(), **kw)


def test_a_pass_needs_the_bar_and_the_margin():
    assert _result(0.5, [0.1, 0.12, 0.11]).verdict == "PASS"
    assert _result(0.14, [0.1, 0.12, 0.11]).verdict == "FAIL"          # over the bar, under margin
    assert _result(0.11, [0.1, 0.12, 0.11]).verdict == "FAIL"
    assert _result(0.5, [0.1], n_hidden=3).verdict == "INCONCLUSIVE"
    assert _result(float("nan"), [0.1]).verdict == "INCONCLUSIVE"
    assert _result(0.5, [0.1], n_hidden=3, min_hidden=3).verdict == "PASS"


def test_the_analytic_null_has_a_bar_and_a_p_value():
    r = _result(0.6, [], analytic=(0.2, 0.05))
    assert r.null_high == pytest.approx(0.2 + 1.645 * 0.05)
    assert r.p_value < 1e-10 and r.passed
    assert _result(0.1, [], analytic=(0.2, 0.0)).p_value == 1.0
    assert _result(0.1, [], analytic=(0.2, 0.05), quantile=80.0).null_high == pytest.approx(
        0.2 + 0.842 * 0.05)


def test_summaries_lead_with_the_verdict_and_serialise():
    for r in (_result(0.5, [0.1, 0.2]), _result(0.5, [0.1], n_hidden=2, note="why")):
        text = r.summary()
        assert text.split(" -- ")[0] in ("PASS", "FAIL", "INCONCLUSIVE")
        json.dumps(r.to_dict())
    assert "why" in _result(0.5, [0.1], n_hidden=2, note="why").summary()
    assert _result(0.5, []).null_sd == 0.0 and np.isnan(_result(0.5, []).null_mean)
    assert _result(0.5, [float("nan")]).to_dict()["null_mean"] is None


# =========================================================================== the context
def test_labels_treat_every_spelling_of_absence_as_absent(planted):
    t = planted.truth("compartment")
    assert "unassigned" not in set(t.dropna())
    assert t.isna().any()
    with pytest.raises(ValueError):
        planted.truth("no_such_column")
    assert planted.values("lineage_specific").isin([0.0, 1.0]).all()


def test_the_guard_bans_the_label_its_near_copies_and_its_edge_layer(planted):
    banned = planted.banned("cellcycle_phase")
    assert "cellcycle_phase" in banned and any(c.startswith("cycle_") for c in banned)
    assert planted.banned(None) == set()
    assert planted.banned("not_a_column") == {"not_a_column"}
    assert "compartment" in planted.banned_layers("compartment")
    assert planted.banned_layers(None) == set()
    assert "compartment" not in planted.measurement_layers("compartment")


def test_a_block_with_any_banned_column_goes_whole(planted):
    """The rule `search` applies; dropping only the banned column once leaked a near-copy."""
    assert "cycle" not in planted.blocks("cellcycle_phase")
    assert "cycle" in planted.blocks(None)
    assert planted.family_of("cycle") == "cycle"


def test_the_real_catalogue_supplies_the_family():
    ctx = S.Context.shipped("Tg")
    assert ctx.family_of("Tg_essentiality_in_a_second_background") == "fitness"
    assert "crispr_gra17ko_phenotype" in ctx.same_kind("fit_invitro_hff")
    assert ctx.same_kind("not_a_column") == set()


def test_gene_lists_are_matched_exactly_and_what_is_missing_is_named(planted):
    first = planted.gene_ids[0]
    pos, missing = planted.resolve_genes(f"{first.lower()}, nope\n{planted.gene_ids[3]}")
    assert list(pos) == [0, 3] and missing == ["nope"]
    pos, missing = planted.resolve_genes([first, first])
    assert list(pos) == [0]


def test_a_sample_keeps_what_it_must_and_favours_the_labelled(planted):
    rows = planted.sample_rows(100, target="compartment", must=[5, 7])
    assert 5 in rows and 7 in rows and len(rows) == 100 and (np.diff(rows) > 0).all()
    assert planted.truth("compartment").iloc[rows].notna().mean() > 0.6
    assert len(planted.sample_rows(None)) == planted.n
    assert len(planted.sample_rows(50, target="not_a_column")) == 50


def test_a_graph_over_a_different_table_is_refused(tmp_path, monkeypatch):
    from starplast import paths
    nodes = pd.DataFrame({"gene_id": ["A", "B", "C"], "x": [1.0, 2.0, 3.0]})
    np.savez(tmp_path / "graph.npz", gene_ids=np.array(["A", "C", "B"]), xyz=np.zeros((3, 3)),
             l__a=np.array([0]), l__b=np.array([1]), l__w=np.array([1.0]))
    monkeypatch.setattr(paths, "cache_file", lambda name: str(tmp_path / name))
    said = []
    ctx = S.Context(nodes, organism="Tg", log=said.append)
    assert ctx.graph == {} and said
    np.savez(tmp_path / "graph.npz", gene_ids=np.array(["A", "B", "C"]), xyz=np.zeros((3, 3)),
             l__a=np.array([0]), l__b=np.array([1]), l__w=np.array([1.0]))
    assert S.Context(nodes, organism="Tg").layers() == ["l"]


def test_a_bound_context_shares_the_caches_and_has_its_own_stop(planted):
    stop = {"now": False}
    bound = planted.bound(log=None, should_stop=lambda: stop["now"])
    assert bound._cache is planted._cache and bound.nodes is planted.nodes
    bound.check()
    stop["now"] = True
    with pytest.raises(Stopped):
        bound.say("anything")


def test_the_other_organism_is_reached_and_its_absence_is_tolerated(monkeypatch, planted):
    assert planted.other() is not None
    lonely = S.Context(planted.nodes, graph={}, other=planted.nodes.head(5))
    assert lonely.other().n == 5
    monkeypatch.setattr(S.Context, "shipped", classmethod(lambda cls, *a, **k: 1 / 0))
    assert S.Context(planted.nodes, graph={}).other() is None


def test_the_organism_is_read_from_the_accessions():
    assert S._guess_organism(["PF3D7_0100100", "PF3D7_0100200"]) == "Pf"
    assert S._guess_organism(["TGME49_200010"]) == "Tg"
    assert S._guess_organism([]) == "Tg"


def test_columns_are_offered_by_what_they_can_be(planted):
    cats = planted.categorical_columns()
    assert "compartment" in cats and "product" not in cats and "gene_id" not in cats
    nums = planted.numeric_columns()
    assert "fit_invitro_hff" in nums and "n_papers_focal" not in nums
    assert "fit_invitro_hff" in planted.numeric_targets()
    assert S.default_category(planted) == "compartment"
    assert S.default_numeric(planted) == "fit_invitro_hff"
    cat, members = S.example_set(planted, "compartment")
    assert cat and len(members) >= 30
    assert S.example_set(planted, None)[0] is None
    assert S.example_set(planted, "compartment", lo=10 ** 6)[0] is None


def test_a_label_free_table_still_offers_defaults():
    ctx = S.Context(pd.DataFrame({"gene_id": list("abc"), "x": [1.0, 2.0, 3.0]}), graph={})
    assert S.default_category(ctx) is None and S.default_numeric(ctx) is None


# =========================================================================== the self-tests
@pytest.mark.parametrize("key", KEYS)
def test_every_strategy_finds_what_was_planted(key, planted, planted_big):
    ctx = _ctx(key, planted, planted_big)
    s = S.get(key)
    result = s.test(ctx, **{**_genes_for(ctx, s), **BIG.get(key, {})})
    assert result.verdict == "PASS", result.summary()


@pytest.mark.parametrize("key", KEYS)
def test_no_strategy_passes_on_tables_with_nothing_in_them(key, nulls):
    s = S.get(key)
    passes = 0
    for small, big in nulls:
        ctx = _ctx(key, small, big)
        passes += s.test(ctx, **{**_genes_for(ctx, s), **BIG.get(key, {})}).passed
    assert passes <= 1, f"{key} passed on {passes} of 3 tables built from noise"


@pytest.mark.parametrize("key", KEYS)
def test_every_strategy_runs_and_says_what_it_found(key, planted, planted_big):
    ctx = _ctx(key, planted, planted_big)
    s = S.get(key)
    result = s.run(ctx, **{**_genes_for(ctx, s), **BIG.get(key, {})})
    assert result.strategy == key and result.summary
    assert result.ok, f"{key} produced no table on the planted organism"
    assert result.seconds >= 0 and result.settings
    if result.labels is not None:
        assert len(result.labels) == ctx.n
    if result.coords is not None:
        assert len(result.coords) == len(result.positions)


def test_gene_list_strategies_refuse_a_list_too_short_to_mean_anything(planted):
    for key in ("geneset_hunt", "positive_unlabeled", "set_enrichment", "seed_expansion"):
        r = S.get(key).run(planted, genes="nope")
        assert not r.ok and "gene" in r.summary.lower()


def test_a_long_user_list_is_tested_on_itself(planted):
    cat, members = S.example_set(planted, "compartment")
    r = S.get("positive_unlabeled").test(planted, genes=list(planted.gene_ids[members]),
                                         exclude="compartment")
    assert "your list" in r.hidden


def test_a_layer_built_from_the_label_is_refused(planted):
    with pytest.raises(ValueError, match="built from"):
        S.get("layer_propagation").run(planted, layer="compartment")
    with pytest.raises(ValueError, match="no edge layer"):
        S.get("layer_propagation").run(planted, layer="nope")


def test_a_missing_column_is_named(planted):
    with pytest.raises(ValueError, match="choose a column"):
        S.get("feature_knn").run(planted, target="nope")
