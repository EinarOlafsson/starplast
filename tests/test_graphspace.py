"""The integrated neighbour space: it finds planted structure, finds nothing in noise, cannot cheat.

Five things are checked here, and each of them is a way this module could have been wrong:

* **it works** -- on the planted organism the space recovers the planted edges and the neighbourhoods
  it lists share the planted compartment far more often than chance;
* **it does not work on nothing** -- on `planted_context(null=True)`, where every label and every edge
  is dealt out at random, the AUROC sits at chance and neither strategy passes. A self-test that
  passes on noise is a self-test of nothing;
* **the nulls differ** -- the degree-matched null is harder than the random one, on a graph built to
  have hubs and on the shipped table (the slow test at the bottom, whose measured numbers are in
  `docs/graphspace.md`);
* **self-exclusion holds** -- a layer is never a feature for predicting itself, neither directly nor
  through a shared partner reached across it;
* **the closure holds** -- holding out a label removes its column closure AND its layer closure from
  the sources, through `Context.banned` and `Context.banned_layers` rather than a copy of them.

Measured on the shipped Toxoplasma table (`docs/graphspace.md` carries the full table): the baseline
reaches AUROC 0.670 against degree-matched non-pairs, 0.678 against uniformly random ones and 0.653
against a configuration model -- a fame gap of +0.008, against +0.024 for the same features trained
the usual way with degree. The learned embedding adds +0.0075 on the degree-matched null, which does
not clear the bar it is held to (0.01, and twice its spread across layers), so the baseline ships.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import graphspace as G  # noqa: E402
from starplast import strategies as S  # noqa: E402


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture(scope="module")
def built(planted):
    """One space, built once: every test that reads a finished space shares this one."""
    return G.space(planted)


@pytest.fixture(scope="module")
def noise():
    return [S.planted_context(null=True, seed=s) for s in (7, 8, 9)]


def _hubbed(n=400, seed=3, hub_share=0.05):
    """A table whose edges follow FAME and not biology: a few hubs, edges by preferential attachment.

    The random null and the degree-matched null can only differ on a graph with a degree spread, and
    the planted organism has almost none (its layers are built with a fixed number of partners per
    gene). This builds the case the degree-matched null exists for, so the difference can be asserted
    without waiting for the shipped table.
    """
    rng = np.random.default_rng(seed)
    nodes = pd.DataFrame({"gene_id": [f"H{i:04d}" for i in range(n)],
                          **{f"m{j}": rng.normal(size=n) for j in range(8)}})
    fame = rng.pareto(1.2, size=n) + 0.05
    fame = fame / fame.sum()
    graph = {}
    for layer in ("xlms", "coexpression", "struct"):
        a = rng.choice(n, size=6 * n, p=fame)
        b = rng.choice(n, size=6 * n, p=fame)
        keep = a != b
        graph[f"{layer}__a"], graph[f"{layer}__b"] = a[keep], b[keep]
        graph[f"{layer}__w"] = rng.uniform(0.5, 1.0, size=int(keep.sum()))
    return S.Context(nodes, graph=graph, other=pd.DataFrame({"gene_id": ["X"]}), organism="Tg")


# =========================================================================== it works
def test_the_space_recovers_the_planted_structure(planted, built):
    """Hidden edges outrank degree-matched non-pairs, and the verdict is a PASS."""
    result = G.hidden_edge_test(planted, "space", layer="xlms")
    assert result.verdict == "PASS", result.summary()
    assert result.observed > 0.7
    assert built.report.auroc("logistic", "degree_matched") > 0.65
    assert 0 < len(built.pairs) <= G.MAX_PAIRS


def test_the_neighbourhoods_it_lists_share_the_planted_compartment(planted, built):
    """The planted organism puts co-expression and crosslinks inside compartments, so a neighbour
    should share one much more often than two genes picked at random do."""
    comp = planted.truth("compartment")
    top = built.to_frame(top=400)
    same = [comp.iloc[built.position(a)] == comp.iloc[built.position(b)]
            for a, b in zip(top["gene_a"], top["gene_b"])]
    share = float(np.mean([s for s in same if s == s]))
    labelled = comp.dropna()
    chance = float((labelled.value_counts(normalize=True) ** 2).sum())
    assert share > 3 * chance, f"{share:.3f} against a chance level of {chance:.3f}"


def test_a_gene_can_be_asked_for_its_neighbours_and_the_evidence_behind_each(planted, built):
    gene = planted.gene_ids[built.pairs[0, 0]]
    table = built.neighbours(gene, 5)
    assert 0 < len(table) <= 5
    assert list(table["probability"]) == sorted(table["probability"], reverse=True)
    assert table["measured"].dtype == bool and "measured_in" in table
    a, b = built.pairs[0]
    ev = built.evidence(a, b)
    assert set(ev["source"]) == set(built.sources)
    assert ev["measured_edge"].any() or not built.measured[0]
    # The contributions and the intercept ARE the logit of the probability: the attribution is the
    # model, not a second story told about it.
    A, B = built.platt
    logit = ev["contribution"].sum() + A * built.intercept + B
    assert 1 / (1 + np.exp(-logit)) == pytest.approx(built.edge_strength(a, b), abs=1e-9)


def test_a_pair_that_was_never_a_candidate_is_still_scored_on_demand(planted, built):
    """The candidate cap decides what is RANKED, not what can be asked about."""
    keys = {(int(a), int(b)) for a, b in built.pairs}
    outside = next((i, j) for i in range(planted.n) for j in range(i + 1, planted.n)
                   if (i, j) not in keys)
    p = built.edge_strength(*outside)
    assert 0.0 <= p <= 1.0 and built.row_of(*outside) == -1
    assert len(built.evidence(*outside)) == len(built.sources)


def test_inferred_edges_are_distinguishable_from_measured_ones_everywhere(built):
    frame = built.to_frame()
    assert set(frame["kind"]) <= {"measured", "inferred"}
    assert (frame.loc[frame["kind"] == "inferred", "measured_in"] == "").all()
    assert (frame.loc[frame["kind"] == "measured", "measured_in"] != "").all()
    gaps = built.gaps(50)
    assert len(gaps) and (gaps["kind"] == "inferred").all()
    assert "probability" in gaps and "top_evidence" in gaps


def test_the_measured_layers_are_never_touched(planted):
    """An integrated edge is a claim. The shipped layers must be byte-for-byte what they were."""
    before = {k: np.array(v, copy=True) for k, v in planted.graph.items()}
    sp = G.space(planted)
    assert len(sp.pairs)
    after = planted.graph
    assert set(after) == set(before)
    for k, v in before.items():
        assert np.array_equal(np.asarray(after[k]), v), k


# =========================================================================== and not on noise
def test_nothing_is_found_in_a_table_with_nothing_in_it(noise):
    for ctx in noise:
        report = G.evaluate(ctx)
        auroc = report.auroc("logistic", "degree_matched")
        assert abs(auroc - 0.5) < 0.12, f"{auroc:.3f} on a table built from noise"


def test_neither_strategy_passes_on_a_table_built_from_noise(noise):
    """The house rule for every strategy, asserted here with no allowance: both must FAIL, on all
    three tables, because a space built from randomly rewired layers has nothing to recover."""
    for ctx in noise:
        for key in ("neighbour_space", "network_training"):
            result = S.get(key).test(ctx)
            assert not result.passed, f"{key} passed on noise: {result.summary()}"


def test_both_strategies_find_the_planted_structure(planted):
    for key in ("neighbour_space", "network_training"):
        result = S.get(key).test(planted)
        assert result.verdict == "PASS", f"{key}: {result.summary()}"
        assert result.numbers["fame_gap"] == result.numbers["fame_gap"]


# =========================================================================== the nulls differ
def test_the_easy_null_flatters_the_usual_recipe_and_the_hard_one_does_not():
    """The whole claim of the module, on a graph whose edges ARE fame: same data, same held-out
    edges, two negative sets, and a gap between them that is the share of the number which is fame.

    Measured on the preferential-attachment table below: the usual recipe (degree as a feature,
    random non-pairs to train against) reaches AUROC 0.92 against random non-pairs and 0.71 against
    degree-matched ones -- a gap of +0.22, and a Brier score that collapses from 0.11 to 0.32 because
    its probabilities are calibrated for a world where every non-pair is two obscure genes. The
    shipped baseline, given no degree and trained against degree-matched non-pairs, scores 0.75 on
    the hard null and BEATS the fame model there while losing to it by 0.38 on the easy one. Its own
    gap is negative, which is what "not exploiting degree" looks like as a number.
    """
    ctx = _hubbed()
    report = G.evaluate(ctx)
    assert report.scored, "the hubbed table produced no scorable fold"
    assert report.fame_recipe_gap > 0.05, report.summary()
    assert report.fame_recipe_gap > report.fame_gap + 0.1, report.summary()
    assert report.auroc("fame", "random") > report.auroc("logistic", "random") + 0.1
    assert report.auroc("fame", "degree_matched") < report.auroc("logistic", "degree_matched")
    assert report.metric("brier", "fame", "degree_matched") > \
        report.metric("brier", "logistic", "degree_matched")
    # And the gate on the learned model is a real gate, not a way of never shipping it: on this graph
    # the embedding gains +0.049 on the degree-matched null, in all three folds, and is offered.
    assert report.learned_earns_its_place, report.summary()
    assert all(g > 0 for g in report.learned_gain_per_layer())


def test_the_degree_matched_negatives_really_are_degree_matched():
    ctx = _hubbed()
    ev = G.Evidence(ctx)
    layer = ev.targets()[0]
    deg = ev.degree(layer)
    pos = ev.edges(layer)
    nodes = np.flatnonzero(deg > 0)
    rng = ctx.rng(5)
    keys = set(G._pack(ev.all_edges()[:, 0], ev.all_edges()[:, 1], ev.n).tolist())
    matched = G.degree_matched_negatives(rng, deg, pos, keys, nodes, ev.n)
    plain = G.random_negatives(rng, np.arange(ctx.n), keys, len(pos), ev.n)
    assert len(matched) > 0.5 * len(pos)
    d = lambda p: np.log1p(deg[p]).mean()
    assert abs(d(matched) - d(pos)) < abs(d(plain) - d(pos))
    # And no negative may be an edge of any permitted layer.
    assert not set(G._pack(matched[:, 0], matched[:, 1], ev.n).tolist()) & keys


def test_the_configuration_null_keeps_the_degrees_and_loses_the_topology(planted):
    ev = G.Evidence(planted)
    pos = ev.edges("xlms")
    rng = planted.rng(1)
    rewired = G.configuration_negatives(rng, pos, set(), ev.n)
    before = np.bincount(pos.ravel(), minlength=ev.n)
    after = np.bincount(rewired.ravel(), minlength=ev.n)
    # Stub pairing loses the self-loops and duplicates it draws, so degrees are preserved on average
    # rather than exactly; what must hold is that the sequence is the positives' and not uniform.
    assert np.corrcoef(before, after)[0, 1] > 0.5
    shared = set(G._pack(pos[:, 0], pos[:, 1], ev.n).tolist()) & set(
        G._pack(rewired[:, 0], rewired[:, 1], ev.n).tolist())
    assert len(shared) < 0.2 * len(pos)


def test_every_number_is_reported_against_all_three_nulls(built):
    frame = built.report.frame()
    for model in G.MODELS:
        got = set(frame.loc[frame["model"] == model, "null"])
        assert got == set(G.NULLS), f"{model}: {got}"
    assert {"auroc", "precision_at_k", "brier", "reliability_gap"} <= set(frame.columns)
    gap, table = G.reliability(built.probability, built.measured)
    assert 0 <= gap <= 1 and len(table) >= 2


def test_the_learned_model_is_only_offered_when_it_beats_the_baseline_where_it_is_hard(built):
    """The gate, asserted as arithmetic: a gain under 0.01 AUROC, or inside twice its own spread
    across folds, is not a gain, and the interpretable model keeps the place."""
    report = built.report
    gains = report.learned_gain_per_layer()
    assert len(gains) == len(report.scored)
    if report.learned_earns_its_place:
        se = np.std(gains, ddof=1) / np.sqrt(len(gains))
        assert report.learned_gain > 0.01 and report.learned_gain > 2 * se
    # And the strategy honours the gate rather than the request.
    result = S.get("network_training").run(built.ctx, model="embedding", top=20)
    assert result.numbers["model_used"] == ("embedding" if report.learned_earns_its_place
                                            else "logistic")
    assert "embedding" in result.summary


# =========================================================================== self-exclusion
def test_a_layer_is_never_a_feature_for_predicting_itself(planted):
    ev = G.Evidence(planted)
    for layer in ev.targets():
        names = ev.names(layer)
        assert layer not in names and f"{layer}_absent" not in names
        assert len(names) == len(ev.names(None)) - 2
    report = G.evaluate(planted)
    for fold in report.folds:
        assert fold.layer not in fold.sources
        assert fold.layer not in fold.coefficients


def test_the_shared_partner_count_is_rebuilt_without_the_held_out_layer(planted):
    """Self-exclusion is not only the layer's own column: a partner shared THROUGH the held-out layer
    is that layer's evidence arriving by a second route, and it has to go too."""
    ev = G.Evidence(planted)
    layer = "xlms"
    assert ev.union(layer).nnz < ev.union(None).nnz
    pairs = ev.edges(layer)[:200]
    a, b = pairs[:, 0], pairs[:, 1]
    # The partners themselves: a subgraph can only have fewer of them, elementwise.
    full_n = np.asarray(ev.union(None)[a].multiply(ev.union(None)[b]).sum(axis=1)).ravel()
    without_n = np.asarray(ev.union(layer)[a].multiply(ev.union(layer)[b]).sum(axis=1)).ravel()
    assert (without_n <= full_n).all() and without_n.sum() < full_n.sum()
    assert int((without_n == 0).sum()) > int((full_n == 0).sum())
    # The feature is Adamic-Adar, not a raw count, so it does NOT fall monotonically: removing a
    # layer also removes degree, and 1/log(degree) rises for the partners that remain. Asserting
    # monotonicity on the feature is the mistake this comment exists to stop being made again.
    full = ev.features(pairs)[:, list(ev.names(None)).index("shared_partners")]
    without = ev.features(pairs, skip=layer)[:, list(ev.names(layer)).index("shared_partners")]
    assert (without[full_n == 0] == 0).all()
    assert not np.allclose(full, without)


def test_a_source_that_could_only_predict_itself_ends_at_weight_zero(built):
    """The integrated weight is the mean over folds counting the excluded fold as zero, so the rule
    is arithmetic rather than intention. Here: every weight is the mean the report says it is."""
    weights = built.report.weights().set_index("source")["weight"].to_dict()
    for source, weight in weights.items():
        assert built.weights[source] == pytest.approx(weight)
    for fold in built.report.scored:
        assert fold.layer in built.weights and f"{fold.layer}_absent" in built.weights


# =========================================================================== the leakage guard
def test_holding_out_a_label_removes_its_column_closure_and_its_layer_closure(planted):
    """Both halves of the guard, through the project's own closure and not a copy of it.

    Two labels, because on this table each demonstrates one half. `compartment` has an edge layer
    built from it and no numeric near-copies, so it tests the LAYER closure; `cellcycle_phase` has no
    layer of its own but its closure removes the six cycle columns, so it tests the COLUMN closure.
    A test that used only the first would pass with the column closure never applied at all.
    """
    open_space = G.Evidence(planted)
    layer_banned = planted.banned_layers("compartment")
    assert "compartment" in layer_banned
    assert "compartment" in open_space.layers, "the layer exists when nothing is held out"
    guarded = G.Evidence(planted, "compartment")
    assert "compartment" not in guarded.layers
    for name in guarded.names(None):
        assert G._group(name) not in layer_banned
    column_banned = planted.banned("cellcycle_phase")
    assert {"cycle_t0", "cycle_t2"} <= column_banned, "the project's closure caught the near-copies"
    phase_guarded = G.Evidence(planted, "cellcycle_phase")
    assert not set(phase_guarded.columns) & column_banned
    assert set(open_space.columns) - set(phase_guarded.columns) == column_banned & set(
        open_space.columns)
    # And the similarity feature really is computed from the guarded matrix, not from ctx.nodes.
    assert phase_guarded.X.shape[1] == len(phase_guarded.columns) < open_space.X.shape[1]


def test_a_space_built_while_holding_out_a_label_carries_the_closure_into_every_source(planted):
    sp = G.space(planted, exclude="compartment")
    assert "compartment" not in sp.layers
    assert not set(sp.columns) & planted.banned("compartment")
    assert all(G._group(s) != "compartment" for s in sp.sources)
    assert sp.report.settings["exclude"] == "compartment"
    # And the strategies pass the same setting through rather than building an unguarded space.
    result = S.get("neighbour_space").run(planted, exclude="compartment", top=20, k=3)
    assert result.ok and "compartment" not in result.tables["how each source is weighted"]["source"]\
        .map(G._group).tolist()


def test_the_literature_layers_are_never_sources(planted):
    """Co-mention follows attention, and this project keeps it out of anything that claims to say
    what the MEASUREMENTS imply."""
    ev = G.Evidence(planted)
    assert "comention" not in ev.layers and "comention_ft" not in ev.layers
    for derived in S.DERIVED_LAYERS:
        assert derived not in ev.layers


# =========================================================================== edges and refusals
def test_a_table_with_no_layers_says_so_rather_than_returning_an_empty_space():
    ctx = S.Context(pd.DataFrame({"gene_id": ["A", "B", "C"], "m": [1.0, 2.0, 3.0]}), graph={})
    with pytest.raises(ValueError, match="no permitted edge layer"):
        G.space(ctx)


def test_a_layer_that_cannot_be_held_out_is_refused_by_name(planted):
    with pytest.raises(ValueError, match="not one that can be held out"):
        G.hidden_edge_test(planted, "k", layer="compartment")
    with pytest.raises(ValueError, match="not a layer this table permits"):
        G.hidden_edge_test(planted, "k", layer="no_such_layer")


def test_too_few_edges_to_hide_is_inconclusive_rather_than_a_pass(planted):
    result = G.hidden_edge_test(planted, "k", layer="xlms", fraction=0.01)
    assert result.verdict in ("INCONCLUSIVE", "FAIL")
    if result.verdict == "INCONCLUSIVE":
        assert "orthogroups" in result.note or "few" in result.note


def test_metrics_and_reliability_survive_a_degenerate_input():
    m = G.metrics([0.1, 0.2], [True, True])
    assert np.isnan(m["auroc"]) and m["positives"] == 2
    gap, table = G.reliability([], [])
    assert np.isnan(gap) and table.empty


def test_the_evaluation_is_repeatable(planted):
    a = G.evaluate(planted).auroc("logistic", "degree_matched")
    b = G.evaluate(planted).auroc("logistic", "degree_matched")
    assert a == b


def test_the_two_strategies_are_registered_with_the_conventions_the_panel_needs():
    catalog = {s.key: s for s in S.catalog()}
    for key, number in (("neighbour_space", 33), ("network_training", 34)):
        s = catalog[key]
        assert s.number == number and s.family in S.families()
        assert len(s.tooltip.split()) >= 25 and s.question.endswith("?")
        assert len(s.walkthrough) >= 4 and len(s.explanation.split()) >= 120
        assert "null" in s.test_description.lower() and "pass" in s.test_description.lower()
        assert "degree-matched" in s.test_description.lower()
        for p in s.params:
            assert p.kind in S.PARAM_KINDS and len(p.tip.split()) >= 15, f"{key}.{p.name}"


# =========================================================================== the shipped table
@pytest.mark.slow
def test_on_the_shipped_table_the_fame_gap_is_positive_and_the_numbers_are_the_documented_ones():
    """The real measurement, and the one the documentation quotes. Minutes, hence the marker.

    The gap is positive and small on this graph, and small is the finding: the shipped co-expression
    and co-fitness layers are built with a cap on partners per gene, so their degree distribution is
    narrow and there is less fame to remove than on a literature-derived network. The null's job is
    to establish that rather than to assume it.
    """
    ctx = S.Context.shipped("Tg")
    sp = G.space(ctx)
    report = sp.report
    assert report.fame_gap > 0, report.summary()
    assert report.fame_recipe_gap > report.fame_gap, report.summary()
    assert report.learned_gain == report.learned_gain
    assert report.auroc("logistic", "degree_matched") > 0.6
    assert len(sp.pairs) > 100000 and (~sp.measured).sum() > 10000
    assert len(sp.gaps(100)) == 100
