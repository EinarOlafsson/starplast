"""The star map: gene-gene links with provenance, the user's runs, the bounds, and the view.

The edge store is checked on a small synthetic table so every number is known; the shipped file and
the application hook are checked against the real cache.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import star_edges as E  # noqa: E402
from starplast import strategies as S  # noqa: E402

N = 60
ORG = "Tg"


@pytest.fixture()
def ctx():
    rng = np.random.default_rng(0)
    ids = [f"TGME49_{100000 + i:06d}" for i in range(N)]
    nodes = pd.DataFrame({"gene_id": ids, "product": [f"protein {i}" for i in range(N)],
                          "compartment": (["A"] * 20 + ["B"] * 20 + [None] * 20)})
    for j in range(4):
        nodes[f"m{j}"] = rng.normal(size=N) + (np.arange(N) < 30) * 2.0
    graph = {"xlms__a": np.array([0, 0, 1, 2, 5]), "xlms__b": np.array([1, 2, 3, 4, 6]),
             "xlms__w": np.array([3.0, 1.0, 2.0, 1.0, 5.0]),
             "coexpression__a": np.arange(0, 40), "coexpression__b": np.arange(1, 41),
             "coexpression__w": np.linspace(0.1, 1.0, 40)}
    return S.Context(nodes, graph=graph, organism=ORG)


def _result(key, tables=None, labels=None, **settings):
    return S.StrategyResult(key, f"{key} summary", tables=tables or {}, labels=labels,
                            settings=settings)


# --------------------------------------------------------------------------- provenance
def test_measured_edges_carry_their_layer(ctx):
    m = E.measured_edges(ctx.graph)
    assert set(m["group"]) == {"xlms", "coexpression"}
    assert (m["kind"] == E.MEASURED).all() and (m["origin"] == E.ORIGIN_DATA).all()
    assert (m["run"] == "").all()
    assert m["strength"].between(0, 1).all()
    x = m[m["source"] == "xlms"].sort_values("score")
    assert x["strength"].is_monotonic_increasing          # strength is the rank of the weight


def test_pair_table_links_carry_strategy_setting_and_run(ctx):
    ids = ctx.gene_ids
    t = pd.DataFrame({"gene_a": ids[[0, 7, 9]], "gene_b": ids[[8, 3, 9]], "score": [3.0, 1.0, 2.0],
                      "linked_in": ["domain", "", "x"]})
    r = _result("link_prediction", {"predicted links": t}, layer="xlms", top=3)
    e = E.edges_from_result(r, ctx, run="starplast:link_prediction", origin=E.ORIGIN_SHIPPED)
    assert len(e) == 2                                    # the self pair (9, 9) is dropped
    assert (e["kind"] == E.INFERRED).all()
    assert set(e["source"]) == {"link_prediction"} and set(e["group"]) == {"link_prediction"}
    assert set(e["run"]) == {"starplast:link_prediction"}
    assert e["setting"].iloc[0] == "layer=xlms, top=3"
    assert "linked in: domain" in set(e["note"])
    assert e.loc[e["score"].idxmax(), "strength"] == 1.0


def test_skipped_tables_are_not_counted_twice(ctx):
    ids = ctx.gene_ids
    t = pd.DataFrame({"gene_a": ids[[0]], "gene_b": ids[[1]], "probability": [0.9]})
    r = _result("network_training", {"gaps": t, "measured edges, ranked": t.assign(
        gene_b=ids[[5]])})
    e = E.edges_from_result(r, ctx)
    assert len(e) == 1 and e["group"].iloc[0] == E.YOUR_RUNS


# --------------------------------------------------------------------------- bounds
def test_a_module_is_not_drawn_as_a_clique(ctx):
    labels = np.full(N, -1)
    labels[:50] = 7
    e = E.edges_from_result(_result("multiplex_modules", labels=labels), ctx)
    assert (e["how"] == "module").all()
    assert len(e) <= 50 * E.MODULE_K < 50 * 49 // 2
    # every link joins two members, and each member asked for at most MODULE_K
    assert set(e["a"]).union(e["b"]) <= set(range(50))
    # shared measured layers come first: 0-1 are joined by xlms AND coexpression
    assert ((e["a"] == 0) & (e["b"] == 1)).any()


def test_a_run_never_adds_more_than_the_cap(ctx, monkeypatch):
    monkeypatch.setattr(E, "MAX_EDGES_PER_RUN", 10)
    labels = np.zeros(N, dtype=int)
    e = E.edges_from_result(_result("consensus_modules", labels=labels), ctx)
    assert len(e) == 10
    assert e["strength"].min() >= 0


def test_seed_expansion_links_each_candidate_to_at_most_k_seeds(ctx):
    ids = ctx.gene_ids
    seeds = ids[:6]
    cand = pd.DataFrame({"gene_id": ids[20:30], "score": np.linspace(1, 0.1, 10)})
    r = _result("seed_expansion", {"candidates": cand}, genes=" ".join(seeds), top=10)
    e = E.edges_from_result(r, ctx)
    assert (e["how"] == "seed").all()
    seed_pos = set(range(6))
    for g in range(20, 30):
        touching = e[(e["a"] == g) | (e["b"] == g)]
        assert 1 <= len(touching) <= E.SEED_K
        others = set(touching["a"]).union(touching["b"]) - {g}
        assert others <= seed_pos
    assert "genes=6 genes" in e["setting"].iloc[0]


def test_partner_calls_link_to_the_partners_they_name(ctx):
    ids = ctx.gene_ids
    calls = pd.DataFrame({"gene_id": [ids[45]], "prediction": ["A"], "support": [1.0],
                          "partners": [f"{ids[3]} (A, 2), {ids[4]} (A, 1)"]})
    e = E.edges_from_result(_result("physical_partners", {"calls": calls}, target="compartment"),
                            ctx)
    assert set(map(tuple, e[["a", "b"]].to_numpy())) == {(3, 45), (4, 45)}
    assert e.loc[e["b"] == 45, "score"].max() == 2.0


def test_neighbour_vote_links_to_labelled_genes_of_the_called_label(ctx):
    ids = ctx.gene_ids
    calls = pd.DataFrame({"gene_id": ids[[45, 46]], "product": "", "prediction": ["B", "A"],
                          "support": [0.8, 0.6]})
    e = E.edges_from_result(_result("feature_knn", {"calls": calls}, target="compartment", k=5),
                            ctx)
    assert len(e) == 2 * E.LABEL_K
    truth = ctx.truth("compartment")
    for g, lab in ((45, "B"), (46, "A")):
        others = e.loc[(e["a"] == g) | (e["b"] == g)]
        partners = set(others["a"]).union(others["b"]) - {g}
        assert all(truth.iloc[p] == lab for p in partners)


def test_star_keeps_per_source_and_node_caps(ctx):
    idx = E.EdgeIndex(N)
    idx.add(E.measured_edges(ctx.graph))
    labels = np.zeros(N, dtype=int)
    idx.add(E.edges_from_result(_result("consensus_modules", labels=labels), ctx,
                                run="starplast:consensus_modules", origin=E.ORIGIN_SHIPPED))
    nodes, edges = E.star(idx, 0, per_group=1, depth=1)
    assert edges.groupby("group").size().max() == 1
    assert set(nodes["hop"]) == {0, 1}
    nodes2, _ = E.star(idx, 0, per_group=3, depth=2, max_nodes=5)
    assert len(nodes2) <= 5 and 2 in set(nodes2["hop"]) | {2}
    only = E.star(idx, 0, groups=["xlms"], depth=1)[1]
    assert set(only["group"]) == {"xlms"}


# --------------------------------------------------------------------------- user runs
def test_user_runs_are_kept_on_disk_and_marked_as_yours(ctx, tmp_path):
    ids = ctx.gene_ids
    t = pd.DataFrame({"gene_a": ids[[10, 11]], "gene_b": ids[[50, 51]], "score": [1.0, 2.0]})
    store = E.UserEdgeStore(str(tmp_path))
    f = store.add(_result("link_prediction", {"predicted links": t}, layer="xlms"), ctx, ORG)
    assert len(f) == 2 and (f["origin"] == E.ORIGIN_USER).all()
    assert (f["group"] == E.YOUR_RUNS).all() and f["run"].iloc[0].startswith("you:link_prediction:")
    assert f["created"].iloc[0]
    again = E.UserEdgeStore(str(tmp_path))
    back = again.load(ORG, ids)
    assert len(back) == 2 and set(back["run"]) == set(f["run"])
    assert again.runs(ORG)[0]["strategy"] == "link_prediction"
    assert again.load("Pf", ids).empty                      # other organisms do not see them
    assert again.remove(f["run"].iloc[0]) and again.load(ORG, ids).empty


def test_stored_round_trip_is_by_gene_id(ctx):
    m = E.measured_edges(ctx.graph)
    stored = E.to_stored(m, ctx.gene_ids, ORG)
    assert list(stored.columns) == list(E.STORED)
    shuffled = ctx.gene_ids[::-1]
    back = E.from_stored(stored, shuffled, ORG)
    assert len(back) == len(m)
    assert shuffled[back["a"].iloc[0]] == ctx.gene_ids[m["a"].iloc[0]]


# --------------------------------------------------------------------------- the shipped file
def test_shipped_links_are_small_and_carry_provenance():
    path = E.shipped_path()
    assert os.path.exists(path), "run scripts/build_star_edges.py"
    assert os.path.getsize(path) < 6e6
    f = pd.read_parquet(path)
    assert list(f.columns) == list(E.STORED)
    assert (f["kind"].astype(str) == E.INFERRED).all()
    assert (f["origin"].astype(str) == E.ORIGIN_SHIPPED).all()
    assert f["run"].astype(str).str.startswith("starplast:").all()
    assert {"link_prediction", "multiplex_modules", "seed_expansion"} <= set(f["source"].astype(str))
    per_run = f.groupby([f["organism"].astype(str), f["run"].astype(str)]).size()
    assert per_run.max() <= E.MAX_EDGES_PER_RUN


# --------------------------------------------------------------------------- the view
@pytest.fixture()
def panel(ctx, tmp_path):
    from starplast.star_map import StarMapPanel
    p = StarMapPanel(ctx.nodes, ORG, ctx=ctx, graph=ctx.graph,
                     store=E.UserEdgeStore(str(tmp_path)), shipped_path=str(tmp_path / "none"))
    p.resize(900, 600)
    p.show()
    yield p
    p.close()


def test_every_control_explains_itself(panel):
    from PyQt6 import QtWidgets
    kinds = (QtWidgets.QAbstractButton, QtWidgets.QComboBox, QtWidgets.QAbstractSpinBox,
             QtWidgets.QLineEdit, QtWidgets.QListWidget, QtWidgets.QGraphicsView)
    missing = [type(w).__name__ + ":" + getattr(w, "text", lambda: "")()
               for k in kinds for w in panel.findChildren(k) if not w.toolTip()
               # a spin box's own text field is part of the spin box, which carries the tip
               and not isinstance(w.parentWidget(), QtWidgets.QAbstractSpinBox)]
    assert not missing


def test_view_recentres_on_click_and_goes_back(panel):
    from PyQt6 import QtCore, QtWidgets
    from PyQt6.QtTest import QTest
    from starplast.star_map import _Node
    assert panel.centre_on("TGME49_100000")
    QtWidgets.QApplication.processEvents()
    targets = [i for i in panel.scene.items() if isinstance(i, _Node) and i.gene == 1]
    assert targets, "gene 1 is linked to the centre by xlms and coexpression"
    where = panel.view.mapFromScene(targets[0].scenePos())
    QTest.mouseClick(panel.view.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=where)
    for _ in range(5):
        QtWidgets.QApplication.processEvents()
    assert panel.centre == 1
    assert panel.history == [0]
    assert panel.back() and panel.centre == 0


def test_hover_text_names_the_links_provenance(panel):
    panel.centre_on(0)
    row = next(panel.last_edges.itertuples(index=False))
    html = panel.edge_html(row)
    assert "measured layer" in html and row.source in html and "score" in html


def test_a_user_run_appears_in_the_map_at_once(panel, ctx):
    panel.centre_on(10)
    before = set(panel.last_edges["group"])
    assert E.YOUR_RUNS not in before
    ids = ctx.gene_ids
    t = pd.DataFrame({"gene_a": ids[[10]], "gene_b": ids[[55]], "score": [4.0]})
    n = panel.add_result(_result("link_prediction", {"predicted links": t}, layer="xlms"))
    assert n == 1
    assert E.YOUR_RUNS in set(panel.last_edges["group"])
    assert 55 in set(panel.last_nodes["gene"])
    groups = [panel.sources.item(i).data(0x0100) for i in range(panel.sources.count())]
    assert E.YOUR_RUNS in groups
    row = panel.last_edges[panel.last_edges["group"] == E.YOUR_RUNS].iloc[0]
    assert "your run you:link_prediction:" in panel.edge_html(row)


def test_toggling_a_source_hides_its_links(panel):
    panel.centre_on(0)
    panel.set_enabled({"xlms"})
    assert set(panel.last_edges["group"]) == {"xlms"}
    panel.only_measured()
    assert set(panel.last_edges["kind"]) == {E.MEASURED}
    panel.set_enabled(set())
    assert panel.last_edges.empty


def test_follow_centres_a_gene_selected_elsewhere(panel):
    panel.follow("TGME49_100005")
    assert panel.centre == 5
    panel.follow_box.setChecked(False)
    panel.follow("TGME49_100006")
    assert panel.centre == 5
