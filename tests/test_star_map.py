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
               # A spin box's or an editable combo's own text field is part of that control, and the
               # control carries the tip -- both are checked here in their own right.
               and not isinstance(w.parentWidget(), (QtWidgets.QAbstractSpinBox,
                                                     QtWidgets.QComboBox))]
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


# --------------------------------------------------------------------------- what counts as an edge
def _debounce_ms() -> int:
    from starplast.star_map import HOVER_DEBOUNCE_MS
    return HOVER_DEBOUNCE_MS


def _links(n, pairs, group, strengths=None):
    """A small edge frame: `pairs` of positions, all from one source."""
    a = np.array([p[0] for p in pairs], dtype=np.int32)
    b = np.array([p[1] for p in pairs], dtype=np.int32)
    s = np.asarray(strengths if strengths is not None else np.linspace(0.1, 1.0, len(pairs)),
                   dtype=np.float32)
    return pd.DataFrame({"a": a, "b": b, "kind": E.MEASURED, "group": group, "source": group,
                         "how": "layer", "setting": "", "score": s, "strength": s, "run": "",
                         "origin": E.ORIGIN_DATA, "created": "", "note": ""})[list(E.COLUMNS)]


def test_a_definition_is_data_and_round_trips():
    d = E.Definition(name="two sources", groups=["xlms", "coexpression"], min_strength=0.4,
                     max_per_gene=5, min_sources=2)
    assert d == E.Definition.from_json(d.to_json())
    assert d.signature() == d.replace(name="other").signature()   # the name is not a rule
    assert d != d.replace(min_sources=1)
    assert "2+ independent sources" in d.describe()
    assert E.Definition(min_strength=5.0).min_strength == 1.0     # clamped, never nonsense


def test_each_rule_removes_what_it_says_and_says_how_many(ctx):
    edges = pd.concat([_links(N, [(0, 1), (0, 2), (0, 3), (1, 2)], "xlms", [0.9, 0.2, 0.8, 0.5]),
                       _links(N, [(0, 1), (4, 5)], "coexpression", [0.7, 0.6])],
                      ignore_index=True)
    kept, c = E.apply_definition(edges, E.Definition(), N)
    assert c["links"] == 6 and c["sources"] == 2 and c["dropped_strength"] == 0
    kept, c = E.apply_definition(edges, E.Definition(groups=["xlms"]), N)
    assert set(kept["group"]) == {"xlms"} and c["dropped_source"] == 2
    kept, c = E.apply_definition(edges, E.Definition(min_strength=0.6), N)
    assert kept["strength"].min() >= 0.6 and c["dropped_strength"] == 2
    # agreement: only 0-1 is stated by two different sources
    kept, c = E.apply_definition(edges, E.Definition(min_sources=2), N)
    assert set(map(tuple, kept[["a", "b"]].to_numpy())) == {(0, 1)}
    assert len(kept) == 2 and c["dropped_agreement"] == 4
    # the per-gene cap keeps the strongest, and keeps the count it dropped
    kept, c = E.apply_definition(edges, E.Definition(max_per_gene=1), N)
    assert c["dropped_cap"] == len(edges) - len(kept)
    for gene in set(kept["a"]).union(kept["b"]):
        assert ((kept["a"] == gene) | (kept["b"] == gene)).sum() <= 1


def test_counts_text_is_the_live_headline():
    assert E.counts_text({"links": 4812, "genes": 3104, "sources": 3}) == \
        "4,812 links between 3,104 genes from 3 sources"
    assert "1 source" in E.counts_text({"links": 1, "genes": 2, "sources": 1})


def test_a_definition_is_saved_and_reloaded_by_name(tmp_path):
    store = E.DefinitionStore(str(tmp_path))
    d = E.Definition(groups=["xlms", "coexpression"], min_sources=2)
    name = store.save(d, "crosslink + co-expression, 2+ sources")
    assert store.names() == [name]
    again = E.DefinitionStore(str(tmp_path))              # a new session
    back = again.get(name)
    assert back == d and back.name == name
    assert again.remove(name) and again.names() == []
    assert not again.remove(name)
    with pytest.raises(ValueError):
        store.save(d, "  ")


# --------------------------------------------------------------------------- the larger network
def test_a_neighbourhood_grows_by_hops_and_stops_at_its_cap(ctx):
    chain = _links(N, [(i, i + 1) for i in range(20)], "xlms")
    one = E.neighbourhood(chain, N, [0], hops=1)
    assert set(one["gene"]) == {0, 1}
    three = E.neighbourhood(chain, N, [0], hops=3)
    assert set(three["gene"]) == {0, 1, 2, 3}
    assert dict(zip(three["gene"], three["hop"]))[3] == 3
    capped = E.neighbourhood(chain, N, [0], hops=E.MAX_HOPS, max_nodes=3)
    assert len(capped) == 3


def test_the_overview_keeps_the_strongest_and_says_what_it_left_out():
    rng = np.random.default_rng(5)
    n = 400
    pairs = list(zip(rng.integers(0, n, 3000), rng.integers(0, n, 3000)))
    edges = _links(n, pairs, "xlms", rng.random(3000))
    nodes, kept, note = E.overview(edges, n, max_nodes=50, max_edges=500)
    assert len(kept) <= 500 and len(nodes) <= 50
    assert "500" in note and "50" in note and "left out" in note
    assert set(kept["a"]).union(kept["b"]) <= set(nodes["gene"])
    # nothing cut, nothing claimed
    small = _links(n, pairs[:10], "xlms")
    assert E.overview(small, n)[2] == ""


def test_clusters_and_layout_are_deterministic_and_bounded():
    n = 40
    # four blocks of ten, every gene joined to every other gene of its own block and to nothing else
    pairs = [(i, j) for block in range(4) for i in range(block * 10, block * 10 + 10)
             for j in range(i + 1, block * 10 + 10)]
    edges = _links(n, pairs, "xlms")
    genes = np.arange(n)
    c1, c2 = E.clusters(genes, edges), E.clusters(genes, edges)
    assert np.array_equal(c1, c2)
    assert len(set(c1.tolist())) == 4                     # the four blocks, found as four clusters
    for block in range(4):
        assert len(set(c1[block * 10:block * 10 + 10].tolist())) == 1
    p1 = E.layout(genes, edges, c1)
    p2 = E.layout(genes, edges, c1)
    assert np.array_equal(p1, p2) and p1.shape == (n, 2)
    assert np.abs(p1).max() < 4000


# --------------------------------------------------------------------------- the view stands still
def _hoverable(panel):
    from starplast.star_map import _Edge, _Node
    return [i for i in panel.scene.items() if isinstance(i, (_Node, _Edge))]


def _geometry(panel):
    """Everything about the picture's geometry, to the bit: the transform, what the view shows,
    every item's position, and the scene's own bounds."""
    from starplast.star_map import _Node
    t = panel.view.transform()
    return (tuple(round(v, 9) for v in (t.m11(), t.m12(), t.m21(), t.m22(), t.dx(), t.dy())),
            panel.view.mapToScene(panel.view.viewport().rect()).boundingRect().getRect(),
            round(panel.scene.itemsBoundingRect().getRect()[2], 9),
            tuple(sorted((i.gene, round(i.pos().x(), 9), round(i.pos().y(), 9))
                         for i in panel.scene.items() if isinstance(i, _Node))))


def test_hovering_every_visible_item_never_moves_the_graph(panel):
    """The twitch the user reported: hovering must change colour and text, and nothing else."""
    from PyQt6 import QtWidgets
    from PyQt6.QtTest import QTest
    panel.centre_on(0)
    QtWidgets.QApplication.processEvents()
    items = _hoverable(panel)
    assert len(items) > 3, "the centre gene has links to hover"
    before = _geometry(panel)
    for item in items:
        QTest.mouseMove(panel.view.viewport(), panel.view.mapFromScene(item.scenePos()))
        for _ in range(3):
            QtWidgets.QApplication.processEvents()
        QTest.qWait(_debounce_ms() + 20)
        assert _geometry(panel) == before, f"hovering {item} moved the graph"
    assert _geometry(panel) == before


def test_the_hover_box_keeps_its_height_whatever_it_says(panel):
    from PyQt6 import QtWidgets
    panel.centre_on(0)
    height = panel.info.height()
    panel.info.setText("x")
    QtWidgets.QApplication.processEvents()
    short = panel.info.height()
    panel.info.setText("<br>".join(["a very long line of provenance text"] * 12))
    QtWidgets.QApplication.processEvents()
    assert panel.info.height() == short == height
    assert panel.info.maximumHeight() == panel.info.minimumHeight()


def test_the_layout_is_computed_once_and_then_held(panel, monkeypatch):
    calls = []
    real = panel.layout_positions
    monkeypatch.setattr(panel, "layout_positions",
                        lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    panel.centre_on(0)
    assert len(calls) == 1
    panel.redraw()
    panel.redraw()
    assert len(calls) == 1, "the same picture must reuse the frozen positions"
    panel.per_source.setValue(panel.per_source.value() + 1)
    assert len(calls) == 2, "a setting that changes the picture does recompute it"


# --------------------------------------------------------------------------- the larger views
def test_the_panel_shows_a_neighbourhood_and_the_whole_network(panel):
    from PyQt6 import QtWidgets
    from starplast.star_map import MODE_NEIGHBOURHOOD, MODE_NETWORK, MODE_STAR
    panel.centre_on(0)
    panel.set_mode(MODE_NEIGHBOURHOOD)
    QtWidgets.QApplication.processEvents()
    assert panel.stack.currentWidget() is panel.canvas
    assert len(panel.last_nodes) > 1 and 0 in set(panel.last_nodes["gene"])
    hops1 = len(panel.last_nodes)
    panel.hops.setValue(3)
    assert len(panel.last_nodes) >= hops1
    panel.set_mode(MODE_NETWORK)
    QtWidgets.QApplication.processEvents()
    assert len(panel.last_nodes) >= hops1
    assert panel.canvas._pos.shape == (len(panel.last_nodes), 2)
    assert "whole network" in panel.headline.text()
    panel.set_mode(MODE_STAR)
    assert panel.stack.currentWidget() is panel.view


def test_the_large_view_caps_what_it_draws_and_says_so(panel):
    from starplast.star_map import MODE_NETWORK
    panel.centre_on(0)
    panel.set_mode(MODE_NETWORK)
    panel.size_cap.setValue(50)
    assert len(panel.last_nodes) <= 50
    assert "left out" in panel.note_label.text() or len(panel.last_nodes) < 50


def test_hovering_the_large_view_moves_nothing(panel):
    from PyQt6 import QtCore, QtWidgets
    from PyQt6.QtTest import QTest
    from starplast.star_map import MODE_NETWORK
    panel.centre_on(0)
    panel.set_mode(MODE_NETWORK)
    QtWidgets.QApplication.processEvents()
    canvas = panel.canvas
    before = (canvas._scale, canvas._eye.x(), canvas._eye.y(), canvas._pos.tobytes())
    for frac in (0.25, 0.5, 0.75):
        QTest.mouseMove(canvas, QtCore.QPoint(int(canvas.width() * frac),
                                              int(canvas.height() * frac)))
        for _ in range(3):
            QtWidgets.QApplication.processEvents()
        QTest.qWait(_debounce_ms() + 20)
        assert (canvas._scale, canvas._eye.x(), canvas._eye.y(), canvas._pos.tobytes()) == before


def test_clicking_a_gene_in_the_large_view_recentres(panel):
    from PyQt6 import QtWidgets
    from starplast.star_map import MODE_NEIGHBOURHOOD
    panel.centre_on(0)
    panel.set_mode(MODE_NEIGHBOURHOOD)
    QtWidgets.QApplication.processEvents()
    other = next(g for g in panel.last_nodes["gene"] if int(g) != 0)
    panel.canvas.picked.emit(int(other))
    assert panel.centre == int(other)
    assert panel.back() and panel.centre == 0


# --------------------------------------------------------------------------- defining connections
def test_the_connection_rules_change_the_network_and_the_counts(panel):
    panel.centre_on(0)
    loose = panel.last_counts["links"]
    panel.min_strength.setValue(90)
    panel._apply_rules()
    assert panel.last_counts["links"] < loose
    assert panel.last_counts["dropped_strength"] > 0
    assert E.counts_text(panel.last_counts) in panel.counts_label.text()
    panel.min_strength.setValue(0)
    panel.min_sources.setValue(2)
    panel._apply_rules()
    assert panel.last_counts["dropped_agreement"] > 0
    panel.min_sources.setValue(1)
    panel.max_per_gene.setValue(1)
    panel._apply_rules()
    assert panel.last_counts["dropped_cap"] > 0


def test_a_definition_is_saved_and_comes_back_by_name(panel, tmp_path):
    panel.definitions = E.DefinitionStore(str(tmp_path))
    panel.set_enabled({"xlms"})
    panel.min_strength.setValue(30)
    panel.min_sources.setValue(2)
    panel.definition_box.setCurrentText("crosslink, 2+ sources")
    assert panel.save_definition() == "crosslink, 2+ sources"
    assert panel.definitions.names() == ["crosslink, 2+ sources"]
    panel.set_definition(E.Definition())
    assert panel.current_definition().min_sources == 1
    panel.set_definition(panel.definitions.get("crosslink, 2+ sources"))
    now = panel.current_definition()
    assert now.groups == ("xlms",) and now.min_sources == 2 and now.min_strength == 0.30
    assert panel.forget_definition() and panel.definitions.names() == []
