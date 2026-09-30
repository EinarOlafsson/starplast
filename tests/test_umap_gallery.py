"""The pregenerated map gallery: its scores on planted answers, the shipped files, and the dock."""
import hashlib
import os

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from starplast import organisms  # noqa: E402
from starplast import umap_gallery as G  # noqa: E402

TG = organisms.TOXOPLASMA


# --------------------------------------------------------------------------- arithmetic
def test_wilson_lower_bound_matches_the_closed_form():
    # 2 of 2: centre (1 + z^2/4) / (1 + z^2/2), margin z*sqrt(z^2/16)/(1 + z^2/2)
    z = G.Z95
    expect = (1 + z * z / 4 - z * z / 4) / (1 + z * z / 2)
    assert G.wilson_lower(2, 2) == pytest.approx(expect, abs=1e-12)
    assert G.wilson_lower(2, 2) == pytest.approx(0.3424, abs=1e-4)
    assert G.wilson_lower(0, 10) == 0.0
    assert G.wilson_lower(5, 0) == 0.0
    assert G.wilson_lower(900, 1000) == pytest.approx(0.8798, abs=1e-3)


def test_labels_that_are_the_clusters_score_one_and_far_above_chance():
    rng = np.random.default_rng(0)
    clusters = rng.integers(0, 5, 800)
    truth = np.array([f"c{k}" for k in clusters], dtype=object)
    s = G.label_scores(truth, clusters)
    assert s["category_f1"] == pytest.approx(1.0)
    assert s["category_precision"] == pytest.approx(1.0)
    assert s["category_recall"] == pytest.approx(1.0)
    assert s["category_f1_chance"] < 0.35
    assert s["category_skill"] == pytest.approx(1.0)
    assert s["best_skill"] > 0.9
    assert s["categories"] == 5 and s["genes_scored"] == 800


def test_labels_unrelated_to_the_clusters_score_at_chance():
    rng = np.random.default_rng(1)
    clusters = rng.integers(0, 6, 3000)
    truth = rng.choice(["a", "b", "c", "d"], 3000).astype(object)
    s = G.label_scores(truth, clusters)
    assert abs(s["category_skill"]) < 0.05
    assert abs(s["best_skill"]) < 0.1


def test_precision_and_recall_are_the_planted_ones():
    # 100 genes of "x": 80 in cluster 0 (with 20 of "y"), 20 in cluster 1 (with 200 of "y").
    truth = np.array(["x"] * 100 + ["y"] * 220, dtype=object)
    clusters = np.array([0] * 80 + [1] * 20 + [0] * 20 + [1] * 200)
    per = G.per_category(truth, clusters).set_index("category")
    assert per.loc["x", "precision"] == pytest.approx(0.8)
    assert per.loc["x", "recall"] == pytest.approx(0.8)
    assert per.loc["y", "precision"] == pytest.approx(200 / 220)
    assert per.loc["y", "recall"] == pytest.approx(200 / 220)
    s = G.label_scores(truth, clusters)
    w = np.array([100, 220]) / 320
    assert s["category_f1"] == pytest.approx(w[0] * 0.8 + w[1] * 200 / 220)


def test_unclustered_genes_count_as_missed_not_as_absent():
    truth = np.array(["x"] * 10 + ["y"] * 90, dtype=object)
    clusters = np.array([0] * 5 + [G.NOISE] * 5 + [1] * 90)
    per = G.per_category(truth, clusters).set_index("category")
    assert per.loc["x", "precision"] == pytest.approx(1.0)
    assert per.loc["x", "recall"] == pytest.approx(0.5), "noise is a miss, not a smaller category"
    assert per.loc["x", "n"] == 10


def test_unlabelled_and_absence_values_are_not_categories():
    truth = np.array(["x"] * 50 + ["unknown"] * 20 + [None] * 20 + ["y"] * 50, dtype=object)
    clusters = np.array([0] * 50 + [0] * 40 + [1] * 50)
    per = G.per_category(truth, clusters)
    assert set(per.category) == {"x", "y"}
    assert G.label_scores(truth, clusters)["genes_scored"] == 100


def test_a_two_gene_category_cannot_win_by_sitting_alone_in_a_cluster():
    """The user's case: a category of 2 genes that maps perfectly must not score 1."""
    rng = np.random.default_rng(3)
    n = 1200
    big = np.zeros(n, bool)
    big[:200] = True
    clusters = rng.integers(1, 8, n)
    # The big category: 90% of its genes in cluster 0, which is 90% pure.
    clusters[:180] = 0
    clusters[200:220] = 0
    truth = rng.choice(["p", "q", "r"], n).astype(object)
    truth[:200] = "big"
    # Two genes alone in a cluster of their own: raw precision 1 and recall 1.
    truth[1000:1002] = "tiny"
    clusters[1000:1002] = 99
    per = G.per_category(truth, clusters).set_index("category")
    assert per.loc["tiny", "precision"] == 1.0 and per.loc["tiny", "recall"] == 1.0
    assert per.loc["tiny", "f1"] == 1.0, "raw F1 would call it perfect"
    assert per.loc["tiny", "f1_lower"] == pytest.approx(G.wilson_lower(2, 2), abs=1e-9)
    assert per.loc["tiny", "f1_lower"] < 0.35
    assert per.loc["big", "f1_lower"] > 0.8
    s = G.label_scores(truth, clusters)
    assert s["best_category"] == "big" and s["best_n"] == 200
    assert s["best_f1_lower"] == pytest.approx(per.loc["big", "f1_lower"])
    assert s["best_precision"] == pytest.approx(0.9) and s["best_recall"] == pytest.approx(0.9)


def test_the_majority_of_a_yes_no_label_does_not_win_on_one_giant_cluster():
    """One giant cluster gives the majority value a high bound by size alone -- and by chance."""
    rng = np.random.default_rng(4)
    clusters = np.zeros(4000, int)
    clusters[:40] = 1
    truth = np.where(rng.random(4000) < 0.9, "False", "True").astype(object)
    s = G.label_scores(truth, clusters)
    assert s["best_f1_lower"] > 0.9, "the raw bound is high"
    assert abs(s["best_skill"]) < 0.05, "and worthless: shuffled labels score the same"
    assert abs(s["category_skill"]) < 0.05


def test_skill_is_zero_at_chance_one_at_perfect():
    assert G.skill(0.2, 0.2) == 0.0
    assert G.skill(1.0, 0.3) == 1.0
    assert np.isnan(G.skill(np.nan, 0.3))
    assert np.isnan(G.skill(0.5, 1.0))


def test_score_map_on_a_planted_table_finds_the_planted_label():
    from starplast import strategies as S
    from starplast.clustering import cluster
    ctx = S.planted_context(n=480, seed=7)
    X, cols = ctx.features()
    comp = ctx.truth("compartment")
    # A "map" whose clusters are the compartments, minus some noise: compartment must lead.
    rows = np.arange(ctx.n)
    codes = pd.factorize(comp)[0]
    table = G.score_map(ctx, rows, codes, labels=["compartment", "cellcycle_phase"]
                        if "cellcycle_phase" in ctx.nodes else ["compartment"],
                        map_id="planted", columns=["expr_stage0"])
    top = table.set_index("label").loc["compartment"]
    assert top.category_f1 == pytest.approx(1.0)
    assert top.category_skill > 0.9
    assert top.coverage == pytest.approx(1.0)
    assert top.circular in (False, np.False_)
    # And an embedding-free check that real clustering of the planted features finds something.
    lab = cluster(X[:, :30], "hdbscan", **G.CLUSTERING)
    assert len(set(lab[lab >= 0])) >= 2


def test_circular_flag_comes_from_the_leakage_closure():
    from starplast import strategies as S
    ctx = S.planted_context(n=300, seed=2)
    rows = np.arange(ctx.n)
    clusters = np.zeros(ctx.n, int)
    t = G.score_map(ctx, rows, clusters, labels=["compartment"], columns=["compartment"])
    assert bool(t.circular.iloc[0]) and "compartment" in t.circular_columns.iloc[0]
    t = G.score_map(ctx, rows, clusters, labels=["compartment"], columns=None)
    assert t.circular.iloc[0] is None


# --------------------------------------------------------------------------- the shipped files
@pytest.fixture(scope="module")
def shipped():
    g = G.Gallery()
    if not g.available():
        pytest.skip("gallery not built")
    return g


def test_the_gallery_covers_every_shipped_parasite_space(shipped):
    assert set(shipped.organisms()) == set(G.space_codes())


def test_the_gallery_files_are_small():
    total = sum(os.path.getsize(os.path.join(G._data_dir(), f))
                for f in (G.COORDS_FILE, G.MANIFEST_FILE, G.SCORES_FILE))
    assert total < 15e6


@pytest.mark.parametrize("code", G.space_codes())
def test_every_map_is_whole_and_over_this_table(shipped, code):
    nodes = pd.read_parquet(organisms.nodes_path(code), columns=["gene_id"])
    ids = nodes.gene_id.astype(str).to_numpy()
    assert (shipped.gene_ids(code) == ids).all(), "built over a different gene table"
    m = shipped.manifest["organisms"][code]
    assert m["gene_ids_sha256"] == hashlib.sha256("\n".join(ids).encode()).hexdigest()
    maps = shipped.maps(code)
    kinds = {r["id"] for r in maps}
    assert {f"all_nn{k}" for k in G.ALL_NEIGHBORS} <= kinds
    assert "all_but_localization" in kinds
    assert sum(r["family"] not in ("all",) for r in maps) >= 4, "one map per evidence family"
    for r in maps:
        d = shipped.load(code, r["id"])
        rows, xyz, lab = d["rows"], d["xyz"], d["clusters"]
        assert xyz.dtype == np.float32 and xyz.shape == (len(rows), 3)
        assert np.isfinite(xyz).all()
        assert len(np.unique(rows)) == len(rows) and rows.min() >= 0 and rows.max() < len(ids)
        assert len(lab) == len(rows)
        assert r["n_genes"] == len(rows)
        assert r["n_clusters"] == len(set(lab[lab >= 0].tolist()))
        assert r["coordinates_sha256"] == hashlib.sha256(xyz.tobytes()).hexdigest()
        assert r["clustering"]["algorithm"] == "hdbscan"
        assert r["recipe"]["n_components"] == 3 and r["recipe"]["method"] == "umap"
        assert r["executed_method"] == "umap"


@pytest.mark.parametrize("code", G.space_codes())
def test_no_map_was_built_from_a_label(shipped, code):
    from starplast import strategies as S
    labels = set(S.shipped(code).categorical_columns())
    for r in shipped.maps(code):
        assert not labels & set(r["columns"]), r["id"]
        if r["id"] == "all_but_localization":
            ctx = S.shipped(code)
            assert all(ctx.family_of(b) != G.LEFT_OUT_FAMILY for b in r["blocks"])


@pytest.mark.parametrize("code", G.space_codes())
def test_every_label_is_scored_on_every_map(shipped, code):
    from starplast import strategies as S
    s = shipped.scores(code)
    labels = set(S.shipped(code).categorical_columns())
    maps = {r["id"] for r in shipped.maps(code)}
    assert set(s.label) == labels and set(s["map"]) == maps
    assert len(s) == len(labels) * len(maps)
    for c in ("category_f1", "category_precision", "category_recall", "category_skill",
              "best_category", "best_n", "best_f1_lower", "best_skill", "category_f1_chance",
              "best_f1_lower_chance", "coverage", "circular"):
        assert c in s.columns
    ok = s.category_f1.notna()
    assert ok.mean() > 0.9
    assert ((s.category_f1[ok] >= 0) & (s.category_f1[ok] <= 1)).all()
    assert (s.best_f1_lower[ok] <= 1).all()


def test_a_localization_label_is_circular_on_the_localization_map(shipped):
    s = shipped.scores(TG).set_index(["label", "map"])
    assert str(s.loc[("compartment", "family_localization"), "circular"]) == "True"
    assert str(s.loc[("compartment", "all_but_localization"), "circular"]) == "False"


def test_a_map_aligns_onto_a_reordered_table(shipped):
    ids = shipped.gene_ids(TG)
    m = shipped.load(TG, "all_nn25")
    rev = ids[::-1]
    a = shipped.aligned(TG, "all_nn25", rev)
    assert (rev[a["rows"]] == ids[m["rows"]]).all()
    assert np.array_equal(a["xyz"], m["xyz"])


def test_every_column_has_a_plain_name_and_a_tooltip():
    from starplast.umap_gallery_panel import TABLE_COLUMNS
    for c in ("map", "label") + TABLE_COLUMNS:
        name, tip = G.SCORE_COLUMNS[c]
        assert name and len(tip.split()) >= 6


# --------------------------------------------------------------------------- the dock
def test_the_panel_lists_the_maps_and_reads_the_table_both_ways(shipped):
    from starplast.umap_gallery_panel import MapGalleryPanel, LABEL_VIEW, MAP_VIEW
    from PyQt6 import QtWidgets
    p = MapGalleryPanel(TG, gallery=shipped)
    try:
        assert p.list.count() == len(shipped.maps(TG))
        assert not p.list.item(0).icon().isNull(), "each map has a thumbnail"
        assert p.label_box.currentText() == "compartment"
        assert p.table.rowCount() == len(shipped.maps(TG))
        headers = [p.table.horizontalHeaderItem(j).text() for j in range(p.table.columnCount())]
        assert headers[0] == "Map" and "Categories → clusters" in headers
        assert all(p.table.horizontalHeaderItem(j).toolTip() for j in range(p.table.columnCount()))
        p.view_box.setCurrentText(MAP_VIEW)
        assert p.table.rowCount() == shipped.scores(TG).label.nunique()
        chosen = []
        p.map_chosen.connect(chosen.append)
        p._on_item(p.list.item(2))
        assert chosen == [shipped.maps(TG)[2]["id"]]
        colored = []
        p.color_label.connect(colored.append)
        p.color_btn.click()
        assert colored == [p.label_box.currentText()]
        p.view_box.setCurrentText(LABEL_VIEW)
        for cls in (QtWidgets.QComboBox, QtWidgets.QPushButton, QtWidgets.QListWidget,
                    QtWidgets.QTableWidget):
            for w in p.findChildren(cls):
                assert w.toolTip().strip(), f"{cls.__name__} without a tooltip"
    finally:
        p.deleteLater()


def test_the_window_shows_a_gallery_map_and_colors_it_in_one_click(shipped):
    from starplast import app as A
    from starplast import umap_gallery_panel as P
    w = A.Window()
    try:
        p = w.maps_panel
        assert w.maps_dock.windowTitle() == "maps"
        rec = next(r for r in shipped.maps(TG) if r["id"] == "family_protein_abundance")
        before = np.array(w.xyz, copy=True)
        P.show_map(w, p, rec["id"])
        assert int(w.placed.sum()) == rec["n_genes"] < w.n
        assert not np.array_equal(before, w.xyz)
        assert w.cluster_labels is not None and len(w.cluster_labels) == w.n
        assert (w.cluster_labels[~w.placed] == -1).all()
        p.label_box.setCurrentText("lopit_unified")
        p.color_btn.click()
        assert w.category == "lopit_unified" and w.color_mode == "compartment"
        p.clusters_btn.click()
        assert w.color_mode == "clusters"
        scores = P.score_on_screen(w, p)
        assert set(scores["map"]) == {P.ON_SCREEN}
        assert "compartment" in set(scores.label)
        assert p.map_box.currentData() == P.ON_SCREEN and p.table.rowCount() == len(scores)
    finally:
        w.console.remove()
        w.close()
