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


# --------------------------------------------------------------------------- structure, label-free
def test_separation_is_high_for_far_apart_clusters_and_a_half_for_one_cloud():
    rng = np.random.default_rng(0)
    far = np.vstack([rng.normal(0, 0.1, (200, 3)), rng.normal(20, 0.1, (200, 3))])
    labels = np.array([0] * 200 + [1] * 200)
    assert G.separation(far, labels) > 0.95
    # The same labels over one cloud: the points of a "cluster" are no nearer each other than they
    # are to the other cluster's, so the silhouette is about 0 and the rescaled value about 0.5.
    one = rng.normal(0, 1, (400, 3))
    assert 0.4 < G.separation(one, labels) < 0.6
    assert not np.isfinite(G.separation(far, np.zeros(400, dtype=int))), "one cluster is not two"
    assert not np.isfinite(G.separation(far, np.full(400, -1)))


def test_structure_needs_both_a_real_partition_and_real_separation():
    rng = np.random.default_rng(1)
    # Even, fully clustered AND separated: good at both, so a high geometric mean.
    good_xyz = np.vstack([rng.normal(k * 20, 0.2, (100, 3)) for k in range(4)])
    good_lab = np.repeat(np.arange(4), 100)
    good = G.structure(good_xyz, good_lab)
    assert good["structure"] > 0.9
    # Even and fully clustered, but the clusters overlap completely: map_quality is still high and
    # the structure score is not, which is the whole reason for the second term.
    flat = G.structure(rng.normal(0, 1, (400, 3)), good_lab)
    assert flat["score"] > 0.9 and flat["structure"] < 0.75
    assert flat["structure"] < good["structure"]
    # Separated, but 95% of the genes unclustered: the silhouette is perfect, the partition is not.
    thin_lab = np.where(np.arange(400) < 20, good_lab, -1)
    thin = G.structure(good_xyz, thin_lab)
    assert thin["separation"] > 0.9 and thin["structure"] < good["structure"]
    # A bisection is refused outright by `clustering.degenerate`, through `search.map_quality`.
    bisect = G.structure(good_xyz, np.array([0] * 200 + [1] * 200))
    assert bisect["structure"] == 0.0


def test_the_structure_score_reads_no_label():
    """`structure` takes coordinates and clusters only: there is no way to pass it a label."""
    import inspect
    assert list(inspect.signature(G.structure).parameters) == ["xyz", "labels"]


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
    # Dozens of maps per organism, so the budget is larger than the first gallery's -- but the
    # coordinates are float32 and the cluster labels int16, and the whole gallery stays well inside
    # what the wheel can carry (`scripts/check_wheel.py`).
    total = sum(os.path.getsize(os.path.join(G._data_dir(), f))
                for f in (G.COORDS_FILE, G.MANIFEST_FILE, G.SCORES_FILE))
    assert total < 25e6, f"{total / 1e6:.1f} MB of gallery"


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
        assert r["group"] in G.GROUPS, r["id"]
        assert 0.0 <= r["structure"] <= 1.0, r["id"]
        assert r["tuned"] == (not r["id"].startswith("all_nn")), r["id"]
        assert (len(r["settings_searched"]) > 0) == r["tuned"], r["id"]
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


@pytest.mark.parametrize("code", G.space_codes())
def test_the_gallery_is_large_and_covers_every_group(shipped, code):
    """Every category on its own and the sensible combinations, not just a handful of maps."""
    maps = shipped.maps(code)
    assert len(maps) >= 40, f"{code} ships only {len(maps)} maps"
    groups = {g: [r for r in maps if r["group"] == g] for g in G.GROUPS}
    assert len(groups[G.GROUPS[0]]) == len(G.ALL_NEIGHBORS)
    for g in G.GROUPS[1:]:
        assert len(groups[g]) >= 2, f"{code}: {g} has {len(groups[g])} map(s)"
    # Single experiments: one map per individual block, each from exactly one block.
    assert all(len(r["blocks"]) == 1 for r in groups[G.GROUPS[3]])
    assert len(groups[G.GROUPS[3]]) >= 10
    # Pairs and triples are built from the stated number of families.
    assert all(len(r["family"].split(" + ")) == 2 for r in groups[G.GROUPS[4]])
    assert all(len(r["family"].split(" + ")) == 3 for r in groups[G.GROUPS[5]])
    assert len({r["id"] for r in maps}) == len(maps), "duplicate map ids"


@pytest.mark.parametrize("code", G.space_codes())
def test_an_all_but_map_never_saw_the_family_it_leaves_out(shipped, code):
    from starplast import strategies as S
    ctx = S.shipped(code)
    left_out = [r for r in shipped.maps(code) if r["id"].startswith("all_but_")]
    assert len(left_out) >= 2
    for r in left_out:
        fam = r["id"][len("all_but_"):]
        assert all(G._slug(ctx.family_of(b)) != fam for b in r["blocks"]), r["id"]


@pytest.mark.parametrize("code", G.space_codes())
def test_the_search_tried_the_grid_and_the_winner_is_the_best_of_it(shipped, code):
    """Each tuned map ships the highest-structure setting the search measured at full size."""
    for r in shipped.maps(code):
        if not r["tuned"]:
            continue
        tried = r["settings_searched"]
        sample = [t for t in tried if t["stage"] == "sample"]
        full = [t for t in tried if t["stage"] == "full"]
        assert len(sample) == len(G.MAP_GRID["n_neighbors"]) * len(G.MAP_GRID["min_dist"]), r["id"]
        assert 1 <= len(full) <= G.SEARCH_KEEP, r["id"]
        best = max(full, key=lambda t: t["structure"])
        assert r["structure"] == best["structure"], r["id"]
        assert r["n_neighbors"] == best["n_neighbors"] and r["min_dist"] == best["min_dist"]
        assert r["clustering"]["min_cluster_size"] == best["min_cluster_size"]
        assert r["clustering"]["min_samples"] == best["min_samples"]
        assert r["clustering"]["cluster_selection_method"] == "leaf"
        for t in tried:
            assert t["n_neighbors"] in G.MAP_GRID["n_neighbors"]
            assert t["min_dist"] in G.MAP_GRID["min_dist"]
            assert t["min_cluster_size"] in G.GALLERY_CLUSTER_GRID["min_cluster_size"]
            assert t["min_samples"] in G.GALLERY_CLUSTER_GRID["min_samples"]


def test_the_search_did_not_simply_pick_one_setting_everywhere(shipped):
    """If the grid always returned the same setting, searching it would be hours wasted."""
    chosen = {(r["n_neighbors"], r["min_dist"],
               r["clustering"]["min_cluster_size"], r["clustering"]["min_samples"])
              for r in shipped.maps(TG) if r["tuned"]}
    assert len(chosen) >= 3, chosen


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
        assert [i for i, _ in p.items()] == [r["id"] for r in shipped.maps(TG)]
        assert p.label_box.currentText() == "compartment"
        assert p.table.rowCount() == len(shipped.maps(TG))
        headers = [p.table.horizontalHeaderItem(j).text() for j in range(p.table.columnCount())]
        assert headers[0] == "Map" and "Categories → clusters" in headers
        assert all(p.table.horizontalHeaderItem(j).toolTip() for j in range(p.table.columnCount()))
        p.view_box.setCurrentText(MAP_VIEW)
        assert p.table.rowCount() == shipped.scores(TG).label.nunique()
        chosen = []
        p.map_chosen.connect(chosen.append)
        third = shipped.maps(TG)[2]["id"]
        p._on_item(third)
        assert chosen == [third]
        colored = []
        p.color_label.connect(colored.append)
        p.color_btn.click()
        assert colored == [p.label_box.currentText()]
        p.view_box.setCurrentText(LABEL_VIEW)
        for cls in (QtWidgets.QComboBox, QtWidgets.QPushButton, QtWidgets.QTreeWidget,
                    QtWidgets.QTableWidget, QtWidgets.QLineEdit, QtWidgets.QCheckBox):
            for w in p.findChildren(cls):
                assert w.toolTip().strip(), f"{cls.__name__} without a tooltip"
    finally:
        p.deleteLater()


def _cells(table) -> list:
    """Every cell of a table as text: what the user can actually read off it."""
    return [[table.item(i, j).text() if table.item(i, j) else None
             for j in range(table.columnCount())] for i in range(table.rowCount())]


def test_choosing_a_second_map_changes_the_score_table(shipped):
    """The 0.49.0 bug: the table changed for the first map chosen and for no map after it.

    The default view lists every map for one label, so its rows are the same whichever map is
    current -- and nothing in it moved when a second map was picked, because `refresh` read only the
    view and the label. Both views now follow the current map, so A, then B, then A again gives
    three different tables and the third equals the first.
    """
    from starplast.umap_gallery_panel import MapGalleryPanel, LABEL_VIEW, MAP_VIEW, SHOWN_MARK
    maps = [r["id"] for r in shipped.maps(TG)]
    a, b = maps[4], maps[7]
    for view in (LABEL_VIEW, MAP_VIEW):
        p = MapGalleryPanel(TG, gallery=shipped)
        try:
            p.view_box.setCurrentText(view)
            p._on_item(a)
            first = _cells(p.table)
            p._on_item(b)
            second = _cells(p.table)
            p._on_item(a)
            third = _cells(p.table)
            assert first != second, f"{view}: choosing a second map left the first map's table"
            assert second != third, f"{view}: choosing a third map left the second map's table"
            assert first == third, f"{view}: the same map gave two different tables"
            assert p.current_map == a
            if view == LABEL_VIEW:
                # Every map is still listed, and exactly one row is marked as the one on screen.
                assert len(first) == len(maps)
                marked = [r[0] for r in third if r[0].startswith(SHOWN_MARK)]
                assert marked == [SHOWN_MARK + p.title(a)]
                assert [r[0] for r in second if r[0].startswith(SHOWN_MARK)] == \
                    [SHOWN_MARK + p.title(b)]
            else:
                # The map view is filtered, so the whole body is that map's scores.
                assert p.map_box.currentData() == a
                assert len(first) == shipped.scores(TG).label.nunique()
        finally:
            p.deleteLater()


def test_the_map_combo_and_a_double_clicked_row_choose_a_map_too(shipped):
    """The other two ways in: the combo in map view, and double-clicking a row in label view."""
    from starplast.umap_gallery_panel import MapGalleryPanel, LABEL_VIEW, MAP_VIEW
    maps = [r["id"] for r in shipped.maps(TG)]
    p = MapGalleryPanel(TG, gallery=shipped)
    try:
        p.view_box.setCurrentText(MAP_VIEW)
        p.map_box.setCurrentIndex(p.map_box.findData(maps[2]))
        one = _cells(p.table)
        assert p.current_map == maps[2]
        p.map_box.setCurrentIndex(p.map_box.findData(maps[5]))
        assert _cells(p.table) != one and p.current_map == maps[5]
        p.view_box.setCurrentText(LABEL_VIEW)
        chosen = []
        p.map_chosen.connect(chosen.append)
        row = next(i for i in range(p.table.rowCount()) if p._row_key(i) == maps[9])
        p._row_activated(row, 0)
        assert chosen == [maps[9]] and p.current_map == maps[9]
    finally:
        p.deleteLater()


def test_the_panel_filters_sorts_and_groups_a_large_gallery(shipped):
    """With dozens of maps the list has to be searchable, sortable and groupable."""
    from starplast.umap_gallery_panel import (MapGalleryPanel, SORT_STRUCTURE, SORT_LABEL,
                                             SORT_GENES, SORT_GROUP)
    from PyQt6 import QtCore
    USER = QtCore.Qt.ItemDataRole.UserRole
    p = MapGalleryPanel(TG, gallery=shipped)
    try:
        everything = [i for i, _ in p.items()]
        assert len(everything) == len(shipped.maps(TG))
        # Grouped by default: the top level is headings, and every heading is a known group.
        tops = [p.tree.topLevelItem(i).text(0) for i in range(p.tree.topLevelItemCount())]
        assert tops and all(t.split(" — ")[0] in G.GROUPS for t in tops)
        assert all(p.tree.topLevelItem(i).data(0, USER) is None
                   for i in range(p.tree.topLevelItemCount())), "a heading is not a map"
        # Filtering: every map kept mentions the word somewhere a user can see.
        p.search.setText("fitness")
        narrowed = [i for i, _ in p.items()]
        assert 0 < len(narrowed) < len(everything)
        assert all("fitness" in " ".join(str(p.record(i)[k]) for k in
                                         ("id", "title", "group", "family", "description")).lower()
                   for i in narrowed)
        p.search.setText("no such evidence anywhere")
        assert p.items() == []
        p.search.setText("")
        assert len([i for i, _ in p.items()]) == len(everything)
        # Sorting, flat so the order is the global one.
        p.group_check.setChecked(False)
        for key, value in ((SORT_STRUCTURE, lambda i: p.record(i)["structure"]),
                           (SORT_GENES, lambda i: p.record(i)["n_genes"]),
                           (SORT_LABEL, lambda i: p.label_skill().get(i, -1e9))):
            p.sort_box.setCurrentText(key)
            got = [value(i) for i, _ in p.items()]
            assert got == sorted(got, reverse=True), key
        p.sort_box.setCurrentText(SORT_GROUP)
        assert [i for i, _ in p.items()] == everything
        # The rows say what the map is and how good it is, and the thumbnails are drawn.
        p.group_check.setChecked(True)
        p.tree.topLevelItem(0).setExpanded(True)
        first_id, first_item = p.items()[0]
        assert not first_item.icon(0).isNull(), "a shown map has a thumbnail"
        assert f"{p.record(first_id)['structure']:.2f}" in first_item.text(0)
        assert "structure" in first_item.toolTip(0) and "n_neighbors" in first_item.toolTip(0)
        assert all(it.toolTip(0).strip() for _, it in p.items())
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
