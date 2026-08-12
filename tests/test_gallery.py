#!/usr/bin/env python3
"""The walk gallery: thumbnails as a walk runs, and what happens when one is opened.

Two failures this file exists to prevent, both of which have already happened elsewhere in this
project in a different form.

The first is a picture that is not the map. A thumbnail is a projection painted by hand, so nothing
guarantees it corresponds to the embedding it claims to show except tests that check the arithmetic:
the same scale on both axes (a map that is genuinely elongated must not be squared up), y pointing
the right way, and near points drawn over far ones rather than summed. The rendering bug that made
every colour mode draw as one white blob passed every array-level test in the suite.

The second is absence drawn as data. A walk embeds a seeded subsample, so most genes have NO position
in one of its maps -- and the code that displays a rebuilt embedding used to leave every one of them
at the origin, where they formed a lump that could be clicked, filtered and counted like structure.
So these tests assert on what is drawn AND on what can be picked, not on the coordinate array alone.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import gallery as G  # noqa: E402
from starplast.embedding import EmbeddingSpec  # noqa: E402
from starplast.tuning import WalkStep  # noqa: E402


@pytest.fixture(scope="module")
def app():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _step(i=1, total=4, n=40, n_genes=100, nn=15, md=0.1, seed=0, scores=True):
    """One walk configuration, shaped exactly as `tuning.walk_umap_iter` yields them."""
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(n, 3)) * 10.0
    genes = np.zeros(n_genes, dtype=bool)
    genes[:n] = True
    row = {"n_neighbors": nn, "min_dist": md, "seed": 42, "sample_size": n}
    if scores:
        row.update(trustworthiness=0.912, continuity_proxy=0.5,
                   n_clusters_hdbscan=7, noise_frac=0.25)
    return WalkStep(index=i, total=total, row=row, coords=coords, genes=genes,
                    spec=EmbeddingSpec(n_neighbors=nn, min_dist=md), name=f"walk_nn{nn}")


# --------------------------------------------------------------------------- projection
def test_both_axes_are_scaled_by_the_same_factor(app):
    """A map that is genuinely elongated must look elongated. Scaling each axis to fill the box
    would make every configuration look equally isotropic -- which is the difference between two
    settings that the gallery exists to show."""
    coords = np.zeros((4, 3))
    coords[:, 0] = [-10, 10, 0, 0]      # wide
    coords[:, 1] = [0, 0, -1, 1]        # narrow
    x, y = G.project(coords, 100)
    assert np.ptp(x) > 10 * np.ptp(y)


def test_the_image_is_not_a_mirror_of_the_map(app):
    """The embedding's y points up and a raster's points down. Without the flip every thumbnail is a
    mirror image of the map it opens into, which is worse than no thumbnail."""
    coords = np.array([[0.0, -5.0, 0.0], [0.0, 5.0, 0.0]])
    _, y = G.project(coords, 100)
    assert y[0] > y[1], "the point that is lower in the map must be lower in the image"


def test_a_cloud_with_no_extent_is_centred_rather_than_dividing_by_zero(app):
    """Every point identical is a real outcome -- a degenerate embedding -- and it must produce a
    picture of a degenerate embedding, not a crash inside a paint call."""
    coords = np.zeros((5, 3))
    x, y = G.project(coords, 100)
    assert np.allclose(x, 50.0) and np.allclose(y, 50.0)


# --------------------------------------------------------------------------- thumbnails
def test_a_thumbnail_is_the_size_asked_for_and_has_points_on_it(app):
    img = G.thumbnail(_step().coords, size=120)
    assert (img.width(), img.height()) == (120, 120)
    px = np.array([[img.pixelColor(x, y).getRgb()[:3] for x in range(0, 120, 3)]
                   for y in range(0, 120, 3)])
    assert len(np.unique(px.reshape(-1, 3), axis=0)) > 1, "nothing was drawn"


def test_an_empty_walk_step_paints_the_background_rather_than_failing(app):
    """A configuration can legitimately produce no coordinates -- every gene dropped by the
    missing-value policy -- and a viewer must show that rather than raise from a paint call."""
    img = G.thumbnail(np.zeros((0, 3)), size=40, background=(1.0, 0.0, 0.0))
    assert img.pixelColor(20, 20).getRgb()[:3] == (255, 0, 0)


def test_coordinates_without_a_third_component_still_draw(app):
    """`n_components=2` is a legal spec, and indexing a missing z is the obvious way to break it."""
    img = G.thumbnail(np.random.default_rng(0).normal(size=(30, 2)), size=60)
    assert img.width() == 60


def test_colours_are_used_when_they_line_up_and_ignored_when_they_do_not(app):
    """A colour array of the wrong length means the caller's colouring does not describe these genes.
    Painting it anyway would colour gene 12 with gene 40's category -- a picture that is confidently
    wrong, which is worse than the neutral default."""
    step = _step(n=30)
    red = np.tile([1.0, 0.0, 0.0], (30, 1))
    lit = G.thumbnail(step.coords, red, size=80)
    wrong = G.thumbnail(step.coords, np.tile([1.0, 0.0, 0.0], (7, 1)), size=80)
    assert lit != wrong


def test_a_thumbnail_takes_rgb_as_well_as_rgba(app):
    """`colours` returns RGBA and a caller's own palette is usually RGB; refusing one of them is a
    papercut in the only path that gives a thumbnail any meaning."""
    step = _step(n=20)
    assert G.thumbnail(step.coords, np.tile([0.2, 0.9, 0.4], (20, 1)), size=50).width() == 50
    assert G.thumbnail(step.coords, np.tile([0.2, 0.9, 0.4, 1.0], (20, 1)), size=50).width() == 50


# --------------------------------------------------------------------------- captions
def test_the_caption_leads_with_the_configuration_then_the_scores(app):
    text = G.caption(_step(nn=25, md=0.5))
    assert text.startswith("nn=25") and "md=0.5" in text
    assert "trust 0.912" in text and "7 clusters" in text and "25% noise" in text


def test_a_walk_without_the_cluster_check_says_nothing_about_clusters(app):
    """"clusters=nan" reads as a failure rather than as a question that was not asked."""
    text = G.caption(_step(scores=False))
    assert "clusters" not in text and "trust" not in text


def test_an_undefined_trustworthiness_is_omitted_rather_than_printed(app):
    """`_quality` returns nan on a degenerate embedding, deliberately, and nan is not a score."""
    step = _step()
    step.row["trustworthiness"] = float("nan")
    assert "trust" not in G.caption(step)


def test_something_that_is_not_a_walk_step_captions_without_raising(app):
    assert G.caption(object()) == "nn=?  ·  md=?"


# --------------------------------------------------------------------------- the panel
@pytest.fixture
def panel(app):
    return G.GalleryPanel()


def test_it_starts_empty_and_says_what_to_do(panel):
    assert panel.steps == [] and panel.list.count() == 0
    assert "no maps yet" in panel.count.text()
    assert panel.hint.isVisibleTo(panel)


def test_each_configuration_appears_as_it_arrives(panel):
    """The whole point: run 12 is on screen while run 13 computes. A gallery that fills at the end
    is a progress bar with pictures."""
    panel.add(_step(i=1, total=3))
    assert panel.list.count() == 1 and "1 of 3 computed" in panel.count.text()
    panel.add(_step(i=2, total=3, seed=1))
    assert panel.list.count() == 2 and "2 of 3 computed" in panel.count.text()
    assert not panel.hint.isVisibleTo(panel)


def test_a_walk_of_unknown_length_still_reports_what_it_has(panel):
    panel.add(_step(total=0))
    assert panel.count.text() == "1 computed"


def test_the_view_follows_the_newest_map_until_the_user_steps_back(panel):
    """The log-tail rule. A walk that yanked the view to run 13 while run 12 was being read would
    make the gallery useless during exactly the period it exists for."""
    for k in range(3):
        panel.add(_step(i=k + 1, seed=k))
    assert panel.slider.value() == 2
    panel.slider.setValue(0)
    panel.add(_step(i=4, seed=9))
    assert panel.slider.value() == 0, "the walk moved the view off what was being looked at"


def test_scroll_mode_renders_one_configuration_large(panel):
    panel.add(_step())
    panel.set_mode("scroll")
    assert panel.stack.currentIndex() == 1
    assert panel.big.pixmap() is not None and not panel.big.pixmap().isNull()
    assert "1 of 1" in panel.caption.text()
    panel.set_mode("grid")
    assert panel.stack.currentIndex() == 0


def test_the_chooser_says_which_mode_is_showing(panel):
    """Rendered and looked at: the panel was showing the scroll page while the control still read
    "grid", because the mode was set in code and only the combo box's own signal updated it."""
    panel.add(_step())
    panel.set_mode("scroll")
    assert panel.mode_box.currentText() == "scroll"
    panel.mode_box.setCurrentText("grid")
    assert panel.stack.currentIndex() == 0


def test_stepping_past_the_end_does_nothing_rather_than_raising(panel):
    panel.show_index(5)
    assert panel.caption.text() == ""


def test_clicking_a_thumbnail_asks_for_that_map(panel):
    """The gallery picks and the main view displays; without this signal a thumbnail is a picture."""
    got = []
    panel.chosen.connect(got.append)
    panel.add(_step(i=1))
    item = panel.add(_step(i=2, seed=3))
    panel._on_item(item)
    assert got and got[0] is panel.steps[1]
    assert panel.slider.value() == 1, "the scroll view should follow the tile that was clicked"


def test_an_item_that_names_no_configuration_is_ignored(panel, app):
    from PyQt6 import QtWidgets
    got = []
    panel.chosen.connect(got.append)
    panel._on_item(QtWidgets.QListWidgetItem("stray"))
    assert got == []


def test_the_open_button_asks_for_the_one_on_screen(panel):
    got = []
    panel.chosen.connect(got.append)
    panel.expand_current()
    assert got == [], "nothing has been computed yet; there is nothing to open"
    panel.add(_step())
    panel.expand_current()
    assert got == [panel.steps[0]]


def test_a_new_walk_clears_the_last_one(panel):
    """Two walks read as one is a comparison between configurations that were never compared."""
    panel.add(_step())
    panel.clear()
    assert panel.steps == [] and panel.list.count() == 0
    assert panel.slider.maximum() == 0 and panel.caption.text() == ""
    assert "no maps yet" in panel.count.text() and panel.hint.isVisibleTo(panel)


def test_current_is_none_before_anything_has_been_computed(panel):
    assert panel.current() is None


# --------------------------------------------------------------------------- colouring
def test_colours_come_from_the_caller_so_the_thumbnail_matches_the_map(app):
    seen = {}

    def colour_fn(mask):
        seen["n"] = int(np.sum(mask))
        return np.tile([1.0, 0.0, 0.0, 1.0], (int(np.sum(mask)), 1))

    p = G.GalleryPanel(colour_fn=colour_fn)
    step = _step(n=40)
    assert p.colours_for(step).shape == (40, 4)
    assert seen["n"] == 40


def test_a_colouring_of_the_wrong_length_is_refused_rather_than_recycled(app):
    p = G.GalleryPanel(colour_fn=lambda mask: np.zeros((3, 4)))
    assert p.colours_for(_step(n=40)) is None


def test_a_colouring_that_fails_costs_the_colour_not_the_picture(app, capsys):
    """A viewer must not raise out of a signal handler because a column was dropped."""
    def boom(mask):
        raise KeyError("compartment")

    p = G.GalleryPanel(colour_fn=boom)
    assert p.colours_for(_step()) is None
    assert "colouring unavailable" in capsys.readouterr().out
    p.add(_step())               # and it still draws
    assert p.list.count() == 1


def test_without_a_colour_function_thumbnails_still_draw(panel):
    assert panel.colours_for(_step()) is None
    panel.add(_step())
    assert panel.list.count() == 1


def test_every_gallery_control_explains_itself(panel):
    """The useful tooltip says WHY. "grid" as a tooltip on a control labelled grid is not one."""
    for attr in ("mode_box", "count", "list", "big", "slider", "open_btn"):
        tip = getattr(panel, attr).toolTip()
        assert tip and len(tip.split()) >= 15, f"{attr}: {tip!r}"


# --------------------------------------------------------------------------- in the window
@pytest.fixture(scope="module")
def win(app):
    from starplast.app import Window
    w = Window()
    w.resize(900, 600)
    yield w
    w.close()


@pytest.fixture
def restored(win):
    """Put the window's map back after a test replaces it with a walk configuration."""
    xyz, placed, sel = win.xyz.copy(), win.placed, win.sel
    yield
    win.xyz = xyz
    win.view.xyz = xyz
    win.placed = placed
    win.view.pickable = placed
    win.sel = sel
    win.redraw()


def test_the_gallery_is_present_and_out_of_the_way_until_a_walk_runs(win):
    assert win.gallery_dock.isHidden()
    assert win.gallery.steps == []


def test_starting_a_walk_clears_the_gallery_and_shows_it(win):
    win.gallery.add(_step(n_genes=win.n))
    win.panel.walk_started.emit()
    assert win.gallery.steps == []
    assert not win.gallery_dock.isHidden()
    win.gallery_dock.hide()


def test_a_configuration_reaches_the_gallery_as_the_walk_emits_it(win):
    win.gallery.clear()
    win.panel.walk_step.emit(_step(n=60, n_genes=win.n))
    assert win.gallery.list.count() == 1


def test_thumbnails_are_coloured_the_way_the_map_is(win):
    """Same function, so "grey means unknown" cannot mean one thing on the map and another here."""
    mask = np.zeros(win.n, bool)
    mask[:50] = True
    c = win.colours_for_genes(mask)
    assert c.shape == (50, 4)
    assert np.allclose(c, win.colours(np.ones(win.n, bool))[:50])


def test_opening_a_configuration_shows_it_in_the_central_view(win, restored):
    """Not a picture of a map: the same 3D view, so a gene can be clicked exactly as before."""
    step = _step(n=60, n_genes=win.n)
    win.gallery.clear()
    win.gallery.add(step)
    win.gallery.expand_current()
    assert win.placed.sum() == 60
    assert np.allclose(win.xyz[step.genes], step.coords.astype(np.float32))
    assert win.view.xyz is win.xyz, "the picker must see the same coordinates as the scatter"
    assert "not in this embedding" in win.status.currentMessage()


def test_genes_outside_the_walk_sample_are_hidden_rather_than_stacked_at_the_origin(win, restored):
    """They have no position in this embedding. Left at the origin they were a dense lump in the
    middle of the map that could be clicked, filtered and counted like real structure."""
    step = _step(n=60, n_genes=win.n)
    win.use_embedding(step.coords, step.genes)
    assert win.visible_mask().sum() == 60
    alpha = np.asarray(win.scatter.color)[:, 3]
    assert np.allclose(alpha[~step.genes], 0.0)
    assert alpha[step.genes].max() > 0.1


def test_an_unplaced_gene_cannot_be_picked(win, restored):
    """Hiding is not enough: an invisible point is still the nearest point to a click where it sits,
    so the evidence panel would open on a gene that is not in the map."""
    from PyQt6 import QtCore, QtGui
    step = _step(n=60, n_genes=win.n)
    win.use_embedding(step.coords, step.genes)
    sx, sy = win.view.project()
    hidden = np.where(~step.genes & np.isfinite(sx) & np.isfinite(sy))[0]
    assert len(hidden), "no unplaced gene projected onto the widget"
    i = int(hidden[0])
    got = []
    win.view.picked.connect(got.append)
    pos = QtCore.QPointF(float(sx[i]), float(sy[i]))
    win.view.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.KeyboardModifier.NoModifier))
    assert all(step.genes[g] for g in got), "clicked through to a gene with no position"


def test_a_selection_that_the_new_map_does_not_cover_is_dropped(win, restored):
    """Otherwise the halo and the evidence panel go on describing a gene that is not on screen."""
    step = _step(n=60, n_genes=win.n)
    win.sel = int(np.where(~step.genes)[0][0])
    win.use_embedding(step.coords, step.genes)
    assert win.sel is None


def test_a_full_embedding_leaves_nothing_hidden(win, restored):
    """The ordinary case -- "build this map" over every gene -- must not acquire a caveat."""
    coords = win.xyz.copy()
    win.use_embedding(coords, np.ones(win.n, bool))
    assert win.placed.all() and "hidden" not in win.status.currentMessage()


def test_an_index_array_says_which_genes_as_clearly_as_a_mask_does(win, restored):
    """`embed` returns a boolean mask and a subsample is naturally a list of positions; a viewer that
    accepted only one of them would put the coordinates on the wrong genes."""
    idx = np.arange(0, 300, 3)
    win.use_embedding(win.xyz[idx].copy(), idx)
    assert win.placed.sum() == len(idx) and win.placed[idx].all()


def test_the_theme_repaints_the_gallery_ground(win):
    """A dark thumbnail on a light page reads as an image that failed to load."""
    win.gallery.clear()
    win.gallery.add(_step(n_genes=win.n))
    win.apply_theme("paper")
    assert max(win.gallery.background) > 0.5
    win.apply_theme("dark")
    assert max(win.gallery.background) < 0.5
