#!/usr/bin/env python3
"""The cell diagram: the map's palette, on a parasite, clickable in both directions.

The point of it is that a list of 27 names becomes a picture of where the genes are, and the risk of
it is that a picture is persuasive whether or not it is true. So these tests are about the joins: the
colours come from the map's own dict rather than a second palette, a shape shared by three classes
says which one it is showing, and a compartment the artwork has no organelle for is NAMED rather
than quietly dropped.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import celldiagram as CD  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def svg():
    return open(CD.icon_path(), encoding="utf8").read()


@pytest.fixture
def diagram(qapp):
    d = CD.CellDiagram()
    d.resize(220, 320)
    return d


# --------------------------------------------------------------------------- the mapping
def test_the_artwork_ships_with_the_package(svg):
    assert CD.available() and len(svg) > 100_000


def test_every_mapped_compartment_names_an_organelle_that_exists(svg):
    """A mapping to a shape the drawing does not have is a silent no-op: the compartment would
    simply never colour anything and nothing would say why."""
    have = CD.groups_in(svg)
    missing = {c: sl for c, sl in CD.COMPARTMENT_SL.items() if sl not in have}
    assert not missing, missing


def test_the_compartments_with_no_organelle_are_named_rather_than_dropped(svg):
    """Between the ribosomes, the proteasome and the apical classes these are a large part of the
    proteome, and a compartment that vanishes from the legend reads as one that does not exist."""
    classes = ["dense granules", "40S ribosome", "19S proteasome", "apical 1",
               "endomembrane vesicles", "rhoptries 1"]
    absent = CD.missing_from_drawing(classes, svg)
    assert "40S ribosome" in absent and "19S proteasome" in absent and "apical 1" in absent
    assert "endomembrane vesicles" in absent, "a mixed class must not be drawn as a COPI vesicle"
    assert "dense granules" not in absent and "rhoptries 1" not in absent


def test_the_classes_that_share_a_shape_are_known(svg):
    """Three plasma-membrane classes on one shape is the case the whole selection dance exists for."""
    assert set(CD.sharing("SL0233")) == {"rhoptries 1", "rhoptries 2"}
    assert set(CD.sharing("SL0039")) == {"PM - integral", "PM - peripheral 1", "PM - peripheral 2"}
    assert CD.sharing("SL0018") == ["apicoplast"]
    assert CD.sharing("SL9999") == []


# --------------------------------------------------------------------------- colouring
def test_the_fill_reaches_the_paths_not_only_the_group(svg):
    """The artwork sets `fill` on the individual paths, so a fill on the group alone is overridden
    by every path in it -- the diagram comes back grey while every assertion on the string passes."""
    out = CD.recolour(svg, {"SL0018": (1.0, 0.0, 0.0)})
    import re
    group = re.search(r'<g[^>]*id="SL0018".*?</g>', out, re.S).group(0)
    assert "#ff0000" in group
    assert not re.search(r'fill\s*:\s*(?!#ff0000)#[0-9a-f]{6}', group), "a path kept its own fill"


def test_a_shape_nobody_selected_is_neutral_not_the_first_class_s_colour(diagram):
    """Filling a shared shape with whichever class sorts first would be a claim nobody made."""
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)}, "")
    assert diagram.fills()["SL0233"] == CD.NEUTRAL


def test_a_shared_shape_takes_the_colour_of_whichever_class_is_selected(diagram):
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 2")
    assert diagram.fills()["SL0233"] == (0.0, 1.0, 0.0)
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 1")
    assert diagram.fills()["SL0233"] == (1.0, 0.0, 0.0)


def test_a_shared_shape_says_which_class_it_is_showing(diagram):
    """A user coming back to the window cannot otherwise tell whether the rhoptry is coloured for 1
    or for 2."""
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 2")
    said = diagram.showing()
    assert "rhoptries 2" in said and "rhoptries 1" in said
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "apicoplast")
    assert diagram.showing() == "", "a shape of its own has nothing to disambiguate"


def test_an_unshared_shape_takes_its_colour_whether_or_not_it_is_selected(diagram):
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    assert diagram.fills()["SL0018"] == (0.0, 0.0, 1.0)


def test_a_compartment_the_palette_does_not_have_colours_nothing(diagram):
    diagram.set_palette({}, "")
    assert diagram.fills() == {}


def test_the_drawing_renders(diagram):
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "apicoplast")
    r = diagram.renderer()
    assert r is not None and r.isValid() and r.defaultSize().width() > 0
    from PyQt6 import QtGui
    img = QtGui.QImage(220, 320, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    p = QtGui.QPainter(img)
    diagram.render(p)
    p.end()
    assert any(QtGui.QColor.fromRgba(img.pixel(x, y)).alpha() > 0
               for x in range(0, 220, 4) for y in range(0, 320, 4)), "nothing was drawn"


def test_a_missing_artwork_hides_the_diagram_rather_than_raising(qapp, tmp_path):
    d = CD.CellDiagram(path=str(tmp_path / "nothing.svg"))
    assert d.svg == "" and d.groups == set() and d.renderer() is None
    d.set_palette({"apicoplast": (0, 0, 1)}, "apicoplast")
    assert d.fills() == {}
    assert d.organelle_at(10, 10) == ""
    assert not CD.available(str(tmp_path / "nothing.svg"))


# --------------------------------------------------------------------------- clicking
def test_clicking_an_organelle_selects_its_compartment(diagram):
    """Both directions, or the diagram is decoration."""
    from PyQt6 import QtCore, QtGui
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    box = diagram.renderer().boundsOnElement("SL0018")
    target, size = diagram._target(), diagram.renderer().defaultSize()
    x = target.x() + box.center().x() * target.width() / size.width()
    y = target.y() + box.center().y() * target.height() / size.height()
    got = []
    diagram.compartment_clicked.connect(got.append)
    pos = QtCore.QPointF(x, y)
    diagram.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert got == ["apicoplast"]


def test_clicking_a_shared_shape_steps_through_the_classes_it_stands_for(diagram):
    """One shape stands for three plasma-membrane classes; without a cycle two of them would be
    unreachable from the drawing."""
    palette = {c: (0.5, 0.5, 0.5) for c in
               ("PM - integral", "PM - peripheral 1", "PM - peripheral 2")}
    diagram.set_palette(palette, "")
    seen = []
    for _ in range(4):
        classes = [c for c in CD.sharing("SL0039") if c in diagram.colour_of]
        i = classes.index(diagram.selected) + 1 if diagram.selected in classes else 0
        nxt = classes[i % len(classes)]
        seen.append(nxt)
        diagram.set_palette(palette, nxt)
    assert seen == ["PM - integral", "PM - peripheral 1", "PM - peripheral 2", "PM - integral"]


def test_clicking_empty_space_selects_nothing(diagram):
    from PyQt6 import QtCore, QtGui
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    got = []
    diagram.compartment_clicked.connect(got.append)
    pos = QtCore.QPointF(-50.0, -50.0)
    diagram.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert got == []


def test_a_right_click_does_not_select(diagram):
    from PyQt6 import QtCore, QtGui
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    got = []
    diagram.compartment_clicked.connect(got.append)
    pos = QtCore.QPointF(110.0, 160.0)
    diagram.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos, QtCore.Qt.MouseButton.RightButton,
        QtCore.Qt.MouseButton.RightButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert got == []


def test_clicking_the_middle_of_a_shape_selects_that_shape(diagram):
    """Two simpler rules are both wrong here. The largest containing box is the cytoplasm every
    time, since the boxes nest; the smallest picks whichever unrelated organelle happens to have a
    tight box over the click, because each group's box includes its text label and a labelled
    organelle is much wider than its drawing."""
    diagram.set_palette({c: (0.5, 0.5, 0.5) for c in CD.COMPARTMENT_SL}, "")
    r = diagram.renderer()
    box = r.boundsOnElement("SL0188")            # nucleolus, inside the nucleus, inside the cell
    target, size = diagram._target(), r.defaultSize()
    x = target.x() + box.center().x() * target.width() / size.width()
    y = target.y() + box.center().y() * target.height() / size.height()
    assert diagram.organelle_at(x, y) == "SL0188"


def test_a_diagram_with_no_drawing_paints_nothing_rather_than_crashing(qapp, tmp_path):
    """The paint handler runs on every repaint; with no artwork it has to return, not raise."""
    from PyQt6 import QtCore, QtGui
    d = CD.CellDiagram(path=str(tmp_path / "nothing.svg"))
    d.resize(50, 50)
    d.paintEvent(QtGui.QPaintEvent(QtCore.QRect(0, 0, 50, 50)))
    assert d._target().width() >= 0


def test_a_widget_with_no_size_yet_does_not_divide_by_zero(diagram):
    """Qt lays widgets out after they are built, so the first paint can arrive at zero size."""
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    diagram.resize(0, 0)
    assert diagram.organelle_at(5, 5) == ""
    diagram.resize(220, 320)


def test_clicking_a_shape_whose_classes_are_not_in_the_palette_still_selects_one(diagram):
    """A shape stands for its compartments whether or not the current palette has them -- the click
    names one and the window decides what to do with it, rather than the drawing going inert."""
    from PyQt6 import QtCore, QtGui
    diagram.set_palette({}, "")
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    diagram.colour_of = {}                      # a different category is showing
    r = diagram.renderer()
    box = r.boundsOnElement("SL0018")
    target, size = diagram._target(), r.defaultSize()
    pos = QtCore.QPointF(target.x() + box.center().x() * target.width() / size.width(),
                         target.y() + box.center().y() * target.height() / size.height())
    got = []
    diagram.compartment_clicked.connect(got.append)
    diagram.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert got == ["apicoplast"]


def test_a_drawing_with_no_size_of_its_own_fills_the_widget(qapp, tmp_path):
    """An SVG that declares no dimensions renders into whatever space there is, rather than being
    scaled by a division by zero."""
    path = tmp_path / "flat.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="0" height="0">'
                    '<g id="SL0018"><rect width="0" height="0" fill="#000"/></g></svg>')
    d = CD.CellDiagram(path=str(path))
    d.resize(80, 200)                 # the widget has a minimum height of its own
    assert d._target().width() == d.width() and d._target().height() == d.height()
    assert d.organelle_at(10, 10) == ""
