#!/usr/bin/env python3
"""The cell diagram: the map's palette, on a parasite, clickable in both directions.

The point of it is that a list of 27 names becomes a picture of where the genes are, and the risk of
it is that a picture is persuasive whether or not it is true. So these tests are about the joins: the
colors come from the map's own dict rather than a second palette, a shape shared by three classes
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

from PyQt6 import QtCore, QtWidgets  # noqa: E402

from starplast import celldiagram as CD  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def svg():
    return open(CD.icon_path(), encoding="utf8").read()


def _ink_pixel(diagram, sl):
    """A widget point that is actually ON the drawn organelle.

    Bounding boxes cannot be used for this: every group in the artwork contains hidden `<text>`
    holding UniProt's description, so the Golgi's box is 4,894 units wide in a 1,190-wide drawing
    and a click anywhere landed on whichever inflated box won. That is why hit-testing renders the
    shapes, and why a test of it has to aim at a rendered pixel.
    """
    from PyQt6 import QtGui
    mask = diagram.masks()[sl]
    pts = [(x, y) for x in range(mask.width()) for y in range(mask.height())
           if QtGui.QColor.fromRgba(mask.pixel(x, y)).alpha() > 200]
    assert pts, f"{sl} drew nothing"
    return pts[len(pts) // 2]


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
    simply never color anything and nothing would say why."""
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


# --------------------------------------------------------------------------- coloring
def test_the_fill_reaches_the_paths_not_only_the_group(svg):
    """The artwork sets `fill` on the individual paths, so a fill on the group alone is overridden
    by every path in it -- the diagram comes back grey while every assertion on the string passes."""
    out = CD.recolor(svg, {"SL0018": (1.0, 0.0, 0.0)})
    import re
    group = re.search(r'<g[^>]*id="SL0018".*?</g>', out, re.S).group(0)
    assert "#ff0000" in group
    assert not re.search(r'fill\s*:\s*(?!#ff0000)#[0-9a-f]{6}', group), "a path kept its own fill"


def test_with_nothing_selected_nothing_is_colored(diagram):
    """The drawing is grey until a compartment is chosen. Coloring a shared shape with whichever
    class sorts first would be a claim nobody made, and coloring ALL of them makes the diagram a
    second legend -- twenty-odd colors to read against twenty-odd names."""
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)}, "")
    assert diagram.fills() == {}


def test_a_shared_shape_takes_the_color_of_whichever_class_is_selected(diagram):
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 2")
    assert diagram.fills()["SL0233"] == (0.0, 1.0, 0.0)
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 1")
    assert diagram.fills()["SL0233"] == (1.0, 0.0, 0.0)


def test_a_shared_shape_says_which_class_it_is_showing(diagram):
    """A user coming back to the window cannot otherwise tell whether the rhoptry is colored for 1
    or for 2."""
    diagram.set_palette({"rhoptries 1": (1.0, 0.0, 0.0), "rhoptries 2": (0.0, 1.0, 0.0)},
                        "rhoptries 2")
    said = diagram.showing()
    assert "rhoptries 2" in said and "rhoptries 1" in said
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "apicoplast")
    assert diagram.showing() == "", "a shape of its own has nothing to disambiguate"


def test_only_the_selected_compartment_is_ever_colored(diagram):
    """One colored organelle at a time, in the color that compartment has in the list beside it."""
    palette = {"apicoplast": (0.0, 0.0, 1.0), "micronemes": (0.0, 1.0, 0.0)}
    diagram.set_palette(palette, "apicoplast")
    assert diagram.fills() == {"SL0018": (0.0, 0.0, 1.0)}
    diagram.set_palette(palette, "micronemes")
    assert diagram.fills() == {"SL0163": (0.0, 1.0, 0.0)}
    diagram.set_palette(palette, "")
    assert diagram.fills() == {}


def test_a_compartment_the_palette_does_not_carry_colors_nothing(diagram):
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "rhoptries 1")
    assert diagram.fills() == {}


def test_a_compartment_the_palette_does_not_have_colors_nothing(diagram):
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
    x, y = _ink_pixel(diagram, "SL0018")
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
        classes = [c for c in CD.sharing("SL0039") if c in diagram.color_of]
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
    x, y = _ink_pixel(diagram, "SL0188")         # nucleolus, inside the nucleus, inside the cell
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
    diagram.color_of = {}                      # a different category is showing
    pos = QtCore.QPointF(*_ink_pixel(diagram, "SL0018"))
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


def test_the_cell_stands_upright_with_the_apex_at_the_top(diagram):
    """Drawn lying down it is a strip a centimetre high in a 260-pixel column -- and apex-first is
    how the parasite is drawn in every paper about invasion, because that is the end that goes in."""
    r = diagram.renderer()
    assert r.defaultSize().height() > r.defaultSize().width(), "the drawing is still landscape"
    # The rhoptries are apical, the nucleus is basal: upright means the rhoptries are ABOVE it.
    rh = r.boundsOnElement("SL0233").center()
    nu = r.boundsOnElement("SL0191").center()
    _, rh_y = diagram.local_to_widget(rh.x(), rh.y())
    _, nu_y = diagram.local_to_widget(nu.x(), nu.y())
    assert rh_y < nu_y, "the apical end is not at the top"


def test_the_artwork_brings_no_background_of_its_own(diagram, svg):
    """A white card behind a parasite on a dark ground reads as an image that failed to load. Two
    things make that card: a white rectangle behind everything, and a single path that traces the
    canvas and then the cell outline, filled even-odd so it paints everything outside the parasite."""
    assert 'id="path_1_"' in svg, "the artwork's background rect is gone from the source file"
    assert 'id="path_1_"' not in diagram.svg
    i = diagram.svg.find("M0,0v")
    assert i > 0, "the surround path is not in the drawing any more"
    tag_start = diagram.svg.rfind("<path", 0, i)
    assert 'fill="none"' in diagram.svg[tag_start:i], "the surround is still painted"


def test_mapping_a_point_to_the_widget_and_back_returns_it(diagram):
    """The two directions are defined together so they cannot drift: an off-by-a-quarter-turn hit
    test is invisible until someone clicks an organelle and selects its neighbour."""
    diagram.set_palette({"apicoplast": (0.0, 0.0, 1.0)}, "")
    x, y = diagram.local_to_widget(300.0, 400.0)
    back = diagram.widget_to_local(x, y)
    assert abs(back[0] - 300.0) < 1e-6 and abs(back[1] - 400.0) < 1e-6


def test_mapping_without_a_drawing_is_the_origin_rather_than_a_crash(qapp, tmp_path):
    d = CD.CellDiagram(path=str(tmp_path / "nothing.svg"))
    assert d.local_to_widget(1, 2) == (0.0, 0.0)
    assert d.widget_to_local(1, 2) == (0.0, 0.0)


def test_mapping_back_from_a_widget_with_no_size_is_the_origin(qapp, tmp_path):
    """Qt lays widgets out after they are built, so a click can arrive before the widget has a
    size -- and dividing by that width is how a hit test becomes an exception."""
    path = tmp_path / "flat.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="0" height="0" '
                    'viewBox="0 0 0 0"><g id="SL0018"><rect width="0" height="0"/></g></svg>')
    d = CD.CellDiagram(path=str(path))
    d.resize(0, 0)
    assert d.widget_to_local(3, 4) == (0.0, 0.0)


def test_the_drawing_is_hollow_and_only_the_selection_has_color(qapp):
    """What the eye actually gets, measured on the painted widget rather than on the source: every
    attempt to neutralise this artwork by rewriting its fills, strokes and gradient stops left it
    rendering in full color anyway, with no error and nothing in the document to explain it. A
    shape with no fill has nothing to render in any color, which is what finally settled it."""
    from PyQt6 import QtGui
    d = CD.CellDiagram()
    d.resize(240, 380)
    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "")

    def saturated(widget):
        img = widget.grab().toImage()
        n = 0
        for x in range(0, img.width(), 3):
            for y in range(0, img.height(), 3):
                c = QtGui.QColor.fromRgba(img.pixel(x, y))
                if max(c.red(), c.green(), c.blue()) - min(c.red(), c.green(), c.blue()) > 40:
                    n += 1
        return n

    assert saturated(d) == 0, "something is colored with nothing selected"
    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "rhoptries 1")
    assert saturated(d) > 0, "the selected compartment was not colored"


def test_the_drawing_keeps_the_panel_s_background(qapp):
    """Transparent where the cell is not, so the diagram sits on whatever the theme paints."""
    from PyQt6 import QtGui
    d = CD.CellDiagram()
    d.resize(200, 320)
    d.set_palette({}, "")
    img = QtGui.QImage(200, 320, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0)
    p = QtGui.QPainter(img)
    d.render(p, QtCore.QPoint(), QtGui.QRegion(d.rect()),
             QtWidgets.QWidget.RenderFlag.DrawChildren)
    p.end()
    corner = QtGui.QColor.fromRgba(img.pixel(2, 2))
    assert corner.alpha() < 40, "the drawing brought a background of its own"


def test_the_cell_is_outlines_with_nothing_filled_behind_them(qapp):
    """Hollow, so the diagram sits on whatever is behind it and the one filled thing in it is
    unmistakably the selection."""
    d = CD.CellDiagram()
    assert 'fill="none"' in d.svg
    assert '#ffffff' in d.svg, "the outlines are not white"
    # The masks come from the SOLID version: a click belongs to the organelle it lands inside, not
    # only to the two pixels of its outline.
    assert d.svg_solid and d.svg_solid != d.svg
    d.resize(240, 380)
    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "rhoptries 1")
    mask = d.masks()["SL0233"]
    from PyQt6 import QtGui
    ink = sum(QtGui.QColor.fromRgba(mask.pixel(x, y)).alpha() > 200
              for x in range(mask.width()) for y in range(mask.height()))
    assert ink > 200, "the hit-test mask is only an outline"


# --------------------------------------------------------------------------- reading the markup
def test_a_group_the_drawing_does_not_have_is_absent_rather_than_a_guess():
    """`group_span` is the one function every other part of this module asks "where is this
    organelle" -- the masks, the isolation, the coloring. Returning a plausible span for a group
    that is not there would put one organelle's ink under another's name."""
    assert CD.group_span('<svg><g id="SL0018"><path/></g></svg>', "SL9999") is None


def test_an_unclosed_group_is_absent_rather_than_running_to_the_end_of_the_file():
    """Malformed markup: the opening tag is found and nothing closes it. Returning everything after
    it would hand the rest of the drawing back as this organelle."""
    assert CD.group_span('<svg><g id="SL0018"><path/>', "SL0018") is None
    assert CD.group_span('<svg><g id="SL0018"><g><path/></g>', "SL0018") is None


def test_recoloring_a_group_that_is_not_in_the_drawing_changes_nothing():
    """A palette can name a compartment the artwork has no organelle for -- 13 of the 27 do. Each
    one has to be a no-op on the markup rather than an exception in the paint handler."""
    svg = '<svg><g id="SL0018"><path fill="#123456"/></g></svg>'
    assert CD.recolor(svg, {"SL9999": (1.0, 0.0, 0.0)}) == svg
    assert "#ff0000" in CD.recolor(svg, {"SL0018": (1.0, 0.0, 0.0)})


def test_isolating_a_group_that_is_not_there_is_empty_rather_than_a_broken_drawing(qapp, tmp_path):
    """The isolated markup is fed to a renderer; a header with no shape in it renders as a blank
    image, which would become a mask covering nothing that still claims to be an organelle."""
    path = tmp_path / "one.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
                    '<g id="SL0018"><rect width="4" height="4" fill="#000"/></g></svg>')
    d = CD.CellDiagram(path=str(path))
    assert d._isolated("SL9999") == ""
    assert d._isolated("SL0018").startswith("<svg")


def test_an_organelle_that_cannot_be_isolated_gets_no_mask(qapp, tmp_path, monkeypatch):
    """A mask that cannot be built is left out, so `organelle_at` reports nothing there. The
    alternative -- an empty mask kept under its name -- is a click target that matches no pixel and
    a compartment that looks drawn but cannot be selected."""
    path = tmp_path / "two.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
                    '<g id="SL0018"><rect width="4" height="4" fill="#000"/></g></svg>')
    d = CD.CellDiagram(path=str(path))
    d.resize(40, 40)
    assert "SL0018" in d.masks()
    d._masks, d._masks_for = {}, None
    monkeypatch.setattr(d, "_isolated", lambda sl: "")
    assert d.masks() == {}


def test_an_organelle_whose_markup_will_not_render_gets_no_mask(qapp, tmp_path, monkeypatch):
    """The other way the same thing fails: markup that is produced but that Qt refuses. Silently
    keeping an unrenderable organelle would give it a mask of nothing."""
    path = tmp_path / "three.svg"
    path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
                    '<g id="SL0018"><rect width="4" height="4" fill="#000"/></g></svg>')
    d = CD.CellDiagram(path=str(path))
    d.resize(40, 40)
    d._masks, d._masks_for = {}, None
    monkeypatch.setattr(d, "_isolated", lambda sl: "<svg>not really svg")
    assert d.masks() == {}


# --------------------------------------------------------------------------- the artwork's credit
def test_the_credit_block_is_not_part_of_the_drawing(qapp):
    """The artwork carries the creator's name and two SIB logos in the corner of its canvas. Hollowed
    and turned upright they render as a solid grey rectangle floating beside the parasite -- a filled
    shape in a drawing whose whole point is that the one filled thing is the selection, and one that
    belongs to no organelle so clicking it does nothing.

    Measured as the longest run of fully opaque pixels across a row, with nothing selected. The
    outlines are thin and anti-aliased and produce none; the logo produced 19."""
    from PyQt6 import QtCore, QtGui
    d = CD.CellDiagram()
    d.resize(240, 380)

    def longest_solid_run(widget):
        img = QtGui.QImage(widget.width(), widget.height(),
                           QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        p = QtGui.QPainter(img)
        widget.render(p, QtCore.QPoint(), QtGui.QRegion(widget.rect()),
                      QtWidgets.QWidget.RenderFlag.DrawChildren)
        p.end()
        worst = 0
        for y in range(img.height()):
            run = 0
            for x in range(img.width()):
                run = run + 1 if QtGui.QColor.fromRgba(img.pixel(x, y)).alpha() > 200 else 0
                worst = max(worst, run)
        return worst

    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "")
    assert longest_solid_run(d) < 6, "something in the drawing is a solid block with nothing selected"
    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "rhoptries 1")
    assert longest_solid_run(d) > 10, "the selected compartment is not filled"

    assert CD.group_span(d.svg, CD.CREDIT_GROUP) is None


def test_the_credit_survives_being_taken_out_of_the_picture(qapp):
    """The artwork is CC BY: taking the logo block out of the drawing without carrying its
    attribution somewhere a person can read would be a licence breach, not a tidy-up. Read out of
    the file rather than written down here, so it cannot quietly stop matching what it credits."""
    d = CD.CellDiagram()
    assert d.credit["creator"] and "creativecommons.org" in d.credit["license"]
    assert d.credit["creator"] in d.toolTip()
    assert d.credit["license"] in d.toolTip()


def test_an_artwork_with_no_credit_block_is_left_alone(qapp, tmp_path):
    """Not every drawing carries one, and removing a group that is not there must not remove
    something else."""
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
           '<g id="SL0018"><rect width="4" height="4" fill="#000"/></g></svg>')
    assert CD.drop_credit(svg) == svg
    assert CD.credit(svg) == {"creator": "", "license": ""}
    path = tmp_path / "plain.svg"
    path.write_text(svg)
    d = CD.CellDiagram(path=str(path))
    assert "Cell drawing by" not in d.toolTip()
