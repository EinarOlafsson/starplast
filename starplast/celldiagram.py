#!/usr/bin/env python3
"""An apicomplexan cell, colored from the same palette as the map.

A list of 27 compartment names is a legend; a parasite with its rhoptries filled in the color the
rhoptries have on the map is a picture of where the genes are. This draws the second from the first,
at draw time, from `Window.colour_of` -- never from a second palette, because two palettes are two
claims about what a color means and one of them will drift.

The artwork is UniProt's subcellular-location diagram (`data/icons/Apicomplexa_cells.svg`): 59
organelles, each a `<g>` carrying a **UniProt SL identifier**. So the drawing already speaks a
controlled vocabulary and the work here is a mapping from hyperLOPIT compartment names onto SL
codes -- `COMPARTMENT_SL` -- rather than any renaming of the artwork.

## Three decisions this module implements, all of them forced by the data

**Several hyperLOPIT classes share one organelle.** Both rhoptry classes, both nucleus classes and
the three plasma-membrane classes map to a single shape, which cannot carry three colors at once.
The shape takes the color of whichever class is SELECTED, and says which one it is showing. With
nothing selected it is neutral: filling it with whichever class sorts first would be a claim nobody
made. Clicking it cycles through the classes that share it, saying each time which is now selected.

**Some compartments have no organelle in the drawing.** The proteasome classes and `apical 1/2` are
not in the artwork, and cytosol has no distinct shape worth clicking. They are listed beside the
diagram rather than silently absent -- between them they are a large part of this proteome, and a
compartment that vanishes from the legend reads as one that does not exist.

**Absence stays gray.** `unassigned` is 4,313 genes. It is never filled with a compartment color
and never looks like a measurement.
"""
from __future__ import annotations

import os
import re

from PyQt6 import QtCore, QtGui, QtSvg, QtWidgets

from . import paths

#: hyperLOPIT compartment -> the UniProt SL group in the artwork. Confirmed against the file rather
#: than assumed: every code here is present in it, and `missing_from_drawing` reports the classes
#: that are not, instead of dropping them.
COMPARTMENT_SL = {
    "rhoptries 1": "SL0233", "rhoptries 2": "SL0233",
    "micronemes": "SL0163",
    # The artwork's granule shape is UniProt's generic "Cytoplasmic granule" (SL0281). Mapped here
    # because it IS the granule in this drawing and dense granules are the compartment the whole
    # annotation workflow is aimed at -- but the substitution is stated rather than quietly made:
    # the shape is a generic granule, not a dense granule drawn from an apicomplexan.
    "dense granules": "SL0281",
    "apicoplast": "SL0018",
    "IMC": "SL0362",
    "mitochondrion - soluble": "SL0173", "mitochondrion - membranes": "SL0171",
    "Golgi": "SL0132",
    "ER": "SL0095", "ER 2": "SL0095",
    "nucleus - chromatin": "SL0191", "nucleus - non-chromatin": "SL0191",
    "nucleolus": "SL0188",
    "cytosol": "SL0091",
    "PM - integral": "SL0039", "PM - peripheral 1": "SL0039", "PM - peripheral 2": "SL0039",
    "tubulin cytoskeleton": "SL0090",
}

#: Deliberately NOT mapped, with the reason, because a wrong mapping is worse than a stated gap:
#:
#:   40S / 60S ribosome     no ribosome in the drawing at all
#:   19S / 20S proteasome   likewise
#:   apical 1 / apical 2    no conoid or apical complex shape
#:   endomembrane vesicles  the drawing has COPI and COPII vesicles, which are specific organelles;
#:                          the hyperLOPIT class is a mixed one and colouring COPI for it would be a
#:                          claim the data does not make
#:   unassigned             not a compartment
#:
#: These are reported by `missing_from_drawing` and named beside the diagram.
UNMAPPED_NOTE = ("no organelle in this drawing — the artwork has no ribosome, proteasome or apical "
                 "complex, and its vesicles are specifically COPI and COPII")

#: Fill for a shared shape when nothing is selected, and for any organelle whose class is not the
#: current one. Neutral on purpose -- see the module docstring.
NEUTRAL = (0.62, 0.62, 0.66)


def icon_path() -> str:
    """Where the artwork lives, resolved the way every other data file is."""
    return os.path.join(paths.data_dir(), "icons", "Apicomplexa_cells.svg")


def available(path: str = "") -> bool:
    """Whether the drawing is there. It ships with the package; a missing file hides the diagram."""
    return os.path.exists(path or icon_path())


def groups_in(svg: str) -> set:
    """Every SL group the artwork actually contains.

    Read from the file rather than trusted from the table above, because a mapping to a shape that
    is not in the drawing is a silent no-op -- the compartment would simply never color anything
    and nothing would say why.
    """
    return set(re.findall(r'<g[^>]*id="(SL\d+)"', svg))


def missing_from_drawing(compartments, svg: str = "") -> list:
    """Compartments with no organelle to color, in the order given.

    Named rather than dropped. Between the proteasome classes, the apical classes and anything the
    artwork lacks, these are a large part of the proteome, and a compartment that vanishes from the
    legend reads as one that does not exist.
    """
    have = groups_in(svg) if svg else None
    out = []
    for c in compartments:
        sl = COMPARTMENT_SL.get(c)
        if sl is None or (have is not None and sl not in have):
            out.append(c)
    return out


def sharing(sl: str) -> list:
    """Every compartment drawn by one organelle, in the table's order.

    What the click cycle steps through, and the reason a shared shape has to say which class it is
    currently showing.
    """
    return [c for c, code in COMPARTMENT_SL.items() if code == sl]


def _hex(colour) -> str:
    r, g, b = (int(max(0.0, min(1.0, float(v))) * 255) for v in tuple(colour)[:3])
    return f"#{r:02x}{g:02x}{b:02x}"


#: Groups that are "outside the cell" and are drawn as a filled sheet across the whole canvas. Made
#: transparent so the panel's own background shows through: the diagram sits on a ground that changes
#: with the theme, and a white -- or grey -- card behind the parasite reads as an image that failed
#: to load.
CANVAS_GROUPS = ("SL0112", "SL0243")

#: The white sheet the artwork is drawn on. Removed rather than recoloured: the diagram sits in a
#: panel whose colour changes with the theme, and a white card behind a parasite on a dark ground
#: reads as an image that failed to load.
BACKGROUND_ID = "path_1_"


def transparent_ground(svg: str) -> str:
    """Drop the artwork's own background so the panel's shows through.

    Two things make the card: a white rectangle behind everything, and the groups that stand for the
    space OUTSIDE the cell, which are drawn as filled sheets across the whole canvas. Both go.
    """
    out = re.sub(r'<rect[^>]*id="%s"[^>]*/>' % BACKGROUND_ID, "", svg, count=1)
    # The surround: a single path that traces the whole canvas and then the cell outline, filled
    # even-odd so it paints everything OUTSIDE the parasite. It is the artwork's way of drawing "not
    # in the cell", and on a themed panel it is a white card with a parasite-shaped hole in it.
    # It carries no fill of its own -- it inherits one -- so the fill is ADDED rather than replaced.
    # Matched on the `M0,0` that starts at the canvas corner: nothing else in the drawing does.
    def blank(m):
        tag = re.sub(r'\sfill="[^"]*"', "", m.group(0))
        return tag[:5] + ' fill="none"' + tag[5:]

    out = re.sub(r'<path[^>]*\sd="M0,0[^"]*"[^>]*>', blank, out, count=1)
    for sl in CANVAS_GROUPS:
        pattern = re.compile(r'(<g[^>]*id="%s".*?</g>)' % sl, re.S)
        out = pattern.sub(lambda m: re.sub(r'fill\s*:\s*[^;"\']+', "fill:none",
                                           re.sub(r'fill="(?!none)[^"]*"', 'fill="none"',
                                                  m.group(0))),
                          out, count=1)
    return out


def view_box(svg: str):
    """The artwork's own (x, y, width, height), or None. See `CellDiagram.viewbox`."""
    m = re.search(r'viewBox="([\d.\s-]+)"', svg or "")
    if not m:
        return None
    parts = [float(v) for v in m.group(1).split()]
    return tuple(parts) if len(parts) == 4 else None


def portrait(svg: str) -> str:
    """Turn the cell upright, apical end -- the rhoptries -- at the top.

    The artwork is drawn lying down, apex to the left, which is the orientation of a textbook figure
    and the wrong one for a tall panel: laid out flat in a 260-pixel column it is a strip a
    centimetre high. Rotated a quarter turn it fills the column, and it matches how the parasite is
    drawn in every paper about invasion -- apex first, because that is the end that goes in.

    Done by wrapping the drawing in a rotation and swapping the viewBox rather than by editing
    coordinates: the artwork has 59 groups and several thousand paths, and touching them would break
    the SL ids that everything else here depends on.
    """
    m = re.search(r'viewBox="([\d.\s-]+)"', svg)
    if not m:
        return svg
    x, y, w, h = (float(v) for v in m.group(1).split())
    # A quarter turn CLOCKWISE, which is the one that sends the left edge -- the apex -- to the top.
    # Anticlockwise puts the nucleus above the rhoptries, i.e. the parasite on its head, and a test
    # asserts the rhoptries end up higher than the nucleus rather than trusting the sign.
    inner = re.sub(r"^.*?<svg[^>]*>", "", svg, count=1, flags=re.S)
    inner = re.sub(r"</svg>\s*$", "", inner, flags=re.S)
    head = svg[:svg.index(">", svg.index("<svg")) + 1]
    head = head.replace(m.group(0), f'viewBox="0 0 {h:g} {w:g}"')
    return (f'{head}<g transform="rotate(90) translate({-x:g} {-(y + h):g})">'
            f'{inner}</g></svg>')


def recolour(svg: str, fills: dict, neutral=NEUTRAL) -> str:
    """Set the fill of each SL group, returning the modified SVG text.

    Done on the text at draw time rather than by shipping a recolored copy: there are four themes
    and the categorical palette changes with each, so any stored copy would be wrong for three of
    them and would drift from the map the moment a palette changed.

    The fill is applied to the group AND its descendants' inline styles, because the artwork sets
    `fill` on the individual paths -- a fill on the group alone is overridden by every path in it and
    the diagram would come back gray while every test on the returned string passed.
    """
    def paint(match, colour: str) -> str:
        body = match.group(0)
        body = re.sub(r'fill\s*:\s*[^;"\']+', f"fill:{colour}", body)
        body = re.sub(r'fill="(?!none)[^"]*"', f'fill="{colour}"', body)
        return body

    out = svg
    # Only the organelles a compartment names. Painting every group -- including the cell body and
    # the space around it -- filled the whole drawing with the neutral grey and produced a grey card
    # with a parasite on it, which is the opposite of taking the container's background.
    for sl in fills:
        colour = _hex(fills.get(sl, neutral))
        # The group and everything inside it, matched non-greedily up to its closing tag.
        pattern = re.compile(r'(<g[^>]*id="%s".*?</g>)' % sl, re.S)
        out = pattern.sub(lambda m: paint(m, colour), out, count=1)
    return out


class CellDiagram(QtWidgets.QWidget):
    """The parasite, filled from the map's own palette, clickable in both directions."""

    #: A compartment the user clicked in the drawing. The window selects it in the list, which is
    #: the other half of "both directions, or the diagram is decoration".
    compartment_clicked = QtCore.pyqtSignal(str)

    def __init__(self, parent=None, path: str = ""):
        super().__init__(parent)
        self.path = path or icon_path()
        raw = open(self.path, encoding="utf8").read() if available(self.path) else ""
        # The artwork's own coordinate system, kept because `boundsOnElement` reports positions in
        # it -- BEFORE the rotation that stands the cell upright. A click has to be mapped back
        # through that rotation or every organelle is hit-tested against the wrong place, which is
        # exactly the kind of "works, but selects the neighbour" bug this project keeps finding.
        self.viewbox = view_box(raw)
        self.svg = portrait(transparent_ground(raw)) if raw else ""
        self.groups = groups_in(self.svg)
        self.colour_of: dict = {}
        self.selected = ""
        self._renderer = None
        self.setMinimumHeight(150)
        self.setToolTip(
            "The same colors as the map, on the organelle each compartment names. Where several "
            "classes share one shape -- both rhoptry classes, both nucleus classes, the three "
            "plasma-membrane classes -- the shape shows the one that is selected and says so; click "
            "it again to step to the next. Gray means nothing is selected for that shape, or that "
            "the class is unassigned, which is not a compartment.")

    # ------------------------------------------------------------------ state
    def set_palette(self, colour_of: dict, selected: str = ""):
        """Take the map's palette and the current selection, and redraw from them."""
        self.colour_of = dict(colour_of or {})
        self.selected = selected or ""
        self._renderer = None
        self.update()

    def fills(self) -> dict:
        """SL group -> color, from the map's palette and the current selection.

        A shape shared by several classes takes the selected one's color, and neutral when none of
        them is selected: filling it with whichever sorts first would be a claim nobody made.
        """
        out = {}
        for sl in self.groups:
            classes = [c for c in sharing(sl) if c in self.colour_of]
            if not classes:
                continue
            if len(classes) == 1:
                out[sl] = self.colour_of[classes[0]]
            elif self.selected in classes:
                out[sl] = self.colour_of[self.selected]
            else:
                out[sl] = NEUTRAL
        return out

    def showing(self) -> str:
        """Which class the shared shape is currently showing, as a sentence, or ""."""
        sl = COMPARTMENT_SL.get(self.selected)
        others = [c for c in sharing(sl or "") if c != self.selected]
        if not self.selected or not others:
            return ""
        return (f"the {self.selected} color is on a shape shared with "
                f"{', '.join(others)} — click it to step through them")

    def renderer(self):
        """The SVG renderer for the current colors, built once per palette change."""
        if self._renderer is None and self.svg:
            data = recolour(self.svg, self.fills()).encode("utf8")
            self._renderer = QtSvg.QSvgRenderer(QtCore.QByteArray(data))
        return self._renderer

    # ------------------------------------------------------------------ drawing and clicking
    def paintEvent(self, ev):
        """Draw the cell, scaled to fit and centerd."""
        r = self.renderer()
        if r is None or not r.isValid():
            return
        p = QtGui.QPainter(self)
        r.render(p, QtCore.QRectF(self._target()))
        p.end()

    def _target(self) -> QtCore.QRectF:
        """Where the drawing goes: as large as fits, keeping its aspect ratio."""
        r = self.renderer()
        size = r.defaultSize() if r is not None else QtCore.QSize(1, 1)
        if size.width() <= 0 or size.height() <= 0:
            return QtCore.QRectF(self.rect())
        scale = min(self.width() / size.width(), self.height() / size.height())
        w, h = size.width() * scale, size.height() * scale
        return QtCore.QRectF((self.width() - w) / 2, (self.height() - h) / 2, w, h)

    def local_to_widget(self, px: float, py: float):
        """A point in the artwork's own frame to widget pixels, through the rotation and the fit.

        Defined beside its inverse so the two cannot drift: an off-by-a-quarter-turn hit test is
        invisible until someone clicks an organelle and selects its neighbour.
        """
        r = self.renderer()
        if r is None or not r.isValid():
            return 0.0, 0.0
        target, size = self._target(), r.defaultSize()
        vx, vy, vw, vh = self.viewbox or (0.0, 0.0, size.width(), size.height())
        ux, uy = (vy + vh) - py, px - vx
        return (target.x() + ux * target.width() / max(size.width(), 1),
                target.y() + uy * target.height() / max(size.height(), 1))

    def widget_to_local(self, x: float, y: float):
        """Widget pixels back to the artwork's own frame -- the inverse of `local_to_widget`."""
        r = self.renderer()
        if r is None or not r.isValid():
            return 0.0, 0.0
        target, size = self._target(), r.defaultSize()
        if target.width() <= 0 or target.height() <= 0:
            return 0.0, 0.0
        ux = (x - target.x()) * size.width() / target.width()
        uy = (y - target.y()) * size.height() / target.height()
        vx, vy, vw, vh = self.viewbox or (0.0, 0.0, size.width(), size.height())
        return uy + vx, (vy + vh) - ux

    def organelle_at(self, x: float, y: float) -> str:
        """The SL group under a widget point, or "".

        Among the boxes containing the point, the one it sits most centrally in wins -- distance
        from the center relative to the box's own size. Two other rules were tried and both are
        wrong here: the largest match is the cytoplasm every time, since the boxes nest, and the
        SMALLEST match picks whichever unrelated organelle happens to have a tight box over the
        click, because each group's box includes its text label and a labeled organelle's box is
        much wider than its drawing. Clicking the middle of a shape should select that shape.
        """
        r = self.renderer()
        if r is None or not r.isValid():
            return ""
        target, size = self._target(), r.defaultSize()
        if target.width() <= 0 or size.width() <= 0:
            return ""
        # Element bounds are reported in the artwork's own frame, and the drawing on screen has been
        # rotated a quarter turn into portrait around it, so the click is mapped back first.
        sx, sy = self.widget_to_local(x, y)
        best, best_score = "", float("inf")
        for sl in self.groups:
            if not sharing(sl):
                continue                       # not a compartment this map knows about
            box = r.boundsOnElement(sl)
            if box.isEmpty() or not box.contains(QtCore.QPointF(sx, sy)):
                continue
            radius = max(box.width(), box.height()) / 2.0 or 1.0
            score = ((sx - box.center().x()) ** 2 + (sy - box.center().y()) ** 2) ** 0.5 / radius
            if score < best_score:
                best, best_score = sl, score
        return best

    def mouseReleaseEvent(self, ev):
        """Click an organelle to select its compartment; click again to step to the next.

        Cycling is what makes a shared shape usable: one shape stands for three plasma-membrane
        classes, and without a cycle two of them would be unreachable from the drawing.
        """
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        p = ev.position()
        sl = self.organelle_at(p.x(), p.y())
        if not sl:
            return
        # Non-empty by construction: `organelle_at` only returns groups that stand for at least one
        # compartment, so there is no "no classes" case to guard here.
        classes = [c for c in sharing(sl) if c in self.colour_of] or sharing(sl)
        i = classes.index(self.selected) + 1 if self.selected in classes else 0
        self.compartment_clicked.emit(classes[i % len(classes)])
