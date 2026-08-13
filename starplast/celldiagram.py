#!/usr/bin/env python3
"""An apicomplexan cell, coloured from the same palette as the map.

A list of 27 compartment names is a legend; a parasite with its rhoptries filled in the colour the
rhoptries have on the map is a picture of where the genes are. This draws the second from the first,
at draw time, from `Window.colour_of` -- never from a second palette, because two palettes are two
claims about what a colour means and one of them will drift.

The artwork is UniProt's subcellular-location diagram (`data/icons/Apicomplexa_cells.svg`): 59
organelles, each a `<g>` carrying a **UniProt SL identifier**. So the drawing already speaks a
controlled vocabulary and the work here is a mapping from hyperLOPIT compartment names onto SL
codes -- `COMPARTMENT_SL` -- rather than any renaming of the artwork.

## Three decisions this module implements, all of them forced by the data

**Several hyperLOPIT classes share one organelle.** Both rhoptry classes, both nucleus classes and
the three plasma-membrane classes map to a single shape, which cannot carry three colours at once.
The shape takes the colour of whichever class is SELECTED, and says which one it is showing. With
nothing selected it is neutral: filling it with whichever class sorts first would be a claim nobody
made. Clicking it cycles through the classes that share it, saying each time which is now selected.

**Some compartments have no organelle in the drawing.** The proteasome classes and `apical 1/2` are
not in the artwork, and cytosol has no distinct shape worth clicking. They are listed beside the
diagram rather than silently absent -- between them they are a large part of this proteome, and a
compartment that vanishes from the legend reads as one that does not exist.

**Absence stays grey.** `unassigned` is 4,313 genes. It is never filled with a compartment colour
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
    is not in the drawing is a silent no-op -- the compartment would simply never colour anything
    and nothing would say why.
    """
    return set(re.findall(r'<g[^>]*id="(SL\d+)"', svg))


def missing_from_drawing(compartments, svg: str = "") -> list:
    """Compartments with no organelle to colour, in the order given.

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


def recolour(svg: str, fills: dict, neutral=NEUTRAL) -> str:
    """Set the fill of each SL group, returning the modified SVG text.

    Done on the text at draw time rather than by shipping a recoloured copy: there are four themes
    and the categorical palette changes with each, so any stored copy would be wrong for three of
    them and would drift from the map the moment a palette changed.

    The fill is applied to the group AND its descendants' inline styles, because the artwork sets
    `fill` on the individual paths -- a fill on the group alone is overridden by every path in it and
    the diagram would come back grey while every test on the returned string passed.
    """
    def paint(match, colour: str) -> str:
        body = match.group(0)
        body = re.sub(r'fill\s*:\s*[^;"\']+', f"fill:{colour}", body)
        body = re.sub(r'fill="(?!none)[^"]*"', f'fill="{colour}"', body)
        return body

    out = svg
    for sl in groups_in(svg):
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
        self.svg = open(self.path, encoding="utf8").read() if available(self.path) else ""
        self.groups = groups_in(self.svg)
        self.colour_of: dict = {}
        self.selected = ""
        self._renderer = None
        self.setMinimumHeight(150)
        self.setToolTip(
            "The same colours as the map, on the organelle each compartment names. Where several "
            "classes share one shape -- both rhoptry classes, both nucleus classes, the three "
            "plasma-membrane classes -- the shape shows the one that is selected and says so; click "
            "it again to step to the next. Grey means nothing is selected for that shape, or that "
            "the class is unassigned, which is not a compartment.")

    # ------------------------------------------------------------------ state
    def set_palette(self, colour_of: dict, selected: str = ""):
        """Take the map's palette and the current selection, and redraw from them."""
        self.colour_of = dict(colour_of or {})
        self.selected = selected or ""
        self._renderer = None
        self.update()

    def fills(self) -> dict:
        """SL group -> colour, from the map's palette and the current selection.

        A shape shared by several classes takes the selected one's colour, and neutral when none of
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
        return (f"the {self.selected} colour is on a shape shared with "
                f"{', '.join(others)} — click it to step through them")

    def renderer(self):
        """The SVG renderer for the current colours, built once per palette change."""
        if self._renderer is None and self.svg:
            data = recolour(self.svg, self.fills()).encode("utf8")
            self._renderer = QtSvg.QSvgRenderer(QtCore.QByteArray(data))
        return self._renderer

    # ------------------------------------------------------------------ drawing and clicking
    def paintEvent(self, ev):
        """Draw the cell, scaled to fit and centred."""
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

    def organelle_at(self, x: float, y: float) -> str:
        """The SL group under a widget point, or "".

        Among the boxes containing the point, the one it sits most centrally in wins -- distance
        from the centre relative to the box's own size. Two other rules were tried and both are
        wrong here: the largest match is the cytoplasm every time, since the boxes nest, and the
        SMALLEST match picks whichever unrelated organelle happens to have a tight box over the
        click, because each group's box includes its text label and a labelled organelle's box is
        much wider than its drawing. Clicking the middle of a shape should select that shape.
        """
        r = self.renderer()
        if r is None or not r.isValid():
            return ""
        target, size = self._target(), r.defaultSize()
        if target.width() <= 0 or size.width() <= 0:
            return ""
        sx = (x - target.x()) * size.width() / target.width()
        sy = (y - target.y()) * size.height() / target.height()
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
