#!/usr/bin/env python3
"""An apicomplexan cell, colored from the same palette as the map.

A list of 27 compartment names is a legend; a parasite with its rhoptries filled in the color the
rhoptries have on the map is a picture of where the genes are. This draws the second from the first,
at draw time, from `Window.color_of` -- never from a second palette, because two palettes are two
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

import numpy as np
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

#: Organelles that are drawn as containers rather than as shapes: they have structure inside them,
#: and a flat fill on top of that structure erases it. Their color goes BEHIND the outlines, so the
#: cristae, the nucleolus and the vesicles stay visible on a colored ground. Everything else -- a
#: rhoptry, a granule, the Golgi -- is a shape, where filling it is the whole message.
BACKGROUND_SL = frozenset({
    "SL0091",   # cytosol
    "SL0191",   # nucleus
    "SL0173",   # mitochondrion, soluble
    "SL0171",   # mitochondrion, membranes
    "SL0018",   # apicoplast
})

#: Deliberately NOT mapped, with the reason, because a wrong mapping is worse than a stated gap:
#:
#:   40S / 60S ribosome     no ribosome in the drawing at all
#:   19S / 20S proteasome   likewise
#:   apical 1 / apical 2    no conoid or apical complex shape
#:   endomembrane vesicles  the drawing has COPI and COPII vesicles, which are specific organelles;
#:                          the hyperLOPIT class is a mixed one and coloring COPI for it would be a
#:                          claim the data does not make
#:   unassigned             not a compartment
#:
#: These are reported by `missing_from_drawing` and named beside the diagram.
UNMAPPED_NOTE = ("no organelle in this drawing — the artwork has no ribosome, proteasome or apical "
                 "complex, and its vesicles are specifically COPI and COPII")

#: Fill for a shared shape when nothing is selected, and for any organelle whose class is not the
#: current one. Neutral on purpose -- see the module docstring.
NEUTRAL = (0.62, 0.62, 0.66)

#: The drawing with nothing selected: light grey fills, mid-grey lines, transparent ground. The
#: artwork ships in color -- a pink cytoplasm, a red-brown nucleus -- and those colors mean nothing
#: here. Worse, they compete with the one color that does: a cytoplasm that is always red says
#: "cytosol is selected" when nothing is. Neutral, and then exactly one compartment takes the color
#: it has in the list beside it.
FILL_GREY = "#d9d9d9"
LINE_GREY = "#8d8d8d"


def strip_metadata(svg: str) -> str:
    """Remove the artwork's hidden text and links.

    UniProt ships each organelle with a `<text>` description, a `<text>` name and an `<a>` link to
    its own page, all `visibility="hidden"` -- and some of them nested INSIDE a `<path>`, where the
    SVG content model does not allow them. Qt's renderer says so, three times per render:

        qt.svg: <input>:42:452: Could not add child element to parent element because the types
        are incorrect.

    It then draws the file correctly, so the warnings are noise -- but noise printed on every repaint
    of every organelle mask is a console nobody can read, and a real warning would be lost in it.
    None of this is drawn: it is metadata for a page that is not this one.
    """
    out = re.sub(r"<text\b[^>]*/>", "", svg)
    out = re.sub(r"<text\b.*?</text>", "", out, flags=re.S)
    out = re.sub(r"<a\b[^>]*/>", "", out)
    out = re.sub(r"<a\b.*?</a>", "", out, flags=re.S)
    # And `<g class="membranes"/>` sitting INSIDE a `<path>`, which is the last of the three: a path
    # takes no children at all. Only groups with no id are removed -- an empty group draws nothing,
    # but an id'd one is an organelle this program addresses by name.
    out = re.sub(r"<g(?![^>]*\bid=)[^>]*/>", "", out)
    return out


def neutralise(svg: str) -> str:
    """Strip the artwork's own colors down to grey fills and grey lines.

    Every fill becomes one light grey and every stroke one mid-grey, so the drawing reads on a dark
    ground and on a light one, and so that ANY color in it is the selection. `fill="none"` is left
    alone: it is not a color, it is the absence of one, and filling those shapes in would turn the
    cell's internal outlines into solid blocks.
    """
    def paint(m):
        attr, value = m.group(1), m.group(2)
        if value.strip().lower() in ("none", "transparent"):
            return m.group(0)
        color = FILL_GREY if attr == "fill" else LINE_GREY
        return f'{attr}{m.group(0)[len(attr)]}{color}' if False else (
            f'{attr}="{color}"' if '="' in m.group(0) else f'{attr}:{color}')

    out = re.sub(r'(fill|stroke)="([^"]*)"', paint, svg)
    out = re.sub(r'(fill|stroke)\s*:\s*([^;"\']+)', paint, out)
    # And the gradients. Most of this artwork's color is not in a fill attribute at all -- it is in
    # 770 gradient stops, and a shape filled with `url(#SVGID_7_)` keeps its pink however many fill
    # attributes have been rewritten. Greying the stops is what actually makes the cell grey.
    out = re.sub(r'stop-color="[^"]*"', f'stop-color="{FILL_GREY}"', out)
    return out


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


def group_span(svg: str, sl: str):
    """(start, end) of one SL group's markup, or None -- counting nested `<g>` as it goes.

    A non-greedy `<g id="SL...".*?</g>` ends at the FIRST inner closing tag, and these groups nest:
    the rhoptry holds its membrane, the Golgi holds four sub-compartments. Matched that way, most
    organelles came back as a fragment -- which rendered as nothing, so eleven of the fourteen
    hit-test masks were silently empty and the coloring painted part of a shape.
    """
    m = re.search(r'<g[^>]*id="%s"[^>]*>' % re.escape(sl), svg)
    if not m:
        return None
    depth, i = 1, m.end()
    # A self-closing `<g/>` opens and closes at once. Counted as an opening tag it leaves the depth
    # permanently short, and six of the fourteen organelles -- the nucleus, the cytosol, the plasma
    # membrane among them -- came back as "no such group" and could neither be colored nor clicked.
    step = re.compile(r"<g\b[^>]*/>|<g\b|</g>")
    while depth and i < len(svg):
        nxt = step.search(svg[i:])
        if not nxt:
            return None
        token = nxt.group(0)
        if token.endswith("/>"):
            pass
        elif token == "</g>":
            depth -= 1
        else:
            depth += 1
        i += nxt.end()
    return (m.start(), i) if depth == 0 else None


def sharing(sl: str) -> list:
    """Every compartment drawn by one organelle, in the table's order.

    What the click cycle steps through, and the reason a shared shape has to say which class it is
    currently showing.
    """
    return [c for c, code in COMPARTMENT_SL.items() if code == sl]


def _hex(color) -> str:
    r, g, b = (int(max(0.0, min(1.0, float(v))) * 255) for v in tuple(color)[:3])
    return f"#{r:02x}{g:02x}{b:02x}"


#: Groups that are "outside the cell" and are drawn as a filled sheet across the whole canvas. Made
#: transparent so the panel's own background shows through: the diagram sits on a ground that changes
#: with the theme, and a white -- or grey -- card behind the parasite reads as an image that failed
#: to load.
CANVAS_GROUPS = ("SL0112", "SL0243")

#: The white sheet the artwork is drawn on. Removed rather than recolored: the diagram sits in a
#: panel whose color changes with the theme, and a white card behind a parasite on a dark ground
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
        tag = re.sub(r'\sstroke="[^"]*"', "", tag)
        # Stroke as well as fill: hollowing the drawing gives every shape an outline, and an outline
        # around the whole canvas is a frame nobody asked for.
        return tag[:5] + ' fill="none" stroke="none" stroke-width="0"' + tag[5:]

    out = re.sub(r'<path[^>]*\sd="M0,0[^"]*"[^>]*>', blank, out, count=1)
    for sl in CANVAS_GROUPS:
        pattern = re.compile(r'(<g[^>]*id="%s".*?</g>)' % sl, re.S)
        out = pattern.sub(lambda m: re.sub(r'fill\s*:\s*[^;"\']+', "fill:none",
                                           re.sub(r'fill="(?!none)[^"]*"', 'fill="none"',
                                                  m.group(0))),
                          out, count=1)
    return out


#: The artwork's own credit block: the creator's name and the SIB logos, sitting in the corner of the
#: canvas. It has to come out of the DRAWING -- hollowed and rotated it renders as a solid grey
#: rectangle beside the parasite, which reads as a compartment nobody can click -- and the credit it
#: carries has to go somewhere a person can actually read it. See `credit`.
CREDIT_GROUP = "sib_copyright"


def credit(svg: str) -> dict:
    """Who drew the cell and under what licence, read out of the artwork itself.

    Read rather than written down here: the file embeds `property="creator"` and a CC licence link,
    and a hand-copied attribution is one that silently stops matching the file it credits. The
    artwork is CC BY, so the attribution is a condition of using it, not a courtesy -- which is why
    dropping the logo block from the drawing has to be paired with showing this.
    """
    name = re.search(r'property="name"[^>]*>([^<]+)<', svg or "")
    lic = re.search(r'property="license"[^>]*href="([^"]+)"', svg or "")
    return {"creator": name.group(1).strip() if name else "",
            "license": lic.group(1).strip() if lic else ""}


def drop_credit(svg: str) -> str:
    """Remove the credit block from the drawing (the credit itself travels via `credit`)."""
    span = group_span(svg, CREDIT_GROUP)
    return svg if span is None else svg[:span[0]] + svg[span[1]:]


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


def recolor(svg: str, fills: dict, neutral=NEUTRAL) -> str:
    """Set the fill of each SL group, returning the modified SVG text.

    Done on the text at draw time rather than by shipping a recolored copy: there are four themes
    and the categorical palette changes with each, so any stored copy would be wrong for three of
    them and would drift from the map the moment a palette changed.

    The fill is applied to the group AND its descendants' inline styles, because the artwork sets
    `fill` on the individual paths -- a fill on the group alone is overridden by every path in it and
    the diagram would come back gray while every test on the returned string passed.
    """
    out = svg
    # Only the organelles a compartment names -- in practice the one that is selected. Painting every
    # group, including the cell body and the space around it, filled the whole drawing with grey.
    for sl in fills:
        span = group_span(out, sl)
        if span is None:
            continue
        a, b = span
        chunk = re.sub(r'fill\s*:\s*(?!none)[^;"\']+', f"fill:{_hex(fills[sl])}", out[a:b])
        chunk = re.sub(r'fill="(?!none)[^"]*"', f'fill="{_hex(fills[sl])}"', chunk)
        out = out[:a] + chunk + out[b:]
    return out


def outline_only(svg: str, stroke: str = "#ffffff", width: float = 1.4) -> str:
    """Hollow every shape out and draw it in white line.

    A cell drawn as outlines sits on whatever is behind it -- the panel, the theme, anything -- and
    the one filled thing in it is then unmistakably the selection. Filled shapes, however pale, are
    a card the diagram sits on.

    Shapes that carried no stroke get one, because a shape whose fill is removed and which has no
    outline is simply gone: the nucleus and the cell body are both drawn as fills in this artwork.
    """
    def hollow(m):
        tag = m.group(0)
        if 'stroke="none"' in tag:
            return tag                      # deliberately blanked -- the canvas surround
        tag = re.sub(r'\sfill="[^"]*"', ' fill="none"', tag)
        if 'fill="none"' not in tag:
            tag = tag[: len(m.group(1)) + 1] + ' fill="none"' + tag[len(m.group(1)) + 1:]
        tag = re.sub(r'\sstroke="[^"]*"', f' stroke="{stroke}"', tag)
        if f'stroke="{stroke}"' not in tag:
            tag = tag[: len(m.group(1)) + 1] + f' stroke="{stroke}"' + tag[len(m.group(1)) + 1:]
        if "stroke-width" not in tag:
            tag = tag[: len(m.group(1)) + 1] + f' stroke-width="{width}"' + tag[len(m.group(1)) + 1:]
        return tag

    return re.sub(r"<(path|rect|circle|ellipse|polygon|polyline|line)\b[^>]*>", hollow, svg)


def _tinted(mask: QtGui.QImage, color) -> QtGui.QImage:
    """The shape of `mask`, painted flat in `color` -- the one colored thing in the diagram."""
    out = QtGui.QImage(mask)
    out = out.convertToFormat(QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    p = QtGui.QPainter(out)
    p.setCompositionMode(QtGui.QPainter.CompositionMode.CompositionMode_SourceIn)
    p.fillRect(out.rect(), QtGui.QColor.fromRgbF(*tuple(color)[:3]))
    p.end()
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
        # The credit is IN the metadata that is about to be stripped -- the creator's name is a
        # `<text property="name">` and the licence an `<a property="license">` -- so it is read
        # first. Stripping it and then looking for it is how the attribution silently disappeared.
        self.credit = credit(raw)
        raw = neutralise(strip_metadata(raw)) if raw else raw
        # The artwork's own coordinate system, kept because `boundsOnElement` reports positions in
        # it -- BEFORE the rotation that stands the cell upright. A click has to be mapped back
        # through that rotation or every organelle is hit-tested against the wrong place, which is
        # exactly the kind of "works, but selects the neighbour" bug this project keeps finding.
        self.viewbox = view_box(raw)
        base = portrait(transparent_ground(drop_credit(raw))) if raw else ""
        # Two versions of the same drawing. The DISPLAYED one is hollow, white line on nothing; the
        # SOLID one is what the hit-test masks are rendered from, because a click belongs to the
        # organelle it lands inside, not only to the two pixels of its outline.
        self.svg_solid = base
        self.svg = outline_only(base) if base else ""
        self.groups = groups_in(self.svg)
        self.color_of: dict = {}
        self.selected = ""
        self._renderer = None
        # Pixel masks for hit-testing, and their ink, cached per widget size -- see `masks`.
        self._masks, self._masks_for, self._ink = {}, None, {}
        self.setMinimumHeight(150)
        self.setToolTip(
            "The same colors as the map, on the organelle each compartment names. Where several "
            "classes share one shape -- both rhoptry classes, both nucleus classes, the three "
            "plasma-membrane classes -- the shape shows the one that is selected and says so; click "
            "it again to step to the next. Gray means nothing is selected for that shape, or that "
            "the class is unassigned, which is not a compartment."
            + (f"\n\nCell drawing by {self.credit['creator']} (SwissBioPics, SIB)"
               + (f", {self.credit['license']}" if self.credit["license"] else "")
               if self.credit["creator"] else ""))

    # ------------------------------------------------------------------ state
    def set_palette(self, color_of: dict, selected: str = ""):
        """Take the map's palette and the current selection, and redraw from them."""
        self.color_of = dict(color_of or {})
        self.selected = selected or ""
        self._renderer = None
        self.update()

    def fills(self) -> dict:
        """SL group -> color: the SELECTED compartment, and nothing else.

        One colored organelle at a time, in the color that compartment has in the list and on the
        map. Coloring every organelle at once makes the diagram a second legend -- twenty-odd
        colors to read against twenty-odd names -- when what a person wants to know is where the
        thing they just clicked is. Everything else is grey, which is also what the map shows: the
        points of one compartment against a grey field.
        """
        sl = COMPARTMENT_SL.get(self.selected)
        if not sl or sl not in self.groups or self.selected not in self.color_of:
            return {}
        return {sl: self.color_of[self.selected]}

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
            data = recolor(self.svg, self.fills()).encode("utf8")
            self._renderer = QtSvg.QSvgRenderer(QtCore.QByteArray(data))
        return self._renderer

    # ------------------------------------------------------------------ drawing and clicking
    def paintEvent(self, ev):
        """Draw the cell grey, then tint the selected organelle.

        The drawing is hollow -- white outlines on nothing -- so it sits on whatever is behind it and
        the one FILLED thing in it is unmistakably the selection, in the color that compartment has
        in the list beside it.

        Hollowing it is also what finally made it neutral. Every fill, stroke and gradient stop in
        this file can be set to grey and it still renders pink: 145 elements carry a `colored` class
        and eleven gradients cross-reference each other, and nothing in the document accounts for it.
        A shape with no fill has nothing to render in any color.
        """
        r = self.renderer()
        if r is None or not r.isValid():
            return
        w, h = max(self.width(), 1), max(self.height(), 1)
        img = QtGui.QImage(w, h, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        # Transparent: the diagram takes the panel's background, whatever the theme makes it.
        img.fill(0)
        q = QtGui.QPainter(img)
        r.render(q, self._target())
        q.end()
        p = QtGui.QPainter(self)
        fills = self.fills()
        # The big compartments are tinted UNDER the drawing, the small ones over it. A rhoptry is a
        # shape and filling it says everything there is to say; the cytosol, the nucleus, the
        # mitochondrion and the apicoplast are containers with structure drawn inside them, and a
        # flat fill on top erases the cristae, the nucleolus, the vesicles -- the picture stops being
        # a cell and becomes a coloured blob where the cell was.
        for sl, color in fills.items():
            mask = self.masks().get(sl)
            if mask is not None and sl in BACKGROUND_SL:
                p.drawImage(0, 0, _tinted(mask, color))
        p.drawImage(0, 0, img)
        for sl, color in fills.items():
            mask = self.masks().get(sl)
            if mask is not None and sl not in BACKGROUND_SL:
                p.drawImage(0, 0, _tinted(mask, color))
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

    def _isolated(self, sl: str) -> str:
        """The drawing with one organelle in it and nothing else, positioned as it is on screen.

        The group's own markup, put back inside the same header and the same rotation, so it lands
        exactly where it lands in the full drawing.
        """
        source = self.svg_solid or self.svg
        span = group_span(source, sl)
        if span is None:
            return ""
        head = source[:source.index(">", source.index("<svg")) + 1]
        turn = re.search(r'<g transform="rotate\([^"]*\)">', source)
        return (head + (turn.group(0) if turn else "") + source[span[0]:span[1]]
                + ("</g>" if turn else "") + "</svg>")

    def masks(self):
        """Which pixels each organelle actually covers, at the current size.

        Hit-testing used bounding boxes, and the boxes are useless here: every group contains hidden
        `<text>` blocks holding UniProt's description of that compartment, so the Golgi's box is
        4,894 units wide in a 1,190-wide drawing. Clicking a rhoptry landed on whichever inflated box
        happened to win -- the proteasome, the plasma membrane, anything.

        So the shapes are rendered. Each organelle is drawn alone into a small image and its ink
        recorded; the hidden text does not render, so the mask is exactly what a person sees. Built
        once per size and cached, because it costs one render per organelle.
        """
        key = (self.width(), self.height())
        if self._masks_for == key and self._masks:
            return self._masks
        out = {}
        for sl in self.groups:
            if not sharing(sl):
                continue                       # not a compartment this map knows about
            markup = self._isolated(sl)
            if not markup:
                continue
            r = QtSvg.QSvgRenderer(QtCore.QByteArray(markup.encode("utf8")))
            if not r.isValid():
                continue
            img = QtGui.QImage(max(self.width(), 1), max(self.height(), 1),
                               QtGui.QImage.Format.Format_ARGB32_Premultiplied)
            img.fill(0)
            p = QtGui.QPainter(img)
            r.render(p, self._target())
            p.end()
            out[sl] = img
        self._masks, self._masks_for = out, key
        return out

    def organelle_at(self, x: float, y: float) -> str:
        """The organelle under a widget point, or "".

        Whichever drawn shape actually contains the pixel; where several do -- the nucleolus sits
        inside the nucleus -- the SMALLEST wins, because the small one is the thing being aimed at.
        """
        x, y = int(x), int(y)
        if not (0 <= x < self.width() and 0 <= y < self.height()):
            return ""
        best, best_ink = "", float("inf")
        for sl, img in self.masks().items():
            if QtGui.QColor.fromRgba(img.pixel(x, y)).alpha() < 40:
                continue
            ink = self._ink.get(sl)
            if ink is None:
                # Counted on a coarse grid: this is a tie-break between shapes, not a measurement.
                ink = sum(QtGui.QColor.fromRgba(img.pixel(a, b)).alpha() > 40
                          for a in range(0, img.width(), 3) for b in range(0, img.height(), 3))
                self._ink[sl] = ink
            if ink < best_ink:
                best, best_ink = sl, ink
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
        classes = [c for c in sharing(sl) if c in self.color_of] or sharing(sl)
        i = classes.index(self.selected) + 1 if self.selected in classes else 0
        self.compartment_clicked.emit(classes[i % len(classes)])
