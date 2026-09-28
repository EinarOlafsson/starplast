"""The strategy card: what a strategy answers and how well it works, before any prose.

Selecting a strategy in the Strategies tab shows this card first. It is built to be read in a few
seconds and to look the same for every strategy, so two cards can be compared by eye:

    header      number and name, the method as a chip, the calibration grade as a coloured badge,
                and one plain line saying what it answers;
    four bars   the same four positions on every card (`scorecard.headline`): Better than chance,
                Reach, and the two metrics of the strategy's task a biologist asks about first --
                each with its 95% interval and a tick at its chance level, and a plain sentence on
                hover;
    buttons     Run, Test and Details (the Guide, the Settings and the Results tables, one click
                away);
    about       four collapsed lines -- what it does, how it is evaluated, what failure looks like,
                what success looks like (`strategy_explainers`) -- each opening to its full text;
    examples    one real failure and one real success, side by side as mini scorecards, and what
                the successful run added about genes nobody has labeled.

The bars are painted, not tables, so the eye reads lengths instead of numbers. Colours come from the
theme the application is showing, never a literal at a call site.
"""
from __future__ import annotations

import math

from PyQt6 import QtCore, QtGui, QtWidgets

from . import scorecard as SC
from . import strategy_explainers as EX
from . import theme as TH

#: The calibration grade's colour role, and a plain sentence for its tooltip.
GRADES = {
    "reliable": ("success", "Reliable: across many held-out tests at its default settings, its "
                            "skill is clearly above zero and it passes most of them."),
    "works when tuned": ("warning", "Works when tuned: weak at its default settings, reliable at "
                                    "the setting calibration chose (Details, Settings, Use tuned "
                                    "settings)."),
    "weak": ("warning", "Weak: some skill on average, but not enough, or not often enough, to rely "
                        "on at any setting tested."),
    "no skill": ("error", "No skill: across its held-out tests it did no better than the same "
                          "procedure on shuffled data."),
    "untestable": ("fg_dim", "Untestable: too few of its tests could be scored on this organism's "
                             "table to judge it."),
}
VERDICTS = {"PASS": "success", "FAIL": "error", "INCONCLUSIVE": "fg_dim", "NOT RUN": "fg_dim"}


def colors() -> dict:
    """The palette of the theme the application is showing now."""
    app = QtWidgets.QApplication.instance()
    name = str(app.property("starplastTheme") or "dark") if app is not None else "dark"
    return TH.palette_for(name)


def _font(widget, delta: float = 0.0, bold: bool = False) -> QtGui.QFont:
    """The widget's font, `delta` points (or the same in pixels) larger, optionally bold."""
    f = QtGui.QFont(widget.font())
    if f.pixelSize() > 0:                       # the application stylesheet sizes in pixels
        f.setPixelSize(max(8, int(round(f.pixelSize() + delta * 1.33))))
    else:
        size = f.pointSizeF() if f.pointSizeF() > 0 else 10.0
        f.setPointSizeF(max(6.0, size + delta))
    f.setBold(bold)
    return f


# --------------------------------------------------------------------------- small painted parts
class Pill(QtWidgets.QWidget):
    """A rounded label: filled for a grade or verdict, outlined for a method chip."""

    def __init__(self, text: str = "", role: str = "accent", filled: bool = True, tip: str = "",
                 parent=None):
        super().__init__(parent)
        self.role, self.filled = role, filled
        self.setText(text)
        if tip:
            self.setToolTip(TH.tip(tip))
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                           QtWidgets.QSizePolicy.Policy.Fixed)

    def setText(self, text: str):
        """Change the label and resize to fit it."""
        self._text = str(text)
        self.updateGeometry()
        self.update()

    def text(self) -> str:
        """The label shown."""
        return self._text

    def sizeHint(self):
        """Wide enough for the label in the small bold face, plus its rounded ends."""
        fm = QtGui.QFontMetrics(_font(self, -1, True))
        return QtCore.QSize(fm.horizontalAdvance(self._text) + 18, fm.height() + 6)

    def paintEvent(self, _event):
        c = colors()
        colour = QtGui.QColor(c.get(self.role, self.role))
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        r = QtCore.QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        radius = r.height() / 2
        if self.filled:
            fill = QtGui.QColor(colour)
            fill.setAlphaF(0.22)
            p.setBrush(fill)
            p.setPen(QtGui.QPen(colour, 1.0))
        else:
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.setPen(QtGui.QPen(QtGui.QColor(c["border"]), 1.0))
        p.drawRoundedRect(r, radius, radius)
        p.setFont(_font(self, -1, self.filled))
        p.setPen(colour if self.filled else QtGui.QColor(c["fg_muted"]))
        p.drawText(r, int(QtCore.Qt.AlignmentFlag.AlignCenter), self._text)
        p.end()


class ScoreBar(QtWidgets.QWidget):
    """One headline bar: plain label and value, a painted bar with its interval and chance tick,
    and the technical name small underneath. `compact` drops the technical line (mini cards)."""

    def __init__(self, compact: bool = False, parent=None):
        super().__init__(parent)
        self.compact = compact
        self.bar: dict = {}
        self.setMinimumWidth(120)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Fixed)

    def set_bar(self, bar: dict):
        """Show one entry of `scorecard.headline_bars`."""
        self.bar = dict(bar or {})
        self.setToolTip(TH.tip(self.bar.get("reading", "")))
        self.update()

    def value_text(self) -> str:
        """The value as the bar prints it."""
        return SC.fmt_value(self.bar.get("scale", "unit"), self.bar.get("value"))

    def sizeHint(self):
        """Tall enough for the label line, the bar and (unless compact) the technical name."""
        fm = QtGui.QFontMetrics(_font(self, -1 if self.compact else 0))
        small = QtGui.QFontMetrics(_font(self, -2))
        h = fm.height() + 12 + (0 if self.compact else small.height() + 2)
        return QtCore.QSize(260, h + 2)

    def _x(self, track: QtCore.QRectF, value) -> float:
        t = SC.position(self.bar.get("scale", "unit"), value)
        return float("nan") if not math.isfinite(t) else track.left() + t * track.width()

    def paintEvent(self, _event):
        c = colors()
        b = self.bar
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        main = _font(self, -1 if self.compact else 0)
        small = _font(self, -2)
        fm, sm = QtGui.QFontMetrics(main), QtGui.QFontMetrics(small)
        w = self.width()
        # Line 1: the plain label, and the value with its interval on the right.
        value = self.value_text()
        low, high = b.get("low"), b.get("high")
        interval = ""
        if low is not None and high is not None and all(
                isinstance(v, (int, float)) and math.isfinite(v) for v in (low, high)):
            interval = (f"  {SC.fmt_value(b.get('scale', 'unit'), low)}"
                        f"-{SC.fmt_value(b.get('scale', 'unit'), high)}")
        p.setFont(small)
        iw = sm.horizontalAdvance(interval)
        p.setFont(main)
        vw = fm.horizontalAdvance(value)
        line = QtCore.QRectF(0, 0, w, fm.height())
        p.setPen(QtGui.QColor(c["fg"]))
        label = fm.elidedText(b.get("label", ""), QtCore.Qt.TextElideMode.ElideRight,
                              max(10, int(w - vw - iw - 8)))
        p.drawText(line, int(QtCore.Qt.AlignmentFlag.AlignLeft
                             | QtCore.Qt.AlignmentFlag.AlignVCenter), label)
        finite = value != "--"
        p.setPen(QtGui.QColor(c["fg"] if finite else c["fg_dim"]))
        p.drawText(QtCore.QRectF(0, 0, w - iw, fm.height()),
                   int(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter),
                   value)
        if interval:
            p.setFont(small)
            p.setPen(QtGui.QColor(c["fg_dim"]))
            p.drawText(QtCore.QRectF(0, 0, w, fm.height()),
                       int(QtCore.Qt.AlignmentFlag.AlignRight
                           | QtCore.Qt.AlignmentFlag.AlignVCenter), interval)
        # Line 2: the bar -- track, value, its interval as a whisker, chance as a marker above.
        y = fm.height() + 5
        thick = 4.0 if self.compact else 6.0
        track = QtCore.QRectF(3, y, w - 6, thick)
        p.setPen(QtCore.Qt.PenStyle.NoPen)
        p.setBrush(QtGui.QColor(c["surface_hi"]))
        p.drawRoundedRect(track, thick / 2, thick / 2)
        xv = self._x(track, b.get("value"))
        below = finite and isinstance(b.get("value"), (int, float)) and b["value"] < 0 \
            and b.get("scale") == "unit"
        if math.isfinite(xv) and xv > track.left() + 0.5:
            p.setBrush(QtGui.QColor(c["accent"]))
            p.drawRoundedRect(QtCore.QRectF(track.left(), y, xv - track.left(), thick),
                              thick / 2, thick / 2)
        if below:
            p.setBrush(QtGui.QColor(c["error"]))
            p.drawEllipse(QtCore.QPointF(track.left() + thick / 2, y + thick / 2), thick / 2,
                          thick / 2)
        xl, xh = self._x(track, low), self._x(track, high)
        if math.isfinite(xl) and math.isfinite(xh) and xh - xl > 1:
            whisker = QtGui.QColor(c["fg"])
            whisker.setAlphaF(0.75)
            p.setPen(QtGui.QPen(whisker, 1.3))
            mid = y + thick / 2
            p.drawLine(QtCore.QPointF(xl, mid), QtCore.QPointF(xh, mid))
            for x in (xl, xh):
                p.drawLine(QtCore.QPointF(x, y - 1), QtCore.QPointF(x, y + thick + 1))
        xc = self._x(track, b.get("chance"))
        if math.isfinite(xc):
            size = 3.0 if self.compact else 4.0
            p.setPen(QtCore.Qt.PenStyle.NoPen)
            p.setBrush(QtGui.QColor(c["fg_muted"]))
            p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(xc - size, y - size - 2),
                                           QtCore.QPointF(xc + size, y - size - 2),
                                           QtCore.QPointF(xc, y - 1)]))
        # Line 3: the technical name, small, and the chance level the marker stands at.
        if not self.compact:
            y3 = y + thick + 4
            p.setFont(small)
            p.setPen(QtGui.QColor(c["fg_dim"]))
            p.drawText(QtCore.QRectF(0, y3, w, sm.height()),
                       int(QtCore.Qt.AlignmentFlag.AlignLeft
                           | QtCore.Qt.AlignmentFlag.AlignVCenter), b.get("technical", ""))
            ch = b.get("chance")
            if isinstance(ch, (int, float)) and math.isfinite(ch):
                p.drawText(QtCore.QRectF(0, y3, w, sm.height()),
                           int(QtCore.Qt.AlignmentFlag.AlignRight
                               | QtCore.Qt.AlignmentFlag.AlignVCenter),
                           f"▾ chance {SC.fmt_value(b.get('scale', 'unit'), ch)}")
        p.end()


class Scorecard(QtWidgets.QWidget):
    """The four headline bars, stacked, in the same positions for every strategy."""

    def __init__(self, compact: bool = False, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4 if compact else 10)
        self.bars = [ScoreBar(compact) for _ in range(4)]
        for b in self.bars:
            lay.addWidget(b)

    def set_bars(self, bars: list):
        """Fill the four bars from `scorecard.headline_bars`."""
        for w, b in zip(self.bars, list(bars) + [{}] * 4):
            w.set_bar(b)


class ElidedLabel(QtWidgets.QLabel):
    """A one-line label that shortens its text with an ellipsis rather than wrapping."""

    def __init__(self, text: str = "", parent=None):
        super().__init__(parent)
        self._full = ""
        self.setText(text)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Fixed)

    def setText(self, text: str):
        """Keep the whole text; it is shortened only when painted."""
        self._full = str(text)
        super().setText(self._full)

    def full_text(self) -> str:
        """The whole text, before eliding."""
        return self._full

    def paintEvent(self, _event):
        p = QtGui.QPainter(self)
        p.setFont(self.font())
        p.setPen(self.palette().color(self.foregroundRole()))
        text = self.fontMetrics().elidedText(self._full, QtCore.Qt.TextElideMode.ElideRight,
                                             self.width())
        p.drawText(self.rect(), int(QtCore.Qt.AlignmentFlag.AlignLeft
                                    | QtCore.Qt.AlignmentFlag.AlignVCenter), text)
        p.end()


def _role(widget: QtWidgets.QWidget, role: str) -> QtWidgets.QWidget:
    """Name a widget by its role, which the card's own stylesheet (`card_style`) sizes and colours.

    A stylesheet, not a palette or a font: the application's stylesheet sets colour and size on
    every widget and beats both, so only a more specific stylesheet reaches these labels.
    """
    widget.setObjectName(f"card_{role}")
    return widget


def card_style(base_px: int = 13) -> str:
    """The card's stylesheet for the theme in force: sizes relative to the body text size."""
    c = colors()
    b = max(8, int(base_px))
    return f"""
    QLabel#card_name {{ font-size: {b + 4}px; font-weight: 600; color: {c['fg']}; }}
    QLabel#card_question {{ font-size: {b}px; color: {c['fg_muted']}; }}
    QLabel#card_section {{ font-size: {b - 3}px; font-weight: 600; color: {c['fg_dim']};
                           letter-spacing: 1px; padding-top: 4px; }}
    QLabel#card_title {{ font-size: {b - 1}px; font-weight: 600; color: {c['fg']}; }}
    QLabel#card_body {{ font-size: {b - 1}px; color: {c['fg_muted']}; }}
    QLabel#card_dim {{ font-size: {b - 2}px; color: {c['fg_dim']}; }}
    QLabel#card_fg {{ font-size: {b - 1}px; color: {c['fg']}; }}
    QToolButton#card_about {{ border: none; background: transparent; color: {c['fg']};
                              font-size: {b - 1}px; padding: 1px 0px; }}
    QToolButton#card_about:hover {{ color: {c['accent_hi']}; }}
    QPushButton#card_gene {{ border: none; background: transparent; color: {c['accent']};
                             font-size: {b - 1}px; text-align: left; padding: 0px 2px; }}
    QPushButton#card_gene:hover {{ color: {c['accent_hi']}; text-decoration: underline; }}
    QFrame#card_rule {{ background: {c['border_soft']}; border: none; }}
    """


# --------------------------------------------------------------------------- the explainer box
class AboutBox(QtWidgets.QWidget):
    """About this test: four fields, each one line until opened."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        title = _role(QtWidgets.QLabel("ABOUT THIS TEST"), "section")
        lay.addWidget(title)
        self.rows = {}
        for field, heading in EX.FIELDS:
            button = _role(QtWidgets.QToolButton(), "about")
            button.setAutoRaise(True)
            button.setCheckable(True)
            button.setArrowType(QtCore.Qt.ArrowType.RightArrow)
            button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
            button.setText(heading)
            button.setToolTip(TH.tip(f"{heading}: click to open or close the full text."))
            line = _role(ElidedLabel(""), "dim")
            full = QtWidgets.QLabel("")
            full.setWordWrap(True)
            full.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
            _role(full, "body")
            full.setContentsMargins(22, 0, 0, 6)
            full.hide()
            head = QtWidgets.QHBoxLayout()
            head.setSpacing(6)
            head.addWidget(button)
            head.addWidget(line, 1)
            lay.addLayout(head)
            lay.addWidget(full)
            button.toggled.connect(lambda on, f=field: self.expand(f, on))
            self.rows[field] = (button, line, full)

    def set_key(self, key: str):
        """Show one strategy's four fields, all collapsed."""
        e = EX.explainer(key)
        for field, (button, line, full) in self.rows.items():
            line.setText(EX.first_sentence(e[field]))
            line.setToolTip(TH.tip(e[field]))
            full.setText(e[field])
            button.setChecked(False)
            self.expand(field, False)

    def expand(self, field: str, on: bool = True):
        """Open (or close) one field to its full text."""
        button, line, full = self.rows[field]
        if button.isChecked() != on:
            button.setChecked(on)
            return
        button.setArrowType(QtCore.Qt.ArrowType.DownArrow if on else
                            QtCore.Qt.ArrowType.RightArrow)
        full.setVisible(on)
        line.setVisible(not on)


# --------------------------------------------------------------------------- worked examples
class ExampleCard(QtWidgets.QFrame):
    """A mini scorecard for one worked example: what was asked, the verdict, four small bars and
    the one or two sentences saying why it failed or what it found."""

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setObjectName("exampleCard")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 8, 10, 10)
        lay.setSpacing(6)
        top = QtWidgets.QHBoxLayout()
        self.title = _role(QtWidgets.QLabel(title), "title")
        self.verdict = Pill("", "fg_dim")
        top.addWidget(self.title)
        top.addStretch(1)
        top.addWidget(self.verdict)
        lay.addLayout(top)
        self.target = _role(ElidedLabel(""), "dim")
        lay.addWidget(self.target)
        self.card = Scorecard(compact=True)
        lay.addWidget(self.card)
        self.text = QtWidgets.QLabel("")
        self.text.setWordWrap(True)
        self.text.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        _role(self.text, "body")
        lay.addWidget(self.text)
        lay.addStretch(1)

    def paintEvent(self, _event):
        c = colors()
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.setPen(QtGui.QPen(QtGui.QColor(c["border_soft"]), 1.0))
        fill = QtGui.QColor(c["surface_alt"])
        fill.setAlphaF(0.55)
        p.setBrush(fill)
        p.drawRoundedRect(QtCore.QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), 10, 10)
        p.end()

    def set_example(self, example: dict, task: str):
        """Fill from one entry of `strategy_examples.json` (a failure or a success)."""
        ex = example or {}
        verdict = str(ex.get("verdict") or "")
        source = ex.get("source")
        if source == "none" or not ex:
            self.verdict.setText("none")
            self.verdict.role = "fg_dim"
            self.target.setText("no run to show")
            self.card.set_bars([])
            self.card.setVisible(False)
            self.text.setText(ex.get("explanation", "No example has been built for this "
                                                    "strategy yet."))
            return
        self.card.setVisible(True)
        self.verdict.setText(verdict + (" - noise table" if source == "noise table" else ""))
        self.verdict.role = VERDICTS.get(verdict, "fg_dim")
        self.verdict.setToolTip(TH.tip(
            "A self-test on the noise table, where every label and edge was dealt out at random: "
            + (ex.get("why_noise") or "no real-data failure to show") + "."
            if source == "noise table" else
            f"A real calibration run on the shipped table (seed {ex.get('seed')}), "
            f"{'at the default settings' if ex.get('at_default') else 'at a tuned setting'}."))
        hidden = ex.get("hidden") or ""
        what = ex.get("target") or ("random labels and edges" if source == "noise table" else "")
        self.target.setText(what or hidden)
        self.target.setToolTip(TH.tip(f"Held out: {hidden}. Compared with: "
                                      f"{ex.get('null_kind') or 'shuffled data'}."))
        self.card.set_bars(SC.headline_bars(task, ex.get("card") or {},
                                            EX.example_verdict(ex)))
        self.text.setText(ex.get("explanation", ""))


class NewInformation(QtWidgets.QWidget):
    """What the successful run added: its top calls for genes without a known label."""

    gene_clicked = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.lay = QtWidgets.QVBoxLayout(self)
        self.lay.setContentsMargins(0, 0, 0, 0)
        self.lay.setSpacing(3)
        self.title = _role(QtWidgets.QLabel("NEW INFORMATION FROM THE RUN THAT WORKED"),
                           "section")
        self.lay.addWidget(self.title)
        self.grid_host = QtWidgets.QWidget()
        self.grid = QtWidgets.QGridLayout(self.grid_host)
        self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.setHorizontalSpacing(10)
        self.grid.setVerticalSpacing(2)
        self.lay.addWidget(self.grid_host)
        self.note = _role(QtWidgets.QLabel(""), "dim")
        self.note.setWordWrap(True)
        self.lay.addWidget(self.note)
        self.buttons: list = []

    def set_new(self, new: dict):
        """Show up to five rows: gene (a button that finds it on the map), product, call, support."""
        while self.grid.count():
            item = self.grid.takeAt(0)
            old = item.widget()
            if old is not None:
                # Hidden and detached now: a deferred delete alone leaves the old row painted
                # under the new one until the event loop next runs its deletions.
                old.hide()
                old.setParent(None)
                old.deleteLater()
        self.buttons = []
        new = new or {}
        rows = new.get("rows") or []
        self.setVisible(bool(rows) or bool(new.get("skipped")))
        # Two lines per row, so a narrow dock loses nothing: the gene and its call, then its
        # product and the support behind the call, small.
        for i, r in enumerate(rows):
            gene = str(r.get("gene_id") or "")
            call_text = str(r.get("call") or "")
            support = str(r.get("support") or "")
            if gene:
                first = gene.split(" + ")[0]
                w = _role(QtWidgets.QPushButton(gene), "gene")
                w.setFlat(True)
                w.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
                w.setToolTip(TH.tip(f"Find {first} on the map and in the search box."))
                w.clicked.connect(lambda _c=False, g=first: self.gene_clicked.emit(g))
                self.buttons.append(w)
                self.grid.addWidget(w, 2 * i, 0)
            call = _role(ElidedLabel(call_text), "fg")
            call.setToolTip(TH.tip(f"{call_text}: {support}"))
            self.grid.addWidget(call, 2 * i, 1 if gene else 0, 1, 1 if gene else 2)
            below = " \u00b7 ".join(x for x in (str(r.get("product") or ""), support) if x)
            sub = _role(ElidedLabel(below), "dim")
            sub.setToolTip(TH.tip(below))
            sub.setContentsMargins(4 if gene else 0, 0, 0, 4)
            self.grid.addWidget(sub, 2 * i + 1, 0, 1, 2)
        self.grid.setColumnStretch(0, 0)
        self.grid.setColumnStretch(1, 1)
        if new.get("skipped"):
            self.note.setText(f"Not re-run: {new['skipped']}.")
        else:
            n = new.get("n_new")
            self.note.setText(f"{n:,} in all; the top {len(rows)} shown." if n else "")


# --------------------------------------------------------------------------- the card
class StrategyCard(QtWidgets.QWidget):
    """The default view of one strategy: header, four bars, buttons, about, examples."""

    run_clicked = QtCore.pyqtSignal()
    test_clicked = QtCore.pyqtSignal()
    details_clicked = QtCore.pyqtSignal()
    gene_clicked = QtCore.pyqtSignal(str)

    def __init__(self, tips: dict | None = None, parent=None):
        super().__init__(parent)
        tips = tips or {}
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(14, 12, 14, 14)
        outer.setSpacing(12)
        # Header: name, chips, question.
        self.name = QtWidgets.QLabel("")
        self.name.setWordWrap(True)
        _role(self.name, "name")
        outer.addWidget(self.name)
        chips = QtWidgets.QHBoxLayout()
        chips.setSpacing(6)
        self.method = Pill("", "fg_muted", filled=False)
        self.grade = Pill("", "fg_dim")
        self.task = Pill("", "fg_muted", filled=False)
        for w in (self.method, self.task, self.grade):
            chips.addWidget(w)
        chips.addStretch(1)
        outer.addLayout(chips)
        self.question = QtWidgets.QLabel("")
        self.question.setWordWrap(True)
        _role(self.question, "question")
        outer.addWidget(self.question)
        # The four bars.
        self.scorecard = Scorecard()
        outer.addWidget(self.scorecard)
        self.basis = _role(QtWidgets.QLabel(""), "dim")
        self.basis.setWordWrap(True)
        outer.addWidget(self.basis)
        # The last self-test run here, in the same four bars (hidden until one runs).
        self.mine = QtWidgets.QWidget()
        ml = QtWidgets.QVBoxLayout(self.mine)
        ml.setContentsMargins(0, 0, 0, 0)
        ml.setSpacing(6)
        mt = QtWidgets.QHBoxLayout()
        mt.addWidget(_role(QtWidgets.QLabel("YOUR LAST TEST"), "section"))
        mt.addStretch(1)
        self.mine_verdict = Pill("", "fg_dim")
        mt.addWidget(self.mine_verdict)
        ml.addLayout(mt)
        self.mine_card = Scorecard(compact=True)
        ml.addWidget(self.mine_card)
        self.mine.hide()
        outer.addWidget(self.mine)
        # Buttons.
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        self.run_btn = QtWidgets.QPushButton("Run")
        self.run_btn.setProperty("primary", True)
        self.test_btn = QtWidgets.QPushButton("Test")
        self.details_btn = QtWidgets.QPushButton("Details ▸")
        for b, key, sig in ((self.run_btn, "run", self.run_clicked),
                            (self.test_btn, "test", self.test_clicked),
                            (self.details_btn, "details", self.details_clicked)):
            b.setToolTip(TH.tip(tips.get(key, b.text())))
            b.clicked.connect(sig.emit)
            row.addWidget(b)
        row.addStretch(1)
        outer.addLayout(row)
        outer.addWidget(self._rule())
        # About this test.
        self.about = AboutBox()
        outer.addWidget(self.about)
        outer.addWidget(self._rule())
        # Worked examples.
        outer.addWidget(_role(QtWidgets.QLabel("WORKED EXAMPLES"), "section"))
        pair = QtWidgets.QHBoxLayout()
        pair.setSpacing(10)
        self.fails = ExampleCard("Fails when…")
        self.works = ExampleCard("Works when…")
        pair.addWidget(self.fails, 1)
        pair.addWidget(self.works, 1)
        outer.addLayout(pair)
        self.new = NewInformation()
        self.new.gene_clicked.connect(self.gene_clicked.emit)
        outer.addWidget(self.new)
        outer.addStretch(1)

    @staticmethod
    def _rule() -> QtWidgets.QFrame:
        line = _role(QtWidgets.QFrame(), "rule")
        line.setFixedHeight(1)
        return line

    def restyle(self) -> bool:
        """Apply the card stylesheet for the theme and text size in force; True if it changed."""
        app = QtWidgets.QApplication.instance()
        theme = str(app.property("starplastTheme") or "dark") if app is not None else "dark"
        own = self.font().pixelSize()             # the application stylesheet's body size
        base = own if own > 0 else 13
        if getattr(self, "_styled", None) == (theme, base):
            return False
        self._styled = (theme, base)
        self.setStyleSheet(card_style(base))
        return True

    def show_strategy(self, s, organism: str, card: dict, verdict: dict, grade: str,
                      basis: str):
        """Fill the card for strategy `s`: `card` and `verdict` as `headline_bars` takes them."""
        self.key, self.task_key = s.key, s.task
        self.restyle()
        self.name.setText(f"{s.number:02d} · {s.title}")
        self.method.setText(s.method)
        self.method.setToolTip(TH.tip(f"Method: {s.method}. Details, Guide lists each technique "
                                      f"it is built from and why."))
        self.task.setText(s.task)
        self.task.setToolTip(TH.tip(f"Task: {SC.TASKS[s.task].description} Every strategy doing "
                                    f"this task shows the same four bars."))
        role, tip = GRADES.get(grade, ("fg_dim", "Not calibrated on this organism: press Test "
                                                 "to measure it on your data."))
        self.grade.setText(grade or "not calibrated")
        self.grade.role = role
        self.grade.setToolTip(TH.tip(tip))
        self.question.setText(s.question)
        self.scorecard.set_bars(SC.headline_bars(s.task, card, verdict))
        self.basis.setText(basis)
        self.mine.hide()
        self.about.set_key(s.key)
        ex = EX.examples(s.key, organism)
        self.fails.set_example(ex.get("failure") or {}, s.task)
        self.works.set_example(ex.get("success") or {}, s.task)
        self.new.set_new((ex.get("success") or {}).get("new") or {})
        for w in (self.method, self.grade, self.task):
            w.updateGeometry()
            w.update()

    def show_test(self, test):
        """Put a finished self-test on the card in the same four bars."""
        self.mine_verdict.setText(test.verdict)
        self.mine_verdict.role = VERDICTS.get(test.verdict, "fg_dim")
        self.mine_verdict.setToolTip(TH.tip(test.summary()))
        task = test.task if test.task in SC.TASKS else getattr(self, "task_key", SC.T_LABEL)
        self.mine_card.set_bars(SC.headline_bars(
            task, dict(test.scorecard or {}),
            {"skill": test.skill, "observed": test.observed, "chance": test.null_mean}))
        self.mine.show()
        self.mine_verdict.update()
