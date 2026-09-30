"""The Start-here tab: one question at a time, ending at the strategies worth running.

Beside the Strategies tab, which is left exactly as it was. The Strategies tab answers "what is
there"; this one answers "where do I begin", for someone who arrives with a gene list from a screen,
their own measurement, or nothing but curiosity.

    question    one per screen, in a few plain words, with large choices and nothing else on the
                screen to read; a gene step searches, a gene-list step pastes, loads, takes the
                genes gated on the map or offers a real example, and a column step shows each
                column's coverage so an almost-empty one is never picked by accident
    trail       every answer so far along the top, each one a button back to the question that
                asked it, so nothing is a one-way door
    the end     three to six recommended strategies, each as the card summary the Strategies tab
                already shows -- name, method, calibration grade and the same four headline bars
                (`strategy_panel.headline`, `strategy_card`, `scorecard.HEADLINE`) -- with one line
                on why it is recommended, and Run it here / Open in Strategies

The questions and the ranking are :mod:`starplast.guided`, which holds no Qt: this module draws
them and nothing more. It also holds no window: it emits what the user asked for (`run_strategy`,
`open_strategy`, `open_view`) and :func:`install` wires those to the application's strategy panel
and docks, which is what lets the whole panel be built and walked in a test with no window at all.
"""
from __future__ import annotations

from PyQt6 import QtCore, QtGui, QtWidgets

from . import guided as G
from . import scorecard as SC
from . import strategies as S
from . import strategy_card as CARD
from . import theme as TH

#: Every button's tooltip, in one place, because the panel builds its controls in several methods
#: and a control without a tooltip is a control nobody can ask about.
TIPS = {
    "back": "Back to the previous question. Your later answers are dropped, because they were "
            "answers to questions that followed from this one.",
    "restart": "Start again from the first question. Nothing is kept.",
    "next": "Take what is entered here and go on to the next question.",
    "trail": "Go back to this question and answer it differently. Everything after it is asked "
             "again.",
    "search": "Type an accession, a symbol, or a word from the product description; the matches "
              "appear below.",
    "matches": "The genes matching what you typed, with their product. Click one to choose it.",
    "paste": "Paste or type identifiers, one per line or separated by commas or spaces. Anything "
             "this table does not hold is reported rather than silently dropped.",
    "file": "Read identifiers from a text or CSV file: the first column of each line, or the whole "
            "file if it is a plain list.",
    "gate": "Use the genes currently gated on the 3D map, so a group you drew around can be handed "
            "straight to the questions.",
    "example": "Fill the list with the members of one known category, so the branch can be tried "
               "before you have a list of your own.",
    "filter": "Type to shorten the list to the columns whose name contains the text.",
    "run": "Run this strategy now, with the settings your answers decided. It runs in the "
           "background like any job, and its tables appear under Strategies, Results.",
    "open": "Open this strategy in the Strategies tab with the settings your answers decided, "
            "where its guide, every setting and its self-test are.",
    "view": "Open this panel: it answers the question you asked better than running anything does.",
    "option": "Choose this and go on to the next question.",
}
#: A column list longer than this gets a filter box above it. Eighteen labels are a list; three
#: hundred and fifty measurements are not.
FILTER_ABOVE = 12
#: How many genes a search offers at once.
MATCH_LIMIT = 40


def _role(widget: QtWidgets.QWidget, role: str) -> QtWidgets.QWidget:
    """Name a widget for `strategy_card.card_style`, which sizes and colours these roles already."""
    widget.setObjectName(f"card_{role}")
    return widget


class Choice(QtWidgets.QAbstractButton):
    """One answer, as a wide flat row: its label, its coverage or kind on the right, and a border.

    Painted rather than styled so the row reads the same on all four themes and so the detail can
    be dimmer than the label, which a single stylesheet on a QPushButton cannot do.
    """

    def __init__(self, option: G.Option, parent=None):
        super().__init__(parent)
        self.option = option
        self.setCheckable(False)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setToolTip(TH.tip(option.hint or TIPS["option"]))
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Fixed)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_Hover, True)

    def sizeHint(self) -> QtCore.QSize:
        """Two lines of breathing room around one line of text."""
        return QtCore.QSize(260, QtGui.QFontMetrics(self.font()).height() + 22)

    def paintEvent(self, _event):
        c = CARD.colors()
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        r = QtCore.QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        hot = self.underMouse() or self.hasFocus()
        fill = QtGui.QColor(c["accent"] if hot else c["surface"])
        fill.setAlphaF(0.16 if hot else 0.5)
        p.setBrush(fill)
        p.setPen(QtGui.QPen(QtGui.QColor(c["accent"] if hot else c["border"]), 1.0))
        p.drawRoundedRect(r, 7, 7)
        inner = r.adjusted(12, 0, -12, 0)
        p.setPen(QtGui.QColor(c["fg"]))
        p.drawText(inner, int(QtCore.Qt.AlignmentFlag.AlignLeft
                              | QtCore.Qt.AlignmentFlag.AlignVCenter), self.option.label)
        if self.option.detail:
            p.setPen(QtGui.QColor(c["fg_dim"]))
            p.drawText(inner, int(QtCore.Qt.AlignmentFlag.AlignRight
                                  | QtCore.Qt.AlignmentFlag.AlignVCenter), self.option.detail)
        p.end()


class Card(QtWidgets.QFrame):
    """A frame whose height follows from its width, which is what wrapped text needs in a scroller.

    A QScrollArea sizes its contents to their minimum height, and a word-wrapped QLabel has no
    minimum height until its width is known: without this the recommendation cards collapsed on top
    of one another and their buttons vanished. Declaring height-for-width and handing the question
    to the card's own layout is Qt's answer to that, and `GuidedPanel._fit_body` asks the same
    question of the column of cards.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        policy = self.sizePolicy()
        policy.setHeightForWidth(True)
        policy.setVerticalPolicy(QtWidgets.QSizePolicy.Policy.Minimum)
        self.setSizePolicy(policy)

    def hasHeightForWidth(self) -> bool:
        """Yes: the height depends on how wide the card is drawn."""
        return True

    def heightForWidth(self, width: int) -> int:
        """How tall this card must be at `width`, as its own layout works it out."""
        layout = self.layout()
        return layout.heightForWidth(width) if layout is not None else super().heightForWidth(width)

    def minimumSizeHint(self) -> QtCore.QSize:
        """Never shorter than the layout needs at the width it currently has."""
        hint = super().minimumSizeHint()
        return QtCore.QSize(hint.width(), max(hint.height(), self.heightForWidth(max(self.width(),
                                                                                    hint.width()))))


class RecommendationCard(Card):
    """One recommended strategy, as the Strategies tab's card summarises it, plus why and Run.

    The four bars are `scorecard.headline_bars` drawn by `strategy_card.Scorecard`, filled from
    `strategy_panel.headline` -- the same numbers the Strategies tab shows for the same strategy,
    computed in one place rather than two.
    """

    run = QtCore.pyqtSignal(str, object)
    open_it = QtCore.pyqtSignal(str, object)

    def __init__(self, rec: G.Recommendation, organism: str, parent=None):
        super().__init__(parent)
        from . import strategy_panel as SP
        self.rec = rec
        self.setObjectName("guided_card")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 12)
        lay.setSpacing(7)
        title = _role(QtWidgets.QLabel(f"{rec.number:02d} · {rec.title}"), "title")
        title.setWordWrap(True)
        lay.addWidget(title)
        chips = QtWidgets.QHBoxLayout()
        chips.setSpacing(6)
        role, tip = CARD.GRADES.get(rec.grade, ("fg_dim", "Not calibrated on this organism: open "
                                                          "it in Strategies and press Test."))
        chips.addWidget(CARD.Pill(rec.task, "fg_muted", filled=False,
                                  tip=SC.TASKS[rec.task].description))
        chips.addWidget(CARD.Pill(rec.grade or "not calibrated", role, tip=tip))
        chips.addStretch(1)
        lay.addLayout(chips)
        why = _role(QtWidgets.QLabel(rec.why), "body")
        why.setWordWrap(True)
        why.setToolTip(TH.tip("Why this strategy is on your list: what it asks, what of yours it "
                              "uses, and what the calibration sweep measured for it here."))
        lay.addWidget(why)
        card, verdict, _grade, basis = SP.headline(rec.key, organism)
        bars = CARD.Scorecard(compact=True)
        bars.set_bars(SC.headline_bars(rec.task, card, verdict))
        bars.setToolTip(TH.tip(basis))
        lay.addWidget(bars)
        if rec.settings:
            filled = ", ".join(f"{k}={_short(v)}" for k, v in sorted(rec.settings.items()))
            note = _role(QtWidgets.QLabel(f"Settings from your answers: {filled}"), "dim")
            note.setWordWrap(True)
            note.setToolTip(TH.tip("These settings are filled in for you from the answers you "
                                   "gave; everything else stays at the strategy's own default."))
            lay.addWidget(note)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        self.run_btn = QtWidgets.QPushButton("Run it here")
        self.run_btn.setProperty("primary", True)
        self.run_btn.setToolTip(TH.tip(TIPS["run"]))
        self.run_btn.clicked.connect(lambda: self.run.emit(rec.key, dict(rec.settings)))
        self.open_btn = QtWidgets.QPushButton("Open in Strategies")
        self.open_btn.setToolTip(TH.tip(TIPS["open"]))
        self.open_btn.clicked.connect(lambda: self.open_it.emit(rec.key, dict(rec.settings)))
        row.addWidget(self.run_btn)
        row.addWidget(self.open_btn)
        row.addStretch(1)
        lay.addLayout(row)

    def paintEvent(self, _event):
        c = CARD.colors()
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        r = QtCore.QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        fill = QtGui.QColor(c["surface"])
        fill.setAlphaF(0.55)
        p.setBrush(fill)
        p.setPen(QtGui.QPen(QtGui.QColor(c["border_soft"]), 1.0))
        p.drawRoundedRect(r, 9, 9)
        p.end()


def _short(value) -> str:
    """One setting's value in a few characters: a long gene list becomes its size."""
    if isinstance(value, (list, tuple)):
        return f"{len(value):,} genes"
    text = str(value)
    return text if len(text) <= 28 else text[:25] + "..."


class GuidedPanel(QtWidgets.QWidget):
    """The guided tab: the current question, the trail of answers, and the recommendations."""

    status = QtCore.pyqtSignal(str)
    #: A strategy the user wants run now, and the settings the answers decided.
    run_strategy = QtCore.pyqtSignal(str, object)
    #: A strategy the user wants to look at in the Strategies tab, with the same settings.
    open_strategy = QtCore.pyqtSignal(str, object)
    #: Another panel that answers the question better: `guided.VIEW_MAPS` or `guided.VIEW_STAR`.
    open_view = QtCore.pyqtSignal(str)

    def __init__(self, ctx, gated=None, parent=None):
        """Build the panel over one organism's context; other spaces are loaded if they are asked
        for. `gated` is a callable returning the gene ids gated on the map, or None."""
        super().__init__(parent)
        self.base, self.gated = ctx, gated
        self._contexts = {ctx.organism: ctx}
        self.answers: dict = {}
        #: The recommendation cards on screen, empty until the last question is answered.
        self.cards: list = []
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(14, 12, 14, 12)
        outer.setSpacing(10)
        self.trail_host = QtWidgets.QWidget()
        self.trail = QtWidgets.QHBoxLayout(self.trail_host)
        self.trail.setContentsMargins(0, 0, 0, 0)
        self.trail.setSpacing(6)
        outer.addWidget(self.trail_host)
        self.question = _role(QtWidgets.QLabel(""), "name")
        self.question.setWordWrap(True)
        outer.addWidget(self.question)
        self.note = _role(QtWidgets.QLabel(""), "dim")
        self.note.setWordWrap(True)
        outer.addWidget(self.note)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.body = QtWidgets.QWidget()
        self.body_layout = QtWidgets.QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 0, 6, 0)
        self.body_layout.setSpacing(8)
        self.scroll.setWidget(self.body)
        outer.addWidget(self.scroll, 1)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        self.back_btn = QtWidgets.QPushButton("‹ Back")
        self.back_btn.setToolTip(TH.tip(TIPS["back"]))
        self.back_btn.clicked.connect(self.back)
        self.next_btn = QtWidgets.QPushButton("Next ›")
        self.next_btn.setProperty("primary", True)
        self.next_btn.setToolTip(TH.tip(TIPS["next"]))
        self.next_btn.clicked.connect(self._accept)
        self.restart_btn = QtWidgets.QPushButton("Start over")
        self.restart_btn.setToolTip(TH.tip(TIPS["restart"]))
        self.restart_btn.clicked.connect(self.restart)
        row.addWidget(self.back_btn)
        row.addWidget(self.next_btn)
        row.addStretch(1)
        row.addWidget(self.restart_btn)
        outer.addLayout(row)
        self.render()

    # ------------------------------------------------------------------ state
    def context(self):
        """The table the questions are asked of: the space chosen, or the one this panel was built
        on until one is."""
        code = self.answers.get("space") or self.base.organism
        if code not in self._contexts:
            try:
                self._contexts[code] = S.shipped(code)
            except Exception as e:      # a space whose table is not installed here
                self.status.emit(f"{code} is not installed here ({e}); staying on "
                                 f"{self.base.organism}")
                self.answers["space"] = self.base.organism
                return self.base
        return self._contexts[code]

    def answer(self, step: str, value):
        """Answer one question and move on. Later answers are dropped, since they followed it."""
        self.drop_after(step)
        self.answers[step] = value
        self.render()

    def drop_after(self, step: str):
        """Forget every answer that came after `step`, which a changed answer has invalidated."""
        chain = G.route(self.answers)
        if step in chain:
            for later in chain[chain.index(step):]:
                self.answers.pop(later, None)

    def back(self):
        """Ask the last answered question again, or do nothing on the first one."""
        answered = [k for k in G.route(self.answers) if G.answered(self.answers, k)]
        if answered:
            self.goto(answered[-1])

    def goto(self, step: str):
        """Ask `step` again, dropping it and everything after it."""
        self.drop_after(step)
        self.render()

    def restart(self):
        """Back to the first question with nothing kept."""
        self.answers = {}
        self.render()

    # ------------------------------------------------------------------ drawing
    def render(self):
        """Draw the question due now, or the recommendations when there is none left."""
        self.restyle()
        _clear(self.trail)
        for step, _value, text in G.trail(self.answers):
            button = QtWidgets.QToolButton()
            button.setText(f"{G.STEPS[step].question.rstrip('?')} · {text}")
            button.setObjectName("guided_crumb")
            button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            button.setToolTip(TH.tip(TIPS["trail"]))
            button.clicked.connect(lambda _c=False, s=step: self.goto(s))
            self.trail.addWidget(button)
        self.trail.addStretch(1)
        self.trail_host.setVisible(bool(self.trail.count() > 1))
        _clear(self.body_layout)
        self.cards = []                 # the cards of the last end screen are gone with the layout
        step = G.next_step(self.answers)
        self.back_btn.setEnabled(bool(self.answers))
        if not step:
            self._draw_end()
            return
        spec = G.STEPS[step]
        self.question.setText(spec.question)
        self.question.setToolTip(TH.tip(spec.tip))
        self.note.setText(spec.note)
        self.note.setVisible(bool(spec.note))
        self.current = step
        self.next_btn.setVisible(spec.kind in ("gene", "genes"))
        if spec.kind == "choice":
            self._draw_choices(step, G.options(step, self.answers, self.context()))
        elif spec.kind == "column":
            self._draw_columns(step)
        elif spec.kind == "gene":
            self._draw_gene()
        else:
            self._draw_genes()
        self.body_layout.addStretch(1)
        self._fit_body()

    def resizeEvent(self, event):
        """Wrapped text is taller in a narrow panel, so the scrolled column is measured again."""
        super().resizeEvent(event)
        self._fit_body()

    def _fit_body(self):
        """Give the scrolled column the height its wrapped text needs at the width it has.

        Deferred as well as done now: a widget added to a layout is not shown until the event loop
        next runs, and a layout does not count a hidden widget -- measured in the same breath as the
        cards are built, the column comes out empty and every card is squashed to a sliver.
        """
        QtCore.QTimer.singleShot(0, self._measure_body)
        self._measure_body()

    def _measure_body(self):
        """The measurement itself: the column's height at the viewport's width, never less than the
        minimum its cards ask for."""
        self.body_layout.activate()
        width = self.scroll.viewport().width()
        height = self.body_layout.heightForWidth(width) if width > 0 else 0
        self.body.setMinimumHeight(max(0, height, self.body_layout.minimumSize().height()))

    def restyle(self) -> bool:
        """Apply the card stylesheet of the theme and text size in force; True if it changed."""
        app = QtWidgets.QApplication.instance()
        theme = str(app.property("starplastTheme") or "dark") if app is not None else "dark"
        own = self.font().pixelSize()
        base = own if own > 0 else 13
        if getattr(self, "_styled", None) == (theme, base):
            return False
        self._styled = (theme, base)
        c = TH.palette_for(theme)
        self.setStyleSheet(CARD.card_style(base) + f"""
        QToolButton#guided_crumb {{ border: 1px solid {c['border']}; border-radius: 9px;
                                    padding: 2px 9px; color: {c['fg_muted']};
                                    background: transparent; font-size: {max(8, base - 2)}px; }}
        QToolButton#guided_crumb:hover {{ color: {c['accent_hi']}; border-color: {c['accent']}; }}
        """)
        return True

    def _draw_choices(self, step: str, options):
        for option in options:
            button = Choice(option)
            button.clicked.connect(lambda _c=False, s=step, v=option.value: self.answer(s, v))
            self.body_layout.addWidget(button)
        if not options:
            self.body_layout.addWidget(_role(QtWidgets.QLabel(
                "Nothing to choose here on this table. Go back and answer the question before "
                "this one differently."), "body"))

    def _draw_columns(self, step: str):
        options = G.options(step, self.answers, self.context())
        rows = QtWidgets.QWidget()
        rows_layout = QtWidgets.QVBoxLayout(rows)
        rows_layout.setContentsMargins(0, 0, 0, 0)
        rows_layout.setSpacing(6)
        buttons = []
        for option in options:
            button = Choice(option)
            button.clicked.connect(lambda _c=False, s=step, v=option.value: self.answer(s, v))
            rows_layout.addWidget(button)
            buttons.append(button)
        if len(options) > FILTER_ABOVE:
            box = QtWidgets.QLineEdit()
            box.setPlaceholderText("filter…")
            box.setToolTip(TH.tip(TIPS["filter"]))
            box.textChanged.connect(lambda text: _filter(buttons, text))
            self.body_layout.addWidget(box)
        self.body_layout.addWidget(rows)

    def _draw_gene(self):
        ctx = self.context()
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText("accession, symbol or product…")
        self.search.setToolTip(TH.tip(TIPS["search"]))
        self.matches = QtWidgets.QListWidget()
        self.matches.setToolTip(TH.tip(TIPS["matches"]))
        self.matches.itemActivated.connect(self._pick_gene)
        self.matches.itemClicked.connect(self._pick_gene)
        self.search.textChanged.connect(lambda text: self._fill_matches(ctx, text))
        self.body_layout.addWidget(self.search)
        self.body_layout.addWidget(self.matches, 1)
        self._fill_matches(ctx, "")

    def _fill_matches(self, ctx, text: str):
        self.matches.clear()
        query = str(text).strip().lower()
        if not query:
            return
        ids = ctx.gene_ids
        products = ctx.product(range(ctx.n))
        for i, gene in enumerate(ids):
            if query in str(gene).lower() or query in str(products[i]).lower():
                item = QtWidgets.QListWidgetItem(f"{gene}  ·  {products[i][:70]}")
                item.setData(QtCore.Qt.ItemDataRole.UserRole, str(gene))
                self.matches.addItem(item)
            if self.matches.count() >= MATCH_LIMIT:
                break

    def _pick_gene(self, item):
        self.answer("gene", str(item.data(QtCore.Qt.ItemDataRole.UserRole)))

    def _draw_genes(self):
        self.genes_box = QtWidgets.QPlainTextEdit()
        self.genes_box.setPlaceholderText("one accession per line, or separated by commas or spaces")
        self.genes_box.setToolTip(TH.tip(TIPS["paste"]))
        self.genes_box.textChanged.connect(self._count_genes)
        self.body_layout.addWidget(self.genes_box, 1)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        load = QtWidgets.QPushButton("Load a file…")
        load.setToolTip(TH.tip(TIPS["file"]))
        load.clicked.connect(lambda: self.load_file())
        gate = QtWidgets.QPushButton("Use the gated genes")
        gate.setToolTip(TH.tip(TIPS["gate"]))
        gate.setEnabled(self.gated is not None)
        gate.clicked.connect(self.use_gated)
        example = QtWidgets.QPushButton("Try an example set")
        example.setToolTip(TH.tip(TIPS["example"]))
        example.clicked.connect(self.use_example)
        for b in (load, gate, example):
            row.addWidget(b)
        row.addStretch(1)
        self.body_layout.addLayout(row)
        self.count = _role(QtWidgets.QLabel(""), "dim")
        self.count.setToolTip(TH.tip("How many of the identifiers you gave this table holds. "
                                     "Anything it does not hold is named, never dropped in "
                                     "silence."))
        self.body_layout.addWidget(self.count)
        self._count_genes()

    def _count_genes(self):
        ctx = self.context()
        found, missing = ctx.resolve_genes(self.genes_box.toPlainText())
        self.count.setText(
            f"{len(found):,} genes found in this table"
            + (f"; {len(missing):,} not found: {', '.join(missing[:6])}"
               f"{'…' if len(missing) > 6 else ''}" if missing else ""))
        self.next_btn.setEnabled(bool(len(found)))

    def set_genes(self, genes) -> int:
        """Put a gene list in the box (and say how many landed); for the buttons and for tests."""
        self.genes_box.setPlainText("\n".join(str(g) for g in genes))
        return len(list(genes))

    def load_file(self, path: str = "") -> int:
        """Fill the gene box from a text or CSV file. Returns how many identifiers were read."""
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Gene list", "", "Text or CSV (*.txt *.csv *.tsv);;All (*)")
        if not path:
            return 0
        tokens = G.read_gene_file(path)
        self.set_genes(tokens)
        self.status.emit(f"read {len(tokens):,} identifiers from {path}")
        return len(tokens)

    def use_gated(self) -> int:
        """Fill the gene box with the genes gated on the 3D map."""
        genes = list(self.gated() or []) if self.gated is not None else []
        self.set_genes(genes)
        self.status.emit(f"{len(genes):,} gated genes taken from the map" if genes else
                         "nothing is gated on the map: drag a selection there first")
        return len(genes)

    def use_example(self) -> int:
        """Fill the gene box with the members of one known category of the default label."""
        category, genes = G.example_genes(self.context(), self.answers)
        self.set_genes(genes)
        self.status.emit(f"example set: the {len(genes):,} genes of {category!r}" if category else
                         "no category of a usable size to take an example from")
        return len(genes)

    def _accept(self):
        """Take what the current step's control holds and go on -- the Next button."""
        step = getattr(self, "current", "")
        if step == "genes":
            found, _missing = self.context().resolve_genes(self.genes_box.toPlainText())
            genes = [str(g) for g in self.context().gene_ids[found]]
            if genes:
                self.answer("genes", genes)
        elif step == "gene":
            found, _missing = self.context().resolve_genes(self.search.text())
            if len(found):
                self.answer("gene", str(self.context().gene_ids[found[0]]))
            else:
                self.status.emit("no gene of that name in this table; pick one from the list")

    # ------------------------------------------------------------------ the end
    def recommendations(self) -> list:
        """The ranked strategies for the answers as they stand."""
        return G.recommend(self.answers, self.context())

    def _draw_end(self):
        ctx = self.context()
        recs = self.recommendations()
        self.question.setText("Start with one of these")
        self.question.setToolTip(TH.tip(
            "Ranked for what you said you have and what you want to know: what each strategy is "
            "for, whether its task matches your goal, and what the calibration sweep measured for "
            "it on this organism."))
        self.note.setText(f"{len(recs)} of {len(S.catalog())} strategies, best first. The four bars "
                          f"are the same on every card, so two can be compared by eye.")
        self.note.setVisible(True)
        self.next_btn.setVisible(False)
        caveat = G.caveat(recs)
        if caveat:
            warn = _role(QtWidgets.QLabel(caveat), "body")
            warn.setWordWrap(True)
            warn.setToolTip(TH.tip("Said rather than hidden: nothing calibrated as reliable for "
                                   "this exists on this table."))
            c = CARD.colors()
            warn.setStyleSheet(f"color: {c['warning']};")
            self.body_layout.addWidget(warn)
        self.cards = []
        for rec in recs:
            card = RecommendationCard(rec, ctx.organism)
            card.run.connect(self.run_strategy.emit)
            card.open_it.connect(self.open_strategy.emit)
            self.body_layout.addWidget(card)
            self.cards.append(card)
        for view in G.views(self.answers, ctx):
            self.body_layout.addWidget(self._view_row(view))
        self.body_layout.addStretch(1)
        self._fit_body()

    def _view_row(self, view: dict) -> QtWidgets.QWidget:
        host = Card()
        lay = QtWidgets.QHBoxLayout(host)
        lay.setContentsMargins(0, 4, 0, 4)
        lay.setSpacing(8)
        text = _role(QtWidgets.QLabel(f"<b>{view['title']}</b> — {view['why']}"), "body")
        text.setWordWrap(True)
        text.setToolTip(TH.tip(view["why"]))
        lay.addWidget(text, 1)
        button = QtWidgets.QPushButton(f"Open the {view['title'].lower()}")
        button.setToolTip(TH.tip(TIPS["view"]))
        button.clicked.connect(lambda _c=False, v=view["view"]: self.open_view.emit(v))
        lay.addWidget(button)
        return host


def _filter(buttons, text: str):
    """Show only the column rows whose name or column name contains `text`."""
    wanted = str(text).strip().lower()
    for button in buttons:
        button.setVisible(wanted in button.option.label.lower()
                          or wanted in str(button.option.value).lower())


def _clear(layout: QtWidgets.QLayout):
    """Empty a layout, deleting what was in it -- every question rebuilds its own controls."""
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)
            widget.deleteLater()
        elif item.layout() is not None:
            _clear(item.layout())


def install(window) -> GuidedPanel | None:
    """Add the Start-here tab to the main window, beside Strategies, and wire what it asks for.

    Everything it emits is carried out by the panels that already do those things: a strategy is
    selected, filled in and run in the Strategies tab, so it appears in the Jobs panel and its
    tables land where every other strategy's tables land; a view request raises that dock.
    """
    from . import organisms
    strategies = getattr(window, "strategy_panel", None)
    ctx = getattr(strategies, "ctx", None)
    if ctx is None:
        code = organisms.detect(window.nodes["gene_id"])
        if code is None:
            return None
        ctx = S.Context(window.nodes, organism=code)
    panel = GuidedPanel(ctx, gated=getattr(strategies, "gated", None))
    panel.status.connect(lambda m: window.statusBar().showMessage(m))

    def apply(key, settings, run=False):
        if strategies is None:
            window.statusBar().showMessage("the Strategies tab is not available in this window")
            return
        strategies.select(key)
        for name, value in (settings or {}).items():
            try:
                strategies.set_setting(name, value)
            except (KeyError, ValueError, TypeError):
                window.statusBar().showMessage(f"{name} could not be filled in; kept its default")
        dock = getattr(window, "strategies_dock", None)
        if run:
            strategies.run_current()
            window.statusBar().showMessage(f"running {key}; watch Jobs, then Strategies ▸ Results")
        elif dock is not None:
            strategies.show_details(1)
            dock.raise_()

    panel.run_strategy.connect(lambda key, settings: apply(key, settings, run=True))
    panel.open_strategy.connect(lambda key, settings: apply(key, settings, run=False))

    def show_view(name):
        dock = getattr(window, "star_map_dock" if name == G.VIEW_STAR else "maps_dock", None)
        if dock is None:
            window.statusBar().showMessage(f"the {name} is not available in this window")
            return
        dock.show()
        dock.raise_()

    panel.open_view.connect(show_view)
    dock = QtWidgets.QDockWidget("start here", window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    dock.setWidget(panel)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = getattr(window, "strategies_dock", None) or getattr(window, "analysis_dock", None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
    if getattr(window, "right_dock", None) is not None:
        window.right_dock.raise_()
    window.guided, window.guided_dock = panel, dock
    return panel
