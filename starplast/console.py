"""A console pane, and the stdout capture that feeds it.

The build and the search both report progress by printing, and in a GUI that output goes nowhere --
the user watches a frozen window while the terminal they cannot see fills up. This tees stdout and
stderr into a widget without taking them away from the terminal, so running from a shell still
behaves normally and running from a launcher is no longer silent.

Capture is deliberately a tee rather than a redirect. Swallowing stderr would hide tracebacks from
anyone debugging from a terminal, which is where they are actually read.

Lines that came through `logging_util` carry their level as a `[LEVEL]` prefix, and the pane reads
it: it colors by level and filters by minimum level, so a warning is findable in a walk that
printed four thousand progress lines. Ordinary prints are left levelless rather than being called
INFO -- most output here is `print`, and folding it into a level would make filtering by that level
useless.

Lines are appended as they arrive rather than by re-rendering the whole pane. Re-rendering was
quadratic in the number of lines, on the one screen whose entire purpose is watching something long
run: 5,000 lines meant 5,000 full rebuilds of a 5,000-line document.
"""
from __future__ import annotations

import sys
from datetime import datetime

from PyQt6 import QtCore, QtGui, QtWidgets

from .logging_util import LEVELS, parse_level

MAX_LINES = 5000        # a runaway loop must not exhaust memory through the log widget

#: One colour per level, chosen to read on both the dark and the light themes rather than to match
#: either. Levelless output -- ordinary prints -- takes the pane's own text colour, so the default
#: case looks exactly as it did.
LEVEL_COLOUR = {"DEBUG": "#7f8c9b", "INFO": "#5aa9e6", "WARNING": "#d79a2b", "ERROR": "#e05561"}


class Tee(QtCore.QObject):
    """Forwards writes to the original stream and emits them as a signal.

    A file-like object handed to `sys.stdout` is written to from worker threads, so the text is
    emitted as a queued signal rather than touching a widget directly. Writing to a QTextEdit off
    the GUI thread is a crash that shows up as an unrelated stack days later.
    """

    text = QtCore.pyqtSignal(str, bool)          # chunk, is_error

    def __init__(self, stream, is_error: bool = False, parent=None):
        super().__init__(parent)
        self.stream = stream
        self.is_error = is_error

    def write(self, chunk):
        """Write to the original stream and emit the text for the pane."""
        if self.stream is not None:
            try:
                self.stream.write(chunk)
            except Exception:
                pass          # a closed original stream must not break the application's printing
        if chunk:
            self.text.emit(str(chunk), self.is_error)
        return len(chunk) if chunk else 0

    def flush(self):
        """Flush the original stream, tolerating one that has been closed."""
        if self.stream is not None:
            try:
                self.stream.flush()
            except Exception:
                pass

    def close(self):
        """Do NOT close the wrapped stream -- this is a tee, not an owner.

        Python's logging shutdown calls close() on whatever it finds at sys.stderr. Without this the
        interpreter reported "'Tee' object has no attribute 'close'" at exit, and closing the real
        stream instead would take stderr down for everything still running.
        """
        self.flush()

    def isatty(self):
        # Progress bars and colour codes key off this. The widget is not a terminal, and claiming
        # otherwise fills the pane with escape sequences.
        """False: the pane is not a terminal, and claiming otherwise fills it with escapes."""
        return False


class ConsolePanel(QtWidgets.QWidget):
    """Read-only log with a timestamp column, a filter, and copy/clear."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._install = None
        L = QtWidgets.QVBoxLayout(self)
        L.setContentsMargins(6, 6, 6, 6)
        L.setSpacing(6)

        bar = QtWidgets.QHBoxLayout()
        self.filter = QtWidgets.QLineEdit(placeholderText="filter…")
        self.filter.setToolTip("Show only lines containing this text. The log is kept whole "
                               "underneath, so clearing the filter brings everything back.")
        self.filter.textChanged.connect(self._apply_filter)
        self.errors_only = QtWidgets.QCheckBox("errors")
        self.errors_only.setToolTip("Show only what was written to stderr.")
        self.errors_only.toggled.connect(self._apply_filter)
        self.level = QtWidgets.QComboBox()
        self.level.addItem("all levels", "")
        for name in LEVELS:
            self.level.addItem(f"{name} and above", name)
        self.level.setToolTip(
            "Hide anything logged below this level. A walk prints thousands of progress lines and "
            "one warning; without this the warning is somewhere in the middle of them. Ordinary "
            "printed output carries no level and is always shown, because most of what this "
            "application writes is a print and folding it into INFO would make INFO useless.")
        self.level.currentIndexChanged.connect(self._apply_filter)
        copy = QtWidgets.QPushButton("copy")
        copy.setToolTip("Copy the visible lines to the clipboard.")
        copy.clicked.connect(self.copy_all)
        clear = QtWidgets.QPushButton("clear")
        clear.setToolTip("Discard everything logged so far. The log is capped anyway, so this is "
                         "for making the next run's output easy to find rather than for memory.")
        clear.clicked.connect(self.clear)
        for w in (self.filter, self.level, self.errors_only, copy, clear):
            bar.addWidget(w)
        bar.setStretch(0, 1)
        L.addLayout(bar)

        # QTextEdit rather than QPlainTextEdit: a level has to be able to colour its own line, and
        # the plain widget has one colour for the whole document.
        self.view = QtWidgets.QTextEdit()
        self.view.setReadOnly(True)
        self.view.document().setMaximumBlockCount(MAX_LINES)
        self.view.setFont(QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.SystemFont.FixedFont))
        self.view.setLineWrapMode(QtWidgets.QTextEdit.LineWrapMode.NoWrap)
        L.addWidget(self.view, 1)

        self._lines: list[tuple[str, str, bool, str]] = []   # time, text, is_error, level ("" = none)
        self._partial = {False: "", True: ""}

    # -------------------------------------------------------------- capture
    def install(self):
        """Tee stdout and stderr into this pane. Returns the pair for later removal."""
        if self._install is not None:
            return self._install
        out, err = Tee(sys.stdout, False), Tee(sys.stderr, True)
        # Queued, because these are written from worker threads.
        out.text.connect(self.append, QtCore.Qt.ConnectionType.QueuedConnection)
        err.text.connect(self.append, QtCore.Qt.ConnectionType.QueuedConnection)
        self._install = (sys.stdout, sys.stderr)
        sys.stdout, sys.stderr = out, err
        return self._install

    def remove(self):
        """Put the original streams back. Called on close, so a crash does not leave a dead tee."""
        if self._install is None:
            return
        sys.stdout, sys.stderr = self._install
        self._install = None

    # -------------------------------------------------------------- content
    def append(self, chunk: str, is_error: bool = False):
        """Add text, splitting on newlines and holding an unterminated tail.

        print() arrives as two writes -- the text, then the newline -- so appending each write as its
        own line would double every line in the pane.
        """
        buf = self._partial[is_error] + chunk
        *lines, self._partial[is_error] = buf.split("\n")
        stamp = datetime.now().strftime("%H:%M:%S")
        at_end = self._at_end()
        for line in lines:
            level = parse_level(line) or ""
            self._lines.append((stamp, line, is_error, level))
            # Appended rather than re-rendered. Rebuilding the whole document per line is quadratic,
            # and this pane exists to be watched during the runs that print the most.
            if self._matches(stamp, line, is_error, level):
                self._write(stamp, line, level)
        if len(self._lines) > MAX_LINES:
            del self._lines[:-MAX_LINES]
        if lines and at_end:
            self._scroll_to_end()

    def _matches(self, t: str, line: str, err: bool, level: str = "") -> bool:
        if self.errors_only.isChecked() and not err:
            return False
        want = self.level.currentData() if hasattr(self, "level") else ""
        # Levelless output is always shown: it is not "below DEBUG", it is output that was never
        # assigned a level, and hiding it would empty the pane for everything this program prints.
        if want and level and LEVELS.index(level) < LEVELS.index(want):
            return False
        f = self.filter.text().strip().lower()
        return not f or f in line.lower()

    def _write(self, stamp: str, line: str, level: str = ""):
        """Put one line in the view, in its level's color."""
        colour = LEVEL_COLOUR.get(level)
        self.view.setTextColor(QtGui.QColor(colour) if colour
                               else self.view.palette().text().color())
        self.view.append(f"{stamp}  {line}")

    def _at_end(self) -> bool:
        bar = self.view.verticalScrollBar()
        return bar.value() >= bar.maximum() - 4

    def _scroll_to_end(self):
        # Only follow the tail if the user was already at the tail. Yanking the view to the bottom
        # while someone is reading further up is the classic log-pane annoyance.
        bar = self.view.verticalScrollBar()
        bar.setValue(bar.maximum())

    def _apply_filter(self):
        """Re-render everything kept. Called when a filter changes, not when a line arrives."""
        self.view.clear()
        for t, s, e, lv in self._lines:
            if self._matches(t, s, e, lv):
                self._write(t, s, lv)
        self._scroll_to_end()

    def text(self) -> str:
        """The currently visible log text."""
        return self.view.toPlainText()

    def copy_all(self):
        """Copy the visible lines to the clipboard."""
        cb = QtWidgets.QApplication.clipboard()
        if cb is not None:
            cb.setText(self.text())

    def clear(self):
        """Discard the log, both the pane and the buffer behind it."""
        self._lines.clear()
        self._partial = {False: "", True: ""}
        self.view.clear()

    def closeEvent(self, ev):
        self.remove()
        super().closeEvent(ev)
