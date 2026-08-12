"""A console pane, and the stdout capture that feeds it.

The build and the search both report progress by printing, and in a GUI that output goes nowhere --
the user watches a frozen window while the terminal they cannot see fills up. This tees stdout and
stderr into a widget without taking them away from the terminal, so running from a shell still
behaves normally and running from a launcher is no longer silent.

Capture is deliberately a tee rather than a redirect. Swallowing stderr would hide tracebacks from
anyone debugging from a terminal, which is where they are actually read.
"""
from __future__ import annotations

import sys
from datetime import datetime

from PyQt6 import QtCore, QtGui, QtWidgets

MAX_LINES = 5000        # a runaway loop must not exhaust memory through the log widget


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
        copy = QtWidgets.QPushButton("copy")
        copy.setToolTip("Copy the visible lines to the clipboard.")
        copy.clicked.connect(self.copy_all)
        clear = QtWidgets.QPushButton("clear")
        clear.setToolTip("Discard everything logged so far. The log is capped anyway, so this is "
                         "for making the next run's output easy to find rather than for memory.")
        clear.clicked.connect(self.clear)
        for w in (self.filter, self.errors_only, copy, clear):
            bar.addWidget(w)
        bar.setStretch(0, 1)
        L.addLayout(bar)

        self.view = QtWidgets.QPlainTextEdit()
        self.view.setReadOnly(True)
        self.view.setMaximumBlockCount(MAX_LINES)
        self.view.setFont(QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.SystemFont.FixedFont))
        self.view.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        L.addWidget(self.view, 1)

        self._lines: list[tuple[str, str, bool]] = []       # time, text, is_error
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
        for line in lines:
            self._lines.append((stamp, line, is_error))
        if len(self._lines) > MAX_LINES:
            del self._lines[:-MAX_LINES]
        if lines:
            self._apply_filter()

    def _matches(self, t: str, line: str, err: bool) -> bool:
        if self.errors_only.isChecked() and not err:
            return False
        f = self.filter.text().strip().lower()
        return not f or f in line.lower()

    def _apply_filter(self):
        keep = [(t, s, e) for t, s, e in self._lines if self._matches(t, s, e)]
        at_end = self.view.verticalScrollBar().value() >= self.view.verticalScrollBar().maximum() - 4
        self.view.setPlainText("\n".join(f"{t}  {s}" for t, s, _ in keep))
        if at_end:
            # Only follow the tail if the user was already at the tail. Yanking the view to the
            # bottom while someone is reading further up is the classic log-pane annoyance.
            self.view.verticalScrollBar().setValue(self.view.verticalScrollBar().maximum())

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
