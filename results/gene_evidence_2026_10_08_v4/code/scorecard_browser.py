"""Desktop host for the same immutable scorecard shown and exported elsewhere.

Definitions, recorded outcomes and exact exports share one supplied evaluation.
Explicit navigation keeps an archived result from opening arbitrary files or
dispatching a different cohort, while the host never fits or recalibrates it.
"""
from __future__ import annotations

from PyQt6 import QtCore, QtGui, QtWidgets

from . import scorecard_view as V


class ScorecardBrowser(QtWidgets.QTextBrowser):
    """Render a supplied scorecard and expose its definitions without refitting."""

    rows_requested = QtCore.pyqtSignal(str)
    outcome_requested = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.view = None
        self.setOpenLinks(False)
        self.setOpenExternalLinks(False)
        self.anchorClicked.connect(self.navigate)

    def set_scorecard(self, view):
        """Display exactly one immutable view, preserving its supplied numbers."""
        if not isinstance(view, V.ScorecardView):
            raise TypeError('An immutable scorecard view is required')
        self.view = view
        self.setHtml(V.render_scorecard_html(view))

    def navigate(self, destination):
        """Dispatch validated definition, expansion and outcome routes."""
        url = destination.toString() if isinstance(destination, QtCore.QUrl) else destination
        try:
            V.validate_scorecard_link(url)
        except (ValueError, TypeError):
            return False
        if QtCore.QUrl(url).scheme().lower()=='https':
            # Only explicitly supplied source links can leave this card.
            if self.view and QtCore.QUrl(url).toString() in {QtCore.QUrl(link.url).toString() for link in self.view.links}:
                return QtGui.QDesktopServices.openUrl(QtCore.QUrl(url))
            return False
        if self.view is None:
            return False
        if url == 'scorecard:expand':
            self.setHtml(V.render_scorecard_html(self.view, expanded=True))
            return True
        route, key = url.removeprefix('scorecard:').split('/', 1)
        if route in {'metric', 'detail'}:
            keys={metric.key for metric in self.view.metrics} if route=='metric' else {detail.key for detail in self.view.details}
            if key not in keys:return False
            try:
                html = V.render_scorecard_detail(self.view, key)
            except KeyError:
                return False
            self.setHtml(html + '<p><a href="scorecard:expand">Full test record</a></p>')
        elif route == 'rows':
            if url not in {link.url for link in self.view.links}:return False
            self.rows_requested.emit(key)
        elif route == 'outcome':
            if url not in {link.url for link in self.view.links}:return False
            self.outcome_requested.emit(key)
        return True

    def export_json(self):
        """Return the exact displayed model and original supplied evaluation."""
        if self.view is None:
            raise ValueError('No scorecard is selected')
        return V.export_scorecard(self.view)

    def clear(self):
        """Clear the old model along with its HTML so it cannot be exported."""
        self.view = None
        super().clear()
