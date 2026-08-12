"""Pin the Qt binding before anything imports pyqtgraph.

starplast imports PyQt6 directly. pytest-qt, if installed, imports PySide6 when the plugin loads, and
pyqtgraph then binds itself to whichever binding it finds in sys.modules first -- PySide6 -- while the app
keeps using PyQt6. The two halves of the GL widget then come from different bindings and the app fails to
import under pytest while working perfectly when launched normally.

Setting PYQTGRAPH_QT_LIB makes pyqtgraph's choice explicit instead of load-order dependent.
"""
import os

os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("PYTEST_QT_API", "PyQt6")
# The GL widget needs a platform plugin; offscreen keeps the suite headless.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


# Preferences persist through QSettings, and a test run must not write into the settings file of
# whoever is running it -- nor read one, since a developer with logging enabled would otherwise have
# every Window in the suite writing to their real log.
import tempfile  # noqa: E402

from PyQt6 import QtCore  # noqa: E402

QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat, QtCore.QSettings.Scope.UserScope,
                         tempfile.mkdtemp(prefix="starplast-test-settings-"))
