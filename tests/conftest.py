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
# And a SOFTWARE GL context, which the suite has always used by accident and now uses on purpose.
# pyqtgraph 0.14 parses the driver's version string as a float, so an NVIDIA context reporting
# `4.6.0 NVIDIA 580.173.02` raises "Requires >= OpenGL 2.1" on a driver that plainly exceeds it --
# and every GL test then fails on a machine that has a card, while passing on one that does not.
# Software GL also makes the render comparisons reproducible, which is what they are for.
os.environ.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
# Mesa reads the variable above; a machine whose GLX vendor library is NVIDIA's does not, so Qt is
# told directly as well. Both, because either one alone leaves the suite passing or failing
# according to which graphics stack the machine happens to have.
os.environ.setdefault("QT_OPENGL", "software")
os.environ.setdefault("__GLX_VENDOR_LIBRARY_NAME", "mesa")


# Preferences persist through QSettings, and a test run must not write into the settings file of
# whoever is running it -- nor read one, since a developer with logging enabled would otherwise have
# every Window in the suite writing to their real log.
#
# The path is set for BOTH formats, which is the part that matters and that a first attempt got
# wrong. On Unix NativeFormat IS the ini format, so a QSettings built from (organisation,
# application) resolves as NativeFormat whatever `setDefaultFormat` says, and setting the path for
# IniFormat alone left every test writing to ~/.config. It did: a test's temporary directory ended
# up in the real settings file and the application then started with a log directory that no longer
# existed. `test_settings_are_isolated_from_the_user_running_the_suite` is the check that this holds.
import tempfile  # noqa: E402

from PyQt6 import QtCore  # noqa: E402

SETTINGS_DIR = tempfile.mkdtemp(prefix="starplast-test-settings-")
QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
for _fmt in (QtCore.QSettings.Format.NativeFormat, QtCore.QSettings.Format.IniFormat):
    QtCore.QSettings.setPath(_fmt, QtCore.QSettings.Scope.UserScope, SETTINGS_DIR)
    QtCore.QSettings.setPath(_fmt, QtCore.QSettings.Scope.SystemScope, SETTINGS_DIR)

# The same rule for the user's own state, and for the same reason. Every `Window` the suite builds
# keeps its clusterings in `paths.user_cache_dir()/runs` and its annotations beside them, so the
# tests were writing into the saved runs of whoever ran them: a real store on this machine had
# `test_run_a`, `test_two_b` and `before` sitting among a user's own work, and any test that keeps a
# run left another one behind. Saved embeddings go the same way, into the package's own directory.
#
# `STARPLAST_STATE` points all three somewhere temporary for the duration of the run.
# `test_the_user_s_saved_runs_are_isolated_from_the_suite` is the check that this holds.
STATE_DIR = tempfile.mkdtemp(prefix="starplast-test-state-")
os.environ["STARPLAST_STATE"] = STATE_DIR
