#!/usr/bin/env python3
"""Levels, the opt-in log file, and the console pane's reading of them.

The two failures this exists to prevent are both silent ones this project has already had: every
AlphaFold fetch returning 404 while reporting the ordinary "no model available", and a search
finishing 40 of 288 configurations with no record of which. Both are something that failed reporting
itself the way something that succeeded does, so these tests are mostly about the failure paths --
that a warning is emitted at all, that it survives with its URL attached, and that nothing is written
to disk unless it was asked for.
"""
from __future__ import annotations

import logging
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import logging_util as L  # noqa: E402


@pytest.fixture(autouse=True)
def restore():
    """Leave the logger as it was found: this configures a process-wide logger."""
    log = logging.getLogger(L.ROOT)
    before, level, prop = list(log.handlers), log.level, log.propagate
    yield
    for h in list(log.handlers):
        log.removeHandler(h)
    for h in before:
        log.addHandler(h)
    log.setLevel(level)
    log.propagate = prop


# --------------------------------------------------------------------------- opt-in
def test_nothing_is_written_to_disk_unless_it_was_asked_for(tmp_path):
    """A tool that starts writing files to someone's disk without being asked is a tool people stop
    trusting -- and this one is meant to be run by someone who did not ask for a logging system."""
    assert L.configure(enabled=False, directory=str(tmp_path)) == ""
    L.get_logger("x").error("something went wrong")
    assert os.listdir(tmp_path) == []


def test_warnings_still_reach_the_console_with_the_file_off(capsys):
    """A warning that is only written to a file nobody enabled is a silent failure with extra
    steps, which is the thing this module exists to end."""
    L.configure(enabled=False)
    L.get_logger("x").warning("the fetch returned 404")
    out = capsys.readouterr().out
    assert "[WARNING] the fetch returned 404" in out


def test_the_console_level_is_independent_of_the_file(tmp_path, capsys):
    """DEBUG is exactly what you want kept during a half-hour walk and exactly what you do not want
    scrolling past while you watch it."""
    path = L.configure(enabled=True, file_level="DEBUG", console_level="ERROR",
                       directory=str(tmp_path))
    log = L.get_logger("x")
    log.debug("configuration 12 of 288")
    log.error("it fell over")
    out = capsys.readouterr().out
    assert "configuration 12" not in out and "[ERROR] it fell over" in out
    kept = open(path).read()
    assert "configuration 12 of 288" in kept and "it fell over" in kept


def test_the_file_records_the_level_the_time_and_the_module(tmp_path):
    path = L.configure(enabled=True, file_level="INFO", directory=str(tmp_path))
    L.get_logger("starplast.jobs").info("job 3 'walk': done in 12.0s")
    line = open(path).read().strip().splitlines()[-1]
    assert "INFO" in line and "jobs" in line and "done in 12.0s" in line
    assert line[:2].isdigit(), "no timestamp"


def test_configuring_twice_does_not_print_everything_twice(tmp_path, capsys):
    """Handlers that accumulate are how a logging system starts being called noisy and then gets
    turned off."""
    for _ in range(3):
        L.configure(enabled=True, console_level="INFO", directory=str(tmp_path))
    capsys.readouterr()
    L.get_logger("x").info("once")
    assert capsys.readouterr().out.count("once") == 1


def test_a_directory_that_cannot_be_written_costs_the_file_not_the_console(tmp_path, capsys):
    """Losing the console because a directory was read-only would be a strange way to handle a
    log."""
    blocked = tmp_path / "wall"
    blocked.write_text("not a directory")
    assert L.configure(enabled=True, directory=str(blocked / "logs")) == ""
    assert "cannot be written" in capsys.readouterr().out
    L.get_logger("x").warning("still audible")
    assert "still audible" in capsys.readouterr().out


def test_an_unknown_level_falls_back_rather_than_raising(tmp_path):
    """A stored preference from a future version must not stop the application from starting."""
    L.configure(enabled=True, file_level="LOUD", console_level="QUIET", directory=str(tmp_path))
    assert L.state()["file_level"] == "DEBUG" and L.state()["console_level"] == "WARNING"


def test_the_state_says_what_logging_is_doing(tmp_path):
    L.configure(enabled=True, file_level="INFO", console_level="ERROR", directory=str(tmp_path))
    s = L.state()
    assert s["enabled"] and s["file_level"] == "INFO" and s["console_level"] == "ERROR"
    assert s["path"] == L.log_file() and os.path.exists(s["path"])


def test_the_file_rotates_so_a_long_session_cannot_fill_a_disk(tmp_path, monkeypatch):
    monkeypatch.setattr(L, "MAX_BYTES", 2000)
    path = L.configure(enabled=True, file_level="DEBUG", directory=str(tmp_path))
    log = L.get_logger("x")
    for i in range(400):
        log.debug("a line that is not especially short, number %d", i)
    files = sorted(os.listdir(tmp_path))
    assert len(files) > 1, files
    assert os.path.getsize(path) <= 4000


def test_every_logger_is_a_child_of_one_root(tmp_path):
    """Configuring one name configures the lot, and no message escapes by being logged under an
    unexpected one."""
    assert L.get_logger("starplast.app").name == "starplast.app"
    assert L.get_logger("some.other.module").name == "starplast.module"
    assert L.get_logger(None).name == L.ROOT
    assert L.get_logger(L.ROOT).name == L.ROOT


def test_records_do_not_reach_the_root_logger(tmp_path, capsys):
    """Another library's basicConfig would otherwise print every one of these a second time, in its
    own format."""
    L.configure(enabled=False, console_level="INFO")
    assert logging.getLogger(L.ROOT).propagate is False


def test_the_default_log_directory_is_the_cache_not_the_working_directory():
    """A log dropped into whatever directory the app was launched from is litter."""
    from starplast import paths
    assert L.log_dir().startswith(paths.user_cache_dir())


# --------------------------------------------------------------------------- timing
def test_a_timed_block_records_how_long_it_took(tmp_path):
    path = L.configure(enabled=True, file_level="INFO", directory=str(tmp_path))
    with L.Timed(L.get_logger("x"), "the search"):
        pass
    text = open(path).read()
    assert "the search: started" in text and "the search: done in" in text


def test_a_timed_block_that_fails_says_so_and_re_raises(tmp_path):
    """It reports; it does not handle. A duration that only appears on success is exactly the case
    where the number matters least."""
    path = L.configure(enabled=True, file_level="INFO", directory=str(tmp_path))
    with pytest.raises(KeyError):
        with L.Timed(L.get_logger("x"), "the fetch"):
            raise KeyError("gone")
    text = open(path).read()
    assert "the fetch: failed after" in text and "KeyError" in text
    assert "Traceback" in text, "the traceback is the reason to keep the file"


# --------------------------------------------------------------------------- parsing for the pane
@pytest.mark.parametrize("line,expected", [
    ("[WARNING] the fetch returned 404", "WARNING"),
    ("[DEBUG] configuration 12", "DEBUG"),
    ("walk_umap: 2,000 of 8,140 genes", None),
    ("[NOTALEVEL] something", None),
    ("", None),
    ("[WARNING", None),
])
def test_a_console_line_reports_the_level_it_was_logged_at(line, expected):
    """Ordinary prints are left levelless rather than being called INFO: most output here is a
    print, and folding it into a level would make filtering by that level useless."""
    assert L.parse_level(line) == expected


# --------------------------------------------------------------------------- the console pane
@pytest.fixture(scope="module")
def qapp():
    """Held for the module's lifetime. Created and dropped inside a test, the QApplication is
    garbage-collected before the widget is built and Qt aborts the interpreter."""
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def pane(qapp):
    from starplast.console import ConsolePanel
    return ConsolePanel()


def test_the_pane_hides_what_is_below_the_chosen_level(pane):
    """A walk prints thousands of progress lines and one warning; without this the warning is
    somewhere in the middle of them."""
    pane.append("[DEBUG] configuration 12\n")
    pane.append("[WARNING] the fetch returned 404\n")
    pane.level.setCurrentIndex(pane.level.findData("WARNING"))
    assert "404" in pane.text() and "configuration 12" not in pane.text()
    pane.level.setCurrentIndex(pane.level.findData(""))
    assert "configuration 12" in pane.text()


def test_ordinary_output_is_never_hidden_by_a_level(pane):
    """It is not "below DEBUG"; it is output that was never assigned a level, and hiding it would
    empty the pane of everything this program prints."""
    pane.append("walk_umap: 2,000 of 8,140 genes\n")
    pane.append("[DEBUG] noise\n")
    pane.level.setCurrentIndex(pane.level.findData("ERROR"))
    assert "walk_umap" in pane.text() and "noise" not in pane.text()


def test_a_level_colors_its_own_line(pane):
    from starplast.console import LEVEL_COLOR
    from PyQt6 import QtGui
    pane.append("[ERROR] it fell over\n")
    doc = pane.view.document()
    block = doc.lastBlock()
    fmt = block.begin().fragment().charFormat()
    assert fmt.foreground().color().name() == QtGui.QColor(LEVEL_COLOR["ERROR"]).name()


def test_levelless_lines_keep_the_pane_s_own_color(pane):
    from starplast.console import LEVEL_COLOR
    pane.append("just a print\n")
    fmt = pane.view.document().lastBlock().begin().fragment().charFormat()
    assert fmt.foreground().color().name() not in {c.lower() for c in LEVEL_COLOR.values()}


def test_the_level_filter_and_the_text_filter_both_apply(pane):
    pane.append("[WARNING] alpha failed\n")
    pane.append("[WARNING] beta failed\n")
    pane.level.setCurrentIndex(pane.level.findData("WARNING"))
    pane.filter.setText("alpha")
    assert "alpha" in pane.text() and "beta" not in pane.text()


def test_a_line_arriving_while_filtered_out_is_kept_underneath(pane):
    """The filter hides; it must not discard, or clearing it would lose the run."""
    pane.level.setCurrentIndex(pane.level.findData("ERROR"))
    pane.append("[INFO] quiet progress\n")
    assert "quiet progress" not in pane.text()
    pane.level.setCurrentIndex(pane.level.findData(""))
    assert "quiet progress" in pane.text()


def test_the_tee_can_be_closed_without_closing_the_stream_it_wraps(qapp):
    """Python's logging shutdown calls close() on whatever it finds at sys.stderr. Closing the real
    stream instead would take stderr down for everything still running -- and with a log file now in
    the picture, shutdown reaches this on every exit rather than occasionally."""
    import io
    from starplast.console import Tee
    buf = io.StringIO()
    t = Tee(buf)
    t.close()
    assert not buf.closed
    t.write("still working\n")
    assert "still working" in buf.getvalue()


def test_settings_are_isolated_from_the_user_running_the_suite(qapp):
    """This failed once, silently and expensively: a test's temporary directory was written into
    the real settings file, and the application then started with a log directory that no longer
    existed and a console level nobody had chosen.

    On Unix NativeFormat IS the ini format, so a QSettings built from (organisation, application)
    resolves as NativeFormat whatever setDefaultFormat says -- setting the path for IniFormat alone
    left every test writing to ~/.config."""
    from PyQt6 import QtCore
    s = QtCore.QSettings("starplast", "starplast")
    assert "starplast-test-settings" in s.fileName(), s.fileName()
    s.setValue("logging/console_level", "DEBUG")
    s.sync()
    assert os.path.exists(s.fileName())
    real = os.path.expanduser("~/.config/starplast/starplast.conf")
    assert not os.path.exists(real) or "test-settings" not in open(real).read()


def test_the_user_s_saved_runs_are_isolated_from_the_suite(qapp):
    """The same failure as the settings one, found the same way -- by looking at a real machine. A
    user's store held `test_run_a`, `test_two_b` and `before` among their own clusterings, because
    every Window the suite builds keeps its runs in `paths.user_cache_dir()`. A test that writes
    into what it is testing has damaged the thing it was checking.

    Saved embeddings travel with it: they land inside the package, so the suite was leaving walk
    output in the shipped directory too.
    """
    import numpy as np
    from starplast import paths
    from starplast.app import Window

    assert "starplast-test-state" in paths.user_cache_dir(), paths.user_cache_dir()
    win = Window()
    try:
        win.keep_run(np.zeros(win.n, dtype=int), name="isolation_check")
        assert "starplast-test-state" in win.runs.root
        assert "starplast-test-state" in win.annotations.path
        assert "starplast-test-state" in win.panel.store.root
        assert os.path.exists(os.path.join(win.runs.root, "isolation_check.npz"))
    finally:
        win.close()
    real = os.path.join(os.path.expanduser("~/.cache/starplast"), "runs")
    assert not os.path.exists(os.path.join(real, "isolation_check.npz")), \
        "the suite wrote a run into the store of whoever ran it"
