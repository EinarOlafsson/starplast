"""The assistant's subprocess path: streaming, failure, and the send/stop cycle.

Driven against a real subprocess rather than a mock, because everything here is about how a child
process actually behaves -- what arrives before it exits, what a non-zero exit means, and whether a
timeout leaves a process behind. A mocked Popen would test the mock.
"""
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtWidgets  # noqa: E402

from starplast import chat as C  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    yield QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def panels():
    """Panels and workers, torn down deterministically.

    A QThread destroyed while it is still running aborts the interpreter, and a panel dropped by the
    garbage collector takes its worker -- and the thread's parent -- with it. Every panel a test makes
    is kept alive here and stopped before the test ends.
    """
    made = []
    yield made
    for p in made:
        p.stop()
        w = p._worker
        if w is not None:
            w.wait(15000)
        p._worker = None
    QtWidgets.QApplication.processEvents()


def _python_provider(code: str, name: str = "claude") -> C.Provider:
    """A provider whose 'CLI' is this interpreter running `code`.

    argv() appends the prompt and system text, which the snippet is free to ignore -- what matters is
    that a real process starts, writes, and exits with a code we choose.
    """
    p = C.Provider(name, "test", sys.executable, "install-hint", "login-command")
    p.argv = lambda prompt, system: [sys.executable, "-c", code]
    return p


def test_stdout_is_streamed_line_by_line(qapp):
    prov = _python_provider("print('one'); print('two')")
    assert "".join(C.stream(prov, "q", "s")).splitlines() == ["one", "two"]


def test_a_non_zero_exit_reports_the_login_command(qapp):
    """A failing vendor CLI almost always means "not logged in", and the fix is a command."""
    prov = _python_provider("import sys; sys.stderr.write('not authenticated'); sys.exit(3)")
    out = "".join(C.stream(prov, "q", "s"))
    assert "exited 3" in out
    assert "not authenticated" in out
    assert "login-command" in out


def test_a_non_zero_exit_with_no_stderr_still_explains_itself(qapp):
    prov = _python_provider("import sys; sys.exit(1)")
    out = "".join(C.stream(prov, "q", "s"))
    assert "exited 1" in out and "login-command" in out


def test_output_written_before_a_failure_is_kept(qapp):
    """Partial output is evidence about what went wrong; discarding it loses the only clue."""
    prov = _python_provider("print('partial answer'); import sys; sys.exit(2)")
    out = "".join(C.stream(prov, "q", "s"))
    assert "partial answer" in out and "exited 2" in out


def test_a_hanging_process_times_out_and_is_killed(qapp):
    """Without the kill, a wedged CLI leaves a process behind for the life of the session."""
    prov = _python_provider("import time; time.sleep(30)")
    out = "".join(C.stream(prov, "q", "s", timeout=1))
    assert "timed out" in out


def test_a_process_that_cannot_start_is_reported_not_raised(qapp):
    """This runs on a worker thread feeding a widget: raising is a dead panel or a crash."""
    prov = C.Provider("x", "x", "/definitely/not/here", "hint-text", "login")
    out = "".join(C.stream(prov, "q", "s"))
    assert "not on PATH" in out and "hint-text" in out


def test_a_directory_as_the_cli_is_reported_not_raised(qapp, tmp_path):
    """FileNotFoundError is not the only way Popen fails; a directory raises PermissionError."""
    prov = C.Provider("x", "x", str(tmp_path), "hint", "login")
    out = "".join(C.stream(prov, "q", "s"))
    assert "could not start" in out or "not on PATH" in out


def test_the_worker_emits_chunks_and_finishes(qapp):
    prov = _python_provider("print('hello from the worker')")
    w = C._Worker(prov, "q", "s")
    got, done = [], []
    w.chunk.connect(got.append)
    w.done.connect(lambda: done.append(True))
    w.start()
    assert w.wait(15000)
    QtWidgets.QApplication.processEvents()
    assert done == [True]
    assert "hello from the worker" in "".join(got)


def test_stopping_the_worker_halts_delivery(qapp):
    prov = _python_provider("for i in range(500): print(i)")
    got = []
    w = C._Worker(prov, "q", "s")
    w.chunk.connect(got.append)
    w.stop()                       # stopped before it ever runs
    w.start()
    assert w.wait(15000), "the worker never finished"
    QtWidgets.QApplication.processEvents()
    assert len(got) <= 1, "a stopped worker went on delivering the whole stream"


def test_send_streams_a_reply_into_the_log(qapp, panels):
    """The whole cycle: send, stream, restore the button."""
    panel = C.ChatPanel(context_provider=lambda: "state")
    panels.append(panel)
    prov = _python_provider("print('the answer')")
    panel.current_provider = lambda: prov
    panel.input.setPlainText("a question")
    panel.send()
    assert panel._worker is not None
    assert panel.send_btn.text() == "stop"
    panel._worker.wait(15000)
    QtWidgets.QApplication.processEvents()
    panel._on_done()
    text = panel.log.toPlainText()
    assert "a question" in text and "the answer" in text
    assert panel.send_btn.text() == "send"
    assert panel._worker is None


def test_an_empty_question_is_not_sent(qapp, panels):
    panel = C.ChatPanel()
    panels.append(panel)
    panel.input.setPlainText("   ")
    panel.send()
    assert panel._worker is None


def test_a_second_question_is_refused_while_one_is_running(qapp, panels):
    """Two workers writing one reply buffer interleave into nonsense."""
    panel = C.ChatPanel()
    panels.append(panel)
    panel.current_provider = lambda: _python_provider("import time; time.sleep(2)")
    panel.input.setPlainText("first")
    panel.send()
    first = panel._worker
    panel.input.setPlainText("second")
    panel.send()
    assert panel._worker is first
    panel.stop()
    first.wait(15000)


def test_a_reply_with_no_output_says_so_rather_than_nothing(qapp, panels):
    panel = C.ChatPanel()
    panels.append(panel)
    panel._reply = ""
    panel._on_done()
    assert "(no output)" in panel.log.toPlainText()


def test_stopping_when_nothing_runs_is_harmless(qapp, panels):
    panel = C.ChatPanel()
    panels.append(panel)
    panel.stop()


def test_closing_the_panel_stops_a_running_worker(qapp, panels):
    panel = C.ChatPanel()
    panels.append(panel)
    panel.current_provider = lambda: _python_provider("import time; time.sleep(5)")
    panel.input.setPlainText("q")
    panel.send()
    w = panel._worker
    panel.close()
    assert w is None or w.isFinished(), "closing left a thread running"


def test_the_accent_falls_back_when_the_theme_is_unavailable(qapp, panels):
    """A chat panel with no window must still colour its two speakers apart."""
    panel = C.ChatPanel()
    panels.append(panel)
    assert panel._accent("you") != panel._accent("assistant")


def test_refresh_providers_disables_the_box_when_none_are_installed(qapp, panels, monkeypatch):
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    panel = C.ChatPanel()
    panels.append(panel)
    panel.refresh_providers()
    assert not panel.provider_box.isEnabled()
    assert panel.current_provider() is None
