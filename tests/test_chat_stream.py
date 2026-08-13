"""The assistant's subprocess path: streaming, failure, and the send/stop cycle.

Driven against a real subprocess rather than a mock, because everything here is about how a child
process actually behaves -- what arrives before it exits, what a non-zero exit means, and whether a
timeout leaves a process behind. A mocked Popen would test the mock.
"""
import os
import subprocess
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


def test_register_receives_the_process_as_soon_as_it_exists(qapp):
    """So a caller can kill it. Without this, stop can only set a flag checked between chunks."""
    got = []
    prov = _python_provider("print('x')")
    list(C.stream(prov, "q", "s", register=got.append))
    assert got and hasattr(got[0], "poll")
    assert got[0].poll() is not None, "the process was left running"


def test_the_worker_body_runs_without_a_thread(qapp):
    """_Worker.run() driven directly.

    Through start() it runs on a Qt-managed thread, which Python's trace hook does not follow, so the
    same body exercised by every other test here is invisible to coverage. This is the identical code.
    """
    prov = _python_provider("print('direct body')")
    w = C._Worker(prov, "q", "s")
    got, done = [], []
    w.chunk.connect(got.append)
    w.done.connect(lambda: done.append(True))
    w.run()
    assert done == [True]
    assert "direct body" in "".join(got)


def test_the_worker_body_stops_delivering_once_stopped(qapp):
    prov = _python_provider("for i in range(200): print(i)")
    w = C._Worker(prov, "q", "s")
    got = []
    w.chunk.connect(got.append)
    w.stop()
    w.run()
    assert len(got) <= 1


def test_registering_after_a_stop_kills_the_process_immediately(qapp):
    """The race that matters: someone closes the window between start() and the process existing.

    Left unhandled, _proc is still None when stop() looks, the kill is skipped, and the thread reads
    a process nobody is waiting for -- then Qt aborts when the panel is destroyed under it.
    """
    prov = _python_provider("import time; time.sleep(30)")
    w = C._Worker(prov, "q", "s")
    w.stop()                                  # stop arrives first
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        w._register(proc)                     # ...then the process appears
        proc.wait(timeout=10)
        assert proc.poll() is not None, "a process registered after a stop was left running"
    finally:
        if proc.poll() is None:
            proc.kill()


def test_killing_an_already_dead_process_is_harmless(qapp):
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait(timeout=10)
    C._Worker._kill(proc)                     # already exited
    C._Worker._kill(None)                     # never started


def test_a_process_that_refuses_to_die_does_not_raise(qapp):
    """kill() can raise if the process is already reaped by something else."""
    class Stubborn:
        def poll(self):
            return None

        def kill(self):
            raise OSError("no such process")

    C._Worker._kill(Stubborn())               # must not raise


def test_a_reader_that_blows_up_still_ends_the_stream(qapp, monkeypatch):
    """The reader thread's guard. A pipe closing under a kill raises mid-read, and if that killed the
    thread without queueing its sentinel the generator would wait on a line that never comes."""
    real_popen = C.subprocess.Popen

    class Exploding:
        def readline(self):
            raise OSError("pipe torn down")

        def close(self):
            pass

    def popen(*a, **k):
        proc = real_popen(*a, **k)
        proc.stdout = Exploding()
        return proc

    monkeypatch.setattr(C.subprocess, "Popen", popen)
    prov = _python_provider("print('never read')")
    out = "".join(C.stream(prov, "q", "s", timeout=10))
    assert "timed out" not in out, "a broken reader was left to hit the deadline instead of ending"


def test_a_process_that_will_not_exit_after_its_output_is_killed(qapp):
    """Covers the wait() timeout inside the normal path: stdout closes but the process lingers.

    Without the kill here the child outlives the answer, which over a session is a pile of stranded
    processes rather than one visible failure.
    """
    # os.close(1), not sys.stdout.close(): closing Python's wrapper leaves the underlying descriptor
    # open until the process exits, so the parent never sees EOF and falls through to the outer
    # deadline instead -- which is a different branch and leaves this one untested.
    prov = _python_provider(
        "import os, sys, time; print('done talking'); sys.stdout.flush(); os.close(1); "
        "time.sleep(30)")
    out = "".join(C.stream(prov, "q", "s", timeout=3))
    assert "done talking" in out
    assert "timed out" not in out, "it hit the read deadline rather than the wait after EOF"


def test_closing_a_stream_that_refuses_to_close_is_survived(qapp, monkeypatch):
    """The final guard. Tearing down the pipes must not raise after a perfectly good answer."""
    real_popen = C.subprocess.Popen

    class Unclosable:
        def __init__(self, inner):
            self._inner = inner

        def __getattr__(self, name):
            return getattr(self._inner, name)

        def close(self):
            raise OSError("cannot close")

    def popen(*a, **k):
        proc = real_popen(*a, **k)
        proc.stdout = Unclosable(proc.stdout)
        proc.stderr = Unclosable(proc.stderr)
        return proc

    monkeypatch.setattr(C.subprocess, "Popen", popen)
    prov = _python_provider("print('answered anyway')")
    assert "answered anyway" in "".join(C.stream(prov, "q", "s"))


def test_the_accent_falls_back_when_the_palette_cannot_be_read(qapp, panels, monkeypatch):
    """A panel with no window, or a theme module that moved, must still color the two speakers apart."""
    panel = C.ChatPanel()
    panels.append(panel)
    monkeypatch.setattr(panel, "window", lambda: (_ for _ in ()).throw(RuntimeError("no window")))
    assert panel._accent("you") == "#7aa2f7"
    assert panel._accent("assistant") == "#9ece6a"


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
    """A chat panel with no window must still color its two speakers apart."""
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
