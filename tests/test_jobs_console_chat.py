"""Background jobs, the console tee, and the assistant panel.

The three exist because a browser that goes unresponsive with no output is indistinguishable from one
that has crashed. What they must not do is introduce a new way to lose information: a job that
vanishes with its traceback, a tee that swallows stderr, or a chat panel that renders a reply as
markup and eats half of it.
"""
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtCore, QtWidgets  # noqa: E402

from starplast import chat as C  # noqa: E402
from starplast import jobs as J  # noqa: E402
from starplast.console import ConsolePanel, Tee  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def _drain(runner, ms=4000):
    """Run the event loop until the pool empties, so queued signals are actually delivered."""
    assert runner.wait(ms)
    QtWidgets.QApplication.processEvents()


# --------------------------------------------------------------------------- jobs
def test_a_job_runs_and_reports_its_result(qapp):
    r = J.JobRunner()
    job = r.submit(lambda: 6 * 7, "answer")
    _drain(r)
    assert job.state == J.DONE and job.result == 42


def test_a_failed_job_keeps_its_traceback(qapp):
    """A failure that vanishes takes the only explanation with it."""
    def boom():
        raise ValueError("no")
    r = J.JobRunner()
    job = r.submit(boom, "boom")
    _drain(r)
    assert job.state == J.FAILED
    assert "ValueError: no" in job.error
    assert "boom" in job.traceback


def test_busy_changes_only_on_the_transitions(qapp):
    """An indicator driven by this must not flicker once per job."""
    r = J.JobRunner(max_threads=1)
    seen = []
    r.busy_changed.connect(seen.append)
    for i in range(3):
        r.submit(lambda: None, f"j{i}")
    _drain(r)
    assert seen == [True, False], seen


def test_a_job_that_takes_an_argument_is_handed_the_job(qapp):
    """So it can report progress and check for cancellation."""
    got = {}

    def work(job):
        got["id"] = job.id
        return "ok"
    r = J.JobRunner()
    j = r.submit(work, "with-arg")
    _drain(r)
    assert got.get("id") == j.id and j.result == "ok"


def test_a_plain_thunk_is_called_without_arguments(qapp):
    r = J.JobRunner()
    j = r.submit(lambda: "bare", "no-arg")
    _drain(r)
    assert j.result == "bare"


def test_an_uninspectable_callable_is_called_bare(qapp, monkeypatch):
    """Some C callables have no signature at all, and guessing would raise inside the worker."""
    import inspect

    def no_signature(*_a, **_k):
        raise ValueError("no signature for builtin")
    monkeypatch.setattr(inspect, "signature", no_signature)
    assert J._takes_arg(lambda job: job) is False, "an unreadable signature must fall back to bare"


def test_takes_arg_distinguishes_required_from_optional_parameters(qapp):
    assert J._takes_arg(lambda job: job) is True        # required positional: hand it the Job
    assert J._takes_arg(lambda *a: None) is True        # *args accepts it
    assert J._takes_arg(lambda x=1: None) is False      # defaulted: it is not asking for one
    assert J._takes_arg(lambda: None) is False
    assert J._takes_arg(lambda *, k=1: None) is False   # keyword-only takes no positional


def test_cancelling_before_the_job_starts_stops_it(qapp):
    r = J.JobRunner(max_threads=1)
    blocker = QtCore.QSemaphore(0)
    r.submit(lambda: blocker.acquire(1), "blocker")
    later = r.submit(lambda: "should not run", "later")
    later.cancel()
    blocker.release(1)
    _drain(r)
    assert later.result is None
    assert later.state == J.CANCELLED


def test_cancel_all_marks_everything_active(qapp):
    r = J.JobRunner(max_threads=1)
    blocker = QtCore.QSemaphore(0)
    r.submit(lambda: blocker.acquire(1), "blocker")
    r.submit(lambda: None, "queued")
    r.cancel_all()
    blocker.release(1)
    _drain(r)
    assert all(j.cancelled for j in r.jobs.values())


def test_the_task_body_runs_a_job_to_completion(qapp):
    """_Task.run() driven directly rather than through the pool.

    The pool runs it on a Qt-managed thread, which Python's trace hook does not follow, so the same
    code exercised through submit() is invisible to coverage. Calling it here traces it and tests the
    identical body -- the state transitions, the result, and the finished signal.
    """
    job = J.Job(id=1, name="direct")
    sig = J._Signals()
    seen = []
    sig.finished.connect(lambda jid, ok: seen.append((jid, ok)))
    J._Task(job, lambda: "value", sig).run()
    assert job.state == J.DONE and job.result == "value"
    assert seen == [(1, True)]


def test_the_task_body_records_a_failure_with_its_traceback(qapp):
    def boom():
        raise KeyError("missing")
    job = J.Job(id=2, name="direct-fail")
    sig = J._Signals()
    seen = []
    sig.finished.connect(lambda jid, ok: seen.append((jid, ok)))
    J._Task(job, boom, sig).run()
    assert job.state == J.FAILED
    assert "KeyError" in job.error and "boom" in job.traceback
    assert seen == [(2, False)]
    # The exception itself, not only its text. A caller that has to tell a deliberate refusal from a
    # crash -- validation refusing a circular target -- needs the type, and reconstructing it from a
    # formatted message is guesswork.
    assert isinstance(job.exception, KeyError)


def test_a_task_cancelled_before_it_starts_never_runs_its_callable(qapp):
    ran = []
    job = J.Job(id=3, name="pre-cancelled")
    job.cancel()
    sig = J._Signals()
    seen = []
    sig.finished.connect(lambda jid, ok: seen.append((jid, ok)))
    J._Task(job, lambda: ran.append(True), sig).run()
    assert ran == [], "a cancelled job still ran its callable"
    assert seen == [(3, False)]


def test_a_job_cancelled_mid_flight_is_marked_cancelled_not_done(qapp):
    """Cooperative: the callable checks the flag and returns, and the result is not a success."""
    def work(job):
        job.cancel()
        return "abandoned"
    job = J.Job(id=4, name="mid-flight")
    J._Task(job, work, J._Signals()).run()
    assert job.state == J.CANCELLED


def test_a_cooperative_stop_arrives_as_cancelled_rather_than_as_a_crash(qapp):
    """A stop unwinds by raising, so it reaches the task looking like a failure. Reported in red
    with a traceback, it would teach people to distrust the failure list.

    Driven directly rather than through the pool: a body that runs on a Qt-managed thread is
    invisible to coverage, and the point of this test is that this branch is exercised."""
    def work(job):
        raise J.Stopped("walk stopped after configuration 12")

    job = J.Job(id=7, name="stoppable")
    sig = J._Signals()
    seen = []
    sig.finished.connect(lambda jid, ok: seen.append((jid, ok)))
    J._Task(job, work, sig).run()
    assert job.state == J.CANCELLED
    assert "configuration 12" in job.note
    assert not job.error and not job.traceback, "a stop is not a failure"
    assert seen == [(7, False)]


def test_a_progress_report_does_not_overwrite_a_newer_one(qapp):
    """The initial "started" is emitted from the worker and applied whenever Qt delivers it, so
    copying it onto the job in the handler made the note travel backwards -- a walk that had reached
    configuration 12 reverted to "started"."""
    r = J.JobRunner()
    job = J.Job(id=8, name="reports")
    r.jobs[8] = job
    job.note = "configuration 12 of 288"
    r._on_progress(8, -1, "started")
    assert job.note == "configuration 12 of 288"


def test_active_lists_only_unfinished_jobs(qapp):
    r = J.JobRunner()
    r.submit(lambda: None, "done-soon")
    _drain(r)
    assert r.active() == [] and r.busy is False


# --------------------------------------------------------------------------- console
def test_the_tee_forwards_to_the_original_stream(qapp):
    """Capturing must not take output away from a terminal, where tracebacks are actually read."""
    import io
    buf = io.StringIO()
    t = Tee(buf)
    seen = []
    t.text.connect(lambda s, e: seen.append(s))
    t.write("hello")
    assert buf.getvalue() == "hello"
    assert seen == ["hello"]


def test_the_tee_survives_a_closed_original_stream(qapp):
    import io
    buf = io.StringIO()
    buf.close()
    t = Tee(buf)
    t.write("x")          # must not raise
    t.flush()


def test_the_tee_is_not_a_terminal(qapp):
    """Progress bars key off isatty and would fill the pane with escape codes."""
    assert Tee(sys.stdout).isatty() is False


def test_print_does_not_produce_two_lines(qapp):
    """print() writes the text and the newline separately."""
    p = ConsolePanel()
    p.append("one")
    p.append("\n")
    p.append("two\n")
    assert p.text().count("one") == 1
    assert "two" in p.text()


def test_an_unterminated_tail_is_held_until_its_newline(qapp):
    p = ConsolePanel()
    p.append("partial")
    assert "partial" not in p.text()
    p.append(" line\n")
    assert "partial line" in p.text()


def test_stdout_and_stderr_tails_do_not_interleave(qapp):
    """Two streams share the panel but not the buffer, or a half line of one splices into the other."""
    p = ConsolePanel()
    p.append("out-part", False)
    p.append("err-part", True)
    p.append("\n", False)
    p.append("\n", True)
    assert "out-part" in p.text() and "err-part" in p.text()
    assert "out-parterr-part" not in p.text()


def test_the_filter_hides_lines_without_discarding_them(qapp):
    p = ConsolePanel()
    for s in ("alpha\n", "beta\n"):
        p.append(s)
    p.filter.setText("alpha")
    assert "alpha" in p.text() and "beta" not in p.text()
    p.filter.setText("")
    assert "beta" in p.text()


def test_errors_only_shows_just_stderr(qapp):
    p = ConsolePanel()
    p.append("normal\n", False)
    p.append("bad\n", True)
    p.errors_only.setChecked(True)
    assert "bad" in p.text() and "normal" not in p.text()


def test_install_and_remove_restore_the_real_streams(qapp):
    p = ConsolePanel()
    real_out, real_err = sys.stdout, sys.stderr
    p.install()
    assert sys.stdout is not real_out
    p.install()                       # idempotent: a second call must not stack tees
    p.remove()
    assert sys.stdout is real_out and sys.stderr is real_err
    p.remove()                        # removing twice is not an error


def test_clear_empties_the_pane_and_the_backing_log(qapp):
    p = ConsolePanel()
    p.append("something\n")
    p.clear()
    assert p.text() == ""
    p.filter.setText("")
    assert p.text() == ""


def test_flushing_a_closed_stream_does_not_raise(qapp):
    """A closed original stream must not break the application's printing."""
    import io
    buf = io.StringIO()
    t = Tee(buf)
    buf.close()
    t.flush()                    # must not raise
    t.write("after close")       # nor this


def test_a_stream_that_raises_on_flush_is_survived(qapp):
    """A closed StringIO happens to tolerate flush(), so it does not reach the guard.

    The case the guard exists for is a stream that genuinely raises -- a pipe whose reader has gone,
    which is the normal way stdout dies when a launcher closes. Printing must not start raising
    because of it.
    """
    class Hostile:
        def write(self, _s):
            raise OSError("broken pipe")

        def flush(self):
            raise OSError("broken pipe")

    t = Tee(Hostile())
    seen = []
    t.text.connect(lambda s, e: seen.append(s))
    t.flush()                    # must not raise
    assert t.write("still emitted") == len("still emitted")
    assert seen == ["still emitted"], "output stopped reaching the pane when the stream broke"


def test_copy_puts_the_visible_lines_on_the_clipboard(qapp):
    p = ConsolePanel()
    p.append("copy me\n")
    p.copy_all()
    cb = QtWidgets.QApplication.clipboard()
    assert "copy me" in cb.text()


def test_closing_the_console_restores_the_streams(qapp):
    """A crash must not leave a dead tee holding stdout for the rest of the session."""
    p = ConsolePanel()
    real = sys.stdout
    p.install()
    assert sys.stdout is not real
    p.close()
    assert sys.stdout is real


def test_the_log_is_capped(qapp):
    """A runaway loop must not exhaust memory through the log."""
    p = ConsolePanel()
    p.append("\n".join(str(i) for i in range(6000)) + "\n")
    assert len(p._lines) <= 5000


# --------------------------------------------------------------------------- chat
def test_the_grounding_forbids_inferring_function_from_position():
    """The assistant is the easiest place in the application to produce the overclaim the rest of the
    project exists to prevent."""
    g = C.GROUNDING.lower()
    assert "do not infer" in g
    assert "negative control" in g
    assert "grey" in g and "unknown" in g


def test_the_system_prompt_carries_the_viewer_state(qapp):
    p = C.ChatPanel(context_provider=lambda: "Selected gene: TGME49_000000")
    s = p.system_prompt()
    assert "TGME49_000000" in s and C.GROUNDING in s


def test_a_broken_context_provider_does_not_block_the_question(qapp):
    def bad():
        raise RuntimeError("nope")
    p = C.ChatPanel(context_provider=bad)
    assert "context unavailable" in p.system_prompt()


def test_a_reply_containing_markup_is_escaped(qapp):
    """Unescaped, a reply mentioning a generic type would swallow the rest of the answer."""
    p = C.ChatPanel()
    p._append("assistant", "use List<int> and 5 < 6")
    assert "List<int>" in p.log.toPlainText()


def test_a_missing_cli_is_reported_with_how_to_install_it(qapp):
    prov = C.Provider("nope", "Nope", "definitely-not-on-path", "pip install nope", "nope login")
    out = "".join(C.stream(prov, "hi", "sys"))
    assert "not on PATH" in out and "pip install nope" in out


def test_each_provider_builds_a_command_carrying_the_system_prompt():
    for p in C.PROVIDERS:
        argv = p.argv("QUESTION", "SYSTEM")
        assert p.cli == argv[0]
        joined = " ".join(argv)
        assert "QUESTION" in joined and "SYSTEM" in joined


def test_available_providers_only_lists_installed_ones(monkeypatch):
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    assert C.available_providers() == []


def test_the_panel_says_what_to_install_when_nothing_is_configured(qapp, monkeypatch):
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    p = C.ChatPanel()
    p.input.setPlainText("hello")
    p.send()
    assert "No assistant CLI" in p.log.toPlainText()


def test_enter_sends_and_shift_enter_does_not(qapp):
    from PyQt6.QtGui import QKeyEvent
    box = C._Input()
    fired = []
    box.submitted.connect(lambda: fired.append(True))
    box.keyPressEvent(QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Return,
                                QtCore.Qt.KeyboardModifier.NoModifier))
    assert fired == [True]
    box.keyPressEvent(QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Return,
                                QtCore.Qt.KeyboardModifier.ShiftModifier))
    assert fired == [True], "shift+enter must insert a newline, not send"
