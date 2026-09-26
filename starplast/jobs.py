"""Background jobs, and a record of them that outlives the job.

Everything slow in this application -- rebuilding an embedding, a hyperparameter walk, a chat reply --
used to run with no visible sign except the window going unresponsive, which is indistinguishable
from a crash. A job here is a named unit of work with a state you can look at after it has finished,
because "what did that run do" is asked more often after the fact than during.

Jobs are kept after they finish rather than discarded. A failed job that vanishes takes its traceback
with it, and the user is left knowing only that something did not work.
"""
from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from PyQt6 import QtCore

from .logging_util import get_logger
from .stopping import Stopped  # noqa: F401 -- re-exported: jobs.Stopped is stopping.Stopped

_log = get_logger(__name__)

PENDING, RUNNING, DONE, FAILED, CANCELLED = "pending", "running", "done", "failed", "cancelled"


@dataclass
class Job:
    """One unit of background work, and its state after it has finished.

    Kept after completion rather than discarded: a failed job that vanishes takes its traceback with
    it, leaving the user knowing only that something did not work.
    """
    id: int
    name: str
    state: str = PENDING
    note: str = ""
    progress: int = -1              # -1 means indeterminate, which is the honest default
    result: Any = None
    error: str = ""
    traceback: str = ""
    #: The exception itself, not only its text. A caller that wants to tell a deliberate refusal
    #: from a crash -- validation refusing a circular target, say -- needs the type, and rebuilding
    #: it from a formatted string is guesswork.
    exception: Optional[BaseException] = None
    _cancel: bool = field(default=False, repr=False)

    @property
    def active(self) -> bool:
        return self.state in (PENDING, RUNNING)

    def cancel(self) -> None:
        """Ask the job to stop. Cooperative: the callable has to check `should_stop`.

        There is no forcible kill. Terminating a thread mid-write leaves half-written caches, and
        this application writes caches.
        """
        self._cancel = True
        if self.state == PENDING:
            self.state = CANCELLED

    @property
    def cancelled(self) -> bool:
        return self._cancel


class _Signals(QtCore.QObject):
    """Signals emitted from the worker thread, delivered to the GUI thread by Qt."""
    progress = QtCore.pyqtSignal(int, int, str)     # job id, percent (-1 unknown), note
    finished = QtCore.pyqtSignal(int, bool)         # job id, ok


class _Task(QtCore.QRunnable):
    """Runs one job on the thread pool, and never lets an exception escape into Qt."""
    def __init__(self, job: Job, fn: Callable, sig: _Signals):
        super().__init__()
        self.job, self.fn, self.sig = job, fn, sig

    def run(self):
        if self.job.cancelled:
            _log.info("job %d %r: canceled before it started", self.job.id, self.job.name)
            self.sig.finished.emit(self.job.id, False)
            return
        t0 = time.monotonic()
        self.job.state = RUNNING
        _log.info("job %d %r: started", self.job.id, self.job.name)
        # Written here, on the worker, rather than left to the queued handler. The signal is
        # delivered on the GUI thread whenever it gets there, so a job that had already reported
        # "configuration 12 of 288" would have its note overwritten with "started" the moment the
        # initial emission arrived -- a progress line that goes backwards.
        self.job.note = "started"
        self.sig.progress.emit(self.job.id, -1, "started")
        try:
            # The callable is handed the job so it can report progress and honour cancellation. Passed
            # positionally only if it accepts an argument, so a plain thunk still works.
            self.job.result = self.fn(self.job) if _takes_arg(self.fn) else self.fn()
            ok = not self.job.cancelled
            self.job.state = DONE if ok else CANCELLED
            _log.info("job %d %r: %s in %.1fs", self.job.id, self.job.name,
                      "done" if ok else "cancelled", time.monotonic() - t0)
        except Stopped as exc:
            # A cooperative stop unwinds by raising, so it arrives here looking like a failure. It is
            # not one: a job the user stopped must not be reported in red with a traceback.
            self.job.state = CANCELLED
            self.job.note = str(exc)
            ok = False
            _log.info("job %d %r: stopped after %.1fs -- %s",
                      self.job.id, self.job.name, time.monotonic() - t0, exc)
        except Exception as exc:                      # a worker thread must never raise into Qt
            self.job.state = FAILED
            self.job.error = f"{type(exc).__name__}: {exc}"
            self.job.traceback = traceback.format_exc()
            self.job.exception = exc
            ok = False
            # ERROR with the traceback: a failed job keeps its traceback in the panel for this
            # session, and in the log for the ones after it.
            _log.error("job %d %r: failed after %.1fs -- %s",
                       self.job.id, self.job.name, time.monotonic() - t0, self.job.error,
                       exc_info=True)
        self.sig.finished.emit(self.job.id, ok)


def _takes_arg(fn: Callable) -> bool:
    """Whether to pass the Job. Anything uninspectable (a builtin, a C callable) is called bare."""
    try:
        import inspect
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return False
    for p in sig.parameters.values():
        if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD):
            return p.default is p.empty
        if p.kind is p.VAR_POSITIONAL:
            return True
    return False


class JobRunner(QtCore.QObject):
    """Runs jobs on a thread pool and reports on them.

    `busy_changed` exists so a caller can drive one indicator without counting jobs itself -- it
    fires only on the transitions, not on every job, so an indicator does not flicker when one job
    ends as another begins.
    """

    started = QtCore.pyqtSignal(int)
    progress = QtCore.pyqtSignal(int, int, str)
    finished = QtCore.pyqtSignal(int, bool)
    busy_changed = QtCore.pyqtSignal(bool)

    def __init__(self, parent=None, max_threads: int = 2):
        super().__init__(parent)
        self.jobs: dict[int, Job] = {}
        self._next = 0
        self._busy = False
        self._pool = QtCore.QThreadPool(self)
        # Deliberately not the core count. These jobs are memory-hungry rather than CPU-bound, and
        # two concurrent UMAP builds on this proteome is already several gigabytes.
        self._pool.setMaxThreadCount(max(1, max_threads))
        self._sig = _Signals(self)
        self._sig.progress.connect(self._on_progress)
        self._sig.finished.connect(self._on_finished)

    def submit(self, fn: Callable, name: str) -> Job:
        """Start a job on the pool and return it, so a caller can watch or cancel it."""
        self._next += 1
        job = Job(id=self._next, name=name)
        self.jobs[job.id] = job
        self.started.emit(job.id)
        self._refresh_busy()
        self._pool.start(_Task(job, fn, self._sig))
        return job

    def _on_progress(self, jid: int, pct: int, note: str):
        """Forward a worker's progress. Nothing is copied back onto the job here.

        The job's own `note` and `progress` are written where they are produced -- by the callable,
        which is handed the job for exactly that -- because this handler runs whenever Qt delivers
        the queued signal, which can be after the job has reported something newer. Assigning the
        note here made a progress line travel backwards: a job that had reached "configuration 12 of
        288" reverted to "started" the moment the initial emission was delivered.
        """
        self.progress.emit(jid, pct, note)

    def _on_finished(self, jid: int, ok: bool):
        self.finished.emit(jid, ok)
        self._refresh_busy()

    def _refresh_busy(self):
        busy = any(j.active for j in self.jobs.values())
        if busy != self._busy:
            self._busy = busy
            self.busy_changed.emit(busy)

    @property
    def busy(self) -> bool:
        return self._busy

    def active(self) -> list[Job]:
        """Jobs that have not finished."""
        return [j for j in self.jobs.values() if j.active]

    def cancel_all(self) -> None:
        """Ask every unfinished job to stop."""
        for j in self.jobs.values():
            if j.active:
                j.cancel()
        self._refresh_busy()

    def wait(self, ms: int = 5000) -> bool:
        """Block until the pool drains. For shutdown and for tests, not for the GUI thread."""
        return self._pool.waitForDone(ms)
