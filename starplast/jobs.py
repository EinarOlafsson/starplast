"""Background jobs, and a record of them that outlives the job.

Everything slow in this application -- rebuilding an embedding, a hyperparameter walk, a chat reply --
used to run with no visible sign except the window going unresponsive, which is indistinguishable
from a crash. A job here is a named unit of work with a state you can look at after it has finished,
because "what did that run do" is asked more often after the fact than during.

Jobs are kept after they finish rather than discarded. A failed job that vanishes takes its traceback
with it, and the user is left knowing only that something did not work.
"""
from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

from PyQt6 import QtCore

PENDING, RUNNING, DONE, FAILED, CANCELLED = "pending", "running", "done", "failed", "cancelled"


class Stopped(Exception):
    """Raised inside a job to unwind it after a cancellation.

    Lives here rather than beside the analysis panel so the runner can tell a deliberate stop from a
    crash without importing the GUI. A stopped job reported as "failed", with a traceback, teaches
    people to distrust the failure list.
    """


@dataclass
class Job:
    id: int
    name: str
    state: str = PENDING
    note: str = ""
    progress: int = -1              # -1 means indeterminate, which is the honest default
    result: Any = None
    error: str = ""
    traceback: str = ""
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
    progress = QtCore.pyqtSignal(int, int, str)     # job id, percent (-1 unknown), note
    finished = QtCore.pyqtSignal(int, bool)         # job id, ok


class _Task(QtCore.QRunnable):
    def __init__(self, job: Job, fn: Callable, sig: _Signals):
        super().__init__()
        self.job, self.fn, self.sig = job, fn, sig

    def run(self):
        if self.job.cancelled:
            self.sig.finished.emit(self.job.id, False)
            return
        self.job.state = RUNNING
        self.sig.progress.emit(self.job.id, -1, "started")
        try:
            # The callable is handed the job so it can report progress and honour cancellation. Passed
            # positionally only if it accepts an argument, so a plain thunk still works.
            self.job.result = self.fn(self.job) if _takes_arg(self.fn) else self.fn()
            ok = not self.job.cancelled
            self.job.state = DONE if ok else CANCELLED
        except Stopped as exc:
            # A cooperative stop unwinds by raising, so it arrives here looking like a failure. It is
            # not one: a job the user stopped must not be reported in red with a traceback.
            self.job.state = CANCELLED
            self.job.note = str(exc)
            ok = False
        except Exception as exc:                      # a worker thread must never raise into Qt
            self.job.state = FAILED
            self.job.error = f"{type(exc).__name__}: {exc}"
            self.job.traceback = traceback.format_exc()
            ok = False
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
        self._next += 1
        job = Job(id=self._next, name=name)
        self.jobs[job.id] = job
        self.started.emit(job.id)
        self._refresh_busy()
        self._pool.start(_Task(job, fn, self._sig))
        return job

    def _on_progress(self, jid: int, pct: int, note: str):
        j = self.jobs.get(jid)
        if j is not None:
            j.progress, j.note = pct, note
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
        return [j for j in self.jobs.values() if j.active]

    def cancel_all(self) -> None:
        for j in self.jobs.values():
            if j.active:
                j.cancel()
        self._refresh_busy()

    def wait(self, ms: int = 5000) -> bool:
        """Block until the pool drains. For shutdown and for tests, not for the GUI thread."""
        return self._pool.waitForDone(ms)
