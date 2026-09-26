"""The exception that unwinds a job after a cancellation, importable without Qt.

`jobs` re-exports it, so ``except jobs.Stopped`` and ``except stopping.Stopped`` catch the same
class. It lives in its own module because the analysis modules raise it too, and they must stay
importable in a plain Python session with no Qt installed or running.
"""


class Stopped(Exception):
    """Raised inside a job to unwind it after a cancellation.

    Kept apart from the analysis panel so the runner can tell a deliberate stop from a crash without
    importing the GUI. A stopped job reported as "failed", with a traceback, teaches people to
    distrust the failure list.
    """
