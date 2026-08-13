#!/usr/bin/env python3
"""Levels, a log file that has to be asked for, and separate control over what reaches the screen.

The console pane already tees stdout, so a walk's progress is visible while it runs. What a tee
cannot do is tell a progress line from a warning from a traceback -- so it cannot filter them, and
it keeps nothing once the window closes. Two failures in this project would have been caught in
minutes by a log and instead ran for months:

* every AlphaFold fetch returning 404, reported each time as the ordinary "no model available";
* a search completing 40 of 288 configurations with no record of which ones, or how long each took.

Both are the same shape: something that failed reporting itself the way something that succeeded
does. That is what levels are for.

## Two switches, because they answer different questions

**Writing a file is opt-in and off by default.** A tool that starts writing to someone's disk
without being asked is a tool people stop trusting, and this one is meant to be run on a laptop by
someone who did not ask for a logging system.

**What reaches the console is a separate choice, and it is never off.** DEBUG is exactly what you
want in the file during a half-hour walk and exactly what you do not want scrolling past on screen,
so the two levels are set independently. With the file disabled the console handler stays, at
WARNING: a warning that is only written to a file nobody enabled is a silent failure with extra
steps, which is the thing this module exists to end.

## The console format is parsed, not just read

Console lines are `[LEVEL] message`, which is what lets the console pane color them and filter by
level. It is a prefix rather than a suffix or a column so that a wrapped or truncated line still
begins with its level.
"""
from __future__ import annotations

import logging
import logging.handlers
import os
import sys
import time

from . import paths

#: The levels offered, coarsest last. Standard names, because a log people have to learn is a log
#: people do not read.
LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")

#: The root of this application's loggers. Everything else is a child of it, so one call configures
#: the lot and a library user of `starplast` can silence it without touching the root logger.
ROOT = "starplast"

#: How large one log file may grow, and how many are kept. A long session must not fill a disk, and
#: three backups is enough to still hold the run before the one that went wrong.
MAX_BYTES = 2_000_000
BACKUPS = 3

CONSOLE_FORMAT = "[%(levelname)s] %(message)s"
FILE_FORMAT = "%(asctime)s  %(levelname)-7s  %(name)-24s  %(message)s"

_state = {"enabled": False, "file_level": "DEBUG", "console_level": "WARNING", "path": ""}


def get_logger(name: str | None = None) -> logging.Logger:
    """The logger for one module: `get_logger(__name__)`.

    Always a child of `starplast`, whatever it is passed, so that configuring one name configures
    everything and no message can escape by being logged under an unexpected one.
    """
    if not name or name == ROOT:
        return logging.getLogger(ROOT)
    short = name.split(".")[-1]
    return logging.getLogger(f"{ROOT}.{short}")


class _LiveStreamHandler(logging.StreamHandler):
    """A stream handler that looks up `sys.stdout` when it emits, not when it is built.

    The console pane installs its tee over `sys.stdout` after this is configured. A handler holding
    the stream it was constructed with would go on writing to the real terminal and never appear in
    the pane -- which is the one place a user of the application can actually see it.
    """

    def __init__(self):
        super().__init__()

    @property
    def stream(self):
        return sys.stdout

    @stream.setter
    def stream(self, value):
        # StreamHandler's __init__ and setStream both assign here. Ignored on purpose: the stream is
        # whatever sys.stdout is at the moment of writing.
        pass


def log_dir() -> str:
    """Where log files go: under the platform cache directory, beside anything else downloaded."""
    return os.path.join(paths.user_cache_dir(), "logs")


def configure(enabled: bool = False, file_level: str = "DEBUG", console_level: str = "WARNING",
              directory: str | None = None) -> str:
    """Set logging up, and return the file being written to, or "" when none is.

    Idempotent: every call replaces the handlers rather than adding to them, so changing a level in
    Preferences twice does not print everything twice -- which is how a logging system starts being
    described as noisy and then gets turned off.

    A file that cannot be opened is reported and the rest of the configuration still applies. Losing
    the console because a directory was read-only would be a strange way to handle a log.
    """
    log = logging.getLogger(ROOT)
    for h in list(log.handlers):
        log.removeHandler(h)
        h.close()
    console_level = console_level if console_level in LEVELS else "WARNING"
    file_level = file_level if file_level in LEVELS else "DEBUG"

    console = _LiveStreamHandler()
    console.setLevel(getattr(logging, console_level))
    console.setFormatter(logging.Formatter(CONSOLE_FORMAT))
    log.addHandler(console)

    path = ""
    if enabled:
        directory = directory or log_dir()
        try:
            os.makedirs(directory, exist_ok=True)
            path = os.path.join(directory, "starplast.log")
            fh = logging.handlers.RotatingFileHandler(
                path, maxBytes=MAX_BYTES, backupCount=BACKUPS, encoding="utf8")
            fh.setLevel(getattr(logging, file_level))
            fh.setFormatter(logging.Formatter(FILE_FORMAT))
            log.addHandler(fh)
        except OSError as exc:
            path = ""
            print(f"[WARNING] logging to file is on but {directory} cannot be written: {exc}")

    # The logger passes everything the loosest handler wants; each handler then applies its own
    # level. Set on the logger rather than on the handlers alone, or a DEBUG record would be
    # dropped before the file handler ever saw it.
    log.setLevel(min(h.level for h in log.handlers) if log.handlers else logging.WARNING)
    # Not propagated to the root logger: another library's basicConfig would otherwise print every
    # one of these a second time, in its own format.
    log.propagate = False
    _state.update(enabled=bool(enabled), file_level=file_level, console_level=console_level,
                  path=path)
    log.debug("logging configured: file=%s console=%s path=%s",
              file_level if enabled else "off", console_level, path or "-")
    return path


def state() -> dict:
    """What logging is currently doing: enabled, both levels, and the file if there is one."""
    return dict(_state)


def log_file() -> str:
    """The file being written to, or "" when the log is not being kept."""
    return _state["path"]


def parse_level(line: str):
    """The level a console line was logged at, or None if it was an ordinary print.

    The console pane uses this to color and filter. Ordinary prints are not forced into a level:
    most of this application's output is `print`, and calling all of it INFO would make an INFO
    filter useless.
    """
    if line.startswith("[") and "]" in line[:12]:
        name = line[1:line.index("]")]
        if name in LEVELS:
            return name
    return None


class Timed:
    """Log how long something took, and what happened, whichever way it ends.

    A duration is the number that turns "the search finished" into something worth keeping: the
    walk that quietly did 40 of 288 configurations looked exactly like the one that did all of them
    until someone timed it.
    """

    def __init__(self, log: logging.Logger, what: str, level: int = logging.INFO):
        self.log, self.what, self.level = log, what, level
        self.t0 = 0.0

    def __enter__(self):
        self.t0 = time.monotonic()
        self.log.log(self.level, "%s: started", self.what)
        return self

    def __exit__(self, exc_type, exc, tb):
        dt = time.monotonic() - self.t0
        if exc is None:
            self.log.log(self.level, "%s: done in %.1fs", self.what, dt)
        else:
            # Not swallowed: this reports, it does not handle. The traceback goes to the file at
            # ERROR, where it can be read after the window has been closed.
            self.log.error("%s: failed after %.1fs -- %s: %s",
                           self.what, dt, type(exc).__name__, exc, exc_info=True)
        return False
