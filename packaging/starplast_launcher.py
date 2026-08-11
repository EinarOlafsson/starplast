"""Cross-platform launcher for the starplast browser.

Every packaging script — Windows, macOS, Debian — wraps this single entry point, so all three
installers behave identically at runtime. Mirrors spaCR's arrangement deliberately: one launcher means
a runtime bug is fixed once rather than three times.
"""
from __future__ import annotations

import multiprocessing
import os
import sys


def main() -> int:
    # A frozen bundle re-executes itself to spawn children; without this the GUI forks copies of itself.
    multiprocessing.freeze_support()
    # pyqtgraph binds to whichever Qt it finds first. A bundle can contain more than one, and the
    # resulting half-PySide half-PyQt GL widget fails at import rather than at use.
    os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
    try:
        from starplast.app import main as run
    except Exception as exc:                                   # pragma: no cover - packaging path
        sys.stderr.write(f"starplast failed to start: {exc}\n")
        return 1
    run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
