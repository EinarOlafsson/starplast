#!/usr/bin/env python3
"""The API reference can actually be generated.

pdoc imports every module to read it, so this catches two things that otherwise surface only in CI: a
module that cannot be imported headlessly, and a docstring whose markup breaks the renderer. Both are
easy to introduce and invisible locally, because nothing else in the suite imports every module at once
purely for its documentation.
"""
from __future__ import annotations

import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

pdoc = pytest.importorskip("pdoc", reason="pdoc is only needed to build the API page")


def test_every_module_carries_a_docstring():
    """These are the reference -- there is no second prose document -- so an undocumented module is a
    hole in the published page, not just an untidy file."""
    import importlib
    import pkgutil

    import starplast
    missing = []
    for m in pkgutil.iter_modules(starplast.__path__):
        if m.name.startswith("_"):
            continue
        mod = importlib.import_module(f"starplast.{m.name}")
        if not (mod.__doc__ or "").strip():
            missing.append(m.name)
    assert not missing, f"modules with no docstring: {missing}"


def test_module_docstrings_say_why_rather_than_only_what():
    """The reasoning is what makes these worth publishing. A one-line summary is a signature restated."""
    import importlib
    import pkgutil

    import starplast
    thin = []
    for m in pkgutil.iter_modules(starplast.__path__):
        if m.name.startswith("_"):
            continue
        doc = (importlib.import_module(f"starplast.{m.name}").__doc__ or "")
        if len(doc.split()) < 40:
            thin.append((m.name, len(doc.split())))
    assert not thin, f"modules whose docstring explains nothing: {thin}"


@pytest.mark.slow
def test_the_api_reference_builds(tmp_path):
    """pdoc imports every module, so this fails on anything that cannot be imported headlessly."""
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", PYQTGRAPH_QT_LIB="PyQt6")
    r = subprocess.run([sys.executable, "-m", "pdoc", "--output-directory", str(tmp_path),
                        "--no-search", "--docformat", "markdown", "starplast"],
                       cwd=ROOT, capture_output=True, env=env, timeout=600)
    assert r.returncode == 0, r.stderr.decode()[-2000:]
    produced = os.listdir(os.path.join(tmp_path, "starplast"))
    assert len(produced) > 20, produced
