#!/usr/bin/env python3
"""The installer specs, checked against where the data actually is.

An installer that ships an application which then cannot find its own cache is the same failure as a
build that cannot find its datasets, wearing a different hat -- and it is worse, because it surfaces on
someone else's machine after a download rather than on the machine that made it.

The spec cannot be executed here (PyInstaller reads it in its own context), so what is asserted is the
thing that actually went wrong: that the paths it names exist, and that they are the paths the resolver
looks in.
"""
from __future__ import annotations

import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SPEC = os.path.join(ROOT, "packaging", "starplast.spec")


@pytest.fixture(scope="module")
def spec():
    return open(SPEC, encoding="utf8").read()


def test_the_spec_bundles_the_cache_from_where_it_actually_lives(spec):
    """It pointed at <repo>/data until the cache moved into the package. After the move that directory
    stopped existing, so PyInstaller would have bundled nothing and produced an application with no
    data and no error."""
    assert 'ROOT / "starplast" / "data"' in spec
    assert '(str(ROOT / "data"), "data")' not in spec, "the pre-move path is back"


def test_the_bundled_cache_lands_where_the_resolver_looks(spec):
    """Bundled to a different relative path, the frozen application would need a special case in
    paths.py -- and a special case for frozen builds is a second resolver to keep in step."""
    assert '"starplast/data"' in spec


def test_the_cache_the_spec_names_is_present():
    from starplast import paths
    d = os.path.join(ROOT, "starplast", "data")
    assert os.path.isdir(d)
    assert os.path.dirname(paths.data_dir()) == os.path.dirname(d) or paths.data_dir() == d


def test_packaging_stops_rather_than_shipping_an_empty_application(spec):
    """A missing cache at build time must fail the build. Silently producing an installer that opens
    to an error is the worst of the options."""
    assert "raise SystemExit" in spec
    assert "build_graph" in spec, "the message has to say how to produce the cache"


def test_the_gui_toolkit_that_is_not_used_is_excluded(spec):
    """pyqtgraph binds to whichever Qt it finds first, so shipping both makes the frozen app depend on
    import order."""
    assert "PySide6" in spec and "excludes" in spec


def test_every_optional_analysis_dependency_is_named(spec):
    """umap and numba are imported lazily, so PyInstaller cannot see them by static analysis and the
    analysis panel would be dead in the installed build."""
    for pkg in ("umap", "numba", "pynndescent", "sklearn.cluster"):
        assert pkg in spec, pkg


# --------------------------------------------------------------------------- what ships
def test_the_shipped_cache_holds_only_what_the_application_reads():
    """Everything under starplast/data is carried by every install. Analysis output belongs in
    results/, where a reader can find it and a user does not have to download it."""
    from starplast import paths
    d = paths.data_dir()
    strays = [f for f in os.listdir(d)
              if f.endswith((".csv", ".log", ".txt")) or f.startswith("search_")]
    assert not strays, f"analysis output inside the shipped cache: {strays}"


def test_the_wheel_carries_the_files_the_application_needs():
    """The three identity tables were silently omitted once, because the package-data glob matched two
    extensions and did not recurse."""
    # tomllib is 3.11+; this project supports 3.10, so the line is read directly rather than adding
    # a dependency to check one setting.
    text = open(os.path.join(ROOT, "pyproject.toml"), encoding="utf8").read()
    line = next(l for l in text.splitlines() if l.strip().startswith("starplast = ["))
    globs = re.findall(r'"([^"]+)"', line)
    joined = " ".join(globs)
    for ext in (".npz", ".parquet", ".tsv"):
        assert ext in joined, f"package-data does not ship {ext}"
    assert "data/icons/*.svg" in globs
    assert "data/hf_release/*.parquet" in globs
    assert "data/*.tsv.gz" in globs
    assert not any("**" in g for g in globs), "recursive globs can ship local saved embeddings"


def test_results_are_kept_out_of_the_package():
    assert os.path.isdir(os.path.join(ROOT, "results"))
    assert os.path.exists(os.path.join(ROOT, "results", "README.md"))
