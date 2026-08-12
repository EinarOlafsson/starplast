#!/usr/bin/env python3
"""The download notebook, which the methods document cites as reproducing every download.

A notebook is the easiest thing in a repository to let rot: nothing imports it, no test runs it, and it
keeps opening fine while quietly referring to an API that changed. This one computed its dataset root as
`Path.cwd().parent.parent` -- the one-machine assumption paths.py was written to remove -- and fetched
with a helper written out inline rather than through the registry, so it did not benefit from the
landing-page check, the HTML-error-page check, or the checksums.

Executing it needs the network, so what is asserted is that it is valid, that it calls the API that
exists, and that it does not reintroduce the assumption.
"""
from __future__ import annotations

import json
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

NB = os.path.join(ROOT, "notebooks", "download_datasets.ipynb")


@pytest.fixture(scope="module")
def notebook():
    with open(NB, encoding="utf8") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def code(notebook):
    return "\n".join("".join(c.get("source", []))
                     for c in notebook["cells"] if c["cell_type"] == "code")


def test_the_notebook_is_valid_and_has_content(notebook):
    assert notebook.get("cells")
    assert any(c["cell_type"] == "code" for c in notebook["cells"])


def test_every_code_cell_parses_as_python(notebook):
    """A notebook that opens fine and cannot run is worse than one that does not open."""
    import ast
    for i, c in enumerate(notebook["cells"]):
        if c["cell_type"] != "code":
            continue
        src = "".join(c.get("source", []))
        if not src.strip():
            continue
        try:
            ast.parse(src)
        except SyntaxError as e:
            pytest.fail(f"cell {i} does not parse: {e}")


def test_it_resolves_the_dataset_root_rather_than_assuming_it(code):
    """The exact assumption paths.py exists to remove.

    Comments are stripped before the check rather than searched: the cell explains what it used to do,
    and a test that cannot tell an explanation from the code it warns about would forbid saying so."""
    executable = "\n".join(l for l in code.splitlines() if not l.lstrip().startswith("#"))
    assert "paths.dataset_root" in executable
    assert "parent.parent" not in executable


def test_it_fetches_through_the_registry_rather_than_its_own_helper(code):
    """An inline downloader does not get the landing-page check, the HTML-error-page check or the
    checksums, and it drifts from the fetcher the rest of the project uses."""
    assert "datasets.ensure" in code


def test_it_reports_what_cannot_be_fetched_rather_than_omitting_it(code):
    """Twelve datasets genuinely have no downloadable form. A notebook that silently skips them implies
    a complete reproduction that is not."""
    assert "datasets.fetchable" in code or "blocked" in code
    assert "datasets.missing" in code


def test_it_surfaces_the_unconfirmed_citations(code):
    """They must not reach a manuscript unchecked, and this notebook is where a reader meets the
    registry."""
    assert "datasets.unresolved" in code


def test_every_registry_function_it_calls_exists(code):
    """The failure this file is really about: an API that moved while the notebook kept opening fine."""
    from starplast import datasets, paths
    for mod, obj in (("datasets", datasets), ("paths", paths)):
        for name in sorted(set(re.findall(rf"{mod}\.(\w+)", code))):
            if name == "py":                      # a mention of paths.py in prose
                continue
            assert hasattr(obj, name), f"notebook calls {mod}.{name}, which does not exist"
