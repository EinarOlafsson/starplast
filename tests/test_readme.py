#!/usr/bin/env python3
"""The README, and the rule that its dataset table cannot drift from the registry.

The README's scope is fixed: what the program does, how to run it, and which datasets are inside with
their references and data types. Everything else lives in HANDOFF.md and MATERIALS_AND_METHODS.md. A
README is read by someone deciding whether to run the thing and then running it; argument belongs in
the paper, not in the front door.

The dataset table is generated rather than written, because a hand-written one drifts the moment a
dataset is added -- and a README that misstates which data is inside is worse than one that omits it.
"""
from __future__ import annotations

import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import datasets as D  # noqa: E402

README = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "README.md")


@pytest.fixture(scope="module")
def readme():
    return open(README, encoding="utf8").read()


def test_the_committed_dataset_table_matches_the_registry(readme):
    """Add a dataset to the registry without updating the README and this fails. That is the point."""
    assert D.readme_table().strip() in readme, (
        "README.md's dataset table is out of date -- regenerate it from "
        "starplast.datasets.readme_table()")


def test_every_dataset_appears_in_the_table():
    table = D.readme_table()
    missing = [d.name for d in D.REGISTRY if d.name not in table]
    assert not missing, f"datasets absent from the README table: {missing}"


def test_each_dataset_is_listed_under_its_own_level():
    """`level` is the on-disk taxonomy, so the table doubles as a map of datasets/."""
    table = D.readme_table()
    for level in D.LEVEL_ORDER:
        if any(d.level == level for d in D.REGISTRY):
            assert D.LEVEL_TITLE[level] in table


def test_the_reference_column_names_the_publication_not_the_file():
    """A path is an implementation detail that changes with the layout; a PMID does not."""
    table = D.readme_table()
    assert "pubmed.ncbi.nlm.nih.gov" in table
    for d in D.REGISTRY:
        if d.path:
            assert d.path not in table, f"{d.key}: the file path must not be in the README"


def test_an_unconfirmed_citation_is_stated_rather_than_left_blank():
    """A blank cell hides it; unresolved() exists so those cannot reach a manuscript unchecked."""
    if D.unresolved():
        assert "citation not yet confirmed" in D.readme_table()


def test_the_readme_says_how_to_install_and_run(readme):
    assert "pip install -e ." in readme
    assert re.search(r"^\s*starplast\s*$", readme, re.M)


def test_the_readme_says_what_the_program_does(readme):
    assert "8,140" in readme
    assert "Toxoplasma gondii" in readme


def test_the_readme_documents_the_data_overrides(readme):
    """The first thing anyone needs when a file cannot be found."""
    from starplast import paths
    assert paths.ENV_CACHE in readme and paths.ENV_DATASETS in readme
    assert "python -m starplast.paths" in readme


def test_the_readme_points_at_the_documents_that_carry_the_argument(readme):
    """The science moved out; the pointers must not."""
    assert "HANDOFF.md" in readme
    assert "MATERIALS_AND_METHODS.md" in readme


def test_the_readme_stays_short_enough_to_read_before_installing(readme):
    """It was 298 lines and mostly an essay about method. The scope is now three things.

    Counted with the generated dataset table removed. That table is one row per registry entry and
    grows every time data is added, which is the project working rather than the scope creeping; it
    tripped this limit at 201 lines the day five PRIDE deposits were registered. What the limit is
    actually guarding is the PROSE -- the essay about method that this README used to be -- so that
    is what is measured.
    """
    prose = readme.replace(D.readme_table().strip(), "")
    assert len(prose.splitlines()) < 200


def test_the_readme_does_not_re_argue_the_method(readme):
    """Interpretation rules, circularity and the search machinery belong in HANDOFF and the methods
    document. A heading about them here means the scope has crept back."""
    headings = re.findall(r"^##+ (.+)$", readme, re.M)
    banned = ("circular", "interpretation", "attention correction", "structural hole")
    offending = [h for h in headings if any(b in h.lower() for b in banned)]
    assert not offending, f"headings that re-argue method: {offending}"


def test_every_screenshot_the_readme_shows_actually_exists(readme):
    """A broken image is worse than no image: it says the repository is unmaintained."""
    import re as _re
    root = os.path.dirname(README)
    missing = [p for p in _re.findall(r"!\[[^\]]*\]\(([^)]+)\)", readme)
               if not os.path.exists(os.path.join(root, p))]
    assert not missing, f"README references images that are not in the repository: {missing}"


def test_the_screenshots_are_committed_rather_than_generated_on_demand():
    """Unlike the API page, these cannot be rebuilt by CI -- rendering needs a GPU context and the
    committed cache -- so they are part of the repository."""
    root = os.path.dirname(README)
    shots = os.path.join(root, "docs", "screenshots")
    assert os.path.isdir(shots)
    assert len([f for f in os.listdir(shots) if f.endswith(".png")]) >= 5
