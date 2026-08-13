"""The per-dataset scripts, and the guarantee that they match the registry.

Generated files rot in a particular way: the registry gains a dataset, nobody regenerates, and the
directory quietly claims to cover 28 sources while covering 27. These tests make that a failing
test rather than something discovered months later.
"""
import os
import runpy
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
GENERATED = os.path.join(SCRIPTS, "datasets")

sys.path.insert(0, SCRIPTS)

from starplast import datasets  # noqa: E402


def test_every_registry_entry_has_a_script():
    for d in datasets.REGISTRY:
        assert os.path.exists(os.path.join(GENERATED, f"{d.key}.py")), \
            f"{d.key} has no script -- run scripts/generate_dataset_scripts.py"


def test_no_script_describes_a_dataset_that_no_longer_exists():
    """A removed dataset leaving its script behind is a file citing a source nothing uses."""
    keys = {d.key for d in datasets.REGISTRY}
    stray = [f[:-3] for f in os.listdir(GENERATED)
             if f.endswith(".py") and not f.startswith("_") and f[:-3] not in keys]
    assert not stray, f"scripts for datasets not in the registry: {stray}"


def test_regenerating_changes_nothing(tmp_path):
    """The check that keeps the two in step: a registry edit without a regeneration fails here."""
    import generate_dataset_scripts as G
    G.main(str(tmp_path), log=lambda *_: None)
    for d in datasets.REGISTRY:
        fresh = (tmp_path / f"{d.key}.py").read_text()
        committed = open(os.path.join(GENERATED, f"{d.key}.py")).read()
        assert fresh == committed, \
            f"{d.key}.py is stale -- run scripts/generate_dataset_scripts.py"


def test_every_special_names_a_real_dataset():
    """A typo in SPECIALS is silent: the entry is ignored and the dataset gets the defaults.

    Six of these were wrong when the generator was first written -- `lopit` for `lopit_tgon`,
    `proteomics` for `proteome_pru` -- so every affected script claimed the wrong normalizing module.
    """
    import generate_dataset_scripts as G
    keys = {d.key for d in datasets.REGISTRY}
    assert not set(G.SPECIALS) - keys, f"SPECIALS names datasets that do not exist: "\
                                       f"{sorted(set(G.SPECIALS) - keys)}"


def test_each_script_carries_its_provenance():
    """A script without the citation is not a record of where the data came from."""
    for d in datasets.REGISTRY:
        text = open(os.path.join(GENERATED, f"{d.key}.py")).read()
        # Several entries genuinely have no citation or no fixed path; those lines are omitted
        # rather than rendered as "None", so they are only required where they exist.
        if d.citation:
            assert d.citation in text, f"{d.key}.py does not cite its source"
        if d.path:
            assert d.path in text, f"{d.key}.py does not say where the file lives"
        assert "None" not in text.split("Fetches the source")[0].replace("normalized_by", ""), \
            f"{d.key}.py renders a missing field as the literal None"
        if d.pmid:
            assert str(d.pmid) in text
        if d.note:
            # The quirks are the expensive knowledge. First few words is enough to prove it carried.
            assert " ".join(str(d.note).split()[:6]) in " ".join(text.split())


def test_each_script_names_the_module_that_normalizes_it():
    """So nobody mistakes these for the code that produced the shipped columns."""
    for d in datasets.REGISTRY:
        text = open(os.path.join(GENERATED, f"{d.key}.py")).read()
        assert "authoritative" in text or "assembled by" in text
        assert "starplast." in text


@pytest.mark.parametrize("key", [d.key for d in datasets.REGISTRY])
def test_every_script_compiles(key):
    """Cheap, and it catches a template change that produces 28 broken files at once."""
    path = os.path.join(GENERATED, f"{key}.py")
    compile(open(path).read(), path, "exec")


def test_the_shared_helpers_import_and_resolve():
    """_common carries the accession resolver every script depends on."""
    mod = runpy.run_path(os.path.join(GENERATED, "_common.py"))
    resolve = mod["resolve"]
    assert resolve("TGME49_308090") == "TGME49_308090"
    assert resolve("not-an-accession") is None


def test_a_script_runs_end_to_end_on_a_dataset_that_is_present(tmp_path):
    """One real run, on the cell-cycle table, which is fetchable and small.

    Also a cross-check on the registry: it claims 873 phased genes, and resolution has to produce
    that many for the coverage figure quoted in the README to mean anything.
    """
    ok, _ = datasets.fetchable("xue_singlecell")
    if not ok:
        pytest.skip("xue_singlecell is not fetchable in this environment")
    proc = subprocess.run([sys.executable, "xue_singlecell.py"], cwd=GENERATED,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "873 distinct genes" in proc.stdout, proc.stdout[-2000:]
