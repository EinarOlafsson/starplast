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


# --------------------------------------------------------------------------- the slot table
def test_the_slot_table_regenerates_identically(tmp_path, monkeypatch):
    """Both files come from one list in the generator, and the coverage numbers are counted from the
    cache on every run. A table checked in by hand goes stale the first time a dataset lands."""
    import runpy
    import shutil
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv_path = os.path.join(root, "instructions", "done", "31_slots.csv")
    md_path = os.path.join(root, "instructions", "done", "31_slots.md")
    if not os.path.exists(csv_path):
        pytest.skip("the slot table has not been generated on this machine")
    before = {p: open(p, encoding="utf8").read() for p in (csv_path, md_path)}
    try:
        runpy.run_path(os.path.join(root, "scripts", "generate_slot_table.py"),
                       run_name="__main__")
        after = {p: open(p, encoding="utf8").read() for p in (csv_path, md_path)}
        assert after == before, "the slot table is stale -- run scripts/generate_slot_table.py"
    finally:
        for p, text in before.items():
            open(p, "w", encoding="utf8").write(text)


def test_every_slot_that_claims_coverage_has_columns_behind_it():
    """A slot graded A with nothing filling it would be a promise the cache does not keep."""
    import csv as _csv
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "instructions", "done", "31_slots.csv")
    if not os.path.exists(path):
        pytest.skip("the slot table has not been generated on this machine")
    rows = list(_csv.DictReader(open(path, encoding="utf8")))
    assert len(rows) > 50, "the taxonomy lost most of its slots"
    for r in rows:
        if r["grade"] != "-":
            assert r["filled_by"], f"{r['slot']} is graded {r['grade']} with nothing filling it"
            assert int(r["genes"]) > 0, r["slot"]
        else:
            assert not r["filled_by"], f"{r['slot']} is graded empty but names {r['filled_by']}"


def test_brain_fpkm_is_assigned_to_transcription_never_fitness():
    """The mixed transcriptome/proteome source was mislabeled as an in-vivo fitness screen.

    Its cache columns are FPKM expression measurements. Keeping the semantic assertion beside the
    generated table prevents a later regeneration from moving them back under fitness by prefix.
    """
    import csv as _csv
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rows = list(_csv.DictReader(open(os.path.join(root, "instructions", "done", "31_slots.csv"),
                                     encoding="utf8")))
    carrying = [r for r in rows
                if any(token.strip().startswith("invivo_")
                       for token in r["filled_by"].split(","))]
    assert carrying
    assert all(r["axis"] == "transcription" for r in carrying)
    assert not [r for r in rows if r["axis"] == "fitness"
                and any(token.strip().startswith("invivo_")
                        for token in r["filled_by"].split(","))]


def test_every_cited_study_carries_its_title():
    """A proposal that names a PMID and nothing else asks the reader to go and look up what is being
    proposed, which is most of the work of reading a list like this."""
    import csv as _csv
    import re
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "instructions", "done", "31_slots.csv")
    if not os.path.exists(path):
        pytest.skip("the slot table has not been generated on this machine")
    sys.path.insert(0, os.path.join(root, "scripts"))
    import generate_slot_table as G

    for r in _csv.DictReader(open(path, encoding="utf8")):
        for pmid in re.findall(r"PMID (\d{7,8})", r["candidates"] or ""):
            assert pmid in G.REFERENCES, f"{r['slot']} cites {pmid} with no reference"
            year, journal, title = G.REFERENCES[pmid]
            assert title and len(title) > 20, f"{pmid} has no usable title"
            assert pmid in r["candidate_titles"], f"{r['slot']} does not render {pmid}'s title"
    # And nothing mangled by the line wrapping that keeps the reference block readable.
    for _, _, title in G.REFERENCES.values():
        assert not [w for w in title.split() if len(w) > 34], f"a wrapped title lost a space: {title}"


def test_candidate_queries_are_derived_from_slot_and_organism():
    import propose_datasets as P
    row = next(r for r in __import__("generate_slot_table").all_slots("Tg")
               if r["axis"] == "translation")
    query = P.query_for(row)
    assert "Toxoplasma" in query and "ribosome profiling" in query
    assert all(word in query for word in ("tachyzoite", "vitro"))
