#!/usr/bin/env python3
"""Publishing a derived release, and fetching the identity tables.

The licensing logic here decides what gets redistributed to strangers, so the property that matters is
that it fails CLOSED: a study whose license cannot be positively confirmed permissive is withheld, and
"no local full text" is treated as unconfirmed rather than as permission. Getting that backwards would
mirror other people's copyrighted files under a CC-BY banner.

Nothing here touches the network or HuggingFace.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import fetch_names as FN  # noqa: E402
from starplast import hf_publish as HF  # noqa: E402


# --------------------------------------------------------------------------- licences
def _cat(pmids_pmcids):
    return pd.DataFrame(pmids_pmcids, columns=["pmid", "pmcid"])


def _jats(d, pmcid, body):
    (d / f"{pmcid}.xml").write_text(body, encoding="utf8")


def test_a_cc_by_href_is_redistributable(tmp_path):
    _jats(tmp_path, "PMC1", '<license xlink:href="https://creativecommons.org/licenses/by/4.0/">x</license>')
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert bool(out.loc[0, "redistributable"])
    assert out.loc[0, "source"] == "href"


def test_a_noncommercial_licence_is_withheld(tmp_path):
    _jats(tmp_path, "PMC1", '<license xlink:href="https://creativecommons.org/licenses/by-nc/4.0/">x</license>')
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert not bool(out.loc[0, "redistributable"])


def test_a_no_derivatives_licence_is_withheld(tmp_path):
    _jats(tmp_path, "PMC1", '<license xlink:href="https://creativecommons.org/licenses/by-nd/4.0/">x</license>')
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert not bool(out.loc[0, "redistributable"])


def test_public_domain_is_redistributable(tmp_path):
    _jats(tmp_path, "PMC1", '<license xlink:href="https://creativecommons.org/publicdomain/zero/1.0/">x</license>')
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert bool(out.loc[0, "redistributable"])


def test_a_licence_stated_only_in_prose_is_accepted_and_marked_as_prose(tmp_path):
    """PERMISSIVE lists `^by \\(text\\)$` explicitly, so a CC-BY stated in a paragraph rather than an
    href does grant redistribution. That is a deliberate choice, and the `source` column records which
    evidence it rested on so a prose-derived permission can be re-audited without re-parsing.

    Worth a human look before the release is made public: a prose match is weaker evidence than a
    machine-readable href, and this is the one place in the project where being wrong redistributes
    somebody else's copyrighted work."""
    _jats(tmp_path, "PMC1", "<license><p>This is a Creative Commons Attribution article.</p></license>")
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert out.loc[0, "source"] == "text"
    assert out.loc[0, "license"] == "by (text)"
    assert bool(out.loc[0, "redistributable"])


def test_noncommercial_in_prose_is_detected(tmp_path):
    _jats(tmp_path, "PMC1", "<license><p>Noncommercial use only.</p></license>")
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert out.loc[0, "license"] == "nc (text)"
    assert not bool(out.loc[0, "redistributable"])


def test_an_unrecognised_prose_licence_is_recorded_as_other(tmp_path):
    _jats(tmp_path, "PMC1", "<license><p>All rights reserved by the publisher.</p></license>")
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert out.loc[0, "license"] == "other (text)"
    assert not bool(out.loc[0, "redistributable"])


def test_no_local_full_text_means_unconfirmed_not_permitted(tmp_path):
    """Failing closed is the whole point: absence of evidence must not become permission."""
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC_missing"]]), str(tmp_path))
    assert out.loc[0, "source"] == "no local full text"
    assert not bool(out.loc[0, "redistributable"])


def test_a_study_with_no_pmcid_is_unconfirmed(tmp_path):
    out = HF.licenses_from_fulltexts(_cat([[1, None]]), str(tmp_path))
    assert not bool(out.loc[0, "redistributable"])


def test_a_full_text_with_no_licence_block_is_unconfirmed(tmp_path):
    _jats(tmp_path, "PMC1", "<article><body>no licence here</body></article>")
    out = HF.licenses_from_fulltexts(_cat([[1, "PMC1"]]), str(tmp_path))
    assert out.loc[0, "license"] == ""
    assert not bool(out.loc[0, "redistributable"])


# --------------------------------------------------------------------------- the plan
def _release_inputs():
    members = pd.DataFrame({"pmid": [1, 1, 2], "gene_id": ["g1", "g2", "g1"]})
    studies = pd.DataFrame({"pmid": [1, 2], "title": ["a", "b"]})
    licenses = pd.DataFrame({"pmid": [1, 2], "pmcid": ["PMC1", "PMC2"],
                             "license": ["by", "nc"], "source": ["href", "href"],
                             "redistributable": [True, False]})
    return members, studies, licenses


def test_the_plan_counts_what_would_be_withheld_and_why():
    """Called before building, so the decision is visible before anything is written."""
    p = HF.plan(*_release_inputs())
    assert p["studies_total"] == 2
    assert p["studies_redistributable"] == 1
    assert p["studies_withheld"] == 1
    assert "not positively confirmed" in p["withheld_reason"]
    assert p["derived_rows"] == 3 and p["derived_genes"] == 2


def test_the_plan_handles_an_empty_membership_table():
    members, studies, licenses = _release_inputs()
    p = HF.plan(members.iloc[0:0], studies, licenses)
    assert p["derived_rows"] == 0 and p["derived_genes"] == 0


def test_a_study_missing_from_the_licence_table_is_withheld():
    members, studies, licenses = _release_inputs()
    p = HF.plan(members, studies, licenses.iloc[0:0])
    assert p["studies_redistributable"] == 0


# --------------------------------------------------------------------------- the release
def test_a_release_writes_the_tables_and_a_card(tmp_path):
    members, studies, licenses = _release_inputs()
    out = HF.build_release(str(tmp_path), members, studies, licenses, log=lambda *_: None)
    assert os.path.exists(os.path.join(out, "study_gene_membership.parquet"))
    assert os.path.exists(os.path.join(out, "studies.parquet"))
    card = open(os.path.join(out, "README.md")).read()
    assert "membership, not interaction" in card


def test_the_card_states_the_membership_caveat_prominently(tmp_path):
    """Someone downloading this will otherwise treat every row as an interaction, which would
    manufacture tens of thousands of false edges."""
    members, studies, licenses = _release_inputs()
    HF.build_release(str(tmp_path), members, studies, licenses, log=lambda *_: None)
    card = open(os.path.join(tmp_path, "README.md")).read()
    assert "complete quantification table" in card
    # The figures are computed from the members passed in, not written into the card by hand. They
    # were constants, and the median had drifted to 754 against an actual 599 -- a stale number in
    # the one section whose whole job is to stop people reading membership as interaction.
    biggest = int(members.groupby("pmid").size().max())
    assert f"{biggest:,}" in card
    assert f"{int(members.groupby('pmid').size().median()):,}" in card


def test_the_card_reports_how_many_licences_were_confirmed(tmp_path):
    members, studies, licenses = _release_inputs()
    HF.build_release(str(tmp_path), members, studies, licenses, log=lambda *_: None)
    card = open(os.path.join(tmp_path, "README.md")).read()
    assert "1 carry a confirmed CC-BY" in card or "Of the source articles, 1" in card


def test_an_empty_licence_table_withholds_everything_instead_of_crashing(tmp_path):
    """It has no columns at all, so merging on "pmid" raised KeyError. Nothing confirmed means nothing
    redistributable -- which is the correct release, not a failure."""
    members, studies, _ = _release_inputs()
    HF.build_release(str(tmp_path), members, studies, pd.DataFrame(), log=lambda *_: None)
    assert os.path.exists(os.path.join(tmp_path, "README.md"))
    out = pd.read_parquet(os.path.join(tmp_path, "studies.parquet"))
    assert not out.redistributable.any()


def test_only_derived_facts_are_written_never_the_source_files(tmp_path):
    """The withheld studies' raw supplements must not be mirrored.

    Asserted as "every file is a derived table or documentation" rather than as a fixed list, so
    adding a document to the release does not read as a licensing regression -- while a mirrored
    .xlsx still does.
    """
    members, studies, licenses = _release_inputs()
    HF.build_release(str(tmp_path), members, studies, licenses, log=lambda *_: None)
    written = sorted(os.listdir(tmp_path))
    assert {"README.md", "studies.parquet", "study_gene_membership.parquet"} <= set(written)
    allowed = {".parquet", ".md", ".ipynb"}
    assert all(os.path.splitext(f)[1] in allowed for f in written), written


def test_the_methods_and_the_notebook_travel_with_the_data(tmp_path):
    """A dataset card is a summary. Anyone deciding whether they may use this needs the membership
    caveat in full and a way to regenerate the tables, and a link to a private repo is neither."""
    members, studies, licenses = _release_inputs()
    HF.build_release(str(tmp_path), members, studies, licenses, log=lambda *_: None)
    written = set(os.listdir(tmp_path))
    assert "METHODS.md" in written
    assert "build_hf_release.ipynb" in written
    methods = open(os.path.join(tmp_path, "METHODS.md")).read()
    # Whitespace-normalized: the source is hard-wrapped, so the phrase spans a newline.
    flat = " ".join(methods.split())
    assert "membership, not interaction" in flat.lower()
    assert "No raw supplementary file is mirrored" in flat


def test_a_release_staged_without_the_docs_says_so(tmp_path, monkeypatch):
    """An installed wheel has no docs/ beside the package, and a release quietly missing its methods
    is exactly the failure the copy exists to prevent."""
    members, studies, licenses = _release_inputs()
    monkeypatch.setattr(HF.os.path, "exists", lambda p: False)
    said = []
    HF.build_release(str(tmp_path), members, studies, licenses, log=said.append)
    assert any("METHODS.md is not in this release" in m for m in said)


# --------------------------------------------------------------------------- upload
def test_upload_without_the_library_reports_rather_than_raising(monkeypatch):
    import builtins
    real = builtins.__import__

    def no_hub(name, *a, **k):
        if name == "huggingface_hub":
            raise ImportError("not installed")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_hub)
    msgs = []
    assert HF.upload("me/ds", "/tmp", log=msgs.append) is None
    assert any("huggingface_hub is not installed" in m for m in msgs)


def test_upload_creates_a_private_repo_by_default(monkeypatch, tmp_path):
    """A dataset made public by accident cannot be made private again in any meaningful sense."""
    seen = {}

    class FakeApi:
        def create_repo(self, repo_id, repo_type=None, private=None, exist_ok=None):
            seen.update(repo_id=repo_id, private=private, repo_type=repo_type)

        def upload_folder(self, folder_path=None, repo_id=None, repo_type=None):
            seen["uploaded"] = folder_path

    import types
    mod = types.ModuleType("huggingface_hub")
    mod.HfApi = FakeApi
    monkeypatch.setitem(sys.modules, "huggingface_hub", mod)
    assert HF.upload("me/ds", str(tmp_path), log=lambda *_: None) == "me/ds"
    assert seen["private"] is True
    assert seen["repo_type"] == "dataset"
    assert seen["uploaded"] == str(tmp_path)


# --------------------------------------------------------------------------- identity tables
def test_toxodb_display_names_are_renamed_to_stable_columns(tmp_path):
    """ToxoDB ships human display names as the header; downstream code keys on stable ones."""
    p = tmp_path / "out.tsv"
    n = FN.write("Gene ID\tGene Name or Symbol\tPrevious ID(s)\n"
                 "TGME49_200010\tGRA16\tTGME49_008830\n", str(p))
    header = open(p).read().splitlines()[0]
    assert header.split("\t") == ["gene_id", "gene_name", "previous_ids"]
    assert n == 1


def test_an_unmapped_column_keeps_its_own_name(tmp_path):
    p = tmp_path / "out.tsv"
    FN.write("Gene ID\tSomething Else\nTGME49_200010\tx\n", str(p))
    assert open(p).read().splitlines()[0].split("\t") == ["gene_id", "Something Else"]


def test_an_empty_response_writes_an_empty_file_rather_than_crashing(tmp_path):
    p = tmp_path / "out.tsv"
    assert FN.write("", str(p)) == -1
    assert os.path.exists(p)


def test_the_request_body_asks_for_the_attributes_it_needs(monkeypatch):
    """gene_previous_ids is the one that matters: without it, papers citing pre-2012 accessions resolve
    to nothing."""
    seen = {}

    class Resp:
        def read(self):
            return b"Gene ID\n"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        seen["body"] = req.data.decode()
        return Resp()

    monkeypatch.setattr(FN.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setenv(FN.API_KEY_ENV, "a-key")
    FN.fetch("Toxoplasma gondii ME49", ["primary_key", "gene_previous_ids"])
    assert "gene_previous_ids" in seen["body"]
    assert "Toxoplasma gondii ME49" in seen["body"]


def test_main_writes_an_identity_table_for_each_arm(monkeypatch, tmp_path):
    """Four tables now: three ToxoDB and one PlasmoDB, each fetched from its own site.

    The Plasmodium one is separate rather than appended, for the reason the whole arm is separate:
    two identifier spaces in one index is the merge this project refuses everywhere else.
    """
    asked = []
    monkeypatch.setattr(FN, "OUT", str(tmp_path))
    monkeypatch.setattr(FN, "fetch", lambda org, attrs, url=FN.URL: (
        asked.append((org, url)) or "Gene ID\tGene Name or Symbol\nTGME49_1\tX\n"))
    FN.main()
    assert sorted(os.listdir(tmp_path)) == ["plasmodb_identity.tsv", "toxodb_identity.tsv",
                                            "toxodb_strain_gt1.tsv", "toxodb_strain_veg.tsv"]
    assert asked[-1] == ("Plasmodium falciparum 3D7", FN.PLASMODB_URL)
    assert all("toxodb" in url for _, url in asked[:-1])


# --------------------------------------------------------------------------- module entry point
def test_the_package_entry_point_calls_main(monkeypatch):
    """`python -m starplast` is a documented way to launch it."""
    called = {}
    import starplast.app as app
    monkeypatch.setattr(app, "main", lambda: called.setdefault("ran", True))
    import runpy
    runpy.run_module("starplast", run_name="__main__")
    assert called.get("ran")


def test_fetch_names_runs_as_a_module(monkeypatch, tmp_path):
    """`python -m starplast.fetch_names` is the documented one-off that produces the identity tables.

    Patched through the environment rather than the imported object: runpy executes a fresh copy of the
    module, so OUT is recomputed from paths.data_dir() and an attribute patch on the already-imported
    module would silently write to the real cache instead."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_CACHE, str(tmp_path))
    monkeypatch.setenv(FN.API_KEY_ENV, "a-key")

    class Resp:
        def read(self):
            return b"Gene ID\tGene Name or Symbol\nTGME49_1\tX\n"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    import urllib.request
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: Resp())
    import runpy
    runpy.run_module("starplast.fetch_names", run_name="__main__")
    assert os.path.exists(os.path.join(tmp_path, "toxodb_identity.tsv"))
    assert open(os.path.join(tmp_path, "toxodb_identity.tsv")).read().startswith("gene_id")
    # And the Plasmodium arm's, which the same one-off now produces: it had no identity layer at all
    # until a deposit keyed on pre-2012 accessions joined zero of its 3,629 genes.
    assert os.path.exists(os.path.join(tmp_path, FN.PLASMODB_IDENTITY))


# --------------------------------------------------------------------------- the GPU install
def _read(*parts):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return open(os.path.join(root, *parts), encoding="utf8").read()


def test_the_checkout_installs_the_gpu_stack_too():
    """`pip install -e .` must be the GPU one, not only `pip install starplast`.

    This reverses what this file used to assert. The code distribution deliberately carried no CUDA so
    that a `starplast-cpu` metapackage could leave it out -- but that made a checkout, which is where
    the person installing knows their own hardware, the one install that came out CPU-only. Chosen
    2026-08-17: the checkout gets the GPU stack, and `starplast-cpu` is gone rather than reduced to a
    name that installs the opposite of what it says (it depended on this distribution, so it would
    have inherited the wheels transitively with no way to decline them).
    """
    main = _read("pyproject.toml")
    assert 'name = "starplast-core"' in main
    defaults = main.split("[project.optional-dependencies]")[0]
    assert "cuml-cu12" in defaults and "cupy-cuda12x" in defaults, \
        "the code distribution must carry CUDA: `pip install -e .` is a GPU install"


def test_the_cuda_wheels_in_the_checkout_are_marked_for_the_platform_that_has_them():
    """Unmarked, `pip install -e .` would FAIL on macOS, Windows and aarch64 -- where RAPIDS and CuPy
    publish nothing -- rather than installing a working CPU program."""
    # Every requirement naming a CUDA wheel, in the defaults AND in the extra -- comments mention the
    # package names too, and a comment is not a requirement.
    for line in _read("pyproject.toml").splitlines():
        stripped = line.strip()
        if not stripped.startswith('"'):
            continue
        if "cuml-cu12" in stripped or "cupy-cuda12x" in stripped:
            assert "sys_platform == 'linux'" in stripped, stripped
            assert "platform_machine == 'x86_64'" in stripped, stripped


def test_the_cpu_name_is_gone_rather_than_lying():
    """Kept as a test because the temptation is to bring the name back. It cannot work: anything that
    depends on `starplast-core` inherits its CUDA wheels, so a `-cpu` metapackage would install the
    two gigabytes it exists to avoid."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    assert not os.path.exists(os.path.join(root, "packaging", "starplast-cpu")), \
        "starplast-cpu is back, and it can only install the CUDA stack it promises to omit"


def test_the_default_install_is_the_gpu_one():
    """`pip install starplast` is the fast one, which is the point of the split."""
    from starplast import __version__
    meta = _read("packaging", "starplast", "pyproject.toml")
    assert 'name = "starplast"' in meta
    assert f'"starplast-core=={__version__}"' in meta
    assert "cuml-cu12" in meta and "cupy-cuda12x" in meta


def test_the_cuda_wheels_are_marked_for_the_platform_that_has_them():
    """RAPIDS and CuPy publish Linux wheels only. Unmarked, `pip install starplast` would FAIL on
    macOS and Windows rather than installing a working program."""
    meta = _read("packaging", "starplast", "pyproject.toml")
    for line in meta.splitlines():
        if "cuml-cu12" in line or "cupy-cuda12x" in line:
            assert "sys_platform == 'linux'" in line and "platform_machine == 'x86_64'" in line, line


def test_the_gpu_name_still_works_as_an_alias():
    """Someone who read the old instructions gets what they expected."""
    from starplast import __version__
    meta = _read("packaging", "starplast-gpu", "pyproject.toml")
    assert f'dependencies = ["starplast=={__version__}"]' in meta


def test_every_distribution_ships_the_same_version():
    """Three files, one release. A metapackage pinned to a version that does not exist is an install
    that fails for a reason nobody can see from the error."""
    from starplast import __version__
    for parts in (("pyproject.toml",), ("packaging", "starplast", "pyproject.toml"),
                  ("packaging", "starplast-gpu", "pyproject.toml")):
        assert f'version = "{__version__}"' in _read(*parts), parts


def test_the_program_says_how_to_install_the_gpu_stack():
    """The switch reports what it found; when it found nothing it has to say what to type."""
    from starplast import gpu
    text = gpu.describe()
    if not any(gpu.available()[k] for k in ("cuml", "cupy", "torch")):
        assert "starplast-gpu" in text or "starplast[gpu]" in text


def test_fetching_without_a_key_says_how_to_get_one(monkeypatch):
    """VEuPathDB closed these reports to anonymous use in August 2026, mid-session.

    The failure to avoid is the one that arrives as `HTTPError: 401` from inside urllib, three
    frames below anything this project wrote, on a machine where the committed tables mean nothing
    is actually broken. The message names the environment variable and says the tables are shipped.
    """
    monkeypatch.delenv(FN.API_KEY_ENV, raising=False)
    with pytest.raises(PermissionError) as caught:
        FN.fetch("Toxoplasma gondii ME49", ["primary_key"])
    assert FN.API_KEY_ENV in str(caught.value)
    assert "COMMITTED" in str(caught.value) or "committed" in str(caught.value)


def test_every_unguarded_import_in_the_package_is_a_declared_dependency():
    """The rule behind three separate incidents, stated once.

    `xlrd` was offered by the importer and never declared; `pypdf` arrived with the interactome;
    `networkx` had been imported by the multiplex community detection all along and reached the
    program only because something else pulled it in -- on the interpreter the user actually runs it
    was simply absent, and four tests failed for a reason that had nothing to do with the code.

    An import inside a `try` is a different statement: it says the feature is optional and the code
    handles its absence. Those are exempt, which is why the GPU backends and the HDBSCAN fallback do
    not have to be declared.
    """
    import ast
    import os
    import sys
    import tomllib

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    declared = tomllib.load(open(os.path.join(root, "pyproject.toml"), "rb"))["project"]["dependencies"]
    names = {d.split(">=")[0].split("==")[0].split("[")[0].strip().lower() for d in declared}
    # Import name -> distribution name, where they differ.
    distribution = {"sklearn": "scikit-learn", "umap": "umap-learn", "opengl": "pyopengl",
                    "pyqt6": "pyqt6", "pil": "pillow", "yaml": "pyyaml"}
    # Deliberately optional, each behind a check the code makes before importing -- an early return
    # on `available()`, or a `try`. Listed by name rather than inferred from syntax, because a guard
    # can be written four ways and a list of four packages is reviewable.
    optional = {"torch": "GPU distances; gpu.available() is asked first",
                "cuml": "GPU UMAP and clustering; same check",
                "cupy": "GPU array work; same check",
                "hdbscan": "the fallback when scikit-learn is too old",
                "pybigwig": "reading coverage tracks, which only the build does",
                "huggingface_hub": "publishing a release, which the application never does"}
    names |= set(optional)
    undeclared = {}
    for folder, _dirs, files in os.walk(os.path.join(root, "starplast")):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(folder, name)
            tree = ast.parse(open(path, encoding="utf8").read(), path)
            guarded = {node for outer in ast.walk(tree)
                       if isinstance(outer, ast.Try) for node in ast.walk(outer)}
            for node in ast.walk(tree):
                if node in guarded:
                    continue
                modules = ([alias.name for alias in node.names] if isinstance(node, ast.Import)
                           else [node.module] if isinstance(node, ast.ImportFrom)
                           and node.level == 0 and node.module else [])
                for module in modules:
                    top = module.split(".")[0]
                    if not top or top == "starplast" or top in sys.stdlib_module_names:
                        continue
                    key = distribution.get(top.lower(), top.lower())
                    if key not in names:
                        undeclared.setdefault(key, os.path.relpath(path, root))
    assert not undeclared, f"imported without a guard and not declared: {undeclared}"
