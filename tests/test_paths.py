#!/usr/bin/env python3
"""Data resolution: the layer that decides where every file comes from.

This is tested harder than its size suggests because the bug it fixes was invisible and total. The
dataset root used to be computed as "the repository's parent directory", which is true on one machine;
from a clean checkout 1 of 16 registry paths resolved, and nothing reported it. The application looked
fine, because the application uses the committed cache; only a rebuild failed, and only on someone
else's computer.

So the properties asserted here are the ones that were silently false before: that resolution never
depends on the current working directory, that an installed layout wins over a source layout, that a
missing raw dataset is a reported state rather than an exception, and that the two roots -- the cache
that must always work and the raw tree that is allowed to be absent -- never get confused.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import paths  # noqa: E402


# --------------------------------------------------------------------------- the cache
def test_the_cache_resolves_and_is_complete():
    ok, msg = paths.check()
    assert ok, msg
    for f in paths.REQUIRED:
        assert os.path.exists(os.path.join(paths.data_dir(), f))


def test_the_cache_lives_inside_the_package_so_a_wheel_carries_it():
    """It used to sit beside the package, so an installed wheel put a generic top-level `data/` into
    site-packages -- collidable with any other package -- and omitted the three identity tables
    entirely, because the glob matched only two extensions and did not recurse."""
    assert os.path.dirname(paths.data_dir()) == os.path.dirname(os.path.abspath(paths.__file__))


def test_resolution_does_not_depend_on_the_working_directory(tmp_path, monkeypatch):
    """The original bug in one line: paths were relative, so the answer changed with the shell's cwd."""
    before = paths.data_dir()
    monkeypatch.chdir(tmp_path)
    assert paths.data_dir() == before


def test_an_explicit_cache_override_wins(tmp_path, monkeypatch):
    (tmp_path / "nodes.parquet").write_bytes(b"")
    monkeypatch.setenv(paths.ENV_CACHE, str(tmp_path))
    assert paths.data_dir() == str(tmp_path)


def test_an_override_pointing_at_nothing_still_names_a_directory_rather_than_raising(monkeypatch):
    """check() is the layer that reports a missing cache; data_dir() must survive to be reported about."""
    monkeypatch.setenv(paths.ENV_CACHE, "/nonexistent/starplast-cache")
    d = paths.data_dir()
    assert isinstance(d, str) and d
    ok, msg = paths.check()
    assert not ok
    assert "build_graph" in msg, "the message has to say what to do, not just what is wrong"


def test_cache_file_joins_onto_the_resolved_cache():
    assert paths.cache_file("nodes.parquet") == os.path.join(paths.data_dir(), "nodes.parquet")


# --------------------------------------------------------------------------- the raw datasets
def test_a_missing_dataset_root_is_a_state_not_an_exception(monkeypatch, tmp_path):
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path / "absent"))
    assert paths.find("anything.csv") is None
    assert isinstance(paths.dataset_root(), str)


def test_dataset_root_can_create_a_writable_download_location(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path / "absent"))
    r = paths.dataset_root(create=True)
    assert os.path.isdir(r) and os.access(r, os.W_OK)


def test_find_searches_every_root_not_only_the_first(monkeypatch, tmp_path):
    """Stopping at the first existing root silently lost the CRISPR screens once already, when half the
    files had moved to datasets/<level>/<type>/<PMID>/ and half had not."""
    first, second = tmp_path / "a", tmp_path / "b"
    first.mkdir()
    second.mkdir()
    (second / "only_here.csv").write_text("x")
    monkeypatch.setattr(paths, "_dataset_candidates", lambda: iter([str(first), str(second)]))
    assert paths.find("only_here.csv") == str(second / "only_here.csv")


def test_roots_are_deduplicated_and_ordered(monkeypatch, tmp_path):
    monkeypatch.setattr(paths, "_dataset_candidates",
                        lambda: iter([str(tmp_path), str(tmp_path), None, "/nonexistent"]))
    assert paths.dataset_roots() == [str(tmp_path)]


def test_find_accepts_path_segments():
    assert paths.find("a", "b", "c.csv") is None      # joins without raising on multiple segments


# --------------------------------------------------------------------------- platform behaviour
def test_user_cache_dir_is_per_user_and_named(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    assert paths.user_cache_dir() == os.path.join(str(tmp_path), "starplast")


def test_user_cache_dir_falls_back_when_xdg_is_unset(monkeypatch):
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
    assert paths.user_cache_dir().startswith(os.path.expanduser("~"))


# --------------------------------------------------------------------------- reporting
def test_describe_names_the_resolved_locations_and_the_overrides():
    """A resolver that picks silently is harder to debug than one that fails, so this is the output a
    user is asked for first."""
    d = paths.describe()
    assert paths.data_dir() in d
    assert paths.ENV_DATASETS in d and paths.ENV_CACHE in d


def test_describe_reports_where_downloads_would_go_when_no_root_exists(monkeypatch, tmp_path):
    monkeypatch.setattr(paths, "dataset_roots", lambda: [])
    assert "downloads would go to" in paths.describe()


def test_describe_marks_an_incomplete_cache(monkeypatch):
    monkeypatch.setattr(paths, "check", lambda: (False, "nope"))
    assert "INCOMPLETE" in paths.describe()


def test_the_cli_runs_and_reports_the_registry(capsys):
    paths._cli()
    out = capsys.readouterr().out
    assert "registered datasets" in out
    assert "cache" in out


def test_the_cli_survives_a_broken_registry(monkeypatch, capsys):
    """Diagnostics must not require the thing being diagnosed to work."""
    import starplast.datasets as D
    monkeypatch.setattr(D, "missing", lambda: (_ for _ in ()).throw(RuntimeError("boom")))
    paths._cli()
    assert "registry unavailable" in capsys.readouterr().out


# --------------------------------------------------------------------------- the other platforms
def test_windows_cache_location(monkeypatch, tmp_path):
    """Packaged as an .exe, so the Windows branch is shipped code, not a hypothetical."""
    monkeypatch.setattr(os, "name", "nt")
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    assert paths.user_cache_dir() == os.path.join(str(tmp_path), "starplast")


def test_windows_cache_falls_back_to_home_when_localappdata_is_unset(monkeypatch):
    monkeypatch.setattr(os, "name", "nt")
    monkeypatch.delenv("LOCALAPPDATA", raising=False)
    assert paths.user_cache_dir().startswith(os.path.expanduser("~"))


def test_macos_cache_location(monkeypatch):
    """Packaged as a .dmg, likewise."""
    monkeypatch.setattr(os, "name", "posix")
    monkeypatch.setattr(os, "uname", lambda: type("u", (), {"sysname": "Darwin"})())
    assert paths.user_cache_dir() == os.path.expanduser("~/Library/Caches/starplast")


# --------------------------------------------------------------------------- the fallback ladder
def test_the_repository_layout_is_still_found(monkeypatch, tmp_path):
    """Caches built before the move sit beside the package rather than inside it. Existing checkouts
    must keep working across the upgrade."""
    pkg, repo = tmp_path / "pkg" / "starplast", tmp_path / "pkg"
    (repo / "data").mkdir(parents=True)
    (repo / "data" / "nodes.parquet").write_bytes(b"")
    pkg.mkdir(exist_ok=True)
    monkeypatch.delenv(paths.ENV_CACHE, raising=False)
    monkeypatch.setattr(paths, "_PKG", str(pkg))
    monkeypatch.setattr(paths, "_REPO", str(repo))
    assert paths.data_dir() == str(repo / "data")


def test_with_nothing_built_it_names_where_the_cache_should_go(monkeypatch, tmp_path):
    """The pre-build state. It must still return a path, because check() needs one to report about."""
    monkeypatch.delenv(paths.ENV_CACHE, raising=False)
    monkeypatch.setattr(paths, "_PKG", str(tmp_path / "nowhere" / "starplast"))
    monkeypatch.setattr(paths, "_REPO", str(tmp_path / "nowhere"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    d = paths.data_dir()
    assert d == os.path.join(str(tmp_path / "nowhere" / "starplast"), "data")
    assert not paths.check()[0]


def test_data_dir_last_resort_is_the_user_cache(monkeypatch, tmp_path):
    """Every candidate empty -- only reachable if the package has no location, but it is the branch that
    keeps data_dir() total, and a resolver that can return None puts a None into an os.path.join."""
    monkeypatch.delenv(paths.ENV_CACHE, raising=False)
    monkeypatch.setattr(paths, "_cache_candidates", lambda: iter([None, None]))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    assert paths.data_dir() == os.path.join(str(tmp_path), "starplast", "data")


def test_a_dataset_override_pointing_at_nothing_is_reported_loudly(monkeypatch, tmp_path):
    """Set-but-wrong looks identical to unset from the outside, and the user believes it took effect."""
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path / "typo"))
    assert "WARNING" in paths.describe()


def test_running_as_a_module_prints_the_report():
    """`python -m starplast.paths` is the documented first step when data cannot be found."""
    import runpy
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        runpy.run_module("starplast.paths", run_name="__main__")
    assert "cache" in buf.getvalue()


def test_an_existing_root_is_used_in_preference_to_the_download_location(monkeypatch, tmp_path):
    """The normal case on a machine that already has the raw data: use it, do not re-fetch it."""
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    assert paths.dataset_root() == str(tmp_path)
    assert paths.dataset_root(create=True) == str(tmp_path)


def test_describe_lists_each_dataset_root_it_found(monkeypatch, tmp_path):
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    out = paths.describe()
    assert str(tmp_path) in out
    assert "none found" not in out
