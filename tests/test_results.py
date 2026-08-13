#!/usr/bin/env python3
"""Saving what a tab computed, and loading it back still clickable.

The point of the format is the last word: a row in these tables is a recipe -- blocks, policy,
scaling, seed, subsample, excluded columns -- which is what lets clicking it rebuild a map. A loader
that produced a table you could only read would be a screenshot with extra steps, so the tests here
are about the round trip keeping the columns that make a row rebuildable, and about a file saying
which table it came from so it cannot be loaded into the wrong one.
"""
from __future__ import annotations

import os
import sys
import zipfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import results as R  # noqa: E402


def _search_rows():
    return pd.DataFrame({"blocks": ["fitness_screens", "expression_summary"],
                         "na_policy": ["median"] * 2, "scaling": ["rank"] * 2,
                         "n_neighbors": [15, 50], "min_dist": [0.1, 0.0],
                         "min_cluster_size": [25, 60], "seed": [42, 42],
                         "sample_size": [3000, 3000], "excluded": ["compartment"] * 2,
                         "mean_f1": [0.31, 0.28], "best_f1": [0.59, 0.44]})


def test_a_saved_table_says_which_table_it_is(tmp_path):
    """Without it, a validation file dropped into the walk tab makes rows nobody can rebuild and an
    error message about the wrong thing."""
    path = R.save_table(str(tmp_path / "s.csv"), "recovery_search", _search_rows())
    assert R.table_kind(path) == "recovery_search"
    assert open(path).readline().startswith(R.MARKER)


def test_a_saved_table_round_trips_with_the_columns_that_rebuild_a_map(tmp_path):
    df = _search_rows()
    got = R.load_table(R.save_table(str(tmp_path / "s.csv"), "recovery_search", df))
    assert list(got.columns) == list(df.columns)
    assert got.equals(df)
    for needed in ("blocks", "seed", "sample_size", "excluded", "min_cluster_size"):
        assert needed in got.columns, needed


def test_an_ordinary_csv_still_loads(tmp_path):
    """Somebody's own table, or one saved before this existed. It loads; the caller decides whether
    a file that does not name its table belongs where it is being put."""
    path = tmp_path / "plain.csv"
    _search_rows().to_csv(path, index=False)
    assert R.table_kind(str(path)) == ""
    assert len(R.load_table(str(path))) == 2


def test_a_file_that_is_not_there_names_no_table(tmp_path):
    assert R.table_kind(str(tmp_path / "nope.csv")) == ""


def test_a_bundle_carries_every_tab_and_says_what_is_in_it(tmp_path):
    tables = {"recovery_search": _search_rows(),
              "validation": pd.DataFrame({"category": ["dense granules"], "precision": [0.62]})}
    path = R.save_bundle(str(tmp_path / "all.starplast"), tables, meta={"target": "compartment"})
    assert zipfile.is_zipfile(path) and R.is_bundle(path)
    back, meta = R.load_bundle(path)
    assert set(back) == set(tables)
    assert back["recovery_search"].equals(tables["recovery_search"])
    assert meta["target"] == "compartment" and meta["tables"]["recovery_search"] == 2
    assert meta["starplast"], "a bundle that will not load in two years has to say what wrote it"


def test_empty_tables_are_left_out_rather_than_written_as_headers(tmp_path):
    """A bundle listing eight tables of which six are empty invites the reader to think six analyses
    returned nothing, when they were never run."""
    path = R.save_bundle(str(tmp_path / "some.starplast"),
                         {"recovery_search": _search_rows(),
                          "validation": pd.DataFrame(),
                          "umap_walk": None})
    back, meta = R.load_bundle(path)
    assert set(back) == {"recovery_search"} and set(meta["tables"]) == {"recovery_search"}


def test_one_unreadable_member_costs_that_table_and_not_the_bundle(tmp_path, capsys):
    """A bundle is a session's work; losing seven tables to one bad row would be its own failure."""
    path = str(tmp_path / "b.starplast")
    R.save_bundle(path, {"recovery_search": _search_rows()})
    with zipfile.ZipFile(path, "a") as z:
        z.writestr("broken.csv", b"\xff\xfe not text at all")
    back, _ = R.load_bundle(path)
    assert "recovery_search" in back
    assert "broken" not in back


def test_a_bundle_reads_a_member_that_never_named_its_table(tmp_path):
    path = str(tmp_path / "b.starplast")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("umap_walk.csv", "a,b\n1,2\n")
    back, meta = R.load_bundle(path)
    assert list(back) == ["umap_walk"] and meta == {}


def test_a_single_table_is_not_mistaken_for_a_bundle(tmp_path):
    path = R.save_table(str(tmp_path / "s.csv"), "umap_walk", _search_rows())
    assert not R.is_bundle(path)


def test_a_bundle_with_a_damaged_manifest_still_gives_up_its_tables(tmp_path):
    """The manifest is a convenience; the tables are the work."""
    path = str(tmp_path / "b.starplast")
    R.save_bundle(path, {"umap_walk": _search_rows()})
    keep = {n: zipfile.ZipFile(path).read(n) for n in zipfile.ZipFile(path).namelist()}
    with zipfile.ZipFile(path, "w") as z:
        for n, data in keep.items():
            z.writestr(n, b"{not json" if n == R.MANIFEST else data)
    back, meta = R.load_bundle(path)
    assert list(back) == ["umap_walk"] and meta == {}


def test_anything_in_a_bundle_that_is_not_a_table_is_ignored(tmp_path):
    """People add things to zips -- a note, a figure -- and a reader that chokes on them is a reader
    that punishes the habit."""
    path = str(tmp_path / "b.starplast")
    R.save_bundle(path, {"umap_walk": _search_rows()})
    with zipfile.ZipFile(path, "a") as z:
        z.writestr("notes.txt", "the second run looked better")
    back, _ = R.load_bundle(path)
    assert list(back) == ["umap_walk"]
