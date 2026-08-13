#!/usr/bin/env python3
"""Kept clusterings: named, dated, saved with their recipe, and knowing which genes they cover.

A clustering used to be computed, drawn and lost, so the second run replaced the first with no way
back -- which makes the one comparison this application exists for impossible to make by looking.
The tests here are mostly about the third property: a run over a subsample has fewer labels than the
map has genes, and lining them up by position would put cluster 3's color on whichever gene happens
to sit at that index.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import runs as R  # noqa: E402


def _run(store, name="r1", n=10, covered=None):
    genes = np.ones(n, bool) if covered is None else covered
    labels = np.array([0, 1, -1] * (int(genes.sum()) // 3 + 1))[:int(genes.sum())]
    return store.add(labels, genes, recipe={"algorithm": "hdbscan"}, name=name)


def test_a_run_knows_which_genes_it_applies_to_not_just_how_many():
    """The subsample case. Placed by position, cluster 3's color lands on whichever gene sits at
    that index -- confidently, and wrongly, for every gene after the first gap."""
    covered = np.zeros(12, bool)
    covered[[1, 4, 7]] = True
    run = R.Run(name="r", labels=np.array([0, 1, 0]), genes=covered)
    v = run.values(12)
    assert list(v[covered]) == ["cluster 0", "cluster 1", "cluster 0"]
    assert (v[~covered] == "").all(), "genes outside the run were given a cluster"


def test_a_gene_outside_a_run_is_absent_not_unclustered():
    """They were not put in no cluster; they were not in the map at all. Conflating the two puts
    thousands of genes into a category that means something else."""
    covered = np.array([True, False])
    v = R.Run(name="r", labels=np.array([-1]), genes=covered).values(2)
    assert v[0] == R.NOISE_NAME and v[1] == ""


def test_a_run_whose_mask_and_labels_disagree_colors_nothing():
    """Better an empty column than 8,140 confident mislabels."""
    v = R.Run(name="r", labels=np.array([0, 1, 2]), genes=np.array([True, False])).values(2)
    assert (v == "").all()


def test_a_run_summarises_itself_for_a_list():
    run = R.Run(name="hdbscan_60", labels=np.array([0, 0, 1, -1]), genes=np.ones(4, bool))
    assert run.n_clusters == 2 and run.noise_frac == pytest.approx(0.25) and run.n_genes == 4
    assert "hdbscan_60" in run.describe() and "2 clusters" in run.describe()


def test_an_empty_run_has_no_clusters_rather_than_dividing_by_zero():
    run = R.Run(name="none", labels=np.array([], dtype=int), genes=np.zeros(3, bool))
    assert run.n_clusters == 0 and run.noise_frac == 0.0


# --------------------------------------------------------------------------- naming
def test_runs_are_named_by_the_clock_to_the_second():
    """Two runs a minute apart have to be distinguishable without anyone typing anything."""
    import datetime
    a = R.timestamp_name(datetime.datetime(2026, 8, 12, 15, 4, 5))
    b = R.timestamp_name(datetime.datetime(2026, 8, 12, 15, 4, 6))
    assert a != b and a.endswith("20260812_150405")
    assert R.timestamp_name().startswith("run_")


def test_a_run_can_be_renamed_and_the_old_file_goes_with_it(tmp_path):
    store = R.RunStore(str(tmp_path))
    _run(store, "run_20260812_150405")
    store.save(store.runs[0])
    assert store.rename("run_20260812_150405", "gras_are_clean")
    assert store.names() == ["gras_are_clean"]
    assert os.path.exists(str(tmp_path / "gras_are_clean.json"))
    assert not os.path.exists(str(tmp_path / "run_20260812_150405.json"))


def test_a_name_already_in_use_is_refused_rather_than_suffixed(tmp_path):
    """Two runs called the same thing in a list you choose from is worse than being told no."""
    store = R.RunStore(str(tmp_path))
    _run(store, "a")
    _run(store, "b")
    assert store.rename("b", "a") is False and sorted(store.names()) == ["a", "b"]
    assert store.rename("nonesuch", "c") is False
    assert store.rename("a", "") is False


# --------------------------------------------------------------------------- on disk
def test_a_run_survives_the_window_with_its_recipe(tmp_path):
    """Without the recipe a name is a label on nothing: the run cannot be rebuilt, and two runs
    cannot be told apart except by their numbers."""
    store = R.RunStore(str(tmp_path))
    run = store.add(np.array([0, 1, -1]), np.ones(3, bool),
                    recipe={"algorithm": "hdbscan", "min_cluster_size": 60}, name="kept")
    store.save(run)
    fresh = R.RunStore(str(tmp_path))
    fresh.load_all()
    back = fresh.get("kept")
    assert back is not None and back.recipe["min_cluster_size"] == 60
    assert list(back.labels) == [0, 1, -1] and back.n_genes == 3


def test_loading_twice_does_not_duplicate(tmp_path):
    store = R.RunStore(str(tmp_path))
    store.save(_run(store, "once"))
    store.load_all()
    store.load_all()
    assert store.names() == ["once"]


def test_one_unreadable_run_does_not_cost_the_others(tmp_path, capsys):
    store = R.RunStore(str(tmp_path))
    store.save(_run(store, "good"))
    (tmp_path / "broken.json").write_text("{not json")
    fresh = R.RunStore(str(tmp_path))
    fresh.load_all()
    assert fresh.names() == ["good"]
    assert "could not read run" in capsys.readouterr().out


def test_a_store_with_nowhere_to_write_still_keeps_runs_in_memory(tmp_path):
    """The panel lists them constantly; the disk copy is so they survive the window."""
    store = R.RunStore(None)
    store.add(np.zeros(3, int), np.ones(3, bool), name="memory only")
    assert store.names() == ["memory only"] and store.save(store.runs[0]) == ""
    assert store.load_all() == store.runs


def test_a_name_with_a_slash_in_it_does_not_write_outside_the_store(tmp_path):
    store = R.RunStore(str(tmp_path))
    store.save(store.add(np.zeros(2, int), np.ones(2, bool), name="../escape"))
    assert sorted(os.listdir(tmp_path)) == [".._escape.json", ".._escape.npz"]


# --------------------------------------------------------------------------- binning
def test_a_quantity_becomes_a_category_with_its_ranges_as_labels():
    v = R.bin_column(pd.Series(np.arange(100.0)), bins=4)
    assert len({x for x in v if x}) == 4
    assert all("," in x for x in v if x), "the labels have to say what range they are"


def test_missing_values_are_absent_rather_than_a_low_bin():
    v = R.bin_column(pd.Series([1.0, 2.0, np.nan, 4.0]), bins=2)
    assert v[2] == "" and v[0] != ""


def test_a_column_that_is_mostly_one_value_gets_fewer_bins_and_says_why():
    """`n_publications` is zero for most of this proteome, so its quartile edges are all zero.
    Splitting the tie by rank or by equal width would draw four colors over a column with one
    level, which is a picture of a distinction that does not exist."""
    said = []
    v = R.bin_column(pd.Series([0.0] * 90 + list(range(1, 11))), bins=4, log=said.append)
    assert len({x for x in v if x}) < 4
    assert any("share one value" in m for m in said)


def test_a_constant_column_is_one_bin_and_says_so():
    said = []
    v = R.bin_column(pd.Series([3.0] * 20), bins=4, log=said.append)
    assert set(v) == {"all one value"} and any("one bin" in m for m in said)


def test_too_little_data_to_bin_leaves_everything_absent():
    assert (R.bin_column(pd.Series([np.nan, np.nan]), bins=4) == "").all()
    assert (R.bin_column(pd.Series([1.0, 2.0]), bins=1) == "").all()


def test_text_that_is_not_a_quantity_bins_to_nothing():
    assert (R.bin_column(pd.Series(["a", "b", "c"]), bins=2) == "").all()


def test_a_value_that_lands_in_no_bin_is_absence_not_a_bin_called_nan():
    """qcut does not raise on a constant column -- it returns NaN categories, and casting those to
    str gives the literal "nan", which would arrive in the interface as a category sitting in the
    legend beside real ones."""
    v = R.bin_column(pd.Series([1.0, 1.0, 1.0, 5.0, 9.0]), bins=3)
    assert "nan" not in set(v)
    assert all(x == "" or "," in x or x == "all one value" for x in v)


def test_qcut_raising_outright_is_also_one_bin(monkeypatch):
    """The other half of the same case: older pandas raises where newer returns NaNs, and the
    interface must not depend on which."""
    def boom(*a, **k):
        raise ValueError("Bin edges must be unique")

    monkeypatch.setattr(pd, "qcut", boom)
    said = []
    v = R.bin_column(pd.Series([3.0] * 10), bins=4, log=said.append)
    assert set(v) == {"all one value"} and any("one bin" in m for m in said)
