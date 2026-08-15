#!/usr/bin/env python3
"""A search survives being put down.

The property under test is a round trip: everything that was in the run before it was written is
readable afterwards, and the parts that cannot be trusted after a reload say so rather than
pretending. The second half matters more than the first -- a run that reloads cleanly against the
wrong node table is worse than one that fails to reload at all.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import searches as S  # noqa: E402


def nodes(n=50, seed=0):
    return pd.DataFrame({"gene_id": [f"TGME49_{i:06d}" for i in range(n)],
                         "product": ["hypothetical protein"] * n,
                         "n_publications": 0})


def finding(cluster=3, **kw):
    row = {"kind": "guilt", "layer": "compartment", "layer_kind": "discrete", "cluster": cluster,
           "category": "IMC", "n_cluster": 20, "n_known": 15, "n_hits": 12, "purity": 0.8,
           "background": 0.05, "lift": 16.0, "p": 1e-12, "q": 1e-10, "n_predicted": 5,
           "circular": False, "genes": [f"TGME49_{i:06d}" for i in range(5)]}
    row.update(kw)
    return row


def climbed(n_configs=3):
    """What `optimize.climb` returns: a table with artefact columns hanging off it."""
    return pd.DataFrame([
        {"restart": 0, "step": i, "score": 10.0 - i, "n_clusters": 20 + i, "mean_auprc": 0.3,
         "algorithm": "kmeans", "blocks": "expression_summary",
         "_labels": np.arange(50) % (3 + i),
         "_findings": pd.DataFrame([finding(cluster=i)])}
        for i in range(n_configs)])


# --------------------------------------------------------------------------- splitting a climb up
def test_a_climb_becomes_a_table_and_its_artefacts():
    run = S.from_climb(climbed(), nodes())
    assert list(run.configs.score) == [10.0, 9.0, 8.0]
    assert not any(c.startswith("_") for c in run.configs.columns), "an artefact reached the table"
    assert len(run.findings) == 3 and set(run.findings.config) == {0, 1, 2}
    assert set(run.labels) == {0, 1, 2}
    assert run.manifest["version"] and run.manifest["created"]
    assert set(("umap", "tsne", "hdbscan", "kmeans", "dbscan")) <= set(run.manifest["backend"])


def test_a_configuration_can_be_asked_for_by_rank():
    """0 is the winner and 1 is the runner-up, which is the comparison this exists to allow: a best
    configuration with nothing behind it cannot be told from a lucky one."""
    run = S.from_climb(climbed(), nodes())
    row, found, labels = run.for_config(1)
    assert row.score == 9.0
    assert list(found.cluster) == [1]
    assert len(labels) == 50


def test_asking_for_a_configuration_that_is_not_there_returns_nothing():
    run = S.from_climb(climbed(), nodes())
    row, found, labels = run.for_config(99)
    assert row is None and found.empty and labels is None
    row, found, labels = run.for_config(-1)
    assert row is None and found.empty and labels is None
    assert S.Search().for_config(0)[0] is None
    assert S.from_climb(None, nodes()).configs.empty


def test_a_climb_that_found_nothing_still_becomes_a_search():
    bare = pd.DataFrame([{"score": 0.0, "_labels": np.zeros(50), "_findings": pd.DataFrame()}])
    run = S.from_climb(bare, nodes())
    assert len(run.configs) == 1 and run.findings.empty
    assert run.for_config(0)[1].empty


# --------------------------------------------------------------------------- the round trip
def test_everything_survives_being_written_and_read_back(tmp_path):
    store = S.SearchStore(str(tmp_path))
    before = S.from_climb(climbed(), nodes(), mode="guilt", layer="compartment", name="mine")
    store.save(before)
    after = store.load("mine")
    pd.testing.assert_frame_equal(before.configs, after.configs, check_dtype=False)
    assert list(after.findings.cluster) == list(before.findings.cluster)
    assert list(after.findings.iloc[0].genes) == list(before.findings.iloc[0].genes)
    for k, v in before.labels.items():
        assert np.array_equal(after.labels[k], v)
    assert after.manifest["mode"] == "guilt" and after.manifest["layer"] == "compartment"


def test_a_saved_search_is_readable_without_this_program(tmp_path):
    """The manifest is JSON and the configurations are CSV on purpose. A run nobody can open in a
    text editor is a run that stops being evaluable the day the reader stops working."""
    store = S.SearchStore(str(tmp_path))
    store.save(S.from_climb(climbed(), nodes(), name="plain", mode="guilt"))
    where = store.path("plain")
    with open(os.path.join(where, "manifest.json")) as fh:
        assert json.load(fh)["mode"] == "guilt"
    assert "score" in open(os.path.join(where, "configs.csv")).readline()


def test_saved_searches_are_listed_newest_first(tmp_path):
    store = S.SearchStore(str(tmp_path))
    for name in ("search_20260101_000000", "search_20260814_120000"):
        store.save(S.from_climb(climbed(1), nodes(), name=name))
    assert [n for n, _m in store.list()] == ["search_20260814_120000", "search_20260101_000000"]


def test_a_directory_that_is_not_a_search_is_skipped_rather_than_fatal(tmp_path):
    store = S.SearchStore(str(tmp_path))
    store.save(S.from_climb(climbed(1), nodes(), name="good"))
    os.makedirs(os.path.join(str(tmp_path), "junk"), exist_ok=True)
    open(os.path.join(str(tmp_path), "junk", "manifest.json"), "w").write("{ not json")
    os.makedirs(os.path.join(str(tmp_path), "empty"), exist_ok=True)
    assert [n for n, _m in store.list()] == ["good"]


def test_a_search_with_pieces_missing_loads_as_far_as_it_can(tmp_path):
    """Half a run is still worth reading. Raising here would make one truncated write cost every
    other thing in the directory."""
    store = S.SearchStore(str(tmp_path))
    store.save(S.from_climb(climbed(), nodes(), name="partial"))
    os.remove(os.path.join(store.path("partial"), "labels.npz"))
    open(os.path.join(store.path("partial"), "findings.json"), "w").write("")
    back = store.load("partial")
    assert len(back.configs) == 3 and back.findings.empty and back.labels == {}
    assert store.load("never_saved").configs.empty


# --------------------------------------------------------------------------- the fingerprint
def test_a_search_knows_which_table_it_was_about():
    run = S.from_climb(climbed(), nodes())
    assert run.matches(nodes())
    assert not run.matches(nodes(n=40)), "a table of a different size passed as the same one"
    shuffled = nodes().iloc[::-1].reset_index(drop=True)
    assert not run.matches(shuffled), "a reordered table passed as the same one"


def test_the_reading_warns_when_the_table_has_moved_on():
    """A search is statements about genes identified by POSITION. Reloaded against a different
    table those positions address different genes, and every claim quietly becomes one about the
    wrong ones."""
    run = S.from_climb(climbed(), nodes())
    assert "different node table" in run.report(nodes(n=40))
    assert "different node table" not in run.report(nodes())
    assert "IMC" in run.report(nodes())


def test_a_search_from_nowhere_has_no_fingerprint_to_check():
    run = S.from_climb(climbed(), None)
    assert run.manifest["fingerprint"]["n_genes"] == 0
    assert not run.matches(nodes())


def test_a_search_describes_itself_in_one_line():
    run = S.from_climb(climbed(), nodes(), mode="disagreement", layer="compartment", name="x")
    line = run.describe()
    assert "disagreement" in line and "3 configurations" in line and "best 10.00" in line
    assert "0 configurations" in S.Search().describe()
