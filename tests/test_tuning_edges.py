#!/usr/bin/env python3
"""The UMAP walk, the embedding store, and CSV import.

Import is the path a user's own data takes into the map, so its identifier handling is where a
newcomer's data silently disappears. A column formatted `TgME49.208830` or `gene|TGME49_208830|v2`
matches nothing and produces an empty join, which looks exactly like a file with no genes in it -- so
the repair is suggested rather than demanded, and the number that resolved is always reported.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import tuning as TU  # noqa: E402
from starplast.embedding import EmbeddingSpec  # noqa: E402


def _nodes(n=150, seed=0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({f"fit_{i}": rng.normal(size=n) for i in range(4)})
    d.insert(0, "gene_id", [f"TGME49_{200000+i}" for i in range(n)])
    return d


# --------------------------------------------------------------------------- quality measures
def test_quality_measures_survive_a_degenerate_embedding():
    """Both measures raise on degenerate input, and reporting a fabricated number would rank that
    setting against real ones."""
    X = np.zeros((40, 3))
    out = TU._quality(X, X, n_neighbors=5)
    assert set(out) == {"trustworthiness", "continuity_proxy"}
    assert all(isinstance(v, float) for v in out.values())


def test_quality_is_reported_for_a_real_embedding():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(120, 5))
    Y = X[:, :2]
    out = TU._quality(X, Y, n_neighbors=10)
    assert 0.0 <= out["trustworthiness"] <= 1.0


# --------------------------------------------------------------------------- the walk
def test_the_walk_is_seeded_and_records_the_seed():
    out = TU.walk_umap(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)),
                       n_neighbors_values=(15,), min_dist_values=(0.1,),
                       sample_size=120, log=lambda *_: None)
    assert len(out) == 1
    assert out.iloc[0].seed == TU.DEFAULT_SEED
    assert out.iloc[0].sample_size == 120


def test_a_neighbourhood_larger_than_the_sample_is_skipped():
    out = TU.walk_umap(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)),
                       n_neighbors_values=(5000,), min_dist_values=(0.1,),
                       sample_size=120, log=lambda *_: None)
    assert out.empty


def test_the_walk_can_report_how_many_clusters_each_setting_produces():
    """The number a hyperparameter walk is usually being run to find out."""
    out = TU.walk_umap(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)),
                       n_neighbors_values=(15,), min_dist_values=(0.0,),
                       cluster_check=True, log=lambda *_: None)
    assert "n_clusters_hdbscan" in out.columns
    assert "noise_frac" in out.columns


def test_without_umap_the_walk_reports_and_returns_empty(monkeypatch):
    import builtins
    real = builtins.__import__

    def no_umap(name, *a, **k):
        if name == "umap":
            raise ImportError("no umap")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_umap)
    msgs = []
    assert TU.walk_umap(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)),
                        log=msgs.append).empty
    assert any("umap-learn not installed" in m for m in msgs)


# --------------------------------------------------------------------------- the store
def test_an_embedding_round_trips_with_its_spec(tmp_path):
    """A result that cannot be reopened with the settings that produced it is not reproducible."""
    store = TU.EmbeddingStore(str(tmp_path))
    spec = EmbeddingSpec(blocks=("fitness_screens",), na_policy="median")
    Y = np.random.default_rng(0).normal(size=(10, 3))
    store.save("run1", Y, spec, gene_ids=[f"g{i}" for i in range(10)], features=["a", "b"])
    got, spec_back, meta = store.load("run1")
    assert np.allclose(got, Y)
    assert spec_back.na_policy == "median"
    assert meta["n_genes"] == 10


def test_loading_something_that_was_never_saved_is_an_explicit_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        TU.EmbeddingStore(str(tmp_path)).load("never_saved")


def test_the_store_lists_what_it_holds(tmp_path):
    store = TU.EmbeddingStore(str(tmp_path))
    Y = np.zeros((5, 3))
    store.save("a", Y, EmbeddingSpec())
    store.save("b", Y, EmbeddingSpec())
    assert sorted(store.list().name) == ["a", "b"]


def test_an_empty_store_lists_nothing(tmp_path):
    assert TU.EmbeddingStore(str(tmp_path)).list().empty


# --------------------------------------------------------------------------- guessing the id column
def test_the_identifier_column_is_recognised_wherever_it_sits():
    df = pd.DataFrame({"description": ["a", "b"], "locus": ["TGME49_200010", "TGME49_200020"]})
    col, frac, needs_repair = TU.suggest_gene_column(df)
    assert col == "locus" and frac > 0.9 and not needs_repair


def test_a_malformed_identifier_column_is_still_suggested_with_a_repair_flag():
    """A malformed column is precisely the case where the user needs a suggestion, so refusing to
    guess would be unhelpful."""
    df = pd.DataFrame({"locus": ["TgME49.208830", "TgME49.208840"]})
    col, frac, needs_repair = TU.suggest_gene_column(df)
    assert col == "locus" and needs_repair


def test_a_table_with_no_identifier_column_suggests_nothing():
    assert TU.suggest_gene_column(pd.DataFrame({"a": [1, 2], "b": ["x", "y"]}))[0] is None


def test_importing_a_table_with_no_identifier_column_is_an_explicit_error():
    with pytest.raises(ValueError, match="no column looks like a gene identifier"):
        TU.import_table(pd.DataFrame({"a": [1, 2]}), log=lambda *_: None)


# --------------------------------------------------------------------------- repairing identifiers
@pytest.mark.parametrize("raw", ["TgME49.208830", "TGME49-208830", "tgme49 208830",
                                 "TGME49208830", "gene|TGME49_208830|v2"])
def test_a_repair_is_suggested_for_each_mangled_form(raw):
    """Every one of these appears in real published tables, and each matches nothing unrepaired."""
    s = pd.Series([raw] * 10)
    sug = TU.suggest_repair(s)
    assert sug is not None, raw
    pat, rep = sug
    assert s.str.replace(pat, rep, regex=True).str.contains("TGME49_208830").all()


def test_no_repair_is_suggested_for_something_that_is_not_an_identifier():
    """Suggesting a pattern for a column of protein names would produce confident nonsense."""
    assert TU.suggest_repair(pd.Series(["alpha", "beta", "gamma"] * 5)) is None


def test_already_correct_identifiers_need_no_repair():
    s = pd.Series(["TGME49_208830"] * 10)
    sug = TU.suggest_repair(s)
    if sug is not None:
        pat, rep = sug
        assert s.str.replace(pat, rep, regex=True).equals(s)


# --------------------------------------------------------------------------- importing
def test_an_imported_table_is_keyed_on_canonical_gene_ids():
    df = pd.DataFrame({"id": ["TGME49_200010", "TGME49_200020"], "score": [1.0, 2.0]})
    out = TU.import_table(df, "id", log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010", "TGME49_200020"]


def test_an_identifier_column_can_be_repaired_on_the_way_in():
    df = pd.DataFrame({"id": ["TgME49.208830"], "score": [1.0]})
    out = TU.import_table(df, "id", pattern=r"(?i)tgme49\s*[.\- ]\s*(\d{5,6})",
                          replacement=r"TGME49_\1", log=lambda *_: None)
    assert list(out.index) == ["TGME49_208830"]


def test_imported_identifiers_go_through_the_identity_layer():
    """Previous and strain accessions map forward -- the failure that cost a published screen every
    one of its rows."""
    df = pd.DataFrame({"id": ["TGGT1_200010"], "score": [1.0]})
    out = TU.import_table(df, "id", resolve=lambda a: a.replace("TGGT1", "TGME49"),
                          log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_a_table_that_matches_nothing_says_so_rather_than_returning_silence():
    """An empty result looks exactly like a file with no genes in it."""
    df = pd.DataFrame({"id": ["not_a_gene", "also_not"], "score": [1.0, 2.0]})
    msgs = []
    out = TU.import_table(df, "id", log=msgs.append)
    assert out.empty
    assert any("nothing matched" in m for m in msgs)


def test_how_many_rows_resolved_is_always_reported():
    df = pd.DataFrame({"id": ["TGME49_200010", "junk"], "score": [1.0, 2.0]})
    msgs = []
    TU.import_table(df, "id", log=msgs.append)
    assert any("1 of 2 rows resolved" in m for m in msgs)


def test_a_csv_path_can_be_imported_directly(tmp_path):
    p = tmp_path / "mine.csv"
    p.write_text("id,score\nTGME49_200010,1.5\n")
    out = TU.import_table(str(p), "id", log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_imported_numeric_columns_are_prefixed_so_they_cannot_collide():
    """A user column called `compartment` must not silently overwrite the measured one."""
    df = pd.DataFrame({"id": ["TGME49_200010"], "score": [1.0]})
    out = TU.import_table(df, "id", prefix="mine_", log=lambda *_: None)
    assert any(c.startswith("mine_") for c in out.columns)


def test_a_trustworthiness_that_cannot_be_computed_is_missing_not_zero(monkeypatch):
    """Reporting a fabricated 0.0 would rank a degenerate setting against real ones."""
    import sklearn.manifold
    monkeypatch.setattr(sklearn.manifold, "trustworthiness",
                        lambda *a, **k: (_ for _ in ()).throw(ValueError("degenerate")))
    out = TU._quality(np.random.default_rng(0).normal(size=(60, 4)),
                      np.random.default_rng(1).normal(size=(60, 2)), n_neighbors=10)
    assert np.isnan(out["trustworthiness"])


def test_a_distance_correlation_that_cannot_be_computed_is_missing_not_zero(monkeypatch):
    import scipy.stats
    monkeypatch.setattr(scipy.stats, "spearmanr",
                        lambda *a, **k: (_ for _ in ()).throw(ValueError("degenerate")))
    out = TU._quality(np.random.default_rng(0).normal(size=(60, 4)),
                      np.random.default_rng(1).normal(size=(60, 2)), n_neighbors=10)
    assert np.isnan(out["continuity_proxy"])


def test_a_walk_where_every_setting_is_skipped_reports_rather_than_crashing():
    """Built from an empty list the frame has no columns, and sorting on one raised KeyError."""
    out = TU.walk_umap(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)),
                       n_neighbors_values=(5000,), min_dist_values=(0.1,),
                       sample_size=120, log=lambda *_: None)
    assert out.empty
    assert "trustworthiness" in out.columns


def test_a_reopened_embedding_knows_which_gene_each_row_is(tmp_path):
    """They were saved into the npz and not returned, so a reopened embedding came back as coordinates
    with no way to say which gene each row is -- useless for anything gene-specific, including turning a
    recovered structure into named predictions, which is the reason for saving embeddings at all."""
    store = TU.EmbeddingStore(str(tmp_path))
    genes = [f"TGME49_{200000+i}" for i in range(10)]
    store.save("run", np.zeros((10, 3)), EmbeddingSpec(), gene_ids=genes, features=["a"])
    _, _, meta = store.load("run")
    assert meta["gene_ids"] == genes


def test_an_embedding_saved_without_gene_ids_still_loads(tmp_path):
    store = TU.EmbeddingStore(str(tmp_path))
    store.save("run", np.zeros((5, 3)), EmbeddingSpec())
    Y, _, meta = store.load("run")
    assert len(Y) == 5
    assert "gene_ids" not in meta
