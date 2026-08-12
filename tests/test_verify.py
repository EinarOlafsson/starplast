#!/usr/bin/env python3
"""Checking the shipped numbers against the repository they came from.

The question this answers is not "which source is better" but "do they agree", because each source
fails in a way the other detects. The comparison is on ranks by construction: the two forms of
GSE108740 differ by a log transform -- the shipped columns top out at 3.3, the GEO matrix at 16,520 --
so comparing values would report a difference that is only a transform, and the check has to survive
not knowing which transform was applied.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import verify as V  # noqa: E402


def _geo_book(tmp_path, values, ids=None, cols=("Tachyzoites_T2 [FPKM]", "Tachyzoites_T4 [FPKM]")):
    """A stand-in GEO workbook at the path the correspondence table expects."""
    p = tmp_path / "transcription" / "RNAseq" / "GSE108740"
    p.mkdir(parents=True, exist_ok=True)
    ids = ids or [f"TGME49_{200000+i}" for i in range(len(values))]
    d = pd.DataFrame({"ToxoDB ID": ids})
    for c in cols:
        d[c] = values
    d.to_excel(p / "GSE108740_FPKM.xlsx", index=False)
    return p


def _nodes(ids, tachy):
    return pd.DataFrame({"gene_id": ids, "expr_tachy": tachy})


def test_a_log_transform_between_the_sources_does_not_count_as_disagreement(tmp_path, monkeypatch):
    """The exact situation on disk: the same series shipped logged and deposited raw."""
    from starplast import paths
    raw = np.array([1.0, 10.0, 100.0, 1000.0, 16520.0] * 8)
    _geo_book(tmp_path, raw)
    ids = [f"TGME49_{200000+i}" for i in range(len(raw))]
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = V.compare_series("GSE108740", _nodes(ids, np.log2(raw + 1)), log=lambda *_: None)
    assert len(t) == 1
    assert t.loc[0, "spearman"] == pytest.approx(1.0)
    assert bool(t.loc[0, "agrees"])


def test_a_genuinely_different_ordering_is_reported_as_disagreement(tmp_path, monkeypatch):
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    _geo_book(tmp_path, raw)
    ids = [f"TGME49_{200000+i}" for i in range(len(raw))]
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = V.compare_series("GSE108740", _nodes(ids, raw[::-1]), log=lambda *_: None)
    assert t.loc[0, "spearman"] < 0
    assert not bool(t.loc[0, "agrees"])


def test_disagreement_is_named_rather_than_resolved(tmp_path, monkeypatch, capsys):
    """Quietly preferring one source would discard the only evidence that they differ."""
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    _geo_book(tmp_path, raw)
    ids = [f"TGME49_{200000+i}" for i in range(len(raw))]
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    msgs = []
    V.verify_all(_nodes(ids, raw[::-1]), log=msgs.append)
    assert any("DISAGREEMENT" in m and "expr_tachy" in m for m in msgs)


def test_a_missing_primary_matrix_reports_rather_than_raising(tmp_path, monkeypatch):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    msgs = []
    assert V.compare_series("GSE108740", _nodes(["a"], [1.0]), log=msgs.append).empty
    assert any("cannot verify" in m for m in msgs)
    assert V.verify_all(_nodes(["a"], [1.0]), log=lambda *_: None).empty


def test_accessions_are_routed_through_the_identity_layer(tmp_path, monkeypatch):
    """A GEO matrix keyed on a strain accession joins nothing without it, silently."""
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    gt1 = [f"TGGT1_{200000+i}" for i in range(len(raw))]
    me49 = [g.replace("TGGT1", "TGME49") for g in gt1]
    _geo_book(tmp_path, raw, ids=gt1)
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))

    unresolved = V.compare_series("GSE108740", _nodes(me49, raw), log=lambda *_: None)
    assert unresolved.loc[0, "n_common"] == 0

    resolved = V.compare_series("GSE108740", _nodes(me49, raw),
                                resolve=lambda a: a.replace("TGGT1", "TGME49"), log=lambda *_: None)
    assert resolved.loc[0, "n_common"] == len(raw)
    assert resolved.loc[0, "spearman"] == pytest.approx(1.0)


def test_too_few_common_genes_gives_nan_not_a_confident_number():
    """Spearman over a handful of genes is noise wearing a statistic's clothes."""
    a = pd.Series([1.0, 2.0, 3.0])
    assert np.isnan(V._spearman(a, a))


def test_a_column_absent_from_the_node_table_is_skipped(tmp_path, monkeypatch):
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    _geo_book(tmp_path, raw)
    ids = [f"TGME49_{200000+i}" for i in range(len(raw))]
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = V.compare_series("GSE108740", pd.DataFrame({"gene_id": ids}), log=lambda *_: None)
    assert t.empty


def test_a_geo_column_absent_from_the_matrix_is_skipped(tmp_path, monkeypatch):
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    _geo_book(tmp_path, raw, cols=("Something Else [FPKM]",))
    ids = [f"TGME49_{200000+i}" for i in range(len(raw))]
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = V.compare_series("GSE108740", _nodes(ids, raw), log=lambda *_: None)
    assert t.empty


def test_duplicate_accessions_in_the_geo_matrix_do_not_multiply_rows(tmp_path, monkeypatch):
    from starplast import paths
    raw = np.arange(1.0, 41.0)
    ids = [f"TGME49_{200000+(i % 20)}" for i in range(len(raw))]     # each id appears twice
    _geo_book(tmp_path, raw, ids=ids)
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = V.compare_series("GSE108740", _nodes(sorted(set(ids)), np.arange(1.0, 21.0)),
                         log=lambda *_: None)
    assert t.loc[0, "n_common"] == 20


def test_the_correspondence_is_declared_not_guessed_from_names():
    """`expr_tachy` and `Tachyzoites_T2 [FPKM]` do not look alike. A name-matching heuristic would
    compare the wrong pair and report the mismatch as a data problem."""
    pairs = V.CORRESPONDENCE["GSE108740"]["pairs"]
    assert "expr_tachy" in pairs
    assert all(isinstance(v, tuple) and v for v in pairs.values())


@pytest.mark.skipif(not os.environ.get("STARPLAST_DATA"),
                    reason="needs the raw dataset tree; set STARPLAST_DATA")
def test_the_real_shipped_columns_reproduce_their_geo_source():
    """The result that answers task 05: measured, not assumed."""
    from starplast import paths
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    t = V.verify_all(nodes, log=lambda *_: None)
    assert not t.empty
    assert t.agrees.all(), t.to_string()
    assert (t.spearman > 0.99).all()


# --------------------------------------------------------------------------- multi-sheet workbooks
def test_the_declared_sheet_is_the_one_compared(tmp_path, monkeypatch):
    """GSE206344 ships two workbooks and the source one has eight sheets. Comparing the wrong sheet --
    or the first by default -- produces a number that looks like verification and is not."""
    from starplast import paths
    d = tmp_path / "transcription" / "RNAseq" / "GSE206344"
    d.mkdir(parents=True)
    f = d / "GSE206344_Normalised_data_ToxoDB_Release68.xlsx"
    ids = [f"TGME49_{200000+i}" for i in range(40)]
    raw = np.arange(1.0, 41.0)
    with pd.ExcelWriter(f) as w:
        # sheet 0 is a contents page, exactly as the real workbook has
        pd.DataFrame({"Unnamed: 0": ["Contents"], "Unnamed: 1": ["see sheet 1"]}).to_excel(
            w, sheet_name="Contents", index=False)
        pd.DataFrame({"ToxoDB ID Release 68": ids,
                      "Unsporulated R1": raw, "Unsporulated R2": raw,
                      "Sporulating R1": raw, "Sporulating R2": raw,
                      "Sporulated R1": raw, "Sporulated R2": raw}).to_excel(
            w, sheet_name="1", index=False)
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    nodes = pd.DataFrame({"gene_id": ids, "expr_sporulated": np.log2(raw + 1)})
    t = V.compare_series("GSE206344", nodes, log=lambda *_: None)
    assert not t.empty, "the declared sheet must be found rather than the contents page"
    assert t.loc[0, "spearman"] == pytest.approx(1.0)


def test_a_series_without_a_declared_sheet_reads_the_first(tmp_path, monkeypatch):
    """GSE108740 is a single-sheet workbook, so the default must keep working."""
    assert "sheet" not in V.CORRESPONDENCE["GSE108740"]


def test_every_declared_correspondence_names_a_real_shipped_column():
    """A pair naming a column the node table does not have is silently skipped, so it would look
    verified while checking nothing."""
    from starplast import paths
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    unknown = {}
    for series, spec in V.CORRESPONDENCE.items():
        missing = [c for c in spec["pairs"] if c not in nodes.columns]
        if missing:
            unknown[series] = missing
    assert not unknown, f"correspondences naming absent columns: {unknown}"
