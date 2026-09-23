#!/usr/bin/env python3
"""Importing a user's own table, with the preprocessing offered rather than assumed.

The rule this module exists to enforce is the one written on the wall of `build_graph`: THE RANGE
DECIDES, NOT THE FILENAME. The same GEO series demonstrates both ways to get it wrong -- its FPKM
file reaches 16,520 and is real FPKM, while the columns derived from it in the node table are still
called `_FPKM` and top out at 9.7 because they were logged upstream. Take the log twice and real
variation compresses to nothing; skip it and one gene dominates every distance. Neither is visible
in the output, which is why the numbers are shown and the guess comes with its reason attached.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import importer as I  # noqa: E402


def _table(n=40, form="fpkm"):
    rng = np.random.default_rng(0)
    ids = [f"TGME49_{200000 + i}" for i in range(n)]
    if form == "fpkm":
        vals = rng.gamma(2, 300, size=(n, 3))
    elif form == "logged":
        vals = rng.uniform(0, 12, size=(n, 3))
    else:
        vals = rng.normal(0, 2, size=(n, 3))
    d = pd.DataFrame(vals, columns=["day3", "day5", "day7"])
    d.insert(0, "gene", ids)
    return d


# --------------------------------------------------------------------------- reading
def test_it_reads_the_formats_this_project_already_reads(tmp_path):
    d = _table()
    for name, write in (("t.csv", lambda p: d.to_csv(p, index=False)),
                        ("t.tsv", lambda p: d.to_csv(p, sep="\t", index=False)),
                        ("t.parquet", lambda p: d.to_parquet(p)),
                        ("t.xlsx", lambda p: d.to_excel(p, index=False))):
        path = tmp_path / name
        write(path)
        got = I.read_any(str(path))
        assert list(got.columns)[-3:] == ["day3", "day5", "day7"], name
        assert len(got) == len(d), name


def test_a_txt_that_is_really_a_tsv_is_read_as_one(tmp_path):
    """Half the supplementary tables in this field are .txt that are really TSV, and guessing wrongly
    gives one column with everything in it -- which looks like a file with no data rather than a file
    read badly."""
    path = tmp_path / "supp.txt"
    _table().to_csv(path, sep="\t", index=False)
    assert I.read_any(str(path)).shape[1] == 4


def test_a_workbook_says_which_sheets_it_has(tmp_path):
    """Published supplements routinely put the table on sheet 3 behind a legend."""
    path = tmp_path / "book.xlsx"
    with pd.ExcelWriter(path) as w:
        pd.DataFrame({"note": ["read me"]}).to_excel(w, sheet_name="legend", index=False)
        _table().to_excel(w, sheet_name="data", index=False)
    assert I.sheets_in(str(path)) == ["legend", "data"]
    assert list(I.read_any(str(path), sheet="data").columns)[0] == "gene"
    assert I.sheets_in(str(tmp_path / "t.csv")) == []


def test_something_that_is_not_a_workbook_lists_no_sheets(tmp_path, capsys):
    path = tmp_path / "broken.xlsx"
    path.write_bytes(b"not a workbook")
    assert I.sheets_in(str(path)) == []


# --------------------------------------------------------------------------- the range decides
def test_the_columns_are_described_before_anything_is_transformed():
    d = I.describe_columns(_table())
    assert list(d.column) == ["day3", "day5", "day7"]
    assert (d["max"] > d["median"]).all() and (d.missing == 0).all()


def test_unlogged_values_are_recognised_by_their_range_not_their_name():
    """`GSE108740_FPKM` reaches 16,520 and is real FPKM. The name is not the evidence."""
    quant, why = I.suggest_quantification(_table(form="fpkm"))
    assert quant == "fpkm" and "has NOT been logged" in why


def test_already_logged_values_are_recognised_too():
    """The columns derived from that same file top out at 9.7 because they were logged upstream.
    Logging them again compresses real variation to nothing."""
    quant, why = I.suggest_quantification(_table(form="logged"))
    assert quant == "log_intensity" and "already logged" in why


def test_a_fold_change_is_recognised_by_its_negatives():
    quant, why = I.suggest_quantification(_table(form="lfc"))
    assert quant == "lfc" and "negative" in why


def test_a_table_with_no_numbers_says_so_rather_than_guessing():
    quant, why = I.suggest_quantification(pd.DataFrame({"gene": ["TGME49_200000"]}))
    assert quant == "none" and "no numeric columns" in why


def test_a_guess_always_comes_with_its_reason():
    """"fpkm" alone is not something a person can agree or disagree with."""
    for form in ("fpkm", "logged", "lfc"):
        quant, why = I.suggest_quantification(_table(form=form))
        assert len(why.split()) > 5, (form, why)


# --------------------------------------------------------------------------- preprocessing
def test_an_import_records_every_choice_that_produced_it():
    """A column whose provenance is a memory of which dropdowns were set is a column nobody can
    defend three weeks later."""
    out, record = I.preprocess(_table(), gene_column="gene", quantification="fpkm",
                               scaling="rank", na_policy="median", log=lambda *_: None)
    assert record["quantification"] == "fpkm" and record["scaling"] == "rank"
    assert record["na_policy"] == "median" and record["duplicates"] == "mean"
    assert record["rows_in"] == 40 and record["rows_resolved"] == 40
    assert record["columns"] == list(out.columns) and record["genes"] == len(out)


def test_imported_columns_are_prefixed_so_they_cannot_pass_as_measurements():
    out, _ = I.preprocess(_table(), gene_column="gene", log=lambda *_: None)
    assert all(c.startswith("imported_") for c in out.columns)
    assert out.index.name == "gene_id" or out.index.str.startswith("TGME49_").all()


def test_the_quantification_transform_is_the_one_the_rest_of_the_program_uses():
    """`sources.normalize`, not a second implementation of log-and-centre."""
    d = _table(form="fpkm")
    out, _ = I.preprocess(d, gene_column="gene", quantification="fpkm", log=lambda *_: None)
    assert out.to_numpy().max() < 20, "the log was not taken"
    plain, _ = I.preprocess(d, gene_column="gene", quantification="none", log=lambda *_: None)
    assert plain.to_numpy().max() > 100


@pytest.mark.parametrize("scaling", ["robust", "zscore", "rank", "none"])
def test_every_scaling_the_embedding_offers_is_offered_here(scaling):
    out, record = I.preprocess(_table(), gene_column="gene", scaling=scaling,
                               log=lambda *_: None)
    assert record["scaling"] == scaling and np.isfinite(out.to_numpy()).any()
    if scaling == "rank":
        assert out.to_numpy().min() >= -0.5001 and out.to_numpy().max() <= 0.5001


def test_a_sign_flip_is_offered_with_the_warning_that_goes_with_it():
    """Screens disagree about which sign means worse, and pooling them without rank-normalizing
    first is how the 64x spread between them bites."""
    said = []
    out, record = I.preprocess(_table(), gene_column="gene", flip=True, log=said.append)
    assert record["flip"] and any("rank-normalize before pooling" in m for m in said)
    plain, _ = I.preprocess(_table(), gene_column="gene", log=lambda *_: None)
    assert np.allclose(out.to_numpy(), -plain.to_numpy())


@pytest.mark.parametrize("how", list(I.DUPLICATES))
def test_several_rows_for_one_gene_are_combined_the_way_the_user_chose(how):
    """A hit table often lists one row per guide or per peptide, and which of them is the claim is
    the user's call rather than the importer's."""
    d = pd.DataFrame({"gene": ["TGME49_200000"] * 3 + ["TGME49_200001"],
                      "score": [1.0, 5.0, 9.0, 2.0]})
    out, record = I.preprocess(d, gene_column="gene", duplicates=how, log=lambda *_: None)
    got = float(out.loc["TGME49_200000", "imported_score"])
    assert got == {"mean": 5.0, "median": 5.0, "max": 9.0, "min": 1.0, "first": 1.0}[how]
    assert record["duplicates"] == how


def test_missing_values_follow_the_policy_that_was_chosen():
    d = _table()
    d.loc[0:9, "day3"] = np.nan
    dropped, _ = I.preprocess(d, gene_column="gene", na_policy="drop_genes", log=lambda *_: None)
    assert len(dropped) == 30
    kept, _ = I.preprocess(d, gene_column="gene", na_policy="median", log=lambda *_: None)
    assert len(kept) == 40 and not kept.isna().any().any()


def test_a_column_that_is_mostly_missing_can_be_dropped(capsys):
    d = _table()
    d.loc[0:35, "day3"] = np.nan
    out, _ = I.preprocess(d, gene_column="gene", na_policy="drop_columns", log=lambda *_: None)
    assert "imported_day3" not in out.columns and "imported_day5" in out.columns


def test_indicator_is_left_to_the_embedding_rather_than_done_twice():
    """The embedding adds missingness indicators itself, weighted below a measurement. Adding them
    here as well would count absence twice."""
    out, _ = I.preprocess(_table(), gene_column="gene", na_policy="indicator",
                          log=lambda *_: None)
    assert not any("missing" in c for c in out.columns)


# --------------------------------------------------------------------------- identifiers
def test_a_malformed_identifier_column_is_repaired_rather_than_refused():
    """`TgME49.208830` and `gene|TGME49_208830|v2` are both real formats from published supplements,
    and both match nothing unrepaired."""
    d = pd.DataFrame({"id": [f"TgME49.{200000 + i}" for i in range(6)],
                      "score": np.arange(6.0)})
    said = []
    out, record = I.preprocess(d, gene_column="id", log=said.append)
    assert len(out) == 6 and record["repair"]
    assert any("reshaped" in m for m in said)


def test_a_table_whose_identifiers_resolve_to_nothing_says_so():
    d = pd.DataFrame({"id": ["not-a-gene"] * 5, "score": np.arange(5.0)})
    with pytest.raises(ValueError, match="no row resolved"):
        I.preprocess(d, gene_column="id", log=lambda *_: None)


def test_a_table_with_no_identifier_column_at_all_says_so():
    with pytest.raises(ValueError, match="no column looks like a gene identifier"):
        I.preprocess(pd.DataFrame({"a": [1, 2], "b": [3, 4]}), log=lambda *_: None)


def test_a_table_with_no_numbers_says_so():
    d = pd.DataFrame({"gene": ["TGME49_200000"], "note": ["interesting"]})
    with pytest.raises(ValueError, match="no numeric columns"):
        I.preprocess(d, gene_column="gene", log=lambda *_: None)


def test_previous_and_strain_accessions_go_through_the_identity_layer():
    """The failure that cost a published screen every one of its rows: a 2019 in vivo screen uses
    pre-2012 ids for EVERY gene and contributed 0 of 8,140 until it was resolved."""
    d = pd.DataFrame({"id": ["TGGT1_208830", "TGGT1_208840"], "score": [1.0, 2.0]})
    out, _ = I.preprocess(d, gene_column="id", resolve=lambda x: x.replace("TGGT1_", "TGME49_"),
                          log=lambda *_: None)
    assert list(out.index) == ["TGME49_208830", "TGME49_208840"]


# --------------------------------------------------------------------------- joining
def test_imported_columns_join_the_table_in_memory_and_never_the_file():
    """The cache is measurement that shipped with the program. An imported column written into the
    same table on disk becomes indistinguishable from one."""
    from starplast import paths
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    imported = pd.DataFrame({"imported_score": [1.0, 2.0]},
                            index=[nodes.gene_id.iloc[0], nodes.gene_id.iloc[5]])
    out = I.merge_into(nodes, imported)
    assert "imported_score" not in nodes.columns, "the node table was modified in place"
    assert out.imported_score.iloc[0] == 1.0 and out.imported_score.iloc[5] == 2.0
    assert np.isnan(out.imported_score.iloc[1]), "a gene not in the file must be missing, not zero"
    after = pd.read_parquet(paths.cache_file("nodes.parquet"))
    assert list(after.columns) == list(nodes.columns), "the cache on disk changed"


def test_numbers_that_fit_no_pattern_are_left_alone_and_say_so():
    """Better "no transform is obviously right" than a confident guess about someone else's data."""
    d = pd.DataFrame({"gene": [f"TGME49_{200000 + i}" for i in range(20)],
                      "x": np.linspace(0.001, 0.05, 20)})
    quant, why = I.suggest_quantification(d)
    assert quant == "none" and "no transform is obviously right" in why


def test_resolving_previous_accessions_survives_a_duplicate_rule_other_than_the_default():
    """The identity layer has to be applied on the path that combines duplicates too, or a table
    keyed on strain accessions resolves under `mean` and vanishes under `max`."""
    d = pd.DataFrame({"id": ["TGGT1_208830", "TGGT1_208830", "TGGT1_208840"],
                      "score": [1.0, 9.0, 4.0]})
    out, _ = I.preprocess(d, gene_column="id", duplicates="max",
                          resolve=lambda x: x.replace("TGGT1_", "TGME49_"),
                          log=lambda *_: None)
    assert list(out.index) == ["TGME49_208830", "TGME49_208840"]
    assert float(out.loc["TGME49_208830", "imported_score"]) == 9.0


def test_every_format_the_importer_offers_has_a_reader_declared():
    """`.xls` was offered and its reader was not declared, which is an ImportError with no author.

    The gap surfaced when a Plasmodium supplement in that format was registered and the loader ran
    on the interpreter the user actually uses. It is the same class of fault either way: the program
    says it can read a format, and whether it can depends on what else happens to be installed.
    """
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    declared = open(os.path.join(root, "pyproject.toml"), encoding="utf8").read()
    readers = {".xlsx": "openpyxl", ".xls": "xlrd", ".parquet": "pyarrow"}
    for ext, package in readers.items():
        assert ext in I.READABLE, f"{ext} stopped being offered; drop it from this test too"
        assert package in declared, f"{ext} is offered but {package} is not a declared dependency"


def test_a_screen_export_can_use_gene_id_as_its_identifier_column():
    """A canonical column name must survive identifier resolution and duplicate aggregation."""
    from starplast.importer import preprocess
    frame = pd.DataFrame({"gene_id": ["TGME49_208830", "tgme49_208830", "TGME49_200000", "unknown"],
                          "score": [2.0, 4.0, 8.0, 10.0]})
    before = frame.copy(deep=True)
    result, record = preprocess(frame, gene_column="gene_id", columns=["score"],
                                prefix="screen_", log=lambda *_: None)
    assert result.loc["TGME49_208830", "screen_score"] == 3.0
    assert result.loc["TGME49_200000", "screen_score"] == 8.0
    assert len(result) == 2
    assert record["gene_column"] == "gene_id"
    pd.testing.assert_frame_equal(frame, before)
