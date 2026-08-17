#!/usr/bin/env python3
"""The arrayed imaging screen, and the letter that had to be earned rather than assumed.

The source table scores each gene with codes like `E3` and `F1 A2` and ships no legend -- it lives in
a figure that is an image. Reading `E` as egress is the whole value of the column, so it is not
asserted here from the initial. The paper names exactly two genes as the egress mutants it went on
to characterise, CGP and SLF, and the test is that those two come out with an egress phenotype. If
the mapping were wrong that assertion fails, which is the point of writing it down.

The other thing worth failing over is missingness. 319 genes were looked at and 99 had a phenotype;
a screened gene with no egress call was examined and found normal, and the other 7,800 genes were
never examined. Both would arrive as an empty cell, and collapsing them would tell the map that
almost every gene in Toxoplasma has been checked for an egress defect and passed.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import phenotype_screen as S  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE = os.path.join(ROOT, "datasets", "DNA", "imaging_screen", "35538310", S.SOURCE)

#: The two genes the paper names as its egress mutants: CGP and SLF.
NAMED_EGRESS = ("TGGT1_240380", "TGGT1_208420")


def _book(tmp_path, library=None, calls=None, sheets=True):
    """A workbook shaped like the real supplementary tables."""
    library = library if library is not None else ["TGGT1_100010", "TGGT1_100020", "TGGT1_100030"]
    calls = calls if calls is not None else [("100010", "a", "x1", "E3", "E₂"),
                                             ("100020", "b", "x2", "F₁ A2 / R3", None)]
    lib = pd.DataFrame({"Gene ID": library, "Annotation in ToxoDB": ["x"] * len(library)})
    call = pd.DataFrame(calls, columns=["Accession number", "Name", "No. of times picked",
                                        "Investigator 1", "Investigator 2"])
    path = tmp_path / "book.xlsx"
    with pd.ExcelWriter(path) as writer:
        name_a = S.LIBRARY[0] if sheets else "other"
        name_b = S.CALLS[0] if sheets else "elsewhere"
        pd.DataFrame([["Supplementary Table 2"], [""]]).to_excel(
            writer, sheet_name=name_a, index=False, header=False)
        lib.to_excel(writer, sheet_name=name_a, index=False, startrow=2)
        pd.DataFrame([[""]] * 4).to_excel(writer, sheet_name=name_b, index=False, header=False)
        call.to_excel(writer, sheet_name=name_b, index=False, startrow=4)
    return str(path)


def _base(tmp_path, **kw):
    """A dataset tree with the workbook where `screen` looks for it."""
    folder = tmp_path / "datasets" / "phenotype"
    folder.mkdir(parents=True)
    path = _book(tmp_path, **kw)
    os.rename(path, folder / S.SOURCE)
    return str(tmp_path)


# --------------------------------------------------------------------------- reading a cell
def test_a_category_letter_is_read_and_its_subscript_is_not():
    assert S.categories("E3") == {"egress"}
    assert S.categories("E₄") == {"egress"}
    assert S.categories("E3") == S.categories("E4"), "the digit must not change the category"


def test_several_clones_in_one_cell_are_all_read():
    assert S.categories("F₁ A2 / R3") == {"actin", "apicoplast", "replication"}


def test_a_cell_with_no_call_is_no_categories():
    assert S.categories(None) == set() and S.categories(float("nan")) == set()
    assert S.categories("") == set() and S.categories("no phenotype") == set()


def test_an_unknown_letter_is_ignored_rather_than_invented():
    assert S.categories("Z9 E1") == {"egress"}


# --------------------------------------------------------------------------- the accession
def test_a_gene_is_recognised_written_either_way():
    assert S._accession("201270") == "TGGT1_201270"
    assert S._accession("TGGT1_201270") == "TGGT1_201270"
    assert S._accession(201270.0) == "TGGT1_201270"


def test_a_cell_that_is_not_a_gene_is_refused():
    assert S._accession(None) is None and S._accession(float("nan")) is None
    assert S._accession("Accession number") is None


# --------------------------------------------------------------------------- missingness
def test_a_screened_gene_with_no_phenotype_is_false_and_not_missing(tmp_path):
    """It was looked at. That is a measurement, and it is the commonest result in any screen."""
    d = S.screen(_base(tmp_path), log=lambda *a: None)
    assert d.loc["TGGT1_100030", "screen_egress_phenotype"] is False or \
        not d.loc["TGGT1_100030", "screen_egress_phenotype"]
    assert not d.loc["TGGT1_100030", "screen_any_phenotype"]
    assert len(d) == 3


def test_a_gene_that_was_never_screened_is_absent_entirely(tmp_path):
    d = S.screen(_base(tmp_path), log=lambda *a: None)
    assert "TGGT1_999999" not in d.index


def test_the_categories_land_in_their_own_columns(tmp_path):
    d = S.screen(_base(tmp_path), log=lambda *a: None)
    assert d.loc["TGGT1_100010", "screen_egress_phenotype"]
    assert not d.loc["TGGT1_100010", "screen_actin_phenotype"]
    assert d.loc["TGGT1_100020", "screen_actin_phenotype"]
    assert d.loc["TGGT1_100020", "screen_replication_phenotype"]


def test_agreement_is_only_defined_where_both_investigators_scored(tmp_path):
    d = S.screen(_base(tmp_path), log=lambda *a: None)
    assert d.loc["TGGT1_100010", "screen_scorers_agree"]        # both said egress
    assert d.loc["TGGT1_100020", "screen_scorers_agree"] is None  # only one scored
    assert d.loc["TGGT1_100030", "screen_scorers_agree"] is None  # neither did


def test_a_library_listing_one_gene_twice_yields_one_row(tmp_path):
    """The real library has 320 gRNA rows for 319 genes."""
    d = S.screen(_base(tmp_path, library=["TGGT1_100010", "TGGT1_100010", "TGGT1_100020"]),
                 log=lambda *a: None)
    assert len(d) == 2 and d.index.is_unique


def test_accessions_are_resolved_when_a_resolver_is_given(tmp_path):
    d = S.screen(_base(tmp_path), log=lambda *a: None,
                 resolve=lambda g: g.replace("TGGT1_", "TGME49_"))
    assert d.index.str.startswith("TGME49_").all()


def test_two_strain_accessions_collapsing_to_one_gene_yield_one_row(tmp_path):
    d = S.screen(_base(tmp_path), log=lambda *a: None, resolve=lambda g: "TGME49_1")
    assert len(d) == 1


# --------------------------------------------------------------------------- refusals
def test_a_missing_source_yields_nothing(tmp_path):
    said = []
    assert S.screen(str(tmp_path), log=said.append).empty
    assert said


def test_a_workbook_without_both_sheets_is_refused(tmp_path):
    assert S.screen(_base(tmp_path, sheets=False), log=lambda *a: None).empty


def test_a_library_sheet_with_no_genes_is_refused(tmp_path):
    assert S.screen(_base(tmp_path, library=["not a gene"]), log=lambda *a: None).empty


def test_a_call_row_that_names_no_gene_is_skipped(tmp_path):
    d = S.screen(_base(tmp_path, calls=[(None, "x", "x1", "E3", None),
                                        ("100010", "a", "x1", "E3", None)]),
                 log=lambda *a: None)
    assert d.loc["TGGT1_100010", "screen_egress_phenotype"]


# --------------------------------------------------------------------------- the real source
@pytest.mark.skipif(not os.path.exists(SOURCE), reason="screen supplement not fetched")
def test_the_egress_code_means_egress():
    """The mapping the whole column rests on, checked against the paper's own named mutants.

    CGP and SLF are the two genes this study identified as egress mutants and went on to
    characterise. Both must carry an egress phenotype here. If `E` meant apicoplast, or if the
    subscripts were being read as the category, this fails.
    """
    d = S.screen(ROOT, log=lambda *a: None)
    for gene in NAMED_EGRESS:
        assert gene in d.index, gene
        assert d.loc[gene, "screen_egress_phenotype"], f"{gene} is a named egress mutant"


@pytest.mark.skipif(not os.path.exists(SOURCE), reason="screen supplement not fetched")
def test_the_real_screen_has_the_shape_the_paper_reports():
    d = S.screen(ROOT, log=lambda *a: None)
    assert len(d) == 319, "320 gRNA rows describe 319 genes"
    assert int(d["screen_any_phenotype"].sum()) == 99, "99 genes were picked with a phenotype"
    # Every gene with a category must also read as having a phenotype, or the summary column and
    # the category columns are telling different stories about the same well.
    per_category = d[[f"screen_{c}_phenotype" for c in sorted(set(S.CATEGORIES.values()))]]
    assert (per_category.any(axis=1) == d["screen_any_phenotype"]).all()


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
                    reason="node table not built")
def test_the_shipped_columns_did_not_swallow_the_unscreened_genes():
    """8,140 genes, 319 of them screened. The rest must stay missing rather than reading normal."""
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"),
                            columns=["gene_id", "screen_egress_phenotype", "screen_any_phenotype"])
    known = nodes["screen_egress_phenotype"].notna().sum()
    assert 300 < known < 400, f"{known} genes carry a screen call; the screen covered 319"
    assert (nodes["screen_any_phenotype"] == True).sum() < known  # noqa: E712
