#!/usr/bin/env python3
"""The metabolite table, and the two things that made it awkward.

The archive's two labelling sheets differ in whether they carry a title row; assuming they did not
silently dropped the glucose arm, and the failure was invisible because the sheet read fine and
simply had no column called `Metabolite`. And a compound has no stable identifier the way a gene
does, so the row key is a normalised name and joining a second study will be lossy.
"""
from __future__ import annotations

import io
import os
import sys
import zipfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import metabolites as M  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _archive(tmp_path, *, title_row_on_labelling=True, levels=True, labelling=True):
    """A supplementary archive shaped like the real one."""
    path = tmp_path / "supp.zip"
    with zipfile.ZipFile(path, "w") as z:
        if levels:
            body = pd.DataFrame({
                "Unnamed: 0": ["Label", "Glycine", "Serine"],
                "Normalised abundance": ["Ctrl_02", 0.1, 0.2],
                "Volcano raw data": ["Fold Change", 0.79, 1.2],
                "Unnamed: 9": ["log2(FC) DFO/Ctrl", -0.33, 0.26],
                "Unnamed: 10": ["p.ajusted", 0.045, 0.5]})
            buf = io.BytesIO()
            with pd.ExcelWriter(buf) as w:
                pd.DataFrame([["Table S3: summary"]]).to_excel(
                    w, sheet_name=M.LEVELS[1], index=False, header=False)
                body.to_excel(w, sheet_name=M.LEVELS[1], index=False, startrow=1)
            z.writestr(M.LEVELS[0], buf.getvalue())
        if labelling:
            frame = pd.DataFrame({"Metabolite": ["Glycine", "Glycine", "Serine"],
                                  "ExpCondition": ["Control", "Control", "Control"],
                                  "M0": [0.4, 0.6, 1.0], "M1": [0.6, 0.4, 0.0]})
            buf = io.BytesIO()
            with pd.ExcelWriter(buf) as w:
                for precursor, (_member, sheet) in M.LABELLING.items():
                    if title_row_on_labelling:
                        pd.DataFrame([["Table S5: labelling"]]).to_excel(
                            w, sheet_name=sheet, index=False, header=False)
                        frame.to_excel(w, sheet_name=sheet, index=False, startrow=1)
                    else:
                        frame.to_excel(w, sheet_name=sheet, index=False)
            z.writestr(next(iter(M.LABELLING.values()))[0], buf.getvalue())
    return str(path)


# --------------------------------------------------------------------------- the row key
def test_a_compound_name_normalises_to_something_two_studies_could_share():
    assert M.norm("L-Glutamic acid") == M.norm("l glutamic  acid")
    assert M.norm("1-18:0-2-18:1-phosphatidylinositol") == "118021 81phosphatidylinositol".replace(" ", "")


def test_an_unusable_name_has_no_key():
    assert M.norm(None) is None and M.norm("") is None and M.norm("   ") is None


# --------------------------------------------------------------------------- building
def test_levels_and_both_labelling_arms_are_read(tmp_path):
    d = M.build(_archive(tmp_path))
    assert "metabolite_level_log2fc_iron_depleted" in d.columns
    assert "labelled_fraction_glucose" in d.columns
    assert "labelled_fraction_glutamine" in d.columns


def test_a_labelling_sheet_without_a_title_row_is_read_too(tmp_path):
    """The two sheets in the real archive differ, and assuming otherwise dropped the glucose arm."""
    d = M.build(_archive(tmp_path, title_row_on_labelling=False))
    assert d["labelled_fraction_glucose"].notna().any()


def test_labelled_fraction_is_one_minus_the_unlabelled_isotopologue(tmp_path):
    """M1..M55 are not comparable between molecules of different carbon number; M0 is."""
    d = M.build(_archive(tmp_path)).set_index("metabolite")
    assert d.loc["Glycine", "labelled_fraction_glucose"] == pytest.approx(0.5)
    assert d.loc["Serine", "labelled_fraction_glucose"] == pytest.approx(0.0)


def test_an_archive_with_only_levels_still_builds(tmp_path):
    d = M.build(_archive(tmp_path, labelling=False))
    assert len(d) and "labelled_fraction_glucose" not in d.columns


def test_an_archive_with_only_labelling_still_builds(tmp_path):
    d = M.build(_archive(tmp_path, levels=False))
    assert len(d) and "metabolite_level_log2fc_iron_depleted" not in d.columns


def test_a_missing_archive_yields_nothing(tmp_path):
    assert M.build(str(tmp_path / "absent.zip")).empty


def test_an_empty_archive_yields_nothing(tmp_path):
    p = tmp_path / "empty.zip"
    zipfile.ZipFile(p, "w").close()
    assert M.build(str(p)).empty


def test_load_returns_nothing_when_the_table_is_not_built(tmp_path):
    assert M.load(str(tmp_path)).empty


# --------------------------------------------------------------------------- the shipped table
@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", M.TABLE)),
                    reason="metabolite table not built")
def test_the_shipped_table_has_rows_that_are_compounds():
    d = M.load(ROOT)
    assert len(d) > 500
    assert not d["metabolite"].astype(str).str.contains("TGME49_").any(), (
        "gene accessions have got into the metabolite table")
    labelled = d["labelled_fraction_glucose"].dropna()
    assert ((labelled >= -0.01) & (labelled <= 1.01)).all(), "a labelled fraction outside 0..1"


def test_a_levels_sheet_with_no_fold_change_column_is_refused(tmp_path):
    """The sheet carries two blocks side by side; without the volcano block there is no summary."""
    p = tmp_path / "supp.zip"
    body = pd.DataFrame({"Unnamed: 0": ["Label", "Glycine"], "Normalised abundance": ["Ctrl", 0.1]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame([["Table S3"]]).to_excel(w, sheet_name=M.LEVELS[1], index=False, header=False)
        body.to_excel(w, sheet_name=M.LEVELS[1], index=False, startrow=1)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(M.LEVELS[0], buf.getvalue())
    assert M.build(str(p)).empty


def test_a_labelling_sheet_with_no_metabolite_column_is_refused(tmp_path):
    p = tmp_path / "supp.zip"
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        for _precursor, (_member, sheet) in M.LABELLING.items():
            pd.DataFrame({"something": [1], "else": [2]}).to_excel(w, sheet_name=sheet, index=False)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(next(iter(M.LABELLING.values()))[0], buf.getvalue())
    assert M.build(str(p)).empty


def test_a_labelling_sheet_with_no_control_samples_is_refused(tmp_path):
    """The labelled fraction is defined in control conditions; a treated-only sheet cannot give it."""
    p = tmp_path / "supp.zip"
    frame = pd.DataFrame({"Metabolite": ["Glycine"], "ExpCondition": ["DFO"], "M0": [0.4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        for _precursor, (_member, sheet) in M.LABELLING.items():
            frame.to_excel(w, sheet_name=sheet, index=False)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(next(iter(M.LABELLING.values()))[0], buf.getvalue())
    assert M.build(str(p)).empty


def test_a_row_seen_in_both_blocks_appears_once(tmp_path):
    d = M.build(_archive(tmp_path))
    assert d["metabolite"].str.lower().duplicated().sum() == 0
