"""Host source versions cannot manufacture equivalence or assign ambiguous groups."""
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from review_rbc_candidate import group_candidates, parse_copy_numbers
from review_rhoptry_version import compare_tables


def _table():
    return pd.DataFrame({"Gene": ["A", "B"], "Combined Rhoptry Score": [1.0, 2.0],
                         "greenN|beta": [0.1, -0.2], "greenN|wald-fdr": [0.0, None]})


def test_exact_version_comparison_keeps_order_independence_and_missingness():
    before = _table()
    assert compare_tables(before, before.iloc[::-1])["assay_values_identical"]
    after = before.copy()
    after.loc[1, "greenN|wald-fdr"] = 0.0
    report = compare_tables(before, after)
    assert not report["assay_values_identical"]
    assert report["fields"]["greenN|wald-fdr"]["different_cells"] == 1
    assert not compare_tables(before, before.iloc[:1])["assay_values_identical"]


def test_version_comparison_rejects_duplicate_genes_and_invalid_numbers():
    with pytest.raises(ValueError, match="unique"):
        compare_tables(_table(), pd.concat([_table(), _table()]))
    invalid = _table().astype({"greenN|beta": object})
    invalid.loc[0, "greenN|beta"] = "not measured"
    with pytest.raises(ValueError):
        compare_tables(_table(), invalid)


def test_group_projection_preserves_all_reviewed_choices_and_excludes_contaminants():
    reviewed = {"P00001", "P00002"}
    assert group_candidates("P00001;P00001-2;unreviewed", reviewed) == ("P00001",)
    assert group_candidates("P00001;P00002;CON__P00001;REV__P00002", reviewed) == ("P00001", "P00002")
    assert group_candidates("CON__P00001;unreviewed", reviewed) == ()


def test_copy_count_parser_refuses_a_plausible_table_with_swapped_fraction_headers(tmp_path):
    raw = pd.DataFrame(None, index=range(4), columns=range(39))
    raw.iloc[1, 0] = "Protein IDs"
    raw.iloc[1, 31] = "White ghosts"
    raw.iloc[1, 35] = "Whole erythrocyte"
    file = tmp_path / "candidate.xlsx"
    raw.to_excel(file, sheet_name="Table S3", header=False, index=False)
    with pytest.raises(ValueError):
        parse_copy_numbers(file)
