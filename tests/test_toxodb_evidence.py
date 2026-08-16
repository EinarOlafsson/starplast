#!/usr/bin/env python3
"""Evidence ToxoDB integrates, and the check that decided one column against two.

Three reports were fetched; two are here and one is in quarantine. The H3K4me1 ChIP-on-chip looked
exactly as usable as the other two until its "marked" genes turned out to have LESS accessible
promoters than unmarked ones. Nothing in its metadata says that. Only computing the column and
testing it against another column in the map does.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import toxodb_evidence as TE  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _report(tmp_path, filename, rows, header="Gene ID\tValue"):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    (data / filename).write_text(header + "\n" + "\n".join(f"{g}\t{v}" for g, v in rows) + "\n")
    return tmp_path


def test_a_plain_count_is_read_as_it_stands(tmp_path):
    _report(tmp_path, "toxodb_epitopes.tsv", [("TGME49_200010", 45)])
    out = TE.evidence(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "iedb_epitope_count"] == 45


def test_a_fold_report_is_read_with_its_sign_convention(tmp_path):
    """-1257.4 means 1257-fold down. Logged directly it is NaN, and every repressed gene vanishes."""
    _report(tmp_path, "toxodb_enteroepithelial.tsv",
            [("TGME49_200010", -4.0), ("TGME49_200020", 8.0)])
    out = TE.evidence(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "ees_vs_tachyzoite_log2"] == pytest.approx(-2.0)
    assert out.loc["TGME49_200020", "ees_vs_tachyzoite_log2"] == pytest.approx(3.0)


def test_the_value_is_taken_by_position_not_by_header(tmp_path):
    """ToxoDB ships display names as the header and they move with the site release."""
    _report(tmp_path, "toxodb_epitopes.tsv", [("TGME49_200010", 3)],
            header="Gene ID\tSome Renamed Column")
    out = TE.evidence(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "iedb_epitope_count"] == 3


def test_accessions_go_through_the_identity_layer(tmp_path):
    _report(tmp_path, "toxodb_epitopes.tsv", [("TGGT1_100010", 3)])
    out = TE.evidence(str(tmp_path), resolve=lambda g: {"TGGT1_100010": "TGME49_200010"}.get(g),
                      log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _report(tmp_path, "toxodb_epitopes.tsv", [("TGME49_999999", 3)])
    out = TE.evidence(str(tmp_path), resolve=lambda g: None, log=lambda *_: None)
    assert list(out.index) == ["TGME49_999999"]


def test_two_rows_for_one_gene_keep_the_stronger(tmp_path):
    _report(tmp_path, "toxodb_epitopes.tsv", [("TGME49_200010", 2), ("TGME49_200010", 9)])
    out = TE.evidence(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "iedb_epitope_count"] == 9


def test_reports_are_joined_across_sources(tmp_path):
    _report(tmp_path, "toxodb_epitopes.tsv", [("A", 1)])
    _report(tmp_path, "toxodb_enteroepithelial.tsv", [("B", 2.0)])
    out = TE.evidence(str(tmp_path), log=lambda *_: None)
    assert set(out.index) == {"A", "B"}
    assert sorted(out.columns) == ["ees_vs_tachyzoite_log2", "iedb_epitope_count"]


def test_a_one_column_report_is_refused(tmp_path):
    """The searches can be asked for gene ids alone; that report carries no measurement."""
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True)
    (data / "toxodb_epitopes.tsv").write_text("Gene ID\nTGME49_200010\n")
    assert TE.evidence(str(tmp_path), log=lambda *_: None).empty


def test_nothing_downloaded_yields_nothing(tmp_path):
    assert TE.evidence(str(tmp_path), log=lambda *_: None).empty


def test_the_rejected_vesicle_report_is_not_a_source():
    """Secreted GRA and MIC proteins come out DEPLETED from the vesicle fraction and ribosomal
    proteins enriched. Whether that is a correct measurement of something else does not matter --
    it is not an answer to what the parasite secretes."""
    assert not any("vesicle" in filename for filename, *_ in TE.SOURCES)
    assert not os.path.exists(os.path.join(ROOT, "starplast", "data", "toxodb_vesicles.tsv"))


def test_the_rejected_histone_report_is_not_a_source():
    """It stays out until somebody explains the backwards correlation, not until it looks tidy."""
    assert not any("h3k4me1" in filename for filename, *_ in TE.SOURCES)
    assert not os.path.exists(os.path.join(ROOT, "starplast", "data", "toxodb_h3k4me1.tsv"))


def test_every_source_records_what_was_asked_for():
    """A search that can be asked five ways produces five different columns."""
    for filename, column, _fold, query in TE.SOURCES:
        assert query and len(query) > 20, column
        assert "Genes" in query, f"{column} does not name the search it came from"


# --------------------------------------------------------------------------- against the real map
@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_the_immunodominant_antigens_carry_the_most_epitopes():
    """SAG1 -- which ToxoDB calls SRS29B -- GRA6, GRA7 and GRA2 are what Toxoplasma serology uses.

    If the epitope join were wrong this ordering is the first thing to go, and a count column would
    otherwise look perfectly reasonable.
    """
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    if "iedb_epitope_count" not in n.columns:
        pytest.skip("epitope column not merged")
    top = n["iedb_epitope_count"].dropna().sort_values(ascending=False).head(5).index
    products = " ".join(n.loc[top, "product"].astype(str)).upper()
    assert "SRS29B" in products, products
    assert sum(p in products for p in ("GRA6", "GRA7", "GRA2")) >= 2, products


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_merozoite_markers_rise_in_the_enteroepithelial_stages():
    """GRA11B is merozoite-specific, and the EES are where merozoites are. It must be strongly up."""
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    column = "ees_vs_tachyzoite_log2"
    if column not in n.columns:
        pytest.skip("enteroepithelial column not merged")
    hit = n["product"].astype(str).str.contains("GRA11", case=False, na=False)
    values = n.loc[hit, column].dropna()
    assert len(values), "GRA11B not measured"
    assert values.median() > 4.0, f"GRA11B at {values.median():+.2f}"
    assert abs(n[column].median()) < 0.5, "the genome as a whole should not have shifted"


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_the_active_histone_mark_behaves_like_one():
    """H4 acetylation must agree with the two independent measures of an active gene already here.

    This is the test the H3K4me1 report failed, and running it against the column that replaced it
    is what keeps that refusal honest: if this one ever drifts the same way, it goes too.
    """
    from scipy.stats import spearmanr
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet"))
    if "h4_acetylation_chip_score" not in n.columns:
        pytest.skip("H4 acetylation column not merged")
    for other, floor in (("expr_tachy", 0.25), ("atac_promoter_ut", 0.30)):
        k = n[["h4_acetylation_chip_score", other]].dropna()
        rho = spearmanr(k["h4_acetylation_chip_score"], k[other]).statistic
        assert rho > floor, f"{other}: rho {rho:+.3f}"


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_novel_isoforms_land_on_the_genes_with_more_exons():
    """More exons, more ways to splice. If the join were wrong this is what would vanish."""
    from scipy.stats import mannwhitneyu
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet"))
    if "novel_transcript_models" not in n.columns or "n_exons" not in n.columns:
        pytest.skip("isoform or exon column not present")
    have = n.loc[n["novel_transcript_models"].notna(), "n_exons"].dropna()
    none = n.loc[n["novel_transcript_models"].isna(), "n_exons"].dropna()
    assert have.median() > none.median()
    assert mannwhitneyu(have, none).pvalue < 1e-10
