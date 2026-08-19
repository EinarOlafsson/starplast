#!/usr/bin/env python3
"""Cell-cycle phase and stage enrichment.

Two labels that must never be confused, so most of what is asserted here is about which is which. One
is a measurement another laboratory published; the other is a restatement of columns already in the
node table. Treating the second as evidence would produce a spectacular, meaningless result -- this
project has already done that once, when a leaked `compartment` column came back as its own top
discovery at V = 0.96 against a truthful mean F1 of 0.362.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import cellcycle as CC  # noqa: E402


# --------------------------------------------------------------------------- reading the study files
@pytest.fixture
def phase_file(tmp_path, monkeypatch):
    """The real file's shape: tab-separated under a .csv name, unnamed first column, quoted phases."""
    from starplast import paths
    d = tmp_path / CC.STUDY
    d.mkdir(parents=True)
    (d / CC.PHASE_FILE).write_text(
        "\tfold\tmax_cluster\tproduct\tcell_cycle_phase\n"
        'TGGT1_300608\t2.0\t"""C"""\tribosomal protein\tC\n'
        'TGGT1_411260\t1.5\t"""M"""\thypothetical\tM\n'
        'TGGT1_408520\t1.2\t"""S"""\ttRNA-Asp\tS\n'
        'TGGT1_410930\t1.1\t"""G1 a"""\thypothetical\tG1a\n'
    )
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    return d


def test_the_file_is_tab_separated_despite_its_extension(phase_file):
    """Read as comma-separated it becomes one column of strings and every join silently yields nothing
    -- which looks like a dataset with no coverage, not like a parse error."""
    t = CC.phase_labels("", log=lambda *_: None)
    assert len(t) == 4
    assert set(t.cellcycle_phase) == {"C", "M", "S", "G1a"}


def test_the_triple_quoting_in_the_export_is_stripped(phase_file):
    t = CC.phase_labels("", log=lambda *_: None)
    assert all('"' not in p for p in t.cellcycle_phase)


def test_only_the_five_real_phases_are_kept(tmp_path, monkeypatch):
    from starplast import paths
    d = tmp_path / CC.STUDY
    d.mkdir(parents=True)
    (d / CC.PHASE_FILE).write_text(
        "\tcell_cycle_phase\nTGGT1_300608\tC\nTGGT1_411260\tnonsense\nTGGT1_408520\t\n")
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = CC.phase_labels("", log=lambda *_: None)
    assert list(t.cellcycle_phase) == ["C"]


def test_a_missing_study_file_reports_and_returns_empty(tmp_path, monkeypatch):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    # Patch the resolver itself, not just the env root. `paths.find` walks EVERY root it can name, and
    # one of them is derived from the package location rather than from the environment -- so an empty
    # `STARPLAST_DATA` never made the file missing. This test claimed a missing study file while a
    # real dataset tree was still reachable, and it only started failing once that tree gained the
    # file. Naming one root is what makes "absent" mean absent.
    monkeypatch.setattr(paths, "dataset_roots", lambda: [str(tmp_path)])
    msgs = []
    assert CC.phase_labels("", log=msgs.append).empty
    assert any("not found" in m for m in msgs)
    assert CC.pseudotime_clusters("", log=msgs.append).empty


# --------------------------------------------------------------------------- the strain trap
def test_without_identity_resolution_tggt1_accessions_join_nothing(phase_file):
    """The documented failure, reproduced. These files use TGGT1_ (type I RH) while the project is
    built on ME49, and an unresolved join produces an empty result rather than an error."""
    t = CC.phase_labels("", resolve=None, log=lambda *_: None)
    assert all(g.startswith("TGGT1_") for g in t.index)
    me49 = pd.Series(["TGME49_300608", "TGME49_411260"])
    assert me49.isin(t.index).sum() == 0


def test_with_identity_resolution_they_land_on_me49(phase_file):
    t = CC.phase_labels("", resolve=lambda a: a.replace("TGGT1_", "TGME49_"), log=lambda *_: None)
    assert all(g.startswith("TGME49_") for g in t.index)


def test_unresolvable_accessions_are_reported_and_dropped_not_kept_as_nan(phase_file):
    msgs = []
    t = CC.phase_labels("", resolve=lambda a: None if a.endswith("608") else a.replace("GT1", "ME49"),
                        log=msgs.append)
    assert len(t) == 3
    assert any("did not resolve" in m for m in msgs)
    assert not t.index.isna().any()


def test_a_gene_with_two_conflicting_phases_is_left_unlabelled(tmp_path, monkeypatch):
    """Two probes disagreeing is a real ambiguity. Keeping whichever sorted first would record a
    confident phase that the data does not support."""
    from starplast import paths
    d = tmp_path / CC.STUDY
    d.mkdir(parents=True)
    (d / CC.PHASE_FILE).write_text(
        "\tcell_cycle_phase\nTGGT1_1\tS\nTGGT1_2\tS\nTGGT1_2\tM\n")
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    msgs = []
    t = CC.phase_labels("", resolve=lambda a: a, log=msgs.append)
    assert list(t.index) == ["TGGT1_1"]
    assert any("conflicting" in m for m in msgs)


def test_a_gene_listed_twice_with_the_same_phase_is_kept_once(tmp_path, monkeypatch):
    from starplast import paths
    d = tmp_path / CC.STUDY
    d.mkdir(parents=True)
    (d / CC.PHASE_FILE).write_text("\tcell_cycle_phase\nTGGT1_1\tS\nTGGT1_1\tS\n")
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = CC.phase_labels("", resolve=lambda a: a, log=lambda *_: None)
    assert list(t.index) == ["TGGT1_1"]


# --------------------------------------------------------------------------- pseudotime
def test_pseudotime_clusters_are_numeric_and_resolved(tmp_path, monkeypatch):
    from starplast import paths
    d = tmp_path / CC.STUDY
    d.mkdir(parents=True)
    (d / CC.PSEUDOTIME_FILE).write_text(
        "\tpseudotime_cluster\tgene_product\nTGGT1_1\t1\thypothetical\n"
        "TGGT1_2\t3\thypothetical\nTGGT1_3\tnot-a-number\thypothetical\n")
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    t = CC.pseudotime_clusters("", resolve=lambda a: a.replace("TGGT1", "TGME49"),
                               log=lambda *_: None)
    assert len(t) == 2
    assert t.cellcycle_pseudotime.tolist() == [1.0, 3.0]
    assert all(g.startswith("TGME49") for g in t.index)


# --------------------------------------------------------------------------- derived stage classes
def _nodes(tachy, cyst, sporo):
    return pd.DataFrame({"gene_id": [f"TGME49_{i:06d}" for i in range(len(tachy))],
                         "expr_tachy": tachy, "expr_cyst": cyst, "expr_sporulated": sporo})


def test_the_derived_column_announces_itself_in_its_name():
    """A reader must not be able to mistake it for a measurement, in any table it appears in."""
    n = _nodes([3.0, 0.0], [0.0, 3.0], [0.0, 0.0])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert all(c.endswith("_derived") for c in out.columns)


def test_the_stage_with_the_highest_expression_wins():
    n = _nodes([5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert list(out.stage_enriched_derived) == ["tachyzoite", "bradyzoite", "oocyst"]


def test_a_gene_with_no_clear_winner_is_left_unlabelled():
    """A label that is really a coin toss is worse than an absent one: it looks like a measurement."""
    n = _nodes([1.0, 1.01, 0.99], [1.0, 1.0, 1.0], [1.0, 0.99, 1.01])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert out.stage_enriched_derived.isna().all()


def test_columns_are_z_scored_so_a_wider_column_cannot_win_on_units_alone():
    """expr_sporulated has a third the spread of expr_tachy. Compared raw, whichever column happens to
    carry the larger numbers wins nearly every gene, which is a fact about units, not biology."""
    n = pd.DataFrame({"gene_id": ["a", "b", "c", "d"],
                      "expr_tachy": [100.0, 200.0, 300.0, 400.0],      # large values, wide
                      "expr_cyst": [1.0, 1.1, 1.2, 1.3],               # small values, narrow
                      "expr_sporulated": [0.5, 0.5, 0.5, 0.5]})
    out = CC.stage_enrichment(n, log=lambda *_: None)
    called = out.stage_enriched_derived.dropna()
    assert set(called) != {"tachyzoite"}, "the largest-unit column must not sweep every gene"


def test_a_constant_column_does_not_divide_by_zero():
    n = _nodes([1.0, 1.0, 1.0], [1.0, 2.0, 3.0], [0.0, 1.0, 2.0])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert np.isfinite(out.stage_margin_derived.dropna()).all()


def test_a_gene_measured_in_only_one_stage_is_not_called():
    """One stage present is not a comparison, and idxmax over a single value would still return it."""
    n = pd.DataFrame({"gene_id": ["a"], "expr_tachy": [5.0],
                      "expr_cyst": [np.nan], "expr_sporulated": [np.nan]})
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert out.stage_enriched_derived.isna().all()


def test_fewer_than_two_stages_present_returns_nothing_and_says_so():
    msgs = []
    out = CC.stage_enrichment(pd.DataFrame({"gene_id": ["a"], "expr_tachy": [1.0]}), log=msgs.append)
    assert out.empty or out.columns.empty
    assert any("not enough" in m for m in msgs)


def test_the_margin_is_recorded_so_a_weak_call_can_be_told_from_a_strong_one():
    n = _nodes([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    called = out.dropna(subset=["stage_enriched_derived"])
    assert (called.stage_margin_derived >= CC.MIN_MARGIN).all()


def test_the_margin_is_absent_wherever_no_call_was_made():
    n = _nodes([1.0, 1.0], [1.0, 1.0], [1.0, 1.0])
    out = CC.stage_enrichment(n, log=lambda *_: None)
    assert out.stage_margin_derived.isna().all()


# --------------------------------------------------------------------------- composition
def test_add_all_attaches_both_kinds_of_column(phase_file, tmp_path):
    n = _nodes([3.0, 0.0], [0.0, 3.0], [0.0, 0.0])
    n["gene_id"] = ["TGME49_300608", "TGME49_411260"]
    out = CC.add_all("", n, resolve=lambda a: a.replace("TGGT1_", "TGME49_"), log=lambda *_: None)
    assert "cellcycle_phase" in out.columns
    assert "stage_enriched_derived" in out.columns
    assert out.loc[0, "cellcycle_phase"] == "C"


def test_the_registry_records_the_derived_column_as_derived():
    """The name is the cheap part; this is what a reader checks.

    Asserted on `derived_from` rather than on `kind`. This test used to require kind == "derived",
    which conflates two axes: kind answers what sort of measurement is underneath, and here that is
    bulk RNA-seq -- the argmax of three expression columns. Pinning the derivation to kind left the
    one entry in the registry whose data type could not be read off its type field. The circularity
    guard reads derived_from, so that is also the field that has to be right for anything to work.
    """
    from starplast import datasets
    d = datasets.get("stage_enriched")
    assert "DERIVED" in d.note
    assert "circular" in d.note.lower()
    assert "DERIVED" in d.name
    assert d.derived_from == ("expr_tachy", "expr_cyst", "expr_sporulated")
    assert d.kind == "RNAseq", "kind should name the assay underneath, not how it was produced"


def test_the_registry_records_the_measured_study_with_its_citation():
    from starplast import datasets
    d = datasets.get("xue_singlecell")
    assert d.pmid == "32065584"
    assert d.citation and "eLife" in d.citation
    assert "cellcycle_phase" in d.columns


# --------------------------------------------------------------------------- phase on a closed loop
def test_a_cycle_has_no_ends_so_the_peak_is_a_phase_and_not_a_column():
    """The construction the Plasmodium time course needed, kept here because it is not about either
    parasite: it is what a time course around a closed loop means.

    Hour 48 is hour 0 of the next cycle, so an argmax splits a peak that straddles the wrap. These
    two genes peak four hours apart across that boundary, and any reading that puts them 44 hours
    apart is reading a circle as a line.
    """
    import numpy as np
    hours = np.arange(3, 49, 3)
    frame = pd.DataFrame(
        {float(h): [10 + 5 * np.cos(2 * np.pi * (h - 46) / 48),
                    10 + 5 * np.cos(2 * np.pi * (h - 2) / 48),
                    7.0] for h in hours},
        index=["just_before", "just_after", "flat"])
    out = CC.cyclic_phase(frame, period=48.0)
    gap = abs(out.loc["just_before", "phase"] - out.loc["just_after", "phase"])
    assert min(gap, 48 - gap) < 6, "the wrap was read as a 44-hour separation"
    assert out.loc["flat", "amplitude"] < 0.01, "a flat profile still has an angle; it must not rank"
    assert out.loc["just_before", "amplitude"] > 0.9


def test_the_phase_of_a_pure_cosine_is_where_it_peaks():
    import numpy as np
    hours = np.arange(0, 48, 4)
    frame = pd.DataFrame({float(h): [np.cos(2 * np.pi * (h - 20) / 48)] for h in hours},
                         index=["peaks_at_20"])
    assert abs(CC.cyclic_phase(frame, period=48.0).loc["peaks_at_20", "phase"] - 20) < 1
