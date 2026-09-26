"""Leakage: the audit finds what it must, and the grouping keeps what it fixed.

Planted tables first -- a monotone copy, a label that is a joint function of two columns, many
categories that explain nothing -- because an audit that cannot find a planted leak says nothing
about a real one. Then the shipped tables: the three leaks the 2026-09-25 audit found stay closed.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import leakage as L  # noqa: E402
from starplast import search  # noqa: E402


def _table(n=600, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    y = rng.normal(size=n)
    return pd.DataFrame({
        "gene_id": [f"TGME49_{i:06d}" for i in range(n)],
        "x": x,
        "x_exp": np.exp(3 * x),                           # a monotone copy of x
        "y": y,
        "label": np.where(x > y, "a", "b"),               # a joint function of x and y
        "noise": rng.normal(size=n),
        "many": rng.choice([f"c{i}" for i in range(30)], size=n),
    })


def test_a_monotone_copy_is_seen_on_ranks_and_missed_on_values():
    ranks, values = L.association_matrix(_table())
    assert ranks.loc["x", "x_exp"] > 0.99
    assert values.loc["x", "x_exp"] < 0.8          # what the closure used to measure


def test_many_categories_that_explain_nothing_score_near_zero():
    ranks, _v = L.association_matrix(_table())
    assert ranks.loc["many", "noise"] < 0.2


def test_the_closure_now_catches_a_monotone_copy():
    t = _table()
    assert "x_exp" in search.excluded_for(t, "x")


def test_calibration_reports_a_quantile_with_an_interval():
    rng = np.random.default_rng(1)
    pairs = pd.DataFrame({"association": rng.uniform(size=400),
                          "rank_association": rng.uniform(size=400),
                          "value_association": rng.uniform(size=400),
                          "relation": "different axis",
                          "slot_a": rng.choice(list("abcdef"), size=400),
                          "slot_b": rng.choice(list("ghijkl"), size=400)})
    c = L.calibrate(pairs, 0.95, n_boot=200)
    assert 0.85 < c["threshold"] < 1.0 and c["low"] <= c["threshold"] <= c["high"]
    assert np.isnan(L.calibrate(pairs.assign(relation="same slot"))["threshold"])


def test_a_gap_is_a_strong_association_the_closure_lets_through(monkeypatch):
    t = _table()
    pairs = pd.DataFrame({"a": ["x"], "b": ["y"], "association": [0.95],
                          "rank_association": [0.95], "value_association": [0.95]})
    gaps = L.closure_gaps(t, pairs, 0.8)
    # Both directions: holding x out lets y through, and holding y out lets x through.
    assert set(zip(gaps["target"], gaps["column"])) == {("x", "y"), ("y", "x")}
    monkeypatch.setattr(search, "excluded_for", lambda nodes, target, **k: {"x", "y"})
    assert L.closure_gaps(t, pairs, 0.8).empty


def test_a_joint_derivation_is_flagged_and_an_honest_predictor_is_not(monkeypatch):
    from starplast import slots

    class Slot:
        def __init__(self, key, axis):
            self.key, self.axis = key, axis

    t = _table(n=1200)
    t["weak"] = t["x"] + 3 * np.random.default_rng(2).normal(size=len(t))
    for i in range(12):                                    # a field of unrelated slots
        t[f"other{i}"] = np.random.default_rng(10 + i).normal(size=len(t))
    table = {"both": ["x", "y"], "weak": ["weak"]} | {f"o{i}": [f"other{i}"] for i in range(12)}
    fake = [Slot(k, "axis") for k in table]
    monkeypatch.setattr(slots, "all_slots", lambda organism=None: tuple(fake))
    monkeypatch.setattr(slots, "declared_columns",
                        lambda nodes, sl, numeric_only=False: tuple(table[sl.key]))
    monkeypatch.setattr(search, "excluded_for", lambda nodes, target, **k: {"label"})
    out = L.residual_leaks(t, "label", organism="Tg")
    flagged = set(out.loc[out["flag"], "slot"])
    assert flagged == {"both"}, out


def test_relations_name_what_the_catalogue_declares():
    assert L.relation("cellcycle_phase", "cellcycle19092_1h_r1", {}, "Tg") == "same quantity"
    assert L.relation("fit_invitro_hff", "crispr_gra17_synthlethal_delta", {}, "Tg") \
        == "same quantity"
    assert L.relation("a", "b", {}) == "undeclared"


# --------------------------------------------------------------------------- the shipped tables
@pytest.fixture(scope="module")
def tg():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("nodes.parquet"))


@pytest.fixture(scope="module")
def pf():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("pf_nodes.parquet"))


def test_holding_out_the_phase_holds_out_the_phase_resolved_expression(tg):
    banned = search.excluded_for(tg, "cellcycle_phase")
    series = [c for c in tg.columns if c.startswith("cellcycle19092_")]
    assert series and set(series) <= banned


def test_holding_out_fibroblast_fitness_holds_out_the_second_screen(tg):
    banned = search.excluded_for(tg, "fit_invitro_hff")
    assert {"crispr_gra17ko_phenotype", "crispr_gra17_synthlethal_delta"} <= banned
    assert "fit_invivo_PE" not in banned               # another condition, another quantity


def test_topology_takes_its_helix_count_with_it(tg):
    assert "n_tm" in search.excluded_for(tg, "dtm_class")


def test_a_small_many_category_overlap_bans_nothing(tg):
    """24 enzyme classes over 43 genes once banned two knockout-screen columns."""
    assert not any(c.startswith("crispr_") for c in search.excluded_for(tg, "ec_number"))


def test_the_plasmodium_stage_label_is_closed_over_its_own_sources(pf):
    banned = search.excluded_for(pf, "stage_enriched_derived")
    assert {"expr_ring", "expr_schizont", "expr_ookinete"} <= banned


def test_a_shared_column_name_resolves_to_the_right_organism():
    from starplast import datasets
    assert datasets.provenance("mean_plddt", "Pf").key == "pf_alphafold_confidence"
    assert datasets.provenance("mean_plddt", "Tg").key == "alphafold"
    assert "expr_ring" in datasets.derived_sources("stage_enriched_derived", "Pf")
    assert "expr_tachy" in datasets.derived_sources("stage_enriched_derived", "Tg")
    assert datasets.organism_of(datasets.provenance("mean_plddt", "Pf")) == "Pf"


@pytest.mark.slow
@pytest.mark.parametrize("target", ["cellcycle_phase", "fit_invitro_hff", "dtm_class",
                                    "stage_enriched_derived"])
def test_no_permitted_slot_predicts_a_held_out_target_like_a_copy(tg, target):
    out = L.residual_leaks(tg, target, organism="Tg")
    assert not out["flag"].any(), out[out["flag"]]
