#!/usr/bin/env python3
"""hyperLOPIT, and the cross-species transfer that is gated on its own measured accuracy.

The interesting property of this module is that it refuses to trust itself. Transfer from Plasmodium
and Cryptosporidium is calibrated against the 826 genes carrying both a measured and a transferable
label, and the result is used to gate: proteasome transfers at 100% and golgi at 6.1%, so one is kept
and the other is not. An average would have reported 71.8% and hidden both facts.

The calibration is not circular, and the tests below pin the reason: transfers are only APPLIED where a
measured label is absent, so the set they are scored on and the set they are used on are disjoint by
construction.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import localization as L  # noqa: E402


# --------------------------------------------------------------------------- fixtures
def _vocab_file(d, rows):
    pd.DataFrame(rows, columns=["species", "native_term", "unified_category"]).to_csv(
        d / "lopit_vocabulary_dictionary.csv", index=False)


def _species_file(d, species, ids, labels):
    fn, idc, labc = L.SPECIES[species]
    pd.DataFrame({idc: ids, labc: labels}).to_csv(d / fn, index=False)


def _master(d, rows):
    pd.DataFrame(rows).to_csv(d / "MASTER_parasite_wide_by_orthogroup.csv", index=False)


def _nodes(ids):
    return pd.DataFrame({"gene_id": list(ids)})


# --------------------------------------------------------------------------- vocabulary
def test_a_missing_vocabulary_is_empty_not_an_error(tmp_path):
    assert L._vocabulary(str(tmp_path)) == {}


def test_the_vocabulary_maps_species_and_native_term_to_one_category(tmp_path):
    _vocab_file(tmp_path, [["tgon", "Apicoplast", "apicoplast"],
                           ["pfal", "apicoplast", "apicoplast"]])
    v = L._vocabulary(str(tmp_path))
    assert v[("tgon", "apicoplast")] == "apicoplast"
    assert v[("pfal", "apicoplast")] == "apicoplast"


def test_native_terms_are_matched_case_insensitively(tmp_path):
    """The three species' tables capitalise their compartment names differently."""
    _vocab_file(tmp_path, [["tgon", "Apicoplast", "apicoplast"]])
    _species_file(tmp_path, "tgon", ["g1"], ["APICOPLAST"])
    out = L._measured(str(tmp_path), "tgon", L._vocabulary(str(tmp_path)))
    assert out["g1"] == "apicoplast"


# --------------------------------------------------------------------------- measured labels
def test_a_missing_species_file_yields_nothing(tmp_path):
    assert L._measured(str(tmp_path), "pfal", {}).empty


def test_a_species_file_without_the_expected_columns_yields_nothing(tmp_path):
    fn, _, _ = L.SPECIES["pfal"]
    pd.DataFrame({"something": [1]}).to_csv(tmp_path / fn, index=False)
    assert L._measured(str(tmp_path), "pfal", {}).empty


def test_unknown_and_blank_labels_are_dropped_rather_than_becoming_categories(tmp_path):
    _vocab_file(tmp_path, [["pfal", "cytosol", "cytosol"]])
    _species_file(tmp_path, "pfal", ["a", "b", "c", "d"], ["cytosol", "unknown", "", np.nan])
    out = L._measured(str(tmp_path), "pfal", L._vocabulary(str(tmp_path)))
    assert list(out.index) == ["a"]


def test_a_term_absent_from_the_vocabulary_is_dropped(tmp_path):
    """An unmapped native term must not silently become its own category."""
    _vocab_file(tmp_path, [["pfal", "cytosol", "cytosol"]])
    _species_file(tmp_path, "pfal", ["a", "b"], ["cytosol", "something novel"])
    out = L._measured(str(tmp_path), "pfal", L._vocabulary(str(tmp_path)))
    assert list(out.index) == ["a"]


# --------------------------------------------------------------------------- the main path
def test_with_no_files_at_all_every_gene_is_unassigned(tmp_path):
    out = L.lopit_labels(str(tmp_path), _nodes(["g1", "g2"]), log=lambda *_: None)
    assert (out.compartment == "unassigned").all()
    assert (out.compartment_source == "unknown").all()


def test_a_missing_call_is_unassigned_not_a_twenty_seventh_compartment(tmp_path):
    """hyperLOPIT assignment tracks abundance, so treating "no call" as a location would create a
    compartment whose members are simply the low-abundance proteins."""
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": ["g1"], "MAP_location": ["apicoplast"],
                  "MCMC_location": ["apicoplast"]}).to_csv(tmp_path / fn, index=False)
    out = L.lopit_labels(str(tmp_path), _nodes(["g1", "g2"]), log=lambda *_: None)
    assert out.set_index("gene_id").loc["g1", "compartment"] == "apicoplast"
    assert out.set_index("gene_id").loc["g2", "compartment"] == "unassigned"


def test_map_and_mcmc_disagreement_is_recorded_as_uncertainty(tmp_path):
    """Where the two inferences differ, that is real inference uncertainty about the same experiment,
    not noise to be averaged away."""
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": ["g1", "g2"],
                  "MAP_location": ["apicoplast", "cytosol"],
                  "MCMC_location": ["apicoplast", "nucleus"]}).to_csv(tmp_path / fn, index=False)
    out = L.lopit_labels(str(tmp_path), _nodes(["g1", "g2"]), log=lambda *_: None).set_index("gene_id")
    assert out.loc["g1", "lopit_methods_agree"] == 1.0
    assert out.loc["g2", "lopit_methods_agree"] == 0.0


def test_agreement_is_missing_rather_than_false_when_only_one_method_called_it(tmp_path):
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": ["g1"], "MAP_location": ["apicoplast"],
                  "MCMC_location": [np.nan]}).to_csv(tmp_path / fn, index=False)
    out = L.lopit_labels(str(tmp_path), _nodes(["g1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "lopit_methods_agree"])


# --------------------------------------------------------------------------- transfer
def _transfer_setup(tmp_path, pf="apicoplast", cp="apicoplast", tgon_measured=None):
    _vocab_file(tmp_path, [["tgon", "apicoplast", "apicoplast"],
                           ["tgon", "cytosol", "cytosol"],
                           ["pfal", "apicoplast", "apicoplast"],
                           ["pfal", "cytosol", "cytosol"],
                           ["cpar", "apicoplast", "apicoplast"],
                           ["cpar", "nucleus", "nucleus"]])
    _species_file(tmp_path, "pfal", ["pf1"], [pf])
    _species_file(tmp_path, "cpar", ["cp1"], [cp])
    if tgon_measured is not None:
        fn, _, _ = L.SPECIES["tgon"]
        pd.DataFrame({"gene_source_id": ["t1"], "MAP_location": [tgon_measured],
                      "MCMC_location": [tgon_measured]}).to_csv(tmp_path / fn, index=False)
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": "pf1", "cpar_gene_id": "cp1"}])


def test_two_donors_agreeing_transfers_a_label(tmp_path):
    _transfer_setup(tmp_path)
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert out.loc[0, "ortholopit_label"] == "apicoplast"
    assert out.loc[0, "ortholopit_donors"] == "cpar+pfal"


def test_donors_that_disagree_transfer_nothing(tmp_path):
    """A conflict is a reason to say nothing, not a reason to pick one."""
    _transfer_setup(tmp_path, pf="cytosol", cp="nucleus")
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "ortholopit_label"])


def test_a_single_donor_can_transfer_but_is_recorded_as_single(tmp_path):
    _vocab_file(tmp_path, [["pfal", "apicoplast", "apicoplast"]])
    _species_file(tmp_path, "pfal", ["pf1"], ["apicoplast"])
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": "pf1", "cpar_gene_id": np.nan}])
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert out.loc[0, "ortholopit_label"] == "apicoplast"
    assert out.loc[0, "ortholopit_donors"] == "pfal"


def test_a_gene_with_no_orthologs_gets_no_transfer(tmp_path):
    _vocab_file(tmp_path, [["pfal", "apicoplast", "apicoplast"]])
    _species_file(tmp_path, "pfal", ["pf1"], ["apicoplast"])
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": np.nan, "cpar_gene_id": np.nan}])
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "ortholopit_label"])


def test_the_unmapped_marker_never_transfers(tmp_path):
    """A lineage-specific compartment that maps to nothing is not a label to hand to another species."""
    _vocab_file(tmp_path, [["pfal", "odd", L.UNMAPPED]])
    _species_file(tmp_path, "pfal", ["pf1"], ["odd"])
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": "pf1", "cpar_gene_id": np.nan}])
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "ortholopit_label"])


def test_transfer_needs_the_vocabulary(tmp_path):
    """Without it the donor labels are in three different naming systems and cannot be compared."""
    _species_file(tmp_path, "pfal", ["pf1"], ["apicoplast"])
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": "pf1"}])
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "ortholopit_label"])


# --------------------------------------------------------------------------- gating
def test_an_uncalibrated_single_donor_transfer_is_not_promoted(tmp_path):
    """Never validated and resting on one species: kept for inspection, not promoted into the column
    the rest of the project reads."""
    _vocab_file(tmp_path, [["pfal", "apicoplast", "apicoplast"]])
    _species_file(tmp_path, "pfal", ["pf1"], ["apicoplast"])
    _master(tmp_path, [{"tgon_gene_id": "t1", "pfal_gene_id": "pf1"}])
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert not bool(out.loc[0, "ortholopit_accepted"])
    assert out.loc[0, "compartment_best"] == "unassigned"


def test_an_uncalibrated_multi_donor_transfer_is_promoted(tmp_path):
    """Two independent species agreeing is worth far more than one: 96.3% against 65.9% measured."""
    _transfer_setup(tmp_path)
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert bool(out.loc[0, "ortholopit_accepted"])
    assert out.loc[0, "compartment_best"] == "apicoplast"
    assert out.loc[0, "compartment_source"] == "orthoLOPIT"


def test_a_measured_label_always_wins_over_a_transferred_one(tmp_path):
    _transfer_setup(tmp_path, tgon_measured="cytosol")
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert out.loc[0, "compartment_source"] == "hyperLOPIT"
    assert out.loc[0, "compartment_best"] == "cytosol"


def test_calibration_and_application_sets_are_disjoint_by_construction(tmp_path):
    """The reason scoring the transfer against measured genes is not circular: a transfer is only ever
    USED where no measured label exists."""
    _transfer_setup(tmp_path, tgon_measured="apicoplast")
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    measured = out.compartment.ne("unassigned")
    used = out.compartment_source.eq("orthoLOPIT")
    assert not (measured & used).any()


def test_a_category_calibrating_below_the_threshold_is_rejected(tmp_path):
    """golgi/secretory transfers at 6.1% and is dropped while proteasome transfers at 100% and is kept.
    Gating on the average would have accepted both."""
    n = L.MIN_CALIBRATION + 5
    _vocab_file(tmp_path, [["tgon", "golgi", "golgi"], ["tgon", "cytosol", "cytosol"],
                           ["pfal", "golgi", "golgi"], ["cpar", "golgi", "golgi"]])
    _species_file(tmp_path, "pfal", [f"pf{i}" for i in range(n + 1)], ["golgi"] * (n + 1))
    _species_file(tmp_path, "cpar", [f"cp{i}" for i in range(n + 1)], ["golgi"] * (n + 1))
    # every dual-labelled gene is really cytosol, so the golgi transfer validates at 0%
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": [f"t{i}" for i in range(n)],
                  "MAP_location": ["cytosol"] * n,
                  "MCMC_location": ["cytosol"] * n}).to_csv(tmp_path / fn, index=False)
    _master(tmp_path, [{"tgon_gene_id": f"t{i}", "pfal_gene_id": f"pf{i}", "cpar_gene_id": f"cp{i}"}
                       for i in range(n + 1)])

    out = L.lopit_labels(str(tmp_path), _nodes([f"t{i}" for i in range(n + 1)]),
                         log=lambda *_: None).set_index("gene_id")
    assert out.loc[f"t{n}", "ortholopit_accuracy"] == pytest.approx(0.0)
    assert not bool(out.loc[f"t{n}", "ortholopit_accepted"])
    assert out.loc[f"t{n}", "compartment_best"] == "unassigned"


def test_a_category_calibrating_above_the_threshold_is_accepted(tmp_path):
    n = L.MIN_CALIBRATION + 5
    _vocab_file(tmp_path, [["tgon", "apicoplast", "apicoplast"],
                           ["pfal", "apicoplast", "apicoplast"],
                           ["cpar", "apicoplast", "apicoplast"]])
    _species_file(tmp_path, "pfal", [f"pf{i}" for i in range(n + 1)], ["apicoplast"] * (n + 1))
    _species_file(tmp_path, "cpar", [f"cp{i}" for i in range(n + 1)], ["apicoplast"] * (n + 1))
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": [f"t{i}" for i in range(n)],
                  "MAP_location": ["apicoplast"] * n,
                  "MCMC_location": ["apicoplast"] * n}).to_csv(tmp_path / fn, index=False)
    _master(tmp_path, [{"tgon_gene_id": f"t{i}", "pfal_gene_id": f"pf{i}", "cpar_gene_id": f"cp{i}"}
                       for i in range(n + 1)])

    out = L.lopit_labels(str(tmp_path), _nodes([f"t{i}" for i in range(n + 1)]),
                         log=lambda *_: None).set_index("gene_id")
    assert out.loc[f"t{n}", "ortholopit_accuracy"] == pytest.approx(1.0)
    assert bool(out.loc[f"t{n}", "ortholopit_accepted"])


def test_a_combination_rarer_than_the_calibration_minimum_is_not_scored(tmp_path):
    """An accuracy computed from three genes is not an accuracy."""
    _transfer_setup(tmp_path, tgon_measured="apicoplast")
    out = L.lopit_labels(str(tmp_path), _nodes(["t1"]), log=lambda *_: None)
    assert pd.isna(out.loc[0, "ortholopit_accuracy"])


def test_the_calibration_report_names_which_categories_were_kept_and_rejected(tmp_path):
    n = L.MIN_CALIBRATION + 2
    _vocab_file(tmp_path, [["tgon", "apicoplast", "apicoplast"],
                           ["pfal", "apicoplast", "apicoplast"],
                           ["cpar", "apicoplast", "apicoplast"]])
    _species_file(tmp_path, "pfal", [f"pf{i}" for i in range(n)], ["apicoplast"] * n)
    _species_file(tmp_path, "cpar", [f"cp{i}" for i in range(n)], ["apicoplast"] * n)
    fn, _, _ = L.SPECIES["tgon"]
    pd.DataFrame({"gene_source_id": [f"t{i}" for i in range(n)],
                  "MAP_location": ["apicoplast"] * n,
                  "MCMC_location": ["apicoplast"] * n}).to_csv(tmp_path / fn, index=False)
    _master(tmp_path, [{"tgon_gene_id": f"t{i}", "pfal_gene_id": f"pf{i}", "cpar_gene_id": f"cp{i}"}
                       for i in range(n)])
    msgs = []
    L.lopit_labels(str(tmp_path), _nodes([f"t{i}" for i in range(n)]), log=msgs.append)
    assert any("calibration on" in m for m in msgs)
    assert any("keep" in m for m in msgs)
