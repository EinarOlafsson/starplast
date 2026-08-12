#!/usr/bin/env python3
"""Physical interactions: crosslink mass spectrometry, the predicted complexes, and the edge lift.

A crosslink is inherently pairwise -- "GRA7 residue 124 joins GRA2 residue 102" -- so it cannot be a
property of one gene. The node table carries only counts; the pairs live in the graph. Most of what is
tested here is that distinction holding, plus the honesty of the model scoring: 60% of the predicted
complexes explain no crosslink at all, and the point of the `model_trustworthy` flag is to say so rather
than to show a confident-looking picture of a pose nothing supports.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import interactions as IN  # noqa: E402


# --------------------------------------------------------------------------- the edge lift
def _edge_table(tmp_path, rows):
    p = tmp_path / "edges.parquet"
    pd.DataFrame(rows).to_parquet(p, index=False)
    return str(p)


def test_a_missing_edge_table_is_reported_and_yields_nothing(tmp_path):
    msgs = []
    edges, detail = IN.load_toxonet(str(tmp_path / "absent.parquet"), ["g1"], log=msgs.append)
    assert edges == {} and detail.empty
    assert any("absent" in m for m in msgs)


def test_edges_referring_to_unknown_genes_are_dropped(tmp_path):
    """An edge to a gene that is not in the node table would index out of range at draw time."""
    label = next(iter(IN.TOXONET_EDGES))
    p = _edge_table(tmp_path, [{"edge_type": label, "src": "g1", "dst": "g2", "weight": 1.0},
                               {"edge_type": label, "src": "g1", "dst": "ghost", "weight": 1.0}])
    edges, _ = IN.load_toxonet(p, ["g1", "g2"], log=lambda *_: None)
    a, b, w, _ = edges[label]
    assert len(a) == 1


def test_self_edges_are_removed(tmp_path):
    label = next(iter(IN.TOXONET_EDGES))
    p = _edge_table(tmp_path, [{"edge_type": label, "src": "g1", "dst": "g1", "weight": 1.0},
                               {"edge_type": label, "src": "g1", "dst": "g2", "weight": 1.0}])
    edges, _ = IN.load_toxonet(p, ["g1", "g2"], log=lambda *_: None)
    a, b, _, _ = edges[label]
    assert not (a == b).any()


def test_a_pair_listed_in_both_directions_becomes_one_undirected_edge(tmp_path):
    """Never trust an upstream table to be free of duplicates: drawn twice, an edge is twice as bright,
    which reads as twice the evidence."""
    label = next(iter(IN.TOXONET_EDGES))
    p = _edge_table(tmp_path, [{"edge_type": label, "src": "g1", "dst": "g2", "weight": 1.0},
                               {"edge_type": label, "src": "g2", "dst": "g1", "weight": 3.0}])
    edges, _ = IN.load_toxonet(p, ["g1", "g2"], log=lambda *_: None)
    a, b, w, _ = edges[label]
    assert len(a) == 1
    assert w[0] == pytest.approx(3.0), "the strongest weight survives the merge"


def test_an_edge_type_with_no_rows_is_absent_rather_than_empty(tmp_path):
    label = next(iter(IN.TOXONET_EDGES))
    p = _edge_table(tmp_path, [{"edge_type": label, "src": "g1", "dst": "g2", "weight": 1.0}])
    edges, _ = IN.load_toxonet(p, ["g1", "g2"], log=lambda *_: None)
    assert set(edges) == {label}


def test_an_edge_type_whose_rows_all_reference_unknown_genes_is_dropped(tmp_path):
    label = next(iter(IN.TOXONET_EDGES))
    p = _edge_table(tmp_path, [{"edge_type": label, "src": "x", "dst": "y", "weight": 1.0}])
    edges, detail = IN.load_toxonet(p, ["g1", "g2"], log=lambda *_: None)
    assert edges == {} and detail.empty


# --------------------------------------------------------------------------- crosslink models
def _starpath(tmp_path, interactions, crosslinks=None, satisfaction=None, cifs=None):
    pd.DataFrame(interactions).to_csv(tmp_path / "starpath_interactions.csv", index=False)
    (tmp_path / "starpath_crosslinks.json").write_text(json.dumps(crosslinks or []))
    if satisfaction is not None:
        d = tmp_path / "starpath_crosslink_mining"
        d.mkdir(exist_ok=True)
        pd.DataFrame(satisfaction).to_csv(d / "crosslink_satisfaction.csv", index=False)
    if cifs is not None:
        d = tmp_path / "starpath_dump" / "cifs"
        d.mkdir(parents=True, exist_ok=True)
        for f in cifs:
            (d / f).write_text("data_")
    return str(tmp_path)


def _inter_row(edge_id=1, a="TGME49_200010", b="TGME49_200020"):
    return {"edge_id": edge_id, "protA_alias": a, "protB_alias": b, "crosslinks_number": 2,
            "identification_score": 0.9, "local_fdr": 0.01, "protA_loc": "dense granule",
            "protB_loc": "dense granule"}


def test_no_starpath_export_is_reported_and_yields_nothing(tmp_path):
    msgs = []
    assert IN.crosslink_models(str(tmp_path), ["g1"], log=msgs.append).empty
    assert any("StarPath export not found" in m for m in msgs)


def test_pairs_are_stored_in_a_canonical_order(tmp_path):
    """A crosslink is pairwise and undirected; storing it twice under two orders would double it."""
    base = _starpath(tmp_path, [_inter_row(a="TGME49_200020", b="TGME49_200010")])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    assert df.loc[0, "gene_a"] == "TGME49_200010"
    assert df.loc[0, "gene_b"] == "TGME49_200020"


def test_a_pair_naming_a_gene_we_do_not_have_is_skipped(tmp_path):
    """StarPath uses RH88 accessions whose numbering does NOT correspond to ME49, so a mapping failure
    shows up here as an unknown gene rather than as a wrong one."""
    base = _starpath(tmp_path, [_inter_row(a="TGRH88_016370")])
    assert IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"],
                               log=lambda *_: None).empty


def test_a_self_pair_is_skipped(tmp_path):
    base = _starpath(tmp_path, [_inter_row(a="TGME49_200010", b="TGME49_200010")])
    assert IN.crosslink_models(base, ["TGME49_200010"], log=lambda *_: None).empty


def test_residue_positions_are_attached_to_the_pair(tmp_path):
    """What was actually joined to what: the residues are the measurement, the model is a picture."""
    base = _starpath(tmp_path, [_inter_row()],
                     crosslinks=[{"id": 1, "crosslinks": [{"pos_a": 124, "pos_b": 102}]}])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    assert json.loads(df.loc[0, "crosslink_positions"]) == [[124, 102]]


def test_positions_for_an_unknown_pair_are_ignored(tmp_path):
    base = _starpath(tmp_path, [_inter_row(edge_id=1)],
                     crosslinks=[{"id": 999, "crosslinks": [{"pos_a": 1, "pos_b": 2}]}])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    assert "crosslink_positions" not in df.columns or df.crosslink_positions.isna().all()


def test_a_model_is_trustworthy_only_when_confident_and_crosslink_consistent(tmp_path):
    """Both conditions, deliberately: a confident pose that explains none of the crosslinks is a
    convincing picture of something the data does not support."""
    base = _starpath(tmp_path,
                     [_inter_row(1), _inter_row(2, "TGME49_200030", "TGME49_200040"),
                      _inter_row(3, "TGME49_200050", "TGME49_200060")],
                     satisfaction=[{"id": 1, "frac_satisfied": 0.9, "min_ca": 8.0,
                                    "chai_iptm": 0.8, "chai_plddt": 90.0},
                                   {"id": 2, "frac_satisfied": 0.0, "min_ca": 40.0,
                                    "chai_iptm": 0.9, "chai_plddt": 95.0},
                                   {"id": 3, "frac_satisfied": 0.9, "min_ca": 8.0,
                                    "chai_iptm": 0.2, "chai_plddt": 40.0}])
    ids = [f"TGME49_2000{i}0" for i in range(1, 7)]
    df = IN.crosslink_models(base, ids, log=lambda *_: None).set_index("starpath_id")
    assert bool(df.loc[1, "model_trustworthy"])
    assert not bool(df.loc[2, "model_trustworthy"]), "confident but explains no crosslink"
    assert not bool(df.loc[3, "model_trustworthy"]), "consistent but not a confident interface"


def test_the_report_says_how_many_models_explain_nothing(tmp_path):
    """60% of them explain no crosslink at all, and that number is the point."""
    base = _starpath(tmp_path, [_inter_row(1)],
                     satisfaction=[{"id": 1, "frac_satisfied": 0.0, "min_ca": 40.0,
                                    "chai_iptm": 0.9, "chai_plddt": 95.0}])
    msgs = []
    IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=msgs.append)
    assert any("explain no" in m for m in msgs)


def test_local_model_files_are_linked_when_present(tmp_path):
    base = _starpath(tmp_path, [_inter_row(1)], cifs=["1_model_a.cif", "1_model_b.cif",
                                                      "2_other.cif", "notanumber.cif"])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    assert df.loc[0, "n_models"] == 2
    assert df.loc[0, "model_files"] == "1_model_a.cif;1_model_b.cif"


def test_absent_model_files_are_reported_rather_than_faked(tmp_path):
    base = _starpath(tmp_path, [_inter_row(1)])
    msgs = []
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=msgs.append)
    assert pd.isna(df.loc[0, "n_models"])
    assert any("not on this machine" in m for m in msgs)


def test_scoring_columns_exist_even_when_nothing_was_scored(tmp_path):
    """Downstream code selects them by name; a missing column fails at the join, far from the cause."""
    base = _starpath(tmp_path, [_inter_row(1)])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    for c in ("frac_satisfied", "min_ca", "chai_iptm", "chai_plddt", "n_models"):
        assert c in df.columns


def test_satisfaction_rows_for_unknown_pairs_are_ignored(tmp_path):
    base = _starpath(tmp_path, [_inter_row(1)],
                     satisfaction=[{"id": 999, "frac_satisfied": 0.9, "min_ca": 8.0,
                                    "chai_iptm": 0.8, "chai_plddt": 90.0}])
    df = IN.crosslink_models(base, ["TGME49_200010", "TGME49_200020"], log=lambda *_: None)
    assert df.frac_satisfied.isna().all()


# --------------------------------------------------------------------------- per-gene attributes
def test_degree_counts_both_ends_of_every_edge():
    nodes = pd.DataFrame({"gene_id": ["g1", "g2", "g3"]})
    edges = {"xlms": (np.array([0, 0]), np.array([1, 2]), np.array([1.0, 1.0]), np.array([1.0, 1.0]))}
    IN.gene_attributes(edges, pd.DataFrame(), nodes)
    assert nodes.n_xlink_partners.tolist() == [2, 1, 1]


def test_an_absent_edge_type_gives_zero_not_missing():
    """Zero crosslink partners is a measurement; missing would read as "not looked at"."""
    nodes = pd.DataFrame({"gene_id": ["g1", "g2"]})
    IN.gene_attributes({}, pd.DataFrame(), nodes)
    assert nodes.n_xlink_partners.tolist() == [0, 0]
    assert nodes.n_struct_similar.tolist() == [0, 0]
    assert nodes.n_ipms_partners.tolist() == [0, 0]


def test_the_best_model_agreement_is_taken_across_a_genes_pairs():
    nodes = pd.DataFrame({"gene_id": ["g1", "g2", "g3"]})
    models = pd.DataFrame({"gene_a": ["g1", "g1"], "gene_b": ["g2", "g3"],
                           "frac_satisfied": [0.2, 0.9]})
    IN.gene_attributes({}, models, nodes)
    assert nodes.set_index("gene_id").best_model_agreement.loc["g1"] == pytest.approx(0.9)
    assert nodes.set_index("gene_id").best_model_agreement.loc["g2"] == pytest.approx(0.2)


def test_a_gene_with_no_model_has_no_agreement_rather_than_zero():
    """Zero would say "no predicted complex explains it", which is a different claim from "none was
    made"."""
    nodes = pd.DataFrame({"gene_id": ["g1"]})
    IN.gene_attributes({}, pd.DataFrame(), nodes)
    assert nodes.best_model_agreement.isna().all()


def test_models_without_any_scored_row_leave_agreement_missing():
    nodes = pd.DataFrame({"gene_id": ["g1"]})
    models = pd.DataFrame({"gene_a": ["g1"], "gene_b": ["g2"], "frac_satisfied": [np.nan]})
    IN.gene_attributes({}, models, nodes)
    assert nodes.best_model_agreement.isna().all()
