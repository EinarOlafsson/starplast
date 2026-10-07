"""Inventory coverage distinguishes stored negatives, unavailable sources and unknown causes."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from starplast import datasets as D, organisms as O, slots as S
from starplast.evidence import Observation
from starplast.inventory import build_inventory


def _source(key="synthetic", columns=("value",), **kwargs):
    return D.Dataset(key, "Synthetic source", "DNA", "synthetic assay", "Synthetic measurements",
                     organism=kwargs.pop("organism", O.TOXOPLASMA), columns=columns, **kwargs)


def _row(frame, key="synthetic", unit=None):
    chosen = frame[frame.source_id == key]
    if unit:
        chosen = chosen[chosen.unit == unit]
    assert len(chosen) == 1
    return chosen.iloc[0]


def test_all_registered_sources_remain_visible_without_installed_data():
    report = build_inventory(sources=D.REGISTRY)
    assert set(report.source_id) == {d.key for d in D.REGISTRY}
    assert not report.duplicated(["source_id", "organism", "unit"]).any()
    assert report.status.eq("unavailable").all()
    assert report.stored_fraction.isna().all()


def test_storage_coverage_retains_zero_false_and_missing_unknowns():
    table = pd.DataFrame({"value": [0., 1., np.nan, np.inf],
                          "detected": pd.Series([False, True, pd.NA, pd.NA], dtype="boolean")})
    source = _source(columns=("value", "detected"))
    result = _row(build_inventory({(O.TOXOPLASMA, "gene"): table}, sources=[source], catalog=[]))
    assert result.table_rows == 4 and result.stored_any_rows == 2
    assert result.stored_fraction == .5 and result.missing_unknown_rows == 2
    assert result.false_cells == 1
    assert result.unmeasured_records is None and result.unmapped_records is None
    assert result.access_status == "not_checked"  # no claim about network access


def test_partial_outputs_and_missing_tables_do_not_become_zero_coverage():
    source = _source(columns=("value", "missing"))
    partial = _row(build_inventory({(O.TOXOPLASMA, "gene"): pd.DataFrame({"value": [1., np.nan]})},
                                  sources=[source], catalog=[]))
    assert partial.status == "partial" and partial.missing_columns == ["missing"]
    absent = _row(build_inventory(sources=[source], catalog=[]))
    assert absent.status == "unavailable" and absent.table_rows is None
    assert absent.stored_fraction is None and absent.missing_unknown_rows is None
    missing_column = _row(build_inventory({(O.TOXOPLASMA, "gene"): pd.DataFrame({"other": [1]})},
                                         sources=[source], catalog=[]))
    assert missing_column.table_rows == 1 and missing_column.stored_any_rows is None


def test_context_and_evidence_facets_come_from_the_existing_slot_matcher():
    slot = S.Slot(O.TOXOPLASMA, "synthetic slot", "DNA", "stage", "gene", ("value",), "one",
                  context_path=("stage", "synthetic"), evidence_path=("measured", "synthetic"))
    result = _row(build_inventory(sources=[_source()], catalog=[slot]))
    assert result.contexts == ["stage / synthetic"]
    assert result.evidence_families == ["measured / synthetic"]


def test_source_only_host_protein_and_metabolite_storage_units_stay_distinct():
    host_source = _source("host", path="starplast/data/" + O.HOST_TABLES[O.HUMAN])
    metabolite = S.Slot(O.TOXOPLASMA, "metabolite", "translation", "culture", "metabolite",
                        ("value",), "one")
    hosts = build_inventory({(O.HUMAN, "protein"): pd.DataFrame({"value": [1., 2.]})},
                            sources=[host_source], catalog=[])
    host = _row(hosts, "host")
    assert host.organism == O.HUMAN and host.unit == "protein" and host.table_rows == 2
    assert host.registry_organism == O.TOXOPLASMA  # legacy exploration address is retained
    assert _row(build_inventory(sources=[_source()], catalog=[metabolite])).unit == "metabolite"
    raw = _row(build_inventory(sources=[_source(columns=())], raw_paths={"synthetic": "/synthetic"}))
    assert raw.status == "source_only" and raw.unit == "source_records" and raw.table_rows is None


def test_known_missingness_records_are_not_invented_from_null_summary_cells():
    records = [Observation("SYNTHETIC_GENE", "value", None, "synthetic", O.TOXOPLASMA,
                           missing_state=state) for state in ("not_assayed", "ambiguous_mapping", "unknown")]
    records.append(Observation("SYNTHETIC_GENE2", "value", False, "synthetic", O.TOXOPLASMA))
    result = _row(build_inventory(sources=[_source()], catalog=[], observations=records))
    assert result.unmeasured_records == 1 and result.unmapped_records == 1
    assert result.unknown_missing_records == 1 and result.observed_records == 1
    assert result.explicit_negative_records == 1
    assert result.failed_qc_records == 0 and result.stored_any_rows is None
    with pytest.raises(ValueError, match="outside"):
        build_inventory(sources=[_source()], observations=[replace(records[0], organism=O.FALCIPARUM)])


def test_pair_counts_do_not_invent_assayed_pair_universes(tmp_path):
    path = tmp_path / "synthetic.npz"
    np.savez(path, layer__a=np.array([0, 1]), layer__b=np.array([1, 2]))
    with np.load(path) as graph:
        rows = build_inventory(graphs={O.TOXOPLASMA: graph},
                               sources=[_source(columns=("edge:layer", "value"))], catalog=[])
    pair = _row(rows, unit="pair")
    assert pair.pair_records == 2 and pair.assayed_pair_denominator is None
    assert pair.stored_fraction is None and pair.status == "installed"
    assert len(rows) == 2  # one source with two declared storage units


def test_graph_inventory_matches_the_actual_shipped_edge_schema():
    with np.load(O.graph_path(O.FALCIPARUM)) as graph:
        source = D.get("pf_crosslink_ms")
        result = _row(build_inventory(graphs={O.FALCIPARUM: graph}, sources=[source]), source.key)
        assert result.status == "installed"
        assert result.pair_records == len(graph["xlms__a"]) == len(graph["xlms__b"])


def test_graph_inventory_refuses_mismatched_endpoint_arrays(tmp_path):
    path = tmp_path / "synthetic_bad.npz"
    np.savez(path, layer__a=np.array([0, 1]), layer__b=np.array([1]))
    with np.load(path) as graph, pytest.raises(ValueError, match="paired"):
        build_inventory(graphs={O.TOXOPLASMA: graph}, sources=[_source(columns=("edge:layer",))])


def test_unattributed_bridges_are_reported_without_assigning_them_to_a_paper():
    table = pd.DataFrame({"bridge": ["host", "host"]})
    source = _source(columns=("bridge:host",))
    result = _row(build_inventory(bridges={O.TOXOPLASMA: table}, sources=[source]))
    assert result.status == "unattributed" and result.pair_records is None
    assert result.unattributed_pair_records == 2
    table["source"] = ["synthetic", "another-source"]
    result = _row(build_inventory(bridges={O.TOXOPLASMA: table}, sources=[source]))
    assert result.status == "installed" and result.pair_records == 1


def test_refusals_are_retained_only_for_the_rejected_question():
    refusal = {"source_id": "candidate:SYNTHETIC", "organism": O.TOXOPLASMA, "unit": "gene",
               "question": "synthetic wrong question", "reason": "Synthetic assay measures another quantity"}
    report = build_inventory(sources=[_source()], refusals=[refusal])
    rejected = _row(report, "candidate:SYNTHETIC")
    assert rejected.status == "rejected" and rejected.rejection_reason == refusal["reason"]
    assert _row(report).status == "unavailable"
    with pytest.raises(ValueError):
        build_inventory(sources=[], refusals=[dict(refusal, reason="")])


def test_wrong_organism_or_duplicate_ids_cannot_supply_a_denominator():
    foreign = pd.DataFrame({"gene_id": ["PF3D7_0100001"], "value": [1.]})
    with pytest.raises(ValueError, match="explicit organism"):
        build_inventory({(O.TOXOPLASMA, "gene"): foreign}, sources=[_source()])
    with pytest.raises(ValueError, match="unique"):
        build_inventory({(O.TOXOPLASMA, "gene"): pd.DataFrame({"gene_id": ["TGME49_100001"] * 2})})
    with pytest.raises(ValueError, match="unique"):
        build_inventory({(O.HUMAN, "protein"): pd.DataFrame({"host_id": [None]})})


def test_empty_and_duplicate_registries_and_conflicting_units():
    assert build_inventory(sources=[]).empty
    with pytest.raises(ValueError, match="unique"):
        build_inventory(sources=[_source(), _source()])
    a = S.Slot(O.TOXOPLASMA, "a", "DNA", "stage", "gene", ("value",), "one")
    with pytest.raises(ValueError, match="Conflicting"):
        build_inventory(sources=[_source()], catalog=[a, replace(a, unit="metabolite")])
