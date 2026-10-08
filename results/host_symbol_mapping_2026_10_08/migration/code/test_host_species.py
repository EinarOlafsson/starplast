"""Host provenance, migration and readers must keep human and mouse evidence separate."""
from dataclasses import replace
import json
from pathlib import Path

import pandas as pd
import pytest

from starplast import datasets, deposits, host, organisms, slots


def legacy():
    """A human zero, a mouse negative and an identity-only row are all real records."""
    return pd.DataFrame({"host_id": ["H", "M", "I"], "host_name": ["human", "mouse", None],
                         "fibroblast_tpm": [0.0, None, None],
                         "bmdm_surface_detected": pd.array([None, False, None], dtype="boolean")})


def test_organism_is_required_and_does_not_come_from_the_key():
    with pytest.raises(TypeError, match="organism"):
        datasets.Dataset("host_unknown", "unknown", "gene", "table", "measurement")
    with pytest.raises(ValueError, match="unknown dataset organism"):
        datasets.Dataset("unknown", "unknown", "gene", "table", "measurement", organism="host")
    for key, code in (("host_gtex_transcriptome", "Hs"),
                      ("host_bmdm_baseline", "Mm"), ("pv_host_uptake", "Hs"),
                      ("pf_spatial_proteome", "Pf")):
        dataset = datasets._BY_KEY[key]
        assert datasets.organism_of(dataset) == code
        assert datasets.organism_of(replace(dataset, key="renamed")) == code
        for column in dataset.columns:
            assert datasets.provenance(column, organism=code) is not None
    assert datasets.provenance("fibroblast_tpm", organism="Tg") is None
    assert datasets.provenance("brain_tpm", organism="Hs") is None
    assert datasets.derived_sources("stage_enriched_derived", organism="Hs") == ()


def test_every_deposit_agrees_with_its_dataset_organism():
    for deposit in deposits.DEPOSITS:
        assert deposit.organism == datasets._BY_KEY[deposit.key].organism


def test_split_keeps_zeros_negatives_and_identity_only_rows():
    original = legacy()
    tables = host.split_legacy_table(original, {"Hs": {"I"}})
    assert list(tables["Hs"].host_id) == ["H", "I"]
    assert list(tables["Mm"].host_id) == ["M"]
    assert tables["Hs"].fibroblast_tpm.iloc[0] == 0
    assert not bool(tables["Mm"].bmdm_surface_detected.iloc[0])
    assert "fibroblast_tpm" not in tables["Mm"]
    assert "bmdm_surface_detected" not in tables["Hs"]
    for table in tables.values():
        expected = original.set_index("host_id").loc[table.host_id, table.columns[1:]]
        pd.testing.assert_frame_equal(table.set_index("host_id"), expected)


def test_split_refuses_to_guess_or_discard_rows():
    with pytest.raises(ValueError, match="exactly one organism"):
        host.split_legacy_table(legacy())
    with pytest.raises(ValueError, match="exactly one organism"):
        host.split_legacy_table(legacy(), {"Hs": {"I", "M"}})
    frame = legacy()
    frame.loc[1, "fibroblast_tpm"] = 2.0
    with pytest.raises(ValueError, match="exactly one organism"):
        host.split_legacy_table(frame, {"Hs": {"I"}})
    with pytest.raises(ValueError, match="without organism provenance"):
        host.split_legacy_table(legacy().assign(unregistered=3))
    with pytest.raises(ValueError, match="unique, nonmissing"):
        host.split_legacy_table(pd.concat([legacy(), legacy()]))


def test_host_deposit_reader_selects_one_species(monkeypatch):
    calls = []
    monkeypatch.setattr(deposits, "DEPOSITS", (
        deposits.Deposit("human", "Hs", None), deposits.Deposit("mouse", "Mm", None)))

    def read(base, key):
        calls.append(key)
        return pd.DataFrame({"host_id": [key], "host_name": [key], "value": [1.]})

    monkeypatch.setattr(deposits, "_read", read)
    assert list(deposits.host_columns("unused", "Mm").host_id) == ["mouse"]
    assert calls == ["mouse"]
    with pytest.raises(ValueError, match="unknown host"):
        deposits.host_columns("unused", "host")
    with pytest.raises(ValueError, match="unknown host"):
        host.tissue_references("unused", "host")


def test_shipped_hosts_preserve_the_migration_totals_and_are_disjoint():
    from starplast import paths
    tables = host.shipped_tables()
    report = json.loads(Path(paths.cache_file("host_species_migration.json")).read_text())
    withdrawals = host.projection_withdrawals()
    assert set(tables) == set(organisms.HOST_TABLES)
    assert sum(len(t) for t in tables.values()) >= report["rows_before"]
    assert not set(tables["Hs"].host_id) & set(tables["Mm"].host_id)
    for code, table in tables.items():
        assert table.host_id.is_unique
        for column, count in report["tables"][code]["coverage"].items():
            recorded = withdrawals[(withdrawals.organism == code) & (withdrawals.field == column)]
            # Entity-assignment corrections reconcile through exact retained gene scores.
            assert int(table[column].notna().sum()) + len(recorded) >= count
        for column in set(table) - {"host_id", "host_name"}:
            assert datasets.provenance(column, code).organism == code
    for name in (host.BRIDGE_TABLE, organisms.get(organisms.FALCIPARUM).host_bridges):
        bridges = pd.read_parquet(paths.cache_file(name))
        assert set(bridges.host_id) <= set(tables["Hs"].host_id)


def test_withdrawal_ledger_cannot_excuse_an_unrelated_or_altered_gene_measurement(tmp_path):
    from starplast import paths
    import shutil
    folder = Path(paths.cache_file(host.HOST_TABLE)).parent
    ledger_path = folder / "host_projection_withdrawals.json"
    if not ledger_path.exists():
        pytest.skip("projection revision not yet installed")
    ledger = json.loads(ledger_path.read_text())
    for filename in (ledger["gene_evidence_file"], ledger["target_file"]):
        shutil.copy2(folder / filename, tmp_path / filename)
    ledger["records"][0]["old_value"] += 1
    (tmp_path / ledger_path.name).write_text(json.dumps(ledger))
    with pytest.raises(ValueError, match="declared original value encoding"):
        host.projection_withdrawals(tmp_path)


def test_slot_reader_uses_one_host_table_and_refuses_ambiguous_sources():
    tables = host.split_legacy_table(legacy(), {"Hs": {"I"}})
    slot = slots.Slot("Tg", "human transcriptome", "a", "c", "host_gene",
                      ("fibroblast_tpm",), "one")
    assert slots.is_filled(slot, tables={"host_gene": tables})
    assert list(slots.unit_table(slot, {"host_gene": tables}).host_id) == ["H", "I"]
    assert not slots.is_filled(slot, tables={"host_gene": {"Mm": tables["Mm"]}})
    with pytest.raises(ValueError, match="more than one organism"):
        slots.unit_table(slot, {"host_gene": {"Hs": tables["Hs"], "Mm": tables["Hs"]}})


def test_host_migration_refuses_to_overwrite_later_acquisitions(tmp_path, monkeypatch):
    from scripts import split_host_tables as migration
    source = tmp_path / "legacy.parquet"
    legacy().to_parquet(source, index=False)
    monkeypatch.setattr(migration, "identity_memberships", lambda *a: ({"Hs": {"I"}}, []))
    migration.migrate(source, tmp_path, tmp_path, write=True)
    path = tmp_path / host.HOST_TABLES["Hs"]
    newer = pd.read_parquet(path).assign(new_measurement=7.0)
    newer.to_parquet(path, index=False)
    with pytest.raises(AssertionError):
        migration.migrate(source, tmp_path, tmp_path, write=True)
    pd.testing.assert_frame_equal(pd.read_parquet(path), newer)


def test_host_merge_cannot_hide_lost_values_behind_new_coverage(tmp_path, monkeypatch):
    from scripts import add_deposits
    path = tmp_path / "human.parquet"
    existing = pd.DataFrame({"host_id": ["old"], "host_name": ["old"], "value": [2.]})
    existing.to_parquet(path, index=False)
    replacement = pd.DataFrame({"host_id": ["new"], "host_name": ["new"], "value": [3.]})
    monkeypatch.setattr(deposits, "host_columns", lambda *args: replacement)
    assert add_deposits.merge_host(path, str(tmp_path), "Hs", log=lambda *a: None) == 2
    pd.testing.assert_frame_equal(pd.read_parquet(path), existing)


def test_slot_atlas_and_audit_read_separate_host_tables(monkeypatch):
    from scripts import generate_slot_table as generator
    from starplast import slot_tree
    tables = host.split_legacy_table(legacy(), {"Hs": {"I"}})
    frame = generator.host_row_space(tables, "human fibroblast", ("fibroblast_tpm",))
    assert list(frame.host_id) == ["H"]
    assert generator.host_row_space({"Mm": tables["Mm"]}, "human fibroblast").empty
    with pytest.raises(ValueError, match="more than one organism"):
        generator.host_row_space({"Hs": tables["Hs"], "Mm": tables["Hs"]}, "human fibroblast")
    monkeypatch.setattr(slots, "all_slots", lambda *a: ())
    tables["Mm"] = tables["Mm"].assign(unclaimed_mouse=1)
    audit = slot_tree.audit("Tg", tables={"host_gene": tables})
    assert audit["orphan"] == 3  # one human measurement, two mouse measurements
