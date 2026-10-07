"""Provenance cannot turn ambiguous mappings or installed summaries into experiments."""
from dataclasses import replace
import json

import pytest

from starplast import datasets as D, organisms as O
from starplast.provenance import (SourceFile, MeasurementTrace, Transform, MappingAudit,
                                  audit_mapping, read_traces, trace_sources, write_traces)


def _source(key, columns, derived=(), organism=O.TOXOPLASMA):
    return D.Dataset(key, key, "reference", "test", "test", columns=columns,
                     derived_from=derived, organism=organism)


def test_file_identity_detects_changed_content_and_refuses_directories(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text("gene\tvalue\na\t1\n")
    record = SourceFile.inspect(path, "processed_input")
    record.verify()
    path.write_text("gene\tvalue\na\t2\n")
    with pytest.raises(ValueError, match="content"):
        record.verify()
    with pytest.raises(ValueError, match="regular"):
        SourceFile.inspect(tmp_path, "raw_input")
    with pytest.raises(ValueError):
        SourceFile(str(path), "not_a_hash", 1, "raw_input")


def test_mapping_counts_unique_entities_ambiguity_loss_and_reverse_cardinality():
    report = audit_mapping(["a", "a", "b", "c", "d", "e", None, ""],
                           [("a", "protein1"), ("a", "protein1"), ("b", "protein1"),
                            ("c", "protein2"), ("c", "protein3"), ("outside", "protein4")],
                           source_organism=O.HUMAN, target_organism=O.HUMAN,
                           source_namespace="Ensembl gene", target_namespace="UniProt protein",
                           source_version="reference_version", target_version="reference_version",
                           reference="synthetic_fixture")
    assert report.source_entities == 5
    assert (report.mapped_entities, report.ambiguous_entities, report.unmapped_entities) == (2, 1, 2)
    assert report.target_entities == report.many_to_one_targets == 1
    assert report.source_identifier_missing_rows == 2
    with pytest.raises(ValueError, match="reconcile"):
        replace(report, unmapped_entities=0)
    with pytest.raises(ValueError, match="explicit"):
        replace(report, target_version="")


def test_ancestry_is_species_qualified_and_keeps_unknown_and_ambiguous_inputs():
    a = _source("a", ("input",))
    unrelated = _source("other_species", ("input",), organism=O.FALCIPARUM)
    b = _source("b", ("output",), ("input", "unknown"))
    report = trace_sources("b", [a, b, unrelated])
    assert report["sources"] == ["a", "b"]
    assert report["gaps"][0]["column"] == "unknown"
    duplicate = _source("duplicate", ("input",))
    assert trace_sources("b", [a, b, duplicate])["gaps"][0]["status"] == "ambiguous"
    with pytest.raises(ValueError, match="Cycle"):
        trace_sources("a", [_source("a", ("input",), ("output",)), b])


def test_trace_roundtrip_preserves_units_steps_gaps_and_tamper_identity(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text("a\t1\n")
    record = SourceFile.inspect(path, "processed_input")
    step = Transform("normalize", "module.normalize", "0" * 64, (record.path,), ("bmdm_tpm",), {"replicates": 3})
    trace = MeasurementTrace("source", O.MOUSE, O.MOUSE, "protein", "bmdm_tpm",
                             "direct_experiment", "TPM", source_files=(record,), transforms=(step,),
                             gaps=("mapping_version_unresolved",))
    target = tmp_path / "traces.json"
    write_traces([trace], target)
    assert read_traces(target) == [trace]
    with pytest.raises(ValueError, match="Duplicate"):
        write_traces([trace, trace], tmp_path / "duplicate.json")
    saved = json.loads(target.read_text())
    saved[0]["quantity_unit"] = "copies_per_cell"
    target.write_text(json.dumps(saved))
    with pytest.raises(ValueError, match="identity"):
        read_traces(target)


def test_transfer_requires_mapping_or_an_explicit_gap():
    with pytest.raises(ValueError, match="transfer gap"):
        MeasurementTrace("source", O.FALCIPARUM, O.FALCIPARUM, "gene", "fitness",
                         evidence_grade="orthology_transfer")
    trace = MeasurementTrace("source", O.FALCIPARUM, O.FALCIPARUM, "gene", "fitness",
                             evidence_grade="orthology_transfer", gaps=("transfer_mapping_unresolved",))
    assert trace.quantity_unit == "unresolved" and trace.redistribution == "unresolved"


def test_registry_scope_does_not_claim_the_original_assay_species():
    trace = MeasurementTrace("legacy_transfer", O.FALCIPARUM, O.FALCIPARUM, "gene", "fitness")
    assert trace.registry_organism == O.FALCIPARUM
    assert trace.source_species == () and trace.source_species_status == "unresolved"
    verified = replace(trace, source_species=("Plasmodium berghei",),
                       source_species_status="primary_metadata_verified", source_species_reference="source artifact hash")
    assert verified.source_species != (O.get(O.FALCIPARUM).species,)
    with pytest.raises(ValueError, match="requires a source reference"):
        replace(trace, source_species_status="primary_metadata_verified")


def test_transforms_refuse_missing_implementation_identity_or_nonjson_parameters():
    with pytest.raises(ValueError):
        Transform("step", "module.function", "wrong", ("input",), ("output",))
    with pytest.raises(ValueError):
        Transform("step", "module.function", "0" * 64, ("input",), ("output",), {"cutoff": float("nan")})
