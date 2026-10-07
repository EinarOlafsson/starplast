"""Freeze traces for every inventory output and reproduce a raw-to-host pilot.

Legacy units, assay lineage, transformation parameters and licenses stay unknown
where no source-specific audit is available. The measured pilot is mouse BMDM
baseline expression, including all mapping losses and a value-by-value comparison
against the installed table. This census changes no data or inference outputs.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import datasets as D, deposits, host, organisms as O, paths  # noqa: E402
from starplast.provenance import (MeasurementTrace, SourceFile, Transform, audit_mapping,
                                  read_traces, trace_sources, write_traces)  # noqa: E402


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def bmdm_pilot(root: Path) -> tuple:
    """Measure every mouse Ensembl-to-reviewed-UniProt loss and reproduce installed TPM."""
    source = D.get("host_bmdm_baseline")
    raw = root.joinpath(*deposits.BMDM)
    mapping, entries = [root / host.UNIPROT_ROOT / name for name in host.UNIPROT_IDMAP["mouse"]]
    data = pd.read_csv(raw, sep="\t")
    source_ids = data.geneID.astype(str).str.split(".").str[0].tolist()
    with gzip.open(entries, "rt") as stream:
        next(stream)
        reviewed = {line.split("\t", 1)[0] for line in stream}
    pairs = []
    with gzip.open(mapping, "rt") as stream:
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 3 and fields[1] == "Ensembl":
                accession = fields[0].split("-")[0]
                if accession in reviewed:
                    pairs.append((fields[2].split(".")[0], accession))
    mapping_version = _sha(mapping) + ":" + _sha(entries)
    audit = audit_mapping(source_ids, pairs, source_organism=O.MOUSE, target_organism=O.MOUSE,
                          source_namespace="Ensembl gene", target_namespace="reviewed UniProt protein",
                          source_version="source file SHA256:" + _sha(raw), target_version="mapping SHA256:" + mapping_version,
                          reference=str(mapping) + " and " + str(entries))
    raw_file = SourceFile.inspect(raw, "processed_input", source.url or "")
    urls, sums = raw.parent / "URLS.txt", raw.parent / "SHA256SUMS.txt"
    if urls.is_file() and sums.is_file() and source.url in urls.read_text().splitlines():
        declared = [line.split()[0] for line in sums.read_text().splitlines() if line.split()[-1] == raw.name]
        if declared != [raw_file.sha256]:
            raise ValueError("BMDM raw source no longer matches its original archive receipt")
        raw_file = replace(raw_file, association="archive_url_and_SHA256SUMS_verified")
    files = (raw_file, SourceFile.inspect(mapping, "mapping_reference"), SourceFile.inspect(entries, "mapping_reference"))
    reproduced = deposits.bmdm_baseline(str(root)).set_index("host_id")
    installed = pd.read_parquet(paths.cache_file(O.HOST_TABLES[O.MOUSE])).set_index("host_id")
    expected = installed.bmdm_tpm.dropna()
    missing = sorted(set(expected.index) - set(reproduced.index))
    if missing or expected.empty:
        raise ValueError("BMDM source reproduction lost installed nonmissing measurements")
    raw_delta = reproduced.bmdm_tpm.reindex(expected.index) - expected
    # derive_all writes %.10g TSV; reproduce that recorded serialization exactly,
    # rather than broadening tolerances to conceal a failed value comparison.
    encoded = reproduced[["bmdm_tpm"]].to_csv(sep="\t", float_format="%.10g")
    serialized = pd.read_csv(io.StringIO(encoded), sep="\t", index_col=0).bmdm_tpm
    delta = serialized.reindex(expected.index) - expected
    if not np.isfinite(delta).all() or not (delta == 0).all():
        raise ValueError("BMDM source reproduction differs from installed measurements")
    m0 = [column for column in data if column.startswith("M0_")]
    if len(m0) != 3:
        raise ValueError("The declared BMDM baseline pilot requires three unstimulated replicates")
    transform = Transform("FPKM to TPM per sample; mean three M0 replicates; map and retain best expressed row",
                          "starplast.deposits.bmdm_baseline", _sha(ROOT / "starplast/deposits.py"),
                          tuple(item.path for item in files), ("unrounded:bmdm_tpm",),
                          {"sample_columns": m0, "sample_normalization": "FPKM / sum(all source FPKM) * 1e6",
                           "ambiguous_Ensembl": "withheld", "duplicate_target": "highest mean TPM",
                           "namespace_normalization": "remove explicitly recorded Ensembl version and UniProt isoform suffix"})
    serialization = Transform("Serialize keyed deposit TSV using the recorded ten-significant-digit format",
                              "starplast.deposits.derive_all", _sha(ROOT / "starplast/deposits.py"),
                              ("unrounded:bmdm_tpm",), ("deposit_host_bmdm_baseline.tsv:bmdm_tpm",),
                              {"float_format": "%.10g", "separator": "tab"})
    merge = Transform("Read keyed deposit TSV and join the installed host protein table",
                      "starplast.host.merge_tissue", _sha(ROOT / "starplast/host.py"),
                      ("deposit_host_bmdm_baseline.tsv:bmdm_tpm",), ("bmdm_tpm",),
                      {"join_key": "reviewed UniProt accession", "source_species": O.MOUSE})
    summary = {"source_rows": len(data), "source_mapping": asdict(audit), "reproduced_rows": len(reproduced),
               "installed_nonmissing_rows": len(expected), "max_abs_difference": float(delta.abs().max()),
               "max_before_recorded_TSV_serialization": float(raw_delta.abs().max()),
               "serialization_policy": "deposits.derive_all float_format=%.10g; exact equality after read",
               "missing_installed_ids": missing, "raw_file_association": raw_file.association,
               "publication_lineage": "registry association remains subject to origin-paper review"}
    return files, audit, (transform, serialization, merge), summary


def build(root: Path, bindings_path: Path, output: Path) -> dict:
    """Account for every inventory output with typed provenance or explicit unresolved fields."""
    if output.exists():
        raise ValueError("Use a new provenance output directory")
    inventory_path = ROOT / "results/information_inventory_2026_10_07/inventory.json"
    publication_path = ROOT / "results/dataset_selection_2026_10_07/source_publications.csv"
    inventory = json.loads(inventory_path.read_text())
    publications = pd.read_csv(publication_path, keep_default_na=False).set_index("dataset_id").to_dict("index")
    bindings = {key: SourceFile(**value) for key, value in json.loads(bindings_path.read_text()).items()}
    for item in bindings.values():
        item.verify()
    inputs = [inventory_path, publication_path, bindings_path, Path(__file__),
              *(ROOT / "starplast" / name for name in ("provenance.py", "datasets.py", "deposits.py", "host.py")),
              Path(paths.cache_file(O.HOST_TABLES[O.MOUSE]))]
    hashes = {str(path): _sha(path) for path in inputs}
    files, mapping, pilot_steps, pilot = bmdm_pilot(root)
    graphs, traces = {}, []
    for row in inventory:
        if row["status"] == "rejected":
            continue  # refusals are proposals, not measurements
        source = D.get(row["source_id"])
        try:
            graph = trace_sources(source.key, D.REGISTRY)
        except ValueError as error:
            graph = {"source_id": source.key, "sources": [source.key], "edges": [],
                     "gaps": [{"status": "invalid_legacy_derivation", "reason": str(error)}]}
        graphs[source.key] = graph
        metadata = publications[source.key]
        source_files = (bindings[source.key],) if source.key in bindings else ()
        gaps = ["license_unresolved", "redistribution_unresolved", "source_publication_association_unverified"]
        if not source_files:
            gaps.append("raw_or_processed_input_not_located")
            if source.path and source.path.startswith("starplast/data/"):
                cache = Path(paths.cache_file(source.path.removeprefix("starplast/data/")))
                if cache.is_file():
                    source_files = (SourceFile.inspect(cache, "installed_cache"),)
        if graph["gaps"]:
            gaps.append("upstream_source_graph_unresolved")
        upstream = tuple(key for key in graph["sources"] if key != source.key)
        for column in row["columns"] or ("source_records",):
            grade = "derived_quantity" if source.derived_from else "unresolved"
            unit, steps, mappings = "unresolved", (), ()
            column_gaps = gaps + ["quantity_unit_unresolved", "mapping_audit_unavailable", "transform_chain_unresolved"]
            deposit = next((entry for entry in deposits.DEPOSITS if entry.key == source.key), None)
            if deposit and source_files and source_files[0].role != "installed_cache":
                steps = (Transform("Registered deposit transformation; historical parameters not fully audited",
                                   "starplast.deposits." + deposit.derive.__name__, _sha(ROOT / "starplast/deposits.py"),
                                   tuple(item.path for item in source_files), (column,),
                                   {"parameters_status": "legacy_unresolved", "historical_run": "notebooks/derive_deposits_2026_09.ipynb"}),)
            if source.key == "host_bmdm_baseline" and column == "bmdm_tpm":
                grade, unit, source_files, steps, mappings = "direct_experiment", "TPM", files, pilot_steps, (mapping,)
                column_gaps = gaps + ["reference_release_label_unavailable; content pinned by mapping hashes"]
            traces.append(MeasurementTrace(source.key, source.organism, row["organism"], row["unit"], column,
                                           grade, unit, tuple(row["contexts"]), source_files, steps, mappings, upstream,
                                           metadata["publication_id"], metadata["publication_status"], gaps=tuple(column_gaps)))
    if {trace.source_id for trace in traces} != {source.key for source in D.REGISTRY}:
        raise ValueError("Provenance census lost registered sources")
    output.mkdir(parents=True)
    write_traces(traces, output / "traces.json")
    if len(read_traces(output / "traces.json")) != len(traces):
        raise ValueError("Trace round-trip failed")
    (output / "source_graphs.json").write_text(json.dumps(graphs, indent=2) + "\n")
    (output / "bmdm_pilot.json").write_text(json.dumps(pilot, indent=2) + "\n")
    summary = {"registered_sources": len(D.REGISTRY), "measurement_addresses": len(traces),
               "source_file_addresses": sum(bool(t.source_files) for t in traces),
               "mapping_audited_addresses": sum(bool(t.mappings) for t in traces),
               "unit_unresolved_addresses": sum(t.quantity_unit == "unresolved" for t in traces),
               "origin_paper_or_license_unresolved": "legacy gaps remain explicit; no permission inferred"}
    for path, digest in hashes.items():
        if _sha(path) != digest:
            raise ValueError("Input changed during provenance audit: " + path)
    manifest = {"schema_version": 1, "summary": summary, "input_sha256": hashes,
                "source_file_sha256": {item.path: item.sha256 for trace in traces for item in trace.source_files},
                "outputs": {name: _sha(output / name) for name in ("traces.json", "source_graphs.json", "bmdm_pilot.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main():
    """Execute the complete offline trace census and raw-to-installed pilot notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Measurement provenance census and raw host-expression reproduction")
    nb.md("Every inventory source/output has a qualified trace. Unknown units, mapping/version audits, "
          "license and origin-paper lineage remain explicit. Existing publication identities are reused, "
          "never invented. The BMDM pilot reproduces raw-file TPM and quantifies mapping ambiguity/loss.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from build_measurement_provenance import build",
            f"build(Path({str(args.root)!r}), Path({str(args.bindings)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "provenance.ipynb"))


if __name__ == "__main__":
    main()
