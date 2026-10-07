"""Freeze an offline evidence inventory and its executed reconciliation notebook.

Run with --out pointing at a new results directory. No downloads or table writes
occur; all registry records and question-specific refusals remain visible.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D, organisms as O, paths, slots  # noqa: E402
from starplast.inventory import build_inventory  # noqa: E402


def _hash(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_snapshot(output: Path) -> pd.DataFrame:
    """Read installed references and persist a new inventory with stable input hashes."""
    from generate_slot_table import REFUSED_CANDIDATES

    if output.exists():
        raise ValueError("Use a new output directory; existing snapshots are immutable")
    inputs = [ROOT / "starplast" / filename for filename in
              ("datasets.py", "inventory.py", "organisms.py", "slots.py", "evidence.py", "data/slots.json")]
    inputs += [Path(__file__), ROOT / "scripts/generate_slot_table.py"]
    tables, bridges, graph_paths, row_counts = {}, {}, {}, {}
    table_paths = {}
    for code in O.codes():
        table_paths[(code, "gene")] = Path(O.nodes_path(code))
        graph_paths[code] = Path(O.graph_path(code))
        bridge = O.get(code).host_bridges
        if bridge:
            path = Path(paths.cache_file(bridge))
            if path.exists():
                bridges[code] = pd.read_parquet(path)
                inputs.append(path)
                row_counts[bridge] = len(bridges[code])
    for code, filename in O.HOST_TABLES.items():
        table_paths[(code, "protein")] = Path(paths.cache_file(filename))
    for key, path in table_paths.items():
        if path.exists():
            tables[key] = pd.read_parquet(path)
            inputs.append(path)
            row_counts[path.name] = len(tables[key])
    metabolites = Path(paths.cache_file("metabolites.parquet"))
    if metabolites.exists():
        frame = pd.read_parquet(metabolites)
        # This is the installed metabolite row universe, not a claim that every
        # row was assayed in both species. Each source selects its own columns.
        for code in O.codes():
            tables[(code, "metabolite")] = frame
        inputs.append(metabolites)
        row_counts[metabolites.name] = len(frame)
    raw_paths = {d.key: D.local_path(d.key) for d in D.REGISTRY}
    slot_map = {s.organism + "_" + s.name: s for s in slots.all_slots()}
    refusals = []
    for (question, accession), reason in sorted(REFUSED_CANDIDATES.items()):
        slot = slot_map[question]
        refusals.append({"source_id": "candidate:" + accession, "organism": slot.organism,
                         "unit": slot.unit, "question": question, "reason": reason})
    inputs += [path for path in graph_paths.values() if path.is_file()]
    hashes = {str(path): _hash(path) for path in sorted(set(inputs))}
    with ExitStack() as stack:
        graphs = {code: stack.enter_context(np.load(path)) for code, path in graph_paths.items() if path.is_file()}
        report = build_inventory(tables, graphs, bridges, raw_paths=raw_paths, refusals=refusals)
    registry_rows = report[report.status != "rejected"]
    if set(registry_rows.source_id) != {d.key for d in D.REGISTRY}:
        raise ValueError("Inventory does not reconcile with every registered source")
    if registry_rows.duplicated(["source_id", "organism", "unit"]).any():
        raise ValueError("Duplicate source/organism/storage-unit address")
    if len(report[report.status == "rejected"]) != len(REFUSED_CANDIDATES):
        raise ValueError("Refusal count changed")
    for path, digest in hashes.items():
        if _hash(path) != digest:
            raise ValueError(f"Input changed during the census: {path}")
    output.mkdir(parents=True)
    report.to_parquet(output / "inventory.parquet", index=False)
    (output / "inventory.json").write_text(report.to_json(orient="records", indent=2) + "\n")
    (output / "registry.json").write_text(json.dumps([asdict(d) for d in D.REGISTRY], indent=2) + "\n")
    manifest = {"schema_version": 1, "base_commit": __import__("subprocess").check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "code_identity": "Exact worktree files are identified by input_sha256; base_commit is the parent baseline",
        "software": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__},
        "input_sha256": hashes, "registry_sources": len(D.REGISTRY), "inventory_rows": len(report),
        "registry_rows": len(registry_rows), "refused_questions": len(refusals),
        "installed_table_rows": row_counts, "remote_access": "not_checked",
        "source_observations": "not_supplied; missingness causes remain unknown",
        "outputs": {name: _hash(output / name) for name in ("inventory.parquet", "inventory.json", "registry.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"sources": len(D.REGISTRY), "rows": len(report), "refusals": len(refusals),
                      "statuses": report.status.value_counts().to_dict(), "table_rows": row_counts}, indent=2))
    return report


def main():
    """Execute the census in a notebook and refuse to overwrite an earlier snapshot."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Registered evidence inventory and installed-table reconciliation")
    notebook.md("# Offline information-space census",
                "Frozen scope: current dataset registry, installed parasite/host reference tables, "
                "declared graph/bridge outputs and retained question-specific refusals. No acquisition "
                "or model fitting. Stored coverage is not biological accuracy or assay reach.")
    notebook.code("from pathlib import Path", "import sys",
                  f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from build_information_inventory import write_snapshot",
                  f"report = write_snapshot(Path({str(args.out.resolve())!r}))",
                  "report.groupby(['organism', 'unit', 'status']).size().rename('rows')")
    notebook.md("## Interpretation limits",
                "Registry organism and output-table organism are retained separately. Host tables "
                "contain protein references, not admitted host gene spaces. Unmeasured and unmapped "
                "counts are unavailable without source-level observations. False and numeric zero "
                "remain stored values; only a measured boolean-negative Observation certifies an "
                "explicit negative. Pair records have no assayed-universe denominator. Bridge rows "
                "without source IDs remain unattributed. Rejections apply only to the named question. "
                "Remote source access was not probed. Reproduce with a new --out directory.")
    notebook.write(str(args.out / "census.ipynb"))


if __name__ == "__main__":
    main()
