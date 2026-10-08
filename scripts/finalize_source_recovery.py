"""Join verified publisher equivalence to the complete source-location census."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import parse_qs, urlparse

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finalize(prior, publisher_review, output):
    """Accept a publisher filename alias only with paper identity and exact cache reproduction."""
    prior, publisher_review, output = Path(prior), Path(publisher_review), Path(output)
    if output.exists():
        raise ValueError("Use a new final source-census directory")
    previous_manifest = json.loads((prior / "manifest.json").read_text())
    for name, digest in previous_manifest["outputs"].items():
        if _sha(prior / name) != digest:
            raise ValueError("Prior source census output changed: " + name)
    rows = json.loads((prior / "source_locations.json").read_text())
    bindings = json.loads((prior / "bindings.json").read_text())
    for row in rows:
        source = D.get(row["source_id"])
        if (row["registry_path"], row["registry_url"], row["registry_organism"]) != (source.path, source.url, source.organism):
            raise ValueError("Source identity changed; do not reuse historical associations")
    for value in bindings.values():
        SourceFile(**value).verify()
    for manifest_name in ("manifest.json", "equivalence_manifest.json"):
        publisher_manifest = json.loads((publisher_review / manifest_name).read_text())
        for name, digest in publisher_manifest["outputs"].items():
            if _sha(publisher_review / name) != digest:
                raise ValueError("Publisher review evidence changed: " + name)
        if manifest_name == "equivalence_manifest.json":
            for path, digest in publisher_manifest["input_sha256"].items():
                if _sha(path) != digest:
                    raise ValueError("Publisher equivalence input changed: " + path)
    equivalence = json.loads((publisher_review / "equivalence.json").read_text())
    table_review = json.loads((publisher_review / "table_review.json").read_text())
    source = D.get(equivalence["source_id"])
    if table_review["source_id"] != source.key or table_review["primary_resource_url"] != source.url:
        raise ValueError("Publisher resource does not match the declared source")
    publications = ROOT / "results/dataset_selection_2026_10_07/source_publications.csv"
    publication = pd.read_csv(publications, dtype=str, keep_default_na=False).set_index("dataset_id").loc[source.key]
    resource_id = parse_qs(urlparse(source.url).query).get("id", [])
    if publication.pmid != source.pmid or not publication.doi or len(resource_id) != 1 or not resource_id[0].startswith(publication.doi + "."):
        raise ValueError("Publisher resource DOI must match the previously verified PMID identity")
    if set(equivalence["installed_reproduction"]) != set(source.columns):
        raise ValueError("Publisher alias must reproduce every installed source field")
    for check in equivalence["installed_reproduction"].values():
        if check["lost_ids"] or check["max_abs_difference"] != 0 or check["installed_nonmissing"] <= 0:
            raise ValueError("Publisher filename alias is not exact or loses installed values")
    if equivalence["file"] != table_review["file"]:
        raise ValueError("Equivalence and publisher review used different input bytes")
    file = SourceFile(**equivalence["file"])
    file.verify()
    file = replace(file, association="verified_publisher_DOI_to_PMID_and_exact_installed_CSPA_fields; original_published_filename_preserved")
    bindings[source.key] = asdict(file)
    row = next(row for row in rows if row["source_id"] == source.key)
    row.update(status="bound_publisher_filename_alias_with_exact_reproduction", file=asdict(file),
               source_publication_association="verified_publisher_DOI_to_previously_resolved_PMID",
               source_observation_semantics=equivalence["negative_semantics"],
               publisher_equivalence=equivalence["installed_reproduction"],
               original_published_filename=equivalence["published_filename"],
               legacy_filename=equivalence["legacy_filename"])
    if len(rows) != len(D.REGISTRY) or {row["source_id"] for row in rows} != {source.key for source in D.REGISTRY}:
        raise ValueError("Final recovery census lost registered source addresses")
    output.mkdir(parents=True)
    for name, value in (("source_locations.json", rows), ("bindings.json", bindings),
                        ("containers.json", json.loads((prior / "containers.json").read_text()))):
        (output / name).write_text(json.dumps(value, indent=2) + "\n")
    summary = {**previous_manifest["summary"], "processed_source_bindings": len(bindings),
               "statuses": dict(Counter(row["status"] for row in rows)),
               "publisher_alias_verified": source.key, "numerical_changes": "none",
               "metadata_correction": "capture non-detection is not proof of biological absence; mouse references are separate from human"}
    inputs = [Path(__file__), prior / "manifest.json", prior / "source_locations.json", prior / "bindings.json", prior / "containers.json",
              publisher_review / "equivalence.json", publisher_review / "table_review.json", publisher_review / "manifest.json",
              publisher_review / "equivalence_manifest.json",
              publications, ROOT / "starplast/datasets.py", ROOT / "starplast/host.py"]
    (output / "manifest.json").write_text(json.dumps({"created_utc": datetime.now(timezone.utc).isoformat(), "summary": summary,
        "input_sha256": {str(path): _sha(path) for path in inputs}, "source_alias_sha256": file.sha256,
        "outputs": {name: _sha(output / name) for name in ("source_locations.json", "bindings.json", "containers.json")}}, indent=2) + "\n")
    return summary


def main():
    """Execute the final source ledger without changing numerical tables or legacy paths."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--publisher-review", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Final source recovery ledger and verified publisher filename equivalence")
    notebook.md("Frozen scope: merge the complete census with one publisher source whose original filename differs from the legacy path. "
                "Require resolved PMID/DOI association, exact input identity, every installed field reproduced and no lost values. "
                "No file rename, numerical dataset change, restricted source query or new biological absence claim.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from finalize_source_recovery import finalize",
                  f"finalize(Path({str(args.prior)!r}),Path({str(args.publisher_review)!r}),Path({str(args.out)!r}))")
    notebook.md("All source/location gaps, unknown table transforms, installed caches and the PMID/PMCID refusal remain visible. "
                "Container membership is not biological admission; remaining source-specific admission and replacements belong to 65.04. "
                "The source census/recovery item can complete while unavailable evidence remains explicitly unresolved.")
    notebook.write(str(args.out / "finalization.ipynb"))


if __name__ == "__main__":
    main()
