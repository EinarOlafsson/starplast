"""Reconcile the complete source queue, project inputs and bounded archive members.

Preserve every failed or unknown association. Direct processed-file URLs can be
retrieved sequentially; no restricted service, raw instrument archive, template
URL or mismatched PMC identity is queried. Runtime paths and data are unchanged.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D, paths  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from recover_source_files import download, validate_file  # noqa: E402

PROCESSED = {".xlsx", ".xls", ".csv", ".tsv", ".txt", ".json", ".jsonl", ".bed"}
RESTRICTED = {"toxodb.org", "plasmodb.org", "cryptodb.org", "vectorbase.org", "veupathdb.org"}


def direct_filename(url):
    """Return only an explicitly named processed resource on a nonrestricted host."""
    if not url or "{" in url or "}" in url:
        return None
    parsed = urlparse(url)
    if parsed.scheme != "https" or any(parsed.hostname == host or (parsed.hostname or "").endswith("." + host) for host in RESTRICTED):
        return None
    name = unquote(Path(parsed.path).name)
    if not name or Path(name).name != name or name in {".", ".."}:
        return None
    plain = name[:-3] if name.endswith(".gz") else name
    if Path(plain).suffix.lower() not in PROCESSED:
        return None
    return name


def member_manifest(directory, allowed_root, *, max_members=500, max_bytes=256 << 20):
    """Hash bounded existing processed members without following external symlinks.

    Membership is a location observation, not assignment of every file to the
    directory's source. Instrument RAW and sequence files are excluded. A cap
    leaves an explicit partial manifest rather than claiming complete recovery.
    """
    directory, allowed_root = Path(directory).resolve(), Path(allowed_root).resolve()
    if not directory.is_relative_to(allowed_root) or not directory.is_dir():
        raise ValueError("Container must be inside the explicitly allowed archive")
    files, total, visited, skipped = [], 0, 0, []
    complete = True
    for current, dirs, names in os.walk(directory, followlinks=False):
        dirs[:] = sorted(name for name in dirs if not (Path(current) / name).is_symlink())
        for name in sorted(names):
            visited += 1
            if visited > max_members * 10:
                complete = False
                break
            file = Path(current) / name
            if file.is_symlink():
                skipped.append({"path": str(file), "reason": "symlink_member_requires_receipt_review"})
                continue
            plain = name[:-3] if name.endswith(".gz") else name
            if Path(plain).suffix.lower() not in PROCESSED or name in {"URLS.txt", "SHA256SUMS.txt"}:
                continue
            size = file.stat().st_size
            if len(files) >= max_members or total + size > max_bytes:
                complete = False
                skipped.append({"path": str(file), "reason": "member_or_byte_cap; not_hashed"})
                continue
            record = SourceFile.inspect(file, "processed_input", association="existing_container_member; biological_source_assignment_unverified")
            files.append(asdict(record))
            total += record.bytes
        if visited > max_members * 10:
            break
    return {"path": str(directory), "members": files, "hashed_bytes": total,
            "complete_within_processed_scope": complete, "visited_files": visited, "skipped": skipped,
            "limits": {"members": max_members, "bytes": max_bytes},
            "scope": "existing processed candidates only; not raw instrument files or proof of article association"}


def project_input(source, root):
    """Locate only the exact registered sibling-project path; never search by basename."""
    if not source.path or not source.path.startswith("toxo_stage_atlas/"):
        return None
    candidate = (root.parent / source.path).resolve()
    project = (root.parent / "toxo_stage_atlas").resolve()
    if not candidate.is_relative_to(project):
        raise ValueError("Project-relative source escapes its declared project")
    return candidate if candidate.exists() else None


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def consolidate(root, output, *, fetch=False):
    """Account for all sources, merge verified recovery and retain every remaining gap."""
    root, output = Path(root).resolve(), Path(output).resolve()
    if not root.is_dir() or output.exists():
        raise ValueError("Use an existing explicit archive root and new result directory")
    prior = ROOT / "results/source_recovery_pmc_2026_10_07_v2"
    initial = ROOT / "results/source_recovery_2026_10_07"
    review_path = ROOT / "results/source_review_input_recovery_2026_10_07/review_inputs.json"
    inputs = [Path(__file__), ROOT / "scripts/recover_source_files.py", ROOT / "starplast/datasets.py",
              prior / "bindings.json", prior / "source_locations.json", initial / "source_locations.json", review_path]
    hashes = {str(path): _sha(path) for path in inputs}
    bindings = json.loads((prior / "bindings.json").read_text())
    for value in bindings.values():
        SourceFile(**value).verify()
    previous = {row["source_id"]: row for row in json.loads((initial / "source_locations.json").read_text())}
    cloud = {row["source_id"]: row for row in json.loads((prior / "source_locations.json").read_text())}
    reviews = json.loads(review_path.read_text())
    output.mkdir(parents=True)
    containers, records, failures = {}, [], []
    for source in D.REGISTRY:
        row = {"source_id": source.key, "registry_path": source.path, "registry_url": source.url,
               "registry_organism": source.organism, "initial_status": previous[source.key]["status"],
               "status": "not_located", "error": "", "cloud_review": cloud.get(source.key),
               "review_inputs": [item for item in reviews if item["source_id"] == source.key],
               "source_publication_association": "requires_source_specific_review"}
        try:
            if source.key in bindings:
                row.update(status="bound_processed_input", file=bindings[source.key])
            elif source.path and source.path.startswith("starplast/data/"):
                file = Path(paths.cache_file(source.path.removeprefix("starplast/data/")))
                row.update(status="installed_cache_not_original_assay", file=asdict(SourceFile.inspect(file, "installed_cache")))
            else:
                project = project_input(source, root)
                declared = (root / source.path.removeprefix("datasets/")).resolve() if source.path and not source.path.startswith("/") else None
                if project and project.is_file():
                    validate_file(project)
                    record = SourceFile.inspect(project, "processed_input", source.url or "", "exact_registry_sibling_project_path; article_file_association_unverified")
                    bindings[source.key] = asdict(record)
                    row.update(status="located_project_input", file=asdict(record))
                elif project and project.is_dir():
                    manifest = member_manifest(project, root.parent / "toxo_stage_atlas")
                    containers[source.key] = manifest
                    row.update(status="project_container_manifest", container=manifest)
                elif declared and declared.is_relative_to(root) and declared.is_dir():
                    manifest = member_manifest(declared, root)
                    containers[source.key] = manifest
                    row.update(status="archive_container_manifest", container=manifest)
                elif declared and declared.is_relative_to(root) and declared.is_file():
                    validate_file(declared)
                    record = SourceFile.inspect(declared, "processed_input", source.url or "", "exact_registry_archive_path; article_file_association_unverified")
                    bindings[source.key] = asdict(record)
                    row.update(status="located_archive_input", file=asdict(record))
                else:
                    name = direct_filename(source.url)
                    if cloud.get(source.key, {}).get("status") == "article_identity_mismatch":
                        row["status"] = "refused_recorded_article_identity_mismatch"
                    elif fetch and name:
                        # The new recovery location never mutates a stale legacy link.
                        target = root / "recovered_direct_sources" / source.key / name
                        download(source.url, target)
                        record = SourceFile.inspect(target, "processed_input", source.url,
                                                    "exact_registry_primary_resource_URL; article_and_transform_unverified")
                        (target.parent / "URLS.txt").write_text(source.url + "\n")
                        (target.parent / "SHA256SUMS.txt").write_text(record.sha256 + "  " + target.name + "\n")
                        bindings[source.key] = asdict(record)
                        row.update(status="retrieved_direct_processed_input", file=asdict(record))
                    elif source.url and any(host in (urlparse(source.url).hostname or "") for host in RESTRICTED):
                        row["status"] = "restricted_service_not_queried"
                    elif source.derived_from or not source.url:
                        row["status"] = "derived_or_unspecified_input_requires_transform_review"
                    elif row["review_inputs"]:
                        row["status"] = "published_review_inputs_available; legacy_transform_unresolved"
                    else:
                        row["status"] = "source_deposit_or_receipt_review_required"
        except Exception as error:
            row.update(status="retrieval_or_location_validation_failed", error=f"{type(error).__name__}: {error}")
            failures.append({"source_id": source.key, "error": row["error"]})
        records.append(row)
        print(source.key, row["status"], row["error"], flush=True)
    if {row["source_id"] for row in records} != {source.key for source in D.REGISTRY}:
        raise ValueError("Source reconciliation lost a registered source")
    for name, value in (("source_locations.json", records), ("bindings.json", bindings), ("containers.json", containers)):
        (output / name).write_text(json.dumps(value, indent=2) + "\n")
    summary = {"registered_sources": len(records), "processed_source_bindings": len(bindings), "container_manifests": len(containers),
               "container_candidate_files": sum(len(value["members"]) for value in containers.values()),
               "partial_containers": sum(not value["complete_within_processed_scope"] for value in containers.values()),
               "statuses": dict(__import__("collections").Counter(row["status"] for row in records)),
               "failures": failures, "runtime_changes": "none; numerical promotion and transform admission remain pending"}
    for path, digest in hashes.items():
        if _sha(path) != digest:
            raise ValueError("Source input changed during reconciliation: " + path)
    (output / "manifest.json").write_text(json.dumps({"created_utc": datetime.now(timezone.utc).isoformat(), "archive_root": str(root),
        "input_sha256": hashes, "summary": summary, "outputs": {name: _sha(output / name) for name in
        ("source_locations.json", "bindings.json", "containers.json")}}, indent=2) + "\n")
    return summary


def main():
    """Execute location reconciliation and sequential downloads in a retained notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    notebook = ExecutedNotebook("Complete source reconciliation, bounded members and missing processed inputs")
    notebook.md("Frozen scope: all 162 registered sources, existing verified recovery bindings, exact sibling-project paths and "
                "bounded processed-container manifests. Direct processed URLs are retrieved sequentially; no raw instrument archives, "
                "restricted services or guessed paper/file identities. Keep missing content as a gap and continue the queue.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from consolidate_source_recovery import consolidate",
                  f"consolidate(Path({str(args.root)!r}), Path({str(args.out)!r}), fetch={args.fetch!r})")
    notebook.md("Located/retrieved inputs are not admission of a new dataset. A resolved PMID/DOI does not establish that a file generated "
                "a legacy matrix. Container members and 43 published review inputs require source/table/transform assignment. Installed "
                "caches remain distinct from original assays. Failed and unavailable inputs remain recorded; runtime data and paths are unchanged.")
    notebook.write(str(args.out / "reconciliation.ipynb"))


if __name__ == "__main__":
    main()
