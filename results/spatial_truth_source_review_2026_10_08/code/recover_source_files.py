"""Recover registered processed inputs without overwriting data or guessing associations.

An explicit archive root is required. Relocated files are accepted only through
an acquisition receipt with the exact registry URL and a matching file hash.
Missing direct files and named members of recorded Europe PMC supplements are
retrieved sequentially, content-checked and recorded in an executed notebook.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import datasets as D, paths  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_file(path: Path):
    """Refuse web/error pages, corrupt archives and spreadsheets of the wrong type."""
    with path.open("rb") as stream:
        start = stream.read(1024).lstrip().lower()
    if not start or b"<html" in start or b"<!doctype html" in start or start.startswith(b"<error"):
        raise ValueError("Empty file or web/error page, not source data")
    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as stream:
            while stream.read(1 << 20):
                pass  # verify the complete compressed payload including CRC
    elif path.suffix.lower() == ".xlsx":
        with zipfile.ZipFile(path) as archive:
            if "xl/workbook.xml" not in archive.namelist() or archive.testzip() is not None:
                raise ValueError("Invalid Excel workbook")
    elif path.suffix.lower() == ".xls":
        import pandas as pd
        pd.read_excel(path, nrows=2)
    elif path.suffix.lower() == ".pdf" and not start.startswith(b"%pdf"):
        raise ValueError("Invalid PDF")
    elif path.suffix.lower() in {".csv", ".tsv", ".txt"}:
        start.decode("utf-8")


def receipt_binding(source, root: Path, receipts) -> tuple[SourceFile | None, str]:
    """Match the recorded URL, source organism, filename and pinned receipt hash."""
    name = Path(source.path or "").name
    choices = {}
    for receipt_path, receipt in receipts:
        if source.url not in receipt.get("urls", []) or receipt.get("space") != source.organism:
            continue
        for item in receipt.get("files", []):
            if item.get("name") != name:
                continue
            candidate = (root / item["path"]).resolve()
            if not candidate.is_relative_to(root.resolve()) or not candidate.is_file():
                continue
            record = SourceFile.inspect(candidate, "processed_input", source.url,
                                        "exact_registry_url_and_acquisition_receipt_hash")
            if record.sha256 != item.get("sha256"):
                raise ValueError("Acquisition receipt hash mismatch: " + str(candidate))
            choices[record.path] = record
    return (next(iter(choices.values())), "receipt_verified") if len(choices) == 1 else \
           (None, "ambiguous_receipt_bindings" if choices else "not_located")


def download(url: str, target: Path, limit=64 << 20):
    """Stream a bounded download into a temporary file; commit only verified content."""
    if target.exists():
        raise ValueError("Refusing to overwrite source file")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(".partial-" + target.name)
    if temporary.exists():
        raise ValueError("A previous partial download needs separate inspection")
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "starplast-source-recovery/1"})
        with urllib.request.urlopen(request, timeout=25) as response, temporary.open("xb") as stream:
            total = 0
            while block := response.read(1 << 20):
                total += len(block)
                if total > limit:
                    raise ValueError("Processed-source download exceeds the 64 MiB per-file scope")
                stream.write(block)
        # A supplemental container is separately validated before named extraction.
        if target.name.endswith(".zip"):
            with zipfile.ZipFile(temporary) as archive:
                if archive.testzip() is not None:
                    raise ValueError("Corrupt supplementary archive")
        else:
            # Validation uses the real filename suffix, without mutating the source.
            validate_file(temporary)
        temporary.rename(target)
    finally:
        if temporary.exists():
            temporary.unlink()


def recover(root: Path, output: Path, keys=None, fetch=False) -> dict:
    """Record every selected source, retrieve only named processed files, retain gaps."""
    root = root.resolve()
    if not root.is_dir() or output.exists():
        raise ValueError("Use an existing explicit archive root and a new result directory")
    output.mkdir(parents=True)
    receipts = []
    receipt_dir = root / "spaces/_acquisition/parts"
    if receipt_dir.is_dir():
        for receipt_path in sorted(receipt_dir.glob("*/*.json")):
            receipts.append((receipt_path, json.loads(receipt_path.read_text())))
    chosen = [d for d in D.REGISTRY if keys is None or d.key in keys]
    if keys is not None and set(keys) != {d.key for d in chosen}:
        raise ValueError("Unknown source key")
    records, bindings = [], {}
    for source in chosen:
        row = {"source_id": source.key, "registry_path": source.path, "url": source.url,
               "organism": source.organism, "status": "not_located", "error": "", "file": None}
        try:
            if not source.path:
                row["status"] = "no_path_declared"
            elif source.path.startswith("starplast/data/"):
                file = Path(paths.cache_file(source.path.removeprefix("starplast/data/")))
                if file.is_file():
                    row.update(status="installed_cache", file=asdict(SourceFile.inspect(file, "installed_cache")))
            else:
                relative = source.path.removeprefix("datasets/")
                declared = (root / relative).resolve()
                if not declared.is_relative_to(root):
                    row["status"] = "outside_declared_archive"
                elif declared.is_file():
                    row.update(status="present", file=asdict(SourceFile.inspect(declared, "processed_input", source.url or "")))
                elif declared.is_dir():
                    row["status"] = "source_container_requires_member_manifest"
                else:
                    bound, row["status"] = receipt_binding(source, root, receipts)
                    if bound:
                        row["file"] = asdict(bound)
                    elif fetch and source.url and row["status"] != "ambiguous_receipt_bindings":
                        url_name = Path(urllib.parse.urlparse(source.url).path).name
                        suffixes = (".xlsx", ".xls", ".csv", ".tsv", ".txt", ".gz", ".pdf")
                        if source.url.endswith("/supplementaryFiles") and declared.name.endswith(suffixes):
                            container = output / "supplements" / (source.key + ".zip")
                            download(source.url, container)
                            with zipfile.ZipFile(container) as archive:
                                matches = [item for item in archive.infolist() if Path(item.filename).name == declared.name]
                                if len(matches) != 1 or matches[0].file_size > 64 << 20:
                                    raise ValueError("Named source member unavailable, ambiguous or outside size scope")
                                declared.parent.mkdir(parents=True, exist_ok=True)
                                temporary = declared.with_name(".partial-" + declared.name)
                                if temporary.exists():
                                    raise ValueError("Partial source member exists")
                                try:
                                    with archive.open(matches[0]) as src, temporary.open("xb") as dst:
                                        while block := src.read(1 << 20):
                                            dst.write(block)
                                    validate_file(temporary)
                                    temporary.rename(declared)
                                finally:
                                    if temporary.exists():
                                        temporary.unlink()
                            row.update(status="retrieved_named_supplement", file=asdict(SourceFile.inspect(declared, "processed_input", source.url)))
                        elif declared.name == url_name and url_name.endswith(suffixes) and "veupathdb" not in source.url and "toxodb.org" not in source.url and "plasmodb.org" not in source.url:
                            download(source.url, declared)
                            row.update(status="retrieved_direct", file=asdict(SourceFile.inspect(declared, "processed_input", source.url)))
                        else:
                            row["status"] = "needs_repository_or_source_review"
            if row["file"] and row["file"]["role"] != "installed_cache":
                bindings[source.key] = row["file"]
        except Exception as error:
            row.update(status="retrieval_or_identity_failed", error=f"{type(error).__name__}: {error}")
        records.append(row)
        print(source.key, row["status"], row["error"], flush=True)
    (output / "source_locations.json").write_text(json.dumps(records, indent=2) + "\n")
    (output / "bindings.json").write_text(json.dumps(bindings, indent=2) + "\n")
    summary = {"sources": len(records), "bound_processed_sources": len(bindings),
               "statuses": {status: sum(r["status"] == status for r in records) for status in sorted({r["status"] for r in records})}}
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "archive_root": str(root),
                "summary": summary, "input_sha256": {str(Path(__file__)): _sha(__file__),
                str(ROOT / "starplast/datasets.py"): _sha(ROOT / "starplast/datasets.py")},
                "receipt_sha256": {str(path): _sha(path) for path, _ in receipts},
                "outputs": {name: _sha(output / name) for name in ("source_locations.json", "bindings.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main():
    """Run explicitly scoped recovery with an executed, annotated download notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--keys", nargs="+")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    notebook = ExecutedNotebook("Registered source locations, receipt verification and bounded downloads")
    notebook.md("Explicit archive root; no scientific data promotion. Exact receipt URL/species/name/hash "
                "are required to accept a relocated source. Downloads are sequential, capped at 64 MiB per file, "
                "content checked and never overwrite an existing input. Missing cases are retained.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from recover_source_files import recover",
                  f"recover(Path({str(args.root)!r}), Path({str(args.out)!r}), keys={args.keys!r}, fetch={args.fetch!r})")
    notebook.write(str(args.out / "recovery.ipynb"))


if __name__ == "__main__":
    main()
