"""Retrieve one exact publisher table under its discovered original filename."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlparse

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from consolidate_source_recovery import direct_filename  # noqa: E402
from recover_source_files import download  # noqa: E402


def retrieve(source_key, root, headers_path, output):
    """Use exact registry URL and primary redirect filename; retain legacy mismatch."""
    root, headers_path, output = Path(root).resolve(), Path(headers_path).resolve(), Path(output)
    if not root.is_dir() or output.exists():
        raise ValueError("Use an existing archive and a new publisher-table output")
    source = D.get(source_key)
    metadata = json.loads(headers_path.read_text())
    if metadata["source_id"] != source.key or metadata["recorded_registry_url"] != source.url or metadata["status"] != 200:
        raise ValueError("Primary resource discovery does not match the declared source URL")
    name = direct_filename(metadata["final_url"])
    if not name or Path(name).suffix.lower() not in {".xlsx", ".xls"}:
        raise ValueError("No primary named spreadsheet; do not guess the legacy filename")
    target = root / "publisher_review_inputs" / source.key / name
    download(source.url, target)
    file = SourceFile.inspect(target, "processed_input", source.url,
                              "exact_registry_resource_URL_and_primary_redirect_filename; legacy_table_transform_unverified")
    workbook = pd.ExcelFile(target)
    sheets = {sheet: pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=5).fillna("").astype(str).values.tolist()
              for sheet in workbook.sheet_names}
    workbook.close()
    output.mkdir(parents=True)
    # Public signed redirect queries expire; primary provenance uses the stable
    # publisher URL. The exact original response is already pinned separately.
    published_path = urlparse(metadata["final_url"]).path
    record = {"source_id": source.key, "registry_organism": source.organism,
              "registry_legacy_path": source.path, "primary_resource_url": source.url,
              "discovered_published_path": published_path, "file": asdict(file), "sheet_headers": sheets,
              "status": "published_table_retrieved; legacy_filename_and_transform_association_requires_review",
              "runtime_changes": "none; no arbitrary renamed source or numerical promotion"}
    (target.parent / "URLS.txt").write_text(source.url + "\n")
    (target.parent / "SHA256SUMS.txt").write_text(file.sha256 + "  " + target.name + "\n")
    (output / "table_review.json").write_text(json.dumps(record, indent=2) + "\n")
    inputs = [Path(__file__), headers_path, ROOT / "starplast/datasets.py", ROOT / "scripts/recover_source_files.py",
              ROOT / "scripts/consolidate_source_recovery.py"]
    (output / "manifest.json").write_text(json.dumps({"input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
        "outputs": {"table_review.json": hashlib.sha256((output / "table_review.json").read_bytes()).hexdigest()},
        "source_sha256": file.sha256}, indent=2) + "\n")
    return {"source_id": source.key, "published_filename": name, "bytes": file.bytes,
            "sheets": list(sheets), "legacy_association": "requires_review"}


def main():
    """Execute one source-specific primary-table retrieval and preserve interpretation limits."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--headers", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    notebook = ExecutedNotebook("Primary publisher table retrieval without renaming a legacy input")
    notebook.md("Exact registry source URL, publisher-discovered filename, bounded validated download. "
                "Retain five header rows per sheet for context review; no assignment to a different legacy filename or runtime data change.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from retrieve_primary_table import retrieve",
                  f"retrieve({args.source!r},Path({str(args.root)!r}),Path({str(args.headers)!r}),Path({str(args.out)!r}))")
    notebook.md("The stable publisher resource and exact bytes establish acquisition. Actual cell/species, units, assay negatives and "
                "legacy transform association require table-specific review. A published filename mismatch is recorded, not papered over.")
    notebook.write(str(args.out / "retrieval.ipynb"))


if __name__ == "__main__":
    main()
