"""Download published processed workbooks for unresolved historical TSV derivations.

Primary PMC article JSON supplies exact URLs and MD5 checksums. These are review
inputs, not automatic bindings to differently named legacy outputs. Figure-data
archives, instrument RAW files and ambiguous article versions stay outside scope.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from recover_pmc_sources import named_media  # noqa: E402
from recover_source_files import download  # noqa: E402


def retrieve(root: Path, prior: Path, output: Path) -> dict:
    """Retrieve reviewed published media while retaining unmatched legacy filenames."""
    if output.exists():
        raise ValueError("Use a new workbook-review snapshot")
    output.mkdir(parents=True)
    inputs = [Path(__file__), prior / "source_locations.json"]
    records = []
    for row in json.loads((prior / "source_locations.json").read_text()):
        if row["status"] != "named_file_or_published_version_unresolved":
            continue
        metadata_files = list(prior.glob(row["pmcid"] + ".*.json"))
        metadata_files = [path for path in metadata_files if json.loads(path.read_text()).get("is_manuscript") is False]
        if len(metadata_files) != 1:
            records.append({"source_id": row["source_id"], "status": "published_article_version_unresolved"})
            continue
        metadata_file = metadata_files[0]
        inputs.append(metadata_file)
        metadata = json.loads(metadata_file.read_text())
        urls = metadata.get("media_urls", [])
        chosen = [url for url in urls if Path(urllib.parse.urlparse(url).path).suffix.lower() in {".xlsx", ".xls"}]
        # This named ZIP is a published supplementary container, not a figure-data archive.
        if row["source_id"] == "pf_isoforms":
            chosen += [url for url in urls if Path(urllib.parse.urlparse(url).path).name.endswith("_MOESM1_ESM.zip")]
        for original in chosen:
            filename = Path(urllib.parse.urlparse(original).path).name
            item = {"source_id": row["source_id"], "source_pmcid": row["pmcid"], "filename": filename,
                    "legacy_requested_filename": row["filename"], "recorded_utc": datetime.now(timezone.utc).isoformat()}
            try:
                matches = named_media(metadata, row["pmcid"], row["pmid"], filename)
                if len(matches) != 1:
                    raise ValueError("Named published review file is not unique")
                url = matches[0]
                target = root / "source_review_inputs" / row["source_id"] / filename
                download(url, target)
                payload = target.read_bytes()
                expected = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["md5"][0].lower()
                if hashlib.md5(payload).hexdigest() != expected:
                    raise ValueError("Published review-input MD5 mismatch")
                item.update(status="retrieved_verified_review_input", path=str(target), url=url,
                            bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                            license=metadata.get("license_code", "unresolved"),
                            association="same published article; legacy transform/table association requires review")
                target.with_name(filename + ".url.txt").write_text(url + "\n")
                target.with_name(filename + ".sha256.txt").write_text(item["sha256"] + "  " + filename + "\n")
            except Exception as error:
                item.update(status="retrieval_or_validation_failed", error=f"{type(error).__name__}: {error}")
            records.append(item)
            print(row["source_id"], filename, item["status"], item.get("error", ""), flush=True)
        if not chosen:
            records.append({"source_id": row["source_id"], "status": "no_workbook_in_published_media"})
    (output / "review_inputs.json").write_text(json.dumps(records, indent=2) + "\n")
    summary = {"processed_review_inputs": sum(item["status"] == "retrieved_verified_review_input" for item in records),
               "failures_or_gaps": sum(item["status"] != "retrieved_verified_review_input" for item in records),
               "new_legacy_output_bindings": 0}
    manifest = {"summary": summary, "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
                "outputs": {"review_inputs.json": hashlib.sha256((output / "review_inputs.json").read_bytes()).hexdigest()}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main():
    """Archive the bounded workbook acquisition with the exact published associations."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--prior", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Recover published processed review inputs for unmatched legacy derivations")
    nb.md("NIH NLM NCBI PubMed Central provides article metadata and MD5 checksums. The frozen scope "
          "is XLS/XLSX supplemental tables plus the named isoform supplemental container from unique "
          "published versions. These inputs retain original names and need transform/table association "
          "review. A workbook from the same paper is not automatically a replacement for a legacy TSV.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from recover_pmc_review_inputs import retrieve",
            f"retrieve(Path({str(args.root)!r}), Path({str(args.prior)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "retrieval.ipynb"))


if __name__ == "__main__":
    main()
