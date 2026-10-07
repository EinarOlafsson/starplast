"""Recover named processed supplements via the current public PMC Cloud Service.

The legacy archive is retired. Article versions are discovered, never guessed;
the JSON metadata must match the source PMCID and PMID. Only an unambiguous
published article version and exact named media file can bind a source. Every
failed retrieval, unavailable version and missing filename remains visible.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import datasets as D  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from recover_source_files import download  # noqa: E402

BASE = "https://pmc-oa-opendata.s3.amazonaws.com/"


def _metadata(url, target):
    request = urllib.request.Request(url, headers={"User-Agent": "starplast-source-recovery/1"})
    with urllib.request.urlopen(request, timeout=25) as response:
        payload = response.read((4 << 20) + 1)
    if len(payload) > 4 << 20:
        raise ValueError("Metadata response exceeds scope")
    target.write_bytes(payload)
    return payload


def named_media(metadata, pmcid, pmid, filename):
    """Match exact source identity and filename without choosing between versions."""
    if str(metadata.get("pmcid")) != pmcid or str(metadata.get("pmid")) != str(pmid):
        raise ValueError("PMC metadata differs from declared PMCID/PMID")
    if str(metadata.get("is_manuscript", "")).lower() != "no":
        return []
    if str(metadata.get("is_retracted", "")).lower() != "no":
        raise ValueError("Retracted or unclassified article cannot recover a source")
    matches = []
    for url in metadata.get("media_urls", []):
        if not isinstance(url, str):
            raise ValueError("Unexpected PMC media URL type")
        parsed = urllib.parse.urlparse(url)
        if parsed.hostname != "pmc-oa-opendata.s3.amazonaws.com" or parsed.scheme != "https":
            raise ValueError("Unexpected PMC media URL origin")
        if Path(urllib.parse.unquote(parsed.path)).name == filename:
            digest = urllib.parse.parse_qs(parsed.query).get("md5", [])
            if len(digest) != 1 or not re.fullmatch(r"[a-fA-F0-9]{32}", digest[0]):
                raise ValueError("PMC media file lacks a published MD5")
            matches.append(url)
    return matches


def recover(root: Path, prior: Path, output: Path) -> dict:
    """Use primary discovered article metadata to retrieve all missing named PMC inputs."""
    if not root.is_dir() or output.exists():
        raise ValueError("Use an existing archive root and a new recovery snapshot")
    output.mkdir(parents=True)
    docs = _metadata(BASE + "README.txt", output / "PMC_README.txt")
    bindings_path = prior / "bindings.json"
    bindings = json.loads(bindings_path.read_text())
    prior_locations = json.loads((prior / "source_locations.json").read_text())
    records, cache = [], {}
    for prior_row in prior_locations:
        source = D.get(prior_row["source_id"])
        if source.key in bindings or not source.path:
            continue
        match = re.search(r"/(PMC\d+)/supplementaryFiles$", source.url or "")
        filename = Path(source.path).name
        if not match or not filename.lower().endswith((".xlsx", ".xls", ".tsv", ".csv", ".pdf", ".txt")):
            continue
        pmcid = match[1]
        row = {"source_id": source.key, "pmcid": pmcid, "pmid": source.pmid, "filename": filename,
               "previous_status": prior_row["status"], "recorded_utc": datetime.now(timezone.utc).isoformat()}
        try:
            if pmcid not in cache:
                query = urllib.parse.urlencode({"list-type": "2", "prefix": pmcid + ".", "delimiter": "/"})
                url = BASE + "?" + query
                payload = _metadata(url, output / (pmcid + "_versions.xml"))
                tree = ET.fromstring(payload)
                if tree.findtext("{*}IsTruncated") != "false":
                    raise ValueError("PMC version listing is truncated or invalid")
                prefixes = [element.text for element in tree.findall("{*}CommonPrefixes/{*}Prefix")]
                versions = []
                for prefix in prefixes:
                    if not re.fullmatch(pmcid + r"\.\d+/", prefix or ""):
                        raise ValueError("Unexpected article prefix")
                    metadata_url = BASE + prefix + prefix.rstrip("/") + ".json"
                    payload = _metadata(metadata_url, output / (prefix.rstrip("/") + ".json"))
                    versions.append(json.loads(payload))
                cache[pmcid] = versions
            candidates = [(metadata, url) for metadata in cache[pmcid]
                          for url in named_media(metadata, pmcid, source.pmid, filename)]
            if len(candidates) != 1:
                row.update(status="named_file_or_published_version_unresolved", candidates=len(candidates),
                           versions=len(cache[pmcid]))
            else:
                metadata, url = candidates[0]
                target = root / "recovered_processed_sources" / source.key / filename
                download(url, target)
                payload = target.read_bytes()
                expected = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["md5"][0].lower()
                if hashlib.md5(payload).hexdigest() != expected:
                    raise ValueError("PMC supplemental checksum differs from published metadata")
                file = SourceFile.inspect(target, "processed_input", url, "PMC_article_PMCID_PMID_exact_filename_and_MD5_verified")
                bindings[source.key] = asdict(file)
                row.update(status="recovered_exact_named_published_supplement", file=asdict(file),
                           license=metadata.get("license_code", "unresolved"),
                           article_version=metadata.get("version"), article_doi=metadata.get("doi"))
                (target.parent / "URLS.txt").write_text(url + "\n")
                (target.parent / "SHA256SUMS.txt").write_text(file.sha256 + "  " + filename + "\n")
        except Exception as error:
            row.update(status="metadata_or_retrieval_failed", error=f"{type(error).__name__}: {error}")
        records.append(row)
        print(source.key, row["status"], row.get("error", ""), flush=True)
    (output / "source_locations.json").write_text(json.dumps(records, indent=2) + "\n")
    (output / "bindings.json").write_text(json.dumps(bindings, indent=2) + "\n")
    summary = {"attempted_named_missing_sources": len(records), "bound_processed_sources": len(bindings),
               "recovered_sources": sum(r["status"].startswith("recovered_exact") for r in records)}
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "summary": summary,
                "primary_service": "NIH NLM NCBI PubMed Central Article Datasets on AWS",
                "service_docs": "https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/",
                "input_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in
                                 (Path(__file__), bindings_path, prior / "source_locations.json")},
                "outputs": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in output.iterdir()}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main():
    """Execute the recovery with archived primary service documentation and metadata."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--prior", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Recover missing named supplements with the current PMC Cloud Service")
    nb.md("NLM retired the legacy archive in August 2026. This run acknowledges NIH NLM NCBI "
          "PubMed Central as the service source, saves current documentation, lists available article "
          "versions, checks PMCID/PMID, manuscript/retraction state and exact media filename, then "
          "verifies the MD5 provided in primary JSON metadata. No higher version is automatically "
          "preferred. Files remain external; article-level license does not authorize wheel redistribution.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from recover_pmc_sources import recover",
            f"recover(Path({str(args.root)!r}), Path({str(args.prior)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "recovery.ipynb"))


if __name__ == "__main__":
    main()
