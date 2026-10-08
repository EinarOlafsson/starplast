"""Retrieve public host candidate evidence from verified deposit metadata.

Downloads remain in the external source archive; the executed notebook retains
primary metadata, URLs, hashes, receipt checks and failures. This does not admit
any measurements to Starplast or change an installed dataset.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from recover_source_files import download  # noqa: E402


def retrieve(root: Path, metadata: Path, output: Path) -> dict:
    """Download processed evidence only, using filenames and URLs supplied by APIs."""
    if output.exists():
        raise ValueError("Use a new candidate retrieval snapshot")
    output.mkdir(parents=True)
    records = []
    source_files = json.loads((metadata / "rbc_project_files.json").read_text())
    project = json.loads((metadata / "rbc_pride.json").read_text())
    accession = project["accession"]
    doi = next(reference["doi"] for reference in project["references"] if reference.get("doi"))
    inputs = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in
              (metadata / "rbc_project_files.json", metadata / "rbc_pride.json", Path(__file__))}
    for item in source_files:
        if item["fileCategory"]["value"] != "SEARCH":
            continue
        if accession not in item["projectAccessions"] or item["fileSizeBytes"] > 64 << 20:
            raise ValueError("Processed archive outside the verified project or size scope")
        name = item["fileName"]
        if Path(name).name != name or not name.endswith(".zip"):
            raise ValueError("Unexpected processed archive filename")
        origin = next(location["value"] for location in item["publicFileLocations"] if location["value"].startswith("ftp://"))
        parsed = urllib.parse.urlparse(origin)
        if parsed.hostname != "ftp.pride.ebi.ac.uk":
            raise ValueError("Unexpected PRIDE file host")
        url = urllib.parse.urlunparse(parsed._replace(scheme="https"))
        record = {"project": accession, "filename": name, "reported_url": origin,
                  "requested_url": url, "transport": "HTTPS for identical PRIDE FTP host/path",
                  "retrieved_utc": datetime.now(timezone.utc).isoformat(), "status": "pending"}
        try:
            target = root / "host/erythrocyte_proteome" / accession / name
            download(url, target)
            payload = target.read_bytes()
            if len(payload) != item["fileSizeBytes"]:
                raise ValueError("Downloaded size differs from PRIDE metadata")
            sha1 = hashlib.sha1(payload).hexdigest()
            record.update(path=str(target), bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                          reported_checksum=item["checksum"], computed_sha1=sha1,
                          reported_checksum_matches_sha1=sha1 == item["checksum"])
            if sha1 != item["checksum"]:
                raise ValueError("Downloaded SHA1 differs from reported PRIDE checksum")
            with zipfile.ZipFile(target) as archive:
                record["archive_members"] = [{"name": member.filename, "bytes": member.file_size}
                                             for member in archive.infolist()]
            (target.parent / "URLS.txt").write_text(url + "\n")
            (target.parent / "SHA256SUMS.txt").write_text(record["sha256"] + "  " + name + "\n")
            record["status"] = "retrieved_validated_processed_archive"
        except Exception as error:
            record.update(status="retrieval_or_validation_failed", error=f"{type(error).__name__}: {error}")
        records.append(record)
        print(record["filename"], record["status"], record.get("error", ""), flush=True)
    query = {"resource_doi": doi, "page_size": 100}
    url = "https://api.figshare.com/v2/articles/search"
    request_record = {"url": url, "request": query, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    try:
        request = urllib.request.Request(url, data=json.dumps(query).encode(),
                                         headers={"Content-Type": "application/json", "User-Agent": "starplast-source-recovery/1"})
        with urllib.request.urlopen(request, timeout=25) as response:
            payload = response.read(4 << 20)
        articles = json.loads(payload)
        if not isinstance(articles, list):
            raise ValueError("Unexpected Figshare search response")
        (output / "figshare_search.json").write_bytes(payload)
        request_record.update(status="retrieved", articles=len(articles), sha256=hashlib.sha256(payload).hexdigest())
        for article in articles:
            details_url = article["url"]
            with urllib.request.urlopen(details_url, timeout=25) as response:
                payload = response.read(4 << 20)
            details = json.loads(payload)
            if details.get("resource_doi", "").lower() != doi.lower():
                raise ValueError("Publisher supplement DOI association differs from query")
            (output / (str(article["id"]) + "_article.json")).write_bytes(payload)
            for file in details["files"]:
                name = file["name"]
                if Path(name).name != name or not name.lower().endswith((".xlsx", ".xls", ".pdf")):
                    continue
                target = root / "host/erythrocyte_proteome" / accession / "publisher" / name
                row = {"filename": name, "url": file["download_url"], "resource_doi": doi,
                       "article_id": article["id"], "license": details.get("license"),
                       "retrieved_utc": datetime.now(timezone.utc).isoformat()}
                try:
                    download(file["download_url"], target)
                    payload = target.read_bytes()
                    if len(payload) != file["size"] or hashlib.md5(payload).hexdigest() != file["computed_md5"]:
                        raise ValueError("Publisher supplemental file size or checksum differs")
                    row.update(status="retrieved_validated_publisher_supplement", path=str(target),
                               bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
                    # Each verified file has its own receipts to avoid replacing earlier receipts.
                    target.with_name(name + ".url.txt").write_text(file["download_url"] + "\n")
                    target.with_name(name + ".sha256.txt").write_text(row["sha256"] + "  " + name + "\n")
                except Exception as error:
                    row.update(status="retrieval_or_validation_failed", error=f"{type(error).__name__}: {error}")
                records.append(row)
                print(name, row["status"], row.get("error", ""), flush=True)
    except Exception as error:
        request_record.update(status="metadata_retrieval_failed", error=f"{type(error).__name__}: {error}")
    (output / "figshare_request.json").write_text(json.dumps(request_record, indent=2) + "\n")
    (output / "downloads.json").write_text(json.dumps(records, indent=2) + "\n")
    manifest = {"input_sha256": inputs, "downloads": len(records),
                "successful": sum(r["status"].startswith("retrieved_validated") for r in records),
                "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob("*.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    """Keep the retrieval and its results in an executed annotated notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook("Public RBC candidate processed archive and publisher supplements")
    nb.md("Project identity, filenames and locations come from the frozen primary PRIDE API response. "
          "Instrument RAW files are outside this candidate review. Search publisher supplements by the "
          "DOI returned in project references, using the documented Figshare resource_doi filter. "
          "Size and deposit checksums must match. Downloaded candidates are not yet admitted measurements.")
    nb.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
            "from retrieve_host_candidates import retrieve",
            f"retrieve(Path({str(args.root)!r}), Path({str(args.metadata)!r}), Path({str(args.out)!r}))")
    nb.write(str(args.out / "retrieval.ipynb"))


if __name__ == "__main__":
    main()
