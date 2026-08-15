#!/usr/bin/env python3
"""Catalog and download processed public Apicomplexan GEO and PRIDE datasets.

The archive is intentionally broader than the current Toxoplasma graph.  It retains source files for
future slot-aware integrations while excluding raw reads, raw mass spectra, alignment tracks, and
other files that would consume disk without being directly usable as gene-level features.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
import zipfile
import gzip


ROOT = Path("/mnt/firecuda2/Claude/toxoplasma_projects/datasets")
OUT = ROOT / "apicomplexa_acquisition_2026_08_15"
GENERA = ("Plasmodium", "Cryptosporidium", "Toxoplasma", "Eimeria", "Babesia",
          "Theileria", "Neospora", "Sarcocystis", "Cyclospora", "Gregarina",
          "Cystoisospora", "Hammondia")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
RAW_SUFFIXES = (".fastq", ".fq", ".bam", ".sam", ".cram", ".sra", ".raw", ".wiff",
                ".d", ".mzml", ".mgf", ".bw", ".bigwig", ".wig", ".bedgraph")
PROCESSED_SUFFIXES = (".csv", ".tsv", ".txt", ".xlsx", ".xls", ".mtx", ".rds", ".rdata",
                      ".h5", ".h5ad", ".loom", ".zip", ".gz", ".mzid", ".pepxml",
                      ".protxml", ".xml", ".json")
HREF = re.compile(r'href="([^"?]+)"')


def arguments(argv=None) -> argparse.Namespace:
    """Parse archive limits and source selection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--source", choices=("all", "geo", "pride"), default="all")
    parser.add_argument("--max-file-mb", type=float, default=512.0)
    parser.add_argument("--max-total-gb", type=float, default=20.0)
    parser.add_argument("--catalog-only", action="store_true")
    parser.add_argument("--delay", type=float, default=0.12)
    return parser.parse_args(argv)


def get(url: str, *, json_data: bool = False, timeout: int = 60):
    """Read an official repository URL with an identifying user agent."""
    request = urllib.request.Request(url, headers={"User-Agent": "starplast-dataset-archive/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read()
    return json.loads(data) if json_data else data


def api(endpoint: str, params: dict) -> dict:
    """Call NCBI E-utilities with provenance-identifying parameters."""
    params = {**params, "tool": "starplast", "email": "einar.olafsson@gmail.com"}
    return get(f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}", json_data=True)


def safe(value: str) -> str:
    """A stable directory component preserving readable taxonomy/accessions."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "unknown"


def sha256(path: Path) -> str:
    """Stream a file checksum without loading large tables into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def existing_names(output: Path) -> dict[str, str]:
    """Files already on disk outside this acquisition, used to avoid obvious duplication."""
    found = {}
    for path in ROOT.rglob("*"):
        if path.is_file() and output not in path.parents:
            found.setdefault(path.name, str(path))
    return found


def downloadable(name: str) -> bool:
    """Whether a repository file is processed/derived rather than raw acquisition data."""
    lower = urllib.parse.unquote(name).lower()
    if "raw.tar" in lower or lower.endswith(RAW_SUFFIXES):
        return False
    return lower.endswith(PROCESSED_SUFFIXES)


def directory_files(url: str) -> list[dict]:
    """List data-file links from an NCBI Apache directory."""
    try:
        html = get(url).decode("utf-8", "replace")
    except (OSError, ValueError):
        return []
    out = []
    for href in HREF.findall(html):
        name = urllib.parse.unquote(href.rsplit("/", 1)[-1])
        if href.endswith("/") or not downloadable(name):
            continue
        full = urllib.parse.urljoin(url, href)
        try:
            request = urllib.request.Request(full, method="HEAD",
                                             headers={"User-Agent": "starplast-dataset-archive/1.0"})
            with urllib.request.urlopen(request, timeout=30) as response:
                size = int(response.headers.get("Content-Length", 0))
        except (OSError, ValueError):
            size = 0
        out.append({"name": name, "url": full, "size": size})
    return out


def geo_catalog(delay: float) -> list[dict]:
    """Discover every GEO series whose indexed sample organism is a requested genus."""
    ids_by_genus: dict[str, set[str]] = {}
    all_ids = set()
    for genus in GENERA:
        result = api("esearch.fcgi", {"db": "gds", "term": f"{genus}[Organism] AND gse[Entry Type]",
                                      "retmode": "json", "retmax": 2000})
        ids = set(result.get("esearchresult", {}).get("idlist", []))
        ids_by_genus[genus] = ids
        all_ids.update(ids)
        print(f"GEO search {genus}: {len(ids)} series", flush=True)
        time.sleep(delay)
    summaries = {}
    ordered = sorted(all_ids, key=int)
    for start in range(0, len(ordered), 150):
        batch = ordered[start:start + 150]
        data = api("esummary.fcgi", {"db": "gds", "id": ",".join(batch), "retmode": "json"})
        summaries.update({uid: data["result"][uid] for uid in batch if uid in data.get("result", {})})
        time.sleep(delay)
    records = []
    for uid, item in summaries.items():
        accession = item.get("accession", "")
        if not accession.startswith("GSE"):
            continue
        genus = next((name for name, ids in ids_by_genus.items() if uid in ids), "Apicomplexa")
        stem = re.sub(r"\d{1,3}$", "nnn", accession)
        base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/"
        files = directory_files(base + "suppl/") + directory_files(base + "matrix/")
        records.append({"repository": "GEO", "accession": accession, "genus": genus,
                        "organism": item.get("taxon", genus), "title": item.get("title", ""),
                        "assay": item.get("gdstype", ""), "pubmed": ";".join(item.get("pubmedids", [])),
                        "landing_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                        "files": files})
        if len(records) % 50 == 0:
            print(f"GEO catalogued {len(records)}/{len(summaries)}", flush=True)
        time.sleep(delay)
    return records


def pride_catalog(delay: float) -> list[dict]:
    """Discover PRIDE projects and retain those whose organism is in the requested genera."""
    projects = {}
    for genus in (*GENERA, "Apicomplexa"):
        for page in range(20):
            url = f"{PRIDE}/search/projects?" + urllib.parse.urlencode(
                {"keyword": genus, "pageSize": 100, "page": page, "sortDirection": "DESC"})
            try:
                rows = get(url, json_data=True)
            except (OSError, ValueError):
                break
            for row in rows:
                projects[row.get("accession", "")] = row
            if len(rows) < 100:
                break
            time.sleep(delay)
        print(f"PRIDE search {genus}: {len(projects)} cumulative candidates", flush=True)
    records = []
    for index, accession in enumerate(sorted(projects), 1):
        if not accession.startswith("PXD"):
            continue
        try:
            detail = get(f"{PRIDE}/projects/{accession}", json_data=True)
        except (OSError, ValueError):
            continue
        organisms = [row.get("name", "") if isinstance(row, dict) else str(row)
                     for row in detail.get("organisms", [])]
        genus = next((g for g in GENERA if any(g.lower() in o.lower() for o in organisms)), None)
        if genus is None:
            continue
        try:
            rows = get(f"{PRIDE}/projects/{accession}/files?pageSize=10000&page=0", json_data=True)
        except (OSError, ValueError):
            rows = []
        files = []
        for row in rows:
            name = row.get("fileName", "")
            category = row.get("fileCategory", {}).get("value", "")
            if category not in ("RESULT", "OTHER") or not downloadable(name):
                continue
            locations = row.get("publicFileLocations", [])
            ftp = next((loc.get("value", "") for loc in locations
                        if loc.get("name") == "FTP Protocol"), "")
            url = ftp.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
            if url:
                files.append({"name": name, "url": url,
                              "size": int(row.get("fileSizeBytes") or 0), "category": category})
        records.append({"repository": "PRIDE", "accession": accession, "genus": genus,
                        "organism": "; ".join(organisms), "title": detail.get("title", ""),
                        "assay": "; ".join(x.get("name", "") for x in detail.get("experimentTypes", [])),
                        "pubmed": ";".join(str(r.get("pubmedId", "")) for r in detail.get("references", [])),
                        "landing_url": f"https://www.ebi.ac.uk/pride/archive/projects/{accession}",
                        "files": files})
        if index % 25 == 0:
            print(f"PRIDE checked {index}/{len(projects)}", flush=True)
        time.sleep(delay)
    return records


def retrieve(record: dict, file: dict, output: Path, max_file: int,
             remaining: int, duplicates: dict[str, str], catalog_only: bool) -> dict:
    """Download one eligible file atomically and return its manifest row."""
    row = {k: record.get(k, "") for k in ("repository", "accession", "genus", "organism",
                                           "title", "assay", "pubmed", "landing_url")}
    row.update({"file": file["name"], "download_url": file["url"], "bytes": file.get("size", 0),
                "sha256": "", "status": "catalogued", "local_path": "", "error": ""})
    size = int(file.get("size") or 0)
    if file["name"] in duplicates:
        row.update(status="already_present", local_path=duplicates[file["name"]])
        return row
    if size > max_file:
        row["status"] = "skipped_file_limit"
        return row
    if size and size > remaining:
        row["status"] = "skipped_total_limit"
        return row
    if catalog_only:
        return row
    directory = output / safe(record["genus"]) / safe(record["organism"]) / record["accession"]
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / safe(file["name"])
    try:
        if not target.exists():
            part = target.with_suffix(target.suffix + ".part")
            request = urllib.request.Request(file["url"],
                                             headers={"User-Agent": "starplast-dataset-archive/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response, part.open("wb") as handle:
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
            os.replace(part, target)
        digest = sha256(target)
        if target.suffix == ".zip":
            with zipfile.ZipFile(target) as archive:
                if archive.testzip() is not None:
                    raise ValueError("ZIP member checksum failed")
        elif target.suffix == ".gz":
            with gzip.open(target, "rb") as handle:
                handle.read(1024)
        row.update(status="downloaded", local_path=str(target), bytes=target.stat().st_size,
                   sha256=digest)
        duplicates[file["name"]] = str(target)
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        part = target.with_suffix(target.suffix + ".part")
        if part.exists():
            part.unlink()
    return row


def write_outputs(records: list[dict], rows: list[dict], output: Path) -> None:
    """Write machine-readable provenance and a concise human inventory."""
    output.mkdir(parents=True, exist_ok=True)
    catalog = [{**{k: v for k, v in record.items() if k != "files"}, "files": record["files"]}
               for record in records]
    (output / "catalog.json").write_text(json.dumps(catalog, indent=2, ensure_ascii=False) + "\n",
                                         encoding="utf-8")
    fields = ("repository", "accession", "genus", "organism", "title", "assay", "pubmed",
              "landing_url", "file", "download_url", "bytes", "sha256", "status", "local_path",
              "error")
    with (output / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(rows)
    statuses = {}
    for row in rows:
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1
    genera = sorted({record["genus"] for record in records})
    downloaded = sum(int(row["bytes"] or 0) for row in rows if row["status"] == "downloaded")
    text = ["# Pan-Apicomplexan processed-data archive — 2026-08-15", "",
            "Official GEO and PRIDE source files for future Starplast versions. Raw sequencing reads,",
            "raw mass spectra, BAM/coverage tracks, and files over the configured limits are catalogued",
            "but not downloaded. No file in this archive is integrated into the current Toxoplasma graph.", "",
            f"- {len(records):,} repository records across {', '.join(genera)}",
            f"- {len(rows):,} eligible processed files considered",
            f"- {downloaded / 1024**3:.2f} GiB downloaded/verified in this run",
            "- status counts: " + ", ".join(f"{k}={v}" for k, v in sorted(statuses.items())), "",
            "`catalog.json` preserves dataset metadata and every candidate file. `manifest.csv` records",
            "the official landing/download URLs, local path, byte count, SHA-256, and failure/skip state.", ""]
    (output / "README.md").write_text("\n".join(text), encoding="utf-8")


def main(argv=None) -> int:
    """Discover, download within limits, verify, and manifest the archive."""
    args = arguments(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    if args.source in ("all", "geo"):
        records.extend(geo_catalog(args.delay))
    if args.source in ("all", "pride"):
        records.extend(pride_catalog(args.delay))
    records.sort(key=lambda row: (row["genus"], row["repository"], row["accession"]))
    duplicates = existing_names(args.output)
    total_limit = int(args.max_total_gb * 1024**3)
    max_file = int(args.max_file_mb * 1024**2)
    used = 0
    rows = []
    for record in records:
        for file in record["files"]:
            row = retrieve(record, file, args.output, max_file, total_limit - used,
                           duplicates, args.catalog_only)
            rows.append(row)
            if row["status"] == "downloaded":
                used += int(row["bytes"])
                print(f"downloaded {record['accession']} {file['name']} "
                      f"({used / 1024**3:.2f}/{args.max_total_gb:.1f} GiB)", flush=True)
    write_outputs(records, rows, args.output)
    print(f"{len(records)} datasets, {len(rows)} files -> {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
