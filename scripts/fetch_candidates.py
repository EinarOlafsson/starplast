#!/usr/bin/env python3
"""Download the proposed datasets, one per slot, with provenance beside each file.

`propose_datasets.py` says what could fill a slot. This fetches it, and the difference between the
two is the difference between a citation and a dataset -- which is the distinction the slot table
exists to keep. A slot is filled when a column exists in the node table; this is the step before
that, putting the raw file where the ingestion scripts can find it.

## What is written beside every file

A manifest: the accession, the URL it came from, the size, a SHA-256, and the date. Without the
checksum a re-download cannot be told from a silently different release, and "the numbers changed
and nobody knows why" is the failure this project has already paid for once in its identifier
handling.

## Organism is checked before the bytes are

A GEO query for a parasite assay returns host studies too -- whole-blood transcriptomes of malaria
patients are *Homo sapiens*, not *Plasmodium*. Those are not junk, they are host-table data that
arrived through the wrong door, so they are recorded under the host arm rather than discarded. What
must not happen is a human blood series landing in a parasite slot, where every later claim would
be about the wrong organism entirely.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

#: Which arm a taxon belongs to. Anything else is a host or an unrelated organism.
ARM_TAXON = {"Tg": re.compile(r"toxoplasma", re.I), "Pf": re.compile(r"plasmodium", re.I)}
HOST_TAXON = re.compile(r"homo sapiens|mus musculus|anopheles", re.I)


def _get(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return fh.read()


def series_metadata(accessions: list, log=print) -> dict:
    """accession -> (taxon, n_samples, title), from GEO's own index."""
    out = {}
    for i in range(0, len(accessions), 20):
        chunk = accessions[i:i + 20]
        term = " OR ".join(f"{a}[Accession]" for a in chunk)
        try:
            found = json.loads(_get(f"{EUTILS}/esearch.fcgi?db=gds&retmode=json&retmax=200"
                                    f"&term={urllib.parse.quote(term)}"))
            ids = found.get("esearchresult", {}).get("idlist", [])
            for j in range(0, len(ids), 40):
                batch = ids[j:j + 40]
                summary = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=gds&retmode=json"
                                          f"&id={','.join(batch)}"))
                for uid in batch:
                    record = summary.get("result", {}).get(uid) or {}
                    accession = record.get("accession", "")
                    if accession.startswith("GSE"):
                        out[accession] = (record.get("taxon", ""),
                                          int(record.get("n_samples", 0) or 0),
                                          record.get("title", ""))
                time.sleep(0.35)
        except (urllib.error.URLError, ValueError) as exc:
            log(f"  metadata lookup failed for {chunk[0]}...: {type(exc).__name__}")
        time.sleep(0.35)
    return out


PRIDE_FILES = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{}/files"


def fetch_pride(accession: str, where: str, log=print, cap_mb: float = 0.0) -> dict:
    """Download a PRIDE project's PROCESSED results, and deliberately not its raw spectra.

    That distinction decides whether the download is worth making at all. A PRIDE project is mostly
    instrument files -- six 600 MB `.raw` in one of these -- and a raw file is a spectrum, not a
    measurement of a gene. Turning one into protein identifications takes a search engine, a
    sequence database and a set of parameters; it is a study of its own and none of it can be
    inferred from the deposit.

    What can fill a slot is the submitter's own search output: the protein and peptide tables PRIDE
    files under SEARCH and OTHER, three orders of magnitude smaller and already saying which protein
    was seen. So a slot filled from PRIDE is filled from what the authors concluded rather than from
    what their instrument recorded -- a weaker claim than reprocessing, and an honest one.
    """
    try:
        listing = json.loads(_get(PRIDE_FILES.format(accession), timeout=90)
                             .decode("utf8", "replace"))
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError) as exc:
        log(f"  {accession}: file listing failed ({type(exc).__name__})")
        return {"accession": accession, "status": "listing failed"}
    wanted = []
    for record in listing:
        category = str(record.get("fileCategory", "")).upper()
        name = str(record.get("fileName", ""))
        if "RAW" in category or "PEAK" in category:
            continue
        if name.lower().endswith((".raw", ".wiff", ".mzml", ".mzxml")):
            continue
        url = next((str(p.get("value")) for p in (record.get("publicFileLocations") or [])
                    if "FTP" in str(p.get("name", "")).upper()), "")
        if url:
            wanted.append((name, url, int(record.get("fileSizeBytes") or 0)))
    if not wanted:
        log(f"  {accession}: raw spectra only, nothing that names a protein")
        return {"accession": accession, "status": "raw only"}
    os.makedirs(where, exist_ok=True)
    files = []
    for name, url, size in wanted:
        target = os.path.join(where, name)
        if os.path.exists(target):
            continue
        if cap_mb and size > cap_mb * 1e6:
            files.append({"file": name, "bytes": size, "status": "over cap"})
            continue
        try:
            blob = _get(url.replace("ftp://", "https://"), timeout=900)
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            log(f"  {name}: {getattr(exc, 'code', type(exc).__name__)}")
            continue
        with open(target, "wb") as fh:
            fh.write(blob)
        files.append({"file": name, "bytes": len(blob),
                      "sha256": hashlib.sha256(blob).hexdigest(), "status": "ok"})
        log(f"  {name}: {len(blob) / 1e6:.1f} MB")
    return {"accession": accession, "url": PRIDE_FILES.format(accession), "status": "ok",
            "fetched": time.strftime("%Y-%m-%d"), "processed_files": files}


def matrix_url(accession: str) -> str:
    """Where GEO keeps a series matrix. The directory is the accession with its last three digits
    replaced by nnn, which is a GEO convention rather than anything derivable."""
    stem = accession[:-3] + "nnn" if len(accession) > 6 else accession + "nnn"
    return f"{FTP}/{stem}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def has_table(blob: bytes) -> bool:
    """Does this series matrix actually carry values, or only metadata?

    The distinction that decides whether a download was worth making. GEO publishes a series matrix
    for every series, but for a SEQUENCING submission it contains the sample descriptions and
    nothing else -- the counts live in supplementary files. A 3 kB matrix looks like a dataset in a
    directory listing and can never fill a slot, which is the same failure the ingestion loaders
    already guard against at the other end.
    """
    import gzip
    try:
        text = gzip.decompress(blob).decode("utf8", "replace")
    except OSError:
        return False
    start = text.find("!series_matrix_table_begin")
    if start < 0:
        return False
    body = text[start:].split("\n")[1:]
    return sum(1 for line in body if line.strip() and not line.startswith("!")) > 1


def supplementary(accession: str, where: str, cap_mb: float = 0.0, log=print) -> list:
    """Fetch a series' supplementary files, which is where sequencing data actually lives."""
    stem = accession[:-3] + "nnn" if len(accession) > 6 else accession + "nnn"
    listing_url = f"{FTP}/{stem}/{accession}/suppl/"
    try:
        listing = _get(listing_url, timeout=90).decode("utf8", "replace")
    except (urllib.error.URLError, urllib.error.HTTPError):
        log(f"  {accession}: no supplementary directory")
        return []
    names = sorted(set(re.findall(r'href="([^"/?][^"]*)"', listing)))
    # Relative names only. The listing page carries site-footer links too, and an absolute URL
    # scraped from one is not a file in this directory -- it 404s against the supplementary path.
    names = [n for n in names if not n.endswith("/") and "://" not in n
             and n.lower() != "filelist.txt"]
    out = []
    for name in names:
        target = os.path.join(where, name)
        if os.path.exists(target):
            continue
        try:
            blob = _get(listing_url + urllib.parse.quote(name), timeout=600)
        except (urllib.error.URLError, urllib.error.HTTPError) as exc:
            log(f"  {name}: {getattr(exc, 'code', type(exc).__name__)}")
            continue
        if cap_mb and len(blob) > cap_mb * 1e6:
            # NO cap by default, and that is the point. A cap silently converts "nobody has measured
            # this" into "we chose not to fetch it", and the slot table cannot tell those apart
            # afterwards. Twelve files were skipped by a 400 MB cap on the first run, including the
            # 2.8 GB mosquito-stage archive that three transcription slots depend on.
            log(f"  {name}: {len(blob) / 1e6:.0f} MB, over the {cap_mb:.0f} MB cap -- skipped")
            out.append({"file": name, "bytes": len(blob), "status": "over cap"})
            continue
        os.makedirs(where, exist_ok=True)
        with open(target, "wb") as fh:
            fh.write(blob)
        out.append({"file": name, "bytes": len(blob),
                    "sha256": hashlib.sha256(blob).hexdigest(), "status": "ok"})
        log(f"  {name}: {len(blob) / 1e6:.1f} MB")
    return out


def fetch(accession: str, where: str, log=print, cap_mb: float = 0.0) -> dict | None:
    """Download one series matrix. Returns its manifest entry, or None if it is not published."""
    url = matrix_url(accession)
    target = os.path.join(where, f"{accession}_series_matrix.txt.gz")
    if os.path.exists(target):
        log(f"  {accession}: already here")
        return None
    try:
        blob = _get(url, timeout=180)
    except (urllib.error.URLError, urllib.error.HTTPError) as exc:
        # A series with no matrix is ordinary: sequencing submissions often publish only
        # supplementary files. Recorded rather than retried, so the slot shows why it is still open.
        log(f"  {accession}: no series matrix ({getattr(exc, 'code', type(exc).__name__)})")
        return {"accession": accession, "url": url, "status": "no series matrix"}
    os.makedirs(where, exist_ok=True)
    with open(target, "wb") as fh:
        fh.write(blob)
    entry = {"accession": accession, "url": url, "bytes": len(blob),
             "sha256": hashlib.sha256(blob).hexdigest(),
             "fetched": time.strftime("%Y-%m-%d"), "status": "ok",
             "path": os.path.relpath(target, start=os.path.dirname(where))}
    if has_table(blob):
        log(f"  {accession}: {len(blob) / 1e6:.1f} MB, values in the matrix")
    else:
        log(f"  {accession}: matrix is metadata only, fetching supplementary files")
        entry["matrix"] = "metadata only"
        entry["supplementary"] = supplementary(accession, where, cap_mb, log=log)
    return entry


def main(argv=None) -> int:
    """`python scripts/fetch_candidates.py --root datasets/acquired --per-slot 1`."""
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--candidates", default="instructions/open/41_candidates.json")
    p.add_argument("--root", default="datasets/quarantine",
                   help="where downloads land. The default is QUARANTINE, not the application's "
                        "dataset folder: nothing is verified on arrival, and a directory nobody "
                        "has checked gets trusted by whoever ingests it next. Verification moves "
                        "a file out of here; it never lands here already trusted")
    p.add_argument("--cap-mb", type=float, default=0.0,
                   help="skip supplementary files above this size. 0 means no cap, which is the "
                        "default")
    p.add_argument("--per-slot", type=int, default=1, help="candidates to fetch per slot")
    p.add_argument("--max-slots", type=int, default=0)
    p.add_argument("--arm", default="", choices=("", "Tg", "Pf"))
    args = p.parse_args(argv)

    with open(args.candidates) as fh:
        proposals = json.load(fh)
    keys = [k for k in sorted(proposals) if not args.arm or k.startswith(f"{args.arm}::")]
    if args.max_slots:
        keys = keys[:args.max_slots]

    wanted, proteomics = [], []
    for key in keys:
        for paper in proposals[key][:args.per_slot]:
            for accession in [a.strip() for a in paper["accession"].split(",")]:
                if accession.startswith("GSE"):
                    wanted.append((key, accession))
                    break
                if accession.startswith("PXD"):
                    # A different repository with a different shape: no series matrix, no taxon
                    # field to check an arm against, and the organism filter was applied at
                    # proposal time instead. Kept in its own list rather than forced through the
                    # GEO path, which would look for a matrix that does not exist.
                    proteomics.append((key, accession))
                    break
    print(f"{len(keys)} slots, {len(wanted)} GEO series and {len(proteomics)} PRIDE projects",
          flush=True)
    meta = series_metadata(sorted({a for _k, a in wanted}), log=lambda m: print(m, flush=True))

    manifest, misrouted = {}, []
    for key, accession in wanted:
        arm = key.split("::")[0]
        taxon, samples, title = meta.get(accession, ("", 0, ""))
        if not ARM_TAXON[arm].search(taxon):
            where = "host" if HOST_TAXON.search(taxon) else "unmatched"
            misrouted.append({"slot": key, "accession": accession, "taxon": taxon,
                              "title": title, "routed_to": where})
            print(f"{key}: {accession} is {taxon or 'unknown'} -> {where}, not fetched here",
                  flush=True)
            continue
        folder = os.path.join(args.root, arm, re.sub(r"[^A-Za-z0-9]+", "_", key.split("::")[1]))
        print(f"{key} <- {accession} ({taxon}, {samples} samples)", flush=True)
        entry = fetch(accession, folder, log=lambda m: print(m, flush=True),
                      cap_mb=args.cap_mb)
        if entry:
            manifest.setdefault(key, []).append({**entry, "taxon": taxon, "samples": samples,
                                                 "title": title})
        time.sleep(0.4)

    for key, accession in proteomics:
        arm = key.split("::")[0]
        folder = os.path.join(args.root, arm, re.sub(r"[^A-Za-z0-9]+", "_", key.split("::")[1]))
        print(f"{key} <- {accession} (PRIDE)", flush=True)
        entry = fetch_pride(accession, folder, log=lambda m: print(m, flush=True),
                            cap_mb=args.cap_mb)
        if entry:
            manifest.setdefault(key, []).append(entry)
        time.sleep(0.4)

    os.makedirs(args.root, exist_ok=True)
    with open(os.path.join(args.root, "manifest.json"), "w") as fh:
        json.dump({"fetched": manifest, "misrouted": misrouted}, fh, indent=1, sort_keys=True)
    got = sum(1 for v in manifest.values() for e in v if e.get("status") == "ok")
    print(f"\n{got} files fetched into {args.root}; {len(misrouted)} series belong to another arm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
