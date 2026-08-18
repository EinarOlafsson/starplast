#!/usr/bin/env python3
"""A provenance-first archive of pan-Apicomplexan source files. Instruction 38.

Deliberately NOT wired into the graph. This builds a future-ready store so that extending the
leakage-aware slot hierarchy to another parasite does not mean repeating dataset discovery, and
integrating it now would quietly make Toxoplasma claims out of Cryptosporidium files.

Three properties, each of which exists because its absence has cost this project time:

**Nothing is overwritten.** A file already on disk is left exactly as it is and recorded as
`present`, because the archive is a store of what was retrieved and when, and a re-run that silently
replaced a 2019 download with a 2026 one would destroy the provenance it exists to keep.

**A failure is recorded, not omitted.** Controlled-access, withdrawn, supplement-less and
plain-404 sources go into the manifest with their status. An archive that lists only what worked is
an archive that sends the next session to rediscover the same dead ends -- which is exactly what the
empty slots' `blocked_by` verdicts prevent in the slot atlas.

**Every retrieved file carries a checksum and a size.** "Downloaded" is not a claim anyone can check;
a sha256 is.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

#: What the manifest can say about a source, and each is a different fact.
#:
#: retrieved  -- fetched now, checksummed, readable.
#: present    -- already on disk and left alone; the archive does not overwrite.
#: unreachable-- the URL failed. The source may be real; the route was not.
#: restricted -- controlled access, and no amount of retrying changes that.
#: absent     -- the study exists and publishes no downloadable table.
STATUSES = ("retrieved", "present", "unreachable", "restricted", "absent")


@dataclass
class Source:
    """One downloadable table, described well enough that somebody else could fetch it.

    `slot_family` is a guess about where this would land IF it were integrated, and it is recorded as
    a guess: the archive does not fill slots, and a file catalogued as `fitness` that turns out to be
    a growth curve should be discoverable as a mistake rather than believed.
    """
    key: str
    species: str
    study: str
    url: str = ""
    accession: str = ""
    doi: str = ""
    repository: str = ""
    licence: str = ""
    scope: str = ""
    slot_family: str = ""
    status: str = ""
    note: str = ""

    def directory(self, root: str) -> str:
        """Where this source's files live: one directory per species, then per study.

        The study name is sanitised rather than trusted: it comes from a catalogue somebody wrote,
        and a slash in it would silently scatter one study across two directories.
        """
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in self.study)[:80]
        return os.path.join(root, self.species.replace(" ", "_"), safe or self.key)


def checksum(path: str, chunk: int = 1 << 20) -> str:
    """sha256 of a file, read in chunks so a multi-gigabyte archive does not become a memory limit."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def readable(path: str) -> str:
    """Whether the file is a table something can open, and what shape. "" when it is not one.

    Instruction 38 asks for verification rather than a size check, and the distinction matters: a
    404 page saved as `.csv` is 4 kB of valid HTML, and every archive that only checked sizes has one.
    """
    try:
        if path.endswith((".fasta", ".fa", ".faa", ".fna")):
            # Counted rather than assumed: the bulk of this archive is annotation FASTA, and
            # reporting it as "not checked" would leave the largest part of the manifest making no
            # claim at all about whether the bytes are what they say.
            with open(path, "r", errors="replace") as fh:
                n = sum(1 for line in fh if line.startswith(">"))
            return f"{n:,} sequences" if n else "UNREADABLE: no FASTA header"
        if path.endswith((".gff", ".gff3")):
            with open(path, "r", errors="replace") as fh:
                rows = sum(1 for line in fh if line.strip() and not line.startswith("#"))
            return f"{rows:,} feature rows" if rows else "UNREADABLE: no features"
        if path.endswith((".csv", ".tsv", ".txt")):
            import pandas as pd
            sep = "," if path.endswith(".csv") else "\t"
            frame = pd.read_csv(path, sep=sep, nrows=50, engine="python", on_bad_lines="skip")
            return f"{len(frame.columns)} columns"
        if path.endswith((".xlsx", ".xls")):
            import pandas as pd
            return f"{len(pd.read_excel(path, nrows=20).columns)} columns"
        if path.endswith((".gz", ".zip", ".tar", ".tgz")):
            import tarfile
            import zipfile
            if zipfile.is_zipfile(path):
                return f"zip, {len(zipfile.ZipFile(path).namelist())} members"
            if tarfile.is_tarfile(path):
                return "tar"
            import gzip
            with gzip.open(path, "rb") as fh:
                fh.read(1024)
            return "gzip"
    except Exception as e:                       # a corrupt download is a finding, not a crash
        return f"UNREADABLE: {type(e).__name__}"
    return "not checked"


@dataclass
class Record:
    """What happened to one source, which is what the manifest is made of."""
    key: str = ""
    species: str = ""
    study: str = ""
    status: str = ""
    path: str = ""
    sha256: str = ""
    bytes: int = 0
    shape: str = ""
    retrieved: str = ""
    url: str = ""
    accession: str = ""
    doi: str = ""
    repository: str = ""
    licence: str = ""
    scope: str = ""
    slot_family: str = ""
    note: str = ""


def fetch(source: Source, root: str, opener=None, timeout: int = 120) -> Record:
    """Retrieve one source into the archive, or record why it was not retrieved.

    `opener` is injected so the whole path is testable without a network: the interesting behaviour
    here is what happens on failure and on a file that already exists, and neither should need a
    live server to exercise.
    """
    record = Record(**{k: v for k, v in asdict(source).items()
                       if k in Record.__dataclass_fields__ and k != "status"})
    if source.status in ("restricted", "absent"):
        # Declared in the catalogue rather than discovered here: some sources are known to be
        # controlled-access or to publish no table, and pretending to try them wastes a request and
        # reports a network failure where the truth is a licence.
        record.status = source.status
        return record
    directory = source.directory(root)
    os.makedirs(directory, exist_ok=True)
    name = os.path.basename(source.url.split("?")[0]) or f"{source.key}.data"
    path = os.path.join(directory, name)
    if os.path.exists(path):
        record.status, record.path = "present", path
        record.sha256, record.bytes = checksum(path), os.path.getsize(path)
        record.shape = readable(path)
        return record
    if opener is None:
        import urllib.request
        opener = lambda url: urllib.request.urlopen(url, timeout=timeout)     # noqa: E731
    try:
        with opener(source.url) as response, open(path, "wb") as out:
            out.write(response.read())
    except Exception as e:
        if os.path.exists(path):
            os.remove(path)                      # a partial file is worse than none
        record.status = "unreachable"
        record.note = f"{type(e).__name__}: {str(e)[:120]}"
        return record
    record.status, record.path = "retrieved", path
    record.sha256, record.bytes = checksum(path), os.path.getsize(path)
    record.shape = readable(path)
    record.retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    return record


def build(sources, root: str, manifest: str, opener=None, log=print) -> list:
    """Fetch every source and write the manifest. Returns the records.

    The manifest is written after EVERY source rather than at the end, so a run interrupted halfway
    keeps what it retrieved -- which for a multi-gigabyte archive over a slow link is the difference
    between a resumable job and a wasted afternoon.
    """
    os.makedirs(root, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(manifest)) or ".", exist_ok=True)
    records = []
    for source in sources:
        record = fetch(source, root, opener=opener)
        records.append(record)
        log(f"{record.status:12} {record.species:22} {record.study[:44]}"
            + (f"  {record.shape}" if record.shape else "")
            + (f"  -- {record.note}" if record.note else ""))
        write_manifest(records, manifest)
    counts = {}
    for record in records:
        counts[record.status] = counts.get(record.status, 0) + 1
    log("  " + ", ".join(f"{n} {status}" for status, n in sorted(counts.items())))
    return records


def write_manifest(records, path: str) -> str:
    """The manifest as CSV and JSON. CSV because a person reads it, JSON because a script does."""
    fields = list(Record.__dataclass_fields__)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))
    with open(os.path.splitext(path)[0] + ".json", "w") as fh:
        json.dump([asdict(r) for r in records], fh, indent=1)
    return path


def load_catalog(path: str) -> list:
    """The checked-in catalogue of sources, which is what makes discovery reproducible."""
    with open(path) as fh:
        return [Source(**row) for row in json.load(fh)]


# --------------------------------------------------------------------------- discovery
#: The VEuPathDB sites and the genera each one serves. This is the only part of discovery that is
#: written down rather than resolved, because it is the one part that is a fact about the WEB rather
#: than about the data: which hostname serves which clade. Everything below it -- the release number,
#: the organisms, the files -- is read from the site's own index at run time.
#:
#: The rule behind that, learned the expensive way: a release number typed into a URL is correct
#: until the next release, and then it is a 404 that looks like a missing dataset.
VEUPATHDB_SITES = {
    "plasmodb.org": "Plasmodium and other haemosporidia",
    "cryptodb.org": "Cryptosporidium, Cyclospora, Gregarina",
    "piroplasmadb.org": "Babesia, Theileria",
    "toxodb.org": "Toxoplasma, Neospora, Eimeria, Sarcocystis, Besnoitia",
}

#: The files worth archiving per organism. Gene-level and slot-compatible, in the sense instruction 38
#: asks for: sequence and annotation are what every other evidence family is keyed against, so an
#: archive without them cannot be joined to anything later.
WANTED = ("AnnotatedProteins.fasta", "AnnotatedTranscripts.fasta", ".gff", "InterproDomains.txt",
          "GeneAliases.txt", "Orthologs.txt")


def _listing(url: str, timeout: int = 40) -> str:
    import urllib.request
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf8", "replace")


def discover_veupathdb(sites=None, wanted=None, limit_organisms: int = 0, fetcher=None,
                       log=print) -> list:
    """Enumerate downloadable gene-level files across the VEuPathDB sites, resolving as it goes.

    Reproducible in the sense instruction 38 requires: nothing here is a remembered accession. The
    site-to-clade table is checked in, the release is whatever `Current_Release` points at today, the
    organisms come from that release's own index, and the files come from each organism's directory.
    Re-run next year and it finds next year's release without an edit.

    `fetcher` is injected so the whole walk can be tested without a network.
    """
    import re
    sites = sites or VEUPATHDB_SITES
    wanted = wanted or WANTED
    get = fetcher or _listing
    out = []
    for site, clade in sites.items():
        base = f"https://{site}/common/downloads/Current_Release/"
        try:
            index = get(base)
        except Exception as e:
            log(f"{site}: release index unreachable ({type(e).__name__})")
            out.append(Source(key=f"{site}-index", species=clade, study=f"{site} release index",
                              url=base, repository=site, status="unreachable",
                              note=f"{type(e).__name__}: {str(e)[:100]}"))
            continue
        organisms = sorted(set(re.findall(r'href="([A-Za-z][A-Za-z0-9_]+)/"', index)))
        if limit_organisms:
            organisms = organisms[:limit_organisms]
        log(f"{site}: {len(organisms)} organisms in the current release")
        for organism in organisms:
            folder = f"{base}{organism}/fasta/data/"
            try:
                files = re.findall(r'href="([^"?/]+)"', get(folder))
            except Exception:
                # A great many organisms publish no fasta directory. That is a fact about the
                # release, not a failure of the walk, so it is not recorded as unreachable.
                continue
            for name in files:
                if not any(token in name for token in wanted):
                    continue
                out.append(Source(
                    key=f"{site}:{organism}:{name}", species=organism, study=f"{organism} annotation",
                    url=folder + name, repository=site, licence="VEuPathDB terms of use",
                    scope=clade, slot_family="sequence and annotation"))
    return out
