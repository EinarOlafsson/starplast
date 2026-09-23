#!/usr/bin/env python3
"""Structures on demand: metadata ships, coordinates do not.

The app is standalone for *data* — every measurement is in the committed cache. Coordinates are the one
exception, deliberately: 6,538 AlphaFold models plus 12,265 crosslink complex CIFs are gigabytes, git is
the wrong place for them, and you do not need any of them until you click a gene.

So this module resolves a structure only when asked, in a fixed order of preference:

1. a local path already on this machine (the crosslink CIFs, or an AlphaFold mirror if one exists),
2. the AlphaFold Database, fetched over the network and cached under ``~/.cache/starplast/structures``.

Everything is best-effort. `find_structure` returns None rather than raising when a machine is offline or
a model does not exist, because a missing structure is a normal state — AlphaFold DB skips the largest
proteins, which in this proteome are disproportionately the secreted effectors people care about.
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

from .logging_util import get_logger

CACHE = os.path.join(os.path.expanduser("~"), ".cache", "starplast", "structures")

# Local mirrors, in preference order. These are machine-specific by nature; absence is not an error.
LOCAL_DIRS = [
    "/home/carruthers/foldseek/proteome_struct/tgon",          # work machine
    "/mnt/wd4tb/af3/af_output",                                 # home AF3 output
    "/mnt/wd4tb/af3_home/af3_structures",
]

# The file endpoint serves only the CURRENT release, so a pinned version 404s the moment AlphaFold
# reissues. This was pinned to v4 and every fetch had been failing silently: v3, v4 and v5 all 404
# today, only v6 answers. Because a 404 is caught below and reported as "no model for this
# accession" -- which is a real and common state -- nothing ever said so.
#
# Pinning v6 would only reset the same clock. The version is resolved per accession from the API,
# which names the current file, and the probe list is the fallback for when the API is unreachable.
_log = get_logger(__name__)

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{v}.cif"
AFDB_VERSIONS = (6, 5, 4)          # newest first; only used if the API cannot be reached


def local_structure(gene_id: str, uniprot: str | None = None) -> str | None:
    """A model already on this machine, or None."""
    # Any release, not just v4: a mirror populated at some earlier date holds whatever was current
    # then, and refusing to read a v5 file that is sitting right there would be perverse.
    names = ([f"AF-{uniprot}-F1.pdb"]
             + [f"AF-{uniprot}-F1-model_v{v}.cif" for v in AFDB_VERSIONS]) if uniprot else []
    names += [f"{gene_id}.pdb", f"{gene_id}.cif"]
    if uniprot:
        names += [f"AF3-{uniprot}.cif"]
    for d in LOCAL_DIRS:
        if not os.path.isdir(d):
            continue
        for n in names:
            p = os.path.join(d, n)
            if os.path.exists(p):
                return p
        # AF3 output uses a directory per job
        p = os.path.join(d, gene_id.lower())
        if os.path.isdir(p):
            for f in sorted(os.listdir(p)):
                if f.endswith((".cif", ".pdb")):
                    return os.path.join(p, f)
    from .structure_catalog import local_model
    try:
        return local_model(gene_id)
    except (OSError, ValueError, KeyError) as exc:
        _log.warning("Cannot read local AF3 index: %s", exc)
        return None


def _afdb_url(uniprot: str, timeout: int = 30) -> str | None:
    """Ask AlphaFold which file it currently serves for this accession.

    One request, and it removes the version from this module's assumptions entirely. Returning None
    means the API could not be reached; it does not mean there is no model, so the caller falls back
    to probing rather than concluding the protein has no structure.
    """
    try:
        req = urllib.request.Request(AFDB_API.format(acc=uniprot),
                                     headers={"User-Agent": "starplast (research)"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            entries = json.loads(r.read().decode("utf-8", "replace"))
    except (urllib.error.URLError, OSError, TimeoutError, ValueError) as exc:
        _log.warning("AlphaFold API unreachable for %s: %s: %s", uniprot, type(exc).__name__, exc)
        return None
    if isinstance(entries, list) and entries and isinstance(entries[0], dict):
        return entries[0].get("cifUrl") or None
    return None


def _cached(uniprot: str) -> str | None:
    """A previously downloaded model for this accession, whatever release it came from."""
    if not os.path.isdir(CACHE):
        return None
    for v in AFDB_VERSIONS:
        p = os.path.join(CACHE, f"AF-{uniprot}-F1-model_v{v}.cif")
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


def fetch_alphafold(uniprot: str, timeout: int = 30) -> str | None:
    """Download one AlphaFold DB model into the cache; return its path, or None if unavailable.

    The version is resolved rather than assumed, because the file endpoint serves only the current
    release. A pinned version turns every fetch into a 404, and a 404 here is indistinguishable from
    the ordinary case of a protein AlphaFold has no model for -- so the failure is completely silent.
    """
    if not uniprot:
        return None
    hit = _cached(uniprot)
    if hit:
        return hit
    os.makedirs(CACHE, exist_ok=True)
    urls = []
    resolved = _afdb_url(uniprot, timeout=timeout)
    if resolved:
        urls.append(resolved)
    # Newest first, so the probe finds the live release before the retired ones.
    urls += [AFDB_URL.format(acc=uniprot, v=v) for v in AFDB_VERSIONS]
    for url in urls:
        name = os.path.basename(url.split("?")[0]) or f"AF-{uniprot}-F1-model.cif"
        dest = os.path.join(CACHE, name)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "starplast (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            # offline, or this release is retired -- try the next candidate before giving up
            _log.debug("AlphaFold %s: %s: %s", url, type(exc).__name__, exc)
            continue
        if not data:
            continue
        with open(dest, "wb") as fh:
            fh.write(data)
        _log.info("AlphaFold %s: %s (%.1f kB)", uniprot, url, len(data) / 1024)
        return dest
    # Every candidate failed. Logged as a WARNING with the URLs tried, because this is precisely the
    # failure that ran for months: a 404 from a pinned version is indistinguishable, to the caller,
    # from the ordinary case of a protein AlphaFold has no model for.
    _log.warning("AlphaFold %s: no model retrieved from any of %d URL(s): %s",
                 uniprot, len(urls), ", ".join(urls))
    return None


def find_structure(gene_id: str, uniprot: str | None = None, allow_network: bool = True):
    """Return (path, origin) for a gene's structure, or (None, reason)."""
    p = local_structure(gene_id, uniprot)
    if p:
        return p, "local"
    if not allow_network:
        return None, "not cached, network disabled"
    if not uniprot:
        return None, "no UniProt accession known for this gene"
    p = fetch_alphafold(uniprot)
    return (p, "AlphaFold DB") if p else (None, "no AlphaFold model available")


def crosslink_model_paths(models_row, base: str) -> list:
    """Absolute paths to the Chai-1 complexes for one crosslinked pair, those that exist locally."""
    d, files = getattr(models_row, "model_dir", None), getattr(models_row, "model_files", None)
    if not d or not isinstance(files, str):
        return []
    out = [os.path.join(base, d, f) for f in files.split(";")]
    return [p for p in out if os.path.exists(p)]
