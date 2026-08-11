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

import os
import urllib.error
import urllib.request

CACHE = os.path.join(os.path.expanduser("~"), ".cache", "starplast", "structures")

# Local mirrors, in preference order. These are machine-specific by nature; absence is not an error.
LOCAL_DIRS = [
    "/home/carruthers/foldseek/proteome_struct/tgon",          # work machine
    "/mnt/wd4tb/af3/af_output",                                 # home AF3 output
    "/mnt/wd4tb/af3_home/af3_structures",
]

AFDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.cif"


def local_structure(gene_id: str, uniprot: str | None = None) -> str | None:
    """A model already on this machine, or None."""
    names = [f"AF-{uniprot}-F1.pdb", f"AF-{uniprot}-F1-model_v4.cif"] if uniprot else []
    names += [f"{gene_id}.pdb", f"{gene_id}.cif"]
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
    return None


def fetch_alphafold(uniprot: str, timeout: int = 30) -> str | None:
    """Download one AlphaFold DB model into the cache; return its path, or None if unavailable."""
    if not uniprot:
        return None
    os.makedirs(CACHE, exist_ok=True)
    dest = os.path.join(CACHE, f"AF-{uniprot}-F1-model_v4.cif")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    try:
        req = urllib.request.Request(AFDB_URL.format(acc=uniprot),
                                     headers={"User-Agent": "starplast (research)"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        if not data:
            return None
        with open(dest, "wb") as fh:
            fh.write(data)
        return dest
    except (urllib.error.URLError, OSError, TimeoutError):
        # offline, or AlphaFold has no model for this accession -- both are ordinary
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
