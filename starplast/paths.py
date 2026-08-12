#!/usr/bin/env python3
"""Where data lives — resolved in one place, so no module ever opens a relative path again.

Before this module existed, `build_graph.py` computed its dataset root as the repository's *parent*
directory. That is true on exactly one machine. Measured from a clean checkout, 1 of 16 registry paths
resolved; the other 15 were two directories away, and nothing said so. Any reader who cloned the
repository and followed the instructions got fifteen missing files and no explanation, which for a paper
whose whole claim is reproducibility is the wrong first impression.

There are two separate roots and confusing them is the mistake this module exists to prevent:

    the CACHE     `starplast/data/` — built artefacts, committed, 17 MB, ships inside the package.
                  The application needs this and nothing else. It must resolve offline, from an
                  installed wheel, with no configuration, forever.

    the DATASETS  `datasets/` — raw published inputs, tens of gigabytes, NOT committed. Only
                  `build_graph.py` needs it. It is allowed to be missing, and when it is missing the
                  answer is to fetch it, not to fail.

So `data_dir()` must never fail and `dataset_root()` is allowed to return a directory that is empty.

## Resolution order

Both roots take an explicit override first, because a user who sets one has a reason:

    $STARPLAST_DATA      the raw datasets tree
    $STARPLAST_CACHE     the built cache

then the layouts that actually occur, in order, and finally the platform cache directory, which is
where anything downloaded lands. `describe()` prints which candidate won, since a resolver that silently
picks the wrong directory is harder to debug than one that fails.
"""
from __future__ import annotations

import os

ENV_DATASETS = "STARPLAST_DATA"
ENV_CACHE = "STARPLAST_CACHE"

_PKG = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_PKG)


def user_cache_dir() -> str:
    """Platform cache directory. Downloads land here, so it must be writable and per-user."""
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    elif os.uname().sysname == "Darwin":
        base = os.path.expanduser("~/Library/Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    return os.path.join(base, "starplast")


# --------------------------------------------------------------------------- the built cache
def _cache_candidates():
    yield os.environ.get(ENV_CACHE)
    yield os.path.join(_PKG, "data")            # inside the package: the installed layout
    yield os.path.join(_REPO, "data")           # beside the package: the pre-0.3 repository layout
    yield os.path.join(user_cache_dir(), "data")


def data_dir() -> str:
    """The built cache. Returns the first candidate that holds a real cache, else the best guess.

    Never raises. A missing cache is reported by `check()` with an instruction, because an exception
    from an import is the least useful way to tell someone a file is absent.
    """
    # An explicit override is honoured absolutely, even when it points at nothing. Falling through to
    # the next candidate would mean a user who mistyped $STARPLAST_CACHE silently gets a DIFFERENT
    # cache than the one they asked for -- which is precisely the "resolver picks the wrong directory
    # in silence" failure this module exists to prevent. Wrong and loud beats wrong and quiet.
    override = os.environ.get(ENV_CACHE)
    if override:
        return override
    for c in _cache_candidates():
        if c and os.path.exists(os.path.join(c, "nodes.parquet")):
            return c
    for c in _cache_candidates():               # nothing built yet: name where it should go
        if c:
            return c
    return os.path.join(user_cache_dir(), "data")


def cache_file(name: str) -> str:
    return os.path.join(data_dir(), name)


REQUIRED = ("nodes.parquet", "graph.npz", "toxodb_identity.tsv")


def check() -> tuple:
    """(ok, message). What the application needs, and what to do when it is not there."""
    d = data_dir()
    missing = [f for f in REQUIRED if not os.path.exists(os.path.join(d, f))]
    if not missing:
        return True, f"cache complete at {d}"
    return False, (f"cache incomplete at {d}: missing {', '.join(missing)}. "
                   f"Run `python -m starplast.build_graph` to build it, or set ${ENV_CACHE} to an "
                   f"existing cache.")


# --------------------------------------------------------------------------- the raw datasets
def _dataset_candidates():
    yield os.environ.get(ENV_DATASETS)
    yield os.path.join(_REPO, "datasets")                     # what a clone implies
    yield os.path.join(os.path.dirname(_REPO), "datasets")    # the historical layout
    yield os.path.join(user_cache_dir(), "datasets")          # where downloads land


def dataset_roots() -> list:
    """Every dataset root that exists, in priority order.

    All of them, not just the first: after the move to `datasets/<level>/<type>/<PMID>/` some files
    live under the new taxonomy and some under the old, and a resolver that stopped at the first
    existing root silently lost the screens once already.
    """
    seen, out = set(), []
    for c in _dataset_candidates():
        if c and os.path.isdir(c) and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def dataset_root(create: bool = False) -> str:
    """The primary dataset root. Falls back to the cache directory, which is always writable."""
    roots = dataset_roots()
    if roots:
        return roots[0]
    d = os.path.join(user_cache_dir(), "datasets")
    if create:
        os.makedirs(d, exist_ok=True)
    return d


def find(*relpath) -> str | None:
    """Search every dataset root for a file. None when it is genuinely absent.

    Returning None rather than raising is deliberate: a missing raw dataset is a normal state that the
    caller decides about -- fetch it, skip it, or report it -- and a resolver is the wrong layer to
    make that decision.
    """
    rel = os.path.join(*relpath)
    for root in dataset_roots():
        p = os.path.join(root, rel)
        if os.path.exists(p):
            return p
    return None


def describe() -> str:
    """What resolved to where. A resolver that picks silently is harder to debug than one that fails."""
    ok, msg = check()
    lines = [f"cache      {data_dir()}   [{'ok' if ok else 'INCOMPLETE'}]"]
    roots = dataset_roots()
    if roots:
        lines += [f"datasets   {r}" for r in roots]
    else:
        lines.append(f"datasets   (none found; downloads would go to {dataset_root()})")
    ds_override = os.environ.get(ENV_DATASETS)
    if ds_override and not os.path.isdir(ds_override):
        # Said out loud rather than quietly skipped: a set-but-wrong override looks identical to an
        # unset one from the outside, and the user believes their configuration took effect.
        lines.append(f"WARNING    ${ENV_DATASETS} is set to {ds_override}, which is not a directory")
    lines.append(f"overrides  ${ENV_DATASETS}={os.environ.get(ENV_DATASETS) or '(unset)'}  "
                 f"${ENV_CACHE}={os.environ.get(ENV_CACHE) or '(unset)'}")
    return "\n".join(lines)


def _cli():
    """`python -m starplast.paths` — what resolved where, and what is missing.

    The first thing to run when something cannot find its data, which is why it prints the overrides
    even when they are unset: the usual fix is to set one.
    """
    print(describe())
    try:
        from . import datasets
        m = datasets.missing()
    except Exception as e:                        # noqa: BLE001 -- diagnostics must not need the registry
        print(f"\n(registry unavailable: {e})")
        return
    print(f"\n{len(datasets.REGISTRY)} registered datasets, {len(m)} not present locally")
    for k in m:
        ok, how = datasets.fetchable(k)
        print(f"  {k:24s} {'fetchable: ' + how if ok else how}")


if __name__ == "__main__":
    _cli()
