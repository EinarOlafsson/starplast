#!/usr/bin/env python3
"""GPU acceleration, offered where it helps and absent everywhere else.

Three backends, in the order they are preferred, because they cover different parts of the work:

    cuml    UMAP and HDBSCAN themselves -- the two things this program spends its minutes on
    cupy    array work: scaling, ranking, pairwise distances
    torch   the same array work, and far more likely to be installed already

None of them is a dependency. `available()` reports what is actually importable, and every function
here falls back to the NumPy path it mirrors, so a machine with no GPU runs the same code it always
did and a machine with one does not have to be configured. The switch in Preferences turns the whole
thing off even where a backend exists, because "it was slow" and "it was wrong" are different
complaints and the way to tell them apart is to rerun on the CPU.

**A GPU result that differs from the CPU result is a bug, not a speed-up** -- with two exceptions,
both stated rather than discovered later:

* **The device path is single precision.** A 3090 does float64 at a thirty-second of its float32
  rate, so distances computed here carry a relative error of about 1e-4 against SciPy's float64.
  That is invisible in what they are used for -- the continuity proxy is a Spearman correlation of
  distances, and a rank does not move for a part in ten thousand -- and the test checks the score
  itself, not only the matrix.
* **cuml's UMAP is not umap-learn.** A map built on the GPU is a different map of the same data, not
  the same map faster. It is announced in the log when it happens, because a walk whose rows came
  from both would be comparing the two libraries rather than the settings.

Ranking is exact either way: it is a permutation, and a permutation has no precision.
"""
from __future__ import annotations

import os

import numpy as np

#: The environment variable that forces the answer, for tests and for a machine where the driver is
#: present but broken. "0" disables, "1" allows what is importable.
ENV_GPU = "STARPLAST_GPU"

#: How large an array has to be before the copy to the device is worth it, measured on an RTX 3090
#: against SciPy for the distance matrix -- the only array job that survived measurement:
#:
#:      n=500    0.8x     n=2,000   1.2x     n=8,000    1.5x
#:      n=1,000  0.1x     n=4,000   1.4x     n=12,000   1.5x
#:
#: So the crossover is around 120,000 elements and the ceiling is about 1.5x. Under the threshold a
#: walk would get SLOWER for having a GPU, which is the failure this constant exists to prevent.
MIN_ELEMENTS = 120_000


def _forced():
    """The override, or None when the environment says nothing."""
    v = os.environ.get(ENV_GPU)
    return None if v is None or v == "" else v not in ("0", "false", "no")


def available() -> dict:
    """What is importable right now: {"cuml": bool, "cupy": bool, "torch": bool, "device": str}."""
    out = {"cuml": False, "cupy": False, "torch": False, "device": ""}
    if _forced() is False:
        return out
    try:
        import cuml  # noqa: F401
        out["cuml"] = True
    except Exception:
        pass
    try:
        import cupy
        out["cupy"] = bool(cupy.cuda.runtime.getDeviceCount())
    except Exception:
        pass
    try:
        import torch
        out["torch"] = bool(torch.cuda.is_available())
        if out["torch"]:
            out["device"] = torch.cuda.get_device_name(0)
    except Exception:
        pass
    return out


def enabled() -> bool:
    """Whether GPU work is both wanted and possible.

    Possible is whether anything imports. Wanted is the switch in Preferences -- and where nobody
    has touched that switch, wanted defaults to **whether a backend is there at all**: installing
    two gigabytes of CUDA wheels is a deliberate act, and making the user then find a setting to get
    what they installed is a way of wasting their afternoon.

    An explicit choice always wins, in both directions, and survives a restart. Read on every call
    rather than cached: someone who installs cuml while the program is open should get it.
    """
    forced = _forced()
    have = any(available()[k] for k in ("cuml", "cupy", "torch"))
    if forced is not None:
        return forced and have
    try:
        from PyQt6 import QtCore
        s = QtCore.QSettings("starplast", "starplast")
        want = s.value("compute/gpu", type=bool) if s.contains("compute/gpu") else have
    except Exception:
        want = False
    return bool(want) and have


def backend() -> dict:
    """Which library will do each job, with its version -- the record a run has to carry.

    A map built by cuml's UMAP is a different map of the same data, so two rows of one table built
    by two implementations would be a comparison of the libraries wearing the look of a comparison
    of settings. Every saved run and every results row records this.
    """
    out = {"umap": "", "cluster": "", "gpu": enabled()}
    if enabled() and available()["cuml"]:
        import cuml
        # getattr, not cuml.__version__: a build without that attribute used to raise here and fall
        # through to the CPU branch, so the log said "umap-learn 0.5.12 on the GPU" while cuml was
        # doing the work. A version nobody can read is a gap in the record; the wrong library name
        # is a false one.
        out["umap"] = out["cluster"] = f"cuml {getattr(cuml, '__version__', 'version unknown')}"
        return out
    try:
        import umap
        out["umap"] = f"umap-learn {umap.__version__}"
    except Exception:
        out["umap"] = "pca (umap-learn not installed)"
    try:
        import sklearn
        out["cluster"] = f"scikit-learn {sklearn.__version__}"
    except Exception:                           # pragma: no cover - sklearn is a hard dependency
        out["cluster"] = "unknown"
    return out


def backend_id() -> str:
    """The backend as one short string, for a table column: "cuml 26.8.0" or "umap-learn 0.5.9"."""
    b = backend()
    return b["umap"] if b["umap"] == b["cluster"] else f"{b['umap']} + {b['cluster']}"


def describe() -> str:
    """One line for the interface and for a saved recipe: what will do the work."""
    have = available()
    if not any(have[k] for k in ("cuml", "cupy", "torch")):
        return ('no GPU backend found -- pip install starplast-gpu  (or, in a checkout, '
                'pip install -e ".[gpu]")')
    parts = [k for k in ("cuml", "cupy", "torch") if have[k]]
    where = f" on {have['device']}" if have["device"] else ""
    state = "on" if enabled() else "off (switch it on in Preferences)"
    # What it will actually buy, because "GPU: torch" reads as a promise about UMAP and is not one.
    gain = ("UMAP and HDBSCAN move to the GPU" if have["cuml"] else
            'array work only, about 1.5x on large distance matrices -- pip install starplast-gpu '
            '(or pip install -e ".[gpu]") to move UMAP and HDBSCAN themselves')
    return f"GPU {state}: {', '.join(parts)}{where}. {gain}."


def worth_it(a) -> bool:
    """Whether an array is large enough that moving it to the device pays for itself."""
    return enabled() and getattr(a, "size", 0) >= MIN_ELEMENTS


def pairwise_distances(X):
    """Euclidean distances, on the device when that is worth doing.

    The walk's own cost centre: trustworthiness and the continuity proxy are computed from the full
    distance matrix of the subsample, which is O(n^2) in both time and memory, and it is recomputed
    for every configuration.
    """
    X = np.asarray(X, dtype=np.float64)
    if not worth_it(X):
        from scipy.spatial.distance import squareform, pdist
        return squareform(pdist(X))
    have = available()
    if have["torch"]:
        import torch
        t = torch.as_tensor(X, device="cuda", dtype=torch.float32)
        # NOT the matrix-multiplication shortcut. |a|^2 + |b|^2 - 2ab in float32 cancels badly for
        # points that are close together, which is most pairs in an embedding: measured on a
        # 1,200-point matrix it costs 5.6e-4 of relative error against SciPy, against 2e-7 for the
        # direct computation. The direct mode is slower and still far faster than the CPU.
        return torch.cdist(t, t, compute_mode="donot_use_mm_for_euclid_dist"
                           ).detach().cpu().numpy().astype(np.float64)
    import cupy
    # float64 here for the same reason, since cupy has no direct mode: the squared-norm identity is
    # the only formulation available, so the precision has to come from the dtype.
    g = cupy.asarray(X, dtype=cupy.float64)
    sq = (g * g).sum(axis=1)
    d2 = cupy.maximum(sq[:, None] + sq[None, :] - 2.0 * (g @ g.T), 0.0)
    return cupy.asnumpy(cupy.sqrt(d2))


#: Deliberately NOT here: rank scaling. Implemented on the device and measured on an 8,140 x 58
#: matrix -- the real one -- it took 244 ms against NumPy's 41 ms, six times slower, because the
#: scatter runs per column and the copy is not amortised by anything. It also disagreed with the CPU
#: by 1.4e-4, since float32 reorders near-ties and a rank is supposed to be exact. Slower AND less
#: exact is not a trade-off, so the CPU keeps that job.


def umap_class():
    """cuml's UMAP when it is there and wanted, else None so the caller uses umap-learn.

    Named rather than wrapped: cuml's UMAP takes the same arguments this project passes, and a
    wrapper that silently mapped them would be one more place for the two paths to disagree.
    """
    if not enabled() or not available()["cuml"]:
        return None
    from cuml.manifold import UMAP
    return UMAP


def hdbscan_class():
    """cuml's HDBSCAN when it is there and wanted, else None."""
    if not enabled() or not available()["cuml"]:
        return None
    from cuml.cluster import HDBSCAN
    return HDBSCAN
