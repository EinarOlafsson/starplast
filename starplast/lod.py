"""Level of detail: what a point stands for at each of the three tiers.

The galaxy tier used to draw one centroid per compartment, and it was wrong in a way that looked
merely ugly. A compartment's genes are scattered across the whole map -- measured on the shipped
embedding, the median compartment's members spread 0.37 of the map radius from their own mean, while
the median compartment centroid sits only 0.17 of the radius from the centre of the map. Averaging
points that are spread out lands you near the middle, so all 27 compartments collapsed into a central
blob. The picture was not a rendering fault; it was an honest average of a quantity that has no
spatial meaning, drawn as though it did.

That the collapse happens at all is the same fact the held-out search reports from the other
direction: `compartment` is the worst-recovered target, below even the study-effort control. The map
is not organised by localisation, so localisation centroids cannot organise the map.

So the galaxy tier is computed from the embedding itself instead, which is what HANDOFF means by
"data-driven, not invented tiers". Coarse spatial structure is found by connected components over an
occupancy grid: bin the cloud, keep cells with enough genes in them, and join cells that touch. No
clustering library is involved, so the browser still opens without scikit-learn, and the result is
deterministic -- a galaxy that moves between runs is not a landmark.
"""
from __future__ import annotations

import numpy as np

# Bins per axis. 10 puts roughly 8 genes in an occupied cell for this proteome, which is coarse
# enough that the tier is genuinely a summary and fine enough that two dense lobes do not merge
# through the thin bridge between them.
GRID = 10
MIN_CELL = 3            # a cell below this is noise, and joining through it welds unrelated lobes
MIN_MEMBERS = 40        # a component below this is not a landmark worth drawing at galaxy scale


def _occupancy(xyz: np.ndarray, grid: int):
    """Bin points into a cubic grid over their own extent. Returns (cell index per point, counts)."""
    lo = xyz.min(0)
    span = xyz.max(0) - lo
    # A degenerate axis (every gene at the same z) would divide by zero and put everything in one
    # cell; treating its span as one bin wide is the same answer without the warning.
    span = np.where(span > 0, span, 1.0)
    ijk = np.clip(((xyz - lo) / span * grid).astype(np.int64), 0, grid - 1)
    flat = (ijk[:, 0] * grid + ijk[:, 1]) * grid + ijk[:, 2]
    return ijk, flat


def galaxies(xyz: np.ndarray, grid: int = GRID, min_cell: int = MIN_CELL,
             min_members: int = MIN_MEMBERS) -> np.ndarray:
    """Label each gene with a coarse spatial component, or -1 for none.

    Labels are assigned in descending order of size, so galaxy 0 is always the largest structure in
    the map. That keeps the numbering stable when the embedding is rebuilt slightly differently,
    which matters because these are the things a user learns the names of.
    """
    xyz = np.asarray(xyz, dtype=np.float64)
    n = len(xyz)
    if n == 0:
        return np.zeros(0, dtype=np.int64)

    ijk, flat = _occupancy(xyz, grid)
    counts = np.bincount(flat, minlength=grid ** 3)
    occupied = counts >= min_cell
    if not occupied.any():
        return np.full(n, -1, dtype=np.int64)

    # Connected components over occupied cells, 6-connected. Diagonal touching would join lobes that
    # meet only at a corner, which is a single point of contact and not a shared structure.
    comp = np.full(grid ** 3, -1, dtype=np.int64)
    cur = 0
    occ_cells = np.flatnonzero(occupied)
    for start in occ_cells:
        if comp[start] != -1:
            continue
        stack = [int(start)]
        comp[start] = cur
        while stack:
            c = stack.pop()
            k = c % grid
            j = (c // grid) % grid
            i = c // (grid * grid)
            for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                a, b, d = i + di, j + dj, k + dk
                if not (0 <= a < grid and 0 <= b < grid and 0 <= d < grid):
                    continue
                nb = (a * grid + b) * grid + d
                if occupied[nb] and comp[nb] == -1:
                    comp[nb] = cur
                    stack.append(nb)
        cur += 1

    raw = comp[flat]                      # -1 for genes in cells too sparse to be occupied
    out = np.full(n, -1, dtype=np.int64)
    if cur == 0:
        return out
    sizes = np.bincount(raw[raw >= 0], minlength=cur)
    # Rank by size so 0 is the largest, and drop components too small to be landmarks.
    order = [c for c in np.argsort(-sizes) if sizes[c] >= min_members]
    for new, old in enumerate(order):
        out[raw == old] = new
    return out


def centroids(xyz: np.ndarray, labels: np.ndarray):
    """Centroid, membership and mean spread per label, ordered by label. Ignores -1.

    `spread` travels with the centroid because it is what says whether the centroid means anything.
    A tight component's mean is a landmark; a diffuse one's mean is the middle of the map, and the
    caller can then decline to draw it rather than drawing a dot that implies a structure.
    """
    labels = np.asarray(labels)
    xyz = np.asarray(xyz, dtype=np.float64)
    keep = labels >= 0
    if not keep.any():
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int64), np.zeros(0)
    ids = np.unique(labels[keep])
    pos = np.stack([xyz[labels == i].mean(0) for i in ids])
    num = np.array([int((labels == i).sum()) for i in ids], dtype=np.int64)
    spread = np.array([float(np.linalg.norm(xyz[labels == i] - xyz[labels == i].mean(0), axis=1).mean())
                       for i in ids])
    return pos, num, spread


def dominant(values, labels: np.ndarray, min_frac: float = 0.0, ignore=None):
    """The most common value within each component, and the fraction of the component holding it.

    Returned as a fraction rather than as a bare name so a caller can refuse to label a component
    that has no majority. "40% cytosol" is a description; "cytosol" would be a claim.

    `ignore` drops absence labels before counting, and without it every galaxy in this map is named
    "unassigned" -- true, since most of the proteome has no measured call, and useless as a landmark.
    The fraction is then of the *named* members, so "cytosol 45%" means 45% of the genes in that
    galaxy that have any call at all, not 45% of the galaxy.
    """
    values = np.asarray(values, dtype=object)
    labels = np.asarray(labels)
    ignore = {str(x).lower() for x in (ignore or ())}
    out = {}
    for i in np.unique(labels[labels >= 0]):
        v = values[labels == i].astype(str)
        if ignore:
            v = v[~np.isin(np.char.lower(v.astype(str)), list(ignore))]
        if len(v) == 0:
            continue
        uniq, cnt = np.unique(v, return_counts=True)
        j = int(np.argmax(cnt))
        frac = float(cnt[j]) / float(len(v))
        if frac >= min_frac:
            out[int(i)] = (str(uniq[j]), frac)
    return out
