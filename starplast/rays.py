#!/usr/bin/env python3
"""Light that can be blocked, and light that bounces -- without a ray tracer.

## What this does

A light in `lighting` reaches every gene that faces it, however much of the map is in the way. That
is what makes a lit cloud look like a painted cloud: the far side of a dense cluster is as bright as
its front, so nothing occludes anything and the eye gets no depth from the light. What a reader
means by "point at it and it lights up" is the other thing -- a gene lights up when there is a clear
line between it and the light, and stays dark when the cluster in front of it is in the way.

The honest way to get that is a shadow ray per gene per light. The expensive way to answer a shadow
ray is to intersect it against the scene, which for 8,140 unconnected points means either a spatial
structure and a per-ray walk, or a GPU pipeline this renderer does not have.

## What it actually is

**Voxel occupancy plus fixed-step marching.** The cloud is binned once into a coarse grid -- how
much of each cell is occupied -- and a shadow ray is answered by sampling that grid at a dozen
points along the segment and multiplying the transmittance through. It is the same trick a volume
renderer uses for a cloud or a candle flame, and a cloud of points is exactly that kind of subject:
there is nothing here with a surface to intersect, so there is nothing to be exact about. Sampling
the density along the line IS the right model, not an approximation of a better one.

The costs are honest and worth stating. It is soft, not sharp: a gene behind one lonely neighbour is
barely shadowed, and one behind a dense cluster goes properly dark, with no hard edge anywhere in
between. It cannot tell you which gene did the blocking. And it is coarse -- at 48 cells across, two
genes in the same cell cannot shadow each other at all.

**Bounce** is the same grid read the other way. A ray leaves the light, marches until it meets an
occupied cell, and that cell becomes a small light of its own. One bounce, no recursion, and the
"colour" it bounces is the light's own -- there is no radiosity here, and nothing here is solving
the rendering equation. What it buys is the thing that was asked for: a light you can aim into the
map that lands somewhere and makes THAT place glow.
"""
from __future__ import annotations

import numpy as np

#: Cells across the longest side of the map. 48 puts a few genes in each occupied cell of a dense
#: cluster, which is the resolution at which "a cluster blocks light" is true and "this gene blocks
#: that gene" is not -- and the second claim would be a lie at any resolution, since a gene is a
#: dimensionless point that would block nothing at all.
RESOLUTION = 48

#: How many genes in a cell before it stops light altogether. Low, because a UMAP cluster is dense:
#: at 6 the core of a cluster is opaque and its edges are haze, which is what a reader sees when
#: they look at one.
OPAQUE = 6.0

#: Samples along a shadow ray, as a range. The count is not fixed, and fixing it was a bug worth
#: keeping the note for: at twelve samples over a long ray the step between them is several cells
#: wide, so a thin occluder falls between two samples and the ray reports a clear line THROUGH a
#: wall. Measured on a dense slab with a gene either side of it, the far gene came back fully lit.
#: The step has to follow the grid, so the count follows the ray's length in cells -- floored, so a
#: short ray is still integrated, and capped, because this runs per gene per light per frame.
SAMPLES = (12, 28)

#: How much of the light a fully occupied cell absorbs at one sample. Under 1 so that a single dense
#: cell dims rather than extinguishes: an all-or-nothing shadow on a scatter this coarse reads as
#: points blinking out rather than as a shadow.
ABSORB = 0.55


class Grid:
    """Where the map is solid, as a coarse box of densities.

    Built once per set of coordinates and reused for every light and every frame -- the cloud does
    not move between redraws, and rebuilding this at 60 Hz would cost more than the shading it is
    there to inform.
    """

    def __init__(self, occupancy, lo, cell):
        self.occupancy = occupancy
        self.lo = lo
        self.cell = cell

    @property
    def resolution(self) -> int:
        """Cells along one side."""
        return int(self.occupancy.shape[0])

    def cells_at(self, points) -> tuple:
        """Grid indices for an array of world points, clamped to the box."""
        idx = np.floor((np.asarray(points, dtype=float) - self.lo) / self.cell).astype(int)
        np.clip(idx, 0, self.resolution - 1, out=idx)
        return idx[..., 0], idx[..., 1], idx[..., 2]

    def transmittance(self, coords, source, samples=None, absorb: float = ABSORB) -> np.ndarray:
        """How much of the light from `source` survives the trip to each of `coords`, in [0, 1].

        Sampled strictly BETWEEN the two ends. A sample at the gene's own position finds the gene
        itself and every gene shadows itself into darkness; a sample at the light's finds whatever
        the light is standing in and shadows the entire map at once.

        The absorption per sample is scaled by how far the step covers, so a ray sampled coarsely
        and one sampled finely report the same shadow -- without it, capping the sample count would
        quietly brighten exactly the long rays that need the samples most.
        """
        coords = np.asarray(coords, dtype=float)
        if coords.size == 0:
            return np.ones((0, 1))
        source = np.asarray(source, dtype=float)
        reach = source - coords
        dist = np.linalg.norm(reach, axis=1)
        n = int(samples or np.clip(round(float(dist.max()) / self.cell), *SAMPLES))
        t = np.linspace(0.0, 1.0, n + 2)[1:-1]
        along = coords[:, None, :] + reach[:, None, :] * t[None, :, None]
        density = self.occupancy[self.cells_at(along)]
        # The step is EACH ray's own length divided by the samples, not the longest ray's. Sharing
        # one step across the map charges a gene two cells from the light the same optical depth as
        # one right across it: measured, the pool around the pointer vanished entirely (1.02x the
        # rest of the map) because the genes nearest the light were being shadowed the hardest.
        step = dist / (n + 1) / self.cell
        return np.exp(-absorb * step * density.sum(axis=1))[:, None]

    def first_hit(self, origin, directions, reach: float, steps: int = 96):
        """Where each ray from `origin` first meets something, or nothing where it meets nothing.

        Fixed-step marching rather than a DDA walk of the grid: a DDA visits every cell a ray
        crosses and is the right algorithm for an exact answer, and there is no exact answer here to
        be right about -- the grid is already a blur of a cloud that has no surfaces.
        """
        origin = np.asarray(origin, dtype=float)
        dirs = np.atleast_2d(np.asarray(directions, dtype=float))
        norms = np.linalg.norm(dirs, axis=1, keepdims=True)
        dirs = dirs / np.where(norms > 1e-9, norms, 1.0)
        # Started clear of the light's own cell. A torch stands INSIDE the cloud, so a march from
        # zero finds the cluster the light is already in, every ray "bounces" where it started, and
        # the bounce is a copy of the source sitting on top of it.
        march = np.linspace(self.cell * 2.0, reach, steps)
        pts = origin[None, None, :] + dirs[:, None, :] * march[None, :, None]
        solid = self.occupancy[self.cells_at(pts)] > 0.0
        # Also drop anything that has marched outside the box: clamping puts those samples on the
        # wall of the grid, where they would report whatever cluster happens to sit against it.
        inside = np.all((pts >= self.lo) & (pts <= self.lo + self.cell * self.resolution), axis=2)
        solid &= inside
        hit = np.argmax(solid, axis=1)
        found = solid[np.arange(len(dirs)), hit]
        return [pts[i, hit[i]] for i in range(len(dirs)) if found[i]]


def build(coords, resolution: int = RESOLUTION, opaque: float = OPAQUE) -> Grid:
    """Bin a cloud into a `Grid`. Cubic cells, so a shadow is the same width in every direction."""
    coords = np.asarray(coords, dtype=float)
    if coords.size == 0:
        return Grid(np.zeros((1, 1, 1)), np.zeros(3), 1.0)
    lo, hi = coords.min(axis=0), coords.max(axis=0)
    # One margin cell all round, so a gene exactly on the boundary is not binned into its neighbour.
    cell = float(max((hi - lo).max(), 1e-6)) / max(resolution - 2, 1)
    lo = lo - cell
    idx = np.floor((coords - lo) / cell).astype(int)
    np.clip(idx, 0, resolution - 1, out=idx)
    counts = np.zeros((resolution, resolution, resolution), dtype=float)
    np.add.at(counts, (idx[:, 0], idx[:, 1], idx[:, 2]), 1.0)
    return Grid(np.clip(counts / opaque, 0.0, 1.0), lo, cell)


def cone(direction, count: int, spread: float = 0.35) -> np.ndarray:
    """`count` directions in a cone around `direction`, spaced evenly rather than at random.

    Evenly, because random directions are a different set on every frame and the bounce would
    flicker between whatever each frame happened to hit. A fixed fan aimed by the pointer moves only
    when the pointer does.
    """
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    d = d / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])
    # Any two axes across the ray. The cross with x fails only when the ray IS x, hence the fallback.
    side = np.cross(d, [1.0, 0.0, 0.0])
    if np.linalg.norm(side) < 1e-6:
        side = np.cross(d, [0.0, 1.0, 0.0])
    side = side / np.linalg.norm(side)
    up = np.cross(d, side)
    out = [d]
    for i in range(1, max(int(count), 1)):
        ang = 2.0 * np.pi * (i - 1) / max(count - 1, 1)
        # The golden ratio spreads the rings out instead of stacking every ray on one circle.
        r = spread * np.sqrt((i * 0.6180339887) % 1.0)
        out.append(d + side * (r * np.cos(ang)) + up * (r * np.sin(ang)))
    return np.array(out[:max(int(count), 1)], dtype=float)
