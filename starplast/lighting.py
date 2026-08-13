#!/usr/bin/env python3
"""Moving lights, and points shaded by where the light is coming from.

## What this is, and what it is not

It is **per-point shading**: a few lights move on their own paths, and every gene's colour is its
own colour modulated by how much light reaches it from those directions. Diffuse, specular and
distance falloff, computed for all 8,140 points as array arithmetic on every frame.

It is **not ray tracing**, and the difference is not a detail. Ray tracing answers "what does this
pixel see, following the light backwards through the scene" -- which buys reflections, refraction,
and shadows cast by one object onto another. None of those exist here, because there are no
surfaces: a point cloud has no geometry to occlude anything. Even given surfaces, this renders
through pyqtgraph's GL scatter, which draws sprites through a fixed pipeline; a ray tracer would
mean replacing the renderer with OptiX or a Vulkan RT pipeline, not adding an option to this one.

**Light Propagation Volumes and Voxel Cone Tracing** are the same answer twice. Both approximate
indirect light by voxelising the scene into a grid and propagating or cone-sampling radiance through
it. Both need geometry to voxelise and a compute pipeline to propagate; a scatter of unconnected
points voxelises to a sparse cloud of isolated cells, and what came back would be a blur of each
point's own colour -- an expensive way to reproduce the ambient term already in the equation below.

So what is offered is what actually changes the picture: direction. A point on the side of the map
facing a light is bright, one facing away is dim, and as the lights move the shape of the cloud
reads in a way that a flat colour per gene does not show. The normal is the direction from the
centre of the map outward, which is the standard trick for a point cloud with no surface -- it
shades the cloud as though it were a solid body, which is exactly the shape a reader is trying to see.
"""
from __future__ import annotations

import math

import numpy as np

#: The modes offered, in the order a menu should list them.
MODES = ("off", "lit", "lit + specular")

#: How many lights, and how fast they travel, at the defaults.
LIGHT_RANGE = (1, 6)
SPEED_RANGE = (0.05, 2.0)
DEFAULT_LIGHTS = 3
DEFAULT_SPEED = 0.35

#: How much of a point's own colour survives where no light reaches it. Not zero: an unlit half of
#: the map that went black would hide half the genes, and this is a data display before it is a
#: rendering. 0.35 keeps every point identifiable while still showing the direction of the light.
AMBIENT = 0.35


def lights(t: float, n: int = DEFAULT_LIGHTS, speed: float = DEFAULT_SPEED, radius: float = 90.0,
           colors=None) -> list:
    """Where each light is at time `t`, and what colour it is.

    Lissajous paths with incommensurable frequencies, so the lights never fall into a repeating
    formation -- three lights orbiting in step read as one light, and the point of several is that
    the shading changes as they pass each other. Position is a function of the clock for the same
    reason the blobs are: any frame can be drawn at any time, and a dropped frame costs nothing.
    """
    n = int(min(max(n, LIGHT_RANGE[0]), LIGHT_RANGE[1]))
    colors = colors or [(1.0, 0.95, 0.85), (0.75, 0.85, 1.0), (1.0, 0.8, 0.9),
                        (0.85, 1.0, 0.9), (1.0, 0.9, 0.7), (0.9, 0.85, 1.0)]
    out = []
    for i in range(n):
        a, b, c = 1.0 + 0.31 * i, 1.37 + 0.19 * i, 0.71 + 0.23 * i
        phase = 2.0 * math.pi * i / n
        tt = t * speed
        out.append({
            "pos": np.array([radius * math.sin(a * tt + phase),
                             radius * math.sin(b * tt + phase * 1.3),
                             radius * math.cos(c * tt + phase * 0.7)], dtype=float),
            "color": np.array(colors[i % len(colors)], dtype=float),
        })
    return out


def shade(coords, colors, lit, specular: bool = False, ambient: float = AMBIENT,
          shininess: float = 24.0) -> np.ndarray:
    """Point colours under a set of lights. Returns RGBA in the shape it was given.

    The normal is the outward direction from the centre of the cloud, because a point has no surface
    of its own. Distance falloff is gentle -- inverse linear rather than inverse square -- since the
    map is not a physical scene and a physically correct falloff makes everything past the first
    light black.
    """
    coords = np.asarray(coords, dtype=float)
    rgba = np.array(colors, dtype=float, copy=True)
    if coords.size == 0 or rgba.size == 0:
        return rgba
    centre = coords.mean(axis=0)
    normal = coords - centre
    length = np.linalg.norm(normal, axis=1, keepdims=True)
    normal = np.divide(normal, np.where(length > 1e-9, length, 1.0))
    scale = float(np.percentile(length, 95)) or 1.0

    total = np.full((len(coords), 3), ambient, dtype=float)
    # The viewer is treated as far away on +Z. A specular term that tracked the real camera would
    # move the highlight when the map is rotated, which reads as the data changing rather than the
    # view -- and rotating to look at a cluster is the commonest thing anyone does here.
    view = np.array([0.0, 0.0, 1.0])
    for light in lit:
        to_light = light["pos"] - coords
        dist = np.linalg.norm(to_light, axis=1, keepdims=True)
        direction = np.divide(to_light, np.where(dist > 1e-9, dist, 1.0))
        falloff = 1.0 / (1.0 + dist / (2.0 * scale))
        diffuse = np.clip((normal * direction).sum(axis=1, keepdims=True), 0.0, None)
        total += (diffuse * falloff) * light["color"]
        if specular:
            half = direction + view
            half /= np.maximum(np.linalg.norm(half, axis=1, keepdims=True), 1e-9)
            spec = np.clip((normal * half).sum(axis=1, keepdims=True), 0.0, None) ** shininess
            total += (spec * falloff) * light["color"]

    rgba[:, :3] = np.clip(rgba[:, :3] * np.clip(total, 0.0, 2.0), 0.0, 1.0)
    return rgba
