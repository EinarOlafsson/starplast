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

#: On or off. What KIND of lit is `FINISHES`, and where the light comes from is `SOURCES` -- three
#: separate questions that used to be one dropdown reading "lit + specular".
MODES = ("off", "lit")

#: How a point's surface answers the light. The parameters are what separate a mineral from a
#: billiard ball: how much of the light scatters (diffuse), how much bounces (specular), how tightly
#: (shininess), how much the silhouette catches the light (rim), and whether the bounce takes the
#: point's own colour (metals) or the light's (dielectrics -- plastic, chalk, skin).
#:
#: The shininess numbers are LOW for a reason worth stating, because the textbook values are not.
#: A Phong lobe of 48 -- ordinary for a glossy solid -- is about four degrees wide, and these normals
#: are not surface normals: they are the outward radial direction, one per gene, on a cloud of a few
#: thousand scattered points. A four-degree lobe on a scatter that sparse lands on one or two points
#: out of thousands, which is invisible; measured, every finish came out matt, which is exactly what
#: was reported. Broad lobes plus a rim term is what reads at this sampling: the rim brightens every
#: point whose normal turns away from the viewer, so the finish shows up along the whole silhouette
#: of the cloud rather than in a highlight too small to find.
#: The diffuse gains stay near 1. Dropping diffuse as gloss rises is what a physical shader does --
#: energy that bounces off did not scatter -- but here it cancelled the effect exactly: glossy lost
#: as much broad light as its highlight added, so the mean barely moved and the finish was invisible.
#: These finishes differ by what they ADD, which is the part a reader can see.
FINISHES = {
    "matt": {"diffuse": 1.0, "specular": 0.0, "shininess": 1.0, "rim": 0.0, "tint": 0.0},
    "satin": {"diffuse": 1.0, "specular": 0.45, "shininess": 5.0, "rim": 0.18, "tint": 0.0},
    "glossy": {"diffuse": 0.95, "specular": 1.1, "shininess": 12.0, "rim": 0.45, "tint": 0.0},
    "metallic": {"diffuse": 0.55, "specular": 1.3, "shininess": 7.0, "rim": 0.70, "tint": 1.0},
}
DEFAULT_FINISH = "satin"

#: Where the light comes from. The default follows the pointer, because the thing a person does with
#: this map is lean into a cluster -- and a light that arrives from wherever they are looking lights
#: the points they are looking at rather than the ones behind them.
SOURCES = ("mouse", "orbiting", "top left", "top right", "bottom left", "bottom right",
           "selected gene", "selected gene and its edges")
DEFAULT_SOURCE = "mouse"

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


def fixed(direction, radius: float, color=(1.0, 0.97, 0.92)) -> list:
    """One light, far away in a given direction -- the corner presets and the pointer both use this.

    Far away rather than at the point: a light placed among the genes lights the near side of a few
    and leaves the rest black, which reads as data missing rather than as a light.
    """
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    d = d / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])
    return [{"pos": d * radius * 2.0, "color": np.array(color, dtype=float)}]


#: The four corner presets, as directions in the view's own frame: x right, y up, z toward the eye.
CORNERS = {
    "top left": (-1.0, 1.0, 0.8), "top right": (1.0, 1.0, 0.8),
    "bottom left": (-1.0, -1.0, 0.8), "bottom right": (1.0, -1.0, 0.8),
}


def at_points(coords, indices, radius: float, color=(1.0, 0.95, 0.85)) -> list:
    """A light sitting on each of the given genes.

    Used for "light the gene I clicked" and "light it and everything it is connected to", where the
    light IS the answer to a question about those genes -- so it goes where they are rather than
    outside the cloud, and its falloff does the rest.
    """
    coords = np.asarray(coords, dtype=float)
    out = []
    for i in np.atleast_1d(np.asarray(indices, dtype=int)):
        if 0 <= int(i) < len(coords):
            out.append({"pos": coords[int(i)].astype(float),
                        "color": np.array(color, dtype=float), "local": True})
    return out or fixed((0.0, 0.0, 1.0), radius)


def shade(coords, colors, lit, specular: bool = False, ambient: float = AMBIENT,
          shininess: float = 24.0, finish: str = None, eye=None) -> np.ndarray:
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

    f = FINISHES.get(finish or "", None)
    if f is not None:
        specular, shininess = f["specular"] > 0, f["shininess"]
    diffuse_gain = f["diffuse"] if f else 1.0
    spec_gain = f["specular"] if f else 1.0
    rim_gain = f["rim"] if f else 0.0
    tint = f["tint"] if f else 0.0

    # Diffuse and specular are kept APART, and this is not a detail. Multiplying the point's colour
    # by everything -- which is what this did first -- means a highlight can never whiten a coloured
    # point: a pure red gene has no green to raise, so a white light glinting off it stayed red and
    # every finish looked the same on saturated colours. Diffuse light multiplies the colour, because
    # that is light the surface absorbed and re-emitted; specular is ADDED, because that is light
    # that bounced off without ever being coloured by it. The exception is metals, which have no
    # separate diffuse colour and tint the bounce itself -- which is the whole difference between a
    # copper bead and a white-glinting plastic one.
    own = np.clip(rgba[:, :3], 0.0, 1.0)
    total = np.full((len(coords), 3), ambient, dtype=float)
    highlight = np.zeros((len(coords), 3), dtype=float)
    # The viewer is treated as far away on +Z. A specular term that tracked the real camera would
    # move the highlight when the map is rotated, which reads as the data changing rather than the
    # view -- and rotating to look at a cluster is the commonest thing anyone does here.
    # The REAL camera where there is one. A specular term computed against a fixed +Z never lands
    # where the viewer is looking once the map has been orbited, so the highlight sits on the far
    # side of the cloud and every finish looks matt from the front -- which is what was reported.
    if eye is not None:
        to_eye = np.asarray(eye, dtype=float) - coords
        view = to_eye / np.maximum(np.linalg.norm(to_eye, axis=1, keepdims=True), 1e-9)
    else:
        view = np.array([0.0, 0.0, 1.0])
    # How far each point's normal has turned away from the viewer. 1 on the silhouette, 0 dead
    # centre. Squared so it stays confined to the edge instead of washing the whole cloud.
    facing = np.abs((normal * view).sum(axis=1, keepdims=True)) if rim_gain else None
    fresnel = (1.0 - np.clip(facing, 0.0, 1.0)) ** 2 if rim_gain else None
    for light in lit:
        to_light = light["pos"] - coords
        dist = np.linalg.norm(to_light, axis=1, keepdims=True)
        direction = np.divide(to_light, np.where(dist > 1e-9, dist, 1.0))
        falloff = 1.0 / (1.0 + dist / (2.0 * scale))
        diffuse = np.clip((normal * direction).sum(axis=1, keepdims=True), 0.0, None)
        total += (diffuse * falloff * diffuse_gain) * light["color"]
        if specular:
            half = direction + view
            half /= np.maximum(np.linalg.norm(half, axis=1, keepdims=True), 1e-9)
            spec = np.clip((normal * half).sum(axis=1, keepdims=True), 0.0, None) ** shininess
            hue = light["color"] * (1.0 - tint) + own * tint
            highlight += (spec * falloff * spec_gain) * hue
            if rim_gain:
                # Tied to the light rather than free-standing: a rim that glows on the side no
                # light reaches would be inventing brightness, and this is a data display. The
                # 0.3 floor keeps the silhouette visible in grazing light instead of switching
                # off the moment a point turns edge-on to the source.
                highlight += (fresnel * falloff * rim_gain * (0.3 + 0.7 * diffuse)) * hue

    rgba[:, :3] = np.clip(own * np.clip(total, 0.0, 2.0) + highlight, 0.0, 1.0)
    return rgba


def on_screen(where, basis, centre, radius: float, color=(1.0, 0.97, 0.92)) -> list:
    """A light placed relative to the CAMERA rather than to the data.

    `where` is (x, y, z) in the screen's own frame: x right, y up, z toward the viewer. `basis` is
    (eye, right, up, forward) from the renderer. This is what makes "top left" mean the top left of
    the picture and keep meaning it while the map is orbited -- placed in world coordinates instead,
    a corner light drifts to the other side of the cloud as soon as anything moves, which looks like
    a light that does not work rather than one that is somewhere else.
    """
    eye, right, up, forward = basis
    x, y, z = (float(v) for v in where)
    direction = right * x + up * y - forward * z
    n = np.linalg.norm(direction)
    direction = direction / n if n > 1e-9 else np.asarray(forward, dtype=float)
    return [{"pos": np.asarray(centre, dtype=float) + direction * radius * 2.0,
             "color": np.array(color, dtype=float)}]


def light_at(coords, source: str, t: float, n: int, speed: float, radius: float,
             pointer=None, basis=None, selected=None, neighbours=None) -> list:
    """The lights for one frame, for whichever source was chosen.

    The screen-relative sources need `basis` -- the camera's own axes -- because "top left" and
    "where the pointer is" are statements about the picture, not about the data.
    """
    coords = np.asarray(coords, dtype=float)
    centre = coords.mean(axis=0) if len(coords) else np.zeros(3)
    if source == "orbiting":
        return lights(t, n, speed, radius=radius)
    if source == "selected gene" and selected is not None:
        return at_points(coords, [selected], radius)
    if source == "selected gene and its edges" and selected is not None:
        return at_points(coords, [selected] + list(neighbours or []), radius)
    if basis is not None:
        if source in CORNERS:
            return on_screen(CORNERS[source], basis, centre, radius)
        if source == "mouse" and pointer is not None:
            # Pushed forward of the screen plane so the light is between the viewer and the cloud
            # rather than in it: at z=0 half the map is behind the light and goes dark.
            x, y, _ = pointer
            return on_screen((x * 1.4, y * 1.4, 1.0), basis, centre, radius)
        return on_screen((0.0, 0.0, 1.0), basis, centre, radius)
    # No camera to speak of -- offscreen, or before the first frame. A light from +Z lights what is
    # facing the reader, which is the honest default rather than a guess.
    return fixed((0.0, 0.0, 1.0), radius)
