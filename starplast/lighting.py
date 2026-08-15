#!/usr/bin/env python3
"""Soft interactive lighting and ray-marched shadows for a point cloud.

## What this is, and what it is not

``soft`` uses broad diffuse/specular illumination while ``ray traced`` additionally marches one
shadow ray from every gene to the active light through the point-density volume. This is real
volumetric shadow-ray tracing: occluding clusters attenuate the ray. It is deliberately not called
hardware path tracing--there are no triangle surfaces, reflections or Vulkan renderer here.

The controls describe independent questions: which interaction places the light, which color mood
it has, how points themselves are drawn, and whether shadow rays are traced. The point modes are
few because a rendering choice that cannot be distinguished on the fixed comparison scene is not a
choice worth presenting.
"""
from __future__ import annotations

import math

import numpy as np

#: Ray tracing is a light-transport mode, not a pretend graphics backend.
MODES = ("off", "soft", "ray traced", "deep ray traced")

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
#:
#: `ambient` is the multiplier on the floor, and it is what makes the finishes tell apart at a
#: glance. Highlights alone will not do it, and the measurement is the argument: with the floor held
#: equal, the MEDIAN point differed between matt and glossy by 1 to 3 values out of 255, which is
#: nothing -- a highlight, however bright, lands on the small share of a scatter whose normal
#: happens to face the light, and everything else is unchanged. What separates chalk from a bead is
#: contrast: chalk is flat and evenly lit everywhere, a bead is dark over most of its body with a
#: few bright places. So matt gets a RAISED floor and the shiny finishes a lowered one, which moves
#: every point in the cloud rather than the lucky ones.
#: "2D" is the flat disc this map has always drawn -- no sphere, no highlight, the point's own
#: colour and nothing else. It is a choice rather than a leftover: shading costs contrast, and a
#: reader comparing colours across clusters is better served by a flat scatter than by a lit one.
#: Every other finish draws each gene as a sphere; there is no separate "3D" entry because they are
#: all 3D, and two names for the same picture is a menu that answers a question nobody asked.
#: Three visibly different point renderings. ``flat`` keeps data colors nearly unmodified;
#: ``glossy 3D`` uses a white broad glint; ``metallic 3D`` has a darker body, colored reflection
#: and stronger rim. These names also drive the GPU PBR sphere shader in :mod:`starplast.sprite`;
#: these coefficients remain the CPU fallback's per-point response.
POINT_MODES = {
    "flat": {"ambient": 1.20, "diffuse": 0.25, "specular": 0.0, "shininess": 1.0,
             "rim": 0.0, "tint": 0.0, "size": 1.0},
    "glossy 3D": {"ambient": 0.88, "diffuse": 0.68, "specular": 1.00,
                  "shininess": 7.0, "rim": 0.30, "tint": 0.0, "size": 1.30},
    "metallic 3D": {"ambient": 0.70, "diffuse": 0.38, "specular": 1.45,
                    "shininess": 11.0, "rim": 0.82, "tint": 1.0, "size": 1.35},
    "brushed metal 3D": {"ambient": 0.76, "diffuse": 0.42, "specular": 1.20,
                         "shininess": 7.0, "rim": 0.65, "tint": 0.85, "size": 1.35},
    "silver 3D": {"ambient": 0.72, "diffuse": 0.28, "specular": 1.65,
                  "shininess": 16.0, "rim": 0.92, "tint": 0.35, "size": 1.35},
    "pearl 3D": {"ambient": 0.84, "diffuse": 0.62, "specular": 1.10,
                 "shininess": 10.0, "rim": 0.48, "tint": 0.20, "size": 1.32},
    "glass 3D": {"ambient": 0.78, "diffuse": 0.34, "specular": 1.45,
                 "shininess": 18.0, "rim": 1.05, "tint": 0.05, "size": 1.38},
    "velvet 3D": {"ambient": 0.92, "diffuse": 0.72, "specular": 0.28,
                  "shininess": 3.0, "rim": 0.72, "tint": 0.30, "size": 1.32},
}

# Kept as a code-level alias for callers written before the UI rename. It contains only the new
# modes; old saved values go through ``normalize_point_mode`` below.
FINISHES = POINT_MODES
LEGACY_POINT_MODES = {
    "2D": "flat", "matt": "flat", "satin": "glossy 3D",
    "glossy": "glossy 3D", "metallic": "metallic 3D",
}

#: However dark a finish is allowed to make the body of the cloud. A metal that reads properly as
#: metal is nearly black away from its highlights, and a gene that is nearly black is a gene the
#: reader cannot find -- this is a data display first, so the floor stops there.
MIN_AMBIENT = 0.32
DEFAULT_POINT_MODE = "glossy 3D"
DEFAULT_FINISH = DEFAULT_POINT_MODE

#: Where the light comes from. The default follows the pointer, because the thing a person does with
#: this map is lean into a cluster -- and a light that arrives from wherever they are looking lights
#: the points they are looking at rather than the ones behind them.
SOURCES = ("mouse flashlight", "selected gene", "selected gene and its edges")
DEFAULT_SOURCE = "mouse flashlight"

# Temporary comparison controls for task 37.  They intentionally cover visibly different optical
# models; after direct use the weak ones can be removed without changing the source semantics.
POINTER_MODES = {
    "broad flashlight": {"inner": 18.0, "outer": 34.0, "gain": 1.00},
    "focused flashlight": {"inner": 7.0, "outer": 17.0, "gain": 1.28},
    "soft flashlight": {"inner": 28.0, "outer": 52.0, "gain": 0.78},
    "parallel wash": None,
}
DEFAULT_POINTER_MODE = "broad flashlight"

RESPONSES = {"direct": 1.0, "smooth": 0.24, "cinematic": 0.09}
DEFAULT_RESPONSE = "smooth"

TARGET_MARKERS = ("none", "halo", "beacon", "pulse")
DEFAULT_TARGET_MARKER = "none"

#: Soft color temperatures. Values are intentionally below pure white so specular light cannot
#: bleach categorical colors. The differences survive multiplication by saturated data colors.
LIGHT_MOODS = {
    "neutral": (0.92, 0.92, 0.90),
    "cool blue": (0.66, 0.82, 1.00),
    "warm": (1.00, 0.70, 0.46),
    "daylight": (0.82, 0.90, 1.00),
    "moonlight": (0.42, 0.58, 1.00),
    "gold": (1.00, 0.82, 0.32),
    "rose": (1.00, 0.54, 0.68),
    "violet": (0.72, 0.52, 1.00),
    "laboratory green": (0.50, 1.00, 0.72),
}
DEFAULT_MOOD = "neutral"

#: How many lights, and how fast they travel, at the defaults.
LIGHT_RANGE = (1, 6)
SPEED_RANGE = (0.05, 2.0)
DEFAULT_LIGHTS = 3
DEFAULT_SPEED = 0.35

#: How brightly a light standing IN the cloud lifts what is around it, regardless of facing, and
#: how wide that pool is as a fraction of the map's own radius. This is the whole of "point at a
#: cluster and it lights up" once the cluster is deep in the map. Both numbers were measured rather
#: than picked: a wide pool (0.5) lifts the entire cloud and reads as the map getting brighter --
#: the neighbourhood came out only 1.3x the rest -- while 0.3 gives 6.8x, which reads as a light
#: with somewhere to be. Turned back DOWN from 2.6 once shadows arrived: with nothing blocking it,
#: a pool had to be bright to be found, and a light that can be blocked reads from the contrast
#: between what it reaches and what it does not.
LOCAL_GLOW = 0.62
LOCAL_WIDTH = 0.46

#: How much of a point's own colour survives where no light reaches it. Not zero: an unlit half of
#: the map that went black would hide half the genes, and this is a data display before it is a
#: rendering. 0.35 keeps every point identifiable while still showing the direction of the light.
AMBIENT = 0.55


def normalize_point_mode(value: str | None) -> str:
    """Translate old saved finish names and refuse unknown render modes safely."""
    value = LEGACY_POINT_MODES.get(str(value), str(value))
    return value if value in POINT_MODES else DEFAULT_POINT_MODE


def normalize_source(value: str | None) -> str:
    """Translate the old ``mouse`` label and discard removed orbit/corner sources."""
    value = "mouse flashlight" if str(value) in ("mouse", "mouse proximity") else str(value)
    return value if value in SOURCES else DEFAULT_SOURCE


def mood_color(mood: str | None) -> np.ndarray:
    """The RGB illumination for a named mood, with neutral as the safe fallback."""
    return np.asarray(LIGHT_MOODS.get(str(mood), LIGHT_MOODS[DEFAULT_MOOD]), dtype=float)


def lights(t: float, n: int = DEFAULT_LIGHTS, speed: float = DEFAULT_SPEED, radius: float = 90.0,
           colors=None) -> list:
    """Where each light is at time `t`, and what colour it is.

    Lissajous paths with incommensurable frequencies, so the lights never fall into a repeating
    formation -- three lights orbiting in step read as one light, and the point of several is that
    the shading changes as they pass each other. Position is a function of the clock for the same
    reason the blobs are: any frame can be drawn at any time, and a dropped frame costs nothing.

    The curve gives the DIRECTION and the radius is fixed, which it was not at first. Taken as a
    position, a Lissajous figure ranges over the whole cube it is inscribed in: measured, the lights
    swung between 0.2 and 1.7 times the intended radius and spent about 4% of frames inside the
    cloud itself. That reads as the map pulsing -- overall brightness pumping up and down, and every
    so often a light surfacing among the genes to blow out its neighbours and leave the rest dark.
    Which is what "the orbit looks strange" was describing. On a shell, only the direction changes.
    """
    n = int(min(max(n, LIGHT_RANGE[0]), LIGHT_RANGE[1]))
    colors = colors or [(1.0, 0.95, 0.85), (0.75, 0.85, 1.0), (1.0, 0.8, 0.9),
                        (0.85, 1.0, 0.9), (1.0, 0.9, 0.7), (0.9, 0.85, 1.0)]
    out = []
    for i in range(n):
        a, b, c = 1.0 + 0.31 * i, 1.37 + 0.19 * i, 0.71 + 0.23 * i
        phase = 2.0 * math.pi * i / n
        tt = t * speed
        where = np.array([math.sin(a * tt + phase),
                          math.sin(b * tt + phase * 1.3),
                          math.cos(c * tt + phase * 0.7)], dtype=float)
        length = float(np.linalg.norm(where))
        # Near the origin the direction is ill-conditioned rather than undefined -- the three sines
        # are all crossing zero at once. Straight up is as good an answer as any and happens for an
        # instant in passing.
        where = where / length if length > 1e-6 else np.array([0.0, 0.0, 1.0])
        out.append({
            "pos": where * radius,
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
          shininess: float = 24.0, finish: str = None, point_mode: str = None,
          eye=None) -> np.ndarray:
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

    requested = point_mode if point_mode is not None else finish
    normalized = normalize_point_mode(requested) if requested is not None else None
    f = POINT_MODES.get(normalized, None) if normalized is not None else None
    if f is not None:
        specular, shininess = f["specular"] > 0, f["shininess"]
    diffuse_gain = f["diffuse"] if f else 1.0
    spec_gain = f["specular"] if f else 1.0
    rim_gain = f["rim"] if f else 0.0
    tint = f["tint"] if f else 0.0
    if f is not None:
        ambient = max(ambient * f["ambient"], MIN_AMBIENT)

    # Diffuse and specular are kept APART, and this is not a detail. Multiplying the point's colour
    # by everything -- which is what this did first -- means a highlight can never whiten a coloured
    # point: a pure red gene has no green to raise, so a white light glinting off it stayed red and
    # every finish looked the same on saturated colours. Diffuse light multiplies the colour, because
    # that is light the surface absorbed and re-emitted; specular is ADDED, because that is light
    # that bounced off without ever being coloured by it. The exception is metals, which have no
    # separate diffuse colour and tint the bounce itself -- which is the whole difference between a
    # copper bead and a white-glinting plastic one.
    own = np.clip(rgba[:, :3], 0.0, 1.0)
    # The mood is a broad fill, not just the few points whose pseudo-normal faces the lamp. With a
    # neutral ambient floor, cool/warm changed only the small highlight share and measured 1.5/255
    # at the median -- a setting that existed in the menu but not in the picture. Tinting the fill
    # preserves the data hue (it still multiplies `own`) while making the temperature readable.
    if lit:
        fill_source = np.mean([np.asarray(light["color"], dtype=float) for light in lit], axis=0)
        fill_color = 0.55 + 0.45 * np.clip(fill_source, 0.0, 1.0)
    else:
        fill_color = np.ones(3)
    total = np.tile(ambient * fill_color, (len(coords), 1))
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
        if light.get("spot"):
            from_light = -direction
            axis = np.asarray(light["direction"], dtype=float)
            axis /= max(float(np.linalg.norm(axis)), 1e-9)
            cosine = (from_light * axis).sum(axis=1, keepdims=True)
            inner = math.cos(math.radians(float(light["inner"])))
            outer = math.cos(math.radians(float(light["outer"])))
            blend = np.clip((cosine - outer) / max(inner - outer, 1e-6), 0.0, 1.0)
            beam = blend * blend * (3.0 - 2.0 * blend)
            falloff *= beam * float(light.get("gain", 1.0))
        # How much of this light survives the trip -- 1 where the line to the gene is clear, less
        # where the map is in the way. Attached by whoever built the light (see `rays`), because
        # working it out needs the whole cloud binned and this function shades one frame.
        falloff = falloff * light.get("shadow", 1.0)
        diffuse = np.clip((normal * direction).sum(axis=1, keepdims=True), 0.0, None)
        total += (diffuse * falloff * diffuse_gain) * light["color"]
        if light.get("local"):
            # A light INSIDE the cloud also brightens whatever is near it, whichever way that point
            # is facing. Without this a torch cannot light the far side of the map at all: the
            # pseudo-normal there points away from the viewer, the light hangs between the viewer
            # and the cluster, and so it lands on the back of every gene it is meant to be lighting.
            # Measured, pointing into a far cluster left it DIMMER than the near face of the map.
            # Physically this is a lamp in a scattering medium rather than a lamp in a vacuum, which
            # is the better model for a cloud of points that have no surfaces anyway.
            glow = LOCAL_GLOW * light.get("gain", 1.0) / (1.0 + (dist / (LOCAL_WIDTH * scale)) ** 2)
            total += (glow * light.get("shadow", 1.0)) * light["color"]
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


def torch(anchor, basis, radius: float, color=(1.0, 0.97, 0.92), stand_off: float = 0.08) -> list:
    """A light hanging just in front of a given gene, between it and the viewer.

    This is what "point at it and it lights up" has to mean once the map has depth. A light placed
    outside the cloud on the viewer's side can only ever light the near face -- it is behind
    everything else -- so pointing into a cluster deep in the map lit the front of the map instead,
    which is what was reported. Anchored on the gene under the pointer, the light is IN the cloud at
    that depth and its falloff does the rest: the genes around it are lit whether they are at the
    front or the back.

    Held off the gene rather than sitting on it, because a light exactly on a point makes that point
    a white dot and its neighbours a hard black shell. The stand-off is a fraction of the map's own
    radius, so the pool of light stays the same size relative to the data at any zoom.
    """
    anchor = np.asarray(anchor, dtype=float)
    eye = np.asarray(basis[0], dtype=float)
    away = eye - anchor
    n = float(np.linalg.norm(away))
    away = away / n if n > 1e-9 else -np.asarray(basis[3], dtype=float)
    return [{"pos": anchor + away * radius * stand_off,
             "color": np.array(color, dtype=float), "local": True}]


def flashlight(pointer, basis, centre, radius: float, cone: dict,
               color=(1.0, 0.97, 0.92), fov: float = 60.0, aspect: float = 1.0) -> list:
    """A camera-origin spotlight following continuous screen coordinates.

    No gene is queried here.  The cursor defines one ray through the camera frustum, so crossing two
    overlapping points at different depths does not change either the origin or direction.  The
    light target is the intersection with the plane through the cloud center, used only for drawing
    an optional marker.
    """
    eye, right, up, forward = (np.asarray(v, dtype=float) for v in basis)
    x, y, _ = pointer if pointer is not None else (0.0, 0.0, 1.0)
    spread = math.tan(math.radians(float(fov)) * 0.5)
    direction = forward + right * float(x) * spread + up * float(y) * spread / max(aspect, 1e-6)
    direction /= max(float(np.linalg.norm(direction)), 1e-9)
    depth = float(np.dot(np.asarray(centre, dtype=float) - eye, forward))
    target = eye + direction * depth / max(float(np.dot(direction, forward)), 1e-6)
    return [{"pos": eye, "color": np.array(color, dtype=float), "spot": True,
             "direction": direction, "inner": float(cone["inner"]),
             "outer": float(cone["outer"]), "gain": float(cone["gain"]),
             "target": target}]


def light_at(coords, source: str, t: float, n: int, speed: float, radius: float,
             pointer=None, basis=None, selected=None, neighbours=None, anchor=None,
             color=None, pointer_mode: str = DEFAULT_POINTER_MODE, fov: float = 60.0,
             aspect: float = 1.0) -> list:
    """The lights for one frame, for whichever source was chosen.

    The screen-relative sources need `basis` -- the camera's own axes -- because "top left" and
    "where the pointer is" are statements about the picture, not about the data. `anchor` is the
    ``anchor`` remains accepted for saved callers but is deliberately ignored: mouse illumination
    is a continuous camera ray, never the identity or depth of the nearest gene.
    """
    coords = np.asarray(coords, dtype=float)
    source = normalize_source(source)
    color = mood_color(DEFAULT_MOOD) if color is None else np.asarray(color, dtype=float)
    centre = coords.mean(axis=0) if len(coords) else np.zeros(3)
    if source == "selected gene" and selected is not None:
        return at_points(coords, [selected], radius, color=color)
    if source == "selected gene and its edges" and selected is not None:
        return at_points(coords, [selected] + list(neighbours or []), radius, color=color)
    if basis is not None:
        if source == "mouse flashlight":
            cone = POINTER_MODES.get(str(pointer_mode), POINTER_MODES[DEFAULT_POINTER_MODE])
            if cone is not None:
                return flashlight(pointer, basis, centre, radius, cone, color=color,
                                  fov=fov, aspect=aspect)
            x, y, _ = pointer if pointer is not None else (0.0, 0.0, 1.0)
            return on_screen((x * 1.2, y * 1.2, 1.0), basis, centre, radius, color=color)
        return on_screen((0.0, 0.0, 1.0), basis, centre, radius, color=color)
    # No camera to speak of -- offscreen, or before the first frame. A light from +Z lights what is
    # facing the reader, which is the honest default rather than a guess.
    return fixed((0.0, 0.0, 1.0), radius, color=color)
