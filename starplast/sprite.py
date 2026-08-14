#!/usr/bin/env python3
"""The picture inside one gene, which is what makes a ball look like a ball.

## Why a texture rather than a shader

pyqtgraph draws a scatter as point sprites: one square of texture per gene, always facing the
viewer, with a vertex shader that does nothing but set the size and pass the colour through. There
is no fragment shader. The texture it ships is `pData[:] = 255` with an alpha disc cut out of it --
every pixel inside a ball carries the identical colour, which is a flat disc, and no amount of
per-point shading can make a flat disc look spherical. That is the whole of "the balls still look
like 2D matt balls": the per-point lighting was working and had nowhere to show up, because a ball
is one colour and a sphere is a gradient.

The fixed-function stage MODULATES: what you see is texture RGBA times the point's own colour. So a
texture holding a shaded unit sphere -- bright where the surface turns toward the light, dark at the
limb, a highlight where it should be -- multiplies each gene's colour into a lit ball of that
colour. That is a real 3D read for the cost of a 64x64 array, and it composes with the per-point
shading already there: the sprite says where the light falls ON a gene, the point colour says how
much light reaches it.

The direction the highlight comes from is the scene's own light, converted to the screen's frame,
so the balls catch the light from wherever the lights are. It is one texture for the whole cloud,
which is what a distant light gives anyway.

## What is not here

Each ball is lit as an isolated sphere. No gene shadows another, and a ball's highlight does not
know what is behind it -- the same limit as everywhere else in this renderer, for the same reason
(see `lighting`). What it buys is the shape of a sphere, which is the thing that was missing.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph.opengl as gl
from OpenGL import GL

#: How each finish answers the light ACROSS ONE BALL. Separate numbers from `lighting.FINISHES`,
#: which says how the cloud answers it, and they are separate questions: a tight highlight is
#: useless spread over a scatter of points and is exactly right inside a single sphere, where there
#: are pixels for it to land on. Shininess here is the textbook range for that reason.
#:
#: "2D" is the plain flat disc -- pyqtgraph's own sprite, and what this map has always drawn. Kept
#: as a choice because a flat disc is the honest way to show a scatter where the reader is comparing
#: colours and a highlight is one more thing in the way.
SPRITES = {
    "2D": None,
    "matt": {"ambient": 0.34, "diffuse": 0.92, "specular": 0.0, "shininess": 1.0, "rim": 0.05},
    "satin": {"ambient": 0.28, "diffuse": 0.88, "specular": 0.35, "shininess": 14.0, "rim": 0.12},
    "glossy": {"ambient": 0.20, "diffuse": 0.80, "specular": 1.0, "shininess": 52.0, "rim": 0.22},
    "metallic": {"ambient": 0.10, "diffuse": 0.52, "specular": 1.25, "shininess": 24.0, "rim": 0.55},
}

#: Texels across one ball. 64 is what pyqtgraph uses, and a gene is drawn at 5 to 12 pixels: the
#: texture is downsampled hard, so a bigger one buys nothing and a smaller one shows its own edge.
WIDTH = 64

#: Where the light comes from when nothing else says -- above and to the left, which is where every
#: reader since the Renaissance has assumed light comes from, and reads as convex rather than
#: hollow. z is toward the viewer.
DEFAULT_LIGHT = (-0.4, 0.6, 0.7)


def disc_alpha(w: int = WIDTH) -> np.ndarray:
    """The round edge, antialiased over the last texel. pyqtgraph's own formula.

    Matched deliberately: a ball with a different edge to the halos and centroids drawn beside it
    reads as a different KIND of thing rather than the same thing lit.
    """
    y, x = np.mgrid[0:w, 0:w]
    r = np.hypot(x - (w - 1) / 2.0, y - (w - 1) / 2.0)
    return np.clip(w / 2 - np.clip(r, w / 2 - 1, w / 2), 0.0, 1.0)


def texture(finish: str, light=DEFAULT_LIGHT, w: int = WIDTH) -> np.ndarray:
    """One ball, lit from `light`, as (w, w, 4) uint8 ready for GL.

    `light` is a direction in the screen's frame: x right, y up, z toward the viewer. The surface
    normal is read off the disc -- a texel at (u, v) from the centre of a unit sphere has normal
    (u, v, sqrt(1 - u^2 - v^2)), which is the whole trick and costs one square root per texel.
    """
    alpha = disc_alpha(w)
    f = SPRITES.get(finish)
    if f is None:                       # "2D", or a name nobody defined: the flat disc
        out = np.empty((w, w, 4), dtype=np.ubyte)
        out[:, :, :3] = 255
        out[:, :, 3] = (alpha * 255).astype(np.ubyte)
        return out

    span = np.linspace(-1.0, 1.0, w)
    u, v = np.meshgrid(span, -span)     # -span: texture rows run down, the sphere's y runs up
    z2 = 1.0 - u * u - v * v
    inside = z2 > 0.0
    nz = np.sqrt(np.where(inside, z2, 0.0))

    d = np.asarray(light, dtype=float)
    n = np.linalg.norm(d)
    if n <= 1e-9:                       # a light with no direction: use the one above and left
        d = np.asarray(DEFAULT_LIGHT, dtype=float)
        n = np.linalg.norm(d)
    d = d / n
    lam = np.clip(u * d[0] + v * d[1] + nz * d[2], 0.0, None)

    # Blinn-Phong: the halfway vector between the light and the viewer, who is straight ahead at
    # (0, 0, 1) because a point sprite always faces them.
    half = d + np.array([0.0, 0.0, 1.0])
    half = half / max(float(np.linalg.norm(half)), 1e-9)
    spec = np.clip(u * half[0] + v * half[1] + nz * half[2], 0.0, None) ** f["shininess"]

    # The limb, where the surface turns away from the viewer. On a real sphere this is where a
    # glancing reflection of everything else in the room shows up; here it is what stops a ball
    # ending in a hard circular edge, and it is most of what a metal looks like.
    rim = (1.0 - nz) ** 3 * np.clip(lam + 0.35, 0.0, 1.0)

    shade = f["ambient"] + f["diffuse"] * lam + f["specular"] * spec + f["rim"] * rim
    shade = np.clip(np.where(inside, shade, f["ambient"]), 0.0, 1.0)

    out = np.empty((w, w, 4), dtype=np.ubyte)
    out[:, :, 0] = out[:, :, 1] = out[:, :, 2] = (shade * 255).astype(np.ubyte)
    out[:, :, 3] = (alpha * 255).astype(np.ubyte)
    return out


def upload(texture_id, data) -> None:
    """Hand a texture to the card. Only ever called with a current context -- see `ShadedScatter`."""
    GL.glBindTexture(GL.GL_TEXTURE_2D, texture_id)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, data.shape[0], data.shape[1], 0,
                    GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, data)


class ShadedScatter(gl.GLScatterPlotItem):
    """A scatter whose sprite is a lit sphere instead of a flat disc.

    The texture is set from outside, on any thread of control -- a finish changing, a light moving --
    and uploaded at the next paint, because a texture can only be uploaded while a GL context is
    current and neither of those events has one. Holding it as `_pending` is what keeps
    "which sprite" a question about the model rather than about the paint clock.
    """

    def __init__(self, **kwds):
        super().__init__(**kwds)
        self._pending = None

    def set_sprite(self, data) -> None:
        """Draw every point with this ball from the next paint onwards."""
        self._pending = data
        self.update()

    def _flush(self) -> bool:
        if self._pending is None or getattr(self, "pointTexture", None) is None:
            return False
        upload(self.pointTexture, self._pending)
        self._pending = None
        return True

    def initializeGL(self):
        """Build the stock sprite, then replace it with ours if one is waiting."""
        super().initializeGL()
        self._flush()

    def paint(self):
        """One frame. The only place with a context current, so the only place a texture can go up."""
        self._flush()
        super().paint()


def to_screen(direction, basis) -> tuple:
    """A world direction as the screen sees it: x right, y up, z toward the viewer.

    The sprite is drawn in the plane of the screen, so a light in world coordinates means nothing to
    it until it is expressed this way -- and doing it every frame is what makes the highlight on the
    balls travel with the lights rather than sitting where it was when they were built.
    """
    _, right, up, forward = basis
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-9:
        return DEFAULT_LIGHT
    d = d / n
    return (float(d @ np.asarray(right, dtype=float)),
            float(d @ np.asarray(up, dtype=float)),
            float(-(d @ np.asarray(forward, dtype=float))))
