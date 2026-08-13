#!/usr/bin/env python3
"""A drifting field of soft colour blobs, behind the panels.

Ported from spaCR's ambient background so the two programs look like one pair of tools, and kept to
the one animation that was asked for -- the blobs -- with the same controls under the same names:
speed, size, density, blur.

The design decision worth carrying across: **position is a pure function of the clock**, two
independent sines per blob rather than a random walk. A walk accumulates error, cannot be rewound,
and gives a different picture every time a frame is dropped; a function of time can be asked for any
instant, which is what makes the field paintable at whatever rate the window can manage and
identical at the same `t` on any machine.

Blobs are seeded on a jittered grid rather than uniformly at random. With a dozen of them, uniform
sampling reliably leaves one corner empty and puts three in the middle -- which reads as a bug in
the drawing rather than as a random arrangement.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

from PyQt6 import QtCore, QtGui, QtWidgets

#: The controls, and the range each is allowed -- the same names and bounds spaCR uses, so a person
#: who has set this up once does not have to learn it twice.
SPEED_RANGE = (0.1, 4.0)
SIZE_RANGE = (0.25, 2.5)
DENSITY_RANGE = (0.25, 3.0)
BLUR_RANGE = (0.0, 3.0)
DEFAULT_SPEED = 1.0
DEFAULT_SIZE = 1.0
DEFAULT_DENSITY = 1.0
DEFAULT_BLUR = 0.0

#: How many blobs at density 1.0, and the grid they are seeded on.
BLOB_COUNT = 12
GRID = (4, 3)


@dataclass
class Blob:
    """One blob: where it sits, how far and how fast it wanders, how it breathes."""
    x: float
    y: float
    drift_x: float
    drift_y: float
    rate_x: float
    rate_y: float
    phase_x: float
    phase_y: float
    radius: float
    pulse: float
    pulse_rate: float
    pulse_phase: float
    color: int


def field(seed: int = 0, density: float = DEFAULT_DENSITY) -> list:
    """A reproducible set of blobs, seeded on a jittered grid.

    The pool is always rolled at the top of the density range and a prefix of it painted, so raising
    the density adds blobs rather than redrawing the ones already there: the random numbers are
    consumed in the same order either way, and blob 3 keeps the numbers blob 3 always had.
    """
    rng = random.Random(seed)
    cols, rows = GRID
    cells = list(range(cols * rows))
    rng.shuffle(cells)
    pool = int(round(BLOB_COUNT * DENSITY_RANGE[1]))
    blobs = []
    for i in range(pool):
        cell = cells[i % len(cells)]
        col, row = cell % cols, cell // cols
        small = (i % 3) == 0
        lo, hi = (0.10, 0.18) if small else (0.22, 0.38)
        blobs.append(Blob(
            x=(col + 0.15 + 0.7 * rng.random()) / cols,
            y=(row + 0.15 + 0.7 * rng.random()) / rows,
            drift_x=rng.uniform(0.02, 0.09), drift_y=rng.uniform(0.02, 0.09),
            rate_x=2 * math.pi / rng.uniform(19.0, 43.0),
            rate_y=2 * math.pi / rng.uniform(23.0, 51.0),
            phase_x=rng.uniform(0.0, 2 * math.pi), phase_y=rng.uniform(0.0, 2 * math.pi),
            radius=rng.uniform(lo, hi), pulse=rng.uniform(0.05, 0.16),
            pulse_rate=2 * math.pi / rng.uniform(11.0, 29.0),
            pulse_phase=rng.uniform(0.0, 2 * math.pi),
            # Round-robin, not a random pick: with three colours and a dozen blobs a random
            # assignment leaves one colour missing about one run in fifty.
            color=i))
    keep = max(1, int(round(BLOB_COUNT * _clamp(density, DENSITY_RANGE))))
    return blobs[:keep]


def _clamp(value: float, bounds) -> float:
    lo, hi = bounds
    return float(min(max(value, lo), hi))


def geometry(blobs, t: float, width: int, height: int, speed: float = DEFAULT_SPEED,
             size: float = DEFAULT_SIZE) -> list:
    """(cx, cy, radius, colour index) per blob, in pixels, at time `t`."""
    short = max(min(width, height), 1)
    out = []
    tt = t * _clamp(speed, SPEED_RANGE)
    for b in blobs:
        cx = (b.x + b.drift_x * math.sin(b.rate_x * tt + b.phase_x)) * width
        cy = (b.y + b.drift_y * math.sin(b.rate_y * tt + b.phase_y)) * height
        r = b.radius * (1.0 + b.pulse * math.sin(b.pulse_rate * tt + b.pulse_phase))
        out.append((cx, cy, max(r * short * _clamp(size, SIZE_RANGE), 1.0), b.color))
    return out


class AmbientWidget(QtWidgets.QWidget):
    """The blob field, painted behind everything else.

    Transparent to the mouse -- it is scenery, and a background that swallowed a click on the panel
    in front of it would be a bug nobody could describe. The timer runs only while the widget is
    visible, because sixty repaints a second of a picture nobody is looking at is sixty repaints a
    second of a picture nobody is looking at.
    """

    def __init__(self, colors=None, background="#0b0d10", parent=None, seed: int = 0):
        super().__init__(parent)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.background = QtGui.QColor(background)
        self.colors = [QtGui.QColor(c) for c in (colors or ["#43c6d8", "#7a5cd8", "#2a9fb0"])]
        # `blob_size`, NOT `size`: QWidget.size() is a method, and an attribute of that name
        # shadows it -- every caller asking this widget how big it is gets a float back and
        # "'float' object is not callable" from somewhere that never mentioned size.
        self.speed, self.blob_size, self.density = DEFAULT_SPEED, DEFAULT_SIZE, DEFAULT_DENSITY
        self.blur = DEFAULT_BLUR
        self.seed = seed
        self.blobs = field(seed, self.density)
        self._t = 0.0
        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(40)                    # 25 fps: scenery, not an animation to watch
        self._timer.timeout.connect(self._tick)

    def configure(self, **kw) -> None:
        """Set any of speed, size, density, blur, colors -- and rebuild if the count changed.

        `size` is accepted under spaCR's name and stored as `blob_size`, since a widget cannot have
        an attribute called `size`.
        """
        density = kw.pop("density", None)
        if "size" in kw:
            kw["blob_size"] = kw.pop("size")
        for k, v in kw.items():
            if k == "colors":
                self.colors = [QtGui.QColor(c) for c in v]
            elif hasattr(self, k):
                setattr(self, k, float(v))
        if density is not None and float(density) != self.density:
            self.density = float(density)
            self.blobs = field(self.seed, self.density)
        self.update()

    def start(self) -> None:
        """Begin animating, if it is not already."""
        if not self._timer.isActive():
            self._timer.start()

    def stop(self) -> None:
        """Stop animating. The field keeps its time, so starting again continues rather than jumps."""
        self._timer.stop()

    def _tick(self) -> None:
        self._t += self._timer.interval() / 1000.0
        self.update()

    def set_time(self, t: float) -> None:
        """Jump to an instant. Position is a function of the clock, so this is exact."""
        self._t = float(t)
        self.update()

    def showEvent(self, ev):
        self.start()
        super().showEvent(ev)

    def hideEvent(self, ev):
        self.stop()
        super().hideEvent(ev)

    def paintEvent(self, _ev):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(), self.background)
        p.setPen(QtCore.Qt.PenStyle.NoPen)
        p.setCompositionMode(QtGui.QPainter.CompositionMode.CompositionMode_Plus)
        for cx, cy, r, ci in geometry(self.blobs, self._t, self.width(), self.height(),
                                      self.speed, self.blob_size):
            base = self.colors[ci % len(self.colors)]
            grad = QtGui.QRadialGradient(cx, cy, r)
            # Alpha to zero at the rim rather than a hard edge: a blob with an edge is a circle, and
            # a circle drifting behind a panel reads as a loading spinner.
            inner = QtGui.QColor(base)
            inner.setAlpha(int(70 * max(0.35, 1.0 - self.blur / (BLUR_RANGE[1] + 1e-9))))
            outer = QtGui.QColor(base)
            outer.setAlpha(0)
            grad.setColorAt(0.0, inner)
            grad.setColorAt(1.0, outer)
            p.setBrush(QtGui.QBrush(grad))
            p.drawEllipse(QtCore.QRectF(cx - r, cy - r, 2 * r, 2 * r))
        p.end()
