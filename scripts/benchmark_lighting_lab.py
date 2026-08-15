#!/usr/bin/env python3
"""Exercise the task-37 lighting lab through the production OpenGL renderer.

One-factor frames are retained for visual comparison.  Every meaningful cross-product is painted
and logged without saving thousands of redundant screenshots, catching shader/material combinations
that compile individually but fail together with a beam, mood, marker, source, or ray mode.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from itertools import product
from pathlib import Path
import statistics
import sys
import time

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
from OpenGL import GL

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from starplast.app import Window
from starplast.lighting import (LIGHT_MOODS, MODES, POINTER_MODES, POINT_MODES, RESPONSES,
                                SOURCES, TARGET_MARKERS)


def arguments() -> argparse.Namespace:
    """Output path and whether to skip the exhaustive combination smoke."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/lighting_lab_2026_08_15"))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def pixels(image: QtGui.QImage) -> np.ndarray:
    """Copy a framebuffer into a stable RGB array."""
    converted = image.convertToFormat(QtGui.QImage.Format.Format_RGBA8888)
    raw = converted.bits().asstring(converted.sizeInBytes())
    return np.frombuffer(raw, np.uint8).reshape(converted.height(), converted.width(), 4)[:, :, :3].copy()


def paint(app, window: Window, ticks: int = 1) -> tuple[QtGui.QImage, float]:
    """Update lighting, force an actual GL paint, and return framebuffer plus elapsed ms."""
    start = time.perf_counter()
    for _ in range(ticks):
        window._light_tick()
    window.view.repaint()
    app.processEvents()
    return window.view.grabFramebuffer(), (time.perf_counter() - start) * 1000.0


def configure(window: Window, *, transport="soft", source="mouse flashlight",
              pointer="broad flashlight", response="smooth", marker="none", mood="neutral",
              material="glossy 3D") -> None:
    """Apply one complete lab configuration through the same public setters as the UI."""
    window.set_lighting(transport)
    window._light_timer.stop()
    for key, value in (("source", source), ("pointer_mode", pointer), ("response", response),
                       ("target_marker", marker), ("mood", mood), ("point_mode", material)):
        window.set_lighting_option(key, value)


def slug(*parts) -> str:
    """Filesystem-safe comparison frame name."""
    return "__".join(str(part).replace(" ", "_") for part in parts) + ".png"


def one_factor(app, window: Window, output: Path) -> list[dict]:
    """Save every value of every independent comparison control against one neutral scene."""
    dimensions = {
        "transport": [m for m in MODES if m != "off"], "source": SOURCES,
        "pointer": POINTER_MODES, "response": RESPONSES, "marker": TARGET_MARKERS,
        "mood": LIGHT_MOODS, "material": POINT_MODES,
    }
    baseline = dict(transport="soft", source="mouse flashlight", pointer="broad flashlight",
                    response="smooth", marker="none", mood="neutral", material="glossy 3D")
    rows = []
    for dimension, values in dimensions.items():
        base_pixels = None
        for value in values:
            settings = {**baseline, dimension: value}
            window._smoothed_pointer = None
            configure(window, **settings)
            ticks = 4
            if dimension == "response":
                # Compare one frame after the same large cursor move. A settled static pointer is
                # deliberately identical for every response and would be a meaningless picture.
                window._smoothed_pointer = np.array([-0.8, 0.0, 1.0], dtype=float)
                window.view.pointer = (0.8, 0.0, 1.0)
                ticks = 1
            image, elapsed = paint(app, window, ticks=ticks)
            path = output / slug(dimension, value)
            image.save(str(path))
            frame = pixels(image)
            if base_pixels is None:
                base_pixels = frame
            delta = np.abs(frame.astype(np.int16) - base_pixels.astype(np.int16))
            rows.append({"dimension": dimension, "value": value, "frame": path.name,
                         "mean_abs_rgb_from_first": round(float(delta.mean()), 4),
                         "pixel_share_over_5": round(float((delta.max(axis=2) > 5).mean()), 6),
                         "paint_ms": round(elapsed, 3), "gpu_failed": window.scatter._gpu_failed})
    return rows


def combinations(app, window: Window, quick: bool) -> list[dict]:
    """Paint the meaningful cross-product and record failures/timing without image bloat."""
    transports = [m for m in MODES if m != "off"]
    mouse = product(transports, ["mouse flashlight"], POINTER_MODES, RESPONSES,
                    TARGET_MARKERS, LIGHT_MOODS, POINT_MODES)
    clicked = product(transports, SOURCES[1:], ["broad flashlight"], ["smooth"],
                      TARGET_MARKERS, LIGHT_MOODS, POINT_MODES)
    settings = list(mouse) + list(clicked)
    if quick:
        settings = settings[::max(len(settings) // 120, 1)]
    rows = []
    for index, (transport, source, pointer, response, marker, mood, material) in enumerate(settings, 1):
        window._smoothed_pointer = None
        configure(window, transport=transport, source=source, pointer=pointer, response=response,
                  marker=marker, mood=mood, material=material)
        image, elapsed = paint(app, window)
        digest = hashlib.sha256(image.bits().asstring(image.sizeInBytes())).hexdigest()
        rows.append({"transport": transport, "source": source, "pointer": pointer,
                     "response": response, "marker": marker, "mood": mood,
                     "material": material, "paint_ms": round(elapsed, 3),
                     "gpu_failed": window.scatter._gpu_failed, "frame_sha256": digest})
        if index % 500 == 0:
            print(f"painted {index}/{len(settings)} combinations", flush=True)
    return rows


def stability(app, window: Window) -> list[dict]:
    """Hash repeated static ray frames for both density-shadow depths."""
    rows = []
    for mode in ("ray traced", "deep ray traced"):
        configure(window, transport=mode, source="selected gene", marker="none",
                  material="metallic 3D")
        hashes, elapsed = [], []
        for _ in range(12):
            image, ms = paint(app, window)
            hashes.append(hashlib.sha256(image.bits().asstring(image.sizeInBytes())).hexdigest())
            elapsed.append(ms)
        rows.append({"transport": mode, "frames": len(hashes),
                     "unique_hashes": len(set(hashes)), "stable": len(set(hashes)) == 1,
                     "median_ms": round(statistics.median(elapsed), 3),
                     "p95_ms": round(float(np.percentile(elapsed, 95)), 3)})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a non-empty list of homogeneous dictionaries."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Build the real GL comparison corpus and exhaustive smoke report."""
    args = arguments()
    args.output.mkdir(parents=True, exist_ok=True)
    state = args.output / "isolated_state"
    state.mkdir(exist_ok=True)
    QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
    QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat,
                             QtCore.QSettings.Scope.UserScope, str(state))
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = Window()
    window.view.setParent(None)
    window.view.resize(1000, 800)
    window.view.show()
    for _ in range(6):
        app.processEvents()
    window.sel = int(np.argmin(np.linalg.norm(window.xyz - window.xyz.mean(0), axis=1)))
    window.view.pointer = (0.38, 0.22, 1.0)
    factors = one_factor(app, window, args.output)
    combos = combinations(app, window, args.quick)
    stable = stability(app, window)
    write_csv(args.output / "one_factor_differences.csv", factors)
    write_csv(args.output / "combinations.csv", combos)
    write_csv(args.output / "stability.csv", stable)
    window.view.makeCurrent()
    renderer = GL.glGetString(GL.GL_RENDERER)
    (args.output / "README.md").write_text(
        "# Lighting lab — 2026-08-15\n\n"
        f"Renderer: `{renderer.decode('utf-8', 'replace') if renderer else 'unknown'}`\n\n"
        f"Painted {len(combos):,} cross-product configurations and retained {len(factors)} "
        "one-factor frames. See the CSV files for timing, hashes, pixel differences, and shader "
        "failure state. The isolated settings directory prevents this benchmark from changing the "
        "user's choices.\n", encoding="utf-8")
    print(f"{len(factors)} frames, {len(combos)} combinations -> {args.output}")
    print(f"GPU failures: {sum(bool(row['gpu_failed']) for row in combos)}")
    print(f"stability: {stable}")
    window.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
