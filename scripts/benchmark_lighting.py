#!/usr/bin/env python3
"""Capture fixed lighting comparisons and measure full-map interactive frame time.

Run under an X server so the real OpenGL scatter path paints rather than a plotting substitute::

    xvfb-run -a -s '-screen 0 1920x1080x24' python scripts/benchmark_lighting.py

The report records the GL renderer because Xvfb commonly selects software rendering. That result is
still a conservative responsiveness check, but it must not be presented as an NVIDIA measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import statistics
import sys
import time

import numpy as np
from PIL import Image, ImageDraw
from PyQt6 import QtCore, QtWidgets
from OpenGL import GL

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from starplast.app import Window
from starplast.lighting import LIGHT_MOODS, POINT_MODES


def arguments() -> argparse.Namespace:
    """Command-line output path and the number of timed frames."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("results/lighting_2026_08_14"))
    parser.add_argument("--frames", type=int, default=20)
    return parser.parse_args()


def render_name(mode: str, point_mode: str, mood: str) -> str:
    """A filesystem-safe stable name for one comparison frame."""
    return "__".join((mode, point_mode, mood)).replace(" ", "_") + ".png"


def capture(app, window: Window, output: Path, mode: str, point_mode: str,
            mood: str) -> Path:
    """Paint one fixed scene through Starplast's real GL view and save its framebuffer."""
    window.set_lighting(mode)
    window._light_timer.stop()  # fixed comparison, not a frame chosen by timer scheduling
    window.set_lighting_option("source", "selected gene")
    window.set_lighting_option("point_mode", point_mode)
    window.set_lighting_option("mood", mood)
    window._light_tick()
    for _ in range(3):
        app.processEvents()
    path = output / render_name(mode, point_mode, mood)
    image = window.view.grabFramebuffer()
    if image.isNull() or not image.save(str(path)):
        raise RuntimeError(f"OpenGL framebuffer could not be saved to {path}")
    return path


def contact_sheet(paths: list[Path], labels: list[str], output: Path) -> None:
    """Join comparison frames with labels while leaving every frame at the same scale."""
    frames = [Image.open(path).convert("RGB") for path in paths]
    thumb = (500, 400)
    header = 34
    sheet = Image.new("RGB", (thumb[0] * len(frames), thumb[1] + header), "#080c10")
    draw = ImageDraw.Draw(sheet)
    for i, (frame, label) in enumerate(zip(frames, labels, strict=True)):
        frame.thumbnail(thumb, Image.Resampling.LANCZOS)
        x = i * thumb[0]
        sheet.paste(frame, (x, header))
        draw.text((x + 12, 10), label, fill="#eef4f8")
    sheet.save(output)


def renderer_name(window: Window) -> str:
    """Return the active OpenGL renderer, or an explicit unknown marker."""
    try:
        window.view.makeCurrent()
        value = GL.glGetString(GL.GL_RENDERER)
        return value.decode("utf-8", "replace") if value else "unknown"
    except Exception as exc:  # renderer reporting must not discard otherwise valid evidence
        return f"unknown ({type(exc).__name__})"


def benchmark(app, window: Window, frames: int, mode: str, point_mode: str,
              mood: str) -> dict:
    """Time update plus actual GL paint for the complete 8,140-point scene."""
    window.set_lighting(mode)
    window._light_timer.stop()
    window.set_lighting_option("point_mode", point_mode)
    window.set_lighting_option("mood", mood)
    for _ in range(3):
        window._light_tick()
        window.view.repaint()
        app.processEvents()
    elapsed = []
    for _ in range(max(frames, 1)):
        start = time.perf_counter()
        window._light_tick()
        window.view.repaint()
        app.processEvents()
        elapsed.append((time.perf_counter() - start) * 1000.0)
    ordered = sorted(elapsed)
    p95 = ordered[min(len(ordered) - 1, int(np.ceil(0.95 * len(ordered))) - 1)]
    return {"light_mode": mode, "point_mode": point_mode, "mood": mood,
            "genes": len(window.xyz), "frames": len(elapsed),
            "median_ms": round(statistics.median(elapsed), 3), "p95_ms": round(p95, 3)}


def image_difference(first: Path, second: Path) -> dict:
    """Visible pixel difference between two fixed-camera framebuffer captures."""
    a = np.asarray(Image.open(first).convert("RGB"), dtype=np.int16)
    b = np.asarray(Image.open(second).convert("RGB"), dtype=np.int16)
    delta = np.abs(a - b)
    return {"first": first.name, "second": second.name,
            "mean_abs_rgb": round(float(delta.mean()), 3),
            "pixel_share_over_5": round(float((delta.max(axis=2) > 5).mean()), 4),
            "max_channel_delta": int(delta.max())}


def static_frame_hashes(app, window: Window, frames: int = 12) -> list[str]:
    """Hashes of repeated static ray-traced frames; more than one means temporal instability."""
    window.set_lighting("ray traced")
    window._light_timer.stop()
    window.set_lighting_option("source", "selected gene")
    window.set_lighting_option("point_mode", "metallic 3D")
    hashes = []
    for _ in range(frames):
        window._light_tick()
        window.view.repaint()
        app.processEvents()
        image = window.view.grabFramebuffer()
        hashes.append(hashlib.sha256(image.bits().asstring(image.sizeInBytes())).hexdigest())
    return hashes


def main() -> int:
    """Create screenshots, contact sheets, and a renderer-labelled timing table."""
    args = arguments()
    args.output.mkdir(parents=True, exist_ok=True)
    state = args.output / "isolated_state"
    state.mkdir(exist_ok=True)
    QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
    QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat,
                             QtCore.QSettings.Scope.UserScope, str(state))
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = Window()
    # Detaching avoids Starplast's analysis panels squeezing the GL view on a headless X screen.
    window.view.setParent(None)
    window.view.resize(1000, 800)
    window.view.show()
    for _ in range(5):
        app.processEvents()
    window.sel = int(np.argmin(np.linalg.norm(window.xyz - window.xyz.mean(axis=0), axis=1)))

    point_paths = [capture(app, window, args.output, "soft", mode, "neutral")
                   for mode in POINT_MODES]
    contact_sheet(point_paths, list(POINT_MODES), args.output / "point_modes.png")
    mood_paths = [capture(app, window, args.output, "soft", "glossy 3D", mood)
                  for mood in LIGHT_MOODS]
    contact_sheet(mood_paths, list(LIGHT_MOODS), args.output / "light_moods.png")
    transport_paths = [capture(app, window, args.output, mode, "glossy 3D", "neutral")
                       for mode in ("soft", "ray traced")]
    contact_sheet(transport_paths, ["soft", "ray traced"],
                  args.output / "light_transport.png")

    comparisons = [
        image_difference(point_paths[0], point_paths[1]),
        image_difference(point_paths[1], point_paths[2]),
        image_difference(transport_paths[0], transport_paths[1]),
        image_difference(mood_paths[0], mood_paths[1]),
        image_difference(mood_paths[0], mood_paths[2]),
    ]
    with (args.output / "visual_differences.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparisons[0].keys())
        writer.writeheader()
        writer.writerows(comparisons)

    hashes = static_frame_hashes(app, window)
    with (args.output / "ray_stability.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("frames", "unique_frame_hashes", "stable"))
        writer.writerow((len(hashes), len(set(hashes)), len(set(hashes)) == 1))

    renderer = renderer_name(window)
    rows = [benchmark(app, window, args.frames, light_mode, point_mode, "neutral")
            for light_mode in ("soft", "ray traced") for point_mode in POINT_MODES]
    report = args.output / "frame_times.csv"
    with report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("renderer", *rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({"renderer": renderer, **row})
    window.view.close()
    window.close()
    print(f"wrote {len(point_paths) + len(mood_paths) + len(transport_paths)} frames and {report}")
    print(f"static ray frames: {len(set(hashes))} unique of {len(hashes)}")
    print(f"OpenGL renderer: {renderer}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
