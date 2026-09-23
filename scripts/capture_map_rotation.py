#!/usr/bin/env python3
"""Capture a seamless rotation of the bundled map using Starplast's OpenGL view.

Requires the application dependencies, an X server, and ffmpeg. For example::

    xvfb-run -a -s '-screen 0 1280x960x24' python scripts/capture_map_rotation.py

Preferences and cache writes are isolated in a temporary directory. Coordinates
and compartment colours come from the app; no embedding is recomputed.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main() -> None:
    """Render one orbit, then encode a looping GIF with a shared colour palette."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/screenshots/map_rotation.gif")
    parser.add_argument("--frames", type=int, default=360)
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--width", type=int, default=1440)
    parser.add_argument("--height", type=int, default=1080)
    parser.add_argument("--gene", default="TGME49_301440", help="Selected gene (CDPK1 by default)")
    parser.add_argument("--edge-layer", default="comention")
    args = parser.parse_args()
    if args.frames < 4 or not 1 <= args.fps <= 50:
        parser.error("use at least four frames and 1–50 frames per second")
    if args.width < 320 or args.height < 240:
        parser.error("use a viewport of at least 320×240")
    from PyQt6 import QtCore, QtGui, QtWidgets
    import numpy as np

    with tempfile.TemporaryDirectory(prefix="starplast-map-") as temporary:
        state = Path(temporary)
        os.environ["STARPLAST_STATE"] = str(state / "state")
        os.environ["XDG_CACHE_HOME"] = str(state / "cache")
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
        QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat,
                                QtCore.QSettings.Scope.UserScope, str(state))
        from starplast.app import Window

        app = QtWidgets.QApplication([])
        window = Window(species="Toxoplasma gondii")
        window.apply_theme("dark")
        window.set_lighting("deep ray traced")
        window.set_lighting_option("source", "selected gene and its edges")
        window.set_lighting_option("point_mode", "glossy 3D")
        window.set_lighting_option("target_marker", "beacon")
        window.set_lighting_option("mood", "neutral")
        window.point_size = 5.0 * args.width / 1440
        window.depth_cue = False
        window.show_ground = False
        for timer in window.findChildren(QtCore.QTimer):
            timer.stop()
        if args.edge_layer not in window.edges:
            parser.error(f"unknown edge layer: {args.edge_layer}")
        window.edge_on = {key: key == args.edge_layer for key in window.edge_on}
        selected = np.flatnonzero(window.nodes.gene_id.eq(args.gene).to_numpy())
        if not len(selected):
            parser.error(f"gene is not in the bundled map: {args.gene}")
        window.on_pick(int(selected[0]))
        window.view.setParent(None)
        window.view.resize(args.width, args.height)
        center = (window.xyz.max(axis=0) + window.xyz.min(axis=0)) / 2
        radius = float(np.linalg.norm(window.xyz - center, axis=1).max())
        window.view.opts["center"] = QtGui.QVector3D(*map(float, center))
        window.view.setCameraPosition(distance=radius * 2.4, elevation=18, azimuth=45)
        window.view.show()
        for _ in range(5):
            app.processEvents()
        for index in range(args.frames):
            window.view.setCameraPosition(azimuth=45 + index * 360 / args.frames)
            # Refresh camera-dependent shading and depth cues at each exact orbit angle.
            window.redraw()
            for edge in window.edge_items:
                edge.setData(width=args.width / 800)
            app.processEvents()
            frame = window.view.grabFramebuffer()
            path = state / f"frame-{index:04}.png"
            if frame.isNull() or not frame.save(str(path)):
                raise RuntimeError(f"Failed to capture {path}")
            if index % 30 == 0:
                print(f"Captured {index + 1}/{args.frames}", flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        (args.output.with_suffix(".png")).write_bytes((state / "frame-0000.png").read_bytes())
        neighbours = window.edge_neighbours(window.sel)
        metadata = {
            "species": window.species, "gene": args.gene,
            "product": str(window.nodes.iloc[window.sel]["product"]),
            "edge_layer": args.edge_layer, "attention_corrected": window.attn_on,
            "neighbours": window.nodes.iloc[neighbours].gene_id.tolist(),
            "lighting": window._lighting, "width": args.width, "height": args.height,
            "point_size": window.point_size, "depth_cue": window.depth_cue,
            "show_ground": window.show_ground,
            "frames": args.frames, "fps": args.fps, "seconds": args.frames / args.fps,
        }
        args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
        subprocess.run([
            "ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
            "-framerate", str(args.fps), "-i", str(state / "frame-%04d.png"),
            "-filter_complex", "[0:v]split[a][b];[a]palettegen=stats_mode=full[p];[b][p]paletteuse=dither=none",
            "-loop", "0", str(args.output),
        ], check=True)
        window.view.close()
        window.close()
        app.processEvents()
    print(f"Saved {args.output} ({args.output.stat().st_size / 1024**2:.2f} MiB)")


if __name__ == "__main__":
    main()
