# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec shared by the Windows and macOS builds.

    pyinstaller packaging/starplast.spec

The committed cache is the whole point of bundling: data/ is 11 MB and makes the installed application
work with no network and no dataset, so it is collected rather than fetched on first run.
"""
from __future__ import annotations

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

ROOT = Path(SPECPATH).resolve().parent
ENTRY = str(ROOT / "packaging" / "starplast_launcher.py")

# PyInstaller's static analysis misses these: pyqtgraph and umap import lazily, and sklearn's
# HDBSCAN and the numba-compiled umap internals are only reachable at run time.
hiddenimports = []
for pkg in ("starplast", "pyqtgraph", "pyqtgraph.opengl", "OpenGL", "OpenGL.platform",
            "sklearn.cluster", "sklearn.manifold", "sklearn.metrics",
            "umap", "numba", "pynndescent", "scipy.stats", "pyarrow"):
    try:
        hiddenimports += collect_submodules(pkg)
    except Exception:
        hiddenimports.append(pkg)

datas = [(str(ROOT / "data"), "data")]          # the committed cache: this is what makes it standalone
for pkg in ("pyqtgraph", "umap", "pynndescent"):
    try:
        datas += collect_data_files(pkg)
    except Exception:
        pass

a = Analysis([ENTRY], pathex=[str(ROOT)], binaries=[], datas=datas,
             hiddenimports=hiddenimports, hookspath=[], runtime_hooks=[],
             excludes=["tkinter", "tensorflow", "torch", "PySide6", "PyQt5"],
             noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="starplast",
          console=False, icon=None)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="starplast")

app = BUNDLE(coll, name="starplast.app", icon=None, bundle_identifier="io.github.einarolafsson.starplast",
             info_plist={"NSHighResolutionCapable": True, "LSMinimumSystemVersion": "11.0"})
