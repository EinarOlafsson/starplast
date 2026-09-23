# Packaging

One launcher, three installers, following spaCR's arrangement so a runtime fix lands once.

    bash packaging/build_debian.sh      # .deb        (Linux)
    bash packaging/build_macos.sh       # .app + .dmg (macOS)
    powershell -File packaging\build_windows.ps1   # .exe (Windows)

All three run `pyinstaller packaging/starplast.spec`, which bundles `starplast/data/`, the built cache and application assets.
That is what makes an installed starplast work with no network, no dataset and no build step.

**Gotchas, learned the hard way in spaCR:**

* `*.spec` is easy to gitignore by accident; this one is committed on purpose.
* Stage the build off the root disk if it is small: `STAGE_DIR=/mnt/big/deb bash packaging/build_debian.sh`.
* The Debian package depends on `libgl1`, `libegl1` and `libxcb-cursor0`. Without the last one Qt 6
  fails with "could not load the xcb platform plugin" and no further explanation.
* `PYQTGRAPH_QT_LIB=PyQt6` is set by the launcher. A bundle can contain more than one Qt binding, and
  pyqtgraph binding to the wrong one fails at import rather than at use.

Python package publishing is documented in [the release guide](../docs/releases.md).
