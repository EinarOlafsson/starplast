#!/usr/bin/env python3
"""pyqtgraph API compatibility, pinned so a version bump cannot break picking again.

`GLViewWidget.projectionMatrix()` has changed signature twice across pyqtgraph releases:

    older     projectionMatrix()
    0.13.x    projectionMatrix(region=None)
    current   projectionMatrix(region, viewport)      -- both required

Calling the old way on a current release raises TypeError *inside* `mouseReleaseEvent`, so every click
printed a traceback and selected nothing. The development machine runs 0.13.7 and never takes the failing
path, which is exactly why this is tested against stubs rather than against whatever happens to be
installed.

The convention is chosen by inspecting the signature rather than by catching TypeError, because an
exception chain cannot distinguish "called wrongly" from "something inside raised TypeError" -- it would
retry a call that failed for an unrelated reason and then report the wrong cause.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


class _Matrix:
    """Stands in for a QMatrix4x4: enough of it for _mvp to multiply and read rows."""

    class _Row:
        def __init__(self, v):
            self._v = v

        def x(self):
            return self._v[0]

        def y(self):
            return self._v[1]

        def z(self):
            return self._v[2]

        def w(self):
            return self._v[3]

    def __init__(self, tag="I"):
        self.tag = tag

    def __mul__(self, other):
        return self

    def row(self, i):
        return _Matrix._Row([1.0 if j == i else 0.0 for j in range(4)])


def _view_with(signature: str):
    """A minimal stand-in for Map3D carrying one of the three signatures."""
    from starplast.app import Map3D
    import numpy as np

    view = Map3D.__new__(Map3D)                 # no Qt widget: only the compat logic is under test
    view.xyz = np.zeros((3, 3), dtype=np.float32)
    view._proj_kind = None
    view.calls = []

    if signature == "none":
        def projectionMatrix():
            view.calls.append(())
            return _Matrix()
    elif signature == "region":
        def projectionMatrix(region=None):
            view.calls.append((region,))
            return _Matrix()
    elif signature == "region_viewport":
        def projectionMatrix(region, viewport):
            view.calls.append((region, viewport))
            return _Matrix()
    else:                                        # a C extension exposing no signature
        projectionMatrix = print                 # builtins refuse inspect.signature

    view.projectionMatrix = projectionMatrix
    view.viewMatrix = lambda: _Matrix()
    view.width = lambda: 800
    view.height = lambda: 600
    view.devicePixelRatioF = lambda: 2.0
    return view


@pytest.mark.parametrize("signature", ["none", "region", "region_viewport"])
def test_every_known_signature_is_called_correctly(signature):
    view = _view_with(signature)
    m = view._projection_matrix()
    assert m is not None
    assert len(view.calls) == 1
    if signature == "none":
        assert view.calls[0] == ()
    elif signature == "region":
        assert view.calls[0] == (None,)
    else:
        # viewport is in device pixels, so the device pixel ratio has to be applied
        assert view.calls[0] == (None, (0, 0, 1600, 1200))


def test_the_convention_is_resolved_once_not_per_frame():
    """project() runs on every click and every spin frame; re-inspecting each time is waste."""
    view = _view_with("region_viewport")
    view._projection_matrix()
    kind = view._proj_kind
    view._projection_matrix()
    assert view._proj_kind == kind == "region_viewport"
    assert len(view.calls) == 2


def test_uninspectable_callable_falls_back_rather_than_raising():
    """A C extension may expose no signature at all; that must degrade, not crash."""
    view = _view_with("builtin")
    view._proj_kind = None
    assert view._proj_kind is None
    try:
        view._projection_matrix()
    except TypeError:
        pass                                     # print() rejects the call, which is fine
    assert view._proj_kind == "none", "an uninspectable callable should resolve to the no-arg form"


def test_mvp_builds_a_4x4_from_the_matrix():
    import numpy as np
    view = _view_with("region")
    mvp = view._mvp()
    assert mvp.shape == (4, 4)
    assert np.isfinite(mvp).all()
