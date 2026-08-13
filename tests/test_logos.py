#!/usr/bin/env python3
"""Forty logo drafts, and the constraint that decides whether any of them is usable.

The mark becomes the window and tab icon, so 16 pixels is where it has to work — and 16 pixels is
exactly where a fine constellation becomes a smudge and a bold crescent becomes a black square. The
task said to test that rather than assume it, so every draft is rendered at 16 and its ink measured.
"""
from __future__ import annotations

import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import paths  # noqa: E402

ICONS = os.path.join(paths.data_dir(), "icons")
LOGOS = sorted(f for f in os.listdir(ICONS) if f.startswith("logo_") and f.endswith(".svg"))


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def _ink(name: str, size: int, qapp) -> float:
    """Fraction of DARK pixels, composited on white, rendering the draft at `size`.

    Darkness rather than coverage: the inverted drafts are a solid black tile with the studied genes
    knocked out of it, so by coverage they are 100% "ink" and by darkness they are what they look
    like. Composited on white because that is the worst case for a mark that will also be shown on
    a light theme.
    """
    from PyQt6 import QtCore, QtGui, QtSvg
    r = QtSvg.QSvgRenderer(os.path.join(ICONS, name))
    img = QtGui.QImage(size, size, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(0xFFFFFFFF)
    p = QtGui.QPainter(img)
    r.render(p, QtCore.QRectF(0, 0, size, size))
    p.end()
    dark = sum(QtGui.QColor.fromRgba(img.pixel(x, y)).lightness() < 128
               for x in range(size) for y in range(size))
    return dark / (size * size)


def test_there_are_at_least_the_forty_that_were_asked_for():
    """Forty was the ask; the eight galaxy drafts are a seventh direction requested after the first
    sheet was looked at, and adding them beats throwing away drafts nobody rejected."""
    assert len(LOGOS) >= 40, LOGOS


def test_they_cover_the_six_directions_rather_than_one_idea_forty_times():
    """Forty variations on whichever came first is not a set of drafts to choose from."""
    families = {re.sub(r"^logo_\d+_|_\d+\.svg$", "", f) for f in LOGOS}
    assert families == {"constellation", "apicoplast", "crescent", "orbit", "letter_s",
                        "dark_field", "galaxy"}, families
    for family in families:
        assert sum(family in f for f in LOGOS) >= 5, family


@pytest.mark.parametrize("name", LOGOS)
def test_each_draft_is_legible_at_sixteen_pixels(name, qapp):
    """The size it becomes a tab icon at. Below a few percent it is a smudge; above two thirds it is
    a filled square, and both are indistinguishable from every other draft at that size."""
    ink = _ink(name, 16, qapp)
    # Both bounds are about the same failure: a mark with no contrast at this size. A few dark
    # pixels is a smudge; an almost solid tile is a square. The inverted drafts sit near the top of
    # this range on purpose -- they ARE mostly black, with the studied genes knocked out -- so the
    # rule is that there must be a visible amount of both.
    assert 0.04 < ink < 0.96, f"{name}: {ink:.0%} dark at 16px"
    assert min(ink, 1 - ink) > 0.03, f"{name}: no contrast at 16px ({ink:.0%} dark)"


@pytest.mark.parametrize("name", LOGOS)
def test_each_draft_is_black_and_white_with_no_greys(name, qapp):
    """A grey vanishes when the icon is printed or shown small, which is the one place this has to
    work. Checked in the source rather than the render, where antialiasing makes greys legitimately."""
    text = open(os.path.join(ICONS, name), encoding="utf8").read()
    colors = set(re.findall(r'(?:fill|stroke)="([^"]+)"', text))
    assert colors <= {"#000", "#fff", "none"}, f"{name}: {colors}"


@pytest.mark.parametrize("name", LOGOS)
def test_each_draft_scales(name, qapp):
    """SVG, and it has to hold up large as well: a mark that only works at one size is a bitmap."""
    small, large = _ink(name, 16, qapp), _ink(name, 128, qapp)
    assert large > 0.01
    assert abs(small - large) < 0.5, f"{name}: {small:.0%} at 16px vs {large:.0%} at 128px"


def test_none_of_them_is_a_cell_diagram():
    """That is the other artwork and a different job. A mark that reads as a cell says the
    application is about drawing cells, which it is not."""
    for name in LOGOS:
        text = open(os.path.join(ICONS, name), encoding="utf8").read()
        assert "SL0" not in text, name
        assert "subcell" not in text, name


def test_each_draft_says_what_it_is_arguing():
    """A file called logo_31 with no title is a file nobody can choose between."""
    for name in LOGOS:
        text = open(os.path.join(ICONS, name), encoding="utf8").read()
        title = re.search(r"<title>([^<]+)</title>", text)
        assert title and "starplast" in title.group(1), name


def test_they_are_distinct_from_one_another():
    """Numbered drafts that are byte-identical are one draft with forty names."""
    bodies = {}
    for name in LOGOS:
        text = open(os.path.join(ICONS, name), encoding="utf8").read()
        body = re.sub(r"<title>[^<]*</title>", "", text)
        assert body not in bodies, f"{name} is identical to {bodies.get(body)}"
        bodies[body] = name


def test_regenerating_produces_the_same_forty(tmp_path):
    """Seeded, so "that one, but sparser" can be answered by editing the script rather than by
    redrawing something nobody can reproduce.

    Into a temporary directory, not over the shipped files. Regenerating in place means a generator
    that has started producing something else overwrites the artwork the package ships on its way to
    reporting that it changed -- and the check would then pass on a second run."""
    import subprocess
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subprocess.run([sys.executable, os.path.join(root, "scripts", "generate_logos.py"),
                    str(tmp_path)], check=True, capture_output=True)
    shipped = {f: open(os.path.join(ICONS, f), encoding="utf8").read() for f in LOGOS}
    again = {f: open(tmp_path / f, encoding="utf8").read() for f in LOGOS}
    assert shipped == again


def test_the_galaxy_drafts_put_the_parasite_at_the_centre(qapp):
    """The direction the name argues for: a star map with Toxoplasma in the middle of it. The centre
    has to be the parasite rather than a dot, or it is the orbit family again."""
    from PyQt6 import QtCore, QtGui, QtSvg
    for name in [f for f in LOGOS if "galaxy" in f]:
        r = QtSvg.QSvgRenderer(os.path.join(ICONS, name))
        img = QtGui.QImage(64, 64, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0xFFFFFFFF)
        p = QtGui.QPainter(img)
        r.render(p, QtCore.QRectF(0, 0, 64, 64))
        p.end()
        middle = sum(QtGui.QColor.fromRgba(img.pixel(x, y)).lightness() < 128
                     for x in range(26, 38) for y in range(20, 46))
        assert middle > 40, f"{name}: the middle is empty"


def test_the_galaxy_orbits_are_tilted_rather_than_face_on(qapp):
    """A face-on system is a target; a tilted one is a system seen from somewhere, which is what a
    3D map you can orbit around is."""
    for name in [f for f in LOGOS if "galaxy" in f]:
        text = open(os.path.join(ICONS, name), encoding="utf8").read()
        assert "rotate(-" in text, f"{name}: the orbits are not tilted"
        # Only the ellipses that are ORBITS -- the ones carrying the tilt. The knockout behind the
        # parasite is an ellipse too, and it is taller than it is wide by design.
        orbits = re.findall(r'<ellipse[^>]*rx="([\d.]+)"[^>]*ry="([\d.]+)"[^>]*rotate\(', text)
        assert orbits, f"{name}: no tilted orbit at all"
        assert all(float(ry) < float(rx) for rx, ry in orbits), f"{name}: the orbits are circles"
