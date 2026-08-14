#!/usr/bin/env python3
"""The ball one gene is drawn as.

The claim under test is the one that was wrong for two rounds: that a gene looks like a sphere. It
cannot be tested by asking whether the shading code ran -- it ran the whole time -- so these tests
ask about the PICTURE. Is there a gradient across one ball, how much of it is highlight, and does
the highlight move when the light does.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import sprite as SP  # noqa: E402


def body(tex):
    """The lit pixels of one ball -- inside the disc, ignoring the antialiased edge."""
    return tex[:, :, 0][tex[:, :, 3] > 200].astype(float)


# --------------------------------------------------------------------------- the texture
def test_a_ball_is_a_gradient_and_not_a_disc_of_one_colour():
    """The whole bug, in one assertion. pyqtgraph's own sprite is `pData[:] = 255` with an alpha
    disc cut out of it: every pixel inside a ball carries the identical colour, so a ball is a flat
    disc and no amount of per-point lighting can make it look round. Reported twice."""
    flat = body(SP.texture("2D"))
    assert flat.max() - flat.min() == 0, "the 2D sprite is not flat"
    for name in ("matt", "satin", "glossy", "metallic"):
        lit = body(SP.texture(name))
        assert lit.max() - lit.min() > 120, f"{name} spans only {lit.max() - lit.min():.0f}/255"


def test_the_finishes_differ_in_how_much_of_the_ball_is_highlight():
    """Chalk is bright nearly everywhere, a metal is dark with a glint. That ordering IS the
    difference between the finishes -- if two of them light the same share of the ball, one of them
    is a duplicate menu entry."""
    share = {name: float((body(SP.texture(name)) > 200).mean())
             for name in ("matt", "satin", "glossy", "metallic")}
    assert share["matt"] > share["satin"] > share["glossy"] > share["metallic"], share
    assert share["metallic"] < 0.25, f"the metal is lit over {share['metallic']:.0%} of its face"


def test_the_highlight_follows_the_light():
    """One texture serves the whole cloud, so the direction baked into it is the scene's own light.
    Rebuilt as the lights move -- otherwise every ball is lit from the top left while the map is lit
    from the right, and the two read as different scenes."""
    left = SP.texture("glossy", (-0.7, 0.2, 0.7))
    right = SP.texture("glossy", (0.7, 0.2, 0.7))
    half = left.shape[1] // 2
    assert left[:, :half, 0].mean() > left[:, half:, 0].mean(), "the light is not on the left"
    assert right[:, half:, 0].mean() > right[:, :half, 0].mean(), "the light is not on the right"


def test_a_light_of_no_direction_falls_back_rather_than_dividing_by_zero():
    assert np.array_equal(SP.texture("glossy", (0.0, 0.0, 0.0)), SP.texture("glossy"))


def test_an_unknown_finish_draws_the_plain_disc():
    """A name from a settings file written by a newer version. A flat ball is the honest answer;
    raising here would take the whole map down over a preference."""
    assert np.array_equal(SP.texture("holographic"), SP.texture("2D"))


def test_the_edge_matches_the_one_pyqtgraph_draws():
    """Halos and centroids are drawn with the stock sprite beside these. A different edge on the
    genes reads as a different kind of object rather than the same object lit."""
    alpha = SP.disc_alpha(64)
    assert alpha[32, 32] == pytest.approx(1.0)
    assert alpha[0, 0] == 0.0
    assert 0.0 < alpha[32, 0] < 1.0, "the edge is not antialiased"


def test_a_ball_is_lit_inside_its_disc_and_flat_outside_it():
    """Outside the circle the texture is invisible anyway, but leaving it at the ambient value keeps
    the mipmap from smearing a bright edge into the alpha as the map is zoomed out."""
    tex = SP.texture("glossy")
    corner = tex[0, 0, 0]
    assert tex[:, :, 3][0, 0] == 0
    assert corner == int(SP.SPRITES["glossy"]["ambient"] * 255)


# --------------------------------------------------------------------------- world to screen
def test_a_world_direction_becomes_a_screen_direction():
    basis = (np.array([0.0, 0.0, 90.0]), np.array([1.0, 0.0, 0.0]),
             np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0]))
    assert SP.to_screen((5.0, 0.0, 0.0), basis) == pytest.approx((1.0, 0.0, 0.0))
    assert SP.to_screen((0.0, 0.0, 4.0), basis) == pytest.approx((0.0, 0.0, 1.0))


def test_a_direction_of_nothing_keeps_the_default():
    basis = (np.zeros(3), np.array([1.0, 0.0, 0.0]),
             np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0]))
    assert SP.to_screen((0.0, 0.0, 0.0), basis) == SP.DEFAULT_LIGHT


# --------------------------------------------------------------------------- getting it to the card
def test_the_sprite_waits_for_a_context_rather_than_uploading_from_wherever(qapp, monkeypatch):
    """A texture can only be uploaded while a GL context is current, and neither "the finish
    changed" nor "a light moved" has one. Holding it until the next paint is what keeps the choice
    of sprite a question about the model rather than about the paint clock."""
    sent = []
    monkeypatch.setattr(SP, "upload", lambda tex, data: sent.append((tex, data.shape)))
    monkeypatch.setattr(SP.gl.GLScatterPlotItem, "paint", lambda self: None)
    item = SP.ShadedScatter(pos=np.zeros((3, 3)), size=5.0)
    item.set_sprite(SP.texture("glossy"))
    assert sent == [], "uploaded with no texture to upload into"
    item.paint()
    assert sent == [], "uploaded before the item had a texture name"
    item.pointTexture = 7
    item.paint()
    assert sent == [(7, (SP.WIDTH, SP.WIDTH, 4))]
    item.paint()
    assert len(sent) == 1, "re-uploaded a texture that had not changed"


def test_the_pending_sprite_survives_until_the_item_is_initialised(qapp, monkeypatch):
    sent = []
    monkeypatch.setattr(SP, "upload", lambda tex, data: sent.append(tex))
    monkeypatch.setattr(SP.gl.GLScatterPlotItem, "initializeGL", lambda self: None)
    item = SP.ShadedScatter(pos=np.zeros((2, 3)), size=5.0)
    item.set_sprite(SP.texture("matt"))
    item.pointTexture = 3
    item.initializeGL()
    assert sent == [3] and item._pending is None


def test_upload_hands_the_array_to_the_card_unchanged(qapp, monkeypatch):
    """The one place that talks to OpenGL. Driven with the calls recorded rather than made, because
    a test suite that needs a graphics card is a test suite that does not run."""
    calls = {}

    class FakeGL:
        GL_TEXTURE_2D, GL_RGBA, GL_UNSIGNED_BYTE = 1, 2, 3

        @staticmethod
        def glBindTexture(target, name):
            calls["bound"] = (target, name)

        @staticmethod
        def glTexImage2D(*args):
            calls["image"] = args

    monkeypatch.setattr(SP, "GL", FakeGL)
    data = SP.texture("satin")
    SP.upload(11, data)
    assert calls["bound"] == (1, 11)
    assert calls["image"][3:5] == (SP.WIDTH, SP.WIDTH)
    assert calls["image"][-1] is data
