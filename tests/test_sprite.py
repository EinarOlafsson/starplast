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
    flat = body(SP.texture("flat"))
    assert flat.max() - flat.min() == 0, "the flat sprite is not flat"
    for name in ("glossy 3D", "metallic 3D"):
        lit = body(SP.texture(name))
        assert lit.max() - lit.min() > 120, f"{name} spans only {lit.max() - lit.min():.0f}/255"


def test_the_point_modes_differ_in_how_much_of_the_ball_is_highlight():
    """A flat marker, glossy bead, and metal bead must be recognizable without reading the menu."""
    share = {name: float((body(SP.texture(name)) > 200).mean())
             for name in ("flat", "glossy 3D", "metallic 3D")}
    assert share["flat"] > share["glossy 3D"] > share["metallic 3D"], share
    assert share["metallic 3D"] < 0.25, \
        f"the metal is lit over {share['metallic 3D']:.0%} of its face"


def test_the_highlight_follows_the_light():
    """One texture serves the whole cloud, so the direction baked into it is the scene's own light.
    Rebuilt as the lights move -- otherwise every ball is lit from the top left while the map is lit
    from the right, and the two read as different scenes."""
    left = SP.texture("glossy 3D", (-0.7, 0.2, 0.7))
    right = SP.texture("glossy 3D", (0.7, 0.2, 0.7))
    half = left.shape[1] // 2
    assert left[:, :half, 0].mean() > left[:, half:, 0].mean(), "the light is not on the left"
    assert right[:, half:, 0].mean() > right[:, :half, 0].mean(), "the light is not on the right"


def test_a_light_of_no_direction_falls_back_rather_than_dividing_by_zero():
    assert np.array_equal(SP.texture("glossy 3D", (0.0, 0.0, 0.0)),
                          SP.texture("glossy 3D"))


def test_an_unknown_finish_draws_the_plain_disc():
    """A name from a settings file written by a newer version. A flat ball is the honest answer;
    raising here would take the whole map down over a preference."""
    assert np.array_equal(SP.texture("holographic"), SP.texture("flat"))


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
    tex = SP.texture("glossy 3D")
    corner = tex[0, 0, 0]
    assert tex[:, :, 3][0, 0] == 0
    assert corner == int(SP.SPRITES["glossy 3D"]["ambient"] * 255)


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
    item.set_sprite(SP.texture("glossy 3D"))
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
    item.set_sprite(SP.texture("flat"))
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
    data = SP.texture("glossy 3D")
    SP.upload(11, data)
    assert calls["bound"] == (1, 11)
    assert calls["image"][3:5] == (SP.WIDTH, SP.WIDTH)
    assert calls["image"][-1] is data


# --------------------------------------------------------------------------- GPU PBR renderer
def test_the_gpu_shader_reconstructs_normals_depth_metal_and_density_rays():
    """Pin the capabilities that distinguish this from another pre-lit flat texture."""
    assert "gl_PointCoord" in SP._SPHERE_FRAGMENT
    assert "gl_FragDepth" in SP._SPHERE_FRAGMENT
    assert "distributionGGX" in SP._SPHERE_FRAGMENT
    assert "reflect(" in SP._SPHERE_FRAGMENT and "metallic" in SP._SPHERE_FRAGMENT
    assert "texture3D" in SP._SPHERE_VERTEX and f"RAY_STEPS = {SP.RAY_STEPS}" in SP._SPHERE_VERTEX


def test_pyqtgraph_014_vbo_shader_uses_attributes_and_core_glsl():
    """0.14 removed the registered pointSprite shader and all fixed-function inputs."""
    vertex, fragment = SP._vbo_shader_sources(core=True)
    assert vertex.startswith("\n#version 140")
    assert "in vec4 a_position" in vertex and "in vec4 a_color" in vertex
    assert "u_mvp * a_position" in vertex and "gl_Vertex" not in vertex
    assert "texture(uDensity" in vertex and "texture3D" not in vertex
    assert "out vec4 fragColor" in fragment and "gl_FragColor" not in fragment


def test_vbo_program_compiles_lazily_binds_all_attributes_and_is_cached(monkeypatch):
    calls = []

    class Format:
        @staticmethod
        def version():
            return (4, 5)

    class Context:
        @staticmethod
        def isOpenGLES():
            return False

        @staticmethod
        def format():
            return Format()

    monkeypatch.setattr(SP.QtGui.QOpenGLContext, "currentContext", lambda: Context())
    monkeypatch.setattr(SP.ogl_shaders, "compileShader",
                        lambda source, kind: calls.append(("shader", kind, source[:20])) or kind)
    monkeypatch.setattr(SP.ogl_shaders, "compileProgram",
                        lambda *compiled: calls.append(("program", compiled)) or 91)
    monkeypatch.setattr(SP.GL, "glBindAttribLocation",
                        lambda *args: calls.append(("attribute", *args)))
    monkeypatch.setattr(SP.GL, "glLinkProgram", lambda program: calls.append(("link", program)))
    program = SP._VboSphereProgram()
    assert program.program() == 91
    assert program.program() == 91, "the unchanged context recompiled the shader"
    assert [call[2:] for call in calls if call[0] == "attribute"] == [
        (0, "a_position"), (1, "a_color"), (2, "a_size")]
    assert len([call for call in calls if call[0] == "program"]) == 1


def test_vbo_program_declines_missing_and_opengles_contexts(monkeypatch):
    monkeypatch.setattr(SP.QtGui.QOpenGLContext, "currentContext", lambda: None)
    with pytest.raises(RuntimeError, match="no current"):
        SP._VboSphereProgram().program()

    class EsContext:
        @staticmethod
        def isOpenGLES():
            return True

    monkeypatch.setattr(SP.QtGui.QOpenGLContext, "currentContext", lambda: EsContext())
    with pytest.raises(RuntimeError, match="desktop OpenGL"):
        SP._VboSphereProgram().program()


def test_a_density_grid_is_transposed_for_opengl_and_uploaded_only_once(qapp):
    class Grid:
        occupancy = np.arange(24, dtype=float).reshape(2, 3, 4)
        lo = np.array([-1.0, -2.0, -3.0])
        cell = 0.5
        resolution = 4

    item = SP.ShadedScatter(pos=np.zeros((2, 3)), size=5.0)
    light = [{"pos": np.ones(3), "color": np.ones(3)}]
    assert item.set_scene("glossy 3D", light, np.ones(3), Grid())
    assert item._density_pending.shape == (4, 3, 2)
    assert item._density_pending[3, 2, 1] == Grid.occupancy[1, 2, 3]
    pending = item._density_pending
    assert item.set_scene("glossy 3D", light, np.ones(3), Grid())
    assert item._density_pending is pending, "the unchanged volume was prepared twice"
    assert not item.set_scene("flat", light, np.ones(3), Grid())
    item._gpu_failed = True
    assert not item.set_scene("metallic 3D", light, np.ones(3), Grid())
    assert not item.gpu_material_enabled("metallic 3D")


def test_density_upload_uses_a_float_3d_texture(qapp, monkeypatch):
    calls = []

    class FakeGL:
        GL_TEXTURE0, GL_TEXTURE1, GL_TEXTURE_3D = 0, 1, 3
        GL_TEXTURE_MIN_FILTER, GL_TEXTURE_MAG_FILTER = 4, 5
        GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_TEXTURE_WRAP_R = 6, 7, 8
        GL_LINEAR, GL_CLAMP_TO_EDGE = 9, 10
        GL_R32F, GL_RED, GL_FLOAT = 11, 12, 13
        glGenTextures = staticmethod(lambda n: 77)
        glActiveTexture = staticmethod(lambda unit: calls.append(("active", unit)))
        glBindTexture = staticmethod(lambda *a: calls.append(("bind", *a)))
        glTexParameteri = staticmethod(lambda *a: calls.append(("parameter", *a)))
        glTexImage3D = staticmethod(lambda *a: calls.append(("image", *a[:-1], a[-1].shape)))

    monkeypatch.setattr(SP, "GL", FakeGL)
    item = SP.ShadedScatter(pos=np.zeros((1, 3)), size=5.0)
    item._density_pending = np.ones((4, 3, 2), np.float32)
    assert item._flush_density()
    image = [call for call in calls if call[0] == "image"][0]
    assert image[4:7] == (2, 3, 4)
    assert calls[-1] == ("active", FakeGL.GL_TEXTURE0)
    assert item._density_pending is None and not item._flush_density()


def test_shader_uniforms_carry_material_lights_mood_and_viewport(qapp, monkeypatch):
    item = SP.ShadedScatter(pos=np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]), size=5.0)

    class View:
        def height(self):
            return 300

        def devicePixelRatioF(self):
            return 2.0

        def update(self):
            pass

    monkeypatch.setattr(item, "view", lambda: View())
    local = {"pos": np.array([1.0, 2.0, 3.0]), "color": np.array([0.4, 0.5, 0.6]),
             "local": True}
    item.set_scene("metallic 3D", [local], np.array([1.0, 0.7, 0.4]))
    values = item._shader_values()
    assert values["material"] == 2 and values["viewport"] == 600
    assert values["positions"][0, 3] == 1.0 and values["colors"][0, 3] == 1.0
    assert values["scale"] > 0 and values["count"] == 1
    item.set_scene("glossy 3D", [], np.ones(3))
    fallback = item._shader_values()
    assert fallback["material"] == 1 and fallback["count"] == 1
    assert np.any(fallback["positions"][0, :3])


def test_typed_shader_program_uploads_every_uniform(monkeypatch):
    calls = []
    monkeypatch.setattr(SP.shaders.ShaderProgram, "__enter__", lambda self: self)
    program = SP._SphereProgram()
    program.uniform = lambda name: name
    assert program.__enter__() is program                         # no values yet
    program.values = {
        "count": 1, "ray": 1, "material": 2, "extent": 48.0, "cell": 1.0,
        "absorb": 0.55, "viewport": 800.0, "scale": 20.0, "lo": np.zeros(3),
            "mood": np.ones(3), "positions": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
            "colors": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
            "directions": np.zeros((SP.MAX_LIGHTS, 3), np.float32),
            "controls": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
    }

    class FakeGL:
        glUniform1i = staticmethod(lambda *a: calls.append(("i", *a)))
        glUniform1f = staticmethod(lambda *a: calls.append(("f", *a)))
        glUniform3fv = staticmethod(lambda *a: calls.append(("3", *a)))
        glUniform4fv = staticmethod(lambda *a: calls.append(("4", *a)))

    monkeypatch.setattr(SP, "GL", FakeGL)
    assert program.__enter__() is program
    assert {call[1] for call in calls} >= {"uLightCount", "uRayEnabled", "uMaterial",
                                           "uDensity", "uScale", "uMood", "uLightPos"}


def test_raw_vbo_program_uploads_every_typed_uniform(monkeypatch):
    calls = []

    class FakeGL:
        glGetUniformLocation = staticmethod(lambda program, name: name)
        glUniform1i = staticmethod(lambda *a: calls.append(("i", *a)))
        glUniform1f = staticmethod(lambda *a: calls.append(("f", *a)))
        glUniform3fv = staticmethod(lambda *a: calls.append(("3", *a)))
        glUniform4fv = staticmethod(lambda *a: calls.append(("4", *a)))

    values = {
        "count": 1, "ray": 1, "material": 2, "extent": 48.0, "cell": 1.0,
        "absorb": 0.55, "viewport": 800.0, "scale": 20.0, "lo": np.zeros(3),
        "mood": np.ones(3), "positions": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
        "colors": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
        "directions": np.zeros((SP.MAX_LIGHTS, 3), np.float32),
        "controls": np.zeros((SP.MAX_LIGHTS, 4), np.float32),
    }
    monkeypatch.setattr(SP, "GL", FakeGL)
    SP._upload_shader_values(17, values)
    assert {call[1] for call in calls} >= {
        "uLightCount", "uRayEnabled", "uMaterial", "uDensity", "uDensityExtent",
        "uCell", "uAbsorb", "uViewportHeight", "uScale", "uDensityLo", "uMood",
        "uLightPos", "uLightColor",
    }


def test_pyqtgraph_014_get_shader_seam_prepares_material_and_keeps_flat_stock(qapp,
                                                                              monkeypatch):
    monkeypatch.setattr(SP.gl.GLScatterPlotItem, "getShaderProgram", lambda self: "stock",
                        raising=False)
    item = SP.ShadedScatter(pos=np.ones((2, 3)), size=5.0)
    assert item._vbo_renderer and item.getShaderProgram() == "stock"

    class Program:
        @staticmethod
        def program():
            return 73

    class Matrix:
        @staticmethod
        def data():
            return np.eye(4, dtype=np.float32).ravel()

    calls = []
    item._pbr = Program()
    item.set_scene("metallic 3D", [], np.ones(3))
    monkeypatch.setattr(item, "projectionMatrix", lambda: Matrix(), raising=False)
    monkeypatch.setattr(item, "_shader_values", lambda: {"ready": True})
    monkeypatch.setattr(SP, "_upload_shader_values",
                        lambda program, values: calls.append((program, values)))
    monkeypatch.setattr(SP.GL, "glUseProgram", lambda program: calls.append(("use", program)))
    monkeypatch.setattr(SP.GL, "glGetUniformLocation", lambda *_: 8)
    monkeypatch.setattr(SP.GL, "glUniformMatrix4fv",
                        lambda *args: calls.append(("matrix", *args)))
    assert item.getShaderProgram() == 73
    assert calls[0] == ("use", 73) and calls[1] == (73, {"ready": True})
    assert calls[-1] == ("use", 0)


def test_paint_uses_pbr_and_falls_back_once_if_the_driver_rejects_it(qapp, monkeypatch, capsys):
    """The 0.13 path, chosen explicitly rather than inherited from whatever is installed.

    `shader` is an attribute of pyqtgraph 0.13's fixed-function scatter and 0.14 removed it, so a
    test that reads it is a test of the older renderer -- and one that only says so by accident
    fails on an upgrade for a reason that has nothing to do with the behaviour being checked. The
    0.14 seam is `getShaderProgram`, tested above.
    """
    painted = []

    def parent_paint(item):
        painted.append(item.shader)
        if item.shader is pbr and len(painted) == 1:
            raise RuntimeError("old driver")

    monkeypatch.setattr(SP.gl.GLScatterPlotItem, "paint", parent_paint)
    monkeypatch.setattr(SP.GL, "glActiveTexture", lambda *_: None)
    monkeypatch.setattr(SP.GL, "glBindTexture", lambda *_: None)
    item = SP.ShadedScatter(pos=np.ones((2, 3)), size=5.0)
    item._vbo_renderer = False
    pbr = SP._SphereProgram()
    stock = object()
    item._pbr = pbr
    item._stock_shader = stock
    item._density_texture = 5
    item.set_scene("metallic 3D", [], np.ones(3))
    monkeypatch.setattr(item, "_shader_values", lambda: {"ready": True})
    item.paint()
    assert painted[0] is pbr and painted[1] is stock
    assert item._gpu_failed and "GPU sphere shader unavailable" in capsys.readouterr().out
    item.paint()
    assert painted[-1] is stock


def test_flat_fallback_never_looks_up_the_removed_point_sprite_name(qapp, monkeypatch):
    """Regression for pyqtgraph 0.14: the name is absent and caused one exception every frame."""
    painted = []
    monkeypatch.setattr(SP.gl.GLScatterPlotItem, "paint", lambda item: painted.append(item.shader))
    monkeypatch.setattr(SP.shaders, "getShaderProgram",
                        lambda *_: pytest.fail("private pointSprite registry was consulted"))
    item = SP.ShadedScatter(pos=np.ones((2, 3)), size=5.0)
    item._vbo_renderer = False          # the renderer whose scatter HAS a `shader` attribute
    stock = object()
    item._stock_shader = stock
    item.paint()
    assert painted == [stock]


def test_a_spotlight_packs_its_cone_into_the_shader_uniforms(qapp):
    """A flashlight is a cone, and the shader needs it as cosines rather than degrees -- comparing
    an angle to a dot product every fragment would be the same arithmetic done a million times a
    frame instead of twice here. The plain lights share the array and must leave the cone off."""
    item = SP.ShadedScatter(pos=np.zeros((3, 3)), size=5.0)
    item.set_scene("glossy 3D", mood=np.ones(3), lights=[
        {"pos": np.array([1.0, 2.0, 3.0]), "color": np.array([1.0, 0.9, 0.8]),
         "spot": True, "direction": np.array([0.0, 0.0, -1.0]),
         "inner": 18.0, "outer": 34.0, "gain": 1.25, "local": True},
        {"pos": np.array([9.0, 0.0, 0.0]), "color": np.ones(3)},
    ])
    v = item._shader_values()
    controls = np.asarray(v["controls"]).reshape(-1, 4)
    assert controls[0, 0] == 1.0, "the spotlight was not flagged as one"
    assert controls[0, 1] == pytest.approx(np.cos(np.radians(18.0)))
    assert controls[0, 2] == pytest.approx(np.cos(np.radians(34.0)))
    assert controls[0, 3] == pytest.approx(1.25)
    assert controls[1, 0] == 0.0, "a plain light was given a cone"
    assert controls[1, 3] == pytest.approx(1.0), "a light with no gain did not default to one"
    assert np.asarray(v["positions"]).reshape(-1, 4)[0, 3] == 1.0, "local flag lost"
    assert v["count"] == 2 and np.allclose(v["directions"][0], (0.0, 0.0, -1.0))
