#!/usr/bin/env python3
"""GPU sphere-impostor materials for the genes in the OpenGL scatter.

## Why a fragment shader

pyqtgraph draws a scatter as point sprites: one square of texture per gene, always facing the
viewer, with a vertex shader that does nothing but set the size and pass the colour through. There
is no fragment shader. The texture it ships is `pData[:] = 255` with an alpha disc cut out of it --
every pixel inside a ball carries the identical colour, which is a flat disc, and no amount of
per-point shading can make a flat disc look spherical. That is the whole of "the balls still look
like 2D matt balls": the per-point lighting was working and had nowhere to show up, because a ball
is one colour and a sphere is a gradient.

Glossy and metallic modes replace that path with a compatibility-profile GLSL program. Each fragment
uses `gl_PointCoord` to reconstruct the visible hemisphere's normal, evaluates a GGX microfacet
response, reflects a soft studio environment, and writes the sphere surface to `gl_FragDepth`.
Unlike a pre-lit texture, the highlight is continuous, follows world-space lights, and changes with
view direction. Metallic reflection is tinted by the gene's data color; glossy reflection remains a
dielectric white, so the two materials differ for a physical reason rather than by brightness alone.

Ray mode uploads the cached occupancy grid as a 3D float texture. The vertex shader marches a fixed
24-sample shadow ray from each gene to each active light and passes one stable visibility value to
all fragments of that sphere. This is GPU-accelerated volumetric ray tracing, not hardware triangle
ray tracing: a gene cloud has density but no triangle surfaces, BVH, reflection rays, or refraction.

The 64×64 CPU textures remain solely as a safe fallback for old drivers and as the flat-disc path.
"""
from __future__ import annotations

import numpy as np
import pyqtgraph.opengl as gl
from pyqtgraph.opengl import shaders
from OpenGL import GL
from OpenGL.GL import shaders as ogl_shaders
from pyqtgraph.Qt import QtGui

#: How each finish answers the light ACROSS ONE BALL. Separate numbers from `lighting.FINISHES`,
#: which says how the cloud answers it, and they are separate questions: a tight highlight is
#: useless spread over a scatter of points and is exactly right inside a single sphere, where there
#: are pixels for it to land on. Shininess here is the textbook range for that reason.
#:
#: "2D" is the plain flat disc -- pyqtgraph's own sprite, and what this map has always drawn. Kept
#: as a choice because a flat disc is the honest way to show a scatter where the reader is comparing
#: colours and a highlight is one more thing in the way.
#: One entry per finish in `lighting.POINT_MODES`, and the same five for the same reason: three of the
#: original eight measured indistinguishable from a neighbour. The note there carries the numbers.
SPRITES = {
    "flat": None,
    "velvet 3D": {"ambient": 0.24, "diffuse": 0.78, "specular": 0.22,
                  "shininess": 5.0, "rim": 0.92},
    "glossy 3D": {"ambient": 0.16, "diffuse": 0.72, "specular": 1.20,
                  "shininess": 64.0, "rim": 0.18},
    "glass 3D": {"ambient": 0.10, "diffuse": 0.34, "specular": 1.90,
                 "shininess": 78.0, "rim": 1.15},
    "metallic 3D": {"ambient": 0.10, "diffuse": 0.50, "specular": 1.65,
                    "shininess": 26.0, "rim": 0.75},
}

LEGACY_POINT_MODES = {
    "2D": "flat", "matt": "flat", "satin": "glossy 3D",
    "glossy": "glossy 3D", "metallic": "metallic 3D",
    "pearl 3D": "glossy 3D", "brushed metal 3D": "metallic 3D", "silver 3D": "metallic 3D",
}

#: Texels across one ball. 64 is what pyqtgraph uses, and a gene is drawn at 5 to 12 pixels: the
#: texture is downsampled hard, so a bigger one buys nothing and a smaller one shows its own edge.
WIDTH = 64

#: Where the light comes from when nothing else says -- above and to the left, which is where every
#: reader since the Renaissance has assumed light comes from, and reads as convex rather than
#: hollow. z is toward the viewer.
DEFAULT_LIGHT = (-0.4, 0.6, 0.7)

# The GL compatibility profile used by pyqtgraph exposes eight varyings comfortably on every
# supported card. Eight local lights is also enough for a selected gene and a readable sample of its
# neighbors; more lights flatten the material into uniform brightness rather than adding detail.
MAX_LIGHTS = 8
RAY_STEPS = 24

#: The id the fragment shader branches on. Renumbered when the duplicate finishes were retired --
#: safe to renumber because nothing stores it: it is looked up from the finish NAME every frame, and
#: the names are what persist in settings and recipes.
MATERIAL_IDS = {
    "glossy 3D": 1, "metallic 3D": 2, "glass 3D": 3, "velvet 3D": 4,
}


_SPHERE_VERTEX = f"""
#version 120
const int MAX_LIGHTS = {MAX_LIGHTS};
const int RAY_STEPS = {RAY_STEPS};
uniform int uLightCount;
uniform int uRayEnabled;
uniform vec4 uLightPos[MAX_LIGHTS];
uniform vec3 uLightDir[MAX_LIGHTS];
uniform vec4 uLightControl[MAX_LIGHTS];
uniform sampler3D uDensity;
uniform vec3 uDensityLo;
uniform float uDensityExtent;
uniform float uCell;
uniform float uAbsorb;
uniform float uViewportHeight;
uniform float uScale;
varying vec4 vColor;
varying vec3 vCenterEye;
varying vec4 vLightEye[MAX_LIGHTS];
varying float vRadiusEye;

void main() {{
    gl_Position = ftransform();
    gl_PointSize = gl_Normal.x;
    gl_FrontColor = gl_Color;
    vColor = gl_Color;
    vec4 eye = gl_ModelViewMatrix * gl_Vertex;
    vCenterEye = eye.xyz;
    vRadiusEye = gl_PointSize * max(-eye.z, 0.001) /
                 max(uViewportHeight * gl_ProjectionMatrix[1][1], 0.001);
    for (int light = 0; light < MAX_LIGHTS; ++light) {{
        vec3 lightEye = (gl_ModelViewMatrix * vec4(uLightPos[light].xyz, 1.0)).xyz;
        float transport = 1.0;
        if (light < uLightCount) {{
            float shadow = 1.0;
            if (uRayEnabled != 0 && uDensityExtent > 0.0) {{
                vec3 reach = uLightPos[light].xyz - gl_Vertex.xyz;
                float distanceCells = length(reach) / max(uCell, 0.000001);
                float density = 0.0;
                for (int rayStep = 1; rayStep <= RAY_STEPS; ++rayStep) {{
                    float t = float(rayStep) / float(RAY_STEPS + 1);
                    vec3 uv = (gl_Vertex.xyz + reach * t - uDensityLo) / uDensityExtent;
                    if (all(greaterThanEqual(uv, vec3(0.0))) && all(lessThanEqual(uv, vec3(1.0))))
                        density += texture3D(uDensity, uv).r;
                }}
                shadow = exp(-uAbsorb * distanceCells * density / float(RAY_STEPS + 1));
            }}
            float beam = 1.0;
            if (uLightControl[light].x > 0.5) {{
                vec3 fromLight = normalize(gl_Vertex.xyz - uLightPos[light].xyz);
                float cosine = dot(fromLight, normalize(uLightDir[light]));
                beam = smoothstep(uLightControl[light].z, uLightControl[light].y, cosine);
            }}
            transport = shadow * beam * uLightControl[light].w;
        }}
        vLightEye[light] = vec4(lightEye, transport);
    }}
}}
"""


_SPHERE_FRAGMENT = f"""
#version 120
const int MAX_LIGHTS = {MAX_LIGHTS};
uniform int uLightCount;
uniform int uMaterial;
uniform vec4 uLightPos[MAX_LIGHTS];
uniform vec4 uLightColor[MAX_LIGHTS];
uniform vec3 uMood;
uniform float uScale;
varying vec4 vColor;
varying vec3 vCenterEye;
varying vec4 vLightEye[MAX_LIGHTS];
varying float vRadiusEye;

float distributionGGX(float nDotH, float roughness) {{
    float a = roughness * roughness;
    float a2 = a * a;
    float d = nDotH * nDotH * (a2 - 1.0) + 1.0;
    return a2 / max(3.14159265 * d * d, 0.0001);
}}

float geometrySchlick(float nDotV, float roughness) {{
    float r = roughness + 1.0;
    float k = r * r / 8.0;
    return nDotV / max(nDotV * (1.0 - k) + k, 0.0001);
}}

vec3 fresnelSchlick(float cosine, vec3 f0) {{
    return f0 + (1.0 - f0) * pow(1.0 - clamp(cosine, 0.0, 1.0), 5.0);
}}

void main() {{
    vec2 q = gl_PointCoord * 2.0 - 1.0;
    q.y = -q.y;
    float radius2 = dot(q, q);
    if (radius2 > 1.0) discard;
    float z = sqrt(max(1.0 - radius2, 0.0));
    vec3 normal = normalize(vec3(q, z));
    vec3 surfaceEye = vCenterEye + normal * vRadiusEye;
    vec4 surfaceClip = gl_ProjectionMatrix * vec4(surfaceEye, 1.0);
    gl_FragDepth = 0.5 * (surfaceClip.z / surfaceClip.w) + 0.5;

    vec3 view = normalize(-surfaceEye);
    float nDotV = max(dot(normal, view), 0.001);
    // Four ids, one per material class. The roughness spread is deliberately wide: the retired
    // `pearl 3D` sat at the SAME roughness as glossy and differed only in base reflectance, which is
    // why the two were indistinguishable on screen. 0.07 / 0.24 / 0.72 is far enough apart that the
    // highlight is a point, a coin and a whole hemisphere respectively.
    bool metal = uMaterial == 2;
    bool glass = uMaterial == 3;
    bool velvet = uMaterial == 4;
    float metallic = metal ? 1.0 : 0.0;
    float roughness = metal ? 0.18 : (velvet ? 0.72 : (glass ? 0.07 : 0.24));
    vec3 albedo = clamp(vColor.rgb, 0.0, 1.0);
    vec3 f0 = mix(vec3(0.045), albedo, metallic);
    vec3 direct = vec3(0.0);

    for (int light = 0; light < MAX_LIGHTS; ++light) {{
        if (light < uLightCount) {{
            vec3 toLight = vLightEye[light].xyz - surfaceEye;
            float lightDistance = length(toLight);
            vec3 lightDir = normalize(toLight);
            vec3 halfDir = normalize(view + lightDir);
            float nDotL = max(dot(normal, lightDir), 0.0);
            float nDotH = max(dot(normal, halfDir), 0.0);
            float hDotV = max(dot(halfDir, view), 0.0);
            float d = distributionGGX(nDotH, roughness);
            float g = geometrySchlick(nDotV, roughness) * geometrySchlick(nDotL, roughness);
            vec3 f = fresnelSchlick(hDotV, f0);
            vec3 specular = d * g * f / max(4.0 * nDotV * nDotL, 0.001);
            vec3 diffuse = (1.0 - f) * (1.0 - metallic) * albedo / 3.14159265;
            float attenuation = uLightPos[light].w > 0.5
                ? (0.18 + 1.35 / (1.0 + pow(lightDistance / max(uScale * 0.46, 0.001), 2.0)))
                : 1.0;
            float localGlow = uLightPos[light].w > 0.5
                ? 0.42 / (1.0 + pow(lightDistance / max(uScale * 0.22, 0.001), 2.0)) : 0.0;
            direct += ((diffuse + specular) * nDotL * attenuation + albedo * localGlow) *
                      uLightColor[light].rgb * vLightEye[light].w;
        }}
    }}

    // A simple studio environment is what gives metal something to reflect. A single point light
    // makes even a correct metal shader look like dark plastic everywhere except one tiny glint.
    vec3 reflected = reflect(-view, normal);
    vec3 sky = vec3(0.32, 0.45, 0.68) * uMood;
    vec3 ground = vec3(0.055, 0.040, 0.030) * uMood;
    vec3 environment = mix(ground, sky, smoothstep(-0.35, 0.75, reflected.y));
    // Two broad softboxes: compact lobes rather than the old infinite vertical strip, which drew
    // an actual line through every metallic point and a dot that followed the cursor.
    // The softbox tightness follows the roughness: a sharp material reflects a small bright box, a
    // rough one smears it across the ball. Glass gets the tightest because that is the whole of what
    // separates it from glossy at this size.
    float boxA = pow(max(dot(reflected, normalize(vec3(-0.48, 0.62, 0.62))), 0.0),
                     glass ? 40.0 : (velvet ? 4.0 : 16.0));
    float boxB = pow(max(dot(reflected, normalize(vec3(0.70, 0.18, 0.69))), 0.0),
                     glass ? 48.0 : (velvet ? 3.0 : 20.0));
    environment += vec3(1.0, 0.88, 0.70) * boxA * (metal ? 1.30 : (glass ? 0.72 : 0.38));
    environment += vec3(0.55, 0.72, 1.0) * boxB * (metal ? 0.82 : (glass ? 0.46 : 0.22));
    vec3 envFresnel = fresnelSchlick(nDotV, f0);
    vec3 moodFill = 0.35 + 0.65 * uMood;
    vec3 diffuseFill = albedo * moodFill * (metal ? 0.09 :
                       (velvet ? 0.46 : (0.30 + 0.18 * max(normal.y, 0.0))));
    vec3 color = diffuseFill + direct * (metal ? 1.35 : 1.0) +
                 environment * envFresnel * (metal ? 1.20 : (glass ? 0.78 : 0.48));
    // A dark body with a bright edge: the reading that says "transparent" at 8 pixels across.
    if (glass) color = mix(albedo * 0.22, color, 0.80) + envFresnel * 0.32;
    // Retroreflective sheen -- brightest where the ball turns away, which is what makes cloth read as
    // cloth. Strengthened when the near-duplicates went, since this is now the only matte 3D finish.
    if (velvet) color += albedo * pow(1.0 - nDotV, 1.8) * 0.62;
    color = vec3(1.0) - exp(-color * (metal ? 1.75 : 2.05));
    color = clamp(color, 0.0, 1.0);

    float edge = 1.0 - smoothstep(0.88, 1.0, radius2);
    gl_FragColor = vec4(color, vColor.a * edge);
}}
"""


class _SphereProgram(shaders.ShaderProgram):
    """Compatibility-profile PBR shader with typed uniforms and a 3D density sampler."""

    def __init__(self):
        super().__init__("starplastPbrSpheres",
                         [shaders.VertexShader(_SPHERE_VERTEX),
                          shaders.FragmentShader(_SPHERE_FRAGMENT)])
        self.values = None

    def __enter__(self):
        super().__enter__()
        values = self.values
        if values is None:
            return self
        one_i = lambda name, value: GL.glUniform1i(self.uniform(name), int(value))
        one_f = lambda name, value: GL.glUniform1f(self.uniform(name), float(value))
        one_i("uLightCount", values["count"])
        one_i("uRayEnabled", values["ray"])
        one_i("uMaterial", values["material"])
        one_i("uDensity", 1)
        one_f("uDensityExtent", values["extent"])
        one_f("uCell", values["cell"])
        one_f("uAbsorb", values["absorb"])
        one_f("uViewportHeight", values["viewport"])
        one_f("uScale", values["scale"])
        GL.glUniform3fv(self.uniform("uDensityLo"), 1, values["lo"])
        GL.glUniform3fv(self.uniform("uMood"), 1, values["mood"])
        GL.glUniform4fv(self.uniform("uLightPos"), MAX_LIGHTS, values["positions"])
        GL.glUniform4fv(self.uniform("uLightColor"), MAX_LIGHTS, values["colors"])
        GL.glUniform3fv(self.uniform("uLightDir"), MAX_LIGHTS, values["directions"])
        GL.glUniform4fv(self.uniform("uLightControl"), MAX_LIGHTS, values["controls"])
        return self


def _vbo_shader_sources(core: bool) -> tuple[str, str]:
    """Translate the compatibility shader for pyqtgraph 0.14's VBO renderer.

    pyqtgraph 0.13 feeds the built-in ``gl_Vertex``, ``gl_Color`` and ``gl_Normal`` values.
    Version 0.14 deliberately stopped doing that and binds three explicit attributes instead.  A
    shader compiled for one API cannot draw the other, even though both classes have the same name.
    Keeping the material body identical here prevents the two supported renderers from slowly
    acquiring different-looking metal and gloss.
    """
    vertex = _SPHERE_VERTEX
    vertex = vertex.replace("gl_Position = ftransform();", "gl_Position = u_mvp * a_position;")
    vertex = vertex.replace("gl_PointSize = gl_Normal.x;", "gl_PointSize = a_size;")
    vertex = vertex.replace("    gl_FrontColor = gl_Color;\n", "")
    vertex = vertex.replace("vColor = gl_Color;", "vColor = a_color;")
    vertex = vertex.replace("gl_ModelViewMatrix", "u_modelview")
    vertex = vertex.replace("gl_ProjectionMatrix", "u_projection")
    vertex = vertex.replace("gl_Vertex", "a_position")
    declarations = """
uniform mat4 u_mvp;
uniform mat4 u_modelview;
uniform mat4 u_projection;
attribute vec4 a_position;
attribute vec4 a_color;
attribute float a_size;
"""
    vertex = vertex.replace("#version 120\n", "#version 120\n" + declarations, 1)

    fragment = _SPHERE_FRAGMENT.replace("gl_ProjectionMatrix", "u_projection")
    fragment = fragment.replace("#version 120\n", "#version 120\nuniform mat4 u_projection;\n", 1)
    if core:
        vertex = vertex.replace("#version 120", "#version 140", 1)
        vertex = vertex.replace("attribute ", "in ").replace("varying ", "out ")
        vertex = vertex.replace("texture3D(", "texture(")
        fragment = fragment.replace("#version 120", "#version 140", 1)
        fragment = fragment.replace("varying ", "in ")
        fragment = fragment.replace("void main() {", "out vec4 fragColor;\n\nvoid main() {", 1)
        fragment = fragment.replace("gl_FragColor", "fragColor")
    return vertex, fragment


class _VboSphereProgram:
    """Lazily compiled raw program for pyqtgraph 0.14 and later."""

    def __init__(self):
        self._program = None

    def program(self):
        if self._program is not None:
            return self._program
        context = QtGui.QOpenGLContext.currentContext()
        if context is None:
            raise RuntimeError("no current OpenGL context")
        if context.isOpenGLES():
            # The desktop application requests desktop OpenGL.  Falling back on unusual ES-only
            # systems is preferable to compiling desktop sampler3D/depth syntax as an ES shader.
            raise RuntimeError("the 3D material needs a desktop OpenGL context")
        vertex, fragment = _vbo_shader_sources(context.format().version() >= (3, 1))
        compiled = [
            ogl_shaders.compileShader(vertex, GL.GL_VERTEX_SHADER),
            ogl_shaders.compileShader(fragment, GL.GL_FRAGMENT_SHADER),
        ]
        program = ogl_shaders.compileProgram(*compiled)
        GL.glBindAttribLocation(program, 0, "a_position")
        GL.glBindAttribLocation(program, 1, "a_color")
        GL.glBindAttribLocation(program, 2, "a_size")
        GL.glLinkProgram(program)
        self._program = program
        return program


def _upload_shader_values(program, values: dict) -> None:
    """Upload the material uniforms to an already-bound raw OpenGL program."""
    uniform = lambda name: GL.glGetUniformLocation(program, name)
    GL.glUniform1i(uniform("uLightCount"), int(values["count"]))
    GL.glUniform1i(uniform("uRayEnabled"), int(values["ray"]))
    GL.glUniform1i(uniform("uMaterial"), int(values["material"]))
    GL.glUniform1i(uniform("uDensity"), 1)
    float_uniforms = {
        "extent": "uDensityExtent", "cell": "uCell", "absorb": "uAbsorb",
        "viewport": "uViewportHeight", "scale": "uScale",
    }
    for key, name in float_uniforms.items():
        GL.glUniform1f(uniform(name), float(values[key]))
    GL.glUniform3fv(uniform("uDensityLo"), 1, values["lo"])
    GL.glUniform3fv(uniform("uMood"), 1, values["mood"])
    GL.glUniform4fv(uniform("uLightPos"), MAX_LIGHTS, values["positions"])
    GL.glUniform4fv(uniform("uLightColor"), MAX_LIGHTS, values["colors"])
    GL.glUniform3fv(uniform("uLightDir"), MAX_LIGHTS, values["directions"])
    GL.glUniform4fv(uniform("uLightControl"), MAX_LIGHTS, values["controls"])


def disc_alpha(w: int = WIDTH) -> np.ndarray:
    """The round edge, antialiased over the last texel. pyqtgraph's own formula.

    Matched deliberately: a ball with a different edge to the halos and centroids drawn beside it
    reads as a different KIND of thing rather than the same thing lit.
    """
    y, x = np.mgrid[0:w, 0:w]
    r = np.hypot(x - (w - 1) / 2.0, y - (w - 1) / 2.0)
    return np.clip(w / 2 - np.clip(r, w / 2 - 1, w / 2), 0.0, 1.0)


def texture(finish: str, light=DEFAULT_LIGHT, w: int = WIDTH) -> np.ndarray:
    """One ball, lit from `light`, as (w, w, 4) uint8 ready for GL.

    `light` is a direction in the screen's frame: x right, y up, z toward the viewer. The surface
    normal is read off the disc -- a texel at (u, v) from the centre of a unit sphere has normal
    (u, v, sqrt(1 - u^2 - v^2)), which is the whole trick and costs one square root per texel.
    """
    alpha = disc_alpha(w)
    finish = LEGACY_POINT_MODES.get(str(finish), str(finish))
    f = SPRITES.get(finish)
    if f is None:                       # flat, or a name nobody defined: the honest fallback
        out = np.empty((w, w, 4), dtype=np.ubyte)
        out[:, :, :3] = 255
        out[:, :, 3] = (alpha * 255).astype(np.ubyte)
        return out

    span = np.linspace(-1.0, 1.0, w)
    u, v = np.meshgrid(span, -span)     # -span: texture rows run down, the sphere's y runs up
    z2 = 1.0 - u * u - v * v
    inside = z2 > 0.0
    nz = np.sqrt(np.where(inside, z2, 0.0))

    d = np.asarray(light, dtype=float)
    n = np.linalg.norm(d)
    if n <= 1e-9:                       # a light with no direction: use the one above and left
        d = np.asarray(DEFAULT_LIGHT, dtype=float)
        n = np.linalg.norm(d)
    d = d / n
    lam = np.clip(u * d[0] + v * d[1] + nz * d[2], 0.0, None)

    # Blinn-Phong: the halfway vector between the light and the viewer, who is straight ahead at
    # (0, 0, 1) because a point sprite always faces them.
    half = d + np.array([0.0, 0.0, 1.0])
    half = half / max(float(np.linalg.norm(half)), 1e-9)
    spec = np.clip(u * half[0] + v * half[1] + nz * half[2], 0.0, None) ** f["shininess"]

    # The limb, where the surface turns away from the viewer. On a real sphere this is where a
    # glancing reflection of everything else in the room shows up; here it is what stops a ball
    # ending in a hard circular edge, and it is most of what a metal looks like.
    rim = (1.0 - nz) ** 3 * np.clip(lam + 0.35, 0.0, 1.0)

    shade = f["ambient"] + f["diffuse"] * lam + f["specular"] * spec + f["rim"] * rim
    shade = np.clip(np.where(inside, shade, f["ambient"]), 0.0, 1.0)

    out = np.empty((w, w, 4), dtype=np.ubyte)
    out[:, :, 0] = out[:, :, 1] = out[:, :, 2] = (shade * 255).astype(np.ubyte)
    out[:, :, 3] = (alpha * 255).astype(np.ubyte)
    return out


def upload(texture_id, data) -> None:
    """Hand a texture to the card. Only ever called with a current context -- see `ShadedScatter`."""
    GL.glBindTexture(GL.GL_TEXTURE_2D, texture_id)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, data.shape[0], data.shape[1], 0,
                    GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, data)


class ShadedScatter(gl.GLScatterPlotItem):
    """A scatter whose glossy modes are GPU-shaded sphere impostors.

    The old CPU sphere texture remains the safe fallback. With a working compatibility-profile
    shader, every fragment reconstructs a sphere normal, computes a GGX glossy/metallic response,
    and writes the sphere surface's depth. Ray mode samples a 3D density texture in the vertex shader,
    so one stable shadow value is computed per gene on the GPU rather than on every CPU timer tick.
    """

    def __init__(self, **kwds):
        super().__init__(**kwds)
        # 0.14 owns a raw VBO shader returned by getShaderProgram(); 0.13 owns a compatibility
        # ShaderProgram in self.shader.  The public class name stayed unchanged across that break.
        self._vbo_renderer = hasattr(gl.GLScatterPlotItem, "getShaderProgram")
        self._pending = None
        self._density_pending = None
        self._density_texture = None
        self._density_key = None
        self._density_meta = (np.zeros(3, np.float32), 1.0, 1.0)
        self._scene = {"mode": "flat", "lights": [], "mood": np.ones(3), "ray": False,
                       "absorb": 0.48}
        self._stock_shader = None
        self._pbr = _VboSphereProgram() if self._vbo_renderer else _SphereProgram()
        self._gpu_failed = False

    def set_sprite(self, data) -> None:
        """Set the CPU texture used by flat mode and by the safe shader fallback."""
        self._pending = data
        self.update()

    def set_scene(self, mode: str, lights, mood, grid=None, absorb: float = 0.48) -> bool:
        """Send material, lights and an optional density volume to the next GPU paint.

        Returns whether ray transport will run on the GPU. Callers use the false result to retain
        the existing CPU shadow fallback on old or failed OpenGL contexts.
        """
        self._scene = {"mode": str(mode), "lights": list(lights or [])[:MAX_LIGHTS],
                       "mood": np.asarray(mood, dtype=np.float32), "ray": grid is not None,
                       "absorb": float(absorb)}
        if grid is not None:
            key = id(grid.occupancy)
            if key != self._density_key:
                self._density_key = key
                # OpenGL's x coordinate is the fastest-moving texel; NumPy's last axis is fastest.
                # Grid stores [x, y, z], so upload [z, y, x] or the shadow volume is x/z-transposed.
                self._density_pending = np.ascontiguousarray(
                    np.transpose(grid.occupancy, (2, 1, 0)), dtype=np.float32)
                extent = float(grid.cell * grid.resolution)
                self._density_meta = (np.asarray(grid.lo, dtype=np.float32), extent,
                                      float(grid.cell))
        self.update()
        return bool(mode != "flat" and grid is not None and not self._gpu_failed)

    def gpu_material_enabled(self, mode: str) -> bool:
        """Whether `mode` is handled by the per-fragment sphere shader on the next paint."""
        return bool(mode != "flat" and not self._gpu_failed)

    def _flush(self) -> bool:
        if self._pending is None or getattr(self, "pointTexture", None) is None:
            return False
        upload(self.pointTexture, self._pending)
        self._pending = None
        return True

    def _flush_density(self) -> bool:
        """Upload the cached occupancy field while the view's GL context is current."""
        if self._density_pending is None:
            return False
        if self._density_texture is None:
            self._density_texture = GL.glGenTextures(1)
        data = self._density_pending
        GL.glActiveTexture(GL.GL_TEXTURE1)
        GL.glBindTexture(GL.GL_TEXTURE_3D, self._density_texture)
        GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_EDGE)
        GL.glTexImage3D(GL.GL_TEXTURE_3D, 0, GL.GL_R32F, data.shape[2], data.shape[1],
                        data.shape[0], 0, GL.GL_RED, GL.GL_FLOAT, data)
        GL.glActiveTexture(GL.GL_TEXTURE0)
        self._density_pending = None
        return True

    def _shader_values(self) -> dict:
        """Fixed-size uniform arrays for the current scene."""
        lights = self._scene["lights"]
        positions = np.zeros((MAX_LIGHTS, 4), dtype=np.float32)
        colors = np.zeros((MAX_LIGHTS, 4), dtype=np.float32)
        directions = np.zeros((MAX_LIGHTS, 3), dtype=np.float32)
        controls = np.zeros((MAX_LIGHTS, 4), dtype=np.float32)
        controls[:, 3] = 1.0
        for i, light in enumerate(lights):
            positions[i, :3] = light["pos"]
            positions[i, 3] = 1.0 if light.get("local") else 0.0
            colors[i, :3], colors[i, 3] = light["color"], 1.0
            if light.get("spot"):
                directions[i] = light["direction"]
                controls[i, 0] = 1.0
                controls[i, 1] = np.cos(np.radians(float(light["inner"])))
                controls[i, 2] = np.cos(np.radians(float(light["outer"])))
            controls[i, 3] = float(light.get("gain", 1.0))
        if not lights:
            positions[0] = (-40.0, 60.0, 90.0, 1.0)
            colors[0] = (0.92, 0.92, 0.90, 1.0)
            count = 1
        else:
            count = len(lights)
        lo, extent, cell = self._density_meta
        view = self.view()
        viewport = max(float(view.height() * view.devicePixelRatioF()), 1.0) if view else 1.0
        points = np.asarray(self.pos, dtype=float)
        scale = float(np.percentile(np.linalg.norm(points - points.mean(0), axis=1), 95)) \
            if len(points) else 1.0
        return {"count": count, "ray": int(self._scene["ray"]),
                "material": MATERIAL_IDS.get(self._scene["mode"], 1),
                "positions": positions, "colors": colors, "mood": self._scene["mood"],
                "directions": directions, "controls": controls,
                "lo": lo, "extent": extent, "cell": cell, "absorb": self._scene["absorb"],
                "viewport": viewport, "scale": max(scale, 1e-6)}

    def initializeGL(self):
        """Let pyqtgraph initialise its renderer and remember the 0.13 fallback shader."""
        super().initializeGL()
        self._flush()
        if not self._vbo_renderer:
            self._stock_shader = getattr(self, "shader", None)

    def getShaderProgram(self):
        """Return and prepare the VBO material program used by pyqtgraph 0.14.

        Its base paint method owns the VBOs and the model/view matrices, so replacing this one
        seam preserves its renderer while supplying Starplast's per-fragment sphere material.
        """
        if (not self._vbo_renderer or self._scene["mode"] == "flat" or self._gpu_failed):
            return super().getShaderProgram()
        program = self._pbr.program()
        GL.glUseProgram(program)
        try:
            _upload_shader_values(program, self._shader_values())
            projection = np.asarray(self.projectionMatrix().data(), dtype=np.float32)
            GL.glUniformMatrix4fv(GL.glGetUniformLocation(program, "u_projection"),
                                  1, False, projection)
        finally:
            GL.glUseProgram(0)
        return program

    def paint(self):
        """Paint flat points normally or glossy/metallic points through the PBR shader."""
        self._flush()
        self._flush_density()
        mode = self._scene["mode"]
        if mode == "flat" or self._gpu_failed:
            if not self._vbo_renderer and self._stock_shader is not None:
                self.shader = self._stock_shader
            super().paint()
            return
        try:
            if not self._vbo_renderer:
                self._pbr.values = self._shader_values()
                self.shader = self._pbr
            if self._density_texture is not None:
                GL.glActiveTexture(GL.GL_TEXTURE1)
                GL.glBindTexture(GL.GL_TEXTURE_3D, self._density_texture)
                GL.glActiveTexture(GL.GL_TEXTURE0)
            super().paint()
        except Exception as exc:
            # One warning and a working map beats a paint-loop traceback on an old driver. The next
            # frame uses the CPU texture and Window notices that ray tracing needs its CPU fallback.
            self._gpu_failed = True
            print(f"starplast: GPU sphere shader unavailable ({type(exc).__name__}: {exc})")
            if not self._vbo_renderer and self._stock_shader is not None:
                self.shader = self._stock_shader
            super().paint()


def to_screen(direction, basis) -> tuple:
    """A world direction as the screen sees it: x right, y up, z toward the viewer.

    The sprite is drawn in the plane of the screen, so a light in world coordinates means nothing to
    it until it is expressed this way -- and doing it every frame is what makes the highlight on the
    balls travel with the lights rather than sitting where it was when they were built.
    """
    _, right, up, forward = basis
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-9:
        return DEFAULT_LIGHT
    d = d / n
    return (float(d @ np.asarray(right, dtype=float)),
            float(d @ np.asarray(up, dtype=float)),
            float(-(d @ np.asarray(forward, dtype=float))))
