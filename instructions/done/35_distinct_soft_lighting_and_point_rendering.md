# Simplify lighting and make the remaining render modes genuinely distinct — done

> Historical first pass. Task 36 supersedes its CPU sphere textures and CPU-first shadow transport
> with per-fragment OpenGL PBR spheres and GPU density-ray marching.

## State

The current lighting controls offer too many choices that do not produce useful visual differences.
The only interaction-lighting states worth keeping are:

- selected gene;
- selected gene and its edges;
- mouse proximity.

Mouse lighting is currently too strong. The old ray controls do not have an understandable visual purpose.
The surface-finish choices (including metallic/glossy variants) look effectively identical on the
UMAP points; none produces convincingly three-dimensional, highly glossy points. Hard white lighting
is also not the desired look.

## Why it matters

A rendering preference is useful only when its alternatives are visually and conceptually distinct.
Near-duplicate controls make the interface harder to understand without improving the map. The points
need enough shape, highlight and depth to read as a 3D cloud, while illumination must stay soft enough
that it does not bleach category colors or hide dense structure.

## What to do

1. Keep exactly the three useful interaction scopes: **selected gene**, **selected gene and its
   edges**, and **mouse proximity**. Reduce the strength and falloff harshness of mouse lighting.
2. Remove the old ray/emitter/bounce controls and the current surface-finish controls. A later
   instruction (2026-08-14) explicitly requested **a ray-tracing option**: retain one only if it
   describes a real, visible transport calculation rather than decorative rays or a renamed preset.
3. Replace them with a small, coherent control set:
   - **rendering engine/path**: expose Vulkan or another path only if the implementation truly uses it
     and produces a meaningful capability or quality difference; do not label a shader preset as an
     engine;
   - **light color/mood**: a few visibly distinct soft-light choices such as neutral white, cool blue,
     and warm light;
   - **point render mode**: a few materially distinct modes, including a simple baseline and a real
     3D glossy mode. A metallic mode belongs only if its response is visibly different from glossy.
4. Implement soft illumination with broad highlights and controlled ambient/fill light rather than a
   hard white spot or decorative rays. Preserve the underlying data color well enough for labels and
   continuous color maps to remain comparable.
5. For 3D glossy points, use an actual sphere/impostor/mesh or shader path with surface normals,
   specular response and depth—not renamed flat sprites. Dense point clouds must remain responsive.
6. Prefer fewer modes. Delete any candidate mode that cannot be distinguished reliably in rendered
   comparisons.

## How to know it worked

- Render the same fixed UMAP, camera and color assignment under every retained mode and inspect the
  output; each choice must be identifiable without reading the control value.
- The 3D glossy mode has visible round form, depth ordering and a soft specular highlight; it must not
  look like the flat baseline.
- Cool, neutral and warm lighting are visibly different but do not replace or wash out the data color.
- Mouse proximity is noticeably softer than the current implementation and does not dominate selected
  gene lighting.
- Selected-gene and selected-gene-plus-edges states remain available and clearly distinct.
- The old rays/emitter/bounce and surface-finish controls are absent. The one retained **ray traced**
  light mode must produce measured occlusion and explain exactly what is traced.
- Check dark and light themes, dense and sparse regions, all three LOD levels, and representative edge
  layers. Record comparison screenshots; tests alone are not sufficient for this rendering task.
- Measure frame time on the full 8,140-point map. The higher-quality mode may cost more, but interaction
  must remain usable and the simple mode must remain a responsive fallback.

## Traps

- Vulkan is a graphics API/backend, not a lighting style. Do not expose it as a cosmetic option unless
  Starplast actually has a Vulkan renderer and can maintain that path.
- pyqtgraph's existing scatter path cannot produce real per-point spherical normals or physically based
  glossy/metallic shading by changing a material label. This may require a custom shader, sphere
  impostors, instanced geometry, or a different renderer.
- Bloom alone makes bright discs glow; it does not make them three-dimensional.
- Hard highlights and additive blending can turn dense regions white and destroy category colors.
- This task explicitly asks for fewer, better choices. Do not satisfy it by adding another layer of
  presets over the controls being removed.

## Outcome — completed 2026-08-14 in v0.32.0

Preferences now asks four independent questions with no decorative controls:

- light render mode: `off`, `soft`, `ray traced`;
- light target: `mouse proximity`, `selected gene`, `selected gene and its edges`;
- light mood: `neutral`, `cool blue`, `warm`;
- point render mode: `flat`, `glossy 3D`, `metallic 3D`.

`ray traced` sends one shadow ray per displayed gene toward the active light through the cached
48-cell point-density volume. Dense clusters attenuate genes behind them. This is genuine volumetric
shadow-ray tracing, but not Vulkan RT or path tracing: Starplast has point sprites rather than
triangles and implements neither reflections nor refraction. The UI tooltip says that explicitly.

The two 3D modes use 64×64 sphere-impostor textures with per-texel normals. Glossy uses a tight white
specular highlight; metallic uses a darker body, colored glint and stronger rim. They draw 30–35%
larger than the flat disc so the sphere form remains visible after texture downsampling. Flat has no
spherical highlight. Mouse local gain was reduced to 0.62 and broadened to 0.46 map radii.

### Verification

- `tests/test_display.py`, `tests/test_sprite.py`, and `tests/test_rays.py`: 114 passed. Tests compare
  median 8-bit differences, point-sprite highlight area, soft mouse contrast, mood color, cached
  occupancy, and front/back transmittance rather than merely checking that settings changed. The
  three rendering modules retain 100% statement coverage (285/285).
- Real OpenGL framebuffer comparisons are in `results/lighting_2026_08_14/`: `point_modes.png`,
  `light_moods.png`, and `light_transport.png`, plus all eight full-resolution source frames.
- `scripts/benchmark_lighting.py` renders the fixed 8,140-gene scene and records the actual GL
  renderer. On conservative headless llvmpipe, soft modes painted in 6.5–7.4 ms median and ray-traced
  modes in 15.3–16.1 ms median (20 frames); NVIDIA should not be inferred from that software result.
- Legacy settings migrate: `lit` → `soft`, `mouse` → `mouse proximity`, and the old 2D/matt/satin/
  glossy/metallic finishes map onto the three retained point modes.
- Complete project regression: 2,315 passed, 5 skipped in 158.56 seconds.
