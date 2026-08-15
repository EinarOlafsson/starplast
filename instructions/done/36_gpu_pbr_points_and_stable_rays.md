# Replace matte sprites with GPU PBR spheres and stabilize ray tracing — done

## State

The first simplified lighting pass made point modes measurably different, but the running map still
reads as flat matte discs. The CPU density-ray mode improves depth but flickers when the mouse moves
between overlapping genes, and its frame cost scales on the CPU.

## What to do

1. Render glossy and metallic genes through the existing OpenGL context rather than another renamed
   CPU texture. Reconstruct sphere normals per fragment, provide genuinely different dielectric and
   metallic responses, and write curved per-fragment depth.
2. Preserve categorical and continuous data colors; reflection may shape a color but must not turn
   every category into white or gray.
3. Move density shadow-ray marching to the GPU where supported. Keep an honest CPU fallback and do
   not label density-volume rays as Vulkan, path tracing, or hardware triangle ray tracing.
4. Remove temporal instability. A static map/light must produce identical frames; mouse targeting
   should ease between discrete gene anchors rather than teleporting the complete shadow field.
5. Judge fixed-camera OpenGL renders directly and time the full 8,140-gene map.

## Done when

- Glossy and metallic points visibly read as curved spheres at the normal full-map camera distance.
- Metallic has colored environment reflection and a Fresnel rim; glossy has a dielectric GGX glint.
- Sphere fragments write their curved depth rather than every point behaving as one flat depth card.
- Static ray-traced frames are pixel-identical across repeated paints.
- GPU ray mode is faster than the CPU fallback on the same headless renderer.
- The full suite and 100% rendering-module coverage remain green.

## Outcome — completed 2026-08-14 in v0.33.0

`sprite.ShadedScatter` now installs a GLSL vertex/fragment program into pyqtgraph's existing OpenGL
context. It supports both the fixed-function scatter API in pyqtgraph 0.13 and the explicit VBO
attributes and matrices introduced in 0.14. Each point fragment reconstructs a hemisphere normal from
`gl_PointCoord`, evaluates a GGX microfacet BRDF, reflects a procedural soft studio environment, and
writes the curved surface to `gl_FragDepth`. Glossy uses a dielectric F0 and white glint; metallic
uses the data color as F0, a darker body, environment strip, and Fresnel rim. The original CPU
textures remain only as a driver fallback.

The 48³ occupancy grid is uploaded once as an `R32F` 3D texture. In glossy/metallic ray mode, the
vertex shader marches 24 density samples from each gene toward up to eight active lights and passes
one visibility value to the sphere's fragments. Flat mode and rejected shader contexts use the
existing NumPy fallback. This is GPU volumetric ray marching, not Vulkan, BVH/path tracing, or use of
hardware ray cores.

Mouse light positions use a 0.28 low-pass step so crossing between overlapping projected genes moves
the ray origin instead of teleporting the entire shadow field. A fixed selected-gene scene produced
one unique framebuffer SHA-256 across 12 repeated paints.

## Verification

- Fixed OpenGL comparisons: `results/pbr_lighting_2026_08_14/point_modes.png`,
  `light_moods.png`, and `light_transport.png`; full-resolution source frames are alongside them.
- On the complete 8,140-gene map under conservative llvmpipe software OpenGL, glossy/metallic soft
  rendering took 6.9–7.0 ms median. GPU ray marching took 7.7–7.8 ms versus 15.3 ms for the flat-mode
  CPU fallback over 30 frames.
- Each retained point mode changed 7.8–8.4% of framebuffer pixels by more than 5/255; ray transport
  changed 7.9%. Cool and warm moods changed 7.9% and 8.1%, respectively.
- `ray_stability.csv`: 12 frames, one unique hash, stable=true.
- The installed production environment (pyqtgraph 0.14.0, OpenGL 4.5 compatibility profile) also
  compiles the VBO shader and renders with `gpu_failed=False`; 12 soft frames and 12 ray-traced
  frames were each pixel-identical under llvmpipe. The flat fallback no longer looks up the private
  `pointSprite` registry name removed in pyqtgraph 0.14.
- Focused renderer suite: 121 passed; `lighting.py`, `rays.py`, and `sprite.py` retain 100% statement
  coverage (390/390). Complete project regression: 2,322 passed, 5 skipped in 158.15 seconds.
