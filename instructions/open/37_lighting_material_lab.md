# Continuous flashlight, material lab, and ray-rendering comparison — open

## State

The GPU sphere renderer now makes glossy points read as three-dimensional, but the mouse light is
still anchored to the nearest projected gene. Crossing between genes on different depth planes
therefore teleports the light through the map and changes the complete illumination. Metallic also
contains a procedural strip reflection that appears as a line through every point, with a bright dot
following the mouse. Ray-traced transport is materially better but needs a broader comparison set.

## What to do

1. Make mouse lighting a true flashlight: its origin belongs to the camera and its continuous
   screen-space direction follows the cursor. It must never use nearest-gene identity or depth.
2. Keep **selected gene** and **selected gene and its edges** as point emitters whose light originates
   at the clicked gene coordinates.
3. Replace the metallic strip artifact with area-like environment reflections and tune metal so it
   remains reflective, three-dimensional, color-preserving, and readable at the full-map distance.
4. Add a deliberately broad temporary comparison matrix for light transport, pointer beam shape,
   target marker, mood, response, location, and point material. These are candidates for the user to
   prune after direct comparison, not a claim that every option should survive.
5. Improve ray-traced density shadows without reintroducing flicker. Static scenes must remain
   pixel-identical and all modes must retain an honest name: volumetric density ray marching is not
   Vulkan path tracing or hardware triangle ray tracing.
6. Render every meaningful cross-product headlessly through the production pyqtgraph 0.14
   environment, record image differences, frame times, failures, and static-frame hashes.

## Done when

- Moving the cursor over overlapping depth planes changes the flashlight continuously and never
  queries the nearest gene.
- Clicked modes emit from the selected coordinates (and visible-edge neighbours where requested).
- No metallic preset contains a line/strip reflection artifact.
- The comparison controls are persistent, documented, visibly distinct, and covered by tests.
- Every tested production render compiles with no paint-loop traceback; static ray frames match.
- The full project regression passes.
