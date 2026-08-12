# Navigate and Select: two modes for the left mouse button — DONE 2026-08-12

## Navigate

Free orbit as before, plus a constraint to x, y or z from **View ▸ Left mouse button**. The camera is
spherical, so each axis is one of its two angles held still: about z is azimuth alone; about x or y
is elevation with the azimuth pinned to that axis. That is the point of it — a free orbit never
returns to the same view twice, and two maps compared from "x" are compared from the same place.

## Select

Drag to draw a gate. Two shapes, and they answer different questions:

- **lasso (2D)** — a polygon in screen space, so it takes everything behind it as well. That is what
  "grab that visual cluster" means: the cluster you can see is a screen-space object.
- **brush (3D)** — a ball in world space around the gene under the press, sized by the drag. Deep in
  the cloud a lasso also catches the far side, which looks like a selection of one structure and is a
  selection of two. The pixels-to-world scale is **measured from the projection itself** rather than
  computed from the camera parameters, for the same reason `_projection_matrix` inspects rather than
  assumes: pyqtgraph has changed that convention twice.

The gate is drawn on a transparent child widget over the GL view, not by painting inside `paintGL` —
mixing QPainter into a QOpenGLWidget's GL painting is driver-dependent, and this project already has
one class of bug that appears only on someone else's machine.

## What a gate produces

`Window.gated` holds the index array. The evidence panel shows the set's **composition by the current
category, leading with how many carry no value for it** — those are the candidates a gate is drawn to
find. The map marks the set by receding everything else and enlarging the gated points; it does
**not** recolour them, because colour here means measurement, inference or absence, and a selection
is none of those.

`File ▸ Export gated selection (CSV)` writes it with coordinates, and it is in the right-click menu
alongside "Clear the gate". That is the annotation input the task asked for: task 17 can take this
set directly.

## Verified

- The gate is a display claim -- "these genes are the ones you drew around" -- so it is tested
  against the projection rather than trusted: every gated gene inside the lasso, every gene inside it
  gated, and for the brush, every gated gene nearer the anchor in WORLD space than any ungated one.
- Genes with no position in the displayed embedding cannot be gated, the same rule picking follows.
- Rendered under Xvfb and looked at: the lasso draws over the live map, and finishing it leaves the
  gated set marked.
- Rotation constraints checked by driving the handler: about z moves the azimuth and leaves the
  elevation exactly where it was.
