# Navigate and Select: two modes for the left mouse button

Independent of everything else and small. Left-click currently rotates and picks at once.

## Navigate
Free rotate, or constrained to x / y / z, plus spin. The constrained axes matter for comparing two
maps: an unconstrained orbit never returns to the same view twice.

## Select
Draw a gate and take the genes inside it.

- **2D** — a lasso in screen space, catching everything behind it. What you want for "grab that
  visual cluster".
- **3D** — a box or brush in world space. What you want when the cloud is deep and a 2D lasso would
  also catch the far side.

## Why it matters beyond convenience

A gated selection is the natural input to annotation: gate a cluster, see its label composition, mark
the unlabelled members. Without it the only way to select a group is one gene at a time.

## Done when

The mode is switchable from the toolbar and the right-click menu, both gate shapes return the correct
gene set, and a gated set can be exported and used as the input to the annotation flow.
