# Tooltips, docstrings, and the scoring explainer in the app — DONE 2026-08-12

Requested: the precision/recall explanation inside the application, tooltips on every option, and
docstrings for the eventual public API.

## Done

**The explainer.** `objectives.EXPLANATION` carries the definitions, the mapping from the four ways
of wanting "good structure" onto precision / recall / F1, and the measured table showing that each
has a degenerate solution that wins it outright. Shown under Help; a test asserts the dialog displays
that exact string rather than a copy, because two copies drift and the one on screen is the one
people act on.

**Tooltips.** An audit found 85 controls with none, almost all in the analysis panel — exactly where
choosing wrongly is expensive. All filled, and a test walks every combo, spin box, checkbox and
button in the window and fails on any without one. They say why the control exists and what choosing
badly costs. Feature-block tooltips generate from the registry so they cannot go stale.

**Docstrings.** An audit found 28 public modules, classes and functions with none — including the
main `Window`, the `Dataset` record and `cluster()` — plus 57 public methods. All 85 written. Three
tests hold the line per module. Qt event handlers and property getters are exempt.

The bar is mechanical: a docstring exists and says something. Whether it says *why* rather than
*what* stays a review judgement.
