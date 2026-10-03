# 58 · The blank panel: a tab switched away from comes back as bare background

**Reported twice, from the desktop, against 0.48.0 and still present in 0.49.0.**

> "if you go to strategy then back to evidence it is not possible to see a gene card again, and then
> going back to strategy, strategy is empty."

> "when I go back and forth between tabs like evidence and maps I can't always see the evidence and I
> only see the background."

## What it is

**Reproduced** on the real display (`DISPLAY=:0`, `QT_QPA_PLATFORM=xcb`) on 2026-09-30, on the second
tab switch: Evidence shows a gene card, Strategies is opened, Evidence is opened again, and the whole
right-hand dock area is nothing but the drifting blob field. It stays that way for every further
switch. Screenshots were taken of each step.

## Why

Stacking, not painting. `Window._apply_ambient` builds the blob background as an ordinary child of
the main window covering the whole of it, and lowers it ONCE. Qt's own dock machinery
(`QMainWindowLayout::tabChanged`) **lowers the dock it has just switched away from** to the bottom of
the sibling stack. After one switch that dock is below the background; when it is brought back the
background is painted over it. The panel is laid out, the right size, `isVisible()` true, its widget
repainting — and invisible. Measured z-order, offscreen, with the background on:

| | background | evidence dock | background above the panel |
|---|---|---|---|
| at startup | 0 | 9 | no |
| after one switch | 1 | 0 | **yes** |
| after switching back | 2 | 1 | **yes** |

It needs no compositor and no translucency. It needs only **Background = blobs**, which is a
Preferences setting the suite's isolated settings leave at `none` — which is exactly why the first
headless check, and every test in the suite, reported each tabified dock visible and painting.

## Ruled out

* `glass.py`: the main window is never `dress`ed, so it carries neither `WA_TranslucentBackground`
  nor a frameless rim; `GlassStyle.polish` only ever dresses `QMenu`, `QTipLabel` and
  `QComboBoxPrivateContainer`, never a dock or its child. With `Background = none` the round trip is
  pixel-identical, with and without compositing.
* The dock title and tab restyle, and `apply_theme` / `set_container_opacity` re-applying the
  stylesheet while a dock is hidden.
* The 0.48 card stack in `strategy_panel.py` / `strategy_card.py` (`old.setParent(None)`): the
  Evidence dock holds a plain `QTextBrowser` and goes blank the same way.
* The `install()` hooks of the newer panels: `star_map`, `umap_gallery_panel` and `guided_panel` only
  tabify and raise, and the bug reproduces between Evidence and Strategies alone, which predate them.
* A tab bar out of step with the dock the layout shows: measured through 120 fuzzed operations
  (tab clicks, `raise_()`, `show()`, `hide()`) under Qt's `vnc` platform — the tab bar and the
  visible dock agreed every time, and no state had nothing visible.

## The fix

* `ambient.AmbientWidget.keep_behind()` — puts the field back at the bottom of its parent's stack,
  and reports whether it had to. Called from `_tick`, so any other reordering (a dock floated and put
  back, a panel hidden and shown) is corrected within a frame.
* `Window._apply_ambient` connects `tabifiedDockWidgetActivated` to it, so the switch itself is
  corrected on the spot rather than up to 40 ms later.

## Checked

`tests/test_display.py`:

* `test_the_background_never_climbs_on_top_of_a_panel` — renders the window and compares the PIXELS
  over the Evidence panel before and after a round trip through another tab, with the blob clock
  stopped so the two renders are comparable, and asserts the background is still below the dock in
  the stacking order. It fails on the unfixed tree and passes on the fixed one; `isVisible()` is true
  either way, which is why it is not asserted.
* `test_the_background_puts_itself_back_at_the_bottom` — `keep_behind` as a unit.

Confirmed on the real display after the fix as well: the panel comes back on every switch.
