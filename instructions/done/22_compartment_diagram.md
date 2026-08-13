# A cell diagram under the compartment list, coloured to match — DONE 2026-08-12

`starplast/celldiagram.py`, under the class list, shown only when the colour source is
`compartment` or `compartment_best` — there is no sensible mapping from cell-cycle phase onto
organelles, and colouring them by one would be a picture of a relationship that does not exist.

## The mapping

The artwork is UniProt's subcellular-location diagram: 59 organelles, each a `<g>` with an SL
identifier, so the drawing already speaks a controlled vocabulary and the work was a table from
hyperLOPIT names onto SL codes. **The table was checked against the file, and two entries in the task
notes were wrong**: the granule shape is `SL0281` (UniProt's generic "Cytoplasmic granule"), not
`SL0246`, and there is no ribosome in the drawing at all. A test asserts every mapped code exists in
the artwork, because a mapping to a shape that is not there is a silent no-op.

Deliberately unmapped, with the reason recorded beside the table: the ribosome and proteasome classes
and `apical 1/2` (no such shape), and `endomembrane vesicles` — the drawing's vesicles are
specifically COPI and COPII, and colouring one for a mixed hyperLOPIT class would be a claim the data
does not make. All of them are **named beside the diagram** rather than dropped.

## The decided behaviour, built

- A shape shared by several classes (both rhoptry classes, both nucleus classes, the three PM
  classes) takes the **selected** class's colour and says which one it is showing. Neutral when none
  is selected: filling it with whichever sorts first would be a claim nobody made.
- Clicking a shared shape **cycles** through the classes it stands for. One shape stands for three PM
  classes; without a cycle two of them would be unreachable from the drawing.
- Clicking selects in the list; selecting in the list recolours the shape. Both directions, or the
  diagram is decoration.
- Fills come from `Window.colour_of` — the same dict the points are drawn from — at draw time, so a
  theme change moves both together and there is never a second palette.
- `unassigned` is never given a compartment colour.

## Two things found by building it

**The fill has to reach the paths, not just the group.** The artwork sets `fill` on each path, so a
fill on the enclosing `<g>` is overridden by every one of them: the diagram comes back grey while
every assertion about the returned string passes.

**Hit-testing by "smallest bounding box" picks the wrong organelle.** Each group's box includes its
text label, so a labelled organelle's box is much wider than its drawing and an unrelated one can be
tighter over the click. Clicking the middle of the nucleolus selected the mitochondrion. The rule is
now the box the click sits most centrally in, relative to that box's own size.

## Verified

- `celldiagram.py` at 100%.
- Rendered under Xvfb and looked at: the parasite draws with rhoptries in the rhoptry colour,
  micronemes in theirs, and the note underneath naming the shared shape.

## Reworked 2026-08-12, after seeing it in use

Three corrections, all of them from the same complaint: the diagram was showing everything at once
and clicking it selected the wrong thing.

**Only the selected compartment is coloured.** Every mapped organelle used to take its colour, which
makes the diagram a second legend — twenty-odd colours to read against twenty-odd names — when the
question a person has is "where is the thing I just clicked". The rest of the drawing is grey.

**The drawing is grey, and transparent where the cell is not.** This is the part that resisted:
every fill, stroke and gradient stop in the file can be set to grey and it still renders in colour,
with no error and nothing in the document to explain it — 145 elements carry a `coloured` class and
eleven gradients cross-reference each other. Chasing that further was archaeology, so the render is
desaturated instead: draw, take the luminance, put the alpha back, then tint the one selected
organelle. A test measures the saturation of the painted widget, because that is the claim.

**Clicks land on the right organelle.** Hit-testing used `boundsOnElement`, and those boxes are
useless here: every group contains hidden `<text>` with UniProt's description of the compartment, so
the Golgi's box is 4,894 units wide in a 1,190-wide drawing. Clicking a rhoptry landed on whatever
inflated box happened to win. Each organelle is now rendered alone into a mask and the click is a
pixel test, with the smallest shape winning where they nest.

Two bugs surfaced under that. `<g id="SL...">` groups were extracted with a non-greedy regex, which
stops at the FIRST nested `</g>` — and these groups nest, so most organelles were being handled as a
fragment: eleven of the fourteen masks were silently empty. And the depth counter that replaced it
treated a self-closing `<g/>` as an opening tag, which left six more groups — the nucleus, the
cytosol, the plasma membrane among them — reported as absent.

**A note for whoever edits this next:** two of the fixes above were written twice, because the first
attempt patched text that the American-spelling pass had already changed (`colour` → `color`), so the
edit silently matched nothing while the code went on doing the old thing. After task 25, patch by
reading the file, not by remembering it.
