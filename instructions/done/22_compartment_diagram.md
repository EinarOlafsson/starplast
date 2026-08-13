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
