# A cell diagram under the compartment list, coloured to match the map

Requested 2026-08-12.

## What to build

When the filter category is `compartment` (or `compartment_best`), show an apicomplexan cell diagram
beneath the class list, with each organelle filled in **the same colour that compartment has in the
3D map**. So the legend stops being a list of names and becomes a picture of a parasite.

Source file: `starplast/data/icons/Apicomplexa_cells.svg`.

Committed 2026-08-12, 206 KB.

## What the file turns out to be, which decides the approach

It is a UniProt-style subcellular-location diagram: 56 distinct organelles, each a `<g>` with a
**UniProt SL identifier** (`id="SL0233"`) and a `<text class="subcell_description">`. So the artwork
already carries a standard controlled vocabulary, and the job is a mapping from hyperLOPIT
compartment names onto SL codes -- not renaming anything in the drawing.

Confirmed present, which covers most of what this proteome measures:

    SL0233  rhoptry            -> rhoptries 1, rhoptries 2
    SL0163  microneme          -> micronemes
    SL0018  apicoplast         -> apicoplast
    SL0362  inner membrane cx  -> IMC
    SL0173  mitochondrion      -> mitochondrion - soluble
    SL0171  mito membrane      -> mitochondrion - membranes
    SL0132  Golgi              -> Golgi
    SL0095  ER                 -> ER, ER 2
    SL0191  nucleus            -> nucleus - chromatin, nucleus - non-chromatin
    SL0188  nucleolus          -> nucleolus
    SL0091  cytosol            -> cytosol
    SL0039  cell membrane      -> PM - integral, PM - peripheral 1, PM - peripheral 2
    SL0090  cytoskeleton       -> tubulin cytoskeleton

Dense granules are in the drawing, which matters: that is the compartment the annotation workflow is
aimed at. Ribosomes and vacuoles are present too.

**Not in the drawing:** conoid, and the proteasome. `19S proteasome`, `20S proteasome` and
`apical 1` / `apical 2` therefore have no organelle to colour and need the same treatment as
cytosol -- named beside the diagram rather than silently absent.

Several hyperLOPIT classes map to ONE organelle (both rhoptry classes, both nucleus classes, the
three PM classes). Those cannot each take their own colour on the same shape.

**DECIDED 2026-08-12: colour the shape for whichever class is selected.** So clicking `rhoptries 2`
fills the single rhoptry shape with the colour `rhoptries 2` has in the map, and clicking
`rhoptries 1` refills the same shape in its colour. The artwork is not split.

What follows from that, and must be built rather than assumed:

- With nothing selected, a shared shape has no single right colour. Draw it neutral -- not the
  colour of whichever class happens to sort first, which would be a claim nobody made.
- The shape must say which class it is currently showing, or a user returning to the window cannot
  tell whether the rhoptry is coloured for 1 or for 2. A label, or the selection visible in the list
  beside it.
- Clicking the shape selects a class, but a shared shape maps to several. Cycle through them on
  repeated clicks, and say which one is now selected.

## Behaviour

- Shown only when the active category is a localisation one; hidden otherwise, because there is no
  sensible mapping from cell-cycle phase onto organelles.
- Organelle fills come from `Window.colour_of`, which is the same dict the points are drawn from, so
  the diagram and the map cannot disagree. Never a second palette.
- Selecting a compartment in the list highlights that organelle; clicking an organelle selects the
  compartment. Both directions, or the diagram is decoration.
- Compartments with no organelle in the drawing (`cytosol`, the two `PM` classes, `unassigned`) need
  a decision rather than silence: label them beside the diagram, or shade the cell body. They must
  not simply vanish, since between them they are most of the proteome.
- **Absence stays grey.** `unassigned` is 4,313 genes and must not be filled with a compartment
  colour or left looking like a measurement.

## How to colour an SVG at runtime

The organelle paths need stable ids or classes matching compartment names. Parse the SVG, rewrite
the `fill` on the matching elements, and render with `QSvgRenderer` from the modified bytes. Do not
ship one recoloured copy per theme -- there are four themes and the categorical palette changes with
each, so it has to be done at draw time from `colour_of`.

If the file's element ids do not match the compartment names, add a mapping table in one place
(`localisation.py` or beside the icon) rather than renaming the artwork, and make a missing mapping
a visible gap rather than a silent no-op.

## Done when

Choosing `compartment` shows the diagram; every organelle that has a colour in the map has the same
colour in the diagram; clicking either side selects on the other; and switching theme recolours
both together.
