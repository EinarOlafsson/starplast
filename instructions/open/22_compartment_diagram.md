# A cell diagram under the compartment list, coloured to match the map

Requested 2026-08-12.

## What to build

When the filter category is `compartment` (or `compartment_best`), show an apicomplexan cell diagram
beneath the class list, with each organelle filled in **the same colour that compartment has in the
3D map**. So the legend stops being a list of names and becomes a picture of a parasite.

Source file: `starplast/data/icons/Apicomplexa_cells.svg`.

**It is not in the repository yet.** It was placed on the work machine at
`/media/carruthers/mnt3/claude/repo/starplast/starplast/data/icons/Apicomplexa_cells.svg` and needs
committing before any of this can be built or tested. First step is `git add` it.

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
