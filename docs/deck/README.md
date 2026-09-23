# Starplast slide deck

[Open the viewer](https://einarolafsson.github.io/starplast/deck/) ·
[Start on GitHub](pages/01.md) · [PDF](starplast_deck.pdf) ·
[Editable PowerPoint](starplast_deck.pptx) · [Transcript](transcript.md)

Forty slides introduce Starplast, its evidence, the desktop controls, guided
workflows, protein features, validation results and reproducible exports. This is
a guide to the 0.43 release. Screenshots are captured from the application; the
screen-comparison values are explicitly illustrative. Benchmark charts are drawn
from the committed aggregate results. Model hypotheses remain distinct from
measured annotations.

The dark layout, Open Sans typography, restrained teal accents and viewer
interaction follow the [spaCR deck](https://einarolafsson.github.io/spacr/_static/deck/).
The approved Starplast logo remains unchanged. Open Sans is included under the
SIL Open Font License in `fonts/OFL.txt`.

## Rebuild

The source is [scripts/deck_content.py](../../scripts/deck_content.py). Text, shapes,
charts and diagrams remain editable in the PowerPoint. The PDF uses native text
and vectors; only application screenshots and the logo are raster images.

```bash
pip install reportlab python-pptx Pillow
# Optional: refresh screenshots under Linux with the app dependencies installed.
QT_QPA_PLATFORM=xcb xvfb-run -a -s '-screen 0 1920x1200x24' python scripts/capture_deck.py
python scripts/build_deck.py
python scripts/build_docs.py
```

The builder also needs Poppler's `pdftoppm`. It checks slide count and text bounds,
then emits PDF, PowerPoint, full-size JPEGs, thumbnails, a transcript, numbered
GitHub pages and an offline-capable viewer. All published forms use the same
content source. The guide version is intentionally labelled; review the content
when updating it for a later release.
