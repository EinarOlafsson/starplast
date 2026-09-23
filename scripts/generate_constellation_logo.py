#!/usr/bin/env python3
"""Draw the Toxoplasma constellation refinement as editable vector assets.

Run with QT_QPA_PLATFORM=offscreen. Requires PyQt6 and Noto Sans Light.
The asymmetric outline follows a tachyzoite's tapered anterior and rounded
posterior, with a deliberately softened apex for the logo. See CDC DPDx:
https://www.cdc.gov/dpdx/toxoplasmosis/index.html
The stars are a symbolic gene network, not anatomical structures.
"""
from pathlib import Path
import zipfile

from PyQt6 import QtCore, QtGui, QtSvg

from generate_logo_options import outline, svg

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/assets/toxoplasma-constellation"

# Each end is a continuous curve: no corner or needle-shaped tip.
BODY = ("M177 29C184 28 190 29 187 37"
        "C178 63 161 89 153 117C145 146 144 174 155 198"
        "C165 219 181 228 185 244C190 265 168 281 145 280"
        "C102 279 73 253 60 220C40 172 50 112 79 77"
        "C105 44 150 27 177 29Z")
STARS = [(151, 65, 7), (117, 103, 8), (89, 147, 9),
         (102, 190, 12), (81, 221, 6), (135, 247, 8), (139, 221, 6)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 5)]


def star_path(x: float, y: float, radius: float) -> str:
    """A four-point star with long rays and a softly curved centre."""
    inset = radius * .23
    return (f"M{x} {y-radius}Q{x+inset} {y-inset} {x+radius} {y}"
            f"Q{x+inset} {y+inset} {x} {y+radius}"
            f"Q{x-inset} {y+inset} {x-radius} {y}"
            f"Q{x-inset} {y-inset} {x} {y-radius}Z")


def icon(inverse: bool = False) -> str:
    """A curved tachyzoite containing a seven-star constellation."""
    ink, paper = ("white", "black") if inverse else ("black", "white")
    links = ''.join(f'<path d="M{STARS[a][0]} {STARS[a][1]}L{STARS[b][0]} {STARS[b][1]}"/>'
                    for a, b in EDGES)
    stars = ''.join(f'<path d="{star_path(x, y, radius)}"/>' for x, y, radius in STARS)
    return (f'<g fill="none" stroke="{ink}" stroke-linecap="round" stroke-linejoin="round">'
            f'<path d="{BODY}" stroke-width="1.8"/>'
            # An apical cap and a short inner-membrane line give the contour biological detail.
            '<g stroke-width="1.05"><path d="M155 42Q165 44 179 51"/>'
            '<path d="M149 47Q160 49 176 57"/>'
            '<path d="M117 49C84 69 62 106 59 143"/>'
            '<ellipse cx="102" cy="190" rx="19" ry="24" transform="rotate(-15 102 190)"/></g>'
            f'<g stroke-width="1.1">{links}</g>'
            f'<g stroke-width="1.4" fill="{paper}">{stars}</g></g>')


def render(source: str, path: Path, width: int, height: int, dark: bool = False) -> None:
    """Rasterize vector artwork for a portable preview without changing the SVG."""
    canvas = QtGui.QImage(width, height, QtGui.QImage.Format.Format_ARGB32)
    canvas.fill(QtCore.Qt.GlobalColor.black if dark else QtCore.Qt.GlobalColor.white)
    painter = QtGui.QPainter(canvas)
    QtSvg.QSvgRenderer(source.encode()).render(painter)
    painter.end()
    canvas.save(str(path))


def main() -> None:
    """Write both colourways, wordmarks, previews, and a standalone comparison page."""
    app = QtGui.QGuiApplication([])
    OUT.mkdir(parents=True, exist_ok=True)
    word = outline("starplast", 86)
    for inverse in (False, True):
        suffix = "-white" if inverse else ""
        ink = "white" if inverse else "black"
        mark = svg(icon(inverse), 240, 310, "Starplast — Toxoplasma constellation")
        lockup = svg(icon(inverse) + f'<g color="{ink}" transform="translate(253 188)">{word}</g>',
                     650, 310, "Starplast — Toxoplasma constellation wordmark")
        (OUT / f"icon{suffix}.svg").write_text(mark)
        (OUT / f"wordmark{suffix}.svg").write_text(lockup)
        render(mark, OUT / f"icon{suffix}.png", 720, 930, inverse)
        render(lockup, OUT / f"wordmark{suffix}.png", 1300, 620, inverse)
    preview = svg('<rect width="1200" height="700" fill="white"/>'
                  '<g transform="translate(58 30) scale(1.8)">' + icon() + '</g>'
                  '<g color="black" transform="translate(475 341)">' + word + '</g>',
                  1200, 700, "Starplast — refined Toxoplasma constellation logo")
    (OUT / "preview.svg").write_text(preview)
    render(preview, OUT / "preview.png", 1200, 700)
    (OUT / "index.html").write_text('''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Starplast · Toxoplasma constellation</title>
<style>*{box-sizing:border-box}body{margin:0;background:white;color:black;font:16px/1.6 system-ui,sans-serif}
main{max-width:1100px;margin:auto;padding:40px 24px}h1{font-size:30px;font-weight:300}p{max-width:750px}
a{color:inherit;text-underline-offset:4px}nav{display:flex;gap:24px;flex-wrap:wrap}.preview img{width:100%;height:auto}
.inverse{display:none}label{display:block;margin-top:24px;cursor:pointer}body:has(#dark:checked){background:black;color:white}
body:has(#dark:checked) .normal{display:none}body:has(#dark:checked) .inverse{display:block}</style>
<main><h1>Toxoplasma constellation</h1><p>An elongated, asymmetric parasite with a softly tapered apex, a rounded rear,
an apical cap, and a nucleus ring. Seven stars form a connected constellation entirely inside its body.</p>
<nav><a href="icon.svg" download>Icon SVG</a><a href="wordmark.svg" download>Wordmark SVG</a>
<a href="starplast-constellation.zip" download>Download both colourways</a><a href="../logo-options/index.html">Earlier options</a></nav>
<label><input id="dark" type="checkbox"> White on black</label><div class="preview">
<img class="normal" src="wordmark.svg" alt="Thin black Toxoplasma outline enclosing a seven-star constellation beside starplast">
<img class="inverse" src="wordmark-white.svg" alt="White Toxoplasma constellation wordmark on black"></div>
<p>The silhouette is informed by <a href="https://www.cdc.gov/dpdx/toxoplasmosis/index.html">CDC DPDx's description of tachyzoite morphology</a>.
The constellation represents a gene network.</p></main></html>\n''')
    with zipfile.ZipFile(OUT / "starplast-constellation.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.glob("*.svg")):
            archive.write(path, path.name)
    print(f"Wrote refinement to {OUT}")


if __name__ == "__main__":
    main()
