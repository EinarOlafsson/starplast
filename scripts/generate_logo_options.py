#!/usr/bin/env python3
"""Draw ten monochrome SVG proposals and a comparison sheet.

Run with QT_QPA_PLATFORM=offscreen. Requires PyQt6 and Noto Sans Light;
lettering is converted to paths so the delivered SVGs need no installed font.
These are proposals, separate from the application's selected logo.
"""
from pathlib import Path
import html
import zipfile

from PyQt6 import QtCore, QtGui, QtSvg

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/assets/logo-options"
OPTIONS = [
    ("crescent", "Crescent", "The curved silhouette of an apicomplexan parasite.",
     '<path d="M79 16C25 15 9 72 29 88C57 103 90 59 79 16Z M79 16C61 28 46 56 29 88"/>'),
    ("apical-star", "Apical star", "A parasite crescent ending in a small four-point star.",
     '<path d="M69 24C27 18 12 64 29 84C51 99 82 68 79 42 M69 24C49 35 43 61 29 84 M80 10Q80 27 65 27Q80 27 80 44Q80 27 95 27Q80 27 80 10Z"/>'),
    ("gene-map", "Gene map", "One gene placed in the context of its neighbours.",
     '<path d="M23 32L50 49L79 25 M50 49L76 77 M50 49L24 80"/><circle cx="23" cy="32" r="4"/><circle cx="79" cy="25" r="4"/><circle cx="76" cy="77" r="4"/><circle cx="24" cy="80" r="4"/><circle cx="50" cy="49" r="8" fill="white"/>'),
    ("plastid", "Plastid", "Nested membranes around a single point of evidence.",
     '<ellipse cx="50" cy="50" rx="28" ry="41" transform="rotate(35 50 50)"/><ellipse cx="50" cy="50" rx="17" ry="29" transform="rotate(35 50 50)"/><circle cx="50" cy="50" r="3" fill="currentColor"/>'),
    ("convergence", "Convergence", "Independent observations meeting at one gene.",
     '<path d="M13 25H30C43 25 43 50 55 50H86 M13 50H86 M13 75H30C43 75 43 50 55 50"/><circle cx="13" cy="25" r="3" fill="white"/><circle cx="13" cy="50" r="3" fill="white"/><circle cx="13" cy="75" r="3" fill="white"/><circle cx="79" cy="50" r="7" fill="white"/>'),
    ("helix", "Helix", "A pared-down gene motif with an open, light outline.",
     '<path d="M29 10C29 38 71 62 71 90 M71 10C71 38 29 62 29 90 M31 22H69 M39 37H61 M39 63H61 M31 78H69"/>'),
    ("pathway", "Pathway S", "The initial drawn as a continuous biological pathway.",
     '<path d="M77 23H42C13 23 13 50 42 50H60C89 50 89 78 60 78H23"/><circle cx="77" cy="23" r="4" fill="white"/><circle cx="23" cy="78" r="4" fill="white"/><circle cx="50" cy="50" r="4" fill="white"/>'),
    ("focus", "Focus", "A chosen gene brought into focus within the map.",
     '<circle cx="44" cy="43" r="28"/><path d="M64 63L87 86 M31 48L44 34L57 47"/><circle cx="31" cy="48" r="3" fill="white"/><circle cx="44" cy="34" r="3" fill="white"/><circle cx="57" cy="47" r="3" fill="white"/>'),
    ("evidence", "Evidence layers", "Separate sources aligned around a common gene.",
     '<path d="M50 13L89 34L50 55L11 34Z M11 50L50 71L89 50 M11 66L50 87L89 66"/><circle cx="50" cy="34" r="3" fill="currentColor"/>'),
    ("star-cell", "Star cell", "A restrained star enclosed by an organic cell contour.",
     '<path d="M78 17C99 43 83 89 49 91C14 94 5 59 20 32C33 9 59 4 78 17Z M51 27Q51 50 30 50Q51 50 51 73Q51 50 72 50Q51 50 51 27Z"/>'),
]


def outline(text: str, size: int) -> str:
    """Return an SVG path for font-independent, lightly weighted lettering."""
    font = QtGui.QFont("Noto Sans")
    font.setWeight(QtGui.QFont.Weight.Light)
    font.setPixelSize(size)
    path = QtGui.QPainterPath()
    path.addText(0, 0, font, text)
    commands = []
    index = 0
    while index < path.elementCount():
        element = path.elementAt(index)
        if element.isMoveTo():
            commands.append(f"M{element.x:.3f},{element.y:.3f}")
        elif element.isLineTo():
            commands.append(f"L{element.x:.3f},{element.y:.3f}")
        elif element.isCurveTo():
            second, third = path.elementAt(index + 1), path.elementAt(index + 2)
            commands.append(f"C{element.x:.3f},{element.y:.3f} {second.x:.3f},{second.y:.3f} {third.x:.3f},{third.y:.3f}")
            index += 2
        index += 1
    return '<path fill="currentColor" stroke="none" d="' + ' '.join(commands) + '"/>'


def svg(body: str, width: int, height: int, title: str) -> str:
    """Wrap artwork in an accessible SVG with a fixed coordinate system."""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img"><title>{html.escape(title)}</title>{body}</svg>\n')


def mark(body: str, inverse: bool = False) -> str:
    """Apply one consistent thin stroke and swap opaque cutouts on dark backgrounds."""
    if inverse:
        body = body.replace('fill="white"', 'fill="black"')
    return '<g fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">' + body + '</g>'


def main() -> None:
    """Write individual marks, wordmarks, a preview, and a downloadable gallery."""
    app = QtGui.QGuiApplication([])
    OUT.mkdir(parents=True, exist_ok=True)
    word = outline("starplast", 61)
    sheet = ['<rect width="1440" height="1560" fill="white"/>',
             '<g color="black" transform="translate(70 82)">' + outline("starplast / ten directions", 34) + '</g>',
             '<g color="black" transform="translate(70 121)">' + outline("Thin lines. Black and white. Gene evidence, in context.", 18) + '</g>']
    cards = []
    for number, (slug, name, description, geometry) in enumerate(OPTIONS, 1):
        stem = f"{number:02}-{slug}"
        composition = '<g transform="translate(12 15)">' + mark(geometry) + '</g><g transform="translate(144 86)">' + word + '</g>'
        (OUT / f"{stem}.svg").write_text(svg('<g color="black">' + composition + '</g>', 420, 130, f"Starplast — {name}"))
        inverse = composition.replace('fill="white"', 'fill="black"')
        (OUT / f"{stem}-white.svg").write_text(svg('<g color="white">' + inverse + '</g>', 420, 130, f"Starplast — {name}, white"))
        (OUT / f"{stem}-icon.svg").write_text(svg('<g color="black">' + mark(geometry) + '</g>', 100, 100, f"Starplast — {name} icon"))
        x, y = 70 + ((number - 1) % 2) * 700, 193 + ((number - 1) // 2) * 266
        sheet.append(f'<g color="black" transform="translate({x} {y})">' + outline(f"{number:02} / {name}", 17) + f'<g transform="translate(12 24)">{composition}</g><g transform="translate(0 197)">' + outline(description, 14) + '</g></g>')
        cards.append(f'<article><h2>{number:02} / {name}</h2><div class="preview"><img class="dark-ink" src="{stem}.svg" alt="Starplast {name} proposal"><img class="light-ink" src="{stem}-white.svg" alt="Starplast {name} in white"></div><p>{description}</p><nav><a href="{stem}.svg" download>Black SVG</a><a href="{stem}-white.svg" download>White SVG</a><a href="{stem}-icon.svg" download>Icon</a></nav></article>')
    sheet_svg = svg(''.join(sheet), 1440, 1560, "Ten monochrome Starplast logo proposals")
    (OUT / "comparison.svg").write_text(sheet_svg)
    image = QtGui.QImage(1440, 1560, QtGui.QImage.Format.Format_ARGB32)
    image.fill(QtCore.Qt.GlobalColor.white)
    painter = QtGui.QPainter(image)
    QtSvg.QSvgRenderer(sheet_svg.encode()).render(painter)
    painter.end()
    image.save(str(OUT / "comparison.png"))
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Starplast — logo options</title><style>
    *{box-sizing:border-box}body{margin:0;background:white;color:black;font:16px/1.6 system-ui,sans-serif}main{max-width:1200px;margin:auto;padding:56px 32px}h1{font-size:36px;font-weight:300;margin:0}header p{max-width:680px}a{color:inherit;text-underline-offset:4px}header a{margin-right:24px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:48px;margin-top:48px}article{border-top:1px solid;padding-top:16px}h2{font-weight:400;font-size:15px;margin:0}.preview{padding:25px 0}.preview img{width:100%;height:auto}.light-ink{display:none}nav{display:flex;gap:24px;font-size:14px}article p{font-size:14px}label{display:block;margin-top:24px;cursor:pointer}body:has(#inverse:checked){background:black;color:white}body:has(#inverse:checked) .dark-ink{display:none}body:has(#inverse:checked) .light-ink{display:block}@media(max-width:700px){.grid{grid-template-columns:1fr}main{padding:32px 20px}}@media(prefers-reduced-motion:reduce){*{scroll-behavior:auto}}
    </style><main><header><h1>starplast / ten directions</h1><p>Thin lines, black and white. Ten different ways to connect the identity to parasites, genes, and evidence. Each includes a vector wordmark and a separate icon.</p><a href="starplast-logo-options.zip" download>Download all SVGs</a><a href="comparison.png">Comparison sheet</a><label><input id="inverse" type="checkbox"> Preview white on black</label></header><section class="grid">'''
    (OUT / "index.html").write_text(page + ''.join(cards) + '</section></main></html>\n')
    with zipfile.ZipFile(OUT / "starplast-logo-options.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.glob("*.svg")):
            archive.write(path, path.name)
    print(f"Wrote 10 proposals to {OUT}")


if __name__ == "__main__":
    main()
