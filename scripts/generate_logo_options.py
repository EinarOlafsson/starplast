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
# All ten marks combine a parasite silhouette, clustered points, and explicit links.
CRESCENT = '<path d="M78 12C43 8 12 36 13 65C14 90 38 99 62 83C38 83 37 62 50 42C59 28 72 21 78 12Z"/>'
SLENDER = '<path d="M71 9C39 17 16 47 21 76C25 97 46 95 61 81C40 79 38 62 47 43C55 28 64 19 71 9Z"/>'


def dots(points, radius=1.0):
    """Draw unconnected observations without turning every dot into a network node."""
    return '<g fill="currentColor" stroke="none">' + ''.join(
        f'<circle cx="{x}" cy="{y}" r="{radius}"/>' for x, y in points) + '</g>'


def cloud(x, y, scale=1.0, rotation=0):
    """A small, irregular UMAP-like cluster; schematic, not a biological dataset."""
    points = [(-8, 0), (-5, -5), (-4, 3), (-1, -2), (0, 6),
              (3, -6), (4, 1), (7, -2), (8, 4), (1, 1)]
    return f'<g transform="translate({x} {y}) rotate({rotation}) scale({scale})">' + dots(points) + '</g>'


def network(points, edges, radius=2.5):
    """Draw sparse connections with outlined nodes overlaid on their endpoints."""
    lines = ''.join(f'<path d="M{points[a][0]} {points[a][1]}L{points[b][0]} {points[b][1]}"/>'
                    for a, b in edges)
    rings = ''.join(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="white"/>' for x, y in points)
    return '<g stroke-width="1.05">' + lines + rings + '</g>'


def star(x, y, radius=5):
    """A four-point star marks a selected gene without adding a filled symbol."""
    return (f'<path d="M{x} {y-radius}Q{x} {y} {x-radius} {y}'
            f'Q{x} {y} {x} {y+radius}Q{x} {y} {x+radius} {y}'
            f'Q{x} {y} {x} {y-radius}Z" fill="white"/>')


OPTIONS = [
    ("crescent-map", "Crescent map", "A parasite crescent opens onto a connected point cloud.",
     CRESCENT + cloud(30, 53, .7, -35) + cloud(49, 28, .65) + cloud(38, 82, .7)
     + cloud(77, 44, .85, -20) + cloud(75, 74, .65)
     + network([(30, 55), (49, 29), (77, 44), (74, 73)], [(0, 1), (0, 2), (2, 3)])
     + star(30, 55, 4)),
    ("apical-atlas", "Apical atlas", "The pointed apex anchors a network across three clusters.",
     SLENDER + '<path d="M58 19L66 25 M55 23L62 29"/>'
     + cloud(30, 66, .7, 65) + cloud(67, 52, .85) + cloud(77, 81, .65)
     + network([(48, 36), (31, 66), (67, 52), (78, 81)], [(0, 1), (0, 2), (2, 3)])),
    ("dissolving-cell", "Cell to constellation", "A curved cell boundary dissolves into mapped genes.",
     '<path d="M76 13C43 8 12 36 13 65C14 87 33 97 53 88 M60 83C38 83 37 62 50 42C59 28 72 21 76 13"/>'
     + dots([(56,87),(60,86),(65,85),(64,90),(70,87),(73,83)])
     + cloud(28, 55, .7, 50) + cloud(69, 45, .8) + cloud(84, 67, .6)
     + network([(28, 55), (50, 33), (69, 45), (84, 66), (64, 84)], [(0, 1), (1, 2), (2, 3), (3, 4)])),
    ("point-cloud", "Point-cloud parasite", "Clustered genes trace the crescent; links reveal its network.",
     '<path d="M76 12C45 10 14 39 15 65 M16 73C23 94 45 94 60 84" stroke-dasharray="2 5"/>'
     + cloud(58, 22, .9, -25) + cloud(41, 35, .9, -40) + cloud(28, 51, 1, -65)
     + cloud(27, 71, .9, 65) + cloud(43, 83, .9, 15)
     + network([(59, 22), (40, 35), (28, 54), (29, 72), (47, 84)], [(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)])),
    ("cluster-bridge", "Cluster bridge", "A parasite connects distinct islands in the gene map.",
     '<g transform="translate(28 14) scale(.65)">' + CRESCENT + '</g>'
     + cloud(18, 33, .85, -20) + cloud(82, 32, .8, 30) + cloud(76, 81, .9)
     + network([(18, 33), (51, 47), (82, 32), (76, 81)], [(0, 1), (1, 2), (1, 3)])
     + star(51, 47, 4)),
    ("orbital-map", "Orbital map", "An apicomplexan sits inside a sparse star-map orbit.",
     '<ellipse cx="50" cy="51" rx="43" ry="28" transform="rotate(-32 50 51)" stroke-width=".8"/>'
     + '<g transform="translate(25 8) scale(.7)">' + SLENDER + '</g>'
     + cloud(21, 63, .65) + cloud(78, 34, .7) + cloud(74, 74, .55)
     + network([(21, 63), (46, 46), (78, 34), (74, 74)], [(0, 1), (1, 2), (1, 3)])),
    ("inner-atlas", "Inner atlas", "A broad parasite outline contains a small UMAP network.",
     '<path d="M81 10C43 6 12 35 13 68C14 90 32 100 55 85C33 69 45 39 81 10Z M66 17L73 24"/>'
     + cloud(48, 27, .7, -30) + cloud(27, 54, .65, 70) + cloud(31, 78, .6)
     + network([(49, 27), (26, 53), (31, 77)], [(0, 1), (1, 2), (0, 2)])),
    ("apical-fan", "Apical fan", "Connections fan from the parasite into separate gene clusters.",
     '<g transform="translate(0 12) scale(.8)">' + SLENDER + '</g>'
     + cloud(79, 24, .7) + cloud(84, 52, .7, 25) + cloud(73, 82, .75, -15)
     + network([(36, 49), (78, 24), (84, 52), (74, 82)], [(0, 1), (0, 2), (0, 3)])
     + star(36, 49, 4)),
    ("membrane-network", "Membrane network", "The parasite contour becomes the backbone of a star map.",
     '<path d="M78 12C43 8 12 36 13 65C14 90 38 99 62 83 M78 12C65 24 54 30 48 44C37 65 40 82 62 83"/>'
     + cloud(27, 51, .65) + cloud(44, 82, .6) + cloud(74, 58, .75, 30)
     + network([(69, 17), (26, 49), (44, 82), (74, 58)], [(0, 1), (1, 2), (1, 3), (2, 3)])),
    ("starplast", "Starplast", "One selected gene ties the crescent to its clustered neighbours.",
     CRESCENT + '<path d="M65 15L70 22"/>'
     + cloud(29, 55, .7, 45) + cloud(44, 82, .55) + cloud(81, 40, .7) + cloud(79, 76, .65)
     + network([(29, 55), (57, 53), (81, 40), (79, 76), (44, 82)], [(0, 1), (1, 2), (1, 3), (1, 4)])
     + star(57, 53, 6)),
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
    return '<g fill="none" stroke="currentColor" stroke-width="1.25" stroke-linecap="round" stroke-linejoin="round">' + body + '</g>'


def main() -> None:
    """Write individual marks, wordmarks, a preview, and a downloadable gallery."""
    app = QtGui.QGuiApplication([])
    OUT.mkdir(parents=True, exist_ok=True)
    # Replace the previous generated proposals; keep exactly ten current options.
    for previous in OUT.glob("[0-9][0-9]-*.svg"):
        previous.unlink()
    word = outline("starplast", 61)
    sheet = ['<rect width="1440" height="1560" fill="white"/>',
             '<g color="black" transform="translate(70 82)">' + outline("starplast / ten directions", 34) + '</g>',
             '<g color="black" transform="translate(70 121)">' + outline("Apicomplexan silhouettes. UMAP clusters. Gene networks.", 18) + '</g>']
    cards = []
    for number, (slug, name, description, geometry) in enumerate(OPTIONS, 1):
        stem = f"{number:02}-{slug}"
        composition = '<g transform="translate(12 10) scale(1.4)">' + mark(geometry) + '</g><g transform="translate(174 100)">' + word + '</g>'
        (OUT / f"{stem}.svg").write_text(svg('<g color="black">' + composition + '</g>', 480, 160, f"Starplast — {name}"))
        inverse = composition.replace('fill="white"', 'fill="black"')
        (OUT / f"{stem}-white.svg").write_text(svg('<g color="white">' + inverse + '</g>', 480, 160, f"Starplast — {name}, white"))
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
    *{box-sizing:border-box}body{margin:0;background:white;color:black;font:16px/1.6 system-ui,sans-serif}main{max-width:1200px;margin:auto;padding:56px 32px}h1{font-size:36px;font-weight:300;margin:0}header p{max-width:680px}a{color:inherit;text-underline-offset:4px}header a{margin-right:24px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:48px;margin-top:48px}article{border-top:1px solid;padding-top:16px}h2{font-weight:400;font-size:15px;margin:0}.preview{padding:25px 0}.preview img{width:100%;height:auto}.light-ink{display:none}nav{display:flex;gap:24px;font-size:14px}article p{font-size:14px}.refinement{margin-top:36px;padding:24px 0;border-block:1px solid}.refinement img{width:100%;max-height:370px}.refinement h2{font-size:23px}.refinement a{display:inline-block;margin-right:24px}label{display:block;margin-top:24px;cursor:pointer}body:has(#inverse:checked){background:black;color:white}body:has(#inverse:checked) .dark-ink{display:none}body:has(#inverse:checked) .light-ink{display:block}@media(max-width:700px){.grid{grid-template-columns:1fr}main{padding:32px 20px}}@media(prefers-reduced-motion:reduce){*{scroll-behavior:auto}}
    </style><main><header><h1>starplast / ten directions</h1><p>Each mark combines an apicomplexan crescent, an irregular UMAP-like point cloud, and a sparse gene network. Black and white, with thin lines. These are schematic symbols, not scientific plots.</p><a href="starplast-logo-options.zip" download>Download all SVGs</a><a href="comparison.png">Comparison sheet</a><label><input id="inverse" type="checkbox"> Preview white on black</label></header><section class="refinement"><h2>Latest refinement / Toxoplasma constellation</h2><a href="../toxoplasma-constellation/index.html"><img class="dark-ink" src="../toxoplasma-constellation/wordmark.svg" alt="Toxoplasma with softened ends enclosing seven connected stars"><img class="light-ink" src="../toxoplasma-constellation/wordmark-white.svg" alt="White Toxoplasma constellation wordmark"></a><p>A more anatomical silhouette with a softly tapered apex, a rounded rear, and an obvious constellation inside.</p><a href="../toxoplasma-constellation/index.html">View the refinement</a><a href="../toxoplasma-constellation/starplast-constellation.zip" download>Download SVGs</a></section><section class="grid">'''
    (OUT / "index.html").write_text(page + ''.join(cards) + '</section></main></html>\n')
    with zipfile.ZipFile(OUT / "starplast-logo-options.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.glob("*.svg")):
            archive.write(path, path.name)
    print(f"Wrote 10 proposals to {OUT}")


if __name__ == "__main__":
    main()
