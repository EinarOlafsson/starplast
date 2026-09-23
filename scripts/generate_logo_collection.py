#!/usr/bin/env python3
"""Draw forty varied monochrome Toxoplasma/network logo proposals as native SVG.

Eight design families explore outline, silhouette, membrane, open contour,
point cloud, motion and seal compositions. The approved application icon is
unchanged. Run with QT_QPA_PLATFORM=offscreen using the project's Qt environment.
"""
from pathlib import Path
import html
import json
import zipfile

from PyQt6 import QtGui
from generate_logo_options import outline, svg
from generate_constellation_logo import BODY, star_path, render
from generate_logo_refinements import BROAD, SLENDER, COMPACT

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/assets/logo-collection-40'
FAMILIES = [
    ('arc', 'Living arc', 'An upright tachyzoite and a clear constellation.', BODY, 0),
    ('glider', 'Glider', 'A long, horizontal parasite in motion.', SLENDER, 68),
    ('vessel', 'Cell vessel', 'A broader, softer cell with room around the stars.', BROAD, -22),
    ('membrane', 'Double membrane', 'Two fine membrane contours and a branching interior.', COMPACT, 25),
    ('open', 'Open contour', 'The cell outline opens where the constellation continues.', BODY, -48),
    ('cloud', 'Star field', 'Sparse point clouds add the feel of a gene embedding.', BROAD, 0),
    ('silhouette', 'Night cell', 'A dark parasite silhouette contains a fine white constellation.', BODY, 42),
    ('seal', 'Celestial seal', 'A fine circular frame joins the cell to a star-map emblem.', COMPACT, -15),
]
PATTERNS = [
    ('path', 'Path', [(148,73,6),(113,106,7),(84,145,8),(102,184,10),(84,215,5),(133,248,7)],
     [(0,1),(1,2),(2,3),(3,4),(4,5)]),
    ('branches', 'Branches', [(150,72,5),(112,105,7),(80,144,6),(117,164,9),(80,193,6),(102,235,7),(149,250,5)],
     [(0,1),(1,2),(2,3),(2,4),(4,5),(5,6),(3,5)]),
    ('islands', 'Islands', [(148,76,6),(111,99,6),(99,126,7),(80,183,7),(105,206,8),(99,232,5),(145,252,6)],
     [(0,1),(1,2),(2,0),(2,3),(3,4),(4,5),(3,5),(5,6)]),
    ('beacon', 'Beacon', [(147,74,4),(106,108,5),(82,143,5),(103,182,14),(76,211,4),(128,247,5),(138,223,4)],
     [(0,1),(1,2),(2,3),(3,4),(3,5),(3,6),(5,6)]),
    ('mesh', 'Mesh', [(147,76,5),(111,106,6),(85,138,6),(120,151,5),(78,176,6),(105,202,8),(91,226,5),(140,250,6)],
     [(0,1),(1,2),(1,3),(2,3),(2,4),(3,5),(4,5),(4,6),(5,6),(6,7),(5,7)]),
]


def mark(family, pattern, inverse=False):
    """One cell/star composition, in a 320-square view box with rounded joins."""
    slug, _, _, body, angle = family
    _, _, stars, edges = pattern
    ink, paper = ('white','black') if inverse else ('black','white')
    inside = paper if slug == 'silhouette' else ink
    starfill = ink if slug == 'silhouette' else paper
    links = ''.join(f'<path d="M{stars[a][0]} {stars[a][1]}L{stars[b][0]} {stars[b][1]}"/>' for a,b in edges)
    nodes = ''.join(f'<path d="{star_path(x,y,r)}"/>' for x,y,r in stars)
    contour = f'<path d="{body}" fill="{ink if slug == "silhouette" else "none"}" stroke-width="1.7"/>'
    if slug == 'open':
        contour = ('<path d="M152 31C120 32 96 52 79 77C50 112 40 172 60 220C73 253 102 279 145 280'
                   'C168 281 190 265 185 244 M180 231C167 216 160 207 155 198C144 174 145 146 153 117'
                   'C161 89 178 63 187 37C190 29 184 28 177 29" stroke-width="1.7"/>')
    details = '<path d="M153 43Q166 44 180 53 M148 49Q162 51 176 60"/>'
    if slug in ('arc','vessel'):
        details += '<path d="M116 50C87 70 64 106 60 138"/>'
    if slug == 'membrane':
        details += ('<path d="M128 52C91 72 63 109 59 150C55 194 73 237 110 254"/>'
                    '<path d="M146 85C128 113 128 163 137 185"/>'
                    '<ellipse cx="103" cy="185" rx="22" ry="27" transform="rotate(-15 103 185)"/>')
    if slug == 'glider':
        details += '<path d="M118 54C86 78 61 119 61 161 M68 212Q81 243 109 255"/>'
    cloud = ''
    if slug == 'cloud':
        for cx,cy in [(119,78),(78,118),(69,157),(79,226),(124,264)]:
            cloud += ''.join(f'<circle cx="{cx+dx}" cy="{cy+dy}" r="{r}" fill="{ink}" stroke="none"/>'
                             for dx,dy,r in [(-5,-3,.8),(0,0,1.1),(4,5,.65),(6,-4,.9),(-3,6,.7)])
    interior = (f'<g stroke="{inside}" stroke-width=".85">{details}</g>'
                f'<g stroke="{inside}" stroke-width="1.0">{links}</g>'
                f'<g stroke="{inside}" fill="{starfill}" stroke-width="1.15">{nodes}</g>'+cloud)
    cell = f'<g transform="translate(40 3) rotate({angle} 120 155)">{contour}{interior}</g>'
    frame = ''
    if slug == 'seal':
        # The cell is a distinct asymmetric organism inside the astronomical frame.
        cell = f'<g transform="translate(24 24) scale(.85)">{cell}</g>'
        frame = ('<circle cx="160" cy="160" r="144" stroke-width=".9"/>'
                 '<path d="M160 9V23 M160 297V311 M9 160H23 M297 160H311" stroke-width=".9"/>')
    return f'<g fill="none" stroke="{ink}" stroke-linecap="round" stroke-linejoin="round">{frame}{cell}</g>'


def main():
    """Write forty icon/wordmark pairs, reverse versions, gallery and contact sheets."""
    app = QtGui.QGuiApplication([])
    OUT.mkdir(parents=True, exist_ok=True)
    word = outline('starplast', 69)
    sheet = ['<rect width="1800" height="3000" fill="white"/>',
             '<g color="black" transform="translate(45 60)">'+outline('starplast / forty directions',32)+'</g>']
    records, cards = [], []
    for f, family in enumerate(FAMILIES):
        for v, pattern in enumerate(PATTERNS):
            number = f*5+v+1
            slug = f'{number:02}-{family[0]}-{pattern[0]}'
            name = family[1] + ' / ' + pattern[1]
            description = family[2] + ' ' + ['A six-star trail.','A branching seven-star network.',
                'Two linked triangular groups.','A larger central star anchors the connections.',
                'An eight-star triangulated network.'][v]
            for inverse in (False, True):
                suffix = '-white' if inverse else ''
                ink = 'white' if inverse else 'black'
                icon = svg(mark(family,pattern,inverse),320,320,name)
                lockup = svg(mark(family,pattern,inverse)+f'<g color="{ink}" transform="translate(340 183)">{word}</g>',750,320,name)
                (OUT/f'{slug}-icon{suffix}.svg').write_text(icon)
                (OUT/f'{slug}-wordmark{suffix}.svg').write_text(lockup)
            render(svg(mark(family,pattern),320,320,name),OUT/f'{slug}.png',640,640)
            x,y=25+355*v,90+360*f
            sheet += [f'<g transform="translate({x+10} {y})">{mark(family,pattern)}</g>',
                      f'<g color="black" transform="translate({x+20} {y+332})">'+outline(f'{number:02} / {family[1]} / {pattern[1]}',14)+'</g>']
            cards.append(f'''<article data-family="{family[0]}"><h2>{number:02} / {html.escape(name)}</h2>
<div class="art"><img class="light" src="{slug}-wordmark.svg" alt="{html.escape(description)}" loading="lazy">
<img class="dark" src="{slug}-wordmark-white.svg" alt="{html.escape(description)}" loading="lazy"></div>
<p>{description}</p><nav><a href="{slug}-icon.svg" download>Icon SVG</a><a href="{slug}-wordmark.svg" download>Wordmark SVG</a>
<a href="{slug}-icon-white.svg" download>White icon</a><a href="{slug}-wordmark-white.svg" download>White wordmark</a></nav></article>''')
            records.append(dict(number=number,name=name,family=family[0],stem=slug,description=description))
    comparison=svg(''.join(sheet),1800,3000,'Forty Starplast logo options')
    (OUT/'comparison.svg').write_text(comparison)
    render(comparison,OUT/'comparison.png',2700,4500)
    (OUT/'manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    choices=''.join(f'<option value="{slug}">{name}</option>' for slug,name,*_ in FAMILIES)
    (OUT/'index.html').write_text('''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Starplast / forty logo directions</title>
<style>*{box-sizing:border-box}body{margin:0;background:white;color:#111;font:16px/1.6 system-ui,sans-serif}main{max-width:1500px;margin:auto;padding:40px 28px}
h1{font-size:40px;line-height:1.15;font-weight:300}h2{font-size:18px;font-weight:450}p{max-width:860px}a{color:inherit;text-underline-offset:4px}
nav,.controls{display:flex;gap:22px;flex-wrap:wrap}.controls{margin:28px 0;align-items:center}select{font:inherit;padding:7px;background:transparent;color:inherit;border:1px solid #888}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px}article{border:1px solid #ddd;padding:24px;border-radius:10px}article[hidden]{display:none}
article nav{font-size:13px;gap:15px}article p{font-size:14px;min-height:45px}.art img{width:100%;height:auto}.dark{display:none}
body.night{background:#080808;color:white}body.night article{border-color:#444}body.night .light{display:none}body.night .dark{display:block}
@media(max-width:780px){.grid{grid-template-columns:1fr}main{padding:24px 16px}h1{font-size:30px}}
</style></head><body><main><h1>Forty directions for Starplast</h1><p>Eight distinct families, each with five constellation layouts.
Every mark combines an asymmetric parasite with stars and gene connections. Thin monochrome lines remain the common thread; the silhouettes,
orientation, cell details and network topology vary. The active app logo is unchanged.</p>
<nav><a href="comparison.png">See all 40 together</a><a href="starplast-40-logos.zip" download>Download all SVGs</a>
<a href="../toxoplasma-constellation/index.html">Current logo</a></nav>
<div class="controls"><label><input type="checkbox" id="theme"> White on black</label><label>Family <select id="family"><option value="all">All forty</option>'''+choices+'''</select></label>
<span id="count" aria-live="polite">40 options</span></div><div class="grid">'''+''.join(cards)+'''</div></main><script>
document.getElementById('theme').addEventListener('change',e=>document.body.classList.toggle('night',e.target.checked));
document.getElementById('family').addEventListener('change',e=>{let n=0;document.querySelectorAll('article').forEach(a=>{a.hidden=e.target.value!=='all'&&a.dataset.family!==e.target.value;if(!a.hidden)n++});document.getElementById('count').textContent=n+' options'});
</script></body></html>\n''')
    with zipfile.ZipFile(OUT/'starplast-40-logos.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.glob('*.svg')):
            archive.write(path,path.name)
        archive.write(OUT/'manifest.json','manifest.json')
    print(f'Wrote forty proposals, 160 SVG assets and gallery: {OUT}')


if __name__=='__main__':
    main()
