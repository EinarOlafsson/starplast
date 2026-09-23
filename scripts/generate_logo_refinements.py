#!/usr/bin/env python3
"""Generate ten refinements of the approved Toxoplasma constellation as native SVGs.

Run with QT_QPA_PLATFORM=offscreen. Uses PyQt6 and Noto Sans Light for outlined
lettering and PNG previews. The active application logo is left to its own generator.
"""
from pathlib import Path
import html
import zipfile

from PyQt6 import QtGui

from generate_logo_options import outline, svg
from generate_constellation_logo import BODY, star_path, render

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/assets/logo-refinements'
SLENDER = ('M174 27C183 26 188 29 184 38C166 76 147 113 142 149'
           'C136 185 145 211 163 231C180 251 169 275 146 277'
           'C111 280 81 253 67 219C45 168 58 113 84 77C109 45 146 28 174 27Z')
BROAD = ('M176 30C185 28 191 31 186 41C165 79 151 113 149 146'
         'C147 183 159 205 180 230C196 251 179 276 155 278'
         'C112 281 77 256 59 220C37 174 45 121 73 83C101 47 148 29 176 30Z')
COMPACT = ('M175 39C185 38 190 42 185 51C163 83 150 118 151 151'
           'C152 181 165 202 183 224C198 245 182 268 156 270'
           'C110 273 76 249 60 217C40 175 47 123 76 87C102 56 146 39 175 39Z')
SWEEP = ('M177 28C186 27 191 31 186 40C161 81 145 119 145 152'
         'C145 191 163 210 182 233C198 254 175 279 150 278'
         'C109 277 76 251 60 216C40 171 50 115 80 77C108 43 149 28 177 28Z')
CAP = '<path d="M153 44Q163 46 178 53 M148 50Q160 52 175 60"/>'
CAP_SMALL = '<path d="M154 43Q167 45 178 53"/>'
MEMBRANE = '<path d="M125 46C88 64 60 102 58 145"/>'
LONG_MEMBRANE = '<path d="M124 49C89 70 62 110 59 150C57 178 64 202 73 218"/>'
A = [(151,70,6),(117,106,7),(85,147,7),(103,185,10),(80,218,5),(132,249,7),(137,217,5)]
B = [(146,75,6),(115,107,7),(91,143,6),(106,180,10),(86,216,6),(127,246,7)]
C = [(152,74,6),(115,99,6),(90,135,8),(117,170,9),(79,190,6),(103,233,7),(151,249,6)]
D = [(152,74,6),(114,103,7),(89,141,6),(121,174,8),(80,193,6),(106,234,8),(151,248,6),(147,216,5)]
CHAIN7 = [(0,1),(1,2),(2,3),(3,4),(4,5),(3,6),(6,5)]
CHAIN6 = [(0,1),(1,2),(2,3),(3,4),(4,5)]
OPTIONS = [
 ('quiet-arc','Quiet arc','A cleaner version of the original, with more space around the stars.',BODY,A,CHAIN7,CAP_SMALL,'',0),
 ('apical-line','Apical line','A slender body and six stars following one uninterrupted path.',SLENDER,B,CHAIN6,CAP,'',0),
 ('nucleus','Nucleus','The central star sits inside a small tilted nucleus ring.',BODY,A,CHAIN7,CAP_SMALL+'<ellipse cx="103" cy="185" rx="17" ry="21" transform="rotate(-14 103 185)"/>','',0),
 ('open-field','Open field','A wider interior gives the constellation a clearer zigzag shape.',BROAD,C,[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6)],CAP_SMALL,'',0),
 ('membrane','Membrane','A fine inner membrane traces the cell beside a branching network.',SWEEP,A,CHAIN7,CAP_SMALL+LONG_MEMBRANE,'',0),
 ('inclined','Inclined','A gentle tilt and a small apical cap give the mark a sense of movement.',SLENDER,B,CHAIN6,CAP_SMALL,'',14),
 ('linked-stars','Linked stars','Two small groups of stars meet across the centre of a rounded cell.',BROAD,D,[(0,1),(1,2),(0,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,3)],CAP_SMALL,'',0),
 ('bright-centre','Bright centre','One larger star anchors a quieter, sparse constellation.',SWEEP,[(151,70,5),(116,108,6),(86,146,5),(104,185,14),(86,220,5),(134,249,6)],[(0,1),(1,2),(2,3),(3,4),(4,5),(3,5)],CAP+MEMBRANE,'',-5),
 ('double-arc','Double arc','A softly rounded silhouette with a short, delicate membrane echo.',BODY,C,[(0,1),(1,2),(2,3),(3,4),(3,5),(5,6)],CAP_SMALL+MEMBRANE,'',-9),
 ('compact','Compact','A fuller body and seven balanced stars for smaller uses.',COMPACT,[(151,79,6),(117,107,7),(90,143,7),(112,180,9),(83,207,6),(113,242,7),(153,239,6)],[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(3,6)],'<path d="M150 56Q163 59 176 67"/>','',0),
]


def mark(option, inverse=False):
    """Draw one rounded cell outline with its constellation entirely inside."""
    _,_,_,body,stars,edges,details,_,angle=option
    ink,paper=('white','black') if inverse else ('black','white')
    links=''.join(f'<path d="M{stars[a][0]} {stars[a][1]}L{stars[b][0]} {stars[b][1]}"/>' for a,b in edges)
    nodes=''.join(f'<path d="{star_path(x,y,r)}"/>' for x,y,r in stars)
    return (f'<g transform="rotate({angle} 120 155)" fill="none" stroke="{ink}" '
            'stroke-linecap="round" stroke-linejoin="round">'
            f'<path d="{body}" stroke-width="1.65"/>'
            f'<g stroke-width=".85">{details}</g>'
            f'<g stroke-width="1.0">{links}</g>'
            f'<g stroke-width="1.2" fill="{paper}">{nodes}</g></g>')


def main():
    """Write the ten SVG sets, a numbered contact sheet, and a light/dark gallery."""
    app=QtGui.QGuiApplication([])
    OUT.mkdir(parents=True,exist_ok=True)
    word=outline('starplast',76)
    sheet=['<rect width="1800" height="1100" fill="white"/>',
           '<g color="black" transform="translate(62 67)">'+outline('starplast / constellation refinements',32)+'</g>',
           '<g color="black" transform="translate(62 102)">'+outline('Ten variations on the approved mark. Thin lines, rounded ends, connected stars.',17)+'</g>']
    cards=[]
    for i,option in enumerate(OPTIONS,1):
        slug,name,description,*_=option
        stem=f'{i:02}-{slug}'
        for inverse in (False,True):
            suffix='-white' if inverse else ''
            ink='white' if inverse else 'black'
            icon=svg(mark(option,inverse),240,310,f'Starplast {i:02} / {name}')
            lockup=svg(mark(option,inverse)+f'<g color="{ink}" transform="translate(240 186)">{word}</g>',620,310,f'Starplast {i:02} / {name} wordmark')
            (OUT/f'{stem}-icon{suffix}.svg').write_text(icon)
            (OUT/f'{stem}-wordmark{suffix}.svg').write_text(lockup)
        render(svg(mark(option),240,310,name),OUT/f'{stem}.png',480,620)
        x,y=45+345*((i-1)%5),157+449*((i-1)//5)
        sheet += [f'<g transform="translate({x+28} {y})">{mark(option)}</g>',
                  f'<g color="black" transform="translate({x+56} {y+343})">'+outline(f'{i:02} / {name}',18)+'</g>',
                  f'<g color="black" transform="translate({x+75} {y+376}) scale(.42)">{word}</g>']
        cards.append(f'''<article><h2>{i:02} / {name}</h2><div class="art">
<img class="light" src="{stem}-wordmark.svg" alt="{html.escape(description)}">
<img class="dark" src="{stem}-wordmark-white.svg" alt="{html.escape(description)}"></div>
<p>{description}</p><nav><a href="{stem}-icon.svg" download>Icon SVG</a>
<a href="{stem}-wordmark.svg" download>Wordmark SVG</a>
<a href="{stem}-icon-white.svg" download>White icon</a>
<a href="{stem}-wordmark-white.svg" download>White wordmark</a></nav></article>''')
    comparison=svg(''.join(sheet),1800,1100,'Starplast / ten Toxoplasma constellation refinements')
    (OUT/'comparison.svg').write_text(comparison)
    render(comparison,OUT/'comparison.png',3600,2200)
    page='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Starplast / ten constellation refinements</title><style>
*{box-sizing:border-box}body{margin:0;color:#111;background:#fff;font:16px/1.6 system-ui,sans-serif}main{max-width:1400px;margin:auto;padding:44px 28px}
h1{font-size:36px;font-weight:300;line-height:1.2}h2{font-size:20px;font-weight:400}p{max-width:840px}a{color:inherit;text-underline-offset:4px}
header nav,article nav{display:flex;gap:20px;flex-wrap:wrap}label{display:block;margin:26px 0;cursor:pointer}.grid{display:grid;grid-template-columns:1fr 1fr;gap:26px}
article{padding:24px;border:1px solid #ddd;border-radius:12px}article p{min-height:52px;font-size:15px}article nav{font-size:13px;gap:14px}.art img{display:block;width:100%;height:auto}.art .dark{display:none}
body:has(#dark:checked){color:#fff;background:#080808}body:has(#dark:checked) article{border-color:#444}body:has(#dark:checked) .light{display:none}body:has(#dark:checked) .art .dark{display:block}
@media(max-width:760px){.grid{grid-template-columns:1fr}main{padding:25px 16px}h1{font-size:28px}}
</style></head><body><main><header><h1>Ten constellation refinements</h1>
<p>Each version keeps the asymmetric Toxoplasma silhouette, softly rounded ends, thin monochrome lines, and an obvious connected constellation. The stars represent a gene network, not cell anatomy.</p>
<nav><a href="starplast-logo-refinements.zip" download>Download all SVGs</a><a href="comparison.png">Numbered comparison</a><a href="../toxoplasma-constellation/index.html">Current logo</a></nav>
<label><input type="checkbox" id="dark"> White on black</label></header><div class="grid">'''+''.join(cards)+'''</div></main></body></html>\n'''
    (OUT/'index.html').write_text(page)
    with zipfile.ZipFile(OUT/'starplast-logo-refinements.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(OUT.glob('*.svg')):archive.write(p,p.name)
    print(f'Wrote ten proposals and gallery to {OUT}')


if __name__=='__main__':
    main()
