#!/usr/bin/env python3
"""Forty black-and-white logo drafts, generated rather than drawn by hand.

Generated for three reasons. Forty hand-drawn variations would be forty variations on whichever one
came first; a script makes the six directions the task asks for actually distinct, and the seeds
mean a draft can be regenerated exactly after someone says "that one, but sparser". It also lets the
16-pixel constraint be a test rather than a hope: `tests/test_logos.py` renders every one of them at
16 and fails on anything that comes out as a smudge or a full square.

The six directions, and what each is arguing:

    constellation  the map IS a point cloud -- 8,140 genes where position means similarity
    apicoplast     the plastid the name refers to, drawn as nested ovals rather than a cell
    crescent       the tachyzoite, the most recognisable shape in the field
    orbit          rings, for the three levels of detail the map zooms through
    letter-S       an S built from points, so the mark reads as a letter and as data at once
    dark-field     the project's actual argument: most of this proteome has never been studied, so
                   most of the mark is empty and the few lit points are the ones anybody has looked
                   at. This is the family worth arguing for.

Black ink on a transparent ground: no greys, because a grey vanishes when the icon is printed or
shown at 16 pixels, and transparency so the same file works on the light and dark themes.
"""
from __future__ import annotations

import math
import os
import random

SIZE = 64                     # the drawing box; every shape is sized in these units
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "starplast", "data", "icons")

HEAD = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE} {SIZE}" '
        f'width="{SIZE}" height="{SIZE}">')
TAIL = "</svg>"


def wrap(body: str, title: str) -> str:
    """One draft, with its idea recorded in the file so the file explains itself."""
    return f"{HEAD}<title>{title}</title>{body}{TAIL}"


def dot(x, y, r, black=True) -> str:
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{"#000" if black else "#fff"}"/>'


# --------------------------------------------------------------------------- 1 constellation
def constellation(seed: int, n: int, links: int, r: float) -> str:
    """Points with a few edges between them: the map, at the scale where you can still count it.

    The edges are the nearest neighbours rather than random pairs, because a constellation whose
    lines cross at random reads as a scribble, and this application's whole point is that proximity
    means something.
    """
    rng = random.Random(seed)
    pts = [(rng.uniform(8, SIZE - 8), rng.uniform(8, SIZE - 8)) for _ in range(n)]
    body = []
    for i, (x, y) in enumerate(pts):
        near = sorted(range(len(pts)), key=lambda j: (pts[j][0] - x) ** 2 + (pts[j][1] - y) ** 2)
        for j in near[1:1 + links]:
            if j > i:
                body.append(f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{pts[j][0]:.2f}" '
                            f'y2="{pts[j][1]:.2f}" stroke="#000" stroke-width="1.1"/>')
    body += [dot(x, y, r) for x, y in pts]
    return "".join(body)


# --------------------------------------------------------------------------- 2 apicoplast
def apicoplast(rings: int, tilt: float, dotted: bool) -> str:
    """Nested ovals: the four membranes of the plastid, which is what the name refers to.

    Deliberately not a cell -- no nucleus, no apex, nothing that reads as a diagram of a parasite.
    """
    # Heavier and wider apart than the first draft of this family, which the 16-pixel test rejected
    # at 3% ink: four hairline ellipses at that size are a grey smear, and the whole point of the
    # constraint is that the mark has to survive being a tab icon.
    body = []
    for k in range(rings):
        rx, ry = 26 - k * 6.0, 17 - k * 3.6
        dash = ' stroke-dasharray="3.2 3.2"' if dotted and k % 2 else ""
        body.append(f'<ellipse cx="32" cy="32" rx="{rx:.1f}" ry="{ry:.1f}" fill="none" '
                    f'stroke="#000" stroke-width="{4.0 - k * 0.5:.2f}"{dash} '
                    f'transform="rotate({tilt:.0f} 32 32)"/>')
    body.append(dot(32, 32, 4.0))
    return "".join(body)


# --------------------------------------------------------------------------- 3 crescent
def crescent(fatness: float, points: int, seed: int) -> str:
    """The tachyzoite: a crescent, apex up, optionally made of points rather than a solid.

    Apex up for the same reason the cell diagram is: that is the end that goes in, and it is how the
    shape is drawn everywhere in the field.
    """
    # Fatter than the first draft, whose outline was a 1.6-unit stroke: at 16 pixels that is a third
    # of a pixel and the crescent disappeared entirely. The 16-pixel test is what said so.
    outer = "M30 5 C 52 18, 54 46, 30 59 C 42 44, 42 20, 30 5 Z"
    if not points:
        return f'<path d="{outer}" fill="#000"/>'
    rng = random.Random(seed)
    body = [f'<path d="{outer}" fill="none" stroke="#000" stroke-width="3.4"/>']
    for _ in range(points):
        t = rng.random()
        x = 32 + 13 * math.sin(t * math.pi) * fatness
        y = 6 + t * 52
        body.append(dot(x, y, 2.4))
    return "".join(body)


# --------------------------------------------------------------------------- 4 orbit
def orbit(rings: int, dots: int, seed: int) -> str:
    """Concentric rings with a point on each: the three levels of detail, as a scale of zoom."""
    # Thicker rings and bigger points than the first draft: three hairline circles came out at 2%
    # ink at 16 pixels, which is a faint ring and nothing else. The test is the whole reason the
    # weights in this file are what they are.
    rng = random.Random(seed)
    body = []
    for k in range(rings):
        r = 10 + k * 9
        body.append(f'<circle cx="32" cy="32" r="{r}" fill="none" stroke="#000" '
                    f'stroke-width="{4.2 - k * 0.6:.2f}"/>')
        for _ in range(dots):
            a = rng.uniform(0, 2 * math.pi)
            body.append(dot(32 + r * math.cos(a), 32 + r * math.sin(a), 3.6))
    body.append(dot(32, 32, 5.0))
    return "".join(body)


# --------------------------------------------------------------------------- 5 letter S
S_PATH = [(44, 14), (36, 10), (26, 11), (20, 17), (22, 25), (30, 30), (40, 34),
          (44, 41), (42, 49), (34, 54), (24, 54), (18, 50)]


def letter_s(radius: float, stroke: float, joined: bool) -> str:
    """An S built from points, so the mark reads as a letter at a glance and as data up close."""
    body = []
    if joined:
        d = " ".join(("M" if i == 0 else "L") + f"{x} {y}" for i, (x, y) in enumerate(S_PATH))
        body.append(f'<path d="{d}" fill="none" stroke="#000" stroke-width="{stroke}" '
                    f'stroke-linecap="round" stroke-linejoin="round"/>')
    body += [dot(x, y, radius) for x, y in S_PATH]
    return "".join(body)


# --------------------------------------------------------------------------- 6 dark field
def dark_field(seed: int, lit: int, total: int, inverted: bool) -> str:
    """Most of the field is dark: the argument the whole application is built on.

    5,574 of 8,140 genes are named in no paper at all. So most of the mark is empty and only a few
    points are lit -- and in the inverted drafts the ground is solid black with the studied genes
    knocked out of it, which is the same claim said the other way round.
    """
    rng = random.Random(seed)
    body = []
    if inverted:
        body.append(f'<rect width="{SIZE}" height="{SIZE}" rx="8" fill="#000"/>')
    for i in range(total):
        x, y = rng.uniform(6, SIZE - 6), rng.uniform(6, SIZE - 6)
        if i < lit:
            body.append(dot(x, y, 4.2, black=not inverted))
        elif not inverted:
            # The unstudied genes are there and they are faint: a ring rather than a filled point.
            # Faint, not invisible -- the first draft used a 0.6-unit stroke, which at 16 pixels is
            # nothing at all, and the mark became six dots on an empty square. The claim is that
            # most of the field is UNSTUDIED, not that most of it is absent.
            body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.2" fill="none" stroke="#000" '
                        f'stroke-width="1.4"/>')
    return "".join(body)


# --------------------------------------------------------------------------- 7 tilted galaxy
def galaxy(orbits: int, tilt: float, planets, seed: int, arms: int = 0, solid: bool = True,
           dotted: bool = False) -> str:
    """A tilted solar system with the parasite at the centre of it.

    The mark the name actually describes: a star map with *Toxoplasma* in the middle. Tilted rather
    than face-on because a face-on system is a target and a tilted one is a system seen from
    somewhere -- which is what a 3D map you can orbit around is.

    The centre is a small tachyzoite, apex up. It is the one place a cell shape belongs in this set:
    not as a diagram of an organism but as the body the whole map turns around.
    """
    rng = random.Random(seed)
    body = []
    for k in range(orbits):
        rx = 17 + k * 9.0
        ry = rx * 0.42                      # the flattening IS the tilt, seen edge-on-ish
        dash = ' stroke-dasharray="3.4 3.4"' if dotted and k % 2 else ""
        body.append(f'<ellipse cx="32" cy="32" rx="{rx:.1f}" ry="{ry:.1f}" fill="none" '
                    f'stroke="#000" stroke-width="{3.4 - k * 0.45:.2f}"{dash} '
                    f'transform="rotate({tilt:.0f} 32 32)"/>')
        # Planets sit ON their orbit, which means the ellipse's own parametric point, rotated with
        # it -- scattered at random they read as dirt rather than as bodies going round something.
        for _ in range(planets[k] if k < len(planets) else 0):
            a = rng.uniform(0, 2 * math.pi)
            px, py = rx * math.cos(a), ry * math.sin(a)
            th = math.radians(tilt)
            x = 32 + px * math.cos(th) - py * math.sin(th)
            y = 32 + px * math.sin(th) + py * math.cos(th)
            body.append(dot(x, y, 3.4))
    for k in range(arms):
        # A spiral arm, for the galaxy reading rather than the solar-system one.
        a0 = 2 * math.pi * k / max(arms, 1)
        pts = []
        for s in range(9):
            a = a0 + s * 0.34
            r = 7 + s * 3.6
            px, py = r * math.cos(a), r * 0.45 * math.sin(a)
            th = math.radians(tilt)
            pts.append((32 + px * math.cos(th) - py * math.sin(th),
                        32 + px * math.sin(th) + py * math.cos(th)))
        d = " ".join(("M" if i == 0 else "L") + f"{x:.1f} {y:.1f}" for i, (x, y) in enumerate(pts))
        body.append(f'<path d="{d}" fill="none" stroke="#000" stroke-width="2.2" '
                    f'stroke-linecap="round"/>')
    # Toxoplasma at the centre, apex up -- with the orbits knocked out behind it. Without the
    # knockout the innermost orbit runs straight through the parasite and the middle of the mark
    # reads as a blob rather than as a body the system turns around, which is the whole idea.
    body.append('<ellipse cx="32" cy="32.5" rx="9.5" ry="14" fill="#fff"/>')
    heart = "M32 17 C 43 25, 44 41, 32 48 C 38 39, 38 26, 32 17 Z"
    body.append(f'<path d="{heart}" fill="#000"/>' if solid else
                f'<path d="{heart}" fill="none" stroke="#000" stroke-width="3.2"/>')
    return "".join(body)


def galaxy_body(orbits: int, tilt: float, planets, seed: int, rhoptries: int = 4) -> str:
    """The parasite IS the innermost orbit, with rhoptries added to it.

    The other galaxy drafts stamp a crescent over the middle of the system, which is two drawings on
    top of each other -- the mark says "a parasite" and "a system" in the same place and the eye has
    to separate them. Here the innermost ellipse is the parasite: same tilt, same family of shapes,
    so the system resolves INTO the organism as it gets smaller instead of being interrupted by it.

    What makes it a parasite rather than a bean is the rhoptries: club-shaped organelles converging
    on the apical end, which is the most recognisable thing an apicomplexan has and the reason the
    phylum is named after it. They are drawn along the long axis of the inner ellipse, narrow ends
    meeting at the apex, so the tilt of the system is also the tilt of the cell.
    """
    rng = random.Random(seed)
    th = math.radians(tilt)

    def place(px, py):
        """A point in the tilted frame of the system."""
        return (32 + px * math.cos(th) - py * math.sin(th),
                32 + px * math.sin(th) + py * math.cos(th))

    body = []
    for k in range(orbits):
        rx = 28 + k * 9.5
        ry = rx * 0.42
        body.append(f'<ellipse cx="32" cy="32" rx="{rx:.1f}" ry="{ry:.1f}" fill="none" '
                    f'stroke="#000" stroke-width="{3.2 - k * 0.5:.2f}" '
                    f'transform="rotate({tilt:.0f} 32 32)"/>')
        for _ in range(planets[k] if k < len(planets) else 0):
            a = rng.uniform(0, 2 * math.pi)
            body.append(dot(*place(rx * math.cos(a), ry * math.sin(a)), 3.4))

    # The innermost "orbit" is the cell: an ellipse of the same family, drawn heavier so it reads as
    # a body rather than as one more ring, and filled white so the orbits behind it do not run
    # through the organelles inside it.
    cell_rx, cell_ry = 20.0, 10.0
    body.append(f'<ellipse cx="32" cy="32" rx="{cell_rx}" ry="{cell_ry}" fill="#fff" '
                f'stroke="#000" stroke-width="3.2" transform="rotate({tilt:.0f} 32 32)"/>')

    # Rhoptries: clubs converging on the apical end, which is the left end of the long axis here.
    apex = (-cell_rx * 0.94, 0.0)
    for i in range(rhoptries):
        # Fanned across the short axis, each tapering from a bulb at the back to the apex.
        spread = (i - (rhoptries - 1) / 2) / max(rhoptries - 1, 1)
        bulb = (cell_rx * 0.34, spread * cell_ry * 0.55)
        waist = (-cell_rx * 0.20, spread * cell_ry * 0.30)
        d = (f"M{place(*apex)[0]:.1f} {place(*apex)[1]:.1f} "
             f"Q{place(*waist)[0]:.1f} {place(*waist)[1]:.1f} "
             f"{place(*bulb)[0]:.1f} {place(*bulb)[1]:.1f}")
        body.append(f'<path d="{d}" fill="none" stroke="#000" stroke-width="1.9" '
                    f'stroke-linecap="round"/>')
        body.append(dot(*place(*bulb), 2.0))
    # The conoid: the apical cap the rhoptries discharge through, and the point the whole cell aims.
    body.append(dot(*place(*apex), 2.4))
    return "".join(body)


def drafts() -> list:
    """The drafts, as (name, title, body). Seven directions, not forty variations on one."""
    out = []
    for i, (n, links, r) in enumerate([(9, 1, 3.2), (14, 1, 2.8), (20, 2, 2.4), (28, 1, 2.0),
                                       (12, 2, 3.0), (18, 3, 2.2), (24, 2, 2.6), (16, 1, 3.4)]):
        out.append((f"constellation_{i + 1}", "starplast — the map is a point cloud",
                    constellation(seed=100 + i, n=n, links=links, r=r)))
    for i, (rings, tilt, dotted) in enumerate([(3, 0, False), (4, 0, False), (3, -20, False),
                                               (4, -20, True), (2, 15, False), (3, 30, True)]):
        out.append((f"apicoplast_{i + 1}", "starplast — the plastid the name refers to",
                    apicoplast(rings, tilt, dotted)))
    for i, (fat, pts) in enumerate([(1.0, 0), (1.0, 14), (0.8, 20), (1.2, 10), (1.0, 26),
                                    (0.9, 34)]):
        out.append((f"crescent_{i + 1}", "starplast — the tachyzoite, apex up",
                    crescent(fat, pts, seed=200 + i)))
    for i, (rings, dots) in enumerate([(3, 1), (3, 2), (2, 3), (4, 1), (3, 4), (2, 1)]):
        out.append((f"orbit_{i + 1}", "starplast — three levels of detail",
                    orbit(rings, dots, seed=300 + i)))
    for i, (rad, stroke, joined) in enumerate([(2.6, 0, False), (2.2, 2.0, True), (3.0, 0, False),
                                               (1.8, 3.0, True), (2.4, 1.4, True), (3.4, 0, False),
                                               (2.0, 2.6, True)]):
        out.append((f"letter_s_{i + 1}", "starplast — an S made of genes",
                    letter_s(rad, stroke, joined)))
    for i, (lit, total, inv) in enumerate([(6, 40, False), (4, 60, False), (8, 30, False),
                                           (5, 45, True), (3, 50, True), (7, 36, True),
                                           (9, 24, False)]):
        out.append((f"dark_field_{i + 1}", "starplast — most of this proteome is unstudied",
                    dark_field(seed=400 + i, lit=lit, total=total, inverted=inv)))
    # The seventh direction, asked for after the first forty were drawn: a tilted system with the
    # parasite at the centre of it. Eight, because it is the direction the name argues for and one
    # variant of it would be a token.
    for i, (orbits, tilt, planets, arms, solid, dotted) in enumerate([
            (2, -22, (1, 1), 0, True, False),
            (3, -22, (1, 1, 1), 0, True, False),
            (3, -18, (1, 2, 1), 0, True, False),
            (2, -30, (2, 1), 0, False, False),
            (3, -22, (1, 1, 2), 0, True, True),
            (3, -14, (1, 1, 1), 2, True, False),
            (2, -26, (1, 2), 3, True, False),
            (4, -20, (1, 1, 1, 1), 0, True, False)]):
        out.append((f"galaxy_{i + 1}", "starplast — a tilted system with Toxoplasma at its centre",
                    galaxy(orbits, tilt, planets, seed=500 + i, arms=arms, solid=solid,
                           dotted=dotted)))
    # The last of the series: the parasite is not stamped over the system, it IS the innermost
    # orbit, with rhoptries on it. The system resolves into the organism rather than being
    # interrupted by it.
    out.append(("galaxy_9", "starplast — the innermost orbit is the parasite",
                galaxy_body(orbits=2, tilt=-22, planets=(1, 1), seed=600, rhoptries=3)))
    return out


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    made = drafts()
    for i, (name, title, body) in enumerate(made, start=1):
        path = os.path.join(OUT, f"logo_{i:02d}_{name}.svg")
        with open(path, "w", encoding="utf8") as fh:
            fh.write(wrap(body, title))
    print(f"{len(made)} logo drafts in {OUT}")
    return len(made)


if __name__ == "__main__":
    main()
