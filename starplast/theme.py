#!/usr/bin/env python3
"""Palettes, color maps and point styles.

Structured after spaCR's theme module so the two tools feel like one toolkit: a palette is a dict of
named roles, every theme defines the same roles, and colors that are a *function* of the theme are
derived rather than typed in. The rule that matters is spaCR's: a hex written into a widget is a hex
that stays dark on the light theme, so nothing here hardcodes a color at a call site.

Four themes. `dark` and `light` are the general pair; `slate` is a low-contrast dark for long sessions
staring at a 3D scatter, and `paper` is a high-contrast light intended for figures and screenshots.

Color maps are separated by what they are *for*, because using the wrong kind is the commonest way to
mislead with a map of this sort:

    sequential   an ordered quantity with a meaningful zero or floor (expression, pLDDT)
    diverging    a quantity with a meaningful midpoint (fitness scores, log ratios, residuals)
    categorical  unordered classes (compartment, cluster id, stage)

A diverging map on an unordered category invents an order; a sequential map on a residual hides its
sign. `CMAPS` therefore records the kind, and the UI offers only the appropriate ones for a column.
"""
from __future__ import annotations

from PyQt6 import QtCore, QtGui, QtWidgets

# --------------------------------------------------------------------------- palettes
DARK = {
    "bg": "#0b0d10", "page": "#101317", "surface": "#141820", "surface_alt": "#1a1f28",
    "surface_hi": "#222833", "border": "#2b323d", "border_soft": "#1c222b",
    "fg": "#eef2f6", "fg_muted": "#9aa6b4", "fg_dim": "#68737f",
    "accent": "#43c6d8", "accent_hi": "#66dcea", "accent_lo": "#2a9fb0", "accent_soft": "#12303a",
    "success": "#4bb96b", "warning": "#d9a13a", "error": "#e5615a", "info": "#43c6d8",
}
LIGHT = {
    "bg": "#f7f9fa", "page": "#ffffff", "surface": "#ffffff", "surface_alt": "#f0f4f6",
    "surface_hi": "#e6ecf0", "border": "#d3dde3", "border_soft": "#e8eef1",
    "fg": "#131a20", "fg_muted": "#53626e", "fg_dim": "#7b8794",
    "accent": "#0f7d8c", "accent_hi": "#0b6b78", "accent_lo": "#095a66", "accent_soft": "#d9eef1",
    "success": "#1f7a3d", "warning": "#a8640a", "error": "#b3261e", "info": "#0f7d8c",
}
SLATE = {
    "bg": "#16191d", "page": "#1b1f24", "surface": "#20252b", "surface_alt": "#262c33",
    "surface_hi": "#2e353d", "border": "#39414a", "border_soft": "#252b32",
    "fg": "#dfe5ea", "fg_muted": "#9aa4ae", "fg_dim": "#6f7983",
    "accent": "#c58f4a", "accent_hi": "#dba765", "accent_lo": "#a4733a", "accent_soft": "#33281a",
    "success": "#5a9e6a", "warning": "#c9a227", "error": "#cf6a63", "info": "#c58f4a",
}
PAPER = {
    "bg": "#ffffff", "page": "#ffffff", "surface": "#ffffff", "surface_alt": "#f4f4f4",
    "surface_hi": "#e9e9e9", "border": "#c9c9c9", "border_soft": "#e2e2e2",
    "fg": "#000000", "fg_muted": "#3c3c3c", "fg_dim": "#6a6a6a",
    "accent": "#0044cc", "accent_hi": "#0036a3", "accent_lo": "#002b80", "accent_soft": "#dde6ff",
    "success": "#0a6b2b", "warning": "#8a5200", "error": "#a30d0d", "info": "#0044cc",
}

_PALETTES = {"dark": DARK, "light": LIGHT, "slate": SLATE, "paper": PAPER}
THEMES = tuple(_PALETTES)
DARK_THEMES = ("dark", "slate")


def palette_for(theme: str = "dark") -> dict:
    """The color palette for a named theme."""
    return _PALETTES.get(theme, DARK)


def is_dark(theme: str) -> bool:
    """Whether a palette is a dark one, which decides contrast choices elsewhere."""
    return theme in DARK_THEMES


def rgbf(hex_color: str, alpha: float = 1.0) -> tuple:
    """Hex to the 0-1 RGBA tuple pyqtgraph's GL items want."""
    h = hex_color.lstrip("#")
    return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, alpha)


def background(theme: str) -> tuple:
    """GL viewport background — the page color, not a literal, so themes stay coherent."""
    return rgbf(palette_for(theme)["bg"])


def unknown_color(theme: str) -> tuple:
    """The gray for 'unknown, not absent'. Derived from fg_dim so it reads on either ground."""
    return rgbf(palette_for(theme)["fg_dim"])


# --------------------------------------------------------------------------- color maps
# name -> (kind, the name pyqtgraph will accept). pyqtgraph ships only five matplotlib maps; its
# diverging maps come from the CET (Peter Kovesi) set, which is perceptually uniform and a better
# choice than coolwarm anyway. Offering a name pyqtgraph cannot build raises FileNotFoundError deep in
# a paint call, which is how this list was wrong the first time.
CMAPS = {
    "viridis": ("sequential", "viridis"), "magma": ("sequential", "magma"),
    "plasma": ("sequential", "plasma"), "inferno": ("sequential", "inferno"),
    "cividis": ("sequential", "cividis"), "CET-L11": ("sequential", "CET-L11"),
    "CET-D1A": ("diverging", "CET-D1A"), "CET-D9": ("diverging", "CET-D9"),
    "CET-D4": ("diverging", "CET-D4"), "CET-D13": ("diverging", "CET-D13"),
    "CET-C7s": ("categorical", "CET-C7s"), "CET-C6": ("categorical", "CET-C6"),
}
DEFAULT_CMAP = {"sequential": "viridis", "diverging": "CET-D1A", "categorical": "CET-C7s"}


def resolve_cmap(name: str):
    """A pyqtgraph colormap, falling back rather than raising inside a paint call."""
    import pyqtgraph as pg
    for candidate in (name, DEFAULT_CMAP["sequential"]):
        try:
            cm = pg.colormap.get(candidate)
            if cm is not None:
                return cm
        except Exception:
            continue
    return pg.colormap.get("viridis")


def cmaps_of(kind: str) -> list:
    """The color maps appropriate to one kind of data: sequential, diverging or categorical."""
    return [n for n, (k, _) in CMAPS.items() if k == kind]


def kind_for_column(values) -> str:
    """Which family of color map a column deserves.

    Diverging only when the data actually straddles a midpoint -- a residual or a score that goes both
    ways. Applying one to a strictly positive quantity invents a midpoint that is not there.
    """
    import numpy as np
    import pandas as pd
    raw = pd.Series(values)
    s = pd.to_numeric(raw, errors="coerce")
    if s.notna().mean() < 0.5:                      # mostly non-numeric -> labels
        return "categorical"
    s = s.dropna()
    if s.empty:
        return "categorical"
    # Few distinct values is not enough on its own: a fitness score can take five values in a small
    # selection and is still a continuous quantity. Categorical means few distinct values that are all
    # whole numbers, which is what an encoded class looks like.
    if s.nunique() <= 12 and float((s % 1 == 0).mean()) == 1.0:
        return "categorical"
    lo, hi = float(s.min()), float(s.max())
    if lo < 0 < hi and abs(lo) / max(hi, 1e-9) > 0.15:
        return "diverging"
    return "sequential"


# --------------------------------------------------------------------------- point styles
POINT_STYLES = {
    "small": {"size": 3.5, "alpha": 0.85},
    "standard": {"size": 5.0, "alpha": 0.95},
    "large": {"size": 8.0, "alpha": 0.95},
    "halo": {"size": 9.0, "alpha": 0.45},        # for density; overlap reads as brightness
    "pinpoint": {"size": 2.0, "alpha": 1.0},     # dense maps, 8,000+ points
}
DEFAULT_POINT_STYLE = "standard"

# How strongly distance from the camera fades a point. 0 disables it. A 3D scatter drawn without depth
# cueing reads as a flat disc of color: near and far points are equally bright, so the eye has nothing
# to build a shape from. This is the cheapest thing that makes the map three-dimensional.
DEPTH_FADE = 0.55        # fraction of alpha the farthest point keeps
DEPTH_SHRINK = 0.45      # fraction of size the farthest point keeps
EDGE_DEPTH_FADE = 0.18   # edges fade harder than points: they are the clutter, and the near
                         # neighbourhood is what the fade is meant to make readable


def is_light(theme_or_palette) -> bool:
    """Whether a theme's ground is light. Decides which direction contrast has to go."""
    p = theme_or_palette if isinstance(theme_or_palette, dict) else palette_for(theme_or_palette)
    r, g, b = rgbf(p["page"])[:3]
    return (0.2126 * r + 0.7152 * g + 0.0722 * b) > 0.5


def depth_t(xyz, camera_pos, rng=None):
    """Distance from the camera, normalized to 0 (nearest) .. 1 (farthest).

    `rng` is the (min, max) distance to normalize against. Pass the WHOLE cloud's range when cueing a
    subset -- edges touching a point must fade by the same amount as that point, and a subset
    normalized against its own extent does not agree with the full set at the same location.
    """
    import numpy as np
    d = np.linalg.norm(np.asarray(xyz, dtype=float) - np.asarray(camera_pos, dtype=float), axis=1)
    lo, hi = (float(d.min()), float(d.max())) if rng is None else rng
    if hi - lo <= 1e-9:                       # camera equidistant from everything: no cue, no divide by zero
        return np.zeros_like(d)
    return np.clip((d - lo) / (hi - lo), 0.0, 1.0)


def depth_cue(xyz, camera_pos, alpha, size, fade=DEPTH_FADE, shrink=DEPTH_SHRINK, rng=None):
    """Scale alpha and size by distance from the camera. Returns (alpha_array, size_array)."""
    t = depth_t(xyz, camera_pos, rng)
    return (alpha * (1.0 - t * (1.0 - fade)),
            size * (1.0 - t * (1.0 - shrink)))


def categorical_colors(n: int, theme: str = "dark", cmap: str | None = None) -> list:
    """`n` distinct colors for unordered classes, legible on this theme's ground.

    A palette tuned for a dark background washes out on a light one -- the same hues that read as
    saturated against near-black read as pastel against near-white. Rather than ship two hand-tuned
    lists, the colors are darkened for light themes and lightened for dark ones, so any chosen
    categorical map stays legible on either ground.
    """
    import numpy as np
    try:
        from matplotlib import colormaps
        name = CMAPS.get(cmap, (None, "tab20"))[1] if cmap else "tab20"
        base = colormaps.get_cmap(name)
        cols = [tuple(base(i / max(n - 1, 1))[:3]) for i in range(n)]
    except Exception:
        rng = np.random.default_rng(0)
        cols = [tuple(rng.uniform(0.25, 0.9, 3)) for _ in range(n)]
    if is_dark(theme):
        return [tuple(min(1.0, c * 1.12 + 0.05) for c in col) for col in cols]
    return [tuple(max(0.0, c * 0.72) for c in col) for col in cols]

# Rendering modes for the scatter itself.
#   occlude    translucent + depth test; nearer points hide farther ones (the correct default)
#   additive   overlapping points sum; reads as density but saturates to white in dense regions
POINT_MODES = ("occlude", "additive")


def gl_options(mode: str) -> str:
    """pyqtgraph GL option set for a point mode.

    `additive` is offered because density is sometimes what you want to see, but it is not the default:
    with 8,140 genes in dense UMAP clusters it saturates every color mode to white and destroys the
    encoding entirely.
    """
    return "additive" if mode == "additive" else "translucent"


# --------------------------------------------------------------------------- glass
# spaCR's look for everything that floats: menus, tooltips, drop-down lists, the help search results
# and Starplast's own windows are rounded panes of translucent black -- translucent white on a light
# theme -- instead of solid grey rectangles. `starplast.glass` decides whether the display can show a
# translucent window at all; these functions only say what colour the pane is either way.

#: How much of the pane is glass, per role. Tooltips and dialogs carry the most reading and are the
#: most opaque; a menu is glanced at, and letting a little of the map show through is what makes it
#: read as floating over it rather than pasted on. Every value keeps body text above 12:1 contrast.
GLASS_ALPHA = {"menu": 0.86, "tooltip": 0.92, "popup": 0.90, "dialog": 0.90, "pane": 0.35,
               "bar": 0.45}

#: Corner radii, spaCR's scale: `sm` for a menu item or a field, `md` for a menu or a list, `lg` for
#: the search results, `card` for a whole window.
RADIUS = {"sm": 4, "md": 8, "lg": 10, "card": 14}

#: spaCR's type. Open Sans where it is installed, the platform sans-serif where it is not.
FONT_FAMILY = '"Open Sans", "Segoe UI", "Helvetica Neue", sans-serif'


def glass_rgba(theme: str, role: str = "menu", translucent: bool = True) -> tuple:
    """The glass colour for one role as an (r, g, b, a) tuple of floats in 0-1.

    Black on a dark theme, white on a light one, at `GLASS_ALPHA[role]`. When the display cannot
    composite (`translucent=False`) the same pane is returned OPAQUE, as it would look laid over the
    theme's own background: a near-black on a dark theme rather than a black square with garbage in
    its corners.
    """
    alpha = GLASS_ALPHA.get(role, GLASS_ALPHA["menu"])
    ink = 0.0 if is_dark(theme) else 1.0
    if translucent:
        return (ink, ink, ink, alpha)
    under = rgbf(palette_for(theme)["bg"])[:3]
    return tuple(ink * alpha + c * (1.0 - alpha) for c in under) + (1.0,)


def glass(theme: str, role: str = "menu", translucent: bool = True) -> str:
    """The glass colour for one role, as a stylesheet colour (`rgba(...)`, or `#rrggbb` when opaque)."""
    r, g, b, a = glass_rgba(theme, role, translucent)
    if a >= 1.0:
        return "#{:02x}{:02x}{:02x}".format(*(int(round(v * 255)) for v in (r, g, b)))
    return f"rgba({int(round(r * 255))}, {int(round(g * 255))}, {int(round(b * 255))}, {a:.3f})"


def rim(theme: str, translucent: bool = True) -> str:
    """The hairline round a glass pane: a faint light edge on black glass, a faint dark one on white.

    A '#aarrggbb' string, which both a stylesheet and `QColor` read. Opaque displays get the
    palette's own border, since a translucent edge over an opaque pane has nothing to show through.
    """
    if not translucent:
        return palette_for(theme)["border"]
    return "#1fffffff" if is_dark(theme) else "#24000000"


# --------------------------------------------------------------------------- stylesheet
def rgba(hex_color: str, alpha: float) -> str:
    """`rgba(r, g, b, a)` from a hex colour, for a panel that has to let the background through."""
    r, g, b = (int(round(v * 255)) for v in rgbf(hex_color)[:3])
    return f"rgba({r}, {g}, {b}, {max(0.0, min(1.0, float(alpha))):.3f})"


def stylesheet(theme: str = "dark", container_opacity: float = 1.0,
               text_scale: float = 1.0, translucent: bool = True) -> str:
    """Qt stylesheet for one theme. Every color comes from the palette, never a literal.

    `container_opacity` below 1 lets whatever is painted behind the window -- the drifting blob
    field -- show through the panels. Only the CONTAINERS take it: a translucent field would put
    moving colour behind text somebody is trying to read, and the point of the background is that it
    is behind things.

    `translucent` says whether the display can composite translucent windows (see
    `starplast.glass.compositing_available`). Menus, tooltips, drop-down lists and glass windows are
    translucent black when it can, and the same colour made opaque when it cannot.
    """
    p = palette_for(theme)
    CONTAINER = rgba(p["page"], container_opacity)
    # The sizes are scaled here as well as on the application font, because a stylesheet font-size
    # BEATS the application font -- which is why the text-size setting appeared to do nothing at all:
    # the font was changing and this rule was overriding it on every widget.
    BODY = max(int(round(13 * text_scale)), 6)
    SMALL = max(int(round(11 * text_scale)), 5)
    MENU = max(int(round(12 * text_scale)), 5)
    MENU_GLASS = glass(theme, "menu", translucent)
    TIP_GLASS = glass(theme, "tooltip", translucent)
    POPUP_GLASS = glass(theme, "popup", translucent)
    PANE_GLASS = glass(theme, "pane", True)       # inside a window: always composited by Qt itself
    BAR_GLASS = glass(theme, "bar", True)
    RIM = rim(theme, translucent)
    INNER_RIM = rim(theme, True)
    ON_ACCENT = p["bg"] if is_dark(theme) else "#ffffff"
    # The field grey is dark in every theme, so its text is light in every theme. With the light
    # themes' dark text on it, what was typed into a field on `light` or `paper` could not be read.
    FIELD_TEXT = p["fg"] if is_dark(theme) else DARK["fg"]
    # The menu bar is the one bar in the window everything else hangs from. ONE FLAT COLOUR, as in
    # spaCR, and never `transparent` for its items: on macOS a transparent item repaints the window's
    # own colour first and shows as a box behind the word being pointed at.
    BAR = p["bg"] if is_dark(theme) else p["surface"]
    R = RADIUS
    return f"""
    QWidget {{ background: {p['bg']}; color: {p['fg']}; font-family: {FONT_FAMILY};
               font-size: {BODY}px; selection-background-color: {p['accent']};
               selection-color: {p['bg'] if is_dark(theme) else p['page']}; }}
    QMainWindow, QDialog {{ background: {p['bg']}; }}
    QDockWidget {{ background: {p['surface']}; color: {p['fg']}; titlebar-close-icon: none; }}
    QDockWidget::title {{ background: {BAR_GLASS}; padding: 7px 10px;
                          border-bottom: 1px solid {INNER_RIM};
                          font-weight: 600; color: {p['fg_muted']}; }}
    QGroupBox {{ border: 1px solid {p['border']}; border-radius: 6px; margin-top: 16px;
                 padding-top: 10px; background: {p['surface']}; }}
    QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 5px;
                        color: {p['fg_muted']}; font-size: {SMALL}px;
                        text-transform: uppercase; letter-spacing: 1px; }}
    QLabel {{ background: transparent; }}
    QPushButton {{ background: {p['surface_hi']}; border: 1px solid {p['border']};
                   border-radius: 5px; padding: 6px 13px; color: {p['fg']}; }}
    QPushButton:hover {{ background: {p['accent_soft']}; border-color: {p['accent']}; }}
    QPushButton:pressed {{ background: {p['accent_lo']}; color: {p['page']}; }}
    QPushButton:disabled {{ color: {p['fg_dim']}; background: {p['surface_alt']}; }}
    QPushButton[primary="true"] {{ background: {p['accent']}; border-color: {p['accent']};
                                   color: {p['bg'] if is_dark(theme) else '#ffffff'};
                                   font-weight: 600; }}
    QPushButton[primary="true"]:hover {{ background: {p['accent_hi']}; }}
    /* One dark grey for every field, rather than the panel colour of wherever it happens to sit:
       asked for directly, and it also makes a field look like a field on a themed background --
       with the ambient blobs behind a translucent panel, a field the colour of its container
       disappears into whatever is drifting past. */
    QDockWidget > QWidget, QTabWidget::pane, QGroupBox {{ background: {CONTAINER}; }}
    QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit, QPlainTextEdit, QTextEdit, QAbstractSpinBox {{
        background: {FIELD_GREY}; border: 1px solid {p['border']};
        border-radius: 5px; padding: 5px 8px; color: {FIELD_TEXT}; }}
    QComboBox:focus, QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
        border-color: {p['accent']}; }}
    QCheckBox, QRadioButton {{ spacing: 7px; background: transparent; }}
    QCheckBox::indicator, QRadioButton::indicator {{ width: 15px; height: 15px;
        border: 1px solid {p['border']}; border-radius: 3px; background: {p['surface_hi']}; }}
    QCheckBox::indicator:checked, QRadioButton::indicator:checked {{
        background: {p['accent']}; border-color: {p['accent']}; }}
    QTabWidget::pane {{ border: 1px solid {p['border']}; border-radius: 6px;
                        background: {p['surface']}; top: -1px; }}
    QTabBar::tab {{ background: transparent; color: {p['fg_muted']};
                    padding: 8px 15px; border-bottom: 2px solid transparent; }}
    QTabBar::tab:selected {{ color: {p['fg']}; border-bottom-color: {p['accent']};
                             font-weight: 600; }}
    QTabBar::tab:hover {{ color: {p['fg']}; }}
    QTableView, QTreeView, QListWidget {{ background: {p['page']};
        alternate-background-color: {p['surface_alt']}; gridline-color: {p['border_soft']};
        border: 1px solid {p['border']}; border-radius: 6px; }}
    QHeaderView::section {{ background: {p['surface_alt']}; color: {p['fg_muted']};
        padding: 6px 8px; border: none; border-bottom: 1px solid {p['border']};
        font-size: {SMALL}px; text-transform: uppercase; letter-spacing: .6px; }}
    QTextBrowser {{ background: {p['page']}; border: 1px solid {p['border']};
                    border-radius: 6px; padding: 4px; }}
    QProgressBar {{ background: {p['surface_hi']}; border: none; border-radius: 4px;
                    height: 6px; text-align: center; color: transparent; }}
    QProgressBar::chunk {{ background: {p['accent']}; border-radius: 4px; }}
    QStatusBar {{ background: {BAR_GLASS}; color: {p['fg_muted']};
                  border-top: 1px solid {INNER_RIM}; }}
    QScrollBar:vertical {{ background: transparent; width: 10px; margin: 0; }}
    QScrollBar::handle:vertical {{ background: {p['border']}; border-radius: 5px;
                                   min-height: 26px; }}
    QScrollBar::handle:vertical:hover {{ background: {p['fg_dim']}; }}
    QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
    QSplitter::handle {{ background: {p['border_soft']}; }}

    /* ---- menu bar and menus, formatted like spaCR's: a flat bar whose words light up rather than
       sitting on a plate, and menus that are rounded panes of glass with rounded item pills. ---- */
    QMenuBar {{ background: {BAR}; color: {p['fg_muted']}; padding: 3px 6px;
                border-bottom: 1px solid {p['border_soft']}; font-size: {MENU}px; }}
    QMenuBar::item {{ background: {BAR}; padding: 5px 10px; border-radius: {R['sm']}px; }}
    QMenuBar::item:selected, QMenuBar::item:pressed {{ background: {BAR}; color: {p['accent']}; }}
    QMenu {{ background: {MENU_GLASS}; color: {p['fg']}; border: 1px solid {RIM};
             border-radius: {R['md']}px; padding: 5px; font-size: {BODY}px; }}
    QMenu::item {{ background: transparent; padding: 5px 24px 5px 12px;
                   border-radius: {R['sm']}px; }}
    QMenu::item:selected {{ background: {p['accent']}; color: {ON_ACCENT}; }}
    QMenu::item:disabled {{ color: {p['fg_dim']}; }}
    QMenu::separator {{ height: 1px; background: {RIM}; margin: 5px 10px; }}
    QMenu::indicator {{ width: 14px; height: 14px; left: 7px; }}
    QMenu::right-arrow {{ right: 8px; }}

    /* ---- tooltips and drop-down lists: the same glass, so nothing that floats is grey. ---- */
    QToolTip {{ background: {TIP_GLASS}; color: {p['fg']}; border: 1px solid {RIM};
                border-radius: {R['md'] - 2}px; padding: 6px 9px; font-size: {MENU}px;
                opacity: 255; }}
    QComboBox QAbstractItemView {{ background: {POPUP_GLASS}; border: 1px solid {RIM};
                                   border-radius: {R['md']}px; padding: 4px; outline: 0;
                                   selection-background-color: {p['accent']};
                                   selection-color: {ON_ACCENT}; color: {p['fg']}; }}
    QComboBox QAbstractItemView::item {{ padding: 4px 8px; border-radius: {R['sm']}px;
                                         min-height: 20px; }}
    QComboBox QAbstractItemView::item:selected {{ background: {p['accent']};
                                                  color: {ON_ACCENT}; }}

    /* ---- the search beside Help. The field keeps an opaque fill like every other field; the list
       under it is a glass pane inside the window, so it is translucent on any display. ---- */
    QLineEdit#HelpSearchField {{ background: {FIELD_GREY}; border: 1px solid {p['border']};
                                 border-radius: {R['md'] - 2}px; padding: 2px 8px;
                                 font-size: {MENU}px; color: {FIELD_TEXT}; }}
    QLineEdit#HelpSearchField:focus {{ border-color: {p['accent']}; }}
    QFrame#HelpSearchResults {{ background: {POPUP_GLASS}; border: 1px solid {INNER_RIM};
                                border-radius: {R['lg']}px; }}
    QLabel#HelpSearchNote {{ background: transparent; color: {p['fg_muted']};
                             padding: 7px 12px; font-size: {MENU}px; }}
    QListWidget#HelpSearchResultList {{ background: transparent; border: none; padding: 4px;
                                        outline: 0; }}
    QListWidget#HelpSearchResultList::item {{ background: transparent; color: {p['fg']};
                                              border-radius: {R['sm'] + 1}px; }}
    QListWidget#HelpSearchResultList::item:selected {{ background: {p['accent_soft']};
                                                       color: {p['fg']}; }}
    *[helpSearchHit="true"] {{ border: 2px solid {p['accent']}; border-radius: {R['sm']}px; }}

    /* ---- Starplast's own windows, dressed as glass cards by `starplast.glass.dress`. The card
       paints the body; the window and its plain containers let it through; controls keep their
       own opaque surfaces, because a field you can see through is one you cannot read. ---- */
    QDialog[glass="true"], QMainWindow[glass="true"] {{ background: transparent; }}
    *[glass="true"] .QWidget, *[glass="true"] QStackedWidget,
    *[glass="true"] QScrollArea, *[glass="true"] QDialogButtonBox {{ background: transparent; }}
    *[glass="true"] QTabWidget::pane {{ background: {PANE_GLASS}; border: 1px solid {INNER_RIM};
                                        border-radius: {R['md']}px; }}
    *[glass="true"] QGroupBox {{ background: {PANE_GLASS}; border-color: {INNER_RIM}; }}
    *[glass="true"] QTabBar {{ background: transparent; }}
    QToolButton#GlassClose {{ background: transparent; border: none; color: {p['fg_muted']};
                              font-size: {BODY}px; padding: 2px 6px;
                              border-radius: {R['sm']}px; }}
    QToolButton#GlassClose:hover {{ background: {p['error']}; color: #ffffff; }}
    QLabel#GlassMessageText {{ padding: 4px 6px 8px 6px; }}
    """


#: The two-state colors of spacr's own switch, so the two programs read as one pair of tools.
#: Every editable field, in every theme. Dark grey rather than the surrounding panel's colour: a
#: field that takes its container's background stops looking like something you can type in, and on
#: a themed background it dissolves into whatever is drifting past behind it.
FIELD_GREY = "#2a2e35"

SWITCH_OFF = "#800080"
SWITCH_ON = "#008080"


class CheckList(QtWidgets.QListWidget):
    """A list of choices carried as TICKS, chosen by dragging across the rows and confirming with the
    right-click menu.

    Selection and choice are separated deliberately, and that separation is the whole point. Both
    category lists in this program used the SELECTION as the choice, which has two costs a reader
    pays constantly: a set assembled over several ctrl-clicks is destroyed by one ordinary click, and
    there is no way to keep a set while clicking somewhere else to look at something. Here a drag
    selects, the context menu turns that selection into ticks, and the ticks survive every later
    click.

    The one-item fast path is not lost: clicking an item's own box ticks it, which is Qt's behaviour
    for a checkable item and needs no menu. Space toggles the selected rows for the same reason.

    `checkedChanged` fires ONCE for a bulk operation rather than once per row. The lists it feeds
    trigger a redraw of 8,140 points, and ticking forty categories should redraw once.
    """

    checkedChanged = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._menu)
        self.itemChanged.connect(lambda _item: self.checkedChanged.emit())

    def add(self, label: str, value=None, **kw):
        """Append one checkable row. `value` is what `checked` returns for it, defaulting to `label`."""
        item = QtWidgets.QListWidgetItem(label)
        item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(QtCore.Qt.CheckState.Unchecked)
        item.setData(QtCore.Qt.ItemDataRole.UserRole, label if value is None else value)
        for key, val in kw.items():
            getattr(item, key)(val)
        self.addItem(item)
        return item

    def items(self) -> list:
        """Every row, in list order. `QListWidget` offers no such accessor of its own."""
        return [self.item(i) for i in range(self.count())]

    def checked(self) -> list:
        """The ticked values, in list order. Empty is a meaningful answer -- both callers read it as
        'all of them' -- so it is never padded into a full list here."""
        return [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in self.items()
                if i.checkState() == QtCore.Qt.CheckState.Checked]

    def set_checked(self, values, emit: bool = True) -> list:
        """Tick exactly these values and nothing else, as one change."""
        wanted = {str(v) for v in values}
        self.blockSignals(True)
        for item in self.items():
            state = (QtCore.Qt.CheckState.Checked
                     if str(item.data(QtCore.Qt.ItemDataRole.UserRole)) in wanted
                     else QtCore.Qt.CheckState.Unchecked)
            item.setCheckState(state)
        self.blockSignals(False)
        if emit:
            self.checkedChanged.emit()
        return self.checked()

    def _apply(self, items, state) -> None:
        """Set a tick state over many rows and report the change once."""
        if not items:
            return
        self.blockSignals(True)
        for item in items:
            item.setCheckState(state)
        self.blockSignals(False)
        self.checkedChanged.emit()

    def keyPressEvent(self, event):
        """Space toggles every selected row, which is what a list of tick boxes should do."""
        if event.key() in (QtCore.Qt.Key.Key_Space, QtCore.Qt.Key.Key_Select):
            chosen = self.selectedItems()
            if chosen:
                on = all(i.checkState() == QtCore.Qt.CheckState.Checked for i in chosen)
                self._apply(chosen, QtCore.Qt.CheckState.Unchecked if on
                            else QtCore.Qt.CheckState.Checked)
                return
        super().keyPressEvent(event)

    def _menu(self, pos) -> QtWidgets.QMenu:
        """The right-click menu: turn what is selected into ticks."""
        return self.build_menu(self.mapToGlobal(pos))

    def build_menu(self, at=None) -> QtWidgets.QMenu:
        """Construct the menu, and show it only when given somewhere to appear.

        Split so a test can read the actions without entering the modal loop `exec` starts -- the same
        split as `build_context_menu` and `build_preferences` in the window.
        """
        chosen = self.selectedItems()
        menu = QtWidgets.QMenu(self)
        act = menu.addAction(f"Check selected ({len(chosen)})")
        act.setEnabled(bool(chosen))
        act.triggered.connect(lambda: self._apply(chosen, QtCore.Qt.CheckState.Checked))
        act = menu.addAction(f"Uncheck selected ({len(chosen)})")
        act.setEnabled(bool(chosen))
        act.triggered.connect(lambda: self._apply(chosen, QtCore.Qt.CheckState.Unchecked))
        act = menu.addAction("Check only selected")
        act.setEnabled(bool(chosen))
        act.triggered.connect(lambda: self.set_checked(
            [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in chosen]))
        menu.addSeparator()
        menu.addAction("Check all").triggered.connect(
            lambda: self._apply(self.items(), QtCore.Qt.CheckState.Checked))
        menu.addAction("Clear all").triggered.connect(
            lambda: self._apply(self.items(), QtCore.Qt.CheckState.Unchecked))
        if at is not None:
            menu.exec(at)
        return menu


class CheckTree(QtWidgets.QTreeWidget):
    """A tree of tick boxes where ticking a CATEGORY ticks everything under it.

    The flat list this replaces had 96 rows, one per feature block, and at that length the list stops
    being a set of choices and becomes a scroll. What a reader actually wants to say is "all the
    transcription evidence", or "everything except the fitness screens" -- and both of those are
    statements about a LEVEL of the slot hierarchy, which already exists and already has the addresses.

    Three properties make it usable for the thing it is for, which is choosing what a map is built on
    before holding something out of it:

    * **Ticking a group ticks its descendants**, and a group whose children are only partly ticked
      shows partially rather than lying in either direction.
    * **`checked()` returns LEAVES only.** A group is a way of saying something about its leaves, not a
      thing that can itself feed a map.
    * **"Everything except this"** is one action, because that is the hold-out question -- build the
      map on all the evidence but this, then ask whether this comes back -- and doing it by hand on 96
      boxes is how a reader ends up holding out something they did not mean to.
    """

    checkedChanged = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setHeaderHidden(True)
        self.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._menu)
        self.itemChanged.connect(self._changed)
        self._quiet = False

    # ------------------------------------------------------------------ building
    def build(self, tree: dict, label=None, enabled=None) -> int:
        """Fill from a nested dict: keys are names, empty values are leaves. Returns the leaf count."""
        self._quiet = True
        self.clear()
        made = self._add(self.invisibleRootItem(), tree, label or (lambda k: k), enabled)
        self._quiet = False
        self.expandToDepth(0)
        return made

    def _add(self, parent, node, label, enabled) -> int:
        made = 0
        for key, child in (node or {}).items():
            item = QtWidgets.QTreeWidgetItem([label(key) if not child else key])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, key)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.CheckState.Unchecked)
            parent.addChild(item)
            if child:
                made += self._add(item, child, label, enabled)
            else:
                made += 1
                if enabled is not None and not enabled(key):
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEnabled)
        return made

    # ------------------------------------------------------------------ state
    def leaves(self, under=None) -> list:
        """Every leaf row, or every leaf beneath one node. Leaves are the blocks; groups are addresses."""
        out, stack = [], [under or self.invisibleRootItem()]
        while stack:
            node = stack.pop()
            for i in range(node.childCount()):
                child = node.child(i)
                (stack if child.childCount() else out).append(child)
        return out

    def checked(self) -> list:
        """The ticked LEAVES, in tree order. A group is not itself an answer."""
        return [i.data(0, QtCore.Qt.ItemDataRole.UserRole) for i in self.leaves()
                if i.checkState(0) == QtCore.Qt.CheckState.Checked]

    def set_checked(self, values, emit: bool = True) -> list:
        """Tick exactly these leaves and nothing else, as one change, skipping any that are disabled."""
        wanted = {str(v) for v in values}
        self._quiet = True
        for item in self.leaves():
            usable = bool(item.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled)
            item.setCheckState(0, QtCore.Qt.CheckState.Checked
                               if usable and str(item.data(0, QtCore.Qt.ItemDataRole.UserRole)) in wanted
                               else QtCore.Qt.CheckState.Unchecked)
        self._refresh_parents()
        self._quiet = False
        if emit:
            self.checkedChanged.emit()
        return self.checked()

    def _changed(self, item, _column=0):
        """Push a group's new state down to its leaves, then recompute every group."""
        if self._quiet:
            return
        self._quiet = True
        if item.childCount():
            state = item.checkState(0)
            for leaf in self.leaves(item):
                if leaf.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled:
                    leaf.setCheckState(0, state)
        self._refresh_parents()
        self._quiet = False
        self.checkedChanged.emit()

    def _refresh_parents(self):
        """A group is checked when all its usable leaves are, partial when some are."""
        def walk(node):
            for i in range(node.childCount()):
                child = node.child(i)
                if not child.childCount():
                    continue
                walk(child)
                usable = [x for x in self.leaves(child)
                          if x.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled]
                on = sum(1 for x in usable
                         if x.checkState(0) == QtCore.Qt.CheckState.Checked)
                child.setCheckState(0, QtCore.Qt.CheckState.Checked if usable and on == len(usable)
                                    else QtCore.Qt.CheckState.PartiallyChecked if on
                                    else QtCore.Qt.CheckState.Unchecked)
        walk(self.invisibleRootItem())

    # ------------------------------------------------------------------ the menu
    def _menu(self, pos):
        return self.build_menu(self.itemAt(pos), self.viewport().mapToGlobal(pos))

    def build_menu(self, item=None, at=None) -> QtWidgets.QMenu:
        """Shown only when given somewhere to appear, so a test can read it without a modal loop."""
        menu = QtWidgets.QMenu(self)
        name = item.text(0) if item is not None else ""
        act = menu.addAction(f"Check everything under “{name}”")
        act.setEnabled(item is not None)
        act.triggered.connect(lambda: self._apply(self.leaves(item), True))
        act = menu.addAction(f"Uncheck everything under “{name}”")
        act.setEnabled(item is not None)
        act.triggered.connect(lambda: self._apply(self.leaves(item), False))
        menu.addSeparator()
        # The hold-out, in one click. Build the map on all the evidence EXCEPT this category, then ask
        # whether this category comes back -- which is the only question a held-out feature answers.
        act = menu.addAction(f"Check everything EXCEPT “{name}”  (hold this out)")
        act.setEnabled(item is not None)
        act.triggered.connect(lambda: self._except(item))
        menu.addSeparator()
        menu.addAction("Check all").triggered.connect(lambda: self._apply(self.leaves(), True))
        menu.addAction("Clear all").triggered.connect(lambda: self._apply(self.leaves(), False))
        if at is not None:
            menu.exec(at)
        return menu

    def _apply(self, items, on: bool):
        self._quiet = True
        for item in items:
            if item.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled:
                item.setCheckState(0, QtCore.Qt.CheckState.Checked if on
                                   else QtCore.Qt.CheckState.Unchecked)
        self._refresh_parents()
        self._quiet = False
        self.checkedChanged.emit()

    def _except(self, item):
        """Everything ticked but this branch -- the hold-out."""
        held = {id(x) for x in self.leaves(item)}
        self._quiet = True
        for leaf in self.leaves():
            if leaf.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled:
                leaf.setCheckState(0, QtCore.Qt.CheckState.Unchecked if id(leaf) in held
                                   else QtCore.Qt.CheckState.Checked)
        self._refresh_parents()
        self._quiet = False
        self.checkedChanged.emit()
        return self.checked()


class Switch(QtWidgets.QWidget):
    """A boolean slider, the shape and colors of `spacr.gui_elements.spacrSwitch`.

    A checkbox would have done the job. This exists because the user runs both programs and a
    setting that looks like a setting in one of them should look like one in the other: 40x20 track,
    a 12-pixel knob that slides, purple for off and teal for on, caption on the left.
    """

    toggled = QtCore.pyqtSignal(bool)

    def __init__(self, text: str = "", checked: bool = False, parent=None):
        super().__init__(parent)
        self._on = bool(checked)
        self._x = 24.0 if self._on else 4.0
        self.text = text
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setMinimumHeight(24)
        self._anim = QtCore.QPropertyAnimation(self, b"knob", self)
        self._anim.setDuration(120)

    # The knob position is a Qt property so QPropertyAnimation can drive it; the slide is the whole
    # point of copying this control rather than using a checkbox.
    def _get_knob(self) -> float:
        return self._x

    def _set_knob(self, value: float) -> None:
        self._x = float(value)
        self.update()

    knob = QtCore.pyqtProperty(float, fget=_get_knob, fset=_set_knob)

    def isChecked(self) -> bool:
        """The current state."""
        return self._on

    def setChecked(self, on: bool) -> None:
        """Set the state, sliding the knob, without emitting."""
        self._on = bool(on)
        self._anim.stop()
        self._anim.setStartValue(self._x)
        self._anim.setEndValue(24.0 if self._on else 4.0)
        self._anim.start()
        self.update()

    def mouseReleaseEvent(self, ev):
        """Click anywhere on the control -- the track or the caption -- to toggle it."""
        if ev.button() == QtCore.Qt.MouseButton.LeftButton:
            self.setChecked(not self._on)
            self.toggled.emit(self._on)
        super().mouseReleaseEvent(ev)

    def sizeHint(self):
        """Wide enough for the caption plus the 40-pixel track and its margin."""
        w = QtGui.QFontMetrics(self.font()).horizontalAdvance(self.text)
        return QtCore.QSize(w + 60, 24)

    def paintEvent(self, _ev):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        fm = QtGui.QFontMetrics(self.font())
        if self.text:
            p.setPen(self.palette().color(QtGui.QPalette.ColorRole.WindowText))
            p.drawText(0, 0, self.width() - 50, self.height(),
                       int(QtCore.Qt.AlignmentFlag.AlignVCenter
                           | QtCore.Qt.AlignmentFlag.AlignLeft), self.text)
        left = self.width() - 42
        p.setPen(QtCore.Qt.PenStyle.NoPen)
        p.setBrush(QtGui.QColor("#ffffff"))
        p.drawRoundedRect(QtCore.QRectF(left, 2, 36, 16), 8, 8)
        p.setBrush(QtGui.QColor(SWITCH_ON if self._on else SWITCH_OFF))
        p.drawEllipse(QtCore.QRectF(left + self._x - 2, 4, 12, 12))
        p.end()


#: How wide a tooltip block is, in CHARACTERS. A tooltip is a paragraph of explanation in this
#: program -- why a control exists and what choosing badly costs -- and 64 characters is about the
#: width of a column of prose, which is what makes a block of it readable.
TOOLTIP_WIDTH = 64


def tip(text: str, width: int = 64) -> str:
    """A tooltip as one block of even lines, wrapped once and left that way.

    Qt lays a plain tooltip out on a single line, so any real explanation becomes a strip wider than
    the screen. Inserting line breaks alone does not fix it: Qt re-wraps rich text at a width of its
    own choosing and strands two words on a row. **A fixed table width does not constrain it
    either** -- that was tried here and measured at 1,588 pixels, because Qt only wraps when it is
    given an explicit text width and a tooltip sets its own.

    `white-space: pre` is what stops the second wrap: 1,188 pixels to 347 on a typical tooltip, with
    the breaks where they were put. Lines are padded to equal length so the block is a rectangle
    rather than a ragged edge; true justification is not in Qt's rich-text subset.

    `width` is in characters, not pixels, because that is what the wrapping is done in.
    """
    import textwrap
    from html import escape
    blocks = []
    for para in [" ".join(p.split()) for p in str(text).split("\n\n") if p.strip()]:
        lines = textwrap.wrap(para, width) or [""]
        longest = max(len(x) for x in lines)
        # Padded with non-breaking spaces, which Qt keeps; ordinary trailing spaces are dropped.
        blocks.append("\n".join(escape(x) + "&#160;" * (longest - len(x)) for x in lines))
    return '<div style="white-space:pre">' + "\n\n".join(blocks) + "</div>"


def wrap_tooltips(root, width: int = TOOLTIP_WIDTH) -> int:
    """Rewrite every tooltip under `root` as a wrapped block, and move it onto its label.

    Two things, because they are the same complaint: a tooltip belongs on the SETTING, not on the
    box you type in. A form row is a label and a field, and hovering the label -- which is what
    names the thing and what the eye goes to first -- said nothing at all. Returns how many were
    rewritten, so a test can tell this ran.
    """
    from PyQt6 import QtWidgets
    done = 0
    # By TYPE, not every QWidget under the root. Walking all of them crashed the interpreter -- this
    # window contains a pyqtgraph GL view, and asking its internals for a tooltip from Python is a
    # segfault rather than an exception. Settings live on these six classes; nothing else in this
    # program carries an explanation worth wrapping.
    kinds = (QtWidgets.QAbstractButton, QtWidgets.QComboBox, QtWidgets.QAbstractSpinBox,
             QtWidgets.QLineEdit, QtWidgets.QLabel, QtWidgets.QAbstractSlider, Switch)
    found = []
    for kind in kinds:
        found.extend(root.findChildren(kind))
    for w in found:
        try:
            text = w.toolTip()
            if not text or "white-space:pre" in text:
                continue
            w.setToolTip(tip(text, width))
            done += 1
            parent = w.parentWidget()
            form = parent.layout() if parent is not None else None
            if isinstance(form, QtWidgets.QFormLayout):
                label = form.labelForField(w)
                if label is not None and not label.toolTip():
                    label.setToolTip(w.toolTip())
        except RuntimeError:
            continue
    return done
