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


def rgbf(hex_colour: str, alpha: float = 1.0) -> tuple:
    """Hex to the 0-1 RGBA tuple pyqtgraph's GL items want."""
    h = hex_colour.lstrip("#")
    return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, alpha)


def background(theme: str) -> tuple:
    """GL viewport background — the page color, not a literal, so themes stay coherent."""
    return rgbf(palette_for(theme)["bg"])


def unknown_colour(theme: str) -> tuple:
    """The gray for 'unknown, not absent'. Derived from fg_dim so it reads on either ground."""
    return rgbf(palette_for(theme)["fg_dim"])


# --------------------------------------------------------------------------- colour maps
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
# cueing reads as a flat disc of colour: near and far points are equally bright, so the eye has nothing
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


def categorical_colours(n: int, theme: str = "dark", cmap: str | None = None) -> list:
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


# --------------------------------------------------------------------------- stylesheet
def stylesheet(theme: str = "dark") -> str:
    """Qt stylesheet for one theme. Every color comes from the palette, never a literal."""
    p = palette_for(theme)
    return f"""
    QWidget {{ background: {p['bg']}; color: {p['fg']};
               font-size: 13px; selection-background-color: {p['accent']};
               selection-color: {p['bg'] if is_dark(theme) else p['page']}; }}
    QMainWindow, QDialog {{ background: {p['bg']}; }}
    QDockWidget {{ background: {p['surface']}; color: {p['fg']}; titlebar-close-icon: none; }}
    QDockWidget::title {{ background: {p['surface_alt']}; padding: 7px 10px;
                          border-bottom: 1px solid {p['border']};
                          font-weight: 600; color: {p['fg_muted']}; }}
    QGroupBox {{ border: 1px solid {p['border']}; border-radius: 6px; margin-top: 16px;
                 padding-top: 10px; background: {p['surface']}; }}
    QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 5px;
                        color: {p['fg_muted']}; font-size: 11px;
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
    QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {{
        background: {p['surface_hi']}; border: 1px solid {p['border']};
        border-radius: 5px; padding: 5px 8px; color: {p['fg']}; }}
    QComboBox:focus, QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
        border-color: {p['accent']}; }}
    QComboBox QAbstractItemView {{ background: {p['surface']}; border: 1px solid {p['border']};
                                   selection-background-color: {p['accent_soft']};
                                   color: {p['fg']}; }}
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
        font-size: 11px; text-transform: uppercase; letter-spacing: .6px; }}
    QTextBrowser {{ background: {p['page']}; border: 1px solid {p['border']};
                    border-radius: 6px; padding: 4px; }}
    QProgressBar {{ background: {p['surface_hi']}; border: none; border-radius: 4px;
                    height: 6px; text-align: center; color: transparent; }}
    QProgressBar::chunk {{ background: {p['accent']}; border-radius: 4px; }}
    QStatusBar {{ background: {p['surface']}; color: {p['fg_muted']};
                  border-top: 1px solid {p['border']}; }}
    QScrollBar:vertical {{ background: transparent; width: 10px; margin: 0; }}
    QScrollBar::handle:vertical {{ background: {p['border']}; border-radius: 5px;
                                   min-height: 26px; }}
    QScrollBar::handle:vertical:hover {{ background: {p['fg_dim']}; }}
    QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
    QSplitter::handle {{ background: {p['border_soft']}; }}
    QToolTip {{ background: {p['surface_hi']}; color: {p['fg']};
                border: 1px solid {p['border']}; padding: 5px; }}
    """
