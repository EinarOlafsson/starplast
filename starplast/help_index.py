"""Everything in Starplast that can be reached by name, and how a query ranks it.

The search field beside the Help menu is only as good as this index: a user who knows what a thing is
CALLED -- a strategy, a slot, a dataset, a menu command, a tab, a setting, a section of the guide --
should reach it without knowing which menu, panel or document it was filed under. The ranking is
spaCR's (`spacr.qt.help_index`), so the two programs answer a query the same way.

THE INDEX IS NEVER A HAND LIST. Every row is derived from the thing that already decides the answer:

==============  ===========================================================================
kind            where the rows come from
==============  ===========================================================================
``action``      every command in the window's menu bar, with its tooltip and status tip
``panel``       every dock, every tab of the analysis, strategy and preference panels, the
                guided workflows and the slot tree
``setting``     every labelled control in the analysis panel and in Preferences, and every
                display choice the map's right-click menu offers
``strategy``    `starplast.strategies.catalog` -- key, number, name, method, family, question
``slot``        `starplast.slots.all_slots` -- the information categories of both organisms
``dataset``     `starplast.datasets.REGISTRY`
``doc``         the headings of ``docs/*.md`` in a source checkout, or the published pages
``tutorial``    the pages and sections of ``docs/tutorial``
==============  ===========================================================================

A command that cannot be found is then a fact about a registry, not about whether somebody remembered
to add a row.

HEADLESS. Importing this module imports no Qt. The providers that read the window (`action`,
`panel`, `setting`) are handed the window and walk its widgets, which works offscreen; the rest read
files and registries and need no application at all, so the ranking and the static rows are tested
without one. What a result DOES when it is chosen is `starplast.help_search`'s business.
"""
from __future__ import annotations

import html
import os
import re
from dataclasses import dataclass, field
from functools import lru_cache

#: The kinds, in the order a tie between two equally good matches is broken: the coarse things a
#: user reaches for first ahead of the ten thousand fine ones.
KIND_ORDER = ("action", "panel", "strategy", "setting", "slot", "dataset", "doc", "tutorial")

#: What each kind is called on a result row.
KIND_LABEL = {"action": "command", "panel": "panel", "strategy": "strategy", "setting": "setting",
              "slot": "slot", "dataset": "dataset", "doc": "guide", "tutorial": "tutorial"}

#: How much a kind is worth when two rows score the same. Small enough that it never outranks a
#: better textual match.
_KIND_BONUS = {kind: (len(KIND_ORDER) - i) * 0.01 for i, kind in enumerate(KIND_ORDER)}

#: The repository's docs directory, when this is a source checkout. An installed wheel has none, and
#: the published pages stand in for it.
DOCS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")

#: Where the documentation is published, and the pages it publishes (scripts/build_docs.py).
SITE = "https://einarolafsson.github.io/starplast/"
REPO = "https://github.com/EinarOlafsson/starplast/blob/main/"
PUBLISHED = {"index": "Starplast documentation", "guide": "User guide", "API": "Python API",
             "datasets": "Dataset catalogue", "structures": "Protein structures",
             "protein-sequences": "Protein sequences", "releases": "Releases",
             "workflows": "Guided workflows", "benchmark-0.43": "Benchmark 0.43",
             "repository-review": "Repository review", "scientific-roadmap": "Scientific roadmap",
             "calibration": "Strategy calibration", "graphspace": "Integrated neighbour space",
             "strategies": "Strategy catalogue"}


@dataclass(frozen=True)
class HelpEntry:
    """One reachable thing, and what it takes to get back to it.

    :ivar kind: which provider produced it, and so which opener acts on it; one of `KIND_ORDER`.
    :ivar title: what a user would type -- a command's label, a strategy's name, a slot's name.
    :ivar subtitle: where it lives, as a path: "View ▸ Color by", "Strategies ▸ Networks".
    :ivar description: what it does, in words. Matched too, which is what lets "density" find the
        additive overlap mode and "HDBSCAN" find every strategy that clusters with it.
    :ivar payload: what the opener needs and nothing more, as sorted (name, value) pairs so the
        entry stays hashable. Read it with `data`.
    """

    kind: str
    title: str
    subtitle: str = ""
    description: str = ""
    payload: tuple = ()

    @property
    def data(self) -> dict:
        """The payload as a dict."""
        return dict(self.payload)

    @property
    def haystack(self) -> str:
        """Title, subtitle and description, lowercased, for matching."""
        return f"{self.title}\n{self.subtitle}\n{self.description}".lower()


def entry(kind: str, title: str, subtitle: str = "", description: str = "", **payload) -> HelpEntry:
    """Build an entry with a payload given as keywords."""
    return HelpEntry(kind, str(title).strip(), str(subtitle).strip(),
                     " ".join(str(description).split()),
                     tuple(sorted((k, str(v)) for k, v in payload.items())))


def plain(text: str) -> str:
    """Tooltip text without its markup: the wrapped tooltips are HTML with padded lines."""
    text = re.sub(r"<[^>]+>", " ", str(text or ""))
    return " ".join(html.unescape(text).replace(" ", " ").split())


def clean_label(text: str) -> str:
    """A menu label as a person reads it: no mnemonic ampersand, no trailing ellipsis."""
    text = str(text or "").replace("&&", "\0").replace("&", "").replace("\0", "&")
    return text.rstrip("…").rstrip(".").strip() if text.endswith(("…", "...")) else text.strip()


# --------------------------------------------------------------------------- ranking
#: How far apart the letters of a fuzzy match may be spread, as a multiple of the term's length,
#: before it stops being a match. spaCR measured it: without a ceiling a long enough name contains
#: almost any short sequence of letters.
_FUZZY_SPAN_LIMIT = 3


def _subsequence_span(term: str, text: str):
    """How wide a window of `text` the letters of `term` fit in, in order; None if they do not."""
    if not term:
        return 0
    start, position = -1, 0
    for letter in term:
        found = text.find(letter, position)
        if found < 0:
            return None
        if start < 0:
            start = found
        position = found + 1
    span = position - start
    return None if span > _FUZZY_SPAN_LIMIT * len(term) else span


@lru_cache(maxsize=None)
def _name_parts(title: str) -> tuple:
    """A title split the three ways a query addresses it: whole, last dotted part, and words."""
    whole = title.lower()
    tail = whole.rpartition(".")[2]
    words = tuple(w for w in re.split(r"[^a-z0-9]+", whole) if w)
    return whole, tail, words


def term_score(term: str, e: HelpEntry):
    """How well one query term matches one entry, or None for no match.

    Five bands, and the order between them is the whole ranking: the name typed in full beats a word
    of it, which beats the start of a word, which beats a fuzzy spelling, which beats a word from the
    description. A description hit must never outrank a name hit, or typing a command's own label
    buries it under every tooltip that mentions it.
    """
    whole, tail, words = _name_parts(e.title)
    if whole == term:
        return 100.0
    if tail == term or term in words:
        return 80.0
    if whole.startswith(term):
        return 60.0 + 10.0 * len(term) / max(len(whole), 1)
    if tail.startswith(term) or any(w.startswith(term) for w in words):
        return 50.0 + 10.0 * len(term) / max(len(tail), 1)
    if term in whole:
        return 40.0
    span = _subsequence_span(term, tail) if len(term) >= 3 else None
    if span is not None:
        return 20.0 + 10.0 * len(term) / max(span, 1)
    if term in f"{e.subtitle}\n{e.description}".lower():
        return 10.0
    return None


def score(e: HelpEntry, terms) -> float | None:
    """The entry's score for a whole query, or None when any term misses: terms are ANDed, so a
    second word narrows the list rather than widening it."""
    total = 0.0
    for term in terms:
        one = term_score(term, e)
        if one is None:
            return None
        total += one
    return total + _KIND_BONUS.get(e.kind, 0.0)


#: How many rows of one kind a result list may hold. There are 290 slots and 150 datasets against a
#: few dozen panels; without a cap the list is one kind, and the point is that one query answers
#: with several.
PER_KIND_LIMIT = 6


def search(index, query: str, limit: int = 24, per_kind: int | None = PER_KIND_LIMIT) -> list:
    """The best matches for `query`, best first; ties broken by kind and then by title."""
    terms = str(query or "").lower().split()
    if not terms:
        return []
    scored = []
    for e in index:
        value = score(e, terms)
        if value is None:
            continue
        rank = KIND_ORDER.index(e.kind) if e.kind in KIND_ORDER else len(KIND_ORDER)
        scored.append((-value, rank, e.title.lower(), e))
    scored.sort(key=lambda row: row[:3])
    out, per = [], {}
    for _v, _r, _t, e in scored:
        if per_kind is not None and per.get(e.kind, 0) >= per_kind:
            continue
        per[e.kind] = per.get(e.kind, 0) + 1
        out.append(e)
        if len(out) >= limit:
            break
    return out


# --------------------------------------------------------------------------- static providers
@lru_cache(maxsize=1)
def strategy_entries() -> tuple:
    """One row per strategy: its number and name, found by key, method, family and question."""
    from . import strategies as S
    return tuple(entry("strategy", f"{s.number:02d} · {s.name}", f"Strategies ▸ {s.family}",
                       " ".join([s.question, plain(s.tooltip), s.method, s.key, s.family]),
                       key=s.key)
                 for s in S.catalog())


@lru_cache(maxsize=1)
def slot_entries() -> tuple:
    """One row per slot of both organisms, found by its axis, context and three addresses."""
    from . import slots as SL
    out = []
    for s in SL.all_slots():
        where = " ▸ ".join(s.evidence_path[:3]) if s.evidence_path else s.axis
        out.append(entry("slot", s.name, f"Slot tree ▸ {s.organism} ▸ {where}",
                         " ".join([s.axis, s.context, s.unit, s.target_family.replace("_", " "),
                                   " ".join(s.biology_path), " ".join(s.context_path), s.key]),
                         organism=s.organism, name=s.name))
    return tuple(out)


@lru_cache(maxsize=1)
def dataset_entries() -> tuple:
    """One row per registered dataset, found by what it provides, its citation and its key."""
    from . import datasets as D
    out = []
    for d in D.REGISTRY:
        text = " ".join(str(x) for x in (d.provides, d.kind, d.level, d.coverage, d.citation,
                                         d.accession or "", d.key) if x)
        out.append(entry("dataset", d.name, f"Datasets ▸ {d.level.replace('_', ' ')}", text,
                         key=d.key))
    return tuple(out)


def _slug(heading: str) -> str:
    """The anchor GitHub and the published site give a heading."""
    s = re.sub(r"[^\w\- ]", "", heading.strip().lower())
    return re.sub(r"\s+", "-", s)


def _heading_text(raw: str) -> str:
    """A markdown heading without its markup: links, code ticks and emphasis."""
    raw = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", raw)
    return re.sub(r"[`*_]", "", raw).strip()


@lru_cache(maxsize=4)
def doc_entries(docs_dir: str = DOCS_DIR) -> tuple:
    """Every heading of every page in the docs directory, or the published pages without one.

    A heading is the unit a reader asks for -- "leakage audit", "held-out recovery" -- and a page is
    too coarse to land on. Code blocks are skipped: a `# comment` inside one is not a heading.
    """
    out = []
    if not os.path.isdir(docs_dir):
        for page, title in PUBLISHED.items():
            out.append(entry("doc", title, "Documentation (online)", page,
                             page=page, url=f"{SITE}{page}.html"))
        return tuple(out)
    for name in sorted(os.listdir(docs_dir)):
        if not name.endswith(".md"):
            continue
        path = os.path.join(docs_dir, name)
        page = name[:-3]
        try:
            lines = open(path, encoding="utf8").read().splitlines()
        except OSError:
            continue
        title, fenced = PUBLISHED.get(page, page.replace("-", " ").capitalize()), False
        for line in lines:
            if line.lstrip().startswith("```"):
                fenced = not fenced
                continue
            m = None if fenced else re.match(r"^(#{1,4})\s+(.+?)\s*#*\s*$", line)
            if not m:
                continue
            text = _heading_text(m.group(2))
            if len(m.group(1)) == 1 and title == page.replace("-", " ").capitalize():
                title = text
            if not text:
                continue
            url = (f"{SITE}{page}.html#{_slug(text)}" if page in PUBLISHED
                   else f"{REPO}docs/{name}#{_slug(text)}")
            out.append(entry("doc", text, f"Guide ▸ {title}", f"{page} {name}",
                             path=path, heading=text, level=len(m.group(1)), url=url))
    return tuple(out)


@lru_cache(maxsize=4)
def tutorial_entries(docs_dir: str = DOCS_DIR) -> tuple:
    """Each tutorial page, and each numbered section of the complete guide."""
    root = os.path.join(docs_dir, "tutorial")
    if not os.path.isdir(root):
        return ()
    out = []
    for name in sorted(os.listdir(root)):
        if not name.endswith(".html"):
            continue
        path = os.path.join(root, name)
        try:
            text = open(path, encoding="utf8").read()
        except OSError:
            continue
        m = re.search(r"<title>([^<]+)</title>", text)
        title = html.unescape(m.group(1)).strip() if m else name
        kind = "notebook" if "notebook" in name else "GUI walkthrough" if "_gui" in name else "page"
        out.append(entry("tutorial", title, f"Tutorials ▸ {kind}", name, path=path, anchor=""))
        for anchor, heading in re.findall(r"<h2 id='([^']+)'>([^<]+)</h2>", text):
            out.append(entry("tutorial", html.unescape(heading).strip(), f"Tutorials ▸ {title}",
                             name, path=path, anchor=anchor))
    return tuple(out)


def static_entries() -> list:
    """Every row that does not need a window: strategies, slots, datasets, docs and tutorials.

    Each provider is guarded on its own: a search field that loses its strategies because the slot
    catalogue could not be read is worse than one that quietly offers fewer kinds.
    """
    out = []
    for provider in (strategy_entries, slot_entries, dataset_entries, doc_entries,
                     tutorial_entries):
        try:
            out.extend(provider())
        except Exception:
            continue
    return out


# --------------------------------------------------------------------------- window providers
def _menu_rows(menu, path: tuple, skip: set, out: list) -> None:
    """Walk one menu, depth first, adding a row for every command and every submenu.

    A command placed in two menus -- the slot tree is in both View and Tools -- is one command, and
    is listed once, at the first place it was found.
    """
    for act in menu.actions():
        if act.isSeparator() or not act.text().strip() or act in skip:
            continue
        # The action itself, not id() of it: a Python wrapper can be collected and its id handed to
        # the next one, which silently dropped "Color by" from the index.
        skip.add(act)
        label = clean_label(act.text())
        where = " ▸ ".join(path)
        tip = plain(act.statusTip()) or plain(act.toolTip())
        tip = "" if tip == act.text() or tip == label else tip
        if act.menu() is not None:
            out.append(entry("action", label, where, f"menu {tip}",
                             path=" ▸ ".join(path + (label,)), menu="1"))
            _menu_rows(act.menu(), path + (label,), skip, out)
            continue
        shortcut = act.shortcut().toString() if not act.shortcut().isEmpty() else ""
        out.append(entry("action", label, where, f"{tip} {shortcut}",
                         path=" ▸ ".join(path + (label,)), shortcut=shortcut))


def action_entries(window) -> list:
    """Every command in the menu bar, including submenus, addressed by its menu path.

    The panels' own show/hide toggles are left to `panel_entries`, which raises a panel rather than
    toggling it -- a search for "console" that hides the console because it was already open would
    be a search that does the opposite of what was asked.
    """
    skip = {d.toggleViewAction() for d in _docks(window)}
    out = []
    for top in window.menuBar().actions():
        if top.menu() is None:
            continue
        _menu_rows(top.menu(), (clean_label(top.text()),), skip, out)
    return out


def _docks(window) -> list:
    from PyQt6 import QtWidgets
    return window.findChildren(QtWidgets.QDockWidget)


#: What each dock is for, so a search for its purpose finds it and the result row says why.
DOCK_HELP = {
    "find": "Find a gene by accession or product, colour the map, and tick the classes to show.",
    "evidence": "Everything recorded about the selected gene, with links to its source records.",
    "analysis": "Data, maps, clusters, inference, search, validation, discovery and questions.",
    "strategies": "Named ways of inferring something, each with a guide, settings and a self-test.",
    "gallery": "Thumbnails of every map a walk or search built; click one to show it.",
    "console": "Everything the program has printed, including errors from background jobs.",
    "jobs": "Every job started this session, with stop buttons and memory use.",
    "assistant": "Ask questions about the map in front of you.",
}


def panel_entries(window) -> list:
    """Every dock, and every tab inside the analysis, strategy, preference and workflow panels."""
    out = []
    for d in _docks(window):
        name = d.windowTitle()
        out.append(entry("panel", name.capitalize(), "Panel", DOCK_HELP.get(name, ""),
                         dock=name))
    panel = getattr(window, "panel", None)
    tabs = getattr(panel, "tabs", None)
    if tabs is not None:
        for i in range(tabs.count()):
            label = tabs.tabText(i)
            out.append(entry("panel", label.split("·")[-1].strip(), "Analysis ▸ tab",
                             f"analysis panel {label} {plain(tabs.tabToolTip(i))}",
                             dock="analysis", tab=label))
    strat = getattr(window, "strategy_panel", None)
    if strat is not None:
        for i in range(strat.tabs.count()):
            label = strat.tabs.tabText(i)
            out.append(entry("panel", f"Strategy {label.lower()}", "Strategies ▸ tab",
                             f"strategies panel {label}", dock="strategies", tab=label))
    for i, label in enumerate(("Explore a gene", "Predict a trait", "Compare a screen")):
        out.append(entry("panel", label, "Tools ▸ Guided workflows",
                         "guided workflow with source evidence, held-out evaluation and exports",
                         workflow=i))
    for label in ("Appearance", "Display", "Window"):
        out.append(entry("panel", f"Preferences: {label}", "File ▸ Preferences",
                         f"preferences settings {label}", preferences=label))
    out.append(entry("panel", "Slot tree", "Tools ▸ Slot tree",
                     "every information category, its address in three hierarchies, and what "
                     "fills it; the empty ones are coloured", slot_tree="1"))
    return out


def _form_rows(root, prefix: str, extra: dict) -> list:
    """Every labelled row of every form layout under `root`, as setting entries."""
    from PyQt6 import QtWidgets
    out, seen = [], set()
    for form in root.findChildren(QtWidgets.QFormLayout):
        for row in range(form.rowCount()):
            label_item = form.itemAt(row, QtWidgets.QFormLayout.ItemRole.LabelRole)
            field_item = form.itemAt(row, QtWidgets.QFormLayout.ItemRole.FieldRole)
            label = label_item.widget() if label_item is not None else None
            fieldw = field_item.widget() if field_item is not None else None
            if not isinstance(label, QtWidgets.QLabel) or fieldw is None:
                continue
            text = clean_label(plain(label.text())).rstrip(":")
            if not text or (text, prefix) in seen:
                continue
            seen.add((text, prefix))
            tip = plain(fieldw.toolTip() or label.toolTip())
            out.append(entry("setting", text, prefix, tip, label=text, **extra))
    return out


def setting_entries(window) -> list:
    """Every labelled control in the analysis panel and Preferences, and every display choice.

    A control is addressed by its tab and its label, and found again from them when it is opened, so
    the index holds no reference to a widget that may since have been rebuilt.
    """
    out = []
    panel = getattr(window, "panel", None)
    tabs = getattr(panel, "tabs", None)
    if tabs is not None:
        for i in range(tabs.count()):
            label = tabs.tabText(i)
            out.extend(_form_rows(tabs.widget(i), f"Analysis ▸ {label}",
                                  {"dock": "analysis", "tab": label}))
    prefs = window.preferences_dialog() if hasattr(window, "preferences_dialog") else None
    from PyQt6 import QtWidgets
    ptabs = prefs.findChild(QtWidgets.QTabWidget) if prefs is not None else None
    if ptabs is not None:
        for i in range(ptabs.count()):
            label = ptabs.tabText(i)
            out.extend(_form_rows(ptabs.widget(i), f"Preferences ▸ {label}",
                                  {"preferences": label}))
    if hasattr(window, "display_choices"):
        from .app import DISPLAY_HELP
        for label, options, _current, _apply in window.display_choices():
            for option in options:
                out.append(entry("setting", f"{label}: {option}", "Map ▸ right-click ▸ Display",
                                 DISPLAY_HELP.get(label, ""), display=label, option=option))
    return out


def window_entries(window) -> list:
    """Every row that is read off a live window: commands, panels and settings."""
    out = []
    for provider in (action_entries, panel_entries, setting_entries):
        try:
            out.extend(provider(window))
        except Exception:
            continue
    return out


def build_index(window=None) -> list:
    """The whole index: the window's rows when a window is given, and every static row.

    Rows that would read identically in the list -- the same heading twice on one page, one per
    organism -- are kept once, since two rows nobody can tell apart are one row and a wasted seat.
    """
    rows = (window_entries(window) if window is not None else []) + static_entries()
    seen, out = set(), []
    for e in rows:
        key = (e.kind, e.title, e.subtitle)
        if key not in seen:
            seen.add(key)
            out.append(e)
    return out
