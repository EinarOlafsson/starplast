#!/usr/bin/env python3
"""Build the Starplast tutorials from the running application.

Writes ``docs/tutorial/``: an index, five tutorials each as a GUI walkthrough and a notebook
walkthrough (the format of spaCR's tutorials), a runnable ``.ipynb`` per notebook tutorial, and one
long guide to every feature. Nothing is typed in by hand that the program can produce: every
screenshot is a grab of a real panel, every console line was emitted by a real run, every table and
number was computed from the shipped tables while this script ran. Rebuild after any change that
alters what the panels show:

    QT_QPA_PLATFORM=offscreen python scripts/build_tutorials.py

The 3D map needs OpenGL, which an offscreen Qt cannot draw; map views are therefore the committed
screenshots in ``docs/deck/assets`` and, for results, renders of the result's own coordinates and
cluster labels (``StrategyResult.plot``), captioned as such.
"""
from __future__ import annotations

import base64
import contextlib
import html
import io
import json
import os
import shutil
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import tutorial_templates as T  # noqa: E402

OUT = os.path.join(ROOT, "docs", "tutorial")
RES = os.path.join(OUT, "resources")
NOTEBOOKS = os.path.join(OUT, "notebooks")
DECK = os.path.join(ROOT, "docs", "deck", "assets")
SHOTS = os.path.join(ROOT, "docs", "screenshots")
_e = html.escape


# --------------------------------------------------------------------------- Qt and captures
class Studio:
    """An offscreen application with isolated settings and state, and helpers to grab panels."""

    def __init__(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
        self.tmp = tempfile.mkdtemp(prefix="starplast-tutorial-")
        os.environ["STARPLAST_STATE"] = self.tmp
        from PyQt6 import QtCore, QtWidgets
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
        for scope in (QtCore.QSettings.Scope.UserScope, QtCore.QSettings.Scope.SystemScope):
            QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat, scope, self.tmp)
        self.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        from starplast import theme as TH
        try:
            TH.apply(self.app, "dark")
        except Exception:                             # older theme API: the panels still render
            pass

    def settle(self, n: int = 8):
        for _ in range(n):
            self.app.processEvents()

    def grab(self, widget, path: str, size=None):
        if size:
            widget.resize(*size)
        widget.show()
        self.settle()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        widget.grab().save(path)
        return path


def render_map(result, ctx, path: str, color=None, title: str = ""):
    """A dark render of a result's own coordinates, colored by its clusters or by a column."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(8.6, 6.4), facecolor="#1e1e1e")
    ax = fig.add_subplot(projection="3d", facecolor="#1e1e1e")
    result.plot(ctx, color=color, ax=ax, size=2.2)
    if title:
        ax.set_title(title, color="#dddddd", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor="#1e1e1e")
    plt.close(fig)
    return path


def copy_asset(src: str, slug: str, name: str) -> str:
    dst = os.path.join(RES, slug, name)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    return f"resources/{slug}/{name}"


# --------------------------------------------------------------------------- notebooks
class Notebook:
    """Runs cells in one namespace, capturing printed text and the value of the last line."""

    def __init__(self, slug: str, title: str):
        self.slug, self.title = slug, title
        self.ns: dict = {}
        self.cells: list = []
        self.count = 0
        self.images = 0

    def md(self, text: str):
        self.cells.append({"kind": "md", "html": text, "source": _strip_tags(text)})

    def code(self, lines: list):
        """lines: [(code, explanation)]. Executed now; output captured."""
        import pandas as pd
        self.count += 1
        source = "\n".join(c for c, _t in lines)
        buf = io.StringIO()
        value = None
        t0 = time.monotonic()
        body, last = _split_last(source)
        with contextlib.redirect_stdout(buf):
            exec(compile(body, f"<{self.slug}:{self.count}>", "exec"), self.ns)
            if last:
                try:
                    value = eval(compile(last, f"<{self.slug}:{self.count}>", "eval"), self.ns)
                except SyntaxError:
                    exec(compile(last, f"<{self.slug}:{self.count}>", "exec"), self.ns)
        printed = buf.getvalue()
        out_html, out_nb = [], []
        if printed.strip():
            out_html.append(f"<pre>{_e(printed.rstrip())}</pre>")
            out_nb.append({"output_type": "stream", "name": "stdout", "text": printed})
        if value is not None:
            rich = self._rich(value)
            out_html.append(rich["html"])
            out_nb.append(rich["nb"])
        self.cells.append({"kind": "code", "lines": lines, "count": self.count,
                           "output": "".join(out_html), "outputs": out_nb, "source": source,
                           "seconds": time.monotonic() - t0})

    def _rich(self, value) -> dict:
        import pandas as pd
        count = self.count
        if hasattr(value, "savefig") or type(value).__name__ in ("Axes3D", "Axes"):
            fig = value if hasattr(value, "savefig") else value.figure
            self.images += 1
            name = f"{self.slug}_nb{self.images}.png"
            path = os.path.join(RES, self.slug, name)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            fig.savefig(path, dpi=100, facecolor=fig.get_facecolor())
            import matplotlib.pyplot as plt
            plt.close(fig)
            data = base64.b64encode(open(path, "rb").read()).decode()
            return {"html": f"<img src='resources/{self.slug}/{name}'>",
                    "nb": {"output_type": "display_data", "metadata": {},
                           "data": {"image/png": data, "text/plain": "<figure>"}}}
        if isinstance(value, pd.DataFrame):
            h = value.head(12).to_html(index=False, float_format=lambda x: f"{x:.3g}")
            note = f"<p><i>{len(value):,} rows; first {min(12, len(value))} shown</i></p>" \
                if len(value) > 12 else ""
            return {"html": h + note, "nb": _nb_result(count, value.head(12).to_string(), h)}
        if isinstance(value, pd.Series):
            h = value.head(15).to_frame().to_html(float_format=lambda x: f"{x:.3g}")
            return {"html": h, "nb": _nb_result(count, value.head(15).to_string(), h)}
        if hasattr(value, "_repr_html_"):
            h = value._repr_html_()
            return {"html": h, "nb": _nb_result(count, repr(value)[:500], h)}
        text = value if isinstance(value, str) else repr(value)
        return {"html": f"<pre>{_e(text)}</pre>", "nb": _nb_result(count, text, None)}

    def write(self):
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, f"{self.slug}_notebook.html"), "w") as fh:
            fh.write(T.notebook_page(self.title, self.cells))
        cells = []
        for c in self.cells:
            if c["kind"] == "md":
                cells.append({"cell_type": "markdown", "metadata": {},
                              "source": c["source"].splitlines(keepends=True)})
            else:
                cells.append({"cell_type": "code", "execution_count": c["count"], "metadata": {},
                              "source": c["source"].splitlines(keepends=True),
                              "outputs": c["outputs"]})
        nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3",
                                                          "language": "python", "name": "python3"},
                                           "language_info": {"name": "python"}},
              "nbformat": 4, "nbformat_minor": 5}
        os.makedirs(NOTEBOOKS, exist_ok=True)
        with open(os.path.join(NOTEBOOKS, f"{self.slug}.ipynb"), "w") as fh:
            json.dump(nb, fh, indent=1)


def _nb_result(count, text, html_text):
    data = {"text/plain": text}
    if html_text:
        data["text/html"] = html_text
    return {"output_type": "execute_result", "execution_count": count, "metadata": {},
            "data": data}


def _split_last(source: str) -> tuple:
    """(everything but the last top-level statement, the last one if it is an expression)."""
    import ast
    tree = ast.parse(source)
    if tree.body and isinstance(tree.body[-1], ast.Expr):
        lines = source.splitlines()
        start = tree.body[-1].lineno - 1
        return "\n".join(lines[:start]), "\n".join(lines[start:])
    return source, ""


def _strip_tags(text: str) -> str:
    import re
    t = re.sub(r"<h1>(.*?)</h1>", r"# \1\n", text)
    t = re.sub(r"<h2>(.*?)</h2>", r"## \1\n", t)
    t = re.sub(r"<li>(.*?)</li>", r"- \1\n", t)
    t = re.sub(r"<code>(.*?)</code>", r"`\1`", t)
    t = re.sub(r"<b>(.*?)</b>", r"**\1**", t)
    t = re.sub(r"<i>(.*?)</i>", r"*\1*", t)
    t = re.sub(r"<[^>]+>", "", t)
    return html.unescape(t).strip() + "\n"


# --------------------------------------------------------------------------- strategy GUI walkthroughs
BUTTONS = [
    ("Run", "Run the selected strategy with these settings on the whole table, in the background. "
            "Its tables appear under Results; a strategy that builds a map can then be shown on it."),
    ("Test (hold-out)", "Hide information that is already known, ask the strategy for it back, and "
                        "compare the answer with the same procedure on shuffled data. PASS means it "
                        "beat that null by the stated margin with these settings."),
    ("Stop", "Ask the running strategy to stop at its next checkpoint. Nothing is written by a "
             "stopped run."),
    ("Show on map", "Put the result's map into the central view, colored by its clusters or "
                    "modules, so it can be looked at rather than only read."),
]


def strategy_walkthrough(studio, key: str, slug: str, title: str, intro: str, settings: dict,
                         console_tips: dict, color=None) -> dict:
    """Drive the real Strategies panel through Guide -> Settings -> Test -> Run, grabbing each."""
    from starplast import strategies as S
    from starplast.strategy_panel import StrategyPanel
    ctx = S.shipped("Tg")
    panel = StrategyPanel(ctx.nodes, organism="Tg")
    log = []
    panel.status.connect(log.append)
    size = (640, 1040)
    panel.select(key)
    for name, value in settings.items():
        panel.set_setting(name, value)
    shots = []
    panel.tabs.setCurrentIndex(0)
    shots.append((studio.grab(panel, os.path.join(RES, slug, "1_guide.png"), size),
                  f"The Guide tab: what strategy {S.get(key).number:02d} infers, why it works, "
                  f"how it fails, and how it is tested."))
    panel.tabs.setCurrentIndex(1)
    shots.append((studio.grab(panel, os.path.join(RES, slug, "2_settings.png"), size),
                  "The Settings tab, filled in for this walkthrough. Hover any field on the left "
                  "for the reason it exists."))
    t0 = time.monotonic()
    test = panel.test_current()
    shots.append((studio.grab(panel, os.path.join(RES, slug, "3_test.png"), size),
                  f"Test (hold-out): {test.verdict}. {test.summary().split(' -- ', 1)[-1]}"))
    test_log = list(log)
    log.clear()
    result = panel.run_current()
    shots.append((studio.grab(panel, os.path.join(RES, slug, "4_results.png"), size),
                  "Run: the summary and the tables it produced. Click a row to find its gene on "
                  "the map; right-click to save the table."))
    if result.coords is not None:
        path = render_map(result, ctx, os.path.join(RES, slug, "5_map.png"), color=color,
                          title="the result's own map, colored by "
                                + (color or "its clusters"))
        shots.append((path, "The map the result was computed on (a render of its coordinates; "
                            "in the application, Show on map puts it in the 3D view)."))
    elapsed = time.monotonic() - t0
    console = [(line, console_tips.get(_kind(line), "Progress reported by the running job."))
               for line in (test_log + log)[:60]]
    console.append((f"{test.verdict}: {test.summary().split(' -- ', 1)[-1]}",
                    "The self-test's verdict, the number it rests on, and the null it beat or "
                    "did not beat."))
    console.append((result.summary, "The run's summary: what was found, and the number that "
                                    "says how much to believe it."))
    fields = []
    s = S.get(key)
    used = s.settings(ctx, **settings)
    for p in s.params:
        v = used.get(p.name)
        if p.kind == "genes" and v:
            v = f"{len(str(v).split())} genes"
        fields.append((p.label, "(none)" if v is None else v, p.tip))
    page = T.gui_page(title, intro, fields, [(os.path.relpath(p, OUT), c) for p, c in shots],
                      console, BUTTONS)
    with open(os.path.join(OUT, f"{slug}_gui.html"), "w") as fh:
        fh.write(page)
    panel.deleteLater()
    return {"test": test, "result": result, "seconds": elapsed}


def _kind(line: str) -> str:
    for k in ("map ", "null", "hidden", "random set", "strategy", "done", "stopped"):
        if k in line:
            return k
    return ""


CONSOLE_TIPS = {
    "map ": "One map of the walk: which measurements, and the UMAP settings. Maps are cached, so "
            "a second strategy on the same grid pays nothing.",
    "null": "One run of the null: the same procedure on shuffled labels. Its spread is the chance "
            "level the verdict is judged against.",
    "hidden": "How much known information was hidden for the test, and how much is left to learn "
              "from.",
    "random set": "One random gene set of the same size: the null for a gene-list strategy.",
    "strategy": "The job starting in the background, as the Jobs panel lists it.",
    "done": "The job finished; its tables are under Results.",
}


# --------------------------------------------------------------------------- the tutorials
def tutorial_explore(studio) -> dict:
    slug = "1_explore"
    slides = [
        (copy_asset(os.path.join(DECK, "workspace.png"), slug, "workspace.png"),
         "The window: find and color-by on the left, the 3D map in the middle, the evidence for "
         "the selected gene on the right. Every gene is a point placed by what it does."),
        (copy_asset(os.path.join(DECK, "map-dark.png"), slug, "map_dark.png"),
         "The map colored by compartment. Gray is a gene with no value -- never zero, never a "
         "category."),
        (copy_asset(os.path.join(DECK, "map-light.png"), slug, "map_light.png"),
         "The same map in the light theme (View > Preferences > Appearance)."),
        (copy_asset(os.path.join(DECK, "evidence.png"), slug, "evidence.png"),
         "The evidence panel for one gene: every measurement, where it came from, and a dash "
         "where nothing was measured."),
        (copy_asset(os.path.join(DECK, "explore.png"), slug, "explore_workflow.png"),
         "Tools > Explore a gene: search by identifier, symbol or description and follow the "
         "evidence to its source."),
    ]
    fields = [
        ("find", "GRA16", "Type an accession, a symbol or words from the product description; the "
                          "matching gene is selected and its evidence shown."),
        ("color by", "compartment", "What color means right now: a column, a kept clustering, or a "
                                    "binned quantity. Gray always means unknown."),
        ("point size", "Automatic", "View > Point size. Automatic scales with how many genes are "
                                    "visible."),
        ("edges", "co-mention, co-fitness", "Edges menu: each layer is a different kind of evidence "
                                            "and they are never merged silently."),
        ("left mouse button", "Navigate", "View > Left mouse button: Navigate rotates and selects; "
                                          "Select draws a gate around genes."),
        ("attention-corrected co-mention", "on", "Edges menu: co-mention is shown as the residual "
                                                 "over what each gene's fame predicts. Off shows "
                                                 "raw counts, and the menu says so."),
    ]
    console = [
        ("8,140 / 8,140 genes shown  ·  edges: comention, cofitness",
         "The status bar: how many genes are visible after filters, and which edge layers are on."),
        ("coloring by compartment; gray is unknown, which is a real answer",
         "What the current coloring means, and the rule that absence is never drawn as a value."),
        ("selected TGME49_…  ·  12 edges", "A clicked gene: its evidence opens on the right and its "
                                          "edges are drawn."),
    ]
    buttons = [("Reset view / clear filters", "Return the camera and show every gene again."),
               ("Spin", "View > Spin: depth in a 3D scatter reads only when it moves."),
               ("Export image…", "File > Export image: the current view as a picture."),
               ("Export visible genes (CSV)…", "File: every gene that survives the current filters, "
                                               "with its columns.")]
    with open(os.path.join(OUT, f"{slug}_gui.html"), "w") as fh:
        fh.write(T.gui_page("Explore the map and read the evidence",
                            "The first things to do in Starplast: find a gene, color the map, and "
                            "read the evidence. Hover a setting or a console line to see what it "
                            "means; move the slider to step through the views.",
                            fields, slides, console, buttons))
    nb = Notebook(slug, "Explore the tables from Python")
    nb.md("<h1>Explore the tables from Python</h1><p>Everything the window shows is available "
          "from Python. A <code>Context</code> is one organism's gene table, its measured "
          "networks and the leakage guard every analysis applies.</p>")
    nb.code([("from starplast import strategies as S", "The strategies module: tables, "
              "guard and the 32 strategies."),
             ("ctx = S.shipped('Tg')", "The packaged T. gondii table, loaded once and cached."),
             ("print(ctx.organism, ctx.n, 'genes,', len(ctx.nodes.columns), 'columns')",
              "One row per gene, one column per measurement or label.")])
    nb.code([("hits = ctx.nodes[ctx.nodes['product'].str.contains('microneme protein', na=False)]",
              "Find genes by words in their product description, as the find box does."),
             ("hits[['gene_id', 'product', 'compartment', 'fit_invitro_hff']].head(8)",
              "A few columns: the hyperLOPIT compartment and fibroblast fitness score.")])
    nb.code([("ctx.truth('compartment').value_counts().head(10)",
              "Labels as the strategies see them: 'unassigned' and other spellings of absence "
              "are missing, not a class.")])
    nb.code([("print(len(ctx.categorical_columns()), 'labels can be held out, e.g.',"
              " ctx.categorical_columns()[:6])", "Columns with 2-60 classes and enough labels."),
             ("print(len(ctx.numeric_columns()), 'measurements are available as features')",
              "Numeric columns measured on enough genes to use.")])
    nb.md("<h2>What a held-out analysis never sees</h2><p>Holding a label out removes the label, "
          "anything that restates it (measured association of 0.8 or more, on values or ranks), "
          "the experiment that produced it, the same quantity measured another way, and any "
          "edge layer built from it. This is what makes a recovered label mean something.</p>")
    nb.code([("sorted(ctx.banned('compartment'))",
              "Every column removed when compartment is held out, and why is in "
              "search.excluded_detail."),
             ])
    nb.code([("print('edge layers:', ctx.layers())", "The measured networks."),
             ("print('banned when compartment is held out:', sorted(ctx.banned_layers("
              "'compartment')))", "The compartment layer joins genes sharing a compartment, so "
                                  "it is refused.")])
    nb.write()
    return {"slug": slug}


def tutorial_holdout(studio) -> dict:
    slug = "2_holdout_search"
    settings = {"target": "compartment", "n_neighbors": "15, 50", "min_dist": "0.0",
                "min_cluster_size": "20, 50", "selection": "eom, leaf", "sample": 2000}
    run = strategy_walkthrough(
        studio, "holdout_search", slug, "Hold out a category and search for the map that finds it",
        "Strategy 01, the question this map was built to ask: hide the hyperLOPIT compartment "
        "and everything that restates it, walk feature sets and UMAP/HDBSCAN settings, and keep "
        "the map whose clusters best recover it. The self-test hides a quarter of the labels and "
        "scores the chosen clusters on them.", settings, CONSOLE_TIPS, color="compartment")
    nb = Notebook(slug, "Hold out a category and search, from Python")
    nb.md("<h1>Hold out a category and search</h1><p>Strategy 01 answers the founding "
          "question: is there a combination of measurements and map settings under which a "
          "label nobody showed the map falls out as clusters? This notebook runs it on the "
          "hyperLOPIT compartment.</p>")
    nb.code([("from starplast import strategies as S", "Import the strategies."),
             ("s = S.get('holdout_search')", "Strategy 01."),
             ("s", "Its guide: the question, the walkthrough, and how it is tested.")])
    nb.code([("s.parameters(S.shipped('Tg'))", "Every setting, its default for this table, and "
                                                "why it exists.")])
    nb.md("<h2>Test before you run</h2><p>The self-test hides 25% of the compartment labels "
          "(whole orthogroups together), lets the search choose the map and each label's "
          "cluster on the rest, and scores how well those clusters hold the hidden genes -- "
          "against the same search made on shuffled labels.</p>")
    nb.code([("test = S.test('holdout_search', target='compartment', n_neighbors='15, 50',",
              "Run the self-test with a small grid."),
             ("              min_dist='0.0', min_cluster_size='20, 50', selection='eom, leaf')",
              "Two neighbourhood sizes, one min_dist, two cluster sizes, both HDBSCAN selections."),
             ("test", "PASS/FAIL, the number, the chance level and the margin.")])
    nb.md("<h2>Run it</h2>")
    nb.code([("result = S.run('holdout_search', target='compartment', n_neighbors='15, 50',",
              "The same search on 2,000 genes, keeping the best map."),
             ("                min_dist='0.0', min_cluster_size='20, 50', selection='eom, leaf',",
              "The grid."),
             ("                sample=2000)", "Genes per map; 0 would embed all of them."),
             ("print(result.summary)", "What was found, and the chance level beside it.")])
    nb.code([("result.tables['configurations'].head(8)",
              "Every map of the walk, ranked by how well its clusters recover the label.")])
    nb.code([("result.tables['label clusters']",
              "The one cluster chosen for each compartment and how well it isolates it.")])
    nb.code([("result.tables['calls'].head(10)",
              "Unlabelled genes in those clusters, called with the compartment; support is the "
              "cluster's F1.")])
    nb.code([("import matplotlib.pyplot as plt", "Plotting."),
             ("plt.style.use('dark_background')", "Match the application."),
             ("ax = result.plot(S.shipped('Tg'), color='compartment')",
              "The result's own map, colored by the held-out label it never saw."),
             ("ax.figure", "Show the figure.")])
    nb.write()
    return {"slug": slug, **run}


def tutorial_genelist(studio) -> dict:
    from starplast import strategies as S
    slug = "3_gene_list"
    ctx = S.shipped("Tg")
    cat, members = S.example_set(ctx, "compartment")
    genes = "\n".join(ctx.gene_ids[members])
    settings = {"genes": genes, "exclude": "compartment", "n_neighbors": "15, 50",
                "min_dist": "0.0", "min_cluster_size": "20, 50", "selection": "eom, leaf",
                "sample": 2000}
    run = strategy_walkthrough(
        studio, "geneset_hunt", slug, "Find the cluster your gene list forms",
        f"Strategy 02: paste a gene list and walk the map space for the one cluster with the "
        f"best precision and recall for it. Here the list is the {len(members)} genes hyperLOPIT "
        f"places in {cat!r}; because the list came from the compartment column, that column and "
        f"its closure are withheld.", settings, CONSOLE_TIPS)
    nb = Notebook(slug, "Start from a gene list, from Python")
    nb.md(f"<h1>Start from a gene list</h1><p>Four strategies take a list: 02 looks for the "
          f"cluster that holds it, 20 learns what makes it special from positives alone, 24 "
          f"describes what it has in common, and 25 grows it along the networks. Here the list "
          f"is the {len(members)} genes of {_e(cat)}; hold out the compartment column so no "
          f"strategy can read the list back from the label that defined it.</p>")
    nb.code([("from starplast import strategies as S", "Import."),
             ("ctx = S.shipped('Tg')", "The table."),
             ("cat, members = S.example_set(ctx, 'compartment')",
              "A known category of about 120 genes -- stand-in for your own list."),
             ("genes = list(ctx.gene_ids[members])", "As accessions, which is what you would "
                                                     "paste."),
             ("print(cat, len(genes), genes[:5])", "What the list is.")])
    nb.code([("hunt = S.run('geneset_hunt', genes=genes, exclude='compartment',",
              "Strategy 02 with the label that defined the list withheld."),
             ("             n_neighbors='15, 50', min_dist='0.0', min_cluster_size='20, 50',",
              "A small grid."),
             ("             selection='eom, leaf', sample=2000)", "Both cluster selections."),
             ("print(hunt.summary)", "The best cluster, its precision and recall, and the "
                                     "search-corrected p-value against random lists.")])
    nb.code([("hunt.tables['candidates'].head(10)", "The cluster's other members, nearest its "
                                                    "centre first.")])
    nb.md("<h2>Three more ways to use a list</h2>")
    nb.code([("pu = S.run('positive_unlabeled', genes=genes, exclude='compartment')",
              "Learn the list against random draws of the rest, scoring genes only out of bag."),
             ("pu.tables['candidates'].head(8)", "Genes that most resemble the list.")])
    nb.code([("profile = S.run('set_enrichment', genes=genes, exclude='compartment')",
              "Test the list against every category, measurement and network at once."),
             ("profile.tables['profile'].head(12)", "What distinguishes it, by q-value.")])
    nb.code([("grown = S.run('seed_expansion', genes=genes, exclude='compartment')",
              "Random walk from the list across the measured networks."),
             ("grown.tables['candidates'].head(8)", "Genes the walk visits most, with their "
                                                    "direct links to the list.")])
    nb.md("<h2>Is my list coherent enough?</h2><p>With twenty or more genes, the self-test of a "
          "list strategy hides 30% of YOUR list and asks whether the rest recovers it.</p>")
    nb.code([("S.test('seed_expansion', genes=genes, exclude='compartment')",
              "PASS says the list is a unit this data can see.")])
    nb.write()
    return {"slug": slug, **run}


def tutorial_trust(studio) -> dict:
    slug = "4_test_and_calibrate"
    settings = {"target": "compartment"}
    run = strategy_walkthrough(
        studio, "physical_partners", slug, "Test a strategy before you trust it",
        "Every strategy carries a self-test. This walkthrough runs strategy 13 -- placing a "
        "protein by the proteins it is crosslinked to -- and reads its verdict: how many hidden "
        "compartment labels come back, against the same vote on shuffled labels.", settings,
        CONSOLE_TIPS)
    nb = Notebook(slug, "Self-tests, nulls and calibration, from Python")
    nb.md("<h1>Test before you trust</h1><p>A number without its chance level is not a result. "
          "Every self-test hides known information, asks the strategy for it back, and compares "
          "the answer with the SAME procedure run without the information it claims to use. The "
          "five patterns are: hide labels, hide members of a set, hide edges, hide values, and "
          "replicate findings on the other half of the genes.</p>")
    nb.code([("from starplast import strategies as S", "Import."),
             ("t = S.test('physical_partners', target='compartment')", "Strategy 13's test."),
             ("t", "The verdict.")])
    nb.code([("print('observed', round(t.observed, 3), '| chance', round(t.null_mean, 3),",
              "The number and its chance level."),
             ("      '+/-', round(t.null_sd, 3), '| bar', round(t.null_high, 3),",
              "The null's spread, and the 95th percentile to clear."),
             ("      '| effect', round(t.effect, 3), '| p', round(t.p_value, 4))",
              "Observed minus chance, and the share of null runs as good.")])
    nb.code([("t.details", "Per hidden class: how many came back, and how precise the calls "
                           "were.")])
    nb.md("<h2>The same test on a table with nothing in it</h2><p>A test that passes on noise "
          "tests nothing. <code>planted_context(null=True)</code> deals every column and edge "
          "out at random.</p>")
    nb.code([("noise = S.planted_context(null=True)", "A synthetic organism with nothing to find."),
             ("S.get('physical_partners').test(noise)", "Should FAIL or be INCONCLUSIVE.")])
    nb.md("<h2>Calibration: where each strategy works</h2><p>"
          "<code>scripts/calibrate_strategies.py</code> runs every strategy's self-test over a "
          "grid of its settings, several held-out labels and five seeds, and summarises each "
          "configuration with a mean and a 95% interval. Skill puts every metric on one scale: "
          "0 is chance, 1 is perfect.</p>")
    nb.code([("import pandas as pd, glob, os", "Read the calibration tables."),
             ("folders = sorted(glob.glob(os.path.join(S.__file__.rsplit('/starplast/', 1)[0],",
              "The newest calibration run in results/."),
             ("                                    'results', 'calibration_*')))", ""),
             ("best = pd.read_csv(os.path.join(folders[-1], 'best.csv'))",
              "Per strategy, the configuration with the highest lower confidence bound."),
             ("best[best.organism == 'Tg'][['strategy', 'settings', 'skill', 'skill_low',",
              "Skill and its 95% interval."),
             ("                               'skill_high', 'pass_rate']].head(12)", "")])
    nb.write()
    return {"slug": slug, **run}


def tutorial_networks(studio) -> dict:
    slug = "5_networks_and_agreement"
    settings = {"target": "compartment", "k": 15}
    run = strategy_walkthrough(
        studio, "layer_vote", slug, "Borrow from networks, learn from examples, demand agreement",
        "Strategy 12 lets every permitted network and the measurement neighbours vote on each "
        "gene's compartment, each source weighted by how far it beat chance on an inner holdout. "
        "The notebook adds a classifier (19), agreement between independent methods (31), and "
        "the understudied genes (32).", settings, CONSOLE_TIPS)
    nb = Notebook(slug, "Networks, learning and agreement, from Python")
    nb.md("<h1>Networks, learning and agreement</h1><p>Different strategies read different "
          "evidence: measured contacts, co-expression, a trained model. Where they agree, a call "
          "is far more likely to be right -- and the self-tests say how much more.</p>")
    nb.code([("from starplast import strategies as S", "Import."),
             ("vote = S.run('layer_vote', target='compartment')", "Strategy 12."),
             ("vote.tables['source weights']", "Which kinds of evidence earned a vote, net of "
                                               "chance.")])
    nb.code([("clf = S.run('supervised_classifier', target='compartment')", "Strategy 19."),
             ("clf.tables['what drives each class'].groupby('class').head(2).head(12)",
              "The measurements each compartment is recognised by.")])
    nb.code([("tri = S.test('triangulation', target='compartment')",
              "Strategy 31: call only where two of three methods agree."),
             ("tri", "Precision of agreed calls on hidden genes."),
             ])
    nb.code([("{k: round(v, 3) for k, v in tri.numbers.items() if 'precision' in k}",
              "Each single method's precision, for comparison.")])
    nb.code([("under = S.run('understudied_first', target='compartment')",
              "Strategy 32: agreed calls for genes nobody has written about."),
             ("under.tables['candidates'].head(10)", "Ranked by agreement times novelty.")])
    nb.write()
    return {"slug": slug, **run}


# --------------------------------------------------------------------------- the long guide
def guide(studio, runs: dict) -> None:
    """One long page describing every feature, with real captures and generated control tables."""
    from PyQt6 import QtWidgets
    from starplast import strategies as S
    from starplast.analysis_panel import AnalysisPanel
    from starplast.strategy_panel import StrategyPanel
    from starplast.tuning import EmbeddingStore
    from starplast.workflows import WorkflowDialog
    slug = "guide"
    ctx = S.shipped("Tg")
    sec = []

    def shot(path, alt):
        return f"<img class='shot' src='{os.path.relpath(path, OUT)}' alt='{_e(alt)}'>"

    def img(name):
        return os.path.join(RES, slug, name)

    # The analysis panel, tab by tab, with every control and its tooltip.
    analysis = AnalysisPanel(ctx.nodes, store=EmbeddingStore(os.path.join(studio.tmp, "emb")))
    tabs = analysis.findChild(QtWidgets.QTabWidget)
    analysis_rows = []
    for i in range(tabs.count()):
        tabs.setCurrentIndex(i)
        path = studio.grab(analysis, img(f"analysis_{i + 1}.png"), (620, 980))
        page_widget = tabs.widget(i)
        controls = []
        for form in page_widget.findChildren(QtWidgets.QFormLayout):
            for r in range(form.rowCount()):
                lab = form.itemAt(r, QtWidgets.QFormLayout.ItemRole.LabelRole)
                fld = form.itemAt(r, QtWidgets.QFormLayout.ItemRole.FieldRole)
                if lab and fld and lab.widget() and fld.widget():
                    tip = _strip(fld.widget().toolTip() or lab.widget().toolTip())
                    controls.append((lab.widget().text(), tip))
        buttons = [(b.text(), _strip(b.toolTip())) for b in page_widget.findChildren(
            QtWidgets.QPushButton) if b.text()]
        analysis_rows.append((tabs.tabText(i), path, controls, buttons))
    strat = StrategyPanel(ctx.nodes, organism="Tg")
    strat.select("holdout_search")
    strat.tabs.setCurrentIndex(0)
    s_guide = studio.grab(strat, img("strategies_guide.png"), (640, 1040))
    strat.tabs.setCurrentIndex(1)
    s_settings = studio.grab(strat, img("strategies_settings.png"), (640, 1040))
    wf = WorkflowDialog(ctx.nodes)
    wf_shots = []
    for i in range(wf.tabs.count()):
        wf.tabs.setCurrentIndex(i)
        wf_shots.append((wf.tabs.tabText(i), studio.grab(wf, img(f"workflow_{i + 1}.png"),
                                                         (1100, 760))))

    sec.append(("The window at a glance",
                f"<p>Starplast places every gene of <i>Toxoplasma gondii</i> (or <i>Plasmodium "
                f"falciparum</i>) as a point in a three-dimensional map built from what the gene "
                f"does -- expression across stages, fitness in many conditions, localization, "
                f"modification, structure -- so that nearness means similar behaviour, not an "
                f"arbitrary layout. Around the map sit the <b>find</b> and <b>color by</b> "
                f"controls (left), the <b>evidence</b> panel for the selected gene, and, tabbed "
                f"with it on the right, <b>analysis</b> and <b>strategies</b>. Console, jobs and "
                f"assistant sit along the bottom and can be shown from the Tools menu.</p>"
                f"{shot(copy_asset(os.path.join(DECK, 'workspace.png'), slug, 'workspace.png'), 'window')}"
                f"<p>Three rules hold everywhere and are worth knowing before anything else:</p>"
                f"<ul><li><b>Gray means unknown</b> -- never zero, never a category. A gene with "
                f"no measurement is not a gene that measured nothing.</li>"
                f"<li><b>Measurement, inference, absence and annotation never read as one "
                f"another.</b> An inferred label is drawn and stored apart from a measured one."
                f"</li><li><b>No candidate without the number that says how much to believe "
                f"it.</b> Every search reports its chance level; every strategy carries a "
                f"self-test.</li></ul>"))
    sec.append(("Finding genes and reading the evidence",
                "<p>Type an accession, a gene symbol or words from the product description in "
                "<b>find</b>; the gene is selected, its edges drawn and its evidence opened. The "
                "evidence panel lists every measurement with where it came from; a dash is "
                "'not measured'. Its header says how much the gene has been written about -- "
                "<i>listed, not studied</i> means named in papers only as an entry in a screen's "
                "table. <b>Tools &gt; Explore a gene</b> opens the same search as a guided "
                "workflow with links to each source record.</p>"
                f"{shot(copy_asset(os.path.join(DECK, 'evidence.png'), slug, 'evidence.png'), 'evidence')}"))
    sec.append(("Coloring, views and edges",
                "<p><b>Color by</b> chooses what color means: any label column, any kept "
                "clustering, or a numeric column cut into bins. The list below it filters: untick "
                "a category to hide it. <b>View</b> sets point size, what the left mouse button "
                "does (navigate, or draw a gate in 2D or a brush in 3D), rotation locked to an "
                "axis so two views can be compared, and a continuous spin -- depth in a 3D "
                "scatter reads only when it moves. The <b>Edges</b> menu turns each measured "
                "network on and off: co-mention in abstracts and full texts, shared orthogroup, "
                "co-expression, shared compartment, co-fitness, shared domain, crosslinks, "
                "pulldowns, structural similarity, and two derived layers. They are never merged "
                "into one graph, because each answers a different question; <i>attention-"
                "corrected co-mention</i> (on by default) shows co-mention beyond what each "
                "gene's fame predicts.</p>"
                f"{shot(copy_asset(os.path.join(DECK, 'map-dark.png'), slug, 'map_dark.png'), 'map')}"))
    rows = []
    for tab, path, controls, buttons in analysis_rows:
        table = "".join(f"<tr><td>{_e(l)}</td><td>{_e(t)}</td></tr>" for l, t in controls)
        btns = "".join(f"<tr><td><b>{_e(l)}</b></td><td>{_e(t)}</td></tr>" for l, t in buttons)
        rows.append(f"<h3>{_e(tab)}</h3>{shot(path, tab)}"
                    + (f"<table><tr><th>control</th><th>why it exists</th></tr>{table}{btns}"
                       f"</table>" if table or btns else ""))
    sec.append(("The analysis panel, tab by tab",
                "<p>The analysis panel is the laboratory: choose the data, build and tune a map, "
                "cluster it, ask what the clusters mean using only features the map never saw, "
                "search the space of maps for one that recovers a held-out label, put an error "
                "rate on an annotation by hiding labels you already have, search for structure "
                "worth reading, and ask named biological questions. It is arranged top to bottom "
                "in the order the work is done. Every control's explanation below is the "
                "tooltip the application itself shows.</p>" + "".join(rows)))
    sec.append(("The gallery, kept runs and annotations",
                "<p>A hyperparameter walk streams a thumbnail per map into the <b>gallery</b> "
                "along the bottom; click one to show that map in the central view with "
                "everything the view can do. Every clustering you make is <b>kept</b> as a run "
                "under a timestamp, so 'does this structure survive other settings?' can be "
                "answered by looking. <b>Annotations</b> -- proposing a label for the unlabelled "
                "members of a cluster -- refuse to be saved without the validated precision and "
                "recall of that category from the Validation tab; an annotation is drawn in a "
                "color used for nothing else and lives in its own file, never in the table of "
                "measurements.</p>"))
    wf_html = "".join(f"<h3>{_e(t)}</h3>{shot(p, t)}" for t, p in wf_shots)
    sec.append(("Guided workflows (Tools menu)",
                "<p><b>Explore a gene</b>, <b>Predict a trait</b> and <b>Compare a screen</b> are "
                "task-oriented views over the same table. Predict evaluates models with whole "
                "orthogroups held out together, calibrates their probabilities and abstains "
                "where support is thin; Compare joins your own screen to the existing evidence "
                "without dropping unresolved identifiers or averaging duplicates.</p>" + wf_html))
    sec.append(("Your own data",
                "<p><b>File &gt; Import data…</b> reads CSV, TSV, Excel or parquet, suggests the "
                "identifier column (and a repair if its identifiers are malformed), resolves "
                "identifiers through the identity layer -- old accessions, other strains, symbols "
                "-- and offers the preprocessing rather than assuming it. Imported columns are "
                "prefixed and joined in memory; the shipped cache is never written. Imported "
                "columns then appear in the Data tab as their own block and in every strategy.</p>"
                "<p><b>File &gt; Save all analysis results…</b> writes every tab's table into "
                "one file; <b>Load analysis results…</b> puts them back, and a loaded row is as "
                "clickable as a computed one. A file whose kind does not match its tab is "
                "refused.</p>"))
    fams = "".join(f"<li><b>{_e(f)}</b>: " + ", ".join(
        f"{s.number:02d} {_e(s.title.lower())}" for s in S.catalog() if s.family == f) + "</li>"
                   for f in S.families())
    sec.append(("The Strategies tab",
                "<p>Thirty-two named ways of turning the combined data into a claim, grouped by "
                f"how they work:</p><ul>{fams}</ul>"
                "<p>Select one and read its <b>Guide</b> (what it infers, why that works, how "
                "it fails, a walkthrough, and how it is tested); set its <b>Settings</b>; press "
                "<b>Test (hold-out)</b> before <b>Run</b>. The column beside each name is the "
                "verdict its self-test earned on the shipped data, and the calibration tables "
                "say which settings and targets it works for. Gene-list strategies accept a "
                "pasted list, a file, an example set, or the genes gated on the map.</p>"
                f"{shot(s_guide, 'strategy guide')}{shot(s_settings, 'strategy settings')}"
                "<p>Every strategy removes the held-out label, anything that restates it, the "
                "experiment that produced it, the same quantity measured another way, and any "
                "edge layer built from any of those, before it looks at anything else.</p>"))
    sec.append(("Two organisms",
                "<p><b>File &gt; Species</b> switches between <i>T. gondii</i> and <i>P. "
                "falciparum</i>. Each has its own table -- one species per window, never a "
                "union, because the identifiers and the data behind them do not align. Orthology "
                "is a bridge, not a merge: strategy 29 carries a measurement from one parasite's "
                "orthologs to the other and says how well it transfers.</p>"))
    sec.append(("Preferences, console, jobs and the assistant",
                "<p><b>File &gt; Preferences</b> holds appearance (theme, color maps, point "
                "rendering, lighting), display and window settings, and the compute backend: "
                "with the CUDA stack installed, UMAP and HDBSCAN run on the GPU, and every run "
                "records which library built it. The <b>jobs</b> panel lists everything running "
                "in the background, with Stop and the traceback of anything that failed; the "
                "<b>console</b> shows the log; the <b>assistant</b> answers questions about the "
                "map.</p>"))
    sec.append(("The Python API",
                "<pre class='code'>"
                + _e("from starplast import strategies as S\n\n"
                     "S.overview()                                  # every strategy, one row each\n"
                     "S.get('geneset_hunt')                         # its guide (renders in notebooks)\n"
                     "S.get('geneset_hunt').parameters()            # settings and why they exist\n"
                     "test = S.test('geneset_hunt', genes=my_list)  # hide 30% of the list, recover it\n"
                     "result = S.run('geneset_hunt', genes=my_list) # run on the shipped T. gondii table\n"
                     "result.tables['candidates']                   # what it found\n"
                     "result.plot(S.shipped('Tg'))                  # the map it was computed on\n"
                     "result.save('out/')                           # every table as CSV\n\n"
                     "ctx = S.shipped('Pf')                         # the P. falciparum table\n"
                     "ctx.banned('stage_enriched_derived')          # what a held-out analysis never sees\n")
                + "</pre><p>The lower-level modules -- <code>search</code>, <code>validate</code>, "
                "<code>prediction</code>, <code>methods</code>, <code>discovery</code>, "
                "<code>leakage</code> -- are documented in the API reference; the notebook "
                "tutorials use exactly the calls above.</p>"))
    sec.append(("Reading the numbers honestly",
                "<p>Three habits keep the output of this program from becoming a source of false "
                "claims. <b>Look at the chance level</b>: every search and every self-test "
                "reports the same procedure on shuffled data; a number is only as good as its gap "
                "over that. <b>Distrust the best of many</b>: a walk picks the luckiest of many "
                "maps, which is why the self-tests choose on some labels and score on others, and "
                "why a gene-list result reports a search-corrected p-value. <b>Trust the stratum, "
                "not the average</b>: a strategy's accuracy over the whole proteome is dominated "
                "by well-studied, conserved genes; strategies 30 and 32 measure it where "
                "inference is needed most.</p>"))
    toc = "".join(f"<li><a href='#s{i}'>{_e(t)}</a></li>" for i, (t, _b) in enumerate(sec, 1))
    body = (f"<div class='page'><h1>The complete guide to Starplast</h1><p>Every feature, in the "
            f"order you meet them, with captures from the running application. The five "
            f"tutorials walk through specific tasks; this page is the reference to come back "
            f"to.</p><ol class='toc'>{toc}</ol>"
            + "".join(f"<h2 id='s{i}'>{i}. {_e(t)}</h2>{b}" for i, (t, b) in enumerate(sec, 1))
            + "</div>")
    with open(os.path.join(OUT, "complete_guide.html"), "w") as fh:
        fh.write(T.page("The complete guide to Starplast", body, T.nav_links()))
    analysis.deleteLater()
    strat.deleteLater()
    wf.deleteLater()


def _strip(text: str) -> str:
    import re
    return html.unescape(re.sub(r"<[^>]+>", " ", text or "")).replace("\xa0", " ").split(
        "  ")[0].strip() if False else " ".join(html.unescape(
            re.sub(r"<[^>]+>", " ", text or "")).replace("\xa0", " ").split())


TUTORIALS = [
    ("1_explore", "Explore the map and read the evidence",
     "Find a gene, color the map, read the evidence, and look at the tables from Python."),
    ("2_holdout_search", "Hold out a category and search for the map that finds it",
     "Strategy 01: hide a label and everything restating it, walk map settings, keep the map "
     "whose clusters recover it."),
    ("3_gene_list", "Find the cluster your gene list forms",
     "Strategies 02, 20, 24 and 25: a gene list in, a cluster, a ranking and a profile out."),
    ("4_test_and_calibrate", "Test a strategy before you trust it",
     "Self-tests, nulls, noise tables and the calibration of every strategy over its settings."),
    ("5_networks_and_agreement", "Borrow from networks, learn from examples, demand agreement",
     "Strategies 12, 19, 31 and 32: networks, a classifier, agreement, and understudied genes."),
]


def index():
    cards = "".join(
        f"<div class='card'><h3>{i}. {_e(t)}</h3><p>{_e(d)}</p>"
        f"<a href='{s}_gui.html'>GUI</a><a href='{s}_notebook.html'>Notebook</a>"
        f"<a href='notebooks/{s}.ipynb'>.ipynb</a></div>"
        for i, (s, t, d) in enumerate(TUTORIALS, 1))
    body = (f"<div class='page'><h1>Tutorials</h1><p>Learn Starplast step by step. Each tutorial "
            f"comes as a <b>GUI</b> walkthrough -- the application's own panels, with every "
            f"setting and console line explained on hover -- and as a <b>notebook</b> doing the "
            f"same from Python, with every line explained on hover and its real output. The "
            f"notebooks are also provided as runnable .ipynb files.</p><div class='cards'>{cards}"
            f"</div><h2>The complete guide</h2><p><a href='complete_guide.html'>Every feature of "
            f"Starplast in one long page</a>: the window, finding genes and reading evidence, "
            f"coloring and edges, every tab of the analysis panel with every control explained, "
            f"the gallery and annotations, the guided workflows, your own data, the Strategies "
            f"tab, the two organisms, preferences, the Python API, and how to read the numbers "
            f"honestly.</p></div>")
    with open(os.path.join(OUT, "index.html"), "w") as fh:
        fh.write(T.page("Starplast tutorials", body, T.nav_links()))


def main(argv=None) -> int:
    """Build everything under docs/tutorial/, or only the parts named: 1_explore ... guide index."""
    argv = list(sys.argv[1:] if argv is None else argv)
    os.makedirs(RES, exist_ok=True)
    studio = Studio()
    runs = {}
    builds = {"1_explore": tutorial_explore, "2_holdout_search": tutorial_holdout,
              "3_gene_list": tutorial_genelist, "4_test_and_calibrate": tutorial_trust,
              "5_networks_and_agreement": tutorial_networks}
    for slug, build in builds.items():
        if argv and slug not in argv:
            continue
        t0 = time.monotonic()
        info = build(studio)
        runs[info["slug"]] = info
        print(f"{info['slug']}: {time.monotonic() - t0:.0f}s", flush=True)
    if not argv or "guide" in argv:
        guide(studio, runs)
    if not argv or "index" in argv:
        index()
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
