"""The Start-here tab: the question tree, the recommendation rules, and the panel that draws them.

The tree is walked exhaustively -- every first answer, every option of every choice it leads to --
because the promise the tab makes is that no answer is a dead end. The recommendations are checked
against the three things they are supposed to respect (what the goal is for, what task the strategy
does, and what the calibration sweep measured for it on THAT organism), and every setting they
prefill is handed to `Strategy.settings`, which is the only judge of whether a setting is real.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6 import QtWidgets  # noqa: E402

from starplast import guided as G  # noqa: E402
from starplast import organisms as ORG  # noqa: E402
from starplast import scorecard as SC  # noqa: E402
from starplast import strategies as S  # noqa: E402

ORGANISM = ORG.codes(available=True)[0]
#: How many options of a long column list a tree walk tries. Three hundred and fifty measurements
#: all lead to the same next question, and walking them all would say nothing more than three do.
COLUMN_SAMPLE = 3


@pytest.fixture(scope="module")
def ctx():
    return S.shipped(ORGANISM)


@pytest.fixture(scope="module")
def label(ctx):
    """The label the Strategies tab itself starts on for this organism."""
    return S.default_category(ctx)


@pytest.fixture(scope="module")
def number(ctx):
    """The numeric screen the Strategies tab itself starts on for this organism."""
    return S.default_numeric(ctx)


@pytest.fixture(scope="module")
def small():
    """A table with no calibration behind it: every grade is unknown, as a user's own table is."""
    rng = np.random.default_rng(0)
    n = 90
    nodes = pd.DataFrame({"gene_id": [f"SYN_{i:05d}" for i in range(n)],
                          "product": [f"protein {i}" for i in range(n)],
                          "compartment": ["A"] * 30 + ["B"] * 30 + [None] * 30})
    for j in range(5):
        nodes[f"m{j}"] = rng.normal(size=n) + (np.arange(n) < 45) * 1.5
    graph = {"xlms__a": np.arange(0, 40), "xlms__b": np.arange(1, 41),
             "xlms__w": np.linspace(0.1, 1.0, 40)}
    return S.Context(nodes, graph=graph, organism=ORGANISM)


def _sample(step: str, answers: dict, ctx) -> list:
    """The answers a walk tries for one step: every choice, a few columns, or one made-up subject."""
    if step == "gene":
        return [str(ctx.gene_ids[0])]
    if step == "genes":
        return [[str(g) for g in ctx.gene_ids[:20]]]
    options = G.options(step, answers, ctx)
    values = [o.value for o in options]
    return values[:COLUMN_SAMPLE] if G.STEPS[step].kind == "column" else values


def _walk(answers: dict, ctx, seen: list, depth: int = 0):
    """Answer every way this branch can be answered, recording each end state reached."""
    assert depth < 10, f"the tree does not end: {answers}"
    step = G.next_step(answers)
    if not step:
        seen.append(dict(answers))
        return
    values = _sample(step, answers, ctx)
    assert values, f"{step} offers nothing to answer it with"
    for value in values:
        nxt = dict(answers)
        nxt[step] = value
        assert G.next_step(nxt) != step, f"{step}={value!r} leads back to itself"
        _walk(nxt, ctx, seen, depth + 1)


# --------------------------------------------------------------------------- the tree
def test_every_first_answer_reaches_an_end_state(ctx):
    for option in G.options("have", {}, ctx):
        seen = []
        _walk({"have": option.value}, ctx, seen)
        assert seen, f"{option.value!r} reaches no end state"
        assert all(G.done(a) for a in seen)


def test_every_option_leads_somewhere_and_nothing_is_a_dead_end(ctx):
    seen = []
    _walk({}, ctx, seen)
    assert len(seen) > 20
    for answers in seen:
        assert G.recommend(answers, ctx), f"no recommendation for {answers}"


def test_the_first_question_is_asked_first_and_the_trail_is_empty():
    assert G.next_step({}) == "have"
    assert G.trail({}) == []
    assert not G.done({})


def test_nothing_specific_skips_the_goal_question(ctx):
    answers = {"have": G.HAVE_NOTHING, "space": ORGANISM}
    assert G.done(answers) and G.goal_of(answers) == G.TOUR
    assert "goal" not in G.route(answers)


def test_comparing_two_conditions_asks_for_the_baseline(ctx):
    answers = {"have": G.HAVE_NUMBER, "space": ORGANISM,
               "measure": G.options("measure", {}, ctx)[0].value, "goal": G.COMPARE}
    assert G.next_step(answers) == "baseline"
    assert answers["measure"] not in [o.value for o in G.options("baseline", answers, ctx)]


def test_one_gene_given_a_label_is_asked_which_label(ctx):
    answers = {"have": G.HAVE_GENE, "space": ORGANISM, "gene": str(ctx.gene_ids[0]),
               "goal": G.PREDICT}
    assert G.next_step(answers) == "label"


def test_changing_an_earlier_answer_drops_what_followed_it(ctx, label):
    answers = {"have": G.HAVE_LABEL, "space": ORGANISM, "label": label, "goal": G.PREDICT}
    assert G.done(answers)
    steps = [step for step, _value, _text in G.trail(answers)]
    assert steps == ["have", "space", "label", "goal"]


def test_every_step_has_a_question_and_a_tooltip():
    for step in G.STEPS.values():
        assert step.question.endswith("?") and len(step.question.split()) <= 6
        assert len(step.tip) > 40


def test_column_options_carry_their_coverage(ctx):
    for option in G.categorical_options(ctx) + G.numeric_options(ctx):
        assert "genes" in option.detail and option.hint
    counts = [int(o.detail.split(" ")[0].replace(",", "")) for o in G.categorical_options(ctx)]
    assert min(counts) > 0


def test_the_space_question_offers_the_installed_spaces():
    assert [o.value for o in G.options("space", {"have": G.HAVE_SET}, None)] == \
        ORG.codes(available=True)


# --------------------------------------------------------------------------- recommendations
def _answers(ctx, **kw):
    out = {"space": ctx.organism}
    out.update(kw)
    return out


def test_a_label_and_a_goal_give_a_short_ranked_list_with_reasons(ctx, label):
    answers = _answers(ctx, have=G.HAVE_LABEL, label=label, goal=G.PREDICT)
    recs = G.recommend(answers, ctx)
    assert G.MIN_RECOMMENDED <= len(recs) <= G.MAX_RECOMMENDED
    assert [r.score for r in recs] == sorted((r.score for r in recs), reverse=True)
    assert all(len(r.why) > 60 and r.why.endswith(".") for r in recs)
    assert all(label.replace("_", " ") in r.why for r in recs)


def test_recommendations_do_the_task_the_goal_asks_for_or_are_what_it_is_for(ctx, label):
    for goal in (G.PREDICT, G.EXPLAIN, G.LEARNABLE):
        answers = _answers(ctx, have=G.HAVE_LABEL, label=label, goal=goal)
        for rec in G.recommend(answers, ctx):
            fits_task = rec.task in G.GOALS[goal].tasks
            is_for_it = rec.key in G.KEYS[(goal, "label")]
            assert fits_task or is_for_it, f"{rec.key} is neither for {goal} nor doing its task"


def test_a_reliable_strategy_outranks_a_weak_one_it_is_otherwise_tied_with(ctx):
    goal = G.GOALS[G.PREDICT]
    strategy = S.get("supervised_classifier")
    assert G._score(strategy, goal, 0, "reliable") > G._score(strategy, goal, 0, "weak")
    assert G._score(strategy, goal, 0, "weak") > G._score(strategy, goal, 0, "no skill")
    assert G._score(strategy, goal, 0, "reliable") > G._score(strategy, goal, 3, "reliable")


def test_the_grade_read_is_the_one_measured_for_that_organism():
    codes = ORG.codes(available=True)
    graded = {code: {s.key: G.grade_of(s.key, code) for s in S.catalog()} for code in codes}
    assert all(any(v for v in per.values()) for per in graded.values())
    if len(codes) > 1:
        assert graded[codes[0]] != graded[codes[1]], "the grades must be read per organism"


def test_a_numeric_screen_is_never_answered_with_a_label_calling_strategy(ctx, number):
    answers = _answers(ctx, have=G.HAVE_NUMBER, measure=number, goal=G.PREDICT)
    recs = G.recommend(answers, ctx)
    assert recs and all(r.task != SC.T_LABEL for r in recs)
    assert all(r.settings for r in recs)


def test_a_gene_list_reaches_the_strategies_that_take_one(ctx):
    genes = [str(g) for g in ctx.gene_ids[:25]]
    answers = _answers(ctx, have=G.HAVE_SET, genes=genes, goal=G.MORE_LIKE)
    recs = G.recommend(answers, ctx)
    assert any(r.settings.get("genes") == genes for r in recs)
    assert all("25 genes" in r.why for r in recs if r.settings.get("genes"))
    # A strategy that takes the very genes chosen leads one that ranks the whole table.
    assert recs[0].settings.get("genes") == genes


def test_every_prefilled_setting_is_accepted_by_the_strategy(ctx):
    seen = []
    _walk({}, ctx, seen)
    for answers in seen:
        for rec in G.recommend(answers, ctx):
            settings = S.get(rec.key).settings(ctx, **rec.settings)
            for name, value in rec.settings.items():
                assert settings[name] is value or settings[name] == value


def test_a_weak_only_list_says_so_and_a_reliable_one_does_not(ctx, label, monkeypatch):
    answers = _answers(ctx, have=G.HAVE_LABEL, label=label, goal=G.PREDICT)
    assert not G.caveat(G.recommend(answers, ctx))              # reliable ones exist here
    monkeypatch.setattr(G, "grade_of", lambda key, organism: "weak")
    weak = G.recommend(answers, ctx)
    assert weak and all(r.weak and r.grade == "weak" for r in weak)
    said = G.caveat(weak)
    assert "reliable" in said and weak[0].title in said and "weak" in said
    assert all("only weak" in r.why for r in weak)
    assert G.caveat([])                                         # and nothing at all is said so too


def test_a_strategy_this_table_cannot_fill_is_not_recommended(small):
    unusable = {s.key for s in S.catalog() if not G.usable(s, small)}
    assert unusable, "this table was supposed to be too thin for some strategies"
    seen = []
    _walk({}, small, seen)
    for answers in seen:
        recommended = {r.key for r in G.recommend(answers, small)}
        assert recommended and not (recommended & unusable)


def test_the_other_panels_are_offered_where_they_answer_it_better(ctx, label):
    maps = G.views(_answers(ctx, have=G.HAVE_LABEL, label=label, goal=G.LEARNABLE), ctx)
    assert [v["view"] for v in maps] == [G.VIEW_MAPS]
    star = G.views(_answers(ctx, have=G.HAVE_GENE, gene=str(ctx.gene_ids[0]), goal=G.PARTNERS), ctx)
    assert G.VIEW_STAR in [v["view"] for v in star]
    assert all(len(v["why"]) > 40 for v in maps + star)


def test_a_gene_file_is_read_the_way_a_pasted_list_is(tmp_path, ctx):
    path = tmp_path / "hits.csv"
    ids = [str(g) for g in ctx.gene_ids[:5]]
    path.write_text("\n".join(f"{g},0.5,notes" for g in ids))
    assert G.read_gene_file(str(path)) == ids
    assert set(G.gene_ids({"genes": ids}, ctx)) == set(ids)


# --------------------------------------------------------------------------- the panel
@pytest.fixture(scope="module")
def panel(ctx):
    from starplast.guided_panel import GuidedPanel
    QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    widget = GuidedPanel(ctx, gated=lambda: [str(g) for g in ctx.gene_ids[:12]])
    widget.resize(620, 820)
    yield widget
    widget.close()


CONTROLS = (QtWidgets.QAbstractButton, QtWidgets.QLineEdit, QtWidgets.QPlainTextEdit,
            QtWidgets.QListWidget)


def _controls(widget) -> list:
    return [w for w in widget.findChildren(QtWidgets.QWidget) if isinstance(w, CONTROLS)]


def _tooltips(widget):
    missing = [w.__class__.__name__ + ":" + (getattr(w, "text", lambda: "")() or w.objectName())
               for w in _controls(widget) if not w.toolTip()]
    assert not missing, f"controls without a tooltip: {missing}"


def test_the_panel_builds_and_every_control_explains_itself(panel):
    panel.restart()
    assert panel.question.text() == G.STEPS["have"].question
    _tooltips(panel)
    for value, step in ((G.HAVE_SET, "genes"), (G.HAVE_LABEL, "label"),
                        (G.HAVE_NUMBER, "measure"), (G.HAVE_GENE, "gene")):
        panel.restart()
        panel.answer("have", value)
        panel.answer("space", ORGANISM)
        assert panel.question.text() == G.STEPS[step].question
        _tooltips(panel)


def test_walking_the_panel_to_the_end_shows_cards_with_their_bars_and_buttons(panel, label):
    panel.restart()
    panel.answer("have", G.HAVE_LABEL)
    panel.answer("space", ORGANISM)
    panel.answer("label", label)
    panel.answer("goal", G.PREDICT)
    assert panel.cards and len(panel.cards) == len(panel.recommendations())
    card = panel.cards[0]
    from starplast import strategy_card as CARD
    bars = card.findChild(CARD.Scorecard)
    assert bars is not None and len(bars.bars) == 4
    assert card.run_btn.isEnabled() and card.open_btn.isEnabled()
    _tooltips(panel)
    sent = []
    card.run.connect(lambda key, settings: sent.append((key, settings)))
    card.run_btn.click()
    assert sent and sent[0][0] == panel.recommendations()[0].key
    assert sent[0][1]["target"] == label


def test_the_trail_goes_back_to_any_earlier_question(panel, label):
    panel.restart()
    panel.answer("have", G.HAVE_LABEL)
    panel.answer("space", ORGANISM)
    panel.answer("label", label)
    panel.goto("label")
    assert panel.question.text() == G.STEPS["label"].question
    assert "label" not in panel.answers
    panel.back()
    assert panel.question.text() == G.STEPS["space"].question
    panel.restart()
    assert panel.answers == {} and panel.question.text() == G.STEPS["have"].question


def test_the_gene_list_step_counts_what_this_table_holds(panel, ctx):
    panel.restart()
    panel.answer("have", G.HAVE_SET)
    panel.answer("space", ORGANISM)
    panel.set_genes([str(ctx.gene_ids[0]), "NOT_A_GENE"])
    assert "1 genes found" in panel.count.text() and "NOT_A_GENE" in panel.count.text()
    assert panel.use_gated() == 12
    assert panel.use_example() > 0
    panel._accept()
    assert G.answered(panel.answers, "genes")


def test_a_gene_is_searched_for_by_product(panel):
    panel.restart()
    panel.answer("have", G.HAVE_GENE)
    panel.answer("space", ORGANISM)
    panel.search.setText("kinase")
    assert panel.matches.count() > 5
    panel._pick_gene(panel.matches.item(0))
    assert G.answered(panel.answers, "gene")
    assert panel.question.text() == G.STEPS["goal"].question


def test_the_panel_asks_for_a_strategy_rather_than_running_it_itself(panel):
    panel.restart()
    panel.answer("have", G.HAVE_NOTHING)
    panel.answer("space", ORGANISM)
    asked = []
    panel.open_strategy.connect(lambda key, settings: asked.append(key))
    panel.cards[0].open_btn.click()
    assert asked == [panel.recommendations()[0].key]
    views = []
    panel.open_view.connect(views.append)
    for row in panel.body.findChildren(QtWidgets.QPushButton):
        if row.text().startswith("Open the"):
            row.click()
    assert views


def test_the_window_gets_the_tab_and_fills_a_strategy_in_from_the_answers(label):
    """The whole hook: the dock is there, beside Strategies, and Open in Strategies fills it in."""
    from starplast import app as A
    window = A.Window()
    try:
        panel = window.guided
        assert panel is not None and window.guided_dock.windowTitle() == "start here"
        assert window.guided_dock in window.tabifiedDockWidgets(window.strategies_dock)
        panel.answer("have", G.HAVE_LABEL)
        panel.answer("space", window.strategy_panel.ctx.organism)
        panel.answer("label", label)
        panel.answer("goal", G.PREDICT)
        best = panel.recommendations()[0]
        panel.cards[0].open_btn.click()
        assert window.strategy_panel.current.key == best.key
        assert window.strategy_panel.settings()["target"] == label
    finally:
        window.console.remove()
        window.close()


def test_the_help_search_can_find_the_tab_and_its_questions():
    from starplast import help_index
    index = help_index.build_index()
    hits = help_index.search(index, "start here")
    assert any(e.title == "Start here" for e in hits)
    assert any(e.data.get("dock") == "start here"
               for e in help_index.search(index, "find more genes like mine"))
