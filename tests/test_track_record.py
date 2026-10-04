"""The track record: every labelled gene held out exactly once, and the arithmetic of the levels.

The properties worth holding are the ones that make the record trustworthy rather than merely
present: complete coverage, folds that never split an orthogroup, abstentions kept apart from
errors, and no percentage printed from too few genes to mean one.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import strategies as S  # noqa: E402
from starplast import track_record as T  # noqa: E402


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture(scope="module")
def ledger(planted):
    return T.evaluate(planted, "feature_knn", "compartment")


def test_every_labelled_gene_is_held_out_exactly_once(planted, ledger):
    truth = planted.truth("compartment")
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    expected = set(np.flatnonzero(kept.notna().to_numpy()).tolist())
    assert set(ledger["gene"]) == expected, "a labelled gene was missed or scored twice"
    assert not ledger["gene"].duplicated().any()


def test_a_fold_never_splits_an_orthogroup(planted):
    """A gene recovered through a paralog that stayed visible is not the inference being tested."""
    truth = planted.truth("compartment")
    labelled = np.flatnonzero(truth.notna().to_numpy())
    folds = T.folds_by_group(planted, labelled)
    groups = np.asarray(planted.groups())[labelled]
    per_group = pd.DataFrame({"group": groups, "fold": folds}).groupby("group")["fold"].nunique()
    assert (per_group == 1).all(), "an orthogroup was split across folds"


def test_the_folds_are_balanced_enough_to_be_useful(planted):
    truth = planted.truth("compartment")
    labelled = np.flatnonzero(truth.notna().to_numpy())
    counts = pd.Series(T.folds_by_group(planted, labelled)).value_counts()
    assert len(counts) == T.FOLDS
    assert counts.max() <= 3 * counts.min(), "one fold swallowed the genes"


def test_an_abstention_is_neither_right_nor_wrong(ledger):
    for _, row in ledger.iterrows():
        if row["abstained"]:
            assert row["prediction"] is None and row["correct"] is None
        else:
            assert isinstance(row["prediction"], str) and row["correct"] in (True, False)
    pooled = T.summary(ledger, "strategy").iloc[0]
    assert pooled["answered"] == pooled["genes"] - pooled["abstained"]
    assert pooled["right"] + pooled["wrong"] == pooled["answered"]


def test_the_rate_is_over_answered_genes_and_carries_its_interval(ledger):
    row = T.summary(ledger, "strategy").iloc[0]
    assert row["rate"] == pytest.approx(row["right"] / row["answered"])
    assert row["rate_low"] <= row["rate"] <= row["rate_high"]
    assert bool(row["enough"]) == bool(row["answered"] >= T.MIN_FOR_RATE)


def test_too_few_answered_genes_is_said_rather_than_shown_as_a_percentage():
    rows = [{"organism": "Tg", "strategy": "s", "target": "t", "setting_key": "", "seed": 0,
             "fold": 0, "mode": "together", "gene": i, "gene_id": f"g{i}", "truth": "a",
             "prediction": "a" if i < 2 else None, "correct": True if i < 2 else None,
             "abstained": i >= 2, "support": 1.0} for i in range(6)]
    row = T.summary(pd.DataFrame(rows), "strategy").iloc[0]
    assert row["answered"] == 2 and not bool(row["enough"])


def test_every_level_pools_the_same_genes(ledger):
    for level in ("gene", "class", "target", "strategy"):
        out = T.summary(ledger, level)
        assert len(out), level
        assert out["genes"].sum() == len(ledger), level
        assert {"rate", "rate_low", "rate_high", "enough"} <= set(out.columns), level


def test_the_class_level_says_what_a_label_is_confused_with(planted):
    led = T.evaluate(planted, "feature_knn", "compartment")
    wrong = led[led["correct"] == False]                        # noqa: E712
    classes = T.summary(led, "class")
    if len(wrong):
        row = classes[classes["truth"] == wrong["truth"].iloc[0]].iloc[0]
        assert row["confused_with"], "a class with errors must name what it was called instead"


def test_a_gene_can_be_asked_about_by_name(planted, ledger):
    gene = ledger["gene_id"].iloc[0]
    one = T.for_gene(ledger, gene)
    assert len(one) == 1 and set(one["strategy"]) == {"feature_knn"}
    line = T.sentence(ledger, gene)
    assert line.startswith("Known ") and "recovered by" in line or "abstained" in line


def test_a_strategy_that_cannot_speak_about_a_target_returns_nothing_rather_than_guessing(planted):
    """A refusal is a fact about the pairing; it must not become an empty-looking success."""
    out = T.evaluate(planted, "holdout_search", "compartment")      # a map walk, not a gene caller
    assert out.empty and list(out.columns) == list(T.COLUMNS)


def test_the_supported_strategies_are_real_and_call_labels():
    keys = T.supported()
    assert len(keys) >= 10
    known = {s.key: s for s in S.catalog()}
    for key in keys:
        assert key in known, key
        assert known[key].task == "label calls", key


# --------------------------------------------------------------------------- sets held out together
def test_a_set_is_hidden_all_at_once(planted):
    """Every member of the set must be invisible while any of them is predicted."""
    sets = T.class_sets(planted, "compartment")
    name, members = next(iter(sets.items()))
    led = T.evaluate_sets(planted, "feature_knn", "compartment", {name: members})
    assert len(led) == len(members)
    assert set(led["mode"]) == {"set"} and set(led["set_name"]) == {name}
    assert set(led["set_size"]) == {len(members)}, "the rows disagree about how many were hidden"


def test_hiding_a_whole_class_leaves_nothing_to_copy(planted):
    """A vote among visible labels cannot name a class when every example of it is gone -- so the
    rate is zero, and what the genes are called INSTEAD is the finding."""
    sets = T.class_sets(planted, "compartment")
    led = T.evaluate_sets(planted, "feature_knn", "compartment", sets)
    summary = T.set_summary(led)
    assert (summary["right"] == 0).all()
    assert summary["placed_at"].str.len().gt(0).any()
    assert ((summary["together"] >= 0) & (summary["together"] <= 1)).all()


def test_the_degradation_curve_covers_the_sizes_asked_for(planted):
    sets = T.random_sets(planted, "compartment", sizes=(1, 5, 20), repeats=3, seed=2)
    assert len(sets) == 9
    assert sorted({len(v) for v in sets.values()}) == [1, 5, 20]
    led = T.evaluate_sets(planted, "feature_knn", "compartment", sets)
    assert len(led) == 3 * (1 + 5 + 20)


def test_set_rows_and_fold_rows_live_in_one_table(planted, ledger):
    sets = T.random_sets(planted, "compartment", sizes=(5,), repeats=2, seed=3)
    both = pd.concat([ledger, T.evaluate_sets(planted, "feature_knn", "compartment", sets)],
                     ignore_index=True)
    assert set(both["mode"]) == {"together", "set"}
    assert len(T.summary(both, "strategy")) == 1, "the levels pool across both kinds of hold-out"
    assert len(T.set_summary(both)) == 2, "and the set level reports only the sets"


# --------------------------------------------------------------------------- the shipped record
@pytest.fixture(scope="module")
def built():
    return T.shipped()


def test_the_shipped_record_covers_both_organisms_and_every_supported_strategy(built):
    if not len(built):
        pytest.skip("the track record has not been built on this machine")
    folds = built[built["mode"] == "together"]
    for organism in built["organism"].unique():
        here = folds[folds["organism"] == organism]
        assert here["strategy"].nunique() >= 10, organism
        # Every labelled gene of each label appears exactly once per strategy.
        per = here.groupby(["target", "strategy"], observed=True)["gene"].agg(["count", "nunique"])
        per = per[per["count"] > 0]
        assert (per["count"] == per["nunique"]).all(), organism


def test_the_shipped_record_kept_the_set_holdouts(built):
    if not len(built):
        pytest.skip("not built here")
    sets = built[built["mode"] == "set"]
    assert len(sets), "no set hold-outs were kept"
    names = set(sets["set_name"].astype(str))
    assert any(n.startswith("random ") for n in names), "the degradation curve is missing"
    assert T.set_summary(built)["together"].notna().any()


def test_a_gene_card_section_is_html_or_nothing(built):
    if not len(built):
        pytest.skip("not built here")
    organism = str(built["organism"].iloc[0])
    gene = str(built[built["mode"] == "together"]["gene_id"].iloc[0])
    html = T.gene_html(gene, organism)
    assert "If this gene were unknown" in html and "<table" in html
    assert T.gene_html("no such gene", organism) == "", "an unknown gene adds no section"


def test_the_class_level_names_the_best_strategy_and_the_confusions(built):
    if not len(built):
        pytest.skip("not built here")
    folds = built[built["mode"] == "together"]
    organism = str(folds["organism"].iloc[0])
    target = T.default_target(organism)
    folds = folds[(folds["organism"] == organism) & (folds["target"].astype(str) == target)]
    label = str(folds["truth"].astype(str).value_counts().index[0])
    html = T.class_html(target, label, organism)
    assert label in html and "<table" in html and "rate [95%]" in html
    assert T.class_html(target, "no such class", organism) == ""


def test_the_category_level_lists_every_class_once_and_links_down(built):
    if not len(built):
        pytest.skip("not built here")
    target = T.default_target("Tg")
    folds = built[(built["mode"] == "together") & (built["organism"] == "Tg")
                  & (built["target"].astype(str) == target)]
    html = T.target_html(target, "Tg")
    classes = set(folds["truth"].astype(str))
    assert html.count("starplast://class/") == len(classes), "a class is missing or repeated"
    assert "best recovered" in html
    assert "per class" in html and "commonest class" in html, "no strategy-side view"
    assert T.target_html("no such category", "Tg") == ""

def test_a_strategy_line_names_where_it_is_weakest_but_only_where_judgeable(built):
    if not len(built):
        pytest.skip("not built here")
    organism = str(built["organism"].iloc[0])
    strategy = str(built[built["mode"] == "together"]["strategy"].iloc[0])
    line = T.weakest(strategy, organism)
    assert "held-out genes" in line
    # Nothing is called weakest on a count too small to judge.
    rows = T.summary(built[(built["mode"] == "together") & (built["strategy"] == strategy)
                           & (built["organism"] == organism)
                           & (built["target"].astype(str) == T.default_target(organism))], "class")
    for row in rows[~rows["enough"].astype(bool)].itertuples():
        assert f"{row.truth} (" not in line, row.truth


def test_the_evidence_panel_drills_from_a_gene_to_its_class_and_back(built):
    """The click path the user asked for: gene, then the class, then back to the gene."""
    if not len(built):
        pytest.skip("not built here")
    from PyQt6 import QtCore, QtWidgets
    _app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast import app as A
    w = A.Window()
    ids = w.nodes.gene_id.astype(str)
    folds = built[(built["mode"] == "together") & (built["organism"] == "Tg")
                  & (built["target"].astype(str) == T.default_target("Tg"))]
    row = int(ids[ids == str(folds["gene_id"].iloc[0])].index[0])
    w.show_detail(row)
    html = w.detail.toHtml()
    assert "If this gene were unknown" in html and "starplast://class/" in html
    link = html.split("starplast://class/", 1)[1].split('"', 1)[0]
    w._detail_link(QtCore.QUrl("starplast://class/" + link))
    assert "rate [95%]" in w.detail.toHtml(), "the class level did not open"
    assert "back to the gene" in w.detail.toHtml()
    html = w.detail.toHtml()
    assert "starplast://target/" in html, "the class page has no way up to its category"
    up = html.split("starplast://target/", 1)[1].split('"', 1)[0]
    w._detail_link(QtCore.QUrl("starplast://target/" + up))
    assert "best recovered" in w.detail.toHtml(), "the category level did not open"
    w._detail_link(QtCore.QUrl(f"starplast://gene/{row}"))
    assert "If this gene were unknown" in w.detail.toHtml(), "back did not return to the gene"
    target = T.default_target("Tg")
    w._record_link(f"starplast://target/{target}")
    assert "best recovered" in w.detail.toHtml(), "a strategy card's link did not open the category"


def test_the_strategy_card_shows_where_it_is_weakest(built):
    if not len(built):
        pytest.skip("not built here")
    from PyQt6 import QtWidgets
    _app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast.strategy_card import StrategyCard
    card = StrategyCard()
    s = S.get("feature_knn")
    card.show_strategy(s, "Tg", {}, {}, "reliable", "")
    assert card.record.isVisibleTo(card) and "weakest on" in card.record.text()
    assert "starplast://target/" in card.record.text(), "the line does not lead to the classes"
    heard = []
    card.record_link.connect(heard.append)
    card.record.linkActivated.emit("starplast://target/compartment")
    assert heard == ["starplast://target/compartment"]
    card.show_strategy(S.get("holdout_search"), "Tg", {}, {}, "weak", "")
    assert not card.record.isVisibleTo(card), "a map walk calls no gene, so it has no record"


def test_your_own_list_is_hidden_together_and_typed_names_work(planted):
    """A pasted list resolves by name, every strategy asked hides all of it, and unknowns are dropped."""
    names = planted.nodes["gene_id"].iloc[:6].tolist()
    led = T.my_list(planted, names + ["NOT_A_GENE"], "compartment",
                    strategies=["feature_knn", "layer_vote"])
    assert set(led["strategy"]) <= {"feature_knn", "layer_vote"} and len(led)
    assert set(led["set_name"]) == {"your list"}
    labelled = led.groupby("strategy")["set_size"].first()
    assert (labelled <= 6).all(), "a name that resolved to nothing was counted as hidden"
    assert set(led["gene_id"]) <= set(names)


def test_a_confusion_is_only_named_if_it_happened(built):
    """The shipped labels are categorical; zero-count categories must never be named as confusions."""
    if not len(built):
        pytest.skip("not built here")
    folds = built[built["mode"] == "together"]
    for row in T.summary(folds, "class").itertuples():
        named = [c for c in str(row.confused_with or "").split(", ") if c]
        if not named:
            continue
        here = folds[(folds["truth"] == row.truth) & (folds["strategy"] == row.strategy)
                     & (folds["organism"] == row.organism) & (folds["correct"] == False)]  # noqa: E712
        assert set(named) <= set(here["prediction"].dropna().astype(str)), (row.truth, named)
    sets = T.set_summary(built)
    for row in sets.itertuples():
        if row.together == 1.0:
            assert ", " not in row.placed_at, f"{row.set_name}: one place, but two named"


def test_one_gene_alone_hides_its_orthogroup_and_is_cached(planted, tmp_path, monkeypatch):
    from starplast import paths
    monkeypatch.setattr(paths, "user_cache_dir", lambda: str(tmp_path))
    truth = planted.truth("compartment")
    gene = int(np.flatnonzero(truth.notna().to_numpy())[0])
    rows = T.alone(planted, gene, "compartment", strategies=("feature_knn", "layer_vote"))
    assert set(rows["mode"]) == {"alone"} and set(rows["gene"]) == {gene}
    assert len(rows) == 2, "one row per strategy asked"
    groups = np.asarray(planted.groups())
    siblings = int((groups == groups[gene]).sum())
    labelled = int(truth.iloc[np.flatnonzero(groups == groups[gene])].notna().sum())
    assert set(rows["set_size"]) == {labelled}, f"the orthogroup ({siblings}) was not hidden with it"
    cached = list((tmp_path / "track_record").glob("alone_*.parquet"))
    assert len(cached) == 1
    again = T.alone(planted, planted.gene_ids[gene], "compartment",
                    strategies=("feature_knn",))
    assert len(again) == 1 and list(again["strategy"]) == ["feature_knn"], "the cache was not reused"


def test_the_gene_card_tests_another_label_on_demand(built, tmp_path, monkeypatch):
    """A link per other label; clicking one runs in the background and shows the answer in place."""
    if not len(built):
        pytest.skip("not built here")
    from PyQt6 import QtCore, QtWidgets
    from starplast import paths
    monkeypatch.setattr(paths, "user_cache_dir", lambda: str(tmp_path))
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast import app as A
    w = A.Window()
    if getattr(w, "strategy_panel", None) is None:
        pytest.skip("no strategy panel here")
    ids = w.nodes.gene_id.astype(str)
    folds = built[(built["mode"] == "together") & (built["organism"] == "Tg")]
    row = int(ids[ids == str(folds["gene_id"].iloc[0])].index[0])
    w.show_detail(row)
    html = w.detail.toHtml()
    if "starplast://record/" in html:
        link = "starplast://record/" + html.split("starplast://record/", 1)[1].split('"', 1)[0]
        w._detail_link(QtCore.QUrl(link))
        shown = w.detail.toHtml()
        assert "were unknown" in shown and "back to the gene" in shown, "a recorded label did not open"
        w.show_detail(row)
        html = w.detail.toHtml()
    if "starplast://alone/" in html:
        link = "starplast://alone/" + html.split("starplast://alone/", 1)[1].split('"', 1)[0]
    else:
        # Every label this gene has is recorded; the on-demand path still serves any other.
        other = next(t for t in T.recorded_targets("Tg") if t != T.default_target("Tg"))
        link = f"starplast://alone/{row}/{other}"
    w._detail_link(QtCore.QUrl(link))
    assert "Hiding" in w.detail.toHtml(), "nothing said while it runs"
    deadline = QtCore.QDeadlineTimer(120_000)
    while w.jobs.busy and not deadline.hasExpired():
        app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 100)
    app.processEvents()
    done = w.detail.toHtml()
    assert "were unknown" in done or "No strategy can speak" in done, done[-400:]
    assert "back to the gene" in done


def test_a_recommendation_quotes_the_record_only_for_its_own_label(built):
    if not len(built):
        pytest.skip("not built here")
    folds = built[(built["mode"] == "together") & (built["organism"] == "Tg")]
    target = T.default_target("Tg")
    phrase = T.record_phrase("feature_knn", "Tg", target)
    assert phrase.startswith(f"with {target.replace('_', ' ')} hidden, right on")
    assert T.record_phrase("feature_knn", "Tg", "some_other_label") == ""
    assert T.record_phrase("holdout_search", "Tg") == ""


def test_start_here_puts_a_strategy_below_its_baseline_last(built):
    if not len(built):
        pytest.skip("not built here")
    from starplast import guided as G
    ctx = S.Context.shipped("Tg")
    target = T.default_target("Tg")
    for goal in G.GOALS:
        recs = G.recommend({"subject": "label", "label": target, "goal": goal}, ctx)
        flags = [T.beats_baseline(r.key, "Tg", target) is False for r in recs]
        assert flags == sorted(flags), f"{goal}: a below-baseline strategy outranks a better one"


def test_start_here_tests_your_own_genes_and_shows_the_grid(built, tmp_path, monkeypatch):
    """From Start here, with a gene list: the button runs the hold-out and the grid opens."""
    if not len(built):
        pytest.skip("not built here")
    from PyQt6 import QtCore, QtWidgets
    from starplast import guided as G
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast import app as A
    window = A.Window()
    try:
        panel = window.guided
        if panel is None:
            pytest.skip("no Start-here tab here")
        folds = built[(built["mode"] == "together") & (built["organism"] == "Tg")]
        genes = folds["gene_id"].astype(str).drop_duplicates().head(6).tolist()
        panel.answer("have", G.HAVE_SET)
        panel.answer("space", window.strategy_panel.ctx.organism)
        panel.set_genes(genes)
        panel._accept()
        panel.answer("goal", G.PREDICT)
        assert getattr(panel, "test_button", None) is not None, "no way to test the list"
        panel.test_button.click()
        assert "Hiding your 6 genes" in window.detail.toHtml()
        deadline = QtCore.QDeadlineTimer(180_000)
        while window.jobs.busy and not deadline.hasExpired():
            app.processEvents(QtCore.QEventLoop.ProcessEventsFlag.AllEvents, 100)
        app.processEvents()
        html = window.detail.toHtml()
        assert "Would they have found your genes?" in html
        assert html.count("starplast://gene/") == 6, "every gene should link to its card"
    finally:
        window.console.remove()
        window.close()


def test_the_list_grid_marks_every_gene_for_every_strategy(planted):
    names = planted.nodes["gene_id"].iloc[:8].tolist()
    rows = T.my_list(planted, names, "compartment", strategies=["feature_knn", "layer_vote"])
    html = T.list_html(rows)
    genes = rows["gene_id"].nunique()
    assert html.count("✓") + html.count("✗") + html.count("<td align='center' style='color:#888888'>·") \
        >= genes * rows["strategy"].nunique()
    assert "None of your genes" in T.list_html(rows.iloc[0:0])


def test_a_strategy_that_favours_small_classes_is_not_called_worse_than_guessing():
    """Losing on accuracy but winning per class is a trade-off, not a failure -- and both are said."""
    big, small = ["big"] * 90, ["small"] * 10
    truth = big + small
    # Right on every small-class gene, wrong on most big ones: accuracy 0.30 < 0.90, per class 0.61.
    prediction = ["big"] * 20 + ["small"] * 70 + ["small"] * 10
    led = pd.DataFrame({"truth": truth, "prediction": prediction,
                        "correct": [t == p for t, p in zip(truth, prediction)],
                        "abstained": False})
    b = T._baselines(led)
    assert b["accuracy"] < b["commonest"] and b["balanced"] > b["chance"]
    # And one that is worse on both: always the wrong class.
    wrong = led.assign(prediction="other", correct=False)
    w = T._baselines(wrong)
    assert w["accuracy"] <= w["commonest"] and w["balanced"] <= w["chance"]


def test_a_gene_list_setting_can_be_tested_from_the_strategies_tab(built):
    """Every gene-list box carries Test these, and it hands the resolved genes to the window."""
    if not len(built):
        pytest.skip("not built here")
    from PyQt6 import QtWidgets
    _app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from starplast import app as A
    window = A.Window()
    try:
        panel = window.strategy_panel
        key = next(s.key for s in S.catalog() if any(p.kind == "genes" for p in s.params))
        panel.select(key)
        host = next(w for _p, w in panel.inputs.values() if hasattr(w, "test"))
        genes = built[built["organism"] == "Tg"]["gene_id"].astype(str).drop_duplicates().head(3)
        host.box.setPlainText("\n".join(genes) + "\nNOT_A_GENE")
        asked = []
        panel.test_genes.disconnect()
        panel.test_genes.connect(asked.append)
        host.test.click()
        assert asked and sorted(asked[0]) == sorted(genes), asked
        host.box.setPlainText("NOT_A_GENE")
        assert panel.test_list(host.box) == 0
    finally:
        window.console.remove()
        window.close()



def test_the_recorded_labels_are_biology_not_bookkeeping():
    ctx = S.Context.shipped("Tg")
    chosen = T.labels(ctx)
    assert chosen[0] == S.default_category(ctx)
    for bookkeeping in ("compartment_source", "ortholopit_donors", "screen_scorers_agree"):
        assert bookkeeping not in chosen



def test_the_record_never_walks_a_layer_the_strategy_would_refuse(planted):
    """A banned default layer is replaced, and the ledger says which layer was walked instead."""
    target = "compartment"
    banned = planted.banned_layers(target)
    if not banned:
        pytest.skip("the planted table bans no layer for this target")
    layer = sorted(banned)[0]
    chosen = T._permitted(planted, {"layer": layer, "restart": 0.5}, target)
    assert chosen["layer"] not in banned
    led = T.evaluate(planted, "layer_propagation", target, settings={"layer": layer})
    if len(led):
        assert f"layer={layer}," not in str(led["setting_key"].iloc[0])
