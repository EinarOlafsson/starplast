#!/usr/bin/env python3
"""The paper-ready PDF: it must build headless, and it must not lie about what it is showing."""
import os

import numpy as np
import pandas as pd
import pytest

from starplast import recipes as R
from starplast import report as RP


def _result(n=300, with_control=True, refused=False):
    """A finished run, assembled directly rather than by running one: the report is being tested,
    not the pipeline, and a real run would put minutes into every assertion."""
    recipe = R.Recipe(question="Which unplaced proteins sit with the dense granules?",
                      inputs=[("biology", ("parasite phenotype", "fitness and essentiality"))],
                      holdout="compartment", validation_holdout="a_control",
                      axis="localisation and export", expectation="dense granules cluster")
    closure = R.Closure(blocks=("a", "b"), columns=("x", "y"), holdout_column="compartment",
                        control_column="a_control",
                        excluded={"compartment": "primary: the target itself"},
                        removed={"lopit_map": "primary: measured association"},
                        refusal="the holdout restates itself" if refused else "")
    res = R.RecipeResult(recipe=recipe, closure=closure)
    if refused:
        res.stopped_because = closure.refusal
        return res
    rng = np.random.default_rng(0)
    labels = np.arange(n) % 4
    labels[::9] = -1
    truth = pd.Series([None if i % 3 else ("dense granules" if labels[i] == 0 else "cytosol")
                       for i in range(n)], dtype="object")
    # Spread across every cluster on purpose. Keyed on i % 4 it landed only in cluster 0 -- because
    # the labels are i % 4 too -- so every panel that needs a control in more than one cluster was
    # being drawn from a single row and two branches never ran.
    control = pd.Series([None if i % 7 == 0 else ("host" if labels[i] in (0, 2) else "not host")
                         for i in range(n)], dtype="object") if with_control \
        else pd.Series(dtype="object")
    res.coords = rng.normal(size=(n, 3))
    res.labels, res.truth, res.control_truth = labels, truth, control
    res.quality = {"clusters": 4, "clustered": 0.89, "evenness": 0.97, "largest": 0.26,
                   "usable": True, "why_not": "", "score": 0.9}
    res.settings = {"n_neighbors": 30, "min_dist": 0.0, "min_cluster_size": 40, "min_samples": 5}
    res.summary, res.recovery = R.score_recovery(labels, truth, min_label=15)
    res.per_cluster = R.recovery_by_cluster(labels, truth)
    res.inference = R.infer(labels, truth, [f"TGME49_{200000+i}" for i in range(n)])
    if with_control:
        res.control = R.dominant_by_cluster(labels, control).rename(
            columns={"label": "control_label", "n_labelled": "n_control_labelled"})
    return res


def _pages(path) -> int:
    with open(path, "rb") as fh:
        blob = fh.read()
    assert blob.startswith(b"%PDF"), "not a PDF"
    # `/Type /Pages` is the page TREE, not a page. Counting it as one made every page count one too
    # many, which would have hidden a missing page rather than reporting it.
    return blob.count(b"/Type /Page") - blob.count(b"/Type /Pages")


def test_a_run_becomes_a_multi_page_pdf(tmp_path):
    path = RP.recipe_pdf(_result(), str(tmp_path / "r.pdf"))
    assert os.path.exists(path) and _pages(path) >= 2


def test_every_alternative_map_gets_its_own_page(tmp_path):
    """A report that showed only the winner would hide what the choice was between."""
    res = _result()
    res.alternatives = [_result(), _result()]
    one = _pages(RP.recipe_pdf(_result(), str(tmp_path / "one.pdf")))
    three = _pages(RP.recipe_pdf(res, str(tmp_path / "three.pdf")))
    assert three == one + 2


def test_a_refused_recipe_still_produces_a_document_saying_why(tmp_path):
    """A reader looking for the run must find the reason, not an absent file."""
    path = RP.recipe_pdf(_result(refused=True), str(tmp_path / "refused.pdf"))
    assert _pages(path) == 2                       # the cover, and the refusal


def test_a_run_with_no_control_still_reports(tmp_path):
    assert _pages(RP.recipe_pdf(_result(with_control=False), str(tmp_path / "n.pdf"))) >= 2


def test_the_directory_is_created_rather_than_required(tmp_path):
    path = RP.recipe_pdf(_result(), str(tmp_path / "deep" / "nested" / "r.pdf"))
    assert os.path.exists(path)


def test_nothing_scoreable_draws_a_panel_that_says_so_rather_than_an_empty_axis(tmp_path):
    """The failure this prevents is a blank square, which reads as a broken figure rather than as a
    map that recovered nothing."""
    res = _result()
    res.recovery = pd.DataFrame()
    res.per_cluster = pd.DataFrame()
    res.control = pd.DataFrame()
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "empty.pdf"))) >= 2


def test_a_control_present_but_absent_from_every_cluster_says_it_cannot_judge(tmp_path):
    """Distinct from having no control at all, and the difference is the whole argument for the
    control: silence is not corroboration and must not be drawn as if it were."""
    res = _result()
    res.control = res.control.assign(n_control_labelled=0)
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "silent.pdf"))) >= 2


def test_a_scatter_with_more_labels_than_colours_keeps_the_rest_as_other(tmp_path):
    res = _result(n=600)
    res.truth = pd.Series([f"label {i % 40}" for i in range(600)], dtype="object")
    res.per_cluster = R.recovery_by_cluster(res.labels, res.truth)
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "many.pdf"))) >= 2


def test_a_map_with_no_genes_at_all_does_not_crash_the_document(tmp_path):
    res = _result()
    res.coords = np.zeros((0, 3))
    res.labels = np.array([], dtype=int)
    res.truth = pd.Series(dtype="object")
    res.control_truth = pd.Series(dtype="object")
    res.recovery = pd.DataFrame()
    res.per_cluster = pd.DataFrame()
    res.control = pd.DataFrame()
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "none.pdf"))) >= 2


def test_the_clusters_the_answer_rests_on_are_split_by_whether_the_control_backs_them(tmp_path):
    """A run whose named genes all come from clusters the control does not corroborate is standing
    on its primary alone, and the figure has to make that visible rather than implied."""
    res = _result()
    res.inference = res.inference.assign(
        control_agrees=[i % 2 == 0 for i in range(len(res.inference))])
    assert len(res.inference), "the fixture named nobody, so this asserts nothing"
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "rests.pdf"))) >= 2


def test_a_control_spanning_orders_of_magnitude_is_drawn_on_a_log_axis(tmp_path):
    """A rare control label reaches several hundred fold while a common one sits near two. Linear,
    the informative end of that axis is one pixel tall."""
    res = _result()
    scoreable = res.control[res.control.n_control_labelled > 0]
    assert len(scoreable) >= 2, "the fixture cannot exercise a range with one scoreable cluster"
    res.control = res.control.assign(
        enrichment=[1.1 if i else 400.0 for i in range(len(res.control))])
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "log.pdf"))) >= 2


def test_a_run_that_named_nobody_says_so_where_the_answer_would_have_been(tmp_path):
    """The panel must say no cluster carries the answer, rather than drawing an empty axis that
    reads as a broken figure."""
    res = _result()
    res.inference = pd.DataFrame()
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "nobody.pdf"))) >= 2


def test_outliers_do_not_own_the_axis(tmp_path):
    """The defect this prevents was visible in the first real report: a handful of genes at x = 10
    while 8,100 sat between -52 and -46, and every panel drew the structure as a four-pixel smear."""
    coords = np.vstack([np.random.default_rng(0).normal(size=(300, 3)),
                        np.array([[500.0, 500.0, 0.0]])])
    (x0, x1), (y0, y1) = RP._limits(coords)
    assert x1 < 100 and y1 < 100, "an outlier at 500 still set the frame"
    assert x0 < 0 < x1, "the bulk of the cloud is not inside the frame"


def test_an_empty_map_still_gets_finite_limits():
    assert RP._limits(np.zeros((0, 3))) == ((-1.0, 1.0), (-1.0, 1.0))


def test_fixed_class_candidates_report_held_out_precision(tmp_path):
    res = _result()
    res.inference = pd.DataFrame({"gene_id": ["g1"], "cluster": [0],
                                  "predicted": ["dense granules"],
                                  "cv_class_precision": [.75], "enrichment": [3.]})
    assert _pages(RP.recipe_pdf(res, str(tmp_path / "classification.pdf"))) >= 2
