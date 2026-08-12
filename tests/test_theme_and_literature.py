#!/usr/bin/env python3
"""Theme resolution, colour-map choice, and the literature layer's remaining edges.

The colour-map tests carry a real constraint: pyqtgraph ships no `coolwarm`, and asking for a name it
does not have raises FileNotFoundError from inside a paint call -- a crash with no useful traceback,
during rendering, from a menu selection. So resolution falls back rather than raising, and the offered
names are only ones the library can actually build.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import literature as LIT  # noqa: E402
from starplast import theme as TH  # noqa: E402


# --------------------------------------------------------------------------- palettes
def test_the_background_is_the_page_colour_not_a_literal():
    """So a new theme changes the viewport with everything else instead of leaving a black hole."""
    for t in TH.THEMES:
        assert TH.background(t)[:3] == TH.rgbf(TH.palette_for(t)["bg"])[:3]
        assert all(0.0 <= c <= 1.0 for c in TH.background(t)[:3])


def test_unknown_has_its_own_grey_derived_from_the_theme():
    """"Unknown" is not "absent" and not a compartment; it needs a colour that reads on either ground."""
    for t in TH.THEMES:
        assert TH.unknown_colour(t)[:3] == TH.rgbf(TH.palette_for(t)["fg_dim"])[:3]


def test_the_stylesheet_is_produced_for_every_theme():
    for t in TH.THEMES:
        css = TH.stylesheet(t)
        assert isinstance(css, str) and css.strip()


def test_dark_and_light_themes_are_distinguished():
    assert TH.is_dark("dark") and not TH.is_dark("light")


def test_rgbf_parses_hex_and_applies_alpha():
    assert TH.rgbf("#ffffff") == pytest.approx((1.0, 1.0, 1.0, 1.0))
    assert TH.rgbf("#000000", 0.5)[3] == pytest.approx(0.5)


def test_every_offered_colour_map_can_actually_be_built():
    """A name pyqtgraph does not have raises FileNotFoundError inside a paint call. Offering
    matplotlib's `coolwarm` -- which pyqtgraph does not ship -- did exactly that."""
    for name in TH.CMAPS:
        assert TH.resolve_cmap(name) is not None, name


def test_an_unknown_colour_map_falls_back_instead_of_raising():
    assert TH.resolve_cmap("no_such_colormap_anywhere") is not None


def test_maps_can_be_listed_by_family():
    for kind in ("sequential", "diverging", "categorical"):
        for n in TH.cmaps_of(kind):
            assert TH.CMAPS[n][0] == kind


# --------------------------------------------------------------------------- choosing a family
def test_labels_get_a_categorical_map():
    assert TH.kind_for_column(pd.Series(["rhoptry", "micronemes", "dense granule"])) == "categorical"


def test_an_all_missing_column_is_treated_as_categorical():
    assert TH.kind_for_column(pd.Series([np.nan, np.nan, np.nan])) == "categorical"


def test_an_encoded_class_is_categorical_but_a_coarse_continuous_score_is_not():
    """Few distinct values is not enough on its own: a fitness score can take five values in a small
    selection and is still continuous. Categorical means few distinct values that are all whole
    numbers."""
    assert TH.kind_for_column(pd.Series([0, 1, 2, 1, 0])) == "categorical"
    assert TH.kind_for_column(pd.Series([-1.5, -0.5, 0.5, 1.5, 2.5])) != "categorical"


def test_a_quantity_straddling_zero_gets_a_diverging_map():
    assert TH.kind_for_column(pd.Series([-2.0, -0.5, 0.3, 1.8])) == "diverging"


def test_a_strictly_positive_quantity_does_not():
    """A diverging ramp on a positive quantity invents a midpoint that is not in the data."""
    assert TH.kind_for_column(pd.Series([0.1, 5.0, 90.5])) == "sequential"


def test_a_barely_negative_tail_does_not_trigger_diverging():
    assert TH.kind_for_column(pd.Series([-0.01, 5.0, 90.5])) == "sequential"


def test_a_mostly_text_column_is_categorical_even_with_some_numbers():
    assert TH.kind_for_column(pd.Series(["a", "b", "c", 1, 2])) == "categorical"


# --------------------------------------------------------------------------- categorical colours
def test_categorical_colours_are_distinct_and_in_range():
    cols = TH.categorical_colours(8, theme="dark")
    assert len(cols) == 8
    assert len({tuple(round(c, 3) for c in col) for col in cols}) > 1
    assert all(0.0 <= c <= 1.0 for col in cols for c in col)


def test_categorical_colours_are_lightened_on_dark_and_darkened_on_light():
    """The same hue at the same value is legible on one ground and not the other."""
    dark = TH.categorical_colours(6, theme="dark")
    light = TH.categorical_colours(6, theme="light")
    assert sum(sum(c) for c in dark) > sum(sum(c) for c in light)


def test_colours_are_produced_even_without_matplotlib(monkeypatch):
    """matplotlib is not a hard dependency of the viewer, and a missing palette must not stop it
    drawing."""
    import builtins
    real = builtins.__import__

    def no_mpl(name, *a, **k):
        if name.startswith("matplotlib"):
            raise ImportError("no matplotlib")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_mpl)
    cols = TH.categorical_colours(5, theme="dark")
    assert len(cols) == 5
    assert all(0.0 <= c <= 1.0 for col in cols for c in col)


def test_a_named_map_can_seed_the_categorical_palette():
    assert len(TH.categorical_colours(4, cmap=next(iter(TH.CMAPS)), theme="dark")) == 4


# --------------------------------------------------------------------------- rendering mode
def test_additive_is_offered_but_is_not_what_anything_else_maps_to():
    """8,140 genes in dense UMAP clusters sum to white under additive blending and every colour mode
    renders as one featureless blob -- while every array-level test passes."""
    assert TH.gl_options("additive") == "additive"
    assert TH.gl_options("occlude") == "translucent"
    assert TH.gl_options("anything else") == "translucent"


def test_both_point_modes_are_declared():
    assert set(TH.POINT_MODES) == {"occlude", "additive"}


# --------------------------------------------------------------------------- literature
def _mentions():
    return pd.DataFrame({
        "gene_id": ["g1", "g1", "g2", "g2", "g3"],
        "pmid": ["1", "2", "1", "3", "4"],
        "doc_id": ["d1", "d2", "d1", "d3", "d4"],
        "section": ["title", "body", "abstract", "body", "caption"],
        "source": ["abstract", "fulltext", "abstract", "fulltext", "fulltext"],
    })


def test_the_strongest_tier_a_paper_affords_a_gene_is_the_one_counted():
    """A paper that names a gene in its title and again in a figure caption has afforded it focal
    attention, not two mentions of different weights."""
    out = LIT.attention_depth(_mentions())
    assert out.loc["g1", "n_papers_focal"] == 1
    assert out.loc["g1", "attention_depth"] == "focal"
    assert out.loc["g2", "n_papers_substantive"] == 1
    assert out.loc["g3", "attention_depth"] == "incidental"


def test_a_gene_named_twice_in_one_paper_counts_once():
    m = pd.DataFrame({"gene_id": ["g1", "g1"], "pmid": ["1", "1"], "doc_id": ["d1", "d1"],
                      "section": ["title", "body"], "source": ["fulltext", "fulltext"]})
    out = LIT.attention_depth(m)
    assert out.loc["g1", "n_papers_focal"] == 1
    assert out.loc["g1", "n_papers_incidental"] == 0


def test_a_document_with_no_pmid_falls_back_to_its_document_id():
    """Preprints and some deposits have no PMID, and keying on an empty string would merge them all
    into one paper."""
    m = pd.DataFrame({"gene_id": ["g1", "g2"], "pmid": ["", ""], "doc_id": ["d1", "d2"],
                      "section": ["title", "title"], "source": ["fulltext", "fulltext"]})
    out = LIT.attention_depth(m)
    assert out.loc["g1", "n_papers_focal"] == 1 and out.loc["g2", "n_papers_focal"] == 1


def test_an_unrecognised_section_is_treated_as_incidental_rather_than_dropped():
    m = pd.DataFrame({"gene_id": ["g1"], "pmid": ["1"], "doc_id": ["d1"],
                      "section": ["some_new_section"], "source": ["fulltext"]})
    assert LIT.attention_depth(m).loc["g1", "attention_depth"] == "incidental"


def test_publication_counts_can_be_restricted_to_one_source():
    """Abstracts and full texts are different populations -- all of PubMed against whichever papers a
    publisher deposited open access -- and support different claims, so they never merge silently."""
    m = _mentions()
    both = LIT.publication_counts(m)
    only_abstract = LIT.publication_counts(m, source="abstract")
    assert both.loc["g1"] == 2
    assert only_abstract.loc["g1"] == 1
    assert "g3" not in only_abstract.index


def test_attention_over_nothing_returns_the_right_empty_shape():
    """An empty frame with the wrong columns fails later, at the join, far from the cause."""
    out = LIT.attention_depth(_mentions().iloc[0:0])
    assert out.empty
    for c in ("n_papers_focal", "n_papers_substantive"):
        assert c in out.columns


def test_a_level_with_no_datasets_is_skipped_from_the_table():
    """The table is grouped by the on-disk taxonomy, and an empty heading would advertise a level of
    the tree that holds nothing."""
    from starplast import datasets as D
    real = D.REGISTRY
    try:
        D.REGISTRY = [d for d in real if d.level == "reference"]
        table = D.readme_table()
        assert D.LEVEL_TITLE["reference"] in table
        assert D.LEVEL_TITLE["DNA"] not in table
    finally:
        D.REGISTRY = real


def test_a_pair_sharing_one_unit_is_not_a_relation():
    """Two genes named once in the same paragraph is a coincidence; the threshold is what stops the
    graph filling with them."""
    from collections import Counter
    co = Counter({("g1", "g2"): 1, ("g1", "g3"): 5})
    meta = {"n_units": {"abstract": 100},
            "unit_hits": {"abstract": {"g1": 10, "g2": 5, "g3": 5}}}
    out = LIT.comention_edges(co, meta, "abstract", min_count=2)
    pairs = {(a, b) for a, b, _, _ in out}
    assert ("g1", "g3") in pairs
    assert ("g1", "g2") not in pairs


def test_an_unknown_name_falls_through_to_the_default(monkeypatch):
    """Resolution happens inside a paint call, so a name pyqtgraph does not ship must fall back rather
    than raise FileNotFoundError from inside rendering."""
    import pyqtgraph as pg
    calls = []
    real = pg.colormap.get

    def only_the_default(name, *a, **k):
        calls.append(name)
        if name == TH.DEFAULT_CMAP["sequential"]:
            return real(name, *a, **k)
        raise FileNotFoundError(name)

    monkeypatch.setattr(pg.colormap, "get", only_the_default)
    assert TH.resolve_cmap("no_such_map") is not None
    assert calls == ["no_such_map", TH.DEFAULT_CMAP["sequential"]]


def test_resolution_ends_rather_than_looping_when_nothing_can_be_built(monkeypatch):
    """pyqtgraph returns None rather than raising for some unknown names, so the loop's guard is
    `is not None` and the function still has to terminate in a value. It cannot invent a colour map,
    but it must not raise from inside a paint call either."""
    import pyqtgraph as pg
    calls = []
    monkeypatch.setattr(pg.colormap, "get", lambda name, *a, **k: calls.append(name))
    assert TH.resolve_cmap("no_such_map") is None
    assert calls[-1] == "viridis", "the last resort is asked for by name"


def test_a_column_of_labels_with_a_few_numbers_is_categorical():
    """Mostly non-numeric means labels, whatever the stray numbers look like."""
    assert TH.kind_for_column(pd.Series(["a", "b", "c", "d", 1, 2])) == "categorical"


def test_the_expectation_is_divided_by_units_not_by_the_sum_of_per_gene_counts():
    """An earlier version used a larger, differently-scaled denominator, which inflated every
    expectation and shrank every residual toward zero -- so the correction did almost nothing while
    appearing to be applied."""
    from collections import Counter
    co = Counter({("g1", "g2"): 10})
    meta = {"n_units": {"abstract": 100},
            "unit_hits": {"abstract": {"g1": 20, "g2": 20}}}
    (_, _, raw, residual), = LIT.comention_edges(co, meta, "abstract")
    assert raw == 10.0
    expected = 20 * 20 / 100                       # 4 units under independence
    import math
    assert residual == pytest.approx(math.log2((10 + 0.5) / (expected + 0.5)))


def test_an_empty_column_is_categorical_rather_than_a_ramp_over_nothing():
    """A filter that matched no genes leaves an empty column, and computing a range over it gives
    NaN bounds -- so the colour bar would be drawn from NaN to NaN."""
    assert TH.kind_for_column(pd.Series([], dtype=float)) == "categorical"
