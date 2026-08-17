#!/usr/bin/env python3
"""Membrane lipid composition, and the one claim the whole slot rests on.

These species are only worth having if they are the PARASITE's lipids. The argument that they are is
not a sentence in a paper but a number in the archive: across four host backgrounds the host cells
differ in over a thousand species and the vesicles the same parasite released in them differ in
almost none. That is checked here against the shipped archive, because if it ever stopped holding the
column would have to come out.

The rest is the arithmetic that makes "composition" mean composition: a vesicle carries far less
total material than a cell, so the difference has to be taken after each sample is centred on itself,
or every lipid looks depleted and the column measures sample size.
"""
from __future__ import annotations

import io
import os
import sys
import zipfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import lipids as L  # noqa: E402
from starplast import metabolites as M  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVE = os.path.join(ROOT, "datasets", "reference", "lipidomics", "41716462",
                       "PMC12913473_supplementary.zip")
SPECIES = ["PI 38:4", "SM 34:1;O2", "CL 72:8", "PC 34:1"]


def _workbooks(tmp_path, *, levels=True, paired=True, loading_only=False, drop_a_replicate=False):
    """An archive shaped like the real one: a vesicle matrix, and vesicle-vs-cell intersected."""
    path = tmp_path / "supp.zip"
    ev = pd.DataFrame({"Lipid_names": SPECIES})
    cell = pd.DataFrame({"Lipid_names": SPECIES})
    base = np.array([10.0, 12.0, 14.0, 16.0])
    for i, host in enumerate(L.HOSTS):
        for rep in (1, 2):
            if drop_a_replicate and host == "Myo" and rep == 2:
                continue
            ev[f"EV_{host}_{rep}"] = base + i
            # A pure loading difference: the cell carries eight log2 units more of everything, and
            # nothing about the proportions differs. Unless `loading_only` is off, in which case one
            # species is genuinely enriched in the vesicle.
            shift = np.zeros(4) if loading_only else np.array([4.0, 0.0, 0.0, 0.0])
            cell[f"Cell_{host}_{rep}"] = base + i + 8.0 - shift
    with zipfile.ZipFile(path, "w") as z:
        if levels:
            buf = io.BytesIO()
            with pd.ExcelWriter(buf) as w:
                ev.to_excel(w, sheet_name=L.LEVELS[1], index=False)
            z.writestr(L.LEVELS[0], buf.getvalue())
        if paired:
            buf = io.BytesIO()
            with pd.ExcelWriter(buf) as w:
                ev.to_excel(w, sheet_name=L.PAIRED[1], index=False)
                cell.to_excel(w, sheet_name=L.PAIRED[2], index=False)
            z.writestr(L.PAIRED[0], buf.getvalue())
    return str(path)


# --------------------------------------------------------------------------- abundance
def test_the_level_is_the_mean_over_the_vesicle_samples(tmp_path):
    d = L.build(_workbooks(tmp_path)).set_index("metabolite")
    # base + host offset 0..3, averaged over hosts, is base + 1.5.
    assert d.loc["PI 38:4", "lipid_ev_level_log2"] == pytest.approx(11.5)
    assert d.loc["PC 34:1", "lipid_ev_level_log2"] == pytest.approx(17.5)


def test_every_species_gets_a_level(tmp_path):
    d = L.build(_workbooks(tmp_path))
    assert len(d) == len(SPECIES)
    assert d["lipid_ev_level_log2"].notna().all()


# --------------------------------------------------------------------------- composition
def test_a_pure_loading_difference_is_not_read_as_composition(tmp_path):
    """The vesicle holds less of everything. That is not a fact about lipids."""
    d = L.build(_workbooks(tmp_path, loading_only=True)).set_index("metabolite")
    assert d["lipid_ev_vs_host_clr"].abs().max() == pytest.approx(0.0, abs=1e-9)


def test_a_species_the_parasite_concentrates_shows_as_enriched(tmp_path):
    d = L.build(_workbooks(tmp_path)).set_index("metabolite")
    assert d.loc["PI 38:4", "lipid_ev_vs_host_clr"] > 0
    # Enriching one species of four pushes the other three down by a quarter of the shift, because
    # centring makes the column a proportion and proportions sum.
    assert d.loc["PI 38:4", "lipid_ev_vs_host_clr"] == pytest.approx(3.0)
    for other in SPECIES[1:]:
        assert d.loc[other, "lipid_ev_vs_host_clr"] == pytest.approx(-1.0)


def test_the_composition_column_is_centred(tmp_path):
    d = L.build(_workbooks(tmp_path))
    assert d["lipid_ev_vs_host_clr"].mean() == pytest.approx(0.0, abs=1e-9)


def test_replicates_are_matched_by_name_and_not_by_position(tmp_path):
    """The real vesicle arm is missing Myo_2; a positional read would take the wrong column."""
    full = L.build(_workbooks(tmp_path)).set_index("metabolite")
    gappy = L.build(_workbooks(tmp_path, drop_a_replicate=True)).set_index("metabolite")
    assert gappy["lipid_ev_vs_host_clr"].equals(full["lipid_ev_vs_host_clr"])


# --------------------------------------------------------------------------- refusals
def test_a_missing_archive_yields_nothing(tmp_path):
    assert L.build(str(tmp_path / "absent.zip")).empty


def test_an_archive_without_the_level_sheet_yields_nothing(tmp_path):
    assert L.build(_workbooks(tmp_path, levels=False)).empty


def test_levels_survive_without_the_paired_sheets(tmp_path):
    """Abundance and composition are separate questions; losing one must not lose the other."""
    d = L.build(_workbooks(tmp_path, paired=False))
    assert d["lipid_ev_level_log2"].notna().all()
    assert d["lipid_ev_vs_host_clr"].isna().all()


def test_a_level_sheet_without_the_name_column_is_refused(tmp_path):
    p = tmp_path / "supp.zip"
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame({"wrong": [1, 2]}).to_excel(w, sheet_name=L.LEVELS[1], index=False)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(L.LEVELS[0], buf.getvalue())
    assert L.build(str(p)).empty


def test_paired_sheets_without_the_name_column_give_no_composition(tmp_path):
    p = tmp_path / "supp.zip"
    ev = pd.DataFrame({"Lipid_names": SPECIES, "EV_Fibro_1": [1.0, 2, 3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        ev.to_excel(w, sheet_name=L.LEVELS[1], index=False)
    other = io.BytesIO()
    with pd.ExcelWriter(other) as w:
        pd.DataFrame({"wrong": [1]}).to_excel(w, sheet_name=L.PAIRED[1], index=False)
        pd.DataFrame({"wrong": [1]}).to_excel(w, sheet_name=L.PAIRED[2], index=False)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(L.LEVELS[0], buf.getvalue())
        z.writestr(L.PAIRED[0], other.getvalue())
    d = L.build(str(p))
    assert len(d) and d["lipid_ev_vs_host_clr"].isna().all()


def test_no_host_background_in_common_gives_no_composition(tmp_path):
    """Nothing named for a host means nothing to difference, not a crash."""
    p = tmp_path / "supp.zip"
    ev = pd.DataFrame({"Lipid_names": SPECIES, "sample_1": [1.0, 2, 3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        ev.to_excel(w, sheet_name=L.LEVELS[1], index=False)
    other = io.BytesIO()
    with pd.ExcelWriter(other) as w:
        ev.to_excel(w, sheet_name=L.PAIRED[1], index=False)
        ev.to_excel(w, sheet_name=L.PAIRED[2], index=False)
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(L.LEVELS[0], buf.getvalue())
        z.writestr(L.PAIRED[0], other.getvalue())
    assert L.build(str(p))["lipid_ev_vs_host_clr"].isna().all()


# --------------------------------------------------------------------------- the real archive
@pytest.mark.skipif(not os.path.exists(ARCHIVE), reason="lipidomics archive not fetched")
def test_the_vesicle_composition_does_not_move_when_the_host_does():
    """The entire case for calling these parasite lipids. Checked, not cited."""
    with zipfile.ZipFile(ARCHIVE) as z:
        host_cells, vesicles = L.conservation(z)
    assert host_cells.min() > 1000, "host cell lipidomes were supposed to differ by background"
    assert vesicles.max() <= 10, (
        "the vesicle lipidome now varies with the host: these can no longer be read as parasite "
        "lipids and the slot must be emptied")
    assert vesicles.max() < host_cells.min() / 50


@pytest.mark.skipif(not os.path.exists(ARCHIVE), reason="lipidomics archive not fetched")
def test_the_real_archive_yields_the_species_the_registry_claims():
    d = L.build(ARCHIVE)
    assert len(d) == 194
    assert d["lipid_ev_level_log2"].notna().all()
    assert d["lipid_ev_vs_host_clr"].notna().sum() == 192


@pytest.mark.skipif(not os.path.exists(ARCHIVE), reason="lipidomics archive not fetched")
def test_lipid_classes_the_parasite_cannot_make_are_depleted():
    """A sanity check on sign that does not depend on the study's own conclusions.

    Toxoplasma scavenges cholesterol and does not build sphingomyelin the way its host does, and
    cardiolipin is mitochondrial inner membrane rather than vesicle. If the sign convention were
    inverted these would read as concentrated in the vesicle.
    """
    d = L.build(ARCHIVE)
    d["klass"] = d["metabolite"].str.split().str[0]
    median = d.groupby("klass")["lipid_ev_vs_host_clr"].median()
    for scavenged in ("CE", "SM", "HexCer", "CL"):
        assert median[scavenged] < 0, scavenged
    assert median["PI"] > 0, "phosphatidylinositol is the parasite's own and should be enriched"


# --------------------------------------------------------------------------- the shipped table
@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", M.TABLE)),
                    reason="metabolite table not built")
def test_lipid_species_and_polar_compounds_do_not_collide_in_the_shipped_table():
    """`metabolites` warns that a second study joined by NAME is lossy. This is that second study.

    It is safe only because the two naming systems are disjoint. If a future source makes them
    overlap, this fails and the table needs a real identifier rather than a normalised name.
    """
    d = M.load(ROOT)
    lipid = d[d["lipid_ev_level_log2"].notna()]
    polar = d[d["metabolite_level_log2fc_iron_depleted"].notna()
              | d["labelled_fraction_glucose"].notna()]
    assert len(lipid) and len(polar)
    assert not set(lipid["metabolite"]) & set(polar["metabolite"])
    assert d["metabolite"].map(M.norm).duplicated().sum() == 0


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", M.TABLE)),
                    reason="metabolite table not built")
def test_the_polar_study_kept_every_row_it_had_before_the_lipids_arrived():
    d = M.load(ROOT)
    assert d["metabolite_level_log2fc_iron_depleted"].notna().sum() == 675
    assert len(d) == 1296
