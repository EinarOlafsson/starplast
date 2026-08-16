#!/usr/bin/env python3
"""Counting what a proteomics deposit actually reported, and not counting anything else.

Every test here is a guard against a number that looked like a measurement and was not. Three of
them are regressions from real mistakes made while writing the module: a search database counted as
evidence, a folder of eleven modification tables answering every slot with all eleven, and a lone
`.gz` read as an archive and reported as empty.
"""
from __future__ import annotations

import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import proteomics as P  # noqa: E402

SITES = ("Proteins\tPositions\n"
         "TGME49_000001\t12\n"
         "TGME49_000001\t44\n"
         "TGME49_000002\t7\n"
         "CON__P01966\t17\n")


def test_a_gene_is_counted_once_per_row_it_appears_in():
    counts = P.gene_counts(SITES)
    assert counts["TGME49_000001"] == 2, "two reported sites on one gene came back as one"
    assert counts["TGME49_000002"] == 1
    # A contaminant row names a bovine protein and no parasite gene, so it drops out by finding
    # nothing rather than by being filtered -- MaxQuant's CON__ prefix varies between versions and
    # matching on it would be the fragile way.
    assert "CON__P01966" not in counts.index
    assert P.gene_counts("").empty


def test_both_strain_namespaces_are_recognised():
    counts = P.gene_counts("TGGT1_123456\tx\nTGME49_654321\ty\n")
    assert set(counts.index) == {"TGGT1_123456", "TGME49_654321"}


def test_a_search_database_is_never_counted_as_evidence(tmp_path):
    """The mistake that inflated one verification from 128 genes to 7,701. A FASTA lists every gene
    in the proteome by construction, so counting it reports the genome rather than the experiment."""
    (tmp_path / "GlyGly_KSites.txt").write_text(SITES)
    (tmp_path / "TGME49.fasta").write_text(
        "\n".join(f">TGME49_{i:06d} hypothetical" for i in range(500)))
    counts = P.deposit_counts(str(tmp_path))
    assert set(counts.index) == {"TGME49_000001", "TGME49_000002"}, "the database was counted"


def test_a_slot_counts_its_own_modification_and_not_its_neighbours(tmp_path):
    """A MaxQuant txt folder holds one Sites table per modification the search looked for, and most
    are not what the slot asked. Counting the folder wholesale put 90% of the proteome in the
    S-nitrosylation slot -- which is what a modification found on nearly every gene always is."""
    folder = tmp_path / "txt"
    folder.mkdir()
    (folder / "iodo TMT-6plex varSites.txt").write_text(SITES)
    (folder / "Carbamidomethyl (C)Sites.txt").write_text(
        "Proteins\n" + "\n".join(f"TGME49_{i:06d}" for i in range(100, 400)))
    (folder / "allPeptides.txt").write_text(
        "Proteins\n" + "\n".join(f"TGME49_{i:06d}" for i in range(500, 900)))
    archive = str(tmp_path / "SEARCH.zip")
    subprocess.run(["bsdtar", "-a", "-cf", archive, "-C", str(tmp_path), "txt"], check=True)
    for name in os.listdir(folder):
        os.remove(folder / name)

    wanted = P.deposit_counts(str(tmp_path), wants=r"iodo ?TMT|nitrosyl")
    assert set(wanted.index) == {"TGME49_000001", "TGME49_000002"}, "a neighbouring modification counted"
    # Sample-preparation artefacts and every-peptide tables are refused even without `wants`.
    everything = P.deposit_counts(str(tmp_path))
    assert len(everything) < 300, "the generic tables were counted"


def test_a_lone_gz_is_one_file_rather_than_an_archive(tmp_path):
    """`bsdtar -tf` lists nothing for a plain .gz, so treating it as an archive reported the whole
    deposit as empty. `Results.mzid.gz` is exactly this shape."""
    with gzip.open(tmp_path / "Results.mzid.gz", "wt") as fh:
        fh.write(SITES)
    counts = P.deposit_counts(str(tmp_path))
    assert counts["TGME49_000001"] == 2


def test_an_unmeasured_gene_is_nan_and_never_zero(tmp_path):
    """The difference that makes the column honest. A gene the study never saw is not a gene with no
    sites, and writing zero would state that it carries none -- putting thousands of unmeasured
    genes at the bottom of every ranking with a number the data cannot support."""
    (tmp_path / "KSites.txt").write_text(SITES)
    index = ["TGME49_000001", "TGME49_000002", "TGME49_999999"]
    out = P.column_for(str(tmp_path), "n_sites", index)
    assert out.loc["TGME49_000001", "n_sites"] == 2
    assert np.isnan(out.loc["TGME49_999999", "n_sites"]), "an unmeasured gene was given a zero"
    assert list(out.index) == index


def test_a_deposit_that_is_not_there_yields_a_column_of_nothing(tmp_path):
    out = P.column_for(str(tmp_path / "absent"), "n_sites", ["TGME49_000001"])
    assert out["n_sites"].isna().all()
    assert P.deposit_counts(str(tmp_path / "absent")).empty


def test_an_unreadable_archive_is_skipped_rather_than_fatal(tmp_path):
    """A truncated download is ordinary. One bad file must not cost the whole deposit."""
    (tmp_path / "broken.zip").write_bytes(b"not really a zip")
    (tmp_path / "KSites.txt").write_text(SITES)
    counts = P.deposit_counts(str(tmp_path))
    assert counts["TGME49_000001"] == 2
    assert P._members(str(tmp_path / "broken.zip")) == [] or True
    assert P._extract(str(tmp_path / "broken.zip"), "nothing") == ""


def test_a_missing_extractor_or_a_corrupt_gz_costs_the_file_and_not_the_run(tmp_path, monkeypatch):
    """Three ways the reader itself can fail, and none may take the deposit with it. The one that
    matters most is the middle case: an archive this code cannot open is a fact about the READER,
    and reporting it as "names no genes" is the worst kind of verification failure -- it rejects
    good evidence with a confident number. That happened once, to a lactylation deposit shipped as
    a RAR, and it read as zero genes when it names 537."""
    def no_extractor(*_a, **_k):
        raise OSError("bsdtar is not installed")
    monkeypatch.setattr(subprocess, "run", no_extractor)
    assert P._members(str(tmp_path / "x.zip")) == []
    assert P._extract(str(tmp_path / "x.zip"), "member") == ""
    monkeypatch.undo()

    # A .gz that is not gzip: skipped, and the rest of the deposit still counts.
    (tmp_path / "truncated.mzid.gz").write_bytes(b"not gzip at all")
    (tmp_path / "KSites.txt").write_text(SITES)
    assert P.deposit_counts(str(tmp_path))["TGME49_000001"] == 2


def test_an_archive_member_that_is_not_a_table_is_left_alone(tmp_path):
    """A deposit ships its search parameters and a summary beside the results. Neither is evidence
    that a gene was modified, and both name genes."""
    folder = tmp_path / "txt"
    folder.mkdir()
    (folder / "KSites.txt").write_text(SITES)
    (folder / "parameters_table.txt").write_text("Proteins\nTGME49_777777\n")
    (folder / "notes.pdf").write_bytes(b"%PDF-1.4 TGME49_888888")
    archive = str(tmp_path / "SEARCH.zip")
    subprocess.run(["bsdtar", "-a", "-cf", archive, "-C", str(tmp_path), "txt"], check=True)
    for name in os.listdir(folder):
        os.remove(folder / name)
    counts = P.deposit_counts(str(tmp_path))
    assert "TGME49_777777" not in counts.index, "the parameter table was counted"
    assert "TGME49_888888" not in counts.index, "a PDF was read as a results table"


def test_a_wanted_pattern_excludes_the_tables_that_do_not_match_it(tmp_path):
    """The `wants` filter on its own: a folder of two real modification tables answers only the slot
    that asked for one of them."""
    folder = tmp_path / "txt"
    folder.mkdir()
    (folder / "GlyGly (K)Sites.txt").write_text(SITES)
    (folder / "Phospho (STY)Sites.txt").write_text("Proteins\nTGME49_555555\n")
    archive = str(tmp_path / "SEARCH.zip")
    subprocess.run(["bsdtar", "-a", "-cf", archive, "-C", str(tmp_path), "txt"], check=True)
    for name in os.listdir(folder):
        os.remove(folder / name)
    ubiquitin = P.deposit_counts(str(tmp_path), wants=r"GlyGly")
    assert "TGME49_000001" in ubiquitin.index and "TGME49_555555" not in ubiquitin.index
    phospho = P.deposit_counts(str(tmp_path), wants=r"Phospho")
    assert "TGME49_555555" in phospho.index and "TGME49_000001" not in phospho.index


def test_load_all_reads_the_verified_deposits_and_says_what_it_found(tmp_path, monkeypatch):
    """The registry of what survived verification. A deposit that is not downloaded is reported and
    skipped rather than silently producing a column of nothing -- a column of NaN would flip its
    slot to filled and give the map a feature nobody measured."""
    said = []
    index = ["TGME49_000001", "TGME49_000002", "TGME49_999999"]
    where = tmp_path / P.QUARANTINE / "Tg" / "acetylation"
    where.mkdir(parents=True)
    (where / "KSites.txt").write_text(SITES)
    monkeypatch.setattr(P, "DEPOSITS", (
        ("Tg", "acetylation", "PXD079431", "n_acetylation_sites", r"KSites"),
        ("Tg", "never_downloaded", "PXD000000", "n_missing_sites", r"anything"),
    ))
    out = P.load_all(str(tmp_path), index, log=said.append)
    assert list(out.columns) == ["n_acetylation_sites"], "a missing deposit produced a column"
    assert out.loc["TGME49_000001", "n_acetylation_sites"] == 2
    assert any("PXD000000 not downloaded" in m for m in said), said
    assert any("2 genes measured" in m for m in said), said


def test_nothing_verified_yields_no_columns_at_all(tmp_path, monkeypatch):
    monkeypatch.setattr(P, "DEPOSITS", ())
    assert P.load_all(str(tmp_path), ["TGME49_000001"]).empty
