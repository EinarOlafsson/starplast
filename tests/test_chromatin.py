#!/usr/bin/env python3
"""Promoter signal from coverage tracks, and the two ways it goes silently wrong.

The first is the strand: a reverse-strand gene's promoter is at its HIGH coordinate, and taking the
low one puts half of every promoter set at the wrong end of the gene while the column still looks
like data. The second is a replicate that cannot be read being dropped without saying so, which makes
the mean quieter than the note beside the column claims.
"""
from __future__ import annotations

import os
import struct
import sys
import tarfile

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import chromatin as CH  # noqa: E402

pyBigWig = pytest.importorskip("pyBigWig")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEADER = "Gene ID\tGenomic Sequence ID\tGenomic Location (Gene)\tGene Strand\tGene Type"


def _locations(tmp_path, rows, header=HEADER):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{g}\tchr1\t{loc}\tforward\tprotein coding gene" for g, loc in rows)
    (data / CH.LOCATION_TABLE).write_text(header + "\n" + body + "\n")
    return tmp_path


def _bigwig(path, intervals, chrom_length=100_000):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr1", chrom_length)])
    for start, end, value in intervals:
        bw.addEntries(["chr1"], [start], ends=[end], values=[float(value)])
    bw.close()
    return str(path)


# --------------------------------------------------------------------------- windows
def test_a_forward_gene_gets_the_window_at_its_low_coordinate(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)")])
    w = CH.gene_windows(str(tmp_path))
    assert (w.loc["TGME49_200010", "start"], w.loc["TGME49_200010", "end"]) == (9000, 11000)


def test_a_reverse_gene_gets_the_window_at_its_high_coordinate(tmp_path):
    """The failure that would look exactly like working data: the promoter at the wrong end."""
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(-)")])
    w = CH.gene_windows(str(tmp_path))
    assert (w.loc["TGME49_200010", "start"], w.loc["TGME49_200010", "end"]) == (19000, 21000)


def test_a_window_is_never_negative(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:100..500(+)")])
    assert CH.gene_windows(str(tmp_path)).loc["TGME49_200010", "start"] == 0


def test_accessions_go_through_the_identity_layer(tmp_path):
    _locations(tmp_path, [("TGGT1_100010", "chr1:10,000..20,000(+)")])
    w = CH.gene_windows(str(tmp_path),
                        resolve=lambda g: {"TGGT1_100010": "TGME49_200010"}.get(g))
    assert list(w.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _locations(tmp_path, [("TGME49_999999", "chr1:10,000..20,000(+)")])
    assert list(CH.gene_windows(str(tmp_path), resolve=lambda g: None).index) == ["TGME49_999999"]


def test_an_unparseable_location_is_dropped(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "unknown"), ("TGME49_200020", "chr1:1..2(+)")])
    assert list(CH.gene_windows(str(tmp_path)).index) == ["TGME49_200020"]


def test_a_table_of_only_unparseable_locations_yields_nothing(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "unknown")])
    assert CH.gene_windows(str(tmp_path)).empty


def test_no_location_table_yields_nothing(tmp_path):
    assert CH.gene_windows(str(tmp_path)).empty


def test_a_table_without_a_location_column_yields_nothing(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:1..2(+)")], header="Gene ID\tA\tB\tC\tD")
    assert CH.gene_windows(str(tmp_path)).empty


def test_a_gene_listed_twice_is_kept_once(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)"),
                          ("TGME49_200010", "chr1:10,000..20,000(+)")])
    assert len(CH.gene_windows(str(tmp_path))) == 1


# --------------------------------------------------------------------------- track summaries
def test_signal_is_relative_to_the_tracks_own_mean(tmp_path):
    """Two tracks differ in depth by whatever the submitters loaded, so raw means are incomparable."""
    windows = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [1000]},
                           index=pd.Index(["g"], name="gene_id"))
    shallow = _bigwig(tmp_path / "a.bw", [(0, 1000, 2.0), (1000, 2000, 1.0)])
    deep = _bigwig(tmp_path / "b.bw", [(0, 1000, 20.0), (1000, 2000, 10.0)])
    assert CH.track_means(shallow, windows)["g"] == pytest.approx(
        CH.track_means(deep, windows)["g"])


def test_a_window_off_the_end_of_a_sequence_is_missing_not_zero(tmp_path):
    """Zero would say 'closed chromatin' about a place the track never covered."""
    windows = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 500_000],
                            "end": [1000, 501_000]},
                           index=pd.Index(["a", "b"], name="gene_id"))
    path = _bigwig(tmp_path / "a.bw", [(0, 2000, 3.0)])
    got = CH.track_means(path, windows)
    assert np.isfinite(got["a"]) and np.isnan(got["b"])


def test_a_window_on_an_unknown_sequence_is_missing(tmp_path):
    windows = pd.DataFrame({"chrom": ["chrZ"], "start": [0], "end": [100]},
                           index=pd.Index(["a"], name="gene_id"))
    path = _bigwig(tmp_path / "a.bw", [(0, 2000, 3.0)])
    assert np.isnan(CH.track_means(path, windows)["a"])


def test_a_region_the_track_does_not_cover_is_missing(tmp_path):
    windows = pd.DataFrame({"chrom": ["chr1"], "start": [50_000], "end": [51_000]},
                           index=pd.Index(["a"], name="gene_id"))
    path = _bigwig(tmp_path / "a.bw", [(0, 2000, 3.0)])
    assert np.isnan(CH.track_means(path, windows)["a"])


# --------------------------------------------------------------------------- readability
def test_a_real_bigwig_is_recognised(tmp_path):
    assert CH.readable(_bigwig(tmp_path / "a.bw", [(0, 100, 1.0)]))


def test_the_corrupt_replicate_shape_is_rejected(tmp_path):
    """GSM8524430_UT_2.bw begins with eight 0xFF bytes, in the tar and at GEO alike."""
    p = tmp_path / "bad.bw"
    p.write_bytes(b"\xff" * 64)
    assert not CH.readable(str(p))


def test_an_unreadable_replicate_is_skipped_out_loud(tmp_path):
    """Silently dropping it makes the mean quieter than the note beside the column claims."""
    good = _bigwig(tmp_path / "GSM1_UT_1.bw", [(0, 100_000, 1.0)])
    bad = tmp_path / "GSM2_UT_2.bw"
    bad.write_bytes(b"\xff" * 64)
    archive = tmp_path / "series.tar"
    with tarfile.open(archive, "w") as tar:
        tar.add(good, arcname="GSM1_UT_1.bw")
        tar.add(bad, arcname="GSM2_UT_2.bw")
    windows = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [1000]},
                           index=pd.Index(["a"], name="gene_id"))
    said = []
    parts = CH._from_tar(str(archive), r"_UT_\d+\.bw$", windows, log=said.append)
    assert len(parts) == 1
    assert any("not a bigWig" in m for m in said), said


def test_members_that_do_not_match_the_arm_are_left_alone(tmp_path):
    good = _bigwig(tmp_path / "GSM1_UT_1.bw", [(0, 100_000, 1.0)])
    other = _bigwig(tmp_path / "GSM2_KO_1.bw", [(0, 100_000, 5.0)])
    archive = tmp_path / "series.tar"
    with tarfile.open(archive, "w") as tar:
        tar.add(good, arcname="GSM1_UT_1.bw")
        tar.add(other, arcname="GSM2_KO_1.bw")
        # A directory whose name matches the arm pattern: extractfile returns None for it.
        info = tarfile.TarInfo("GSM9_UT_9.bw")
        info.type = tarfile.DIRTYPE
        tar.addfile(info)
    windows = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [1000]},
                           index=pd.Index(["a"], name="gene_id"))
    assert len(CH._from_tar(str(archive), r"_UT_\d+\.bw$", windows, log=lambda *_: None)) == 1


# --------------------------------------------------------------------------- the loader
def test_nothing_downloaded_yields_no_columns(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)")])
    assert CH.chromatin_signals(str(tmp_path), log=lambda *_: None).empty


def test_no_locations_means_no_signal_at_all(tmp_path):
    assert CH.chromatin_signals(str(tmp_path), log=lambda *_: None).empty


def test_both_tracks_become_columns(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)")])
    root = tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg"
    (root / "acetylation").mkdir(parents=True)
    (root / "chromatin_accessibility").mkdir(parents=True)
    _bigwig(root / "acetylation" / "GSE313048_ATACseq_GCN5b-KD_UT.bw", [(0, 100_000, 4.0)])
    one = _bigwig(tmp_path / "GSM1_UT_1.bw", [(0, 100_000, 2.0)])
    with tarfile.open(root / "chromatin_accessibility" / "GSE277553_RAW.tar", "w") as tar:
        tar.add(one, arcname="GSM1_UT_1.bw")
    out = CH.chromatin_signals(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["atac_promoter_ut", "cuttag_hdac3_promoter_ut"]


def test_a_corrupt_atac_track_produces_no_column(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)")])
    root = tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg" / "acetylation"
    root.mkdir(parents=True)
    (root / "GSE313048_ATACseq_GCN5b-KD_UT.bw").write_bytes(b"\xff" * 64)
    assert CH.chromatin_signals(str(tmp_path), log=lambda *_: None).empty


def test_a_tar_of_only_unreadable_replicates_produces_no_column(tmp_path):
    _locations(tmp_path, [("TGME49_200010", "chr1:10,000..20,000(+)")])
    root = (tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg"
            / "chromatin_accessibility")
    root.mkdir(parents=True)
    bad = tmp_path / "GSM1_UT_1.bw"
    bad.write_bytes(b"\xff" * 64)
    with tarfile.open(root / "GSE277553_RAW.tar", "w") as tar:
        tar.add(bad, arcname="GSM1_UT_1.bw")
    assert CH.chromatin_signals(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- against the real map
@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_accessible_promoters_are_the_expressed_ones():
    """The check that says the strand and the join are both right.

    Open chromatin at a promoter is the single most reliable predictor of that gene being
    transcribed. If the window were at the wrong end of the gene, or the sequence names did not line
    up, this relationship is the first thing to vanish -- and the column would otherwise look fine.
    """
    from scipy.stats import spearmanr
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet"))
    if "atac_promoter_ut" not in n.columns:
        pytest.skip("ATAC column not merged")
    k = n[["atac_promoter_ut", "expr_tachy"]].dropna()
    assert spearmanr(k["atac_promoter_ut"], k["expr_tachy"]).statistic > 0.35
    q = k["expr_tachy"].quantile([0.1, 0.9])
    quiet = k.loc[k.expr_tachy <= q[0.1], "atac_promoter_ut"].median()
    loud = k.loc[k.expr_tachy >= q[0.9], "atac_promoter_ut"].median()
    assert loud - quiet > 1.0, f"silent {quiet:.2f} against expressed {loud:.2f}"
