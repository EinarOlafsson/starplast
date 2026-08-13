#!/usr/bin/env python3
"""Fetching from GEO and PRIDE, and normalizing by quantification type.

The quantification tests are the important ones, and they pin a pair of failures the same GEO series
produced in both directions. `GSE108740_FPKM.xlsx` reaches 16,520 and is genuine FPKM: skip the log and
one gene dominates every distance in the matrix. The node-table columns derived from it, still named
`rna108740_*_FPKM`, top out at 9.13 because they were logged upstream: log them again and real variation
compresses into nothing. Neither error raises anything, and neither is visible in the output.

So the range decides, not the filename, and both directions are regression cases here. An earlier
version of `infer_quant` returned "log_intensity" for anything non-integer, which classified real FPKM
as already-logged -- the exact failure the module exists to prevent.

No test touches the network. Both fetchers are exercised against captured responses.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import sources as S  # noqa: E402


# --------------------------------------------------------------------------- quantification type
def test_real_fpkm_is_recognised_as_needing_a_log():
    """16,520 in a column named FPKM is genuine FPKM. Returning "log_intensity" here -- as an earlier
    version did for anything non-integer -- skips the log and leaves one gene dominating every
    distance."""
    v = pd.Series([0.3949, 212.8, 6955.0, 16520.0, 1.5])
    assert S.infer_quant(v, "Tissue_cysts_A [FPKM]") == "fpkm"


def test_already_logged_values_are_not_logged_a_second_time():
    """The node table's rna108740_*_FPKM columns top out at 9.13 and were logged upstream. The name
    still says FPKM, so only the range can tell."""
    v = pd.Series([0.0, 2.5, 6.7, 9.13])
    assert S.infer_quant(v, "rna108740_Tachyzoites_T2_FPKM") == "log_intensity"


def test_a_log_fold_change_is_recognised_by_having_negatives():
    assert S.infer_quant(pd.Series([-3.2, 0.1, 2.8]), "MEDIAN_L2FC_IN_VIVO") == "lfc"


def test_a_name_saying_lfc_is_believed():
    assert S.infer_quant(pd.Series([0.5, 1.5, 2.5]), "gene_log2FC") == "lfc"


def test_a_ratio_is_left_alone():
    assert S.infer_quant(pd.Series([0.5, 1.5, 47.8]), "iTRAQ 115/113 ratio") == "ratio"


def test_integer_counts_are_recognised():
    assert S.infer_quant(pd.Series([0, 5, 100, 4000]), "raw_counts") == "counts"


def test_tpm_and_ibaq_are_recognised_by_name_when_the_range_allows():
    assert S.infer_quant(pd.Series([0.5, 100.0, 9000.5]), "gene_TPM") == "tpm"
    assert S.infer_quant(pd.Series([0.5, 100.0, 9000.5]), "protein_iBAQ") == "ibaq"


def test_an_empty_column_is_unknown_rather_than_guessed():
    assert S.infer_quant(pd.Series([np.nan, np.nan]), "whatever") == "unknown"


def test_non_numeric_values_do_not_crash_the_inference():
    assert S.infer_quant(pd.Series(["a", "b", None]), "x") == "unknown"


# --------------------------------------------------------------------------- normalization
def test_counts_are_logged_and_centred_per_sample():
    d = pd.DataFrame({"s1": [1.0, 3.0, 7.0], "s2": [10.0, 30.0, 70.0]})
    out = S.normalize(d, "counts", log=lambda *_: None)
    assert out.median().abs().max() == pytest.approx(0.0, abs=1e-9)
    assert out.s1.iloc[0] < out.s1.iloc[2]


def test_already_logged_data_is_centred_but_not_logged_again():
    d = pd.DataFrame({"s1": [1.0, 2.0, 3.0]})
    out = S.normalize(d, "log_intensity", log=lambda *_: None)
    assert out.s1.tolist() == pytest.approx([-1.0, 0.0, 1.0])


def test_a_ratio_is_returned_untouched():
    """Centring a ratio moves the zero, which IS the reference condition."""
    d = pd.DataFrame({"s1": [0.5, 1.0, 2.0]})
    assert S.normalize(d, "ratio", log=lambda *_: None).s1.tolist() == pytest.approx([0.5, 1.0, 2.0])


def test_a_log_fold_change_is_returned_untouched():
    d = pd.DataFrame({"s1": [-2.0, 0.0, 2.0]})
    assert S.normalize(d, "lfc", log=lambda *_: None).s1.tolist() == pytest.approx([-2.0, 0.0, 2.0])


def test_centring_is_per_sample_not_per_gene():
    """Per-sample corrects for how much material was loaded in each run, which is the systematic
    difference between columns. Per-gene would erase the between-gene differences that are the signal."""
    d = pd.DataFrame({"s1": [1.0, 2.0, 3.0], "s2": [11.0, 12.0, 13.0]})
    out = S.normalize(d, "log_intensity", log=lambda *_: None)
    assert out.s1.tolist() == pytest.approx(out.s2.tolist())
    assert out.iloc[0].tolist() != pytest.approx([0.0, 0.0])


def test_negative_values_are_clipped_before_the_log_rather_than_producing_nan():
    d = pd.DataFrame({"s1": [-5.0, 0.0, 3.0]})
    assert np.isfinite(S.normalize(d, "counts", log=lambda *_: None).to_numpy()).all()


def test_an_unknown_quantification_is_reported_and_left_alone():
    """Silently applying the wrong transform is how this module's failures happen."""
    msgs = []
    d = pd.DataFrame({"s1": [1.0, 2.0]})
    out = S.normalize(d, "something_new", log=msgs.append)
    assert out.s1.tolist() == [1.0, 2.0]
    assert any("unknown quantification" in m for m in msgs)


def test_every_declared_quantification_type_has_a_transform():
    """A type in QUANT_TYPES with no TRANSFORM entry would silently take the untouched path."""
    assert set(S.QUANT_TYPES) - {"unknown"} <= set(S.TRANSFORM)


# --------------------------------------------------------------------------- ranks
def test_ranks_are_bounded_and_order_preserving():
    d = pd.DataFrame({"s1": [10.0, 20.0, 30.0, 40.0]})
    out = S.rank_normalize(d)
    assert out.s1.min() >= -0.5 and out.s1.max() <= 0.5
    assert out.s1.is_monotonic_increasing


def test_missing_values_stay_missing_rather_than_becoming_a_rank():
    d = pd.DataFrame({"s1": [1.0, np.nan, 3.0]})
    assert pd.isna(S.rank_normalize(d).s1.iloc[1])


def test_ranks_put_incomparable_units_on_one_axis():
    """The only defensible way to compare an FPKM column with an iBAQ one: keep the ordering, which
    both support, and discard the units, which are not comparable."""
    fpkm = pd.DataFrame({"a": [1.0, 100.0, 16520.0]})
    ibaq = pd.DataFrame({"b": [0.001, 0.5, 0.9]})
    assert S.rank_normalize(fpkm).a.tolist() == pytest.approx(S.rank_normalize(ibaq).b.tolist())


# --------------------------------------------------------------------------- harmonising
def test_harmonise_joins_datasets_on_gene_and_prefixes_their_columns():
    a = pd.DataFrame({"x": [1.0, 2.0, 3.0]}, index=["g1", "g2", "g3"])
    b = pd.DataFrame({"y": [3.0, 2.0, 1.0]}, index=["g1", "g2", "g3"])
    out = S.harmonise({"first": (a, "log_intensity"), "second": (b, "log_intensity")},
                      log=lambda *_: None)
    assert list(out.columns) == ["first__x", "second__y"]
    assert len(out) == 3


def test_harmonise_infers_a_missing_quantification_type():
    a = pd.DataFrame({"x": [0.0, 2.5, 9.1]}, index=["g1", "g2", "g3"])
    msgs = []
    S.harmonise({"first": (a, None)}, log=msgs.append)
    assert any("as log_intensity" in m for m in msgs)


def test_harmonising_nothing_returns_an_empty_frame():
    assert S.harmonise({}, log=lambda *_: None).empty


# --------------------------------------------------------------------------- GEO
_GEO_LISTING = """<html><body>
<a href="../">Parent</a>
<a href="/geo/">GEO</a>
<a href="GSE108740_FPKM.xlsx">GSE108740_FPKM.xlsx</a>
<a href="filelist.txt">filelist.txt</a>
</body></html>"""


def test_a_geo_listing_is_parsed_and_the_files_downloaded(tmp_path, monkeypatch):
    grabbed = []

    def fake_get(url, timeout=300):
        if url.endswith("/"):
            return _GEO_LISTING.encode()
        grabbed.append(url)
        return b"payload"

    monkeypatch.setattr(S, "_get", fake_get)
    got = S.geo_supplementary("GSE108740", str(tmp_path), log=lambda *_: None)
    assert len(got) == 2
    assert all(os.path.exists(p) for p in got)
    assert any("GSE108740_FPKM.xlsx" in u for u in grabbed)


def test_parent_and_absolute_links_are_not_mistaken_for_files(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "_get",
                        lambda url, timeout=300: _GEO_LISTING.encode() if url.endswith("/") else b"x")
    got = S.geo_supplementary("GSE108740", str(tmp_path), log=lambda *_: None)
    assert not any(p.endswith("GEO") or p.endswith("Parent") for p in got)


def test_an_already_downloaded_file_is_not_fetched_again(tmp_path, monkeypatch):
    (tmp_path / "GSE108740_FPKM.xlsx").write_bytes(b"already here")
    (tmp_path / "filelist.txt").write_bytes(b"already here")
    calls = []

    def fake_get(url, timeout=300):
        calls.append(url)
        return _GEO_LISTING.encode() if url.endswith("/") else b"x"

    monkeypatch.setattr(S, "_get", fake_get)
    got = S.geo_supplementary("GSE108740", str(tmp_path), log=lambda *_: None)
    assert len(got) == 2
    assert len(calls) == 1, "only the listing should have been requested"


def test_a_failed_listing_returns_nothing_and_says_so(tmp_path, monkeypatch):
    def boom(url, timeout=300):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(S, "_get", boom)
    msgs = []
    assert S.geo_supplementary("GSE108740", str(tmp_path), log=msgs.append) == []
    assert any("listing failed" in m for m in msgs)


def test_one_failed_file_does_not_lose_the_others(tmp_path, monkeypatch):
    def fake_get(url, timeout=300):
        if url.endswith("/"):
            return _GEO_LISTING.encode()
        if "filelist" in url:
            raise OSError("connection reset")
        return b"payload"

    monkeypatch.setattr(S, "_get", fake_get)
    got = S.geo_supplementary("GSE108740", str(tmp_path), log=lambda *_: None)
    assert len(got) == 1 and got[0].endswith("FPKM.xlsx")


# --------------------------------------------------------------------------- PRIDE
def _pride_payload(categories):
    return json.dumps([
        {"fileCategory": {"value": c}, "fileName": f"f{i}.raw", "fileSizeBytes": 1000,
         "publicFileLocations": [{"name": "FTP Protocol", "value": f"ftp://x/f{i}.raw"}]}
        for i, c in enumerate(categories)]).encode()


def test_pride_categories_are_reported(monkeypatch):
    monkeypatch.setattr(S, "_get", lambda url, timeout=300: _pride_payload(["RAW", "RAW", "RESULT"]))
    t = S.pride_files("PXD004083", log=lambda *_: None)
    assert len(t) == 3
    assert set(t.category) == {"RAW", "RESULT"}
    assert t.url.notna().all()


def test_a_deposit_with_no_result_files_says_so_plainly(monkeypatch):
    """Checked across seven accessions now: PRIDE deposits instrument output, so per-protein numbers
    have to come from the paper's supplementary tables. Saying nothing would leave a caller believing
    a download was possible."""
    monkeypatch.setattr(S, "_get", lambda url, timeout=300: _pride_payload(["RAW", "SEARCH"]))
    msgs = []
    S.pride_files("PXD019729", log=msgs.append)
    assert any("no RESULT files" in m for m in msgs)


def test_a_pride_response_wrapped_in_embedded_is_handled(monkeypatch):
    """The v2 and v3 APIs disagree about whether the list is at the top level."""
    payload = json.dumps({"_embedded": {"files": [
        {"fileCategory": {"value": "RESULT"}, "fileName": "a.mzid", "fileSizeBytes": 5,
         "publicFileLocations": []}]}}).encode()
    monkeypatch.setattr(S, "_get", lambda url, timeout=300: payload)
    t = S.pride_files("PXD000001", log=lambda *_: None)
    assert len(t) == 1 and t.loc[0, "category"] == "RESULT"


def test_a_file_with_no_category_is_marked_unknown_not_dropped(monkeypatch):
    monkeypatch.setattr(S, "_get", lambda url, timeout=300: json.dumps(
        [{"fileName": "x", "fileSizeBytes": None, "publicFileLocations": None}]).encode())
    t = S.pride_files("PXD000001", log=lambda *_: None)
    assert t.loc[0, "category"] == "?"
    assert t.loc[0, "size_mb"] == 0.0


def test_a_pride_failure_returns_an_empty_frame(monkeypatch):
    def boom(url, timeout=300):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(S, "_get", boom)
    msgs = []
    assert S.pride_files("PXD000001", log=msgs.append).empty
    assert msgs


def test_an_empty_deposit_reports_nothing_rather_than_raising(monkeypatch):
    monkeypatch.setattr(S, "_get", lambda url, timeout=300: b"[]")
    assert S.pride_files("PXD000001", log=lambda *_: None).empty


# --------------------------------------------------------------------------- the transport
def test_the_fetcher_identifies_itself(monkeypatch):
    """Both archives ask for a contactable user agent, and an anonymous scraper gets rate-limited."""
    seen = {}

    class Resp:
        def read(self):
            return b"ok"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        seen["ua"] = req.headers
        return Resp()

    monkeypatch.setattr(S.urllib.request, "urlopen", fake_urlopen)
    assert S._get("https://example.org/x") == b"ok"
    assert any("starplast" in str(v) for v in seen["ua"].values())


def test_a_failed_fetch_is_logged_loudly_and_re_raised(monkeypatch):
    """WARNING rather than DEBUG, and re-raised rather than swallowed. This project shipped a
    version-pinned URL that returned 404 for months while the interface said "no model available" --
    a download that fails quietly is indistinguishable from a source that has nothing to give.

    The record is read off a handler on the module's own logger rather than through caplog: this
    project's loggers do not propagate to the root, so caplog sees them only when nothing else has
    configured logging yet, which makes it pass alone and fail in the suite.
    """
    import logging

    seen = []

    class Collect(logging.Handler):
        def emit(self, record):
            seen.append(record)

    def boom(req, timeout=None):
        raise OSError("connection reset")

    monkeypatch.setattr(S.urllib.request, "urlopen", boom)
    log = logging.getLogger("starplast.sources")
    handler = Collect(level=logging.WARNING)
    log.addHandler(handler)
    try:
        with pytest.raises(OSError, match="connection reset"):
            S._get("https://example.org/gone")
    finally:
        log.removeHandler(handler)
    assert any(r.levelno >= logging.WARNING and "example.org/gone" in r.getMessage()
               for r in seen), "the URL that failed has to be in the record"
