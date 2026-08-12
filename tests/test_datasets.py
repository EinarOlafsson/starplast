#!/usr/bin/env python3
"""The provenance registry, and the three situations a dataset can be in.

The registry is the single source of truth for where every number came from: the download notebook
reads it, the methods section is generated from it, and the application shows it when asked. So the
properties that matter are that it is internally consistent, and that it tells the truth about what can
and cannot be obtained -- including when the answer is "this data does not exist in a repository".
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import datasets as D  # noqa: E402


# --------------------------------------------------------------------------- shape
def test_keys_are_unique():
    keys = [d.key for d in D.REGISTRY]
    assert len(keys) == len(set(keys))


def test_every_level_is_one_of_the_five_on_disk():
    """`level` is the directory layout, not a label: datasets/<level>/<type>/<PMID>/."""
    allowed = {"DNA", "transcription", "translation", "post_translation", "reference"}
    bad = {d.key: d.level for d in D.REGISTRY if d.level not in allowed}
    assert not bad, f"levels outside the on-disk taxonomy: {bad}"


def test_registry_filters_by_level_and_kind():
    assert D.registry(level="reference")
    assert all(d.level == "reference" for d in D.registry(level="reference"))
    k = D.REGISTRY[0].kind
    assert all(d.kind == k for d in D.registry(kind=k))


def test_get_raises_for_an_unknown_key():
    with pytest.raises((KeyError, StopIteration, ValueError)):
        D.get("no_such_dataset")


def test_provenance_maps_a_node_column_back_to_its_dataset():
    """The application's "where did this number come from" answer."""
    with_cols = [d for d in D.REGISTRY if d.columns]
    assert with_cols, "no dataset declares the columns it produces"
    d = with_cols[0]
    assert D.provenance(d.columns[0]).key == d.key


def test_provenance_of_an_unknown_column_is_none_not_a_guess():
    assert D.provenance("column_that_does_not_exist") is None


def test_unresolved_lists_exactly_the_entries_without_a_confirmed_citation():
    """None means 'confirm before citing'. These must not reach a manuscript unchecked, so the list is
    derived rather than maintained by hand."""
    assert {d.key for d in D.unresolved()} == {d.key for d in D.REGISTRY if d.citation is None}


def test_as_table_exposes_every_field():
    t = D.as_table()
    assert len(t) == len(D.REGISTRY)
    for c in ("key", "name", "level", "kind", "coverage", "pmid", "accession", "citation", "url"):
        assert c in t.columns


# --------------------------------------------------------------------------- what can be fetched
def test_fetchable_rejects_a_landing_page():
    """GEO's acc.cgi and ProteomeXchange's GetDataset return HTML describing the data. Saving that HTML
    as though it were the dataset stays invisible until something tries to parse it."""
    page = [k for k in (d.key for d in D.REGISTRY)
            if D.get(k).url and "proteomecentral" in D.get(k).url]
    assert page, "expected at least one landing-page URL in the registry"
    ok, how = D.fetchable(page[0])
    assert not ok and "landing page" in how


def test_fetchable_says_no_url_rather_than_guessing():
    no_url = [d.key for d in D.REGISTRY if not d.url]
    assert no_url
    ok, how = D.fetchable(no_url[0])
    assert not ok and "no download URL" in how


def test_a_geo_series_is_fetchable_through_the_ftp_listing_not_the_landing_page():
    gse = [d.key for d in D.REGISTRY if d.accession and d.accession.startswith("GSE")]
    assert gse
    ok, how = D.fetchable(gse[0])
    assert ok and how == "geo"


def test_direct_supplementary_urls_are_recognised():
    direct = [d.key for d in D.REGISTRY
              if d.url and ("MediaObjects" in d.url or "type=supplementary" in d.url)]
    assert direct
    assert D.fetchable(direct[0])[0]


def test_missing_reports_what_is_absent_here_without_raising():
    m = D.missing()
    assert isinstance(m, list)
    assert all(isinstance(k, str) for k in m)
    assert set(m) <= {d.key for d in D.REGISTRY}


def test_local_path_is_none_for_an_entry_with_no_path_recorded():
    no_path = [d.key for d in D.REGISTRY if not d.path]
    if no_path:
        assert D.local_path(no_path[0]) is None


def test_local_path_resolves_the_committed_cache_entries():
    """Registry paths starting starplast/data/ are part of the shipped cache, not raw inputs, and must
    resolve with no configuration at all."""
    cached = [d.key for d in D.REGISTRY if d.path and d.path.startswith("starplast/data/")]
    assert cached, "expected the identity table to be registered as part of the cache"
    for k in cached:
        assert D.local_path(k), f"{k} is part of the shipped cache and must always resolve"


def test_ensure_returns_none_and_says_why_when_it_cannot_fetch():
    """A path to something that is not the dataset would be worse than None."""
    msgs = []
    unfetchable = [d.key for d in D.REGISTRY
                   if not D.fetchable(d.key)[0] and not D.local_path(d.key)]
    if not unfetchable:
        pytest.skip("everything registered is either present or fetchable on this machine")
    assert D.ensure(unfetchable[0], log=msgs.append) is None
    assert msgs and "cannot fetch automatically" in msgs[0]


def test_ensure_returns_the_local_file_without_downloading_when_it_is_already_here(monkeypatch):
    present = [d.key for d in D.REGISTRY if D.local_path(d.key)]
    assert present
    monkeypatch.setattr("starplast.sources._get",
                        lambda *a, **k: pytest.fail("must not download what is already local"))
    assert D.ensure(present[0]) == D.local_path(present[0])


def test_ensure_refuses_an_html_error_page(monkeypatch, tmp_path):
    """A web server answering an error with a styled 200 page is the normal case, not the odd one, so
    size alone never proves a file arrived."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr("starplast.sources._get", lambda *a, **k: b"<!DOCTYPE html><html>nope")
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[0] and D.fetchable(d.key)[1] == "direct"][0]
    monkeypatch.setattr(D, "local_path", lambda k: None)
    msgs = []
    assert D.ensure(key, log=msgs.append) is None
    assert any("web page" in m for m in msgs)


def test_ensure_reports_a_transport_failure_rather_than_raising(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)

    def boom(*a, **k):
        raise OSError("network is down")

    monkeypatch.setattr("starplast.sources._get", boom)
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[1] == "direct"][0]
    msgs = []
    assert D.ensure(key, log=msgs.append) is None
    assert any("download failed" in m for m in msgs)


def test_ensure_writes_a_real_file_to_the_download_cache(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.sources._get", lambda *a, **k: b"gene_id,value\nTGME49_200010,1\n")
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[1] == "direct"][0]
    out = D.ensure(key, log=lambda *_: None)
    assert out and os.path.exists(out)
    assert open(out).read().startswith("gene_id")


def test_downloadable_lists_entries_with_a_url_in_notebook_order():
    """The reproduction notebook iterates this, so it is the list a reader actually runs."""
    d = D.downloadable()
    assert d == [x for x in D.REGISTRY if x.url]
    assert all(x.url for x in d)


def test_a_geo_series_is_fetched_through_the_ftp_supplementary_listing(monkeypatch, tmp_path):
    """Not the acc.cgi landing page: that returns HTML about the data rather than the data."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    seen = {}

    def fake_geo(acc, out_dir, log=print):
        seen["acc"], seen["dir"] = acc, out_dir
        p = os.path.join(out_dir, f"{acc}_matrix.tsv")
        open(p, "w").write("gene\tvalue\n")
        return [p]

    monkeypatch.setattr("starplast.sources.geo_supplementary", fake_geo)
    key = [x.key for x in D.REGISTRY if x.accession and x.accession.startswith("GSE")][0]
    out = D.ensure(key, log=lambda *_: None)
    assert out and out.endswith("_matrix.tsv")
    assert seen["acc"].startswith("GSE")


def test_a_geo_series_with_no_supplementary_files_returns_none(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.sources.geo_supplementary", lambda *a, **k: [])
    key = [x.key for x in D.REGISTRY if x.accession and x.accession.startswith("GSE")][0]
    assert D.ensure(key, log=lambda *_: None) is None


# --------------------------------------------------------------------------- checksums
def test_a_first_fetch_pins_what_arrived(monkeypatch, tmp_path):
    """A publisher reissuing a supplement under the same URL is the failure this exists for: the file
    changes, the build re-runs, every number moves a little, and nothing says why."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.sources._get", lambda *a, **k: b"gene_id,value\nTGME49_1,1\n")
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[1] == "direct"][0]
    out = D.ensure(key, log=lambda *_: None)
    assert out
    pinned = D.recorded_checksums()
    assert key in pinned
    assert pinned[key]["sha256"] == D.digest(out)
    assert pinned[key]["bytes"] == os.path.getsize(out)


def test_a_matching_file_verifies_silently(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    f = tmp_path / "thing.csv"
    f.write_text("a,b\n1,2\n")
    D.record_checksum("k", str(f))
    msgs = []
    assert D.verify_checksum("k", str(f), log=msgs.append)
    assert msgs == []


def test_a_changed_file_is_reported_rather_than_absorbed(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    f = tmp_path / "thing.csv"
    f.write_text("a,b\n1,2\n")
    D.record_checksum("k", str(f))
    f.write_text("a,b\n1,999\n")            # the publisher reissued it
    msgs = []
    assert not D.verify_checksum("k", str(f), log=msgs.append)
    joined = " ".join(msgs)
    assert "CHECKSUM MISMATCH" in joined
    assert "re-fetch" in joined, "the message has to say what to do, not only what is wrong"


def test_an_unpinned_file_is_not_a_failure_and_gets_pinned(monkeypatch, tmp_path):
    """Most of this tree arrived before checksums existed, so the honest answer for it is 'no claim' --
    and the digest is recorded so the NEXT fetch has something to check against."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    f = tmp_path / "thing.csv"
    f.write_text("a,b\n1,2\n")
    assert D.recorded_checksums() == {}
    assert D.verify_checksum("k", str(f), log=lambda *_: None)
    assert "k" in D.recorded_checksums()


def test_a_rejected_download_is_never_pinned(monkeypatch, tmp_path):
    """Otherwise the first bad download becomes the truth every good one is measured against."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.sources._get", lambda *a, **k: b"<!DOCTYPE html><html>nope")
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[1] == "direct"][0]
    assert D.ensure(key, log=lambda *_: None) is None
    assert key not in D.recorded_checksums()


def test_a_corrupt_checksum_file_does_not_stop_the_build(monkeypatch, tmp_path):
    """It is a cache of claims, not the data. Losing it costs the check, not the fetch."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    (tmp_path / D.CHECKSUMS).write_text("{not json")
    assert D.recorded_checksums() == {}


def test_the_digest_streams_rather_than_loading_the_file(tmp_path):
    """Some of these are gigabytes; reading one into memory to hash it is how a build gets killed."""
    f = tmp_path / "big.bin"
    f.write_bytes(b"x" * (3 << 20))
    import hashlib
    assert D.digest(str(f)) == hashlib.sha256(b"x" * (3 << 20)).hexdigest()


def test_an_already_downloaded_file_is_re_verified_rather_than_trusted(monkeypatch, tmp_path):
    """The check has to run on the cached copy, not only at the moment of download -- otherwise a file
    that changed on disk after being fetched is never noticed."""
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.sources._get", lambda *a, **k: b"gene_id,value\nTGME49_1,1\n")
    key = [d.key for d in D.REGISTRY if D.fetchable(d.key)[1] == "direct"][0]
    out = D.ensure(key, log=lambda *_: None)

    with open(out, "w") as fh:                      # something edited it after the fact
        fh.write("gene_id,value\nTGME49_1,999\n")
    msgs = []
    again = D.ensure(key, log=msgs.append)
    assert again == out, "the path is still returned; the caller decides what to do"
    assert any("CHECKSUM MISMATCH" in m for m in msgs)


def test_the_toxodb_identity_table_is_fetchable_despite_its_url_looking_like_a_page():
    """Its tabular report is a POST with a JSON body, so the URL alone looks like a landing page while
    the data is entirely fetchable -- fetch_names has done it all along. Reported as unfetchable, it
    understated what a clean machine can rebuild."""
    ok, how = D.fetchable("toxodb_identity")
    assert ok and how == "toxodb"


def test_the_toxodb_route_writes_the_identity_table(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)
    monkeypatch.setattr("starplast.fetch_names.fetch",
                        lambda org, attrs: "Gene ID\tGene Name or Symbol\nTGME49_1\tGRA16\n")
    out = D.ensure("toxodb_identity", log=lambda *_: None)
    assert out and os.path.exists(out)
    assert open(out).read().startswith("gene_id\tgene_name")
    assert "toxodb_identity" in D.recorded_checksums()


def test_kind_names_an_assay_never_a_provenance_word():
    """`kind` answers "what kind of measurement is this", not "how did we get it".

    stage_enriched was registered with kind "derived", which is the one entry in the registry whose
    type could not be read off its type field -- it is bulk RNA-seq underneath. Provenance has its
    own machine-readable home in derived_from, so mixing the two axes only loses the assay.
    """
    provenance_words = {"derived", "computed", "inferred", "transferred", "predicted", "internal"}
    for d in D.REGISTRY:
        assert d.kind.lower() not in provenance_words, (
            f"{d.key} has kind {d.kind!r}, which describes how it was produced rather than what it "
            f"is; put the assay in kind and mark the derivation with derived_from")


def test_a_derivation_is_marked_by_derived_from_not_by_its_kind():
    """The circularity guard reads derived_from, so a derivation that only says so in prose is one
    the guard cannot see."""
    derived = [d for d in D.REGISTRY if d.derived_from]
    assert derived, "no registry entry declares a derivation"
    for d in derived:
        assert d.kind.lower() not in ("derived",)
        # Whatever it was computed from has to be a real column of some other dataset, or the
        # declared provenance points at nothing.
        known = {c for other in D.REGISTRY for c in other.columns}
        for src in d.derived_from:
            assert src in known, f"{d.key} says it derives from {src!r}, which no dataset provides"


def test_every_dataset_says_how_to_obtain_it():
    """A dataset with neither a URL nor a declared derivation is unreproducible and unexplained.

    Nineteen entries had no URL, because the registry was written by describing files already on
    disk rather than by recording where they came from. Every one either has a download URL now, or
    is computed by this project and says so through derived_from.
    """
    for d in D.REGISTRY:
        assert d.url or d.derived_from, (
            f"{d.key} has no URL and declares no derivation, so nothing says how to get it; give it "
            f"a download URL, or set derived_from if this project computes it")


def test_a_url_template_says_what_to_substitute():
    """Some sources are per-accession endpoints rather than one file. A template is fine; a template
    nobody can tell is a template is not."""
    for d in D.REGISTRY:
        if d.url and "{" in d.url:
            assert d.note, f"{d.key} has a templated URL but no note explaining what to substitute"


def test_every_entry_states_what_kind_of_data_it_is():
    for d in D.REGISTRY:
        assert d.kind and d.kind.strip(), f"{d.key} has no kind"
        assert d.level and d.level.strip(), f"{d.key} has no level"


def test_a_failed_toxodb_request_reports_rather_than_raising(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_DATASETS, str(tmp_path))
    monkeypatch.setattr(D, "local_path", lambda k: None)

    def boom(org, attrs):
        raise OSError("toxodb is down")

    monkeypatch.setattr("starplast.fetch_names.fetch", boom)
    msgs = []
    assert D.ensure("toxodb_identity", log=msgs.append) is None
    assert any("ToxoDB request failed" in m for m in msgs)
