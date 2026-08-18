#!/usr/bin/env python3
"""The pan-Apicomplexan archive: provenance first, and failures recorded rather than omitted."""
import json
import os

import pytest

from starplast import archive as A


class _Response:
    def __init__(self, payload): self.payload = payload
    def read(self): return self.payload
    def __enter__(self): return self
    def __exit__(self, *a): return False


def _opener(payload=b"gene,value\nTGME49_1,3\n"):
    return lambda url: _Response(payload)


def _source(**kw):
    base = dict(key="k", species="Toxoplasma gondii", study="a study",
                url="https://example.org/data/table.csv")
    base.update(kw)
    return A.Source(**base)


def test_a_retrieved_file_carries_a_checksum_and_a_shape(tmp_path):
    """"Downloaded" is not a claim anyone can check; a sha256 is."""
    record = A.fetch(_source(), str(tmp_path), opener=_opener())
    assert record.status == "retrieved"
    assert len(record.sha256) == 64 and record.bytes > 0
    assert "columns" in record.shape and record.retrieved


def test_an_existing_file_is_left_alone_and_recorded_as_present(tmp_path):
    """The archive is a store of what was retrieved and when. A re-run that silently replaced a 2019
    download with a 2026 one would destroy the provenance it exists to keep."""
    first = A.fetch(_source(), str(tmp_path), opener=_opener(b"original\n"))
    second = A.fetch(_source(), str(tmp_path), opener=_opener(b"REPLACEMENT\n"))
    assert second.status == "present"
    assert second.sha256 == first.sha256
    assert open(first.path).read() == "original\n"


def test_an_unreachable_source_is_recorded_with_its_reason_and_leaves_no_file(tmp_path):
    """An archive listing only what worked sends the next session to rediscover the same dead ends."""
    def broken(url):
        raise TimeoutError("no route")
    record = A.fetch(_source(), str(tmp_path), opener=broken)
    assert record.status == "unreachable" and "TimeoutError" in record.note
    assert not record.path
    assert not any(os.path.isfile(os.path.join(dp, f))
                   for dp, _dn, fn in os.walk(tmp_path) for f in fn)


def test_a_restricted_source_is_not_even_attempted(tmp_path):
    """Retrying a licence does not change it, and reporting a network failure where the truth is
    controlled access sends somebody to debug the wrong thing."""
    tried = []
    A.fetch(_source(status="restricted"), str(tmp_path),
            opener=lambda url: tried.append(url) or _Response(b""))
    assert not tried


def test_a_source_that_publishes_no_table_is_recorded_as_absent(tmp_path):
    assert A.fetch(_source(status="absent"), str(tmp_path), opener=_opener()).status == "absent"


def test_a_corrupt_download_is_reported_rather_than_counted_as_a_table(tmp_path):
    """A 404 page saved as .csv is a few kilobytes of valid HTML, and every archive that only checked
    file sizes has one."""
    record = A.fetch(_source(url="https://example.org/x.gz"), str(tmp_path),
                     opener=_opener(b"this is not gzip"))
    assert record.status == "retrieved" and record.shape.startswith("UNREADABLE")


def test_the_manifest_is_written_after_every_source_so_an_interrupted_run_keeps_what_it_got(tmp_path):
    manifest = str(tmp_path / "m" / "manifest.csv")
    sources = [_source(key=f"k{i}", study=f"study {i}") for i in range(3)]
    written = []

    def watching(url):
        written.append(os.path.exists(manifest))
        return _Response(b"a,b\n1,2\n")
    A.build(sources, str(tmp_path / "root"), manifest, opener=watching, log=lambda *a: None)
    assert written[1:] == [True] * (len(written) - 1), "the manifest waited until the end"
    assert os.path.exists(manifest) and os.path.exists(manifest.replace(".csv", ".json"))
    assert len(json.load(open(manifest.replace(".csv", ".json")))) == 3


def test_the_manifest_records_every_status_not_only_the_successes(tmp_path):
    sources = [_source(key="a", study="ok"), _source(key="b", study="gated", status="restricted"),
               _source(key="c", study="empty", status="absent")]
    records = A.build(sources, str(tmp_path / "root"), str(tmp_path / "m.csv"),
                      opener=_opener(), log=lambda *a: None)
    assert [r.status for r in records] == ["retrieved", "restricted", "absent"]
    assert set(A.STATUSES) >= {r.status for r in records}


def test_two_studies_do_not_collide_in_one_directory(tmp_path):
    a = A.fetch(_source(study="study one"), str(tmp_path), opener=_opener())
    b = A.fetch(_source(study="study/two: with punctuation"), str(tmp_path), opener=_opener())
    assert os.path.dirname(a.path) != os.path.dirname(b.path)
    assert "/" not in os.path.basename(os.path.dirname(b.path))


def test_a_query_string_does_not_become_part_of_the_filename(tmp_path):
    record = A.fetch(_source(key="bare", url="https://example.org/download?id=7"), str(tmp_path),
                     opener=_opener())
    assert os.path.basename(record.path) == "download", "the query string reached the filesystem"


def test_a_url_with_no_filename_at_all_still_lands_somewhere(tmp_path):
    record = A.fetch(_source(key="bare", url="https://example.org/files/"), str(tmp_path),
                     opener=_opener())
    assert record.status == "retrieved" and os.path.basename(record.path) == "bare.data"


@pytest.mark.parametrize("name,payload,expected", [
    ("t.zip", None, "zip"),
    ("t.gz", None, "gzip"),
    ("t.xlsx", None, "columns"),
])
def test_every_archive_kind_is_opened_rather_than_sized(tmp_path, name, payload, expected):
    """Verification, not a size check: the point is that the bytes are what they claim to be."""
    import gzip
    import io
    import zipfile
    if name.endswith(".zip"):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("inner.txt", "hello")
        payload = buf.getvalue()
    elif name.endswith(".gz"):
        payload = gzip.compress(b"gene,value\n1,2\n")
    else:
        import pandas as pd
        target = tmp_path / name
        pd.DataFrame({"a": [1], "b": [2]}).to_excel(target, index=False)
        payload = target.read_bytes()
        target.unlink()
    record = A.fetch(_source(url=f"https://example.org/{name}"), str(tmp_path),
                     opener=_opener(payload))
    assert expected in record.shape, f"{name} came back as {record.shape!r}"


def test_a_tar_is_recognised(tmp_path):
    import io
    import tarfile
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        info = tarfile.TarInfo("inner.txt")
        data = b"hello"
        info.size = len(data)
        tf.addfile(info, io.BytesIO(data))
    record = A.fetch(_source(url="https://example.org/t.tar"), str(tmp_path),
                     opener=_opener(buf.getvalue()))
    assert record.shape == "tar"


def test_an_unknown_extension_is_stored_without_a_false_claim(tmp_path):
    record = A.fetch(_source(url="https://example.org/t.bin"), str(tmp_path), opener=_opener())
    assert record.status == "retrieved" and record.shape == "not checked"


def test_a_download_that_dies_midway_leaves_no_partial_file(tmp_path):
    """A partial file is worse than none: it has a size, a checksum and no content."""
    class Dying:
        def read(self): raise ConnectionResetError("dropped")
        def __enter__(self): return self
        def __exit__(self, *a): return False
    record = A.fetch(_source(), str(tmp_path), opener=lambda url: Dying())
    assert record.status == "unreachable"
    assert not [f for dp, _dn, fn in os.walk(tmp_path) for f in fn], "a partial file was kept"


def test_the_catalogue_round_trips(tmp_path):
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps([{"key": "k", "species": "Pf", "study": "s",
                                 "url": "https://example.org/a.csv"}]))
    assert A.load_catalog(str(path))[0].species == "Pf"


# --------------------------------------------------------------------------- discovery
def _site_pages():
    release = '<a href="Pfalciparum3D7/">Pfalciparum3D7/</a><a href="Pberghei/">Pberghei/</a>'
    files = ('<a href="PlasmoDB-68_Pfalciparum3D7_AnnotatedProteins.fasta">p</a>'
             '<a href="PlasmoDB-68_Pfalciparum3D7_Genome.fasta">g</a>')

    def fetch(url):
        if url.endswith("Current_Release/"):
            return release
        if url.endswith("/fasta/data/"):
            return files
        raise FileNotFoundError(url)
    return fetch


def test_discovery_resolves_the_release_rather_than_remembering_it():
    """A release number typed into a URL is correct until the next release, and then it is a 404 that
    looks like a missing dataset."""
    found = A.discover_veupathdb(sites={"plasmodb.org": "Plasmodium"}, fetcher=_site_pages(),
                                 log=lambda *a: None)
    assert found, "nothing was discovered"
    assert all("Current_Release" in s.url for s in found)
    assert {s.species for s in found} == {"Pfalciparum3D7", "Pberghei"}


def test_discovery_keeps_only_the_gene_level_files():
    found = A.discover_veupathdb(sites={"plasmodb.org": "Plasmodium"}, fetcher=_site_pages(),
                                 log=lambda *a: None)
    assert all("AnnotatedProteins" in s.url for s in found), "a genome fasta was archived as gene-level"


def test_an_unreachable_site_is_recorded_rather_than_skipped_silently():
    def broken(url):
        raise ConnectionError("down")
    found = A.discover_veupathdb(sites={"plasmodb.org": "Plasmodium"}, fetcher=broken,
                                 log=lambda *a: None)
    assert len(found) == 1 and found[0].status == "unreachable"


def test_an_organism_with_no_fasta_directory_is_not_a_failure():
    """A great many organisms publish no fasta directory; that is a fact about the release rather
    than a broken walk, and recording it as unreachable would fill the manifest with noise."""
    def partial(url):
        if url.endswith("Current_Release/"):
            return '<a href="Ponly/">Ponly/</a>'
        raise FileNotFoundError(url)
    found = A.discover_veupathdb(sites={"plasmodb.org": "Plasmodium"}, fetcher=partial,
                                 log=lambda *a: None)
    assert found == []


def test_the_walk_can_be_bounded_for_a_quick_check():
    found = A.discover_veupathdb(sites={"plasmodb.org": "Plasmodium"}, fetcher=_site_pages(),
                                 limit_organisms=1, log=lambda *a: None)
    assert {s.species for s in found} == {"Pberghei"}


# --------------------------------------------------------------------------- the real urllib path
@pytest.fixture
def server(tmp_path):
    """A real HTTP server on localhost, so the production download path is exercised for real.

    Injecting an opener everywhere else keeps the suite off the network, but it also means the code
    that actually runs -- `urllib.request.urlopen` -- would never be executed by any test. Serving a
    directory locally covers it without depending on anything outside this machine.
    """
    import functools
    import http.server
    import threading
    root = tmp_path / "served"
    root.mkdir()
    (root / "table.csv").write_text("gene,value\nTGME49_1,3\n")
    (root / "index").mkdir()
    (root / "index" / "Current_Release").mkdir(parents=True, exist_ok=True)
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}", root
    httpd.shutdown()
    httpd.server_close()


def test_the_real_download_path_works_against_a_live_server(server, tmp_path):
    base, _root = server
    record = A.fetch(_source(url=f"{base}/table.csv"), str(tmp_path / "archive"))
    assert record.status == "retrieved" and record.bytes > 0
    assert "columns" in record.shape and len(record.sha256) == 64


def test_the_real_listing_path_reads_a_directory_index(server):
    base, _root = server
    assert "table.csv" in A._listing(f"{base}/")


def test_a_real_404_is_recorded_as_unreachable(server, tmp_path):
    base, _root = server
    record = A.fetch(_source(url=f"{base}/nothing-here.csv"), str(tmp_path / "archive"))
    assert record.status == "unreachable" and "HTTPError" in record.note


def test_a_fasta_is_counted_rather_than_left_unchecked(tmp_path):
    """The bulk of this archive is annotation FASTA. Reporting it as "not checked" would leave the
    largest part of the manifest making no claim about whether the bytes are what they say."""
    record = A.fetch(_source(url="https://example.org/p.fasta"), str(tmp_path),
                     opener=_opener(b">g1 one\nMKV\n>g2 two\nMKA\n"))
    assert record.shape == "2 sequences"


def test_a_fasta_with_no_header_is_reported_unreadable(tmp_path):
    record = A.fetch(_source(url="https://example.org/p.fasta"), str(tmp_path),
                     opener=_opener(b"<html>404 not found</html>"))
    assert record.shape.startswith("UNREADABLE")


def test_a_gff_reports_its_feature_rows_and_ignores_comments(tmp_path):
    record = A.fetch(_source(url="https://example.org/a.gff"), str(tmp_path),
                     opener=_opener(b"##gff-version 3\nchr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1\n"))
    assert record.shape == "1 feature rows"


def test_an_empty_gff_is_reported_unreadable(tmp_path):
    record = A.fetch(_source(url="https://example.org/a.gff"), str(tmp_path),
                     opener=_opener(b"##gff-version 3\n"))
    assert record.shape.startswith("UNREADABLE")
