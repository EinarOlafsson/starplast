"""Data packs must round-trip, reject tampering and never destroy an installed version."""
import io
import json
import stat
import zipfile

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, packs, paths
from starplast.spaces import build_pack, validate_graph, validate_nodes


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path / "state"))
    space = O.Space("Xx", "Example species", "1", O.HOST, r"EX\d+", ("EX",), "", "",
                    "x_nodes.parquet", graph="x_graph.npz", distribution="pack")
    monkeypatch.setitem(O.SPACES, space.code, space)
    data = tmp_path / "data"
    data.mkdir()
    nodes = pd.DataFrame({"gene_id": ["EX1", "EX2"], "value": [0., 2.]})
    nodes.to_parquet(data / space.nodes, index=False)
    np.savez(data / space.graph, xyz=np.zeros((2, 3)), gene_ids=nodes.gene_id.to_numpy(dtype=str),
             example__a=np.array([0]), example__b=np.array([1]), example__w=np.array([1.]))
    archive = tmp_path / "pack.zip"
    sha = build_pack(space, data, "1.0", {c: "O" for c in nodes}, archive)
    return space, data, nodes, archive, sha


def rewrite(archive, mutate):
    with zipfile.ZipFile(archive) as source:
        members = {i.filename: source.read(i) for i in source.infolist()}
    mutate(members)
    with zipfile.ZipFile(archive, "w") as target:
        for name, contents in members.items():
            target.writestr(name, contents)


def test_roundtrip_active_paths_and_idempotent_install(prepared):
    space, data, nodes, archive, sha = prepared
    assert packs.installed(space.code) is None
    assert space.code not in O.codes(available=True)
    with pytest.raises(FileNotFoundError, match="no installed"):
        O.nodes_path(space.code)
    installed = packs.install(archive, sha, space.code)
    assert installed == packs.install(archive, sha, space.code)
    assert packs.installed(space.code) == installed
    assert O.nodes_path(space.code) == str(installed / space.nodes)
    assert O.graph_path(space.code) == str(installed / space.graph)
    assert space.code in O.codes(available=True)
    pd.testing.assert_frame_equal(pd.read_parquet(O.nodes_path(space.code)), nodes)
    assert packs.digest(installed / space.nodes) == packs.digest(data / space.nodes)


def test_tampered_archive_and_wrong_species_leave_active_pack_intact(prepared):
    space, data, nodes, archive, sha = prepared
    active = packs.install(archive, sha)
    with pytest.raises(ValueError, match="organism mismatch"):
        packs.install(archive, sha, "Other")
    rewrite(archive, lambda m: m.__setitem__(space.nodes, b"tampered"))
    with pytest.raises(ValueError, match="SHA256"):
        packs.install(archive, sha)
    with pytest.raises(ValueError, match="size mismatch"):
        packs.install(archive)
    assert packs.installed(space.code) == active
    pd.testing.assert_frame_equal(pd.read_parquet(active / space.nodes), nodes)


def test_same_version_cannot_replace_different_data(prepared):
    space, data, nodes, archive, sha = prepared
    original = packs.install(archive, sha)
    nodes.assign(value=3.).to_parquet(data / space.nodes, index=False)
    build_pack(space, data, "1.0", {c: "O" for c in nodes}, archive)
    with pytest.raises(ValueError, match="choose a new version"):
        packs.install(archive)
    pd.testing.assert_frame_equal(pd.read_parquet(original / space.nodes), nodes)
    build_pack(space, data, "2.0", {c: "O" for c in nodes}, archive)
    newer = packs.install(archive)
    assert newer.name == "2.0" and packs.installed(space.code, "1.0") == original


@pytest.mark.parametrize("name", ["../outside.txt", "/absolute.txt", "C:drive.txt", "a\\b.txt", "code.py"])
def test_unsafe_member_names_are_refused(prepared, name):
    space, data, nodes, archive, sha = prepared

    def mutate(members):
        manifest = json.loads(members[packs.MANIFEST])
        manifest["files"][name] = manifest["files"].pop(space.nodes)
        manifest["nodes"] = name
        members[name] = members.pop(space.nodes)
        members[packs.MANIFEST] = json.dumps(manifest)

    rewrite(archive, mutate)
    with pytest.raises(ValueError):
        packs.install(archive)
    assert packs.installed(space.code) is None


def test_symlinks_duplicate_and_undeclared_members_are_refused(prepared):
    space, data, nodes, archive, sha = prepared
    clean = archive.read_bytes()
    with zipfile.ZipFile(archive, "a") as out:
        out.writestr("extra.txt", "unlisted")
    with pytest.raises(ValueError, match="match manifest"):
        packs.verify(archive)
    archive.write_bytes(clean)
    with zipfile.ZipFile(archive, "a") as out, pytest.warns(UserWarning):
        out.writestr(packs.MANIFEST, "{}")
    with pytest.raises(ValueError, match="duplicate"):
        packs.verify(archive)
    with zipfile.ZipFile(io.BytesIO(clean)) as src, zipfile.ZipFile(archive, "w") as out:
        for item in src.infolist():
            item.external_attr = (stat.S_IFLNK | 0o777) << 16
            out.writestr(item, src.read(item))
    with pytest.raises(ValueError, match="regular files"):
        packs.verify(archive)


def test_size_and_checksum_limits(prepared, monkeypatch):
    space, data, nodes, archive, sha = prepared
    with pytest.raises(ValueError, match="64 hexadecimal"):
        packs.verify(archive, "invalid")
    monkeypatch.setattr(packs, "MAX_UNPACKED_BYTES", 1)
    with pytest.raises(ValueError, match="unpacked size"):
        packs.verify(archive)
    monkeypatch.setattr(packs, "MAX_ARCHIVE_BYTES", 1)
    with pytest.raises(ValueError, match="archive size"):
        packs.install(archive)


def test_build_refuses_lost_values_columns_ids_and_unreviewed_licenses(prepared):
    space, data, nodes, archive, sha = prepared
    for broken in (nodes.iloc[:1], nodes.drop(columns="value"), nodes.assign(value=np.nan),
                   nodes.assign(gene_id=["EX1", "EX1"]), nodes.assign(gene_id=["WRONG", "EX2"])):
        with pytest.raises(ValueError):
            validate_nodes(space, broken, nodes)
    validate_nodes(space, nodes.iloc[::-1], nodes)
    for licenses in ({"gene_id": "O"}, {"gene_id": "O", "value": "X"},
                     {"gene_id": "O", "value": "V"}):
        with pytest.raises(ValueError, match="license"):
            build_pack(space, data, "2", licenses, archive)
    assert packs.digest(archive) == sha


def test_graph_order_and_indices_are_checked(prepared):
    space, data, nodes, archive, sha = prepared
    graph = data / space.graph
    np.savez(graph, xyz=np.zeros((2, 3)), gene_ids=nodes.gene_id.iloc[::-1].to_numpy(dtype=str))
    with pytest.raises(ValueError, match="gene order"):
        validate_graph(nodes, graph)
    np.savez(graph, xyz=np.zeros((2, 3)), gene_ids=nodes.gene_id.to_numpy(dtype=str),
             layer__a=np.array([3]), layer__b=np.array([0]), layer__w=np.array([1.]))
    with pytest.raises(ValueError, match="indices"):
        validate_graph(nodes, graph)


def test_failed_activation_preserves_previous_version(prepared, monkeypatch):
    space, data, nodes, archive, sha = prepared
    original = packs.install(archive)
    build_pack(space, data, "2.0", {c: "O" for c in nodes}, archive)

    def fail(*args):
        raise OSError("disk failure")

    monkeypatch.setattr(packs.os, "replace", fail)
    with pytest.raises(OSError, match="disk failure"):
        packs.install(archive)
    assert packs.installed(space.code) == original


def test_download_requires_checksum_and_streams_into_verified_install(prepared, monkeypatch):
    space, data, nodes, archive, sha = prepared

    class Response(io.BytesIO):
        def geturl(self):
            return "https://example.invalid/pack.zip"

    monkeypatch.setattr(packs, "urlopen", lambda *a, **k: Response(archive.read_bytes()))
    for url, checksum in (("http://example.invalid/pack.zip", sha),
                          ("https://example.invalid/pack.zip", "")):
        with pytest.raises(ValueError, match="HTTPS"):
            packs.download(url, checksum)
    assert packs.download("https://example.invalid/pack.zip", sha, space.code).name == "1.0"


def test_commands_build_against_the_installed_baseline(prepared, capsys):
    from scripts import build_space
    space, data, nodes, archive, sha = prepared
    licenses = data / "licenses.json"
    licenses.write_text(json.dumps({c: "O" for c in nodes}))
    args = ["--organism", space.code, "--data", str(data), "--version", "2.0",
            "--licenses", str(licenses), "--output", str(archive)]
    assert packs.main(["installed", space.code]) == 1
    assert build_space.main(args) == 0  # a first build has no installed baseline
    assert packs.main(["verify", str(archive)]) == 0
    assert packs.main(["install", str(archive), "--organism", space.code]) == 0
    assert packs.main(["installed", space.code]) == 0
    before = archive.read_bytes()
    nodes.drop(columns="value").to_parquet(data / space.nodes, index=False)
    with pytest.raises(ValueError, match="loses"):
        build_space.main(args)
    assert archive.read_bytes() == before
