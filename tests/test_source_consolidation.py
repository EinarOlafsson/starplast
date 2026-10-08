"""Recovery location observations must not become guessed biological associations."""
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from consolidate_source_recovery import direct_filename, member_manifest, project_input
from starplast import datasets as D, organisms as O


@pytest.mark.parametrize("url", [None, "https://toxodb.org/source.tsv", "https://sub.plasmodb.org/source.csv",
                                "https://veupathdb.org/source.xlsx", "http://publisher.example/table.xlsx",
                                "https://publisher.example/{accession}.xlsx", "https://publisher.example/GSE_RAW.tar",
                                "https://publisher.example/reads.fastq.gz", "https://publisher.example/foo%2Fbar.xlsx"])
def test_restricted_templates_raw_archives_and_guessed_names_are_not_requested(url):
    assert direct_filename(url) is None


def test_explicit_processed_resources_preserve_the_recorded_filename():
    assert direct_filename("https://ftp.ncbi.nlm.nih.gov/data/table.xls.gz") == "table.xls.gz"
    assert direct_filename("https://publisher.example/data%20table.xlsx") == "data table.xlsx"
    assert direct_filename("https://publisher.example/data.tsv?md5=123") == "data.tsv"


def test_container_members_are_bounded_and_raw_or_external_symlinks_are_not_bound(tmp_path):
    archive = tmp_path / "archive"
    directory = archive / "source"
    directory.mkdir(parents=True)
    (directory / "one.tsv").write_text("id\tvalue\na\t1\n")
    (directory / "two.csv").write_text("id,value\nb,2\n")
    (directory / "instrument.raw").write_bytes(b"raw")
    outside = tmp_path / "outside.tsv"
    outside.write_text("outside")
    (directory / "alias.tsv").symlink_to(outside)
    manifest = member_manifest(directory, archive)
    assert {Path(row["path"]).name for row in manifest["members"]} == {"one.tsv", "two.csv"}
    assert manifest["complete_within_processed_scope"]
    assert manifest["skipped"][0]["reason"] == "symlink_member_requires_receipt_review"
    partial = member_manifest(directory, archive, max_members=1)
    assert len(partial["members"]) == 1 and not partial["complete_within_processed_scope"]
    bytes_capped = member_manifest(directory, archive, max_bytes=1)
    assert not bytes_capped["members"] and not bytes_capped["complete_within_processed_scope"]
    with pytest.raises(ValueError, match="allowed archive"):
        member_manifest(tmp_path, archive)


def test_project_recovery_requires_the_exact_registered_path_and_refuses_escape(tmp_path):
    root = tmp_path / "datasets"
    root.mkdir()
    project = tmp_path / "toxo_stage_atlas"
    project.mkdir()
    (project / "table.tsv").write_text("id\tvalue\na\t1\n")
    source = D.Dataset("fixture", "fixture", "transcription", "test", "test", organism=O.TOXOPLASMA,
                       path="toxo_stage_atlas/table.tsv")
    assert project_input(source, root) == project / "table.tsv"
    missing = D.Dataset("missing", "missing", "transcription", "test", "test", organism=O.TOXOPLASMA,
                        path="toxo_stage_atlas/missing.tsv")
    assert project_input(missing, root) is None
    escape = D.Dataset("escape", "escape", "transcription", "test", "test", organism=O.TOXOPLASMA,
                       path="toxo_stage_atlas/../outside.tsv")
    with pytest.raises(ValueError, match="escapes"):
        project_input(escape, root)
