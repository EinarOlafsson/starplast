"""Release guards prevent mismatched dependency pins and unintended publications."""
from pathlib import Path
import shutil
import subprocess

import pytest

from scripts import release


@pytest.fixture
def checkout(tmp_path):
    for relative in (*release.PROJECTS, Path("starplast/__init__.py")):
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(release.ROOT / relative, target)
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "."], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                    "commit", "-qm", "baseline"], cwd=tmp_path, check=True)
    return tmp_path


def revision(root):
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()


def test_unchanged_version_does_not_publish(checkout):
    assert not release.changed(revision(checkout), checkout)


def test_bump_updates_every_distribution_and_enables_publication(checkout):
    previous = revision(checkout)
    assert release.bump("0.99.0", checkout) == "0.99.0"
    assert release.changed(previous, checkout)
    assert release.check(checkout) == "0.99.0"


def test_mismatched_metapackage_cannot_publish(checkout):
    p = checkout / release.PROJECTS[1]
    p.write_text(p.read_text().replace(f'version = "{release.check(checkout)}"', 'version = "0.1.0"'))
    with pytest.raises(ValueError, match="expected"):
        release.check(checkout)


def test_stale_extra_pin_cannot_publish(checkout):
    p = checkout / release.PROJECTS[1]
    p.write_text(p.read_text().replace(f'starplast-core[gpu]=={release.check(checkout)}', 'starplast-core[gpu]==0.1.0'))
    with pytest.raises(ValueError, match="dependency pins"):
        release.check(checkout)


@pytest.mark.parametrize("new", ["0.1.0", "0.42.0", "0.43.0;echo bad", "v0.43.0"])
def test_invalid_bump_leaves_files_untouched(checkout, new):
    before = {p: (checkout / p).read_bytes() for p in release.PROJECTS}
    with pytest.raises(ValueError):
        release.bump(new, checkout)
    assert before == {p: (checkout / p).read_bytes() for p in release.PROJECTS}


def test_publishing_rejects_a_consistent_but_older_version(checkout):
    release.bump("0.99.0", checkout)
    subprocess.run(["git", "add", "."], cwd=checkout, check=True)
    subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                    "commit", "-qm", "newer"], cwd=checkout, check=True)
    previous = revision(checkout)
    subprocess.run(["git", "checkout", "HEAD~1", "--", "."], cwd=checkout, check=True)
    with pytest.raises(ValueError, match="downgrade"):
        release.changed(previous, checkout)


def test_only_object_ids_can_be_used_as_previous_revision(checkout):
    with pytest.raises(ValueError, match="object ID"):
        release.changed("HEAD:unrelated-file", checkout)


def test_first_push_can_publish(checkout):
    assert release.changed("0" * 40, checkout)
