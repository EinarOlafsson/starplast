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


@pytest.mark.parametrize("version, prerelease", [
    ("0.43.0", "false"), ("0.43.0rc1", "true"), ("0.43.0.dev1", "true"),
])
def test_detection_marks_github_prereleases(tmp_path, monkeypatch, version, prerelease):
    """Stable releases, release candidates, and development versions reach GitHub correctly."""
    output = tmp_path / "outputs"
    monkeypatch.setattr(release, "check", lambda: version)
    monkeypatch.setattr("sys.argv", ["release.py", "detect", "--output", str(output)])
    release.main()
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert values == {"version": version, "publish": "true", "prerelease": prerelease}


def test_nightly_version_bump_publishes_after_merge_to_main(checkout):
    """Compare the merged version with main before the push, not the last nightly commit."""
    def git(*args):
        return subprocess.run(["git", *args], cwd=checkout, check=True, capture_output=True)

    git("branch", "-M", "main")
    previous_main = revision(checkout)
    git("switch", "-c", "nightly")
    release.bump("0.99.0", checkout)
    git("add", ".")
    git("-c", "user.name=Test", "-c", "user.email=test@example.invalid",
        "commit", "-qm", "Prepare release on nightly")
    git("switch", "main")
    git("merge", "--ff-only", "nightly")
    assert release.changed(previous_main, checkout)
    assert not release.changed(revision(checkout), checkout)
