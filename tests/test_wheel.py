"""Reject auxiliary distributions and unexpected upload artifacts before publishing."""
from pathlib import Path
import sys
import zipfile

import pytest

from scripts.check_wheel import check_wheel, main


def test_rejects_an_auxiliary_distribution(tmp_path):
    wheel = tmp_path / 'starplast_core-0.42.0-py3-none-any.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('starplast_core-0.42.0.dist-info/METADATA', 'Name: starplast-core\nVersion: 0.42.0\n')
    with pytest.raises(ValueError, match='Unexpected distribution'):
        check_wheel(wheel)


def test_rejects_missing_bundled_data(tmp_path):
    wheel = tmp_path / 'starplast-0.42.0-py3-none-any.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        archive.writestr('starplast-0.42.0.dist-info/METADATA', 'Name: starplast\nVersion: 0.42.0\n')
    with pytest.raises(ValueError, match='Incomplete wheel'):
        check_wheel(wheel)


@pytest.mark.parametrize('extra', ['starplast_gpu-0.42.0.tar.gz', 'starplast-0.41.0.tar.gz', 'notes.txt'])
def test_rejects_extra_upload_artifacts(tmp_path, monkeypatch, extra):
    for name in ('starplast-0.42.0-py3-none-any.whl', 'starplast-0.42.0.tar.gz', extra):
        (tmp_path / name).touch()
    monkeypatch.setattr(sys, 'argv', ['check_wheel.py', str(tmp_path)])
    with pytest.raises(ValueError, match='only the matching'):
        main()
