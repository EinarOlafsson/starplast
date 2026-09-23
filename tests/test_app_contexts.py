"""Exercise shader reuse across species windows with real shared OpenGL contexts."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.skipif(not os.environ.get('DISPLAY') or os.environ.get('QT_QPA_PLATFORM') == 'offscreen',
                    reason='needs a real GL display; run under xvfb-run with QT_QPA_PLATFORM=xcb')
def test_species_switch_renders_without_invalid_shader_programs(tmp_path):
    # Isolate QApplication creation from pytest-qt: the application must configure
    # sharing before Qt creates its global OpenGL context.
    code = '''
from starplast.app import Window
from PyQt6 import QtWidgets, QtGui
app = QtWidgets.QApplication([])
first = Window(species="Toxoplasma gondii")
first.show()
for _ in range(3): app.processEvents()
assert first.view.isValid()
assert not first.view.grabFramebuffer().isNull()
second = first.open_species("Plasmodium falciparum")
for _ in range(3): app.processEvents()
assert second.view.isValid()
assert QtGui.QOpenGLContext.areSharing(first.view.context(), second.view.context())
assert not second.view.grabFramebuffer().isNull()
third = second.open_species("Toxoplasma gondii")
for _ in range(3): app.processEvents()
assert third.view.isValid()
assert not third.view.grabFramebuffer().isNull()
third.close()
app.processEvents()
'''
    env = dict(os.environ, STARPLAST_STATE=str(tmp_path), XDG_CONFIG_HOME=str(tmp_path))
    result = subprocess.run([sys.executable, '-c', code], env=env,
                            cwd=Path(__file__).resolve().parents[1],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert 'Error while drawing item' not in result.stdout + result.stderr
    assert 'GLError' not in result.stderr
