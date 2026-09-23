"""Capture current application widgets for the getting-started slide deck.

Run under xvfb-run with QT_QPA_PLATFORM=xcb. Only temporary preferences and an
explicitly illustrative screen table are written; bundled gene evidence is read.
"""
from pathlib import Path
import os
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    """Render real controls, evidence and both map themes with isolated settings."""
    import numpy as np
    import pandas as pd
    from PyQt6 import QtCore, QtGui, QtWidgets, QtSvg
    from starplast.app import Window
    output = ROOT / 'docs/deck/assets'
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='starplast-deck-') as temporary:
        os.environ['STARPLAST_STATE'] = temporary
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
        QtCore.QSettings.setPath(QtCore.QSettings.Format.IniFormat,
                                QtCore.QSettings.Scope.UserScope, temporary)
        app = QtWidgets.QApplication([])
        window = Window(species='Toxoplasma gondii')
        window.resize(1600, 940)
        window.show()

        def save(widget, name, height=None):
            for _ in range(5): app.processEvents()
            area = QtCore.QRect(0, 0, widget.width(), height or widget.height())
            if not widget.grab(area).save(str(output / (name + '.png'))):
                raise RuntimeError('capture failed: ' + name)

        gene = 'TGME49_316400'
        window.on_pick(int(np.flatnonzero(window.nodes.gene_id.eq(gene))[0]))
        window.analysis_dock.hide()
        window.right_dock.show()
        save(window, 'workspace')
        save(window.detail, 'evidence')
        window._open_workflows(0)
        dialog = window.workflows_dialog
        dialog.gene.setText(gene)
        dialog._explore()
        save(dialog, 'explore')
        dialog.tabs.setCurrentIndex(1)
        save(dialog, 'predict', height=410)
        path = Path(temporary) / 'illustrative_screen.csv'
        pd.DataFrame({'gene_id': [gene, window.nodes.gene_id.iloc[0], 'unresolved_example'],
                      'example_score': [0., .4, np.nan]}).to_csv(path, index=False)
        dialog.load_screen(path)
        dialog._compare()
        dialog.tabs.setCurrentIndex(2)
        save(dialog, 'compare', height=390)
        dialog.hide()
        window.view.setParent(None)
        window.view.resize(1400, 850)
        window.view.show()
        center = np.median(window.xyz, axis=0)
        radius = float(np.quantile(np.linalg.norm(window.xyz - center, axis=1), .97))
        window.view.opts['center'] = QtGui.QVector3D(*map(float, center))
        window.view.setCameraPosition(distance=radius * 2.0, elevation=18, azimuth=45)
        window.show_ground = False
        window.depth_cue = False
        for theme in ('dark', 'light'):
            window.apply_theme(theme)
            window.redraw()
            save(window.view, 'map-' + theme)
        renderer = QtSvg.QSvgRenderer(str(ROOT / 'starplast/data/icons/starplast.svg'))
        picture = QtGui.QImage(930, 930, QtGui.QImage.Format.Format_ARGB32)
        picture.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(picture)
        renderer.render(painter)
        painter.end()
        picture.save(str(output / 'logo.png'))
        window.view.close()
        window.close()
        app.processEvents()
        print('Captured deck screenshots:', output)


if __name__ == '__main__':
    main()
