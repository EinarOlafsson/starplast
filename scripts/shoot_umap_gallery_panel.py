"""Offscreen screenshots of the maps dock: the large gallery, and the table following the map.

Run offscreen; the images go to `results/umap_gallery_2026_09_30/`.

    QT_QPA_PLATFORM=offscreen python scripts/shoot_umap_gallery_panel.py
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("STARPLAST_GPU", "0")

from PyQt6 import QtWidgets  # noqa: E402

OUT = os.path.join(ROOT, "results", "umap_gallery_2026_09_30")


def shoot(widget, name: str, app):
    """One PNG of a widget, after letting Qt lay it out and paint."""
    os.makedirs(OUT, exist_ok=True)
    for _ in range(4):
        app.processEvents()
    path = os.path.join(OUT, name)
    widget.grab().save(path)
    print(f"{path}  {os.path.getsize(path) / 1e3:.0f} kB")
    return path


def main():
    from starplast import organisms, umap_gallery as G
    from starplast.umap_gallery_panel import (MapGalleryPanel, LABEL_VIEW, MAP_VIEW, SORT_STRUCTURE,
                                             SORT_LABEL)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    code = organisms.TOXOPLASMA
    gallery = G.shipped()
    if not gallery.available():
        raise SystemExit("no gallery is built")
    p = MapGalleryPanel(code, gallery=gallery)
    p.resize(760, 1000)
    p.show()
    maps = [r["id"] for r in gallery.maps(code)]
    print(f"{len(maps)} maps for {code}")

    # 1. the gallery as it opens: grouped, every heading named and counted
    for i in range(p.tree.topLevelItemCount()):
        p.tree.topLevelItem(i).setExpanded(False)
    shoot(p, "01_groups_collapsed.png", app)

    # 2. a group opened, with thumbnails, recipes and structure scores on the rows
    for name in ("Single experiments", "Pairs of families"):
        for i in range(p.tree.topLevelItemCount()):
            it = p.tree.topLevelItem(i)
            it.setExpanded(it.text(0).startswith(name))
        shoot(p, f"02_{name.split()[0].lower()}_open.png", app)

    # 3. flat, best-structured first: how to find the single best map
    p.group_check.setChecked(False)
    p.sort_box.setCurrentText(SORT_STRUCTURE)
    shoot(p, "03_flat_by_structure.png", app)

    # 4. sorted by how well the chosen label maps
    p.sort_box.setCurrentText(SORT_LABEL)
    shoot(p, "04_flat_by_label_skill.png", app)

    # 5. filtered
    p.search.setText("fitness")
    shoot(p, "05_filtered_fitness.png", app)
    p.search.setText("")
    p.group_check.setChecked(True)

    # 6-9. THE BUG: the table has to change for every map chosen, not only the first
    a, b = maps[4], maps[7]
    for view, tag in ((LABEL_VIEW, "label"), (MAP_VIEW, "map")):
        p.view_box.setCurrentText(view)
        for which, map_id in (("a", a), ("b", b)):
            p._on_item(map_id)
            shoot(p, f"0{6 if tag == 'label' else 8}{which}_table_{tag}_view_{map_id}.png", app)

    # 10. the whole dock in the real window, to show it fits
    from starplast import app as A
    w = A.Window()
    w.resize(1700, 1000)
    w.maps_dock.raise_()
    w.maps_panel._on_item(maps[6])
    shoot(w, "10_window_with_the_dock.png", app)
    w.console.remove()
    w.close()
    p.deleteLater()


if __name__ == "__main__":
    main()
