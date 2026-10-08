"""Execute dataset-browser acceptance against installed local evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _scope(row):
    return tuple(row[key] for key in ('source_id', 'organism', 'unit', 'question'))


def audit(output):
    """Check filters, every source, original values, paging and application routes."""
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    os.environ.setdefault('PYQTGRAPH_QT_LIB', 'PyQt6')
    from PyQt6 import QtCore, QtWidgets
    import pandas as pd
    import numpy as np
    from starplast import organisms as O, paths
    from starplast.dataset_browser import DatasetBrowser
    from starplast.dataset_space_inputs import from_sources
    from starplast.slot_tree import shipped_sources
    from starplast.scorecard_view import export_scorecard

    output = Path(output)
    if output.exists():
        raise ValueError('Use a new immutable dataset-browser audit directory')
    output.mkdir(parents=True)
    started = time.monotonic()
    code_files = ['starplast/dataset_space.py', 'starplast/dataset_space_inputs.py',
        'starplast/dataset_browser.py', 'starplast/source_refusals.py',
        'starplast/inventory.py', 'starplast/datasets.py', 'starplast/slots.py',
        'starplast/provenance.py', 'starplast/evidence.py', 'starplast/scorecard.py',
        'starplast/organisms.py', 'starplast/scorecard_view.py',
        'starplast/scorecard_browser.py', 'starplast/app.py', 'starplast/slot_tree.py',
        'scripts/audit_dataset_browser.py', 'scripts/notebook_runner.py']
    inputs = [ROOT / p for p in code_files] + [Path(O.nodes_path(c)) for c in O.codes()] + [Path(O.graph_path(c)) for c in O.codes()]
    inputs.append(ROOT / 'starplast/data/slots.json')
    inputs += [Path(paths.cache_file(p)) for p in (*O.HOST_TABLES.values(), 'metabolites.parquet',
        *(O.get(c).host_bridges for c in O.codes() if O.get(c).host_bridges))]
    receipts = {str(p): _sha(p) for p in inputs if p.is_file()}
    sources = shipped_sources()
    raw_tables = {}
    for code, source in sources.items():
        if source.get('nodes') is not None:
            raw_tables[(code, 'gene')] = source['nodes']
        side = source.get('tables') or {}
        if side.get('metabolite') is not None:
            raw_tables[(code, 'metabolite')] = side['metabolite']
        for host, frame in (side.get('host_gene') or {}).items():
            raw_tables[(host, 'protein')] = frame
    space = from_sources(sources)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    browser.resize(1450, 950)
    checks, reports, exports = [], [], {}

    def check(name, condition):
        assert condition, name
        checks.append(name)

    def select(row):
        index = next(i for i, shown in enumerate(browser.shown_rows) if _scope(shown) == _scope(row))
        browser.source_table.setCurrentCell(index, 0)
        app.processEvents()

    try:
        browser.show()
        app.processEvents()
        base = space.rows
        check('full registered inventory population', len(base) == 180)
        check('all question-specific refusals retained', sum(r['status'] == 'rejected' for r in base) == 16)
        for combo in browser.filters.values():
            combo.setCurrentIndex(0)
        for axis, combo in browser.filters.items():
            for i in range(1, combo.count()):
                value = combo.itemData(i)
                combo.setCurrentIndex(i)
                key = {'family': 'evidence_families', 'context': 'contexts'}.get(axis, axis)
                expected = [r for r in base if value in r[key]] if axis in {'family', 'context'} else [r for r in base if r[key] == value]
                check(axis + ': ' + str(value), [_scope(r) for r in browser.shown_rows] == [_scope(r) for r in expected])
                check(axis + ': visible count', browser.source_table.rowCount() == len(expected))
            combo.setCurrentIndex(0)
        for row in base:
            select(row)
            card = browser.card.view
            check(str(_scope(row)) + ': evidence only', card.kind == 'evidence' and card.metrics == ())
            payload = json.loads(browser.card.export_json())
            check(str(_scope(row)) + ': exact source scope', payload['snapshot']['card']['evidence']['inventory'] == row)
            frame = browser.entity_frame
            if row['unit'] in {'gene', 'protein', 'metabolite'} and frame.attrs['total_rows'] is not None:
                raw = raw_tables[(row['organism'], row['unit'])]
                columns = [frame.attrs['entity_column'], *[c for c in row['columns'] if c in raw and c != frame.attrs['entity_column']]]
                pd.testing.assert_frame_equal(frame.reset_index(drop=True), raw.iloc[:200][columns].reset_index(drop=True), check_exact=True)
                check(str(_scope(row)) + ': original values exact', True)
                mask = pd.Series(False, index=raw.index)
                for column in columns[1:]:
                    values = raw[column]
                    present = values.notna()
                    if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
                        present &= np.isfinite(values.fillna(0).astype(float))
                    elif not pd.api.types.is_bool_dtype(values):
                        present &= values.map(lambda value: not isinstance(value, str) or bool(value.strip()))
                    mask |= present
                check(str(_scope(row)) + ': independently counted storage coverage', row['stored_any_rows'] == int(mask.sum()) and row['table_rows'] == len(raw))
            elif frame.attrs['status'] == 'unavailable':
                check(str(_scope(row)) + ': unavailable gap', bool(frame.attrs['gap']) and bool(browser.value_note.text()))
            check(str(_scope(row)) + ': bounded page', len(frame) <= 200)
            for detail in card.details:
                check(str(_scope(row)) + ': detail ' + detail.key, browser.card.navigate('scorecard:detail/' + detail.key))
            exports['|'.join(_scope(row))] = payload
            reports.append({'scope': _scope(row), 'status': row['status'], 'origin': row['origin'],
                'stored_any_rows': row['stored_any_rows'], 'table_rows': row['table_rows'],
                'pair_records': row['pair_records'], 'displayed_rows': len(frame), 'total_rows': frame.attrs['total_rows']})
        gene_row = next(r for r in base if r['organism'] == O.TOXOPLASMA and r['unit'] == 'gene' and r['status'] == 'installed' and r['table_rows'] > 400)
        select(gene_row)
        count = 0
        while True:
            frame = browser.entity_frame
            raw = raw_tables[(gene_row['organism'], 'gene')]
            pd.testing.assert_frame_equal(frame.reset_index(drop=True), raw.iloc[count:count + 200][frame.columns].reset_index(drop=True), check_exact=True)
            count += len(frame)
            if not browser.page_next.isEnabled():
                break
            browser.page_next.click()
        check('every original gene row reachable through pages', count == len(raw))
        browser.page_previous.click()
        check('previous page navigation', browser.page_start < count - len(frame))
        select(gene_row)
        browser.card.navigate('scorecard:expand')
        app.processEvents()
        check('source screenshot', browser.grab().save(str(output / 'dataset_source.png')))
        browser.filters['organism'].setCurrentIndex(browser.filters['organism'].findData(O.HUMAN))
        app.processEvents()
        check('host filter screenshot', browser.grab().save(str(output / 'host_sources.png')))
    finally:
        browser.close()

    # Real menu/catalog routing and imported values; isolate user settings/state.
    state = output / 'temporary_state'
    os.environ['STARPLAST_STATE'] = str(state)
    QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
    for fmt in (QtCore.QSettings.Format.NativeFormat, QtCore.QSettings.Format.IniFormat):
        QtCore.QSettings.setPath(fmt, QtCore.QSettings.Scope.UserScope, str(state))
    from starplast.app import Window
    window = Window(species=O.get(O.TOXOPLASMA).species)
    try:
        window.dataset_browser_act.trigger()
        panel = window._slot_tree.datasets
        check('menu opens existing catalog dataset tab', panel is not None and window._slot_tree.tabs.currentIndex() == 1)
        row = next(r for r in panel.shown_rows if r['unit'] == 'gene' and r['status'] == 'installed')
        index = next(i for i, r in enumerate(panel.shown_rows) if _scope(r) == _scope(row))
        panel.source_table.setCurrentCell(index, 0)
        entity = panel.entity_frame.iloc[0]['gene_id']
        panel._open_entity(panel.entity_table.model().index(0, 0))
        check('real application gene search route', window.search.text() == entity)
        column = 'imported_browser_acceptance'
        window.nodes = window.nodes.assign(**{column: False})
        record = {'columns': [column], 'genes': len(window.nodes), 'quantification': 'authored software check'}
        window.imports.append(record)
        current = dict(window._slot_tree.view._sources)
        current[O.TOXOPLASMA] = {**current[O.TOXOPLASMA], 'nodes': window.nodes}
        window._slot_tree.refresh_sources(current, imports=window.imports, imported_organism=O.TOXOPLASMA)
        imported = next(r for r in panel.shown_rows if r['origin'] == 'user_imported')
        index = next(i for i, r in enumerate(panel.shown_rows) if _scope(r) == _scope(imported))
        panel.source_table.setCurrentCell(index, 0)
        check('session import origin remains distinct', panel.selected_row['source_id'].startswith('session_import:'))
        check('False import cells retained', panel.entity_frame[column].eq(False).all())
    finally:
        window._slot_tree.close()
        window.close()
    for p, h in receipts.items():
        check('input unchanged: ' + p, _sha(p) == h)
    snapshot = output / 'code'
    snapshot.mkdir()
    for relative in code_files:
        (snapshot / Path(relative).name).write_bytes((ROOT / relative).read_bytes())
    summary = {'checks': len(checks), 'inventory_rows': len(reports), 'sources': reports,
        'runtime_seconds': time.monotonic() - started,
        'process_peak_rss_bytes': next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:')),
        'input_sha256': receipts, 'limitations': ['Coverage counts stored values, not biological accuracy',
            'Missing lineage, assay units and raw provenance remain unresolved',
            'Session import is an authored UI control; no publication or biological admission',
            'No source downloads, strategy fits or calibrated deployment']}
    for name, data in [('summary.json', summary), ('source_cards.json', exports), ('checks.json', checks)]:
        (output / name).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return summary


def main():
    """Retain executed notebook evidence or a diagnostic in a new directory."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a new immutable audit directory')
    nb = ExecutedNotebook('Dataset browser acceptance')
    nb.md('Replay the full installed inventory, scoped filters, source cards, original values and application navigation. Missing provenance stays unresolved; storage coverage is not biological accuracy.')
    try:
        nb.code('from scripts.audit_dataset_browser import audit', f'summary = audit({str(args.out)!r})', "{'checks': summary['checks'], 'inventory_rows': summary['inventory_rows']}")
    except Exception as exc:
        args.out.mkdir(parents=True, exist_ok=True)
        nb.md('Diagnostic: ' + type(exc).__name__ + ': ' + str(exc))
        nb.write(str(args.out / 'diagnostic.ipynb'))
        raise
    nb.write(str(args.out / 'audit.ipynb'))
    (args.out / 'SHA256SUMS').write_text(''.join(_sha(p) + '  ' + p.relative_to(args.out).as_posix() + '\n' for p in sorted(args.out.rglob('*')) if p.is_file()))
    print(json.dumps({'checks': nb.ns['summary']['checks'], 'inventory_rows': nb.ns['summary']['inventory_rows']}))


if __name__ == '__main__':
    main()
