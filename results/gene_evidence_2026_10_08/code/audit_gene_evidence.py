"""Replay gene lookup, complete stored evidence and existing scorecard navigation."""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    """Hash a local artifact without loading it into RAM."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def audit(output):
    """Check both installed gene spaces without fitting or admitting biological truth."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    os.environ['STARPLAST_STATE'] = str(output / 'state')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    os.environ.setdefault('PYQTGRAPH_QT_LIB', 'PyQt6')
    from PyQt6 import QtCore, QtWidgets
    import pandas as pd
    from starplast import app as A, organisms as O, paths, track_record

    started = time.monotonic()
    code = ['starplast/gene_evidence.py', 'starplast/gene_evidence_browser.py',
        'starplast/app.py', 'starplast/query.py', 'starplast/identity.py',
        'starplast/slots.py', 'starplast/datasets.py', 'starplast/inventory.py',
        'starplast/dataset_space.py', 'starplast/dataset_space_inputs.py',
        'starplast/scorecard_browser.py', 'starplast/scorecard_view.py',
        'starplast/track_record.py', 'scripts/audit_gene_evidence.py', 'scripts/notebook_runner.py']
    inputs = [ROOT / name for name in code] + [Path(O.nodes_path(c)) for c in O.codes()]
    inputs += [Path(paths.cache_file(name)) for name in ('toxodb_identity.tsv', 'plasmodb_identity.tsv',
        'toxodb_strain_gt1.tsv', 'toxodb_strain_veg.tsv', 'track_record.parquet')]
    receipts = {str(path): sha(path) for path in inputs if path.is_file()}
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    checks, populations = [], []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    for organism, alias, expected in ((O.TOXOPLASMA, 'GRA16', 'TGME49_208830'),
            (O.FALCIPARUM, 'AMA1', 'PF3D7_1133400')):
        window = A.Window(species=O.get(organism).species)
        try:
            resolver = window._gene_resolver()
            result = resolver.resolve(alias)
            check('recorded alias in ' + organism, result.status == 'resolved' and result.entity.identifier == expected)
            check('alias source retained', all(choice.source for choice in result.choices))
            window.search.setText(alias)
            window.do_search()
            check('alias selects correct organism gene', window.nodes.gene_id.iloc[window.sel] == expected)
            check('gene page exposes full evidence entry', 'All evidence by biological question' in window.detail.toPlainText())
            panel = window.open_gene_evidence(window.sel)
            check('condensed evidence initial view', len(panel.shown_rows) <= 12)
            panel.show_all.setChecked(True)
            check('every original stored column reachable', {row['column'] for row in panel.shown_rows} == set(window.nodes.columns) - {'gene_id'})
            original = window.nodes.iloc[window.sel].drop('gene_id')
            reconstructed = pd.Series({row['column']: row['value'] for row in panel.rows}, dtype=original.dtype).reindex(original.index)
            pd.testing.assert_series_equal(original, reconstructed, check_exact=True, check_names=False)
            check('every stored value exact', True)
            for record in panel.rows:
                check('question/context/unit/gaps accessible', all(key in record for key in ('question', 'contexts', 'quantity_unit', 'gaps')))
                for source in record['sources']:
                    check('source remains organism/gene/column qualified', source['organism'] == organism and source['unit'] == 'gene' and record['column'] in source['columns'])
            index = next(i for i, record in enumerate(panel.shown_rows) if record['column'] == 'compartment')
            panel.table.setCurrentCell(index, 0)
            check('known label and class routes enabled', panel.label_open.isEnabled() and panel.class_open.isEnabled())
            label = str(panel.selected_row['value'])
            panel.class_open.click()
            check('class route reaches recorded result', 'Nothing recorded' not in window.detail.toPlainText())
            check('class retains actual ledger result', track_record.class_html('compartment', label, organism) in window.detail.toHtml() or label in window.detail.toPlainText())
            panel = window.open_gene_evidence(window.sel)
            panel.show_all.setChecked(True)
            index = next(i for i, record in enumerate(panel.shown_rows) if record['column'] == 'compartment')
            panel.table.setCurrentCell(index, 0)
            panel.label_open.click()
            check('label route reaches recorded result', 'Nothing recorded' not in window.detail.toPlainText())
            panel = window.open_gene_evidence(window.sel)
            panel.show_all.setChecked(True)
            index = next(i for i, record in enumerate(panel.shown_rows) if record['sources'])
            panel.table.setCurrentCell(index, 0)
            source = panel.source_choice.currentData()
            check('original source card reachable', panel.card.view == panel.space.source_card(source))
            panel.source_open.click()
            address = ('source_id', 'organism', 'unit', 'question')
            check('source route opens exact address', all(window._slot_tree.datasets.selected_row[key] == source[key] for key in address))
            panel.resize(1200, 850)
            panel.show()
            application.processEvents()
            panel.grab().save(str(output / (organism + '_evidence.png')))
            window.search.setText('kinase')
            window.do_search()
            check('product ambiguity requires choice', window.sel is None and 'Choose a gene' in window.detail.toPlainText())
            populations.append({'organism': organism, 'gene': expected, 'columns': len(panel.rows)})
            window._slot_tree.close()
        finally:
            window.close()
            window.deleteLater()
            QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete)
            application.processEvents()
            del window
            gc.collect()
    for path, digest in receipts.items():
        check('frozen input unchanged: ' + path, sha(path) == digest)
    snapshot = output / 'code'
    snapshot.mkdir()
    for name in code:
        (snapshot / Path(name).name).write_bytes((ROOT / name).read_bytes())
    summary = {'checks': len(checks), 'populations': populations, 'input_sha256': receipts,
        'runtime_seconds': time.monotonic() - started,
        'process_peak_rss_bytes': next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:')),
        'limitations': ['Stored evidence and existing annotation-recovery scorecards only',
            'No new independent biological truth, strategy fits, host installation or source promotion']}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (output / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    return summary


def main():
    """Execute acceptance in a notebook and preserve diagnostics on failure."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a fresh immutable directory')
    notebook = ExecutedNotebook('Gene evidence acceptance')
    try:
        notebook.code('from scripts.audit_gene_evidence import audit', f'summary = audit({str(args.out)!r})', "summary['checks']")
    except Exception as exc:
        args.out.mkdir(parents=True, exist_ok=True)
        notebook.md(type(exc).__name__ + ': ' + str(exc))
        notebook.write(str(args.out / 'diagnostic.ipynb'))
        raise
    notebook.write(str(args.out / 'audit.ipynb'))
    files = sorted(path for path in args.out.rglob('*') if path.is_file())
    (args.out / 'SHA256SUMS').write_text(''.join(sha(path) + '  ' + str(path.relative_to(args.out)) + '\n' for path in files))
    print(json.dumps({'checks': notebook.ns['summary']['checks']}))


if __name__ == '__main__':
    main()
