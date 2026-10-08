"""Verify frozen task cards in two real desktop hosts without fitting models."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/scorecard_task_views_2026_10_08'


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(output):
    """Check native exports and navigation, preserving immutable audit evidence."""
    started = time.monotonic()
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6 import QtCore, QtGui, QtWidgets
    from starplast import host_foundation_scorecard as HF, scorecard_view as V
    from starplast.scorecard_browser import ScorecardBrowser

    output = Path(output)
    if output.exists():
        raise ValueError('Use a new immutable desktop audit directory')
    output.mkdir(parents=True)
    receipts = {}
    for line in (SOURCE / 'SHA256SUMS').read_text().splitlines():
        digest, relative = line.split(maxsplit=1)
        assert _sha(SOURCE / relative) == digest, relative
        receipts[relative] = digest
    exports = json.loads((SOURCE / 'exports.json').read_text())
    (output / 'source_exports.json').write_bytes((SOURCE / 'exports.json').read_bytes())
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    checks, reports, tasks = [], [], set()

    def check(name, condition):
        assert condition, name
        checks.append(name)

    views = {}
    for name, original in exports.items():
        snap = original['snapshot']
        view = V.build_scorecard_view(snap['card'], task=original['task'],
            kind=original['kind'], title=original['title'], status=original['status'],
            details=snap['details'], links=tuple(V.ScorecardLink(**link) for link in original['links']))
        current = json.loads(V.export_scorecard(view))
        check(name + ': original snapshot exact', current['snapshot'] == snap)
        # JSON normalizes dictionary order; supplemental metric order is not identity.
        check(name + ': original metrics exact',
            {m['key']: m for m in current['metrics']} == {m['key']: m for m in original['metrics']})
        check(name + ': original provenance exact', current['source'] == original['source'])
        check(name + ': dedicated calibration section', any(d.key == 'calibration' for d in view.details))
        views[name] = view
    host, warning = HF.shipped()
    check('human foundation source verified', host is not None and not warning)
    views['human_foundation'] = host

    for name, view in views.items():
        if view.task is not None:
            tasks.add(view.task)
        browsers = [ScorecardBrowser(), ScorecardBrowser()]
        try:
            for browser in browsers:
                browser.set_scorecard(view)
                browser.resize(1100, 800)
                browser.show()
            app.processEvents()
            check(name + ': compact text identical', browsers[0].toPlainText() == browsers[1].toPlainText())
            check(name + ': desktop exports identical', browsers[0].export_json() == browsers[1].export_json() == V.export_scorecard(view))
            expected = QtGui.QTextDocument()
            for metric in view.metrics:
                route = 'scorecard:metric/' + metric.key
                expected.setHtml(V.render_scorecard_detail(view, metric.key) + '<p><a href="scorecard:expand">Full test record</a></p>')
                for browser in browsers:
                    check(name + ': metric ' + metric.key, browser.navigate(route))
                    check(name + ': metric text ' + metric.key, browser.toPlainText() == expected.toPlainText())
            for detail in view.details:
                expected.setHtml(V.render_scorecard_detail(view, detail.key) + '<p><a href="scorecard:expand">Full test record</a></p>')
                for browser in browsers:
                    check(name + ': detail ' + detail.key, browser.navigate('scorecard:detail/' + detail.key))
                    check(name + ': detail text ' + detail.key, browser.toPlainText() == expected.toPlainText())
            for browser in browsers:
                check(name + ': expansion', browser.navigate('scorecard:expand'))
                check(name + ': export unchanged after navigation', browser.export_json() == V.export_scorecard(view))
                check(name + ': unknown outcome rejected', not browser.navigate('scorecard:outcome/not_supplied'))
                check(name + ': unrelated source rejected', not browser.navigate('https://example.invalid/unregistered'))
            check(name + ': expanded text identical', browsers[0].toPlainText() == browsers[1].toPlainText())
            screenshot = output / (name + '.png')
            check(name + ': screenshot saved', browsers[0].grab().save(str(screenshot)))
            (output / (name + '.json')).write_text(browsers[0].export_json())
            reports.append({'name': name, 'kind': view.kind, 'task': view.task,
                'metrics': len(view.metrics), 'details': [d.key for d in view.details]})
            for browser in browsers:
                browser.clear()
                check(name + ': clear rejects old route', not browser.navigate('scorecard:expand'))
                try:
                    browser.export_json()
                except ValueError:
                    check(name + ': clear rejects old export', True)
                else:
                    check(name + ': clear rejects old export', False)
        finally:
            for browser in browsers:
                browser.close()

    from starplast import scorecard as SC
    check('all six task contracts represented', tasks == set(SC.TASKS))
    route_view = V.build_scorecard_view(exports['labels']['snapshot']['card'], links=(
        V.ScorecardLink('Fixture rows', 'scorecard:rows/fixture'),
        V.ScorecardLink('Fixture outcome', 'scorecard:outcome/fixture'),
        V.ScorecardLink('Authored source URL', 'https://example.invalid/source')))
    browser = ScorecardBrowser()
    source_calls, rows, outcomes = [], [], []
    original_open = QtGui.QDesktopServices.openUrl
    try:
        QtGui.QDesktopServices.openUrl = lambda url: source_calls.append(url.toString()) or True
        browser.set_scorecard(route_view)
        browser.rows_requested.connect(rows.append)
        browser.outcome_requested.connect(outcomes.append)
        for route in ('scorecard:rows/fixture', 'scorecard:outcome/fixture', 'https://example.invalid/source'):
            browser.anchorClicked.emit(QtCore.QUrl(route))
        check('real anchor signal dispatches registered rows', rows == ['fixture'])
        check('real anchor signal dispatches registered outcome', outcomes == ['fixture'])
        check('registered HTTPS dispatch (browser opening mocked)', source_calls == ['https://example.invalid/source'])
    finally:
        QtGui.QDesktopServices.openUrl = original_open
        browser.close()
    code = output / 'code'
    code.mkdir()
    paths = ['scripts/audit_desktop_scorecards.py', 'scripts/notebook_runner.py',
        'starplast/scorecard_view.py', 'starplast/scorecard_browser.py',
        'starplast/scorecard.py', 'starplast/host_foundation_scorecard.py', 'starplast/organisms.py']
    code_receipts = {}
    for relative in paths:
        path = ROOT / relative
        (code / path.name).write_bytes(path.read_bytes())
        code_receipts[relative] = _sha(path)
    summary = {'checks': len(checks), 'cards': reports, 'source_output_receipts': receipts,
        'code_sha256': code_receipts, 'qt_version': QtCore.QT_VERSION_STR,
        'runtime_seconds': time.monotonic() - started,
        'process_peak_rss_bytes': next(int(line.split()[1]) * 1024
            for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:')),
        'limitations': ['Authored task fixtures test presentation, not biological accuracy',
            'Registered external URL dispatch is mocked; no source was downloaded',
            'Human foundation is source evidence, not an installed host inference space']}
    (output / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def main():
    """Execute the desktop audit in a notebook and retain failed attempts."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a new immutable desktop audit directory')
    nb = ExecutedNotebook('Desktop scorecard acceptance')
    nb.md('Verify all six frozen task fixtures and the reviewed human evidence card in two actual desktop hosts. No model fitting or biological admission.')
    try:
        nb.code('from scripts.audit_desktop_scorecards import audit',
            f'summary = audit({str(args.out)!r})', "{'checks': summary['checks'], 'cards': len(summary['cards'])}")
    except Exception as exc:
        args.out.mkdir(parents=True, exist_ok=True)
        nb.md('Diagnostic: ' + type(exc).__name__ + ': ' + str(exc))
        nb.write(str(args.out / 'diagnostic.ipynb'))
        raise
    nb.write(str(args.out / 'audit.ipynb'))
    (args.out / 'README.md').write_text('Desktop presentation acceptance: all six task contracts, exact exports, definitions, full details, registered navigation and reviewed human source evidence. See summary.json for checks and limitations; external browser dispatch was mocked. Missing biology and calibration remain unavailable.\n')
    files = sorted(path for path in args.out.rglob('*') if path.is_file())
    (args.out / 'SHA256SUMS').write_text(''.join(_sha(path) + '  ' + path.relative_to(args.out).as_posix() + '\n' for path in files))
    print(json.dumps({'checks': nb.ns['summary']['checks'], 'cards': len(nb.ns['summary']['cards'])}))


if __name__ == '__main__':
    main()
