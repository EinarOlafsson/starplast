"""Verify packaged enzyme and domain recovery through actual Discoveries routes."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def audit(output, organism='Tg'):
    from PyQt6 import QtWidgets
    from starplast import functional_results as F, strategies as S
    from starplast.discoveries_panel import DiscoveriesPanel
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    code = [Path(__file__), ROOT / 'scripts/notebook_runner.py',
        *[ROOT / 'starplast' / name for name in ('functional_results.py', 'functional_class_view.py',
            'discoveries_panel.py', 'functional_coverage.py', 'scorecard_view.py', 'functional_domains.py')]]
    sources = [ROOT / 'starplast/data' / name for name in ('nodes.parquet', 'pf_nodes.parquet',
        'claims.parquet', 'claim_recipes.parquet', 'track_record.parquet',
        'functional_domain_names.json', 'functional_results.json')]
    before = {str(path): sha(path) for path in code + sources}
    (output / 'code').mkdir()
    for path in code:
        (output / 'code' / path.name).write_bytes(path.read_bytes())
    (output / 'input_manifest.json').write_text(json.dumps(before, indent=2) + '\n')
    benchmarks, reason = F.shipped(organism)
    expected_namespaces={'ec_major','pfam'} if organism=='Tg' else {'ec_major'}
    assert not reason and {item.namespace for item in benchmarks} == expected_namespaces
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = DiscoveriesPanel(organism, context=S.Context.shipped(organism))
    summaries = []
    chosen = []
    panel.gene_chosen.connect(chosen.append)
    try:
        assert not panel.coverage_warning and len(panel.functional_benchmarks) == len(benchmarks)
        for benchmark in benchmarks:
            source = benchmark.metadata['source_targets'][0]
            panel.label.setCurrentIndex(panel.label.findData(source))
            assert panel.annotation_functional.isEnabled()
            panel.annotation_functional.click()
            assert panel.tabs.currentWidget() == panel.functional_page
            assert panel.functional_choice.currentData().metadata == benchmark.metadata
            assert panel.functional_shown_rows == benchmark.rows
            assert dict(panel.functional_card.view.counts) == benchmark.card['counts']
            panel.functional_view.setCurrentIndex(panel.functional_view.findData('strategy:'))
            assert panel.functional_shown_rows == benchmark.rows
            for name, card in benchmark.baseline_cards.items():
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('baseline:' + name))
                exported = json.loads(panel.functional_card.export_json())['snapshot']
                assert exported['card'] == card
                assert panel.functional_rows.rowCount() == 0
            selected = []
            for predicate in (lambda card: card['true_positive'] > 0,
                              lambda card: card['false_positive'] > 0,
                              lambda card: card['precision'] is None and card['false_negative'] > 0):
                matches = [card for card in benchmark.major_class_cards if predicate(card)]
                if matches and matches[0] not in selected:
                    selected.append(matches[0])
            assert selected
            for card in selected:
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('major:' + card['class']))
                values = [panel.functional_rows.item(i, 4).text() for i in range(panel.functional_rows.rowCount())]
                assert values.count('recovered reference class') == card['true_positive']
                assert values.count('unexpected reference-class call') == card['false_positive']
                assert sum('class missed' in item for item in values) == card['false_negative']
                exported = json.loads(panel.functional_card.export_json())['snapshot']
                assert exported['details']['native_class_card'] == card
                metrics = {metric.key: metric.value for metric in panel.functional_card.view.metrics}
                for key in ('precision', 'recall', 'f1'):
                    assert metrics['class_' + key] == card[key]
                assert metrics['coverage'] == card['coverage']
                if benchmark.namespace == 'pfam':
                    assert 'domain' in panel.functional_card.view.title.lower()
            profile = next(card for card in benchmark.profile_class_cards if card['counts']['abstained'] > 0)
            panel.functional_view.setCurrentIndex(panel.functional_view.findData('profile:' + profile['class']))
            assert json.loads(panel.functional_card.export_json())['snapshot']['card'] == profile
            assert any(row['prediction'] is None for row in panel.functional_shown_rows)
            panel.functional_view.setCurrentIndex(0)
            panel._functional_gene_clicked(0, 0)
            assert chosen[-1] == benchmark.rows[0]['entity']
            panel.annotation_coverage.click()
            index = next(i for i, row in enumerate(panel.coverage_shown)
                if row['reference_recovery'].get('source_artifact_identity') == benchmark.metadata['source_artifact_identity'])
            panel._coverage_selected(index, 0)
            assert panel.coverage_open.isEnabled()
            panel.coverage_open.click()
            assert panel.functional_choice.currentData().metadata == benchmark.metadata
            assert panel.functional_shown_rows == benchmark.rows
            panel.resize(1400, 950)
            panel.show()
            app.processEvents()
            assert panel.grab().save(str(output / (benchmark.namespace + '.png')))
            summaries.append({'namespace': benchmark.namespace, 'source': source,
                'counts': benchmark.card['counts'], 'selected_membership_cards': len(selected),
                'controls': list(benchmark.baseline_cards), 'artifact_identity': benchmark.metadata['source_artifact_identity'],
                'routes': ['annotation', 'strategy', 'control', 'membership', 'profile', 'gene', 'coverage'],
                'independent_biology': None, 'calibration': None})
        assert all(not F.shipped(host)[0] and F.shipped(host)[1] for host in ('Hs','Mm'))
        changed = panel.context.nodes.copy()
        changed.loc[changed.index[0], 'gene_id'] = 'altered_audit_gene'
        invalid = S.Context(changed, graph={}, organism=organism)
        for benchmark in benchmarks:
            try:
                F.require_context(benchmark, invalid)
            except ValueError:
                pass
            else:
                raise AssertionError('Altered source borrowed archived accuracy')
            other=S.Context(panel.context.nodes,graph={},organism='Pf' if organism=='Tg' else 'Tg')
            try:F.require_context(benchmark,other)
            except ValueError:pass
            else:raise AssertionError('Other organism borrowed archived accuracy')
    finally:
        panel.close()
    assert before == {str(path): sha(path) for path in code + sources}
    summary = {'benchmarks': summaries, 'elapsed_seconds': time.monotonic() - start,
        'process_peak_bytes': next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines()
            if line.startswith('VmHWM:')), 'original_inputs_unchanged': True,
        'organism':organism,'hosts_unavailable':True,'cross_organism_refused':True,'altered_context_refused':True}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def main():
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--organism', choices=('Tg','Pf'),default='Tg')
    args=parser.parse_args();output=args.out
    nb = ExecutedNotebook('Exact domain and enzyme Discoveries scorecard navigation')
    nb.ns.update(audit=audit, output=output,organism=args.organism)
    nb.md('Use original held-out artifacts and unchanged installed node tables; no fit, source promotion or biological admission. Verify class populations, controls, abstentions and all displayed routes.')
    try:
        nb.code('summary = audit(output,organism)', 'summary')
        print(json.dumps(nb.ns['summary'], indent=2))
    except Exception as exc:
        output.mkdir(parents=True, exist_ok=True)
        (output / 'failure.json').write_text(json.dumps({'error': type(exc).__name__, 'detail': str(exc)}) + '\n')
        nb.md('Execution failed: ' + type(exc).__name__ + ': ' + str(exc))
        raise
    finally:
        output.mkdir(parents=True, exist_ok=True)
        nb.write(str(output / 'executed.ipynb'))
        (output / 'SHA256SUMS.txt').write_text('\n'.join(f'{sha(path)}  {path.relative_to(output)}'
            for path in sorted(output.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
