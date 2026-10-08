"""Execute scoped coverage, shared scorecard and host status navigation checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(output):
    """Retain exact source/card replay and real UI routes in a fresh snapshot."""
    from PyQt6 import QtWidgets
    from scripts.audit_functional_discoveries import audit as source_audit
    from starplast import functional_coverage as FC, functional_results as F, host_foundation_scorecard as HF
    from starplast import organisms as O, strategies as S
    from starplast.discoveries_panel import DiscoveriesPanel
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable coverage UI directory')
    output.mkdir(parents=True)
    source_summary=source_audit(output/'source_replay')
    app=QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    summaries=[]
    for organism in (O.TOXOPLASMA,O.FALCIPARUM):
        panel=DiscoveriesPanel(organism,context=S.Context.shipped(organism))
        try:
            assert not panel.coverage_warning
            expected=FC.build({organism:panel.inventory},benchmarks=F.shipped(organism)[0],organisms=(organism,))
            assert panel.functional_coverage==expected
            panel.label.setCurrentIndex(panel.label.findData('ec_number'))
            panel.annotation_coverage.click()
            assert panel.tabs.currentWidget()==panel.coverage_page
            assert all(row['source_label']=='ec_number' and row['organism']==organism for row in panel.coverage_shown)
            available=[i for i,row in enumerate(panel.coverage_shown)
                if row['reference_recovery']['status']=='verified_reference_recovery']
            assert len(available)==(1 if organism==O.TOXOPLASMA else 0)
            panel._coverage_selected(available[0] if available else 0,0)
            assert panel.coverage_open.isEnabled()==bool(available)
            panel.resize(1400,1000);panel.show();app.processEvents()
            assert panel.grab().save(str(output/(organism+'_coverage.png')))
            if available:
                panel.coverage_open.click()
                assert panel.tabs.currentWidget()==panel.functional_page
                benchmark=panel.functional_choice.currentData()
                view=panel.functional_card.view
                snapshot=json.loads(view.snapshot_json)
                assert snapshot['card']==benchmark.card
                supplied={metric.key:metric.value for metric in view.metrics}
                assert all(supplied[key]==value for key,value in benchmark.card['metrics'].items())
                assert dict(view.counts)==benchmark.card['counts']
                assert panel.functional_shown_rows==benchmark.rows
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('major:2'))
                assert panel.functional_card.view is not None
                assert 'Class precision, whole cohort' in panel.functional_card.toPlainText()
                assert panel.functional_card.navigate('scorecard:metric/class_precision')
                assert 'Whole-cohort calls of this class' in panel.functional_card.toPlainText()
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('strategy:'))
                panel._functional_outcome(0,0)
                panel.functional_outcome_dialog.close()
                assert panel.functional_card.navigate('scorecard:metric/precision_of_calls')
                assert 'Calls made: 171' in panel.functional_card.toPlainText()
                assert panel.functional_card.navigate('scorecard:expand')
                assert 'Calibrated confidence: unavailable' in panel.functional_card.toPlainText()
                exported=json.loads(panel.functional_card.export_json())
                assert exported['snapshot']['card']==benchmark.card
                (output/'functional_card_export.json').write_text(panel.functional_card.export_json())
                app.processEvents();assert panel.grab().save(str(output/'functional_full_record.png'))
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('baseline:majority'))
                assert panel.functional_shown_rows==[]
                assert json.loads(panel.functional_card.view.snapshot_json)['card']==benchmark.baseline_cards['majority']
            panel.tabs.setCurrentWidget(panel.host_page)
            assert panel.host_card.view.kind=='evidence' and not panel.host_card.view.metrics
            host_snapshot=json.loads(panel.host_card.view.snapshot_json)
            assert dict(panel.host_card.view.counts)['canonical_genes']==58988
            assert 'not validated HFF' in panel.host_card.toPlainText()
            assert _sha(ROOT/'starplast/data/host_foundation_summary.json')==HF.SUMMARY_SHA256
            app.processEvents();assert panel.grab().save(str(output/(organism+'_host_status.png')))
            if organism==O.TOXOPLASMA:(output/'host_card_export.json').write_text(panel.host_card.export_json())
            (output/(organism+'_coverage.json')).write_text(json.dumps(expected,indent=2,sort_keys=True)+'\n')
            summaries.append({'organism':organism,'coverage_addresses':len(expected['rows']),
                'reference_recovery_addresses':len(available),'independent_biological_admissions':0,
                'calibrated_deployments':0,'host_view_organism':O.HUMAN,'host_registered':False,
                'shared_card_values':'exact native snapshot','scope_navigation':'verified'})
        finally:panel.close()
    code=[Path(__file__),ROOT/'scripts/audit_functional_discoveries.py',
        *(ROOT/'starplast'/name for name in ('discoveries_panel.py','functional_coverage.py',
            'functional_exclusions.py','functional_domain_profiles.py','functional_results.py',
            'scorecard.py','scorecard_view.py','scorecard_browser.py','host_foundation_scorecard.py',
            'functional_class_view.py'))]
    (output/'code').mkdir()
    for path in code:(output/'code'/path.name).write_bytes(path.read_bytes())
    (output/'summary.json').write_text(json.dumps({'ui':summaries,'source_replay':source_summary},indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps({'input_sha256':{str(path):_sha(path) for path in code},
        'output_sha256':{str(path.relative_to(output)):_sha(path) for path in output.rglob('*') if path.is_file()},
        'interpretation':'Source recovery and evidence quality only; no independent biological admission or host space registration'},indent=2)+'\n')
    return summaries


def main():
    """Write executed evidence or preserve a diagnostic without overwriting earlier work."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    nb=ExecutedNotebook('Functional strategy coverage and shared scorecards')
    nb.md('Replay every original source membership and functional result. Inspect actual scoped coverage routes, exact card values, definitions, full test record and export. A separate human host evidence card retains source and mapping counts, with no performance metrics or admission. Original snapshots stay immutable.')
    try:
        nb.code('from scripts.audit_functional_coverage_ui import audit',f'summary=audit({str(args.out)!r})','summary')
    except Exception as exc:
        nb.md('Diagnostic: '+type(exc).__name__+': '+str(exc))
        nb.write(str(args.out/'diagnostic.ipynb'));raise
    nb.write(str(args.out/'audit.ipynb'))
    print(json.dumps(nb.ns['summary'],indent=2))


if __name__=='__main__':main()
