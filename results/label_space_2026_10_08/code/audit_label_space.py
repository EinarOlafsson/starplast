"""Execute native annotation and held-out label/class navigation acceptance."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def audit(output):
    """Replay two actual spaces and independently count native class errors."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    os.environ['STARPLAST_STATE'] = str(output / 'state')
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PyQt6 import QtCore, QtWidgets
    import numpy as np
    from starplast import app as A, organisms as O, paths, track_record as T, strategies as S
    from starplast.label_space import LabelSpace

    started = time.monotonic()
    code = ['starplast/label_space.py', 'starplast/label_space_browser.py', 'starplast/app.py',
        'starplast/discovery_labels.py', 'starplast/domain_profiles.py', 'starplast/query.py',
        'starplast/record_scorecards.py', 'starplast/scorecard.py', 'starplast/scorecard_view.py',
        'starplast/scorecard_browser.py', 'starplast/track_record.py', 'starplast/strategy_capabilities.py',
        'scripts/audit_label_space.py', 'scripts/notebook_runner.py']
    inputs = [ROOT / name for name in code] + [Path(O.nodes_path(c)) for c in O.codes()]
    inputs.append(Path(paths.cache_file('track_record.parquet')))
    hashes = {str(p): sha(p) for p in inputs if p.is_file()}
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    checks, populations = [], []

    def check(name, condition):
        assert condition, name
        checks.append(name)

    for organism in O.codes():
        window = A.Window(species=O.get(organism).species)
        try:
            panel = window.open_label_browser()
            model = panel.model
            check('all 39 mechanisms visible', panel.strategy_choice.count() == len(S.catalog()) == 39)
            for label in model.labels.to_dict('records'):
                target = label['target']
                members, classes = model.members(target), model.classes(target)
                check('native membership stays organism scoped', members.empty or set(members.organism) == {organism})
                check('every native class balance exact', all(int(row.annotated_genes) == members.loc[members.value.eq(row.value), 'gene_id'].nunique() for row in classes.itertuples()))
                check('whole-universe unknown count exact', int(label['genes']) == len(window.nodes) and int(label['unannotated_genes']) == len(window.nodes) - members.gene_id.nunique())
                saved = json.loads(model.evidence_card(target).snapshot_json)['card']['evidence']
                check('unknown hierarchy and unmeasured classes explicit', saved['hierarchy']['status'] == saved['unmeasured_classes']['status'] == 'unavailable')
            target = T.default_target(organism)
            ledger = T.shipped(organism)
            observed = ledger.loc[(ledger.target == target) & (ledger['mode'] == 'together'), 'truth'].astype(str)
            value = next(value for value in model.classes(target).value if value in set(observed))
            panel.select(target, value)
            records = model.mechanisms(target, value)
            evaluations = 0
            for record in records:
                for evaluation in record['evaluations']:
                    scope = evaluation['scope']
                    mask = ledger.organism.eq(organism) & ledger.target.eq(target) & ledger.strategy.eq(record['strategy'])
                    mask &= ledger.setting_key.astype(str).eq(scope['settings']) & ledger.seed.eq(scope['seed'])
                    partition = ledger['mode'].astype(str) + ':' + ledger.set_name.fillna('').astype(str)
                    native = ledger[mask & partition.eq(scope['partition'])]
                    called = ~native.abstained.astype(bool) & native.prediction.astype(str).eq(value)
                    positive = native.truth.astype(str).eq(value)
                    tp, fp, fn = int((called & positive).sum()), int((called & ~positive).sum()), int((~called & positive).sum())
                    metrics = evaluation['class_metrics']
                    check('native class TP/FP/FN exact', (metrics['true_positive'], metrics['false_positive'], metrics['false_negative']) == (tp, fp, fn))
                    precision = tp / (tp + fp) if tp + fp else None
                    recall = tp / (tp + fn) if tp + fn else None
                    check('native precision and recall exact', metrics['precision'] == precision and metrics['recall'] == recall)
                    evaluations += 1
            check('actual held-out evaluations reached', evaluations > 0)
            index = next(i for i in range(panel.strategy_choice.count()) if panel.strategy_choice.itemData(i)['evaluations'])
            panel.strategy_choice.setCurrentIndex(index)
            check('immutable exact cohort card shown', panel.strategy_card.view == panel.evaluations[panel.evaluation_choice.currentData()]['card'])
            check('confusions accessible', bool(panel.evaluation_details.toPlainText()))
            changed = window.nodes.copy(deep=True)
            changed.loc[changed.index[0], target] = 'authored changed context'
            altered = LabelSpace(changed, organism, ledger=ledger, ledger_nodes=window.nodes)
            check('changed table cannot borrow accuracy', not any(row['evaluations'] for row in altered.mechanisms(target)))
            if 'is_exported' in set(model.labels.target):
                row = model.classes('is_exported')
                check('native False identity retained', any(type(value) in (bool, np.bool_) and not value for value in row.native_value))
            gene = str(panel.member_frame.gene_id.iloc[0])
            panel._open_gene(0, 0)
            check('class to gene keeps organism', window.nodes.gene_id.iloc[window.sel] == gene)
            panel.resize(1400, 950)
            panel.show()
            application.processEvents()
            panel.grab().save(str(output / (organism + '_class.png')))
            populations.append({'organism': organism, 'labels': len(model.labels), 'target': target, 'class': value, 'evaluations': evaluations})
        finally:
            window.close()
            window.deleteLater()
            QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete)
            application.processEvents()
    for path, digest in hashes.items():
        check('input unchanged: ' + path, sha(path) == digest)
    snapshot = output / 'code'
    snapshot.mkdir()
    for name in code:
        if (ROOT / name).is_file():
            (snapshot / Path(name).name).write_bytes((ROOT / name).read_bytes())
    summary = {'checks': len(checks), 'populations': populations, 'input_sha256': hashes,
        'runtime_seconds': time.monotonic() - started,
        'process_peak_rss_bytes': next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:')),
        'limitations': ['Legacy held-out arithmetic is source-annotation recovery; independent biology remains unknown',
            'Derived EC benchmark scopes remain separate from raw EC classes', 'Unknown ontology hierarchy/unmeasured classes remain unavailable',
            'No fitting, downloads, source promotion or host installation']}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (output / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    return summary


def main():
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a fresh immutable directory')
    notebook = ExecutedNotebook('Label/class native acceptance')
    try:
        notebook.code('from scripts.audit_label_space import audit', f'summary = audit({str(args.out)!r})', "summary['checks']")
    except Exception as error:
        args.out.mkdir(parents=True, exist_ok=True)
        notebook.md(type(error).__name__ + ': ' + str(error))
        notebook.write(str(args.out / 'diagnostic.ipynb'))
        raise
    notebook.write(str(args.out / 'audit.ipynb'))
    files = sorted(path for path in args.out.rglob('*') if path.is_file())
    (args.out / 'SHA256SUMS').write_text(''.join(sha(path) + '  ' + str(path.relative_to(args.out)) + '\n' for path in files))
    print(json.dumps({'checks': notebook.ns['summary']['checks']}))


if __name__ == '__main__':
    main()
