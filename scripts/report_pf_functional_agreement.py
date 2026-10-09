"""Freeze a descriptive comparison of the accepted Pf kNN and forest tests."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import organisms as _ORGANISMS


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run(output):
    from starplast import functional_agreement as A, functional_results as F
    from starplast.scorecard_view import export_scorecard, render_scorecard_html
    benchmarks, reason = F.shipped(_ORGANISMS.FALCIPARUM)
    assert not reason
    pair = [next(b for b in benchmarks if b.metadata['strategy'] == strategy)
            for strategy in ('feature_knn', 'random_forest')]
    inputs = {str(ROOT / 'starplast/data/functional_results.json'): F.FUNCTIONAL_RESULTS_SHA256}
    for path in (Path(__file__), ROOT / 'starplast/functional_agreement.py',
                 ROOT / 'starplast/functional_results.py', ROOT / 'starplast/scorecard_view.py',
                 ROOT / 'tests/test_functional_agreement.py', ROOT / 'scripts/notebook_runner.py'):
        inputs[str(path)] = sha(path)
        destination = output / 'code' / path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    models = []; source_outputs = 0; artifact_files = 0
    for benchmark in pair:
        source = ROOT / benchmark.metadata['source_directory']
        manifest_path = source / 'held_out/manifest.json'
        manifest = json.loads(manifest_path.read_text())
        assert manifest == benchmark.metadata['source_manifest']
        inputs[str(manifest_path)] = sha(manifest_path)
        sums = source / 'SHA256SUMS.txt'
        inputs[str(sums)] = sha(sums)
        for line in sums.read_text().splitlines():
            digest, name = line.split('  ', 1); path = source / name
            assert not path.is_symlink() and path.resolve().is_relative_to(source.resolve())
            assert sha(path) == digest, name
            source_outputs += 1
        for name, receipt in manifest['files'].items():
            path = source / 'held_out' / name
            assert sha(path) == receipt['sha256'] and path.stat().st_size == receipt['bytes']
            artifact_files += 1
        model_path = source / 'held_out/model_state.json'
        inputs[str(model_path)] = manifest['files']['model_state.json']['sha256']
        models.append(json.loads(model_path.read_text()))
        for dependency in manifest['spec']['dependencies']:
            if dependency['kind'] == 'table':
                path = ROOT / 'starplast/data/pf_nodes.parquet'
                assert dependency['name'] == 'installed_nodes'
                inputs[str(path)] = dependency['sha256']
    report = A.compare(*pair)
    for name in ('input_columns', 'kept_columns'):
        columns = [model[name] for model in models]
        assert all(len(value) == len(set(value)) for value in columns)
        report['overlap'][name] = {'left_count': len(columns[0]), 'right_count': len(columns[1]),
            'shared': sorted(set(columns[0]) & set(columns[1])),
            'left_only': sorted(set(columns[0]) - set(columns[1])),
            'right_only': sorted(set(columns[1]) - set(columns[0]))}
    report['overlap']['selected_feature_overlap'] = 'Exact native model input/retained column sets, verified against model receipts'
    assert report['counts']['eligible'] == 152 and report['counts']['unsupported'] == 4
    assert report['counts']['joint_answered'] == 143
    assert all(sha(path) == digest for path, digest in inputs.items())
    view = A.build_scorecard(report)
    (output / 'report.json').write_text(F._canonical(report) + '\n')
    (output / 'scorecard.json').write_text(export_scorecard(view))
    (output / 'scorecard.html').write_text(render_scorecard_html(view, expanded=True))
    assert json.loads((output / 'scorecard.json').read_text())['snapshot'] == report
    (output / 'input_manifest.json').write_text(json.dumps(inputs, indent=2) + '\n')
    result = {'source_outputs_verified': source_outputs, 'artifact_files_verified': artifact_files,
              'inputs_verified': len(inputs), 'counts': report['counts'], 'rates': report['rates'],
              'biological_admission': False, 'calibrated_confidence': None,
              'group_uncertainty': 'unavailable', 'fits_performed': 0}
    (output / 'acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=False)
    scope = {'organism': _ORGANISMS.FALCIPARUM, 'strategies': ['feature_knn', 'random_forest'],
             'target': 'ec_direct_complete_major_profile', 'test_genes': 152,
             'unsupported_test_genes': 4, 'fits_performed': 0,
             'biological_admission': False, 'group_uncertainty': 'unavailable'}
    (args.out / 'predeclared_scope.json').write_text(json.dumps(scope, indent=2) + '\n')
    nb = ExecutedNotebook('Paired Pf source-profile agreements and errors')
    nb.ns.update(run=run, output=args.out)
    nb.md('Compare the previously frozen 152 genes, preserving complete profiles and abstentions.',
          'Shared labels, training genes and inputs prevent an independence claim. '
          'Source recovery is not biological accuracy. No group-aware interval, '
          'confidence calibration, model selection, deployment or fitting is performed.')
    try:
        nb.code('result=run(output)', 'result')
        print(json.dumps(nb.ns['result'], indent=2))
    except Exception as exc:
        (args.out / 'failure.json').write_text(json.dumps({'error': type(exc).__name__, 'detail': str(exc)}) + '\n')
        raise
    finally:
        nb.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text(''.join(
            f'{sha(path)}  {path.relative_to(args.out)}\n' for path in sorted(args.out.rglob('*'))
            if path.is_file() and path.name != 'SHA256SUMS.txt'))


if __name__ == '__main__': main()
