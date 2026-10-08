"""Replay three bounded view partitions from the verified frozen EC V4 artifact.

No source acquisition, fitting or prediction is performed. These partitions retain
the existing unresolved biological grade and uncalibrated method support.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from starplast import artifacts as A, functional_ontology as F, precompute_jobs as P, record_scorecards as R
from starplast.splits import read_split

SOURCE = ROOT / 'results/functional_ec_knn_pilot_2026_10_08_v4'
SOURCE_ID = '349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162'
BUDGET = P.Budget(max_jobs=3, max_seconds=120, max_storage_bytes=4 * 1024 * 1024,
    max_package_bytes=2 * 1024 * 1024, memory_limit_bytes=400 * 1024 * 1024)


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def _plain(value):
    return json.loads(json.dumps(value, allow_nan=False,
        default=lambda scalar: scalar.item() if isinstance(scalar, np.generic) else _unsupported(scalar)))


def _unsupported(value):
    raise ValueError('Unsupported scorecard scalar: ' + type(value).__name__)


def _source():
    manifest = json.loads((SOURCE / 'manifest.json').read_text())
    receipts = []
    for relative, expected_hash in sorted(manifest['output_sha256'].items()):
        path = SOURCE / relative
        if path.is_symlink() or not path.is_file() or _sha(path) != expected_hash:
            raise ValueError('Frozen source output changed: ' + relative)
        receipts.append({'path': str(path), 'sha256': expected_hash, 'bytes': path.stat().st_size})
    split, exclusions = read_split(SOURCE / 'split.json')
    spec = A.spec_from_dict(json.loads((SOURCE / 'held_out/manifest.json').read_text())['spec'])
    source = A.read_artifact(SOURCE / 'held_out', expected=spec, split=split)
    if source.identity != SOURCE_ID:
        raise ValueError('Frozen held-out identity differs from the verified V4 pilot')
    source.verify_contents()
    excluded = hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest()
    if A.Dependency('exclusions', 'training_source_closure', excluded) not in spec.dependencies:
        raise ValueError('Source exclusions differ from the artifact lineage')
    if len(source.payloads['rows.json']) != 182 or tuple(row['entity'] for row in source.payloads['rows.json']) != split.entities('test'):
        raise ValueError('The full ordered 182-row test population must be retained')
    if spec.fit_entities != split.entities('train') or len(spec.fit_entities) != 679:
        raise ValueError('Source fitting population differs from its frozen split')
    evaluation = json.loads(spec.evaluation_scope_json)
    if evaluation['truth_grade'] != 'unresolved' or evaluation['biological_admission'] is not False:
        raise ValueError('This replay cannot promote biological evidence grade')
    # Verify, but never fit or duplicate, the 7.3 MB native model state. Keep only
    # the small view payloads required for exact downstream comparisons.
    selected = {name: source.payloads[name] for name in ('rows.json', 'row_columns.json', 'class_scores.json',
        'card.json', 'profile_class_cards.json', 'major_class_cards.json')}
    return spec, split, exclusions, selected, receipts


def prepare(output):
    """Freeze replay scope, byte snapshots, code and budget before invoking builders."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new immutable operational replay directory')
    source_spec, split, exclusions, source_payloads, receipts = _source()
    for name in ('record_scorecards.py', 'scorecard.py', 'functional_ontology.py', 'artifacts.py', 'splits.py'):
        if _sha(ROOT / 'starplast' / name) != _sha(SOURCE / 'code' / name):
            raise ValueError('Replay code differs from the source pilot: ' + name)
    output.mkdir(parents=True)
    code_directory = output / 'code'
    code_directory.mkdir()
    code_paths = [Path(__file__), ROOT / 'scripts/notebook_runner.py', *(ROOT / 'starplast' / name for name in
        ('precompute_jobs.py', 'artifacts.py', 'splits.py', 'record_scorecards.py', 'scorecard.py',
         'functional_ontology.py', 'capabilities.py', 'strategies.py', 'query.py', 'provenance.py'))]
    snapshots = []
    for code in code_paths:
        frozen = code_directory / code.name
        frozen.write_bytes(code.read_bytes())
        snapshots.append(P.Snapshot(A.Dependency('code', 'replay:' + code.stem, _sha(frozen)), str(frozen.resolve())))
    source_manifest = SOURCE / 'manifest.json'
    receipts.append({'path': str(source_manifest), 'sha256': _sha(source_manifest), 'bytes': source_manifest.stat().st_size})
    for receipt in receipts:
        relative = str(Path(receipt['path']).relative_to(SOURCE))
        snapshots.append(P.Snapshot(A.Dependency('source', 'frozen_v4:' + relative, receipt['sha256']), receipt['path']))
    dependencies = source_spec.dependencies + (A.Dependency('artifact', 'frozen_ec_v4', SOURCE_ID),) + tuple(s.dependency for s in snapshots)
    templates = {}
    for name, role in (('held_out_rows', 'held_out'), ('scorecards', 'scorecard'), ('member_cards', 'scorecard')):
        settings = {**json.loads(source_spec.settings_json), 'precompute_view_operation': name,
            'replay_source_identity': SOURCE_ID}
        templates[name] = replace(source_spec, role=role, partition_id='precompute_ec_v4:' + name,
            settings_json=A.canonical_object(settings), dependencies=dependencies,
            code_version='frozen-ec-v4-view-replay-v1')
    score_spec = replace(templates['scorecards'], dependencies=dependencies + (
        A.Dependency('artifact', 'retained_held_out', SOURCE_ID),))
    member_spec = replace(templates['member_cards'], dependencies=dependencies + (
        A.Dependency('artifact', 'retained_held_out', SOURCE_ID), A.Dependency('artifact', 'matching_scorecard', SOURCE_ID)))
    plan = P.Plan((P.Job('held_out_rows', templates['held_out_rows'], split=split),
        P.Job('scorecards', score_spec, (P.Upstream('held_out_rows', 'artifact', 'retained_held_out'),), split),
        P.Job('member_cards', member_spec, (P.Upstream('held_out_rows', 'artifact', 'retained_held_out'),
            P.Upstream('scorecards', 'artifact', 'matching_scorecard')), split)), tuple(snapshots))
    scope = {'source_directory': str(SOURCE), 'source_artifact_identity': SOURCE_ID,
        'source_spec_identity': source_spec.identity, 'plan_identity': plan.identity,
        'job_fingerprints': dict(plan.fingerprints()), 'test_rows': 182, 'fit_entities': 679,
        'split_identity': split.identity, 'exclusion_identity': hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest(),
        'budget': asdict(BUDGET), 'first_run_max_jobs': 1, 'jobs': [j.name for j in plan.ordered()],
        'biological_admission': False, 'truth_grade': 'unresolved', 'confidence': source_spec.confidence_kind,
        'source_output_files_verified': len(receipts) - 1, 'snapshots': [asdict(s) for s in snapshots],
        'memory_measurement': 'Current process RSS at checkpoints; peak is the executed process lifetime, not per-job peak',
        'budget_scope': 'Runner byte limits cover typed artifacts; final entire evidence directory also checked against package limit',
        'runtime': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__},
        'limitations': ['Operational replay only; no fitting or new source acquisition',
            'Production feature/map/model/deployment builders and full coverage remain missing',
            'Cooperative budgets require the external MemoryMax400M process cap',
            *source_spec.gaps]}
    _write(output / 'scope.json', scope)
    _write(output / 'source_receipts.json', receipts)
    (output / 'split.json').write_bytes((SOURCE / 'split.json').read_bytes())
    return {'output': output, 'plan': plan, 'scope': scope, 'payloads': source_payloads, 'calls': []}


def _builders(session):
    def retained(context):
        session['calls'].append('held_out_rows')
        context.checkpoint()
        spec, split, exclusions, payloads, receipts = _source()
        return {name: payloads[name] for name in ('rows.json', 'row_columns.json', 'class_scores.json')} | {
            'scope.json': payloads['card.json']['scope']}

    def scorecards(context):
        session['calls'].append('scorecards')
        payloads = context.parents['held_out_rows'].payloads
        rows = pd.DataFrame(payloads['rows.json'], columns=payloads['row_columns.json'])
        scope_data = {**payloads['scope.json'], 'gaps': tuple(payloads['scope.json']['gaps'])}
        scope = R.RecordScope(**scope_data)
        scores = pd.DataFrame(payloads['class_scores.json']['scores'], columns=payloads['class_scores.json']['classes'])
        card = _plain(R.aggregate(rows, scope, parameters={'class_scores': scores}))
        profiles = _plain(R.class_cards(rows, scope))
        if card != session['payloads']['card.json'] or profiles != session['payloads']['profile_class_cards.json']:
            raise ValueError('Recomputed profile scorecards do not exactly match the source artifact')
        context.checkpoint()
        return {'card.json': card, 'profile_class_cards.json': profiles}

    def members(context):
        session['calls'].append('member_cards')
        rows_payload = context.parents['held_out_rows'].payloads
        rows = pd.DataFrame(rows_payload['rows.json'], columns=rows_payload['row_columns.json'])
        card = context.parents['scorecards'].payloads['card.json']
        if card['counts']['eligible'] != len(rows):
            raise ValueError('Member and profile cards must share the complete test cohort')
        member_cards = F.member_class_cards(rows)
        if member_cards != session['payloads']['major_class_cards.json']:
            raise ValueError('Recomputed member cards do not exactly match the source artifact')
        context.checkpoint()
        return {'major_class_cards.json': member_cards}

    return {'held_out_rows': retained, 'scorecards': scorecards, 'member_cards': members}


def stage(session, name):
    """Execute the stop/resume/no-builder sequence, asserting actual builder calls."""
    session['calls'].clear()
    budget = replace(BUDGET, max_jobs=1) if name == 'stop_after_one' else BUDGET
    result = P.run(session['plan'], {} if name == 'verified_replay' else _builders(session),
        session['output'] / 'run', budget)
    expected = {'stop_after_one': ['held_out_rows'], 'resume_remaining': ['scorecards', 'member_cards'], 'verified_replay': []}
    if session['calls'] != expected[name]:
        raise ValueError('Unexpected builder calls: ' + repr(session['calls']))
    if name == 'stop_after_one':
        assert dict(result.states) == {'held_out_rows': 'succeeded', 'scorecards': 'pending', 'member_cards': 'pending'}
        assert result.stop_reason == 'Run job attempt budget exhausted'
    else:
        assert set(result.states.values()) == {'succeeded'} and not result.stop_reason
    log = json.loads(result.journal.read_text())['journal']['runs'][-1]
    record = {'stage': name, 'states': dict(result.states), 'builder_calls': list(session['calls']),
        'artifact_identities': {key: artifact.identity for key, artifact in result.artifacts.items()}, 'measured_run': log}
    _write(session['output'] / (name + '.json'), record)
    session['result'] = result
    return record


def verify(session):
    """Compare retained rows and recomputed cards exactly, preserving frozen lineage."""
    outputs, source = session['result'].artifacts, session['payloads']
    exact = {'rows': outputs['held_out_rows'].payloads['rows.json'] == source['rows.json'],
        'row_columns': outputs['held_out_rows'].payloads['row_columns.json'] == source['row_columns.json'],
        'class_scores': outputs['held_out_rows'].payloads['class_scores.json'] == source['class_scores.json'],
        'aggregate_card': outputs['scorecards'].payloads['card.json'] == source['card.json'],
        'profile_cards': outputs['scorecards'].payloads['profile_class_cards.json'] == source['profile_class_cards.json'],
        'member_cards': outputs['member_cards'].payloads['major_class_cards.json'] == source['major_class_cards.json']}
    assert all(exact.values())
    for job in session['plan'].jobs:
        artifact = outputs[job.name]
        artifact.verify_contents()
        assert artifact.spec.fit_entities == session['plan'].jobs[0].spec.fit_entities
        assert artifact.spec.evaluation_scope_json == session['plan'].jobs[0].spec.evaluation_scope_json
        assert A.Dependency('artifact', 'frozen_ec_v4', SOURCE_ID) in artifact.spec.dependencies
        for upstream in job.upstream:
            assert A.Dependency(upstream.kind, upstream.name, outputs[upstream.job].identity) in artifact.spec.dependencies
    summary = {'source_artifact_identity': SOURCE_ID, 'test_rows_retained': len(source['rows.json']),
        'exact_matches': exact, 'member_classes': len(outputs['member_cards'].payloads['major_class_cards.json']),
        'counts': source['card.json']['counts'], 'artifact_identities': {key: artifact.identity for key, artifact in outputs.items()},
        'biological_admission': False, 'instruction_64_25_complete': False,
        'remaining': 'Production feature/map/model/deployment builders and full target/strategy coverage'}
    _write(session['output'] / 'summary.json', summary)
    return summary


def main():
    """Write an actually executed bounded replay notebook and checksummed evidence."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('EC V4 operational precomputation: bounded stop, resume and verified replay')
    nb.md('Freeze the already verified EC V4 held-out artifact and all consumed file/code byte identities before any job runs. '
        'The original 182 test rows, 679 fitting entities, split, exclusions, unresolved evidence grade and uncalibrated support remain recorded. '
        'No fitting, new source acquisition or biological promotion occurs. External process cap: MemoryMax400M. '
        'Declared per-run budget: 120 seconds, three job attempts, 4 MiB typed-artifact storage, 2 MiB package, 400 MiB current RSS.')
    try:
        nb.code('from scripts.freeze_precompute_ec_views import prepare, stage, verify',
            f'session = prepare({str(args.out.resolve())!r})', "session['scope']")
        nb.md('Run with max_jobs=1. Only the independently reusable held-out row partition should complete; descendants remain pending.')
        nb.code("first = stage(session, 'stop_after_one')", 'first')
        nb.md('Resume verified completed rows; build only scorecard and member-card partitions. Each downstream dependency binds complete parent content identity.')
        nb.code("second = stage(session, 'resume_remaining')", 'second')
        nb.md('Replay with an empty builder mapping. Every partition must be read and verified from its current spec/split/payload/journal identity.')
        nb.code("third = stage(session, 'verified_replay')", 'third')
        nb.md('Exact comparison to existing source cards; retained abstentions and unknown biological accuracy. Instruction 64.25 remains partial.')
        nb.code('verification = verify(session)', 'verification')
    except Exception as error:
        nb.md('Diagnostic: ' + type(error).__name__ + ': ' + str(error))
        nb.write(str(args.out / 'diagnostic.ipynb'))
        raise
    nb.write(str(args.out / 'pilot.ipynb'))
    files = {str(path.relative_to(args.out)): {'sha256': _sha(path), 'bytes': path.stat().st_size}
        for path in args.out.rglob('*') if path.is_file()}
    package_bytes = sum(record['bytes'] for record in files.values())
    if package_bytes >= BUDGET.max_package_bytes:
        raise ValueError('Entire operational evidence directory exceeds the declared package budget')
    _write(args.out / 'manifest.json', {'files': files, 'evidence_bytes_before_manifest': package_bytes,
        'package_budget_bytes': BUDGET.max_package_bytes})
    assert sum(path.stat().st_size for path in args.out.rglob('*') if path.is_file()) <= BUDGET.max_package_bytes
    print(json.dumps(nb.ns['verification'], indent=2))
    print('Evidence bytes including manifest:', sum(path.stat().st_size for path in args.out.rglob('*') if path.is_file()))


if __name__ == '__main__':
    main()
