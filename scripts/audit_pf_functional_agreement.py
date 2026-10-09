"""Independently check the frozen paired report against original source records."""
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import organisms as _ORGANISMS
PACKET = ROOT / 'results/pf_functional_agreement_2026_10_08_v2'
OUTPUT = ROOT / 'results/pf_functional_agreement_acceptance_2026_10_08_v2'
REPORT_SHA = '940da900b3d5ba356c6e360f220787b849d7382138132ccc96f2d1684178e98b'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def audit():
    from starplast import functional_results as F
    assert sha(PACKET / 'report.json') == REPORT_SHA
    receipts = 0
    for line in (PACKET / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1); path = PACKET / name
        assert not path.is_symlink() and path.resolve().is_relative_to(PACKET.resolve())
        assert sha(path) == digest, name
        receipts += 1
    inputs = json.loads((PACKET / 'input_manifest.json').read_text())
    for name, digest in inputs.items(): assert sha(name) == digest, name
    entries, reason = F.shipped(_ORGANISMS.FALCIPARUM); assert not reason
    pair = [next(b for b in entries if b.metadata['strategy'] == strategy)
            for strategy in ('feature_knn', 'random_forest')]
    report = json.loads((PACKET / 'report.json').read_text())
    a, b = (p.rows for p in pair)
    assert len(a) == len(b) == len(report['rows']) == 152
    for left, right, retained in zip(a, b, report['rows']):
        for key in ('entity', 'truth', 'group', 'training_supported'):
            assert left[key] == right[key] == retained[key]
        assert left['prediction'] == retained['left_prediction']
        assert right['prediction'] == retained['right_prediction']
    calls = [(x['prediction'], y['prediction'], x['truth']) for x, y in zip(a, b)]
    expected = {
        'both_abstain': sum(x is None and y is None for x, y, t in calls),
        'left_only_correct': sum(x == t and y is None for x, y, t in calls),
        'left_only_wrong': sum(x is not None and x != t and y is None for x, y, t in calls),
        'right_only_correct': sum(x is None and y == t for x, y, t in calls),
        'right_only_wrong': sum(x is None and y is not None and y != t for x, y, t in calls),
        'agree_correct': sum(x == y == t for x, y, t in calls),
        'agree_wrong': sum(x is not None and x == y and x != t for x, y, t in calls),
        'conflict_left_correct': sum(x == t and y is not None and y != x for x, y, t in calls),
        'conflict_right_correct': sum(y == t and x is not None and x != y for x, y, t in calls),
        'conflict_both_wrong': sum(x is not None and y is not None and x != y and x != t and y != t
                                   for x, y, t in calls),
        'eligible': len(calls), 'joint_answered': sum(x is not None and y is not None for x, y, t in calls),
        'agreements': sum(x is not None and x == y for x, y, t in calls),
        'conflicts': sum(x is not None and y is not None and x != y for x, y, t in calls),
        'shared_wrong': sum(x is not None and y is not None and x != t and y != t for x, y, t in calls),
        'unsupported': sum(not row['training_supported'] for row in a),
        'groups': len({row['group'] for row in a})}
    assert report['counts'] == expected
    for name, n, d in (
        ('joint_coverage', expected['joint_answered'], 152),
        ('agreement_among_joint_calls', expected['agreements'], expected['joint_answered']),
        ('source_recovery_among_agreements', expected['agree_correct'], expected['agreements']),
        ('correct_agreements_all_eligible', expected['agree_correct'], 152),
        ('shared_wrong_among_joint_calls', expected['shared_wrong'], expected['joint_answered'])):
        assert report['rates'][name] == {'numerator': n, 'denominator': d, 'value': n / d if d else None}
    for method, original in zip(report['methods'], pair):
        assert method['counts'] == original.card['counts']
        assert method['artifact_identity'] == original.metadata['source_artifact_identity']
        assert method['settings'] == original.metadata['source_manifest']['spec']['settings_json']
    assert json.loads((PACKET / 'scorecard.json').read_text())['snapshot'] == report
    assert report['biological_accuracy'] is None and report['calibrated_confidence'] is None
    assert report['overlap']['independent_evidence'] is False and report['biological_admission'] is False
    return {'packet_outputs_verified': receipts, 'current_inputs_verified': len(inputs),
            'all_counts_and_rates_reproduced_exactly': True, 'test_genes': 152,
            'unsupported_test_genes': 4, 'agreements': 94, 'correct_agreements': 33,
            'same_wrong_agreements': 61, 'both_wrong': 78, 'group_uncertainty': 'unavailable'}


if __name__ == '__main__':
    from scripts.notebook_runner import ExecutedNotebook
    assert OUTPUT.is_dir() and not any(OUTPUT.iterdir()), 'Use an untouched acceptance directory'
    shutil.copyfile(__file__, OUTPUT / Path(__file__).name)
    nb = ExecutedNotebook('Independent frozen paired-record verification')
    nb.ns.update(audit=audit)
    nb.md('Verify every packet receipt and current input. Recount directly from the pinned '
          'original method rows, without invoking the comparison implementation. No fitting.')
    try:
        nb.code('result=audit()', 'result')
        (OUTPUT / 'acceptance.json').write_text(json.dumps(nb.ns['result'], indent=2) + '\n')
        print(json.dumps(nb.ns['result']))
    finally:
        nb.write(str(OUTPUT / 'executed.ipynb'))
        (OUTPUT / 'SHA256SUMS.txt').write_text(''.join(
            f'{sha(path)}  {path.relative_to(OUTPUT)}\n' for path in sorted(OUTPUT.rglob('*'))
            if path.is_file() and path.name != 'SHA256SUMS.txt'))
