"""Freeze legacy row-versus-scorecard arithmetic, without claiming new biological validation."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import record_scorecards as R, scorecard as SC, track_record as T  # noqa: E402


def reconcile(output):
    """Check every legacy cohort and preserve the two accuracy denominators separately."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new reconciliation snapshot directory')
    output.mkdir(parents=True)
    source = ROOT / 'starplast/data/track_record.parquet'
    ledger = pd.read_parquet(source)
    reports, cards, seen = [], [], set()
    for scope, rows in R.legacy_cohorts(ledger):
        card = R.aggregate(rows, scope)
        standard = SC.label_calls(rows.prediction.astype(object).where(rows.prediction.notna(), None), rows.truth, range(len(rows)))
        assert card['metrics'] == R._finite(standard)
        legacy = T.summary(rows.rename(columns={'entity': 'gene_id'}), 'strategy').iloc[0]
        counts = card['counts']
        assert counts['correct'] == int(legacy.right) and counts['answered'] == int(legacy.answered)
        expected = card['extra']['accuracy_among_calls']
        assert (expected is None and not np.isfinite(legacy.rate)) or np.isclose(expected, legacy.rate)
        assert np.isclose(card['extra']['accuracy_all_hidden'], card['metrics']['accuracy'])
        cards.append(card)
        seen.update((scope.organism, str(entity)) for entity in rows.entity)
        reports.append({'organism': scope.organism, 'strategy': scope.strategy, 'target': scope.target,
            'settings': scope.settings, 'seed': scope.seed, 'partition': scope.partition,
            **counts, 'accuracy_all_hidden': card['extra']['accuracy_all_hidden'],
            'accuracy_among_calls': expected, 'legacy_rate_matches_call_precision': True,
            'standard_accuracy_matches_all_hidden': True, 'independent_truth': 'unresolved'})
    assert sum(c['counts']['eligible'] for c in cards) == len(ledger)
    frame = pd.DataFrame(reports)
    frame.to_csv(output / 'cohort_reconciliation.csv', index=False)
    (output / 'cards.json').write_text(json.dumps(cards, indent=2, allow_nan=False) + '\n')
    by = frame.groupby(['organism', 'partition'], observed=True).agg(cohorts=('eligible', 'size'), evaluation_rows=('eligible', 'sum')).reset_index()
    by.to_csv(output / 'scope_coverage.csv', index=False)
    differences = frame.dropna(subset=['accuracy_among_calls'])
    summary = {'legacy_rows': len(ledger), 'distinct_cohorts': len(cards), 'unique_organism_gene_addresses': len(seen),
        'cohorts_with_abstentions': int((frame.abstained > 0).sum()),
        'different_all_hidden_vs_among_calls': int((~np.isclose(differences.accuracy_all_hidden, differences.accuracy_among_calls)).sum()),
        'standard_metric_mismatches': 0, 'legacy_precision_mismatches': 0,
        'biological_admission': 'not established; truth/exclusion/nested-fit lineage unresolved',
        'uncertainty': 'unavailable without supplied biological groups; repeats not independent samples',
        'runtime_inference_changes': 'none'}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    inputs = [Path(__file__), source, ROOT / 'starplast/record_scorecards.py', ROOT / 'starplast/scorecard.py', ROOT / 'starplast/track_record.py', ROOT / 'starplast/capabilities.py']
    (output / 'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(), 'summary': summary,
        'input_sha256': {str(p): sha(p) for p in inputs},
        'outputs': {p.name: sha(p) for p in output.iterdir() if p.suffix in {'.json', '.csv'}}}, indent=2) + '\n')
    return summary


def main():
    """Execute and annotate the frozen all-record reconciliation census."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Reconcile individual held-out outcomes and scorecard denominators')
    nb.md('Frozen scope: every shipped track-record row, grouped by organism/strategy/target/settings/seed/hold-out mode/named set. The existing file lacks independent truth and nested-fit lineage: recomputing arithmetic does not admit its biological benchmarks or calibrate its support.')
    nb.code('from scripts.reconcile_record_scorecards import reconcile', f'summary = reconcile({str(args.out)!r})', 'summary')
    nb.md('The legacy rate measures correct / answered; the standard accuracy measures correct / all eligible hidden entities. Abstentions lower coverage and all-hidden accuracy, but are neither wrong calls nor missing records. Settings/seeds/protocols stay separate; unique gene addresses are not independent studies. Class false positives use the full cohort and component performance requires isolated tests.')
    nb.write(str(args.out / 'reconciliation.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
