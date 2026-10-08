"""Execute one frozen ortholog-transfer pilot on existing spatial predictions."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import artifacts as A, baselines as B, capabilities as C, organisms as O, record_scorecards as R, scorecard as SC, strategies as S, strategy_catalog as N, transfer_records as T  # noqa: E402
from starplast.ground_truth import read_registry  # noqa: E402
from starplast.query import Query  # noqa: E402
from starplast.splits import read_split, write_split  # noqa: E402


def _records(frame):
    return frame.astype(object).where(frame.notna(), None).to_dict(orient='records')


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def freeze(output):
    """Freeze native donor aggregation and receiver-training-only label mapping.

    The receiver population, groups and nested roles are exactly those from the
    first kNN pilot. Both spatial targets are predictions and stage/context
    equivalence and biological independence remain unresolved.
    """
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new transfer-pilot directory')
    registry = ROOT / 'results/ground_truth_registry_2026_10_07_v2/registry.json'
    entries, _ = read_registry(registry, relocate_snapshots_to=registry.parent)
    receiver = next(e for e in entries if e.organism == O.TOXOPLASMA and e.target == 'compartment' and e.task == SC.T_LABEL)
    donor = next(e for e in entries if e.organism == O.FALCIPARUM and e.target == 'lopit_pf_location' and e.task == SC.T_LABEL)
    split_path = ROOT / 'results/label_knn_pilot_2026_10_08_v5/split.json'
    split, exclusions = read_split(split_path)
    if split.benchmark_id != receiver.benchmark_id or exclusions.benchmark_id != receiver.benchmark_id:
        raise ValueError('Pilot receiver truth and split differ')
    paths = {e.organism: Path(O.nodes_path(e.organism)) for e in (receiver, donor)}
    tables = {organism: pd.read_parquet(path).set_index('gene_id', drop=False) for organism, path in paths.items()}
    for entry in (receiver, donor):
        universe = json.loads(Path(entry.universe_file.path).read_text())
        if tables[entry.organism].index.tolist() != universe:
            raise ValueError('Installed and registered entity orders differ')
        stored = pd.read_parquet(entry.truth_file.path).set_index('gene_id')[entry.target]
        pd.testing.assert_series_equal(tables[entry.organism][entry.target], stored, check_dtype=False)
    ids = tuple(row.entity for row in split.assignments)
    universe = tables[receiver.organism].index.tolist()
    mask = receiver.mask(universe)
    if list(mask.index[mask]) != list(ids):
        raise ValueError('Frozen receiver split omits/reorders eligible candidate genes')
    donor_nodes = tables[donor.organism][['gene_id', 'orthogroup', donor.target]].reset_index(drop=True)
    donor_context = S.Context(donor_nodes, graph={}, organism=donor.organism)
    # The projection context has no receiver target column or label values.
    receiver_nodes = tables[receiver.organism].loc[list(ids), ['gene_id', 'orthogroup']].reset_index(drop=True)
    context = S.Context(receiver_nodes, graph={}, other=donor_context, organism=receiver.organism)
    mapped, numeric = N._through_orthologs(context, donor.target)
    mapped.index = list(ids)
    projection = {'organism': receiver.organism, 'donor_organism': donor.organism, 'donor_column': donor.target,
        'receiver_entities': list(ids), 'receiver_orthogroups': receiver_nodes.orthogroup.astype(object).where(receiver_nodes.orthogroup.notna(), None).tolist(),
        'source_values': mapped.astype(object).where(mapped.notna(), None).tolist(),
        'donor_table_sha256': _sha(paths[donor.organism]), 'projection': 'native orthogroup median/mode; missing groups excluded'}
    projection = json.loads(json.dumps(projection, allow_nan=False))
    mapping_hash = hashlib.sha256(A.canonical_object(projection).encode()).hexdigest()
    source = T.TransferSource(donor.organism, donor.target, _sha(paths[donor.organism]), mapping_hash,
        donor.evidence_grade, 'Existing spatial-classifier outputs in different organisms/stages. '
        'Training receiver modes learn a statistical association; semantic equivalence and independent biology are unresolved.',
        receiver_truth_dependency='unknown')
    truth = pd.read_parquet(receiver.truth_file.path).set_index('gene_id')[receiver.target]
    labels = truth.loc[list(split.entities('train'))].astype(str)
    batch = T.ortholog_transfer(mapped, labels, split=split, exclusions=exclusions, source=source, numeric_source=numeric)
    # Outer truth becomes visible only after transfer fitting and prediction.
    rows = batch.rows.copy()
    rows['truth'] = truth.loc[rows.entity].astype(str).to_numpy()
    rows['correct'] = pd.array([p == t if isinstance(p, str) else None for p, t in zip(rows.prediction, rows.truth)], dtype='boolean')
    groups = {assignment.entity: assignment.group for assignment in split.assignments}
    rows['group'] = rows.entity.map(groups)
    gaps = tuple(sorted(set((*receiver.gaps, *donor.gaps, *batch.gaps,
        'Both targets are stored classifier predictions; recovery is surrogate agreement, not independent biological accuracy',
        'Source-stage semantic compatibility and shared study/annotation dependence require review'))))
    settings = {'source_organism': donor.organism, 'source_column': donor.target, 'numeric_source': numeric,
        'donor_aggregation': 'native orthogroup mode', 'receiver_mapping': 'training-only native mode'}
    scope = R.RecordScope(receiver.organism, 'ortholog_transfer', receiver.target, SC.T_LABEL,
        A.canonical_object(settings), split.seed, split.identity, 'outer_test', receiver.benchmark_id,
        receiver.evidence_grade, 'gene', receiver.negative_semantics, gaps)
    card = R.aggregate(rows, scope)
    classes = R.class_cards(rows, scope)
    baseline = B.label_baselines(labels, split, seed=split.seed)
    baseline_cards = {}
    for column in ('majority', 'prevalence_call'):
        baseline_rows = rows.drop(columns=['prediction', 'correct', 'abstained']).assign(prediction=baseline[column].to_numpy())
        baseline_cards[column] = R.aggregate(baseline_rows, scope)
    exclusion_hash = hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest()
    code = {'transfer_adapter': ROOT/'starplast/transfer_records.py', 'label_batch': ROOT/'starplast/label_records.py',
        'native_transfer': ROOT/'starplast/strategy_catalog.py', 'native_context': ROOT/'starplast/strategies.py',
        'native_text_mapping': ROOT/'starplast/embedding.py', 'scorecard_aggregation': ROOT/'starplast/record_scorecards.py',
        'standard_metrics': ROOT/'starplast/scorecard.py', 'baseline_recipes': ROOT/'starplast/baselines.py',
        'pilot_assembler': Path(__file__)}
    dependencies = [A.Dependency('code', name, _sha(path)) for name, path in code.items()]
    dependencies.extend(A.Dependency(kind, name, digest) for kind, name, digest in (
        ('table', 'receiver_gene_universe', _sha(paths[receiver.organism])), ('table', 'donor_gene_universe', _sha(paths[donor.organism])),
        ('truth', receiver.benchmark_id, receiver.truth_file.sha256), ('source', donor.benchmark_id, donor.truth_file.sha256),
        ('mapping', 'native_orthogroup_projection', mapping_hash), ('split', 'nested', split.identity),
        ('exclusions', 'training_source_closure', exclusion_hash)))
    code_version = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() + '+transfer-adapter-worktree'
    spec = A.ArtifactSpec(Query(receiver.organism, 'label', target=receiver.target).to_json(), 'ortholog_transfer',
        SC.T_LABEL, next(o for o in C.get('ortholog_transfer').outputs if o.kind == 'label_calls'), receiver.target, tuple(rows.entity), 'held_out',
        'outer:test:'+split.identity, tuple(dependencies), A.canonical_object(settings),
        A.canonical_object({'unit': 'gene', 'eligible_population': len(rows),
            'truth_grade': receiver.evidence_grade, 'donor_grade': donor.evidence_grade,
            'negative_semantics': receiver.negative_semantics,
            'context': {'receiver': receiver.contexts, 'donor': donor.contexts}, 'limitations': gaps,
            'benchmark_status': receiver.status, 'benchmark_admitted': False}), split.seed, code_version,
        fit_entities=split.entities('train'), fit_role='train', benchmark_id=receiver.benchmark_id,
        confidence_kind='none', evaluation_partition='test', gaps=gaps)
    payloads = {'rows.json': _records(rows), 'row_columns.json': rows.columns.tolist(),
        'model_state.json': batch.model_state, 'projection.json': projection,
        'card.json': card, 'class_cards.json': classes, 'baseline_cards.json': baseline_cards}
    payloads = json.loads(json.dumps(payloads, allow_nan=False))
    output.mkdir(parents=True)
    write_split(output/'split.json', split, exclusions=exclusions)
    identity = A.write_artifact(output/'held_out', spec, payloads, split=split)
    assert A.read_artifact(output/'held_out', expected=spec, split=split).payloads == payloads
    summary = {'strategy': 'ortholog_transfer', 'receiver': receiver.benchmark_id, 'donor': donor.benchmark_id,
        'receiver_grade': receiver.evidence_grade, 'donor_grade': donor.evidence_grade, 'benchmark_admitted': False,
        'eligible_candidate_genes': len(ids), 'retained_test_rows': len(rows), 'mapped_receiver_genes': int(mapped.notna().sum()),
        'mapped_training_genes': len(batch.model_state['training_mapped_entities']), 'mapped_test_genes': int(rows.donor_mapped.sum()),
        'counts': card['counts'], 'metrics': card['metrics'], 'baseline_metrics': baseline_cards['majority']['metrics'],
        'native_class_scores': 'unavailable', 'native_support': 'unavailable', 'artifact_identity': identity,
        'mapping_identity': mapping_hash, 'interpretation': 'stored-prediction surrogate recovery only', 'runtime_changes': 'none'}
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    inputs = [registry, split_path, *paths.values(), *(Path(e.truth_file.path) for e in (receiver, donor)),
        *(Path(e.universe_file.path) for e in (receiver, donor)), *code.values(),
        ROOT/'starplast/splits.py', ROOT/'starplast/artifacts.py', ROOT/'starplast/capabilities.py', ROOT/'starplast/ground_truth.py']
    (output/'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256': {str(path): _sha(path) for path in inputs}, 'summary': summary,
        'outputs': {str(path.relative_to(output)): _sha(path) for path in output.rglob('*.json')}}, indent=2)+'\n')
    archive = output/'code'
    archive.mkdir()
    for path in inputs:
        if path.suffix == '.py':
            (archive/path.name).write_bytes(path.read_bytes())
    return summary


def main():
    """Execute the scoped native transfer pilot and preserve its annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Frozen ortholog-transfer records with receiver-training-only mapping')
    nb.md('Fixed pilot: existing Plasmodium schizont spatial-prediction labels projected through native orthogroups onto the same 3,827 eligible Toxoplasma compartment genes and seed-17 nested split as the first kNN pilot. The projection context has no receiver target column. Receiver mapping fits only on training labels. Stage/context equivalence and independent biological truth remain unresolved; no data is added or replaced.')
    nb.code('from scripts.freeze_transfer_pilot import freeze', f'summary=freeze({str(args.out)!r})', 'summary')
    nb.md('All 560 test genes are retained, including missing projections and unsupported donor categories. Native transfer provides no confidence or class-score matrix. Counts, per-class outcomes and baselines use the identical receiver test cohort. Numerical source values, installed strategies and calibration are unchanged. Full outer coverage and the other missing label-capable adapters remain open under 64.10.')
    nb.write(str(args.out/'pilot.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
