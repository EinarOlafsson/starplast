"""Freeze declared strategy/technique capabilities and synthetic adapter routing."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import capabilities as C, organisms as O, strategies as S  # noqa: E402
from starplast.query import EntityRef, Query  # noqa: E402


def freeze(output):
    """Write one immutable schema pilot; no strategy predictions are executed."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new capability snapshot directory')
    output.mkdir(parents=True)
    caps = C.catalog()
    declarations = [c.to_dict() for c in caps]
    rows, native = [], []
    for organism in (O.TOXOPLASMA, O.FALCIPARUM):
        ids = tuple(f'TGME49_{i:06d}' if organism == O.TOXOPLASMA else f'PF3D7_{i:07d}' for i in range(100001, 100007))
        entities = tuple(EntityRef(organism, 'gene', i) for i in ids)
        other = O.FALCIPARUM if organism == O.TOXOPLASMA else O.TOXOPLASMA
        state = C.EvidenceState(organism, ids,
            columns=(('location', 'categorical'), ('growth', 'numeric'), ('baseline', 'numeric')),
            layers=(('contacts', 'physical'), ('folds', 'structural'), ('papers', 'literature')),
            permitted_feature_count=10, orthogroups=True, attention=True,
            other_organism=other, other_columns=('other_growth',), other_orthogroups=True,
            class_values=(('location', ('rhoptry', 'nucleus')),), admitted_gene_space=True)
        for cap in caps:
            settings = {p: 'other_growth' if kind == 'other' else 'growth' if kind == 'numeric' else 'location' for p, kind in cap.column_parameters}
            if cap.strategy == 'condition_shift':
                settings['baseline'] = 'baseline'
            if 'layer' in cap.required:
                settings['layer'] = 'folds' if cap.strategy == 'structural_homology' else 'papers' if cap.strategy == 'attention_correction' else 'contacts'
            preferred = 'gene_set' if cap.minimum_seeds or cap.strategy == 'neighbour_space' else 'pair' if cap.outputs[0].kind == 'pair_ranking' else 'trait' if cap.column_parameters and cap.column_parameters[0][1] == 'numeric' else 'label' if cap.column_parameters else 'gene'
            for kind in ('gene', 'protein', 'gene_set', 'label', 'class', 'trait', 'pair'):
                query = Query(organism, kind,
                    entities=entities if kind == 'gene_set' else entities[:2] if kind == 'pair' else (EntityRef(organism, 'protein', 'P00001'),) if kind == 'protein' else entities[:1] if kind == 'gene' else (),
                    target='growth' if kind == 'trait' else 'location' if kind in {'label', 'class'} else '',
                    values=('rhoptry',) if kind == 'class' else ())
                plan = C.resolve(query, cap.strategy, state, settings=settings)
                row = {'organism': organism, 'strategy': cap.strategy, 'query_kind': kind,
                       'applicable': plan.applicable, 'reason': plan.reason, 'detail': plan.detail,
                       'outputs': [asdict(o) for o in plan.outputs], 'validation_gaps': plan.validation_gaps}
                rows.append(row)
                if kind == preferred:
                    assert plan.applicable, row
                    native.append(asdict(plan))
    for organism in (O.HUMAN, O.MOUSE):
        query = Query(organism, 'gene', (EntityRef(organism, 'gene', 'synthetic_host_gene'),))
        state = C.EvidenceState(organism, ('synthetic_host_gene',), admitted_gene_space=False)
        for cap in caps:
            plan = C.resolve(query, cap.strategy, state)
            assert not plan.applicable and plan.reason == 'organism_adapter_unavailable'
            rows.append({'organism': organism, 'strategy': cap.strategy, 'query_kind': 'gene',
                         'applicable': False, 'reason': plan.reason, 'detail': plan.detail, 'outputs': [], 'validation_gaps': []})
    parameters = {s.key: s.parameters().to_dict(orient='records') for s in S.catalog()}
    for name, payload in (('capabilities.json', declarations), ('routing.json', rows),
                          ('native_plans.json', native), ('parameters.json', parameters)):
        (output / name).write_text(json.dumps(payload, indent=2) + '\n')
    inputs = [Path(__file__), ROOT / 'starplast/capabilities.py', ROOT / 'starplast/query.py',
              ROOT / 'starplast/strategies.py', ROOT / 'starplast/strategy_catalog.py',
              ROOT / 'starplast/strategy_graph.py', ROOT / 'starplast/strategy_learning.py',
              ROOT / 'starplast/techniques.py']
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    summary = {'strategies': len(caps), 'techniques': len({t.technique for c in caps for t in c.techniques}),
               'query_strategy_pairs': len(rows), 'native_plans': len(native),
               'scope': 'synthetic adapter contracts; no biological predictions or fitted models',
               'component_validation': 'declared_not_measured', 'numerical_changes': 'none'}
    (output / 'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
        'schema_version': C.SCHEMA_VERSION, 'summary': summary,
        'input_sha256': {str(p): sha(p) for p in inputs},
        'outputs': {p.name: sha(p) for p in output.glob('*.json')}}, indent=2) + '\n')
    return summary


def main():
    """Execute and record the bounded catalogue/adapter declaration census."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Strategy capability and technique-validation contracts')
    nb.md('Frozen scope: all 39 current catalogue keys, all 40 underlying techniques, seven typed query kinds on two parasite fixtures and explicit unavailable host adapters. This is software routing evidence, not biological accuracy or inference capacity. No inference runner is invoked.')
    nb.code('from scripts.freeze_capability_catalog import freeze', f'summary = freeze({str(args.out)!r})', 'summary')
    nb.md('Score meanings remain distinct. Class recovery and evidence ablation are not gene calls; relationship support ranks pairs; imputation emits transformed percentiles; condition shifts are rank-scale residuals. Independent truth admission and executed component tests remain open under subsequent benchmark cards.')
    nb.write(str(args.out / 'catalogue.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
