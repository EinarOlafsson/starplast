"""Execute bounded shared-scorecard presentation checks as a recorded notebook.

The fixtures exercise software contracts for all six task types, explicit
controls, prediction sets, numeric intervals, unavailable metadata and individual
outcomes. No estimator is fitted, no biological source is acquired, and no
independent biological admission or calibrated probability is inferred. Every
output comes from the executed notebook namespace and a new output directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook


def main():
    """Execute and retain one immutable presentation audit with bounded fixtures.

    The notebook creates the input rows and computes standard record cards, then
    independently renders two hosts from identical supplied cards. This function
    writes only verified notebook outputs and byte receipts; an existing output
    directory is refused rather than overwritten, and unchanged source identities
    are checked before the receipt and completion summary are written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/scorecard_task_views_2026_10_08')
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError('Audit outputs are immutable; choose a new output directory')
    nb = ExecutedNotebook('Shared scorecard task views: synthetic software acceptance')
    nb.ns['root_path'] = str(ROOT)
    nb.md('Scope: ordinary presentation software only.',
        'These tiny hand-authored records exercise all six existing task contracts and explicit '
        'controls. They are not biological benchmarks. No downloads, estimator fitting, Qt host '
        'creation, activity truth or organism benchmark acquisition occurs.',
        'The two hosts are identical pure HTML/model/export consumers. Desktop navigation and '
        'scientific source, split, calibration and task-admission acceptance remain separate gates.')
    nb.code("""
from dataclasses import asdict
from html import escape
import hashlib, json, platform, re
from pathlib import Path
import numpy as np
import pandas as pd
from starplast import organisms as O, record_scorecards as R, scorecard as SC
from starplast import scorecard_view as V

def canonical(value):
    return json.dumps(V._plain(value), sort_keys=True, ensure_ascii=False,
                      separators=(',', ':'), allow_nan=False)

def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()

code_files = ['scripts/audit_scorecard_task_views.py', 'scripts/notebook_runner.py',
    'starplast/scorecard_view.py', 'starplast/scorecard.py', 'starplast/record_scorecards.py',
    'starplast/strategies.py', 'starplast/capabilities.py', 'starplast/provenance.py',
    'starplast/organisms.py']
code_receipts = {name: hashlib.sha256((Path(root_path) / name).read_bytes()).hexdigest()
                 for name in code_files}
checks = []

def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append({'check': name, 'passed': True})

def scope(strategy, task, unit='gene'):
    return R.RecordScope(O.TOXOPLASMA, strategy, 'synthetic software target', task,
        'fixed authored fixture; no fitting', 3, 'synthetic software fixture protocol',
        'software_test', 'SCORECARD-SOFTWARE-ONLY', 'synthetic_control', unit,
        'Artificial fixture truth only; no biological negative inference',
        ('Synthetic software behavior does not admit biological performance',))

labels = pd.DataFrame({'entity': ['g1','g2','g3','g4'], 'truth': ['A','A','B','B'],
                      'prediction': ['A',None,'A','B']})
binary = pd.DataFrame({'entity': labels.entity, 'positive': [True,False,True,False],
                      'score': [4.,3.,2.,1.], 'returned': [True,True,False,False]})
numeric = pd.DataFrame({'entity': labels.entity, 'truth': [1.,2.,3.,4.],
                       'prediction': [1.2,2.1,3.5,np.nan], 'baseline_prediction': [2.,2.,2.,2.]})
cluster = labels.assign(cluster=[0,0,1,-1])
findings = pd.DataFrame({'entity': ['finding1','finding2','finding3','finding4'],
                         'replicated': [True,False,True,False]})
fixtures = {
    'labels': (labels, scope('feature_knn', SC.T_LABEL), {}),
    'ranking': (binary, scope('positive_unlabeled', SC.T_RANK), {'verified_negatives': True}),
    'retrieval': (binary, scope('geneset_hunt', SC.T_SET), {'verified_negatives': True}),
    'clusters': (cluster, scope('holdout_search', SC.T_CLUSTER), {'chosen_clusters': {'A':0,'B':1}}),
    'values': (numeric, scope('trait_regression', SC.T_VALUES),
               {'quantity_unit': 'synthetic assay units', 'baseline_name': 'authored constant-two control'}),
    'replication': (findings, scope('blind_battery', SC.T_REPL, 'finding'), {'null_rates':[.25,.5]}),
}
controls = {
    'labels': labels.assign(prediction='A'),
    'ranking': binary.assign(score=[1.,2.,3.,4.]),
    'retrieval': binary.assign(returned=True),
    'clusters': cluster.assign(cluster=0),
    'values': numeric.assign(prediction=2.),
    'replication': findings.assign(replicated=False),
}
cards = {name: R.aggregate(rows, address, parameters=params)
         for name, (rows, address, params) in fixtures.items()}
control_cards = {name: R.aggregate(controls[name], address, parameters=params)
                 for name, (_, address, params) in fixtures.items()}
check('all six distinct task contracts', {card['scope']['task'] for card in cards.values()} == set(SC.TASKS))
pd.DataFrame([{'fixture': name, 'task': card['scope']['task'], **card['counts']}
              for name, card in cards.items()])
""")
    nb.md('Set and interval variants preserve supplied parameters.',
        'Prediction sets and bounds are hand-authored fixture outputs. Nominal coverage is '
        'declared before this software check; the nominal set field is explicitly attached '
        'presentation metadata because the row aggregator does not emit it. It is not '
        'estimated confidence or a claim of '
        'biological exchangeability. Empty, unbounded and unavailable cases remain explicit.')
    nb.code("""
set_rows = labels.copy()
set_rows['prediction_set'] = [['A'],['A','B'],[],['B']]
set_rows['prediction'] = ['A',None,None,'B']
set_params = {'prediction_classes':['A','B'], 'nominal_set_coverage':.9}
set_card = R.aggregate(set_rows, scope('conformal_calls', SC.T_LABEL), parameters=set_params)
set_card['metrics']['promised_coverage'] = set_params['nominal_set_coverage']
set_card['metric_status']['promised_coverage'] = 'available'
set_control = set_rows.copy()
set_control['prediction_set'] = [['A','B'] for _ in range(4)]
set_control['prediction'] = None
set_control_card = R.aggregate(set_control, scope('conformal_calls', SC.T_LABEL), parameters=set_params)

interval_rows = numeric.copy()
interval_rows['lower'] = [0.,1.,None,None]
interval_rows['upper'] = [2.,2.5,None,None]
interval_rows['interval_status'] = ['finite','finite','unbounded','unavailable']
interval_params = {'quantity_unit':'synthetic assay units',
    'baseline_name':'authored constant-two control','nominal_interval_coverage':.9}
interval_card = R.aggregate(interval_rows, scope('trait_regression', SC.T_VALUES), parameters=interval_params)
check('sets include coverage and size together', set_card['metrics']['set_coverage'] == .75
      and set_card['metrics']['mean_set_size'] == 1.)
check('all-class control is broad without efficiency', set_control_card['metrics']['set_coverage'] == 1.
      and set_control_card['metrics']['mean_set_size'] == 2. and set_control_card['metrics']['set_efficiency'] == 0.)
check('unbounded interval does not invent finite width', interval_card['metrics']['interval_coverage'] == .75
      and interval_card['metrics']['mean_interval_width'] is None)
cards.update(prediction_sets=set_card, numeric_intervals=interval_card)
fixtures.update(prediction_sets=(set_rows, scope('conformal_calls',SC.T_LABEL), set_params),
                numeric_intervals=(interval_rows, scope('trait_regression',SC.T_VALUES), interval_params))
control_cards['prediction_sets'] = set_control_card
control_cards['numeric_intervals'] = control_cards['values']

positive_only = binary.iloc[:3].copy()
positive_only['positive'] = [True,None,True]
for name, strategy, task, params in [
    ('positive_only_ranking','positive_unlabeled',SC.T_RANK,{'verified_negatives':False,'depth':2}),
    ('positive_only_retrieval','geneset_hunt',SC.T_SET,{'verified_negatives':False})]:
    address = scope(strategy,task)
    cards[name] = R.aggregate(positive_only,address,parameters=params)
    fixtures[name] = (positive_only,address,params)
    control_cards[name] = None
check('unknown negative precision remains unavailable', cards['positive_only_retrieval']['metrics']['precision'] is None)
check('unknown negative ranking rates remain unavailable', cards['positive_only_ranking']['metrics']['auroc'] is None
      and cards['positive_only_ranking']['metrics']['auprc'] is None)
{'set_metrics':set_card['metrics'], 'interval_metrics':interval_card['metrics']}
""")
    nb.md('Two pure consumers reuse identical supplied values and definitions.',
        'The audit compares independent model construction, HTML and exact JSON export; '
        'each displayed metric destination resolves to its task authority. Input row digests '
        'are fixture-source identities, distinct from record cohort/card/code identities. '
        'No QObject or desktop window is instantiated.')
    nb.code("""
inputs = {}
exports = {}
html_outputs = {}
reports = []
for name, card in cards.items():
    frame, address, params = fixtures[name]
    records = V._plain(frame.astype(object).where(frame.notna(),None).to_dict('records'))
    inputs[name] = {'scope':asdict(address), 'rows':records, 'parameters':params,
                    'controls':control_cards[name]}
    source_identity = digest(inputs[name])
    supplied_details = {
        'source':{'name':'Authored synthetic '+name+' fixture','grade':'synthetic_control',
            'sha256':source_identity,'context':'Software presentation only',
            'lineage':'Hand-authored records; no fit or biological source',
            'negative_semantics':address.negative_semantics},
        'split':{'protocol':address.protocol,'partition':address.partition,
                 'fit_population':None,'reason':'No fitting in this software fixture'},
        'controls':{'card':control_cards[name], 'interpretation':'Explicit authored software control; not biological chance'},
        'rows':records,'freshness':{'fixture_source_identity':source_identity},
    }
    first = V.build_scorecard_view(card,title='Identical pure host',details=supplied_details)
    second = V.build_scorecard_view(card,title='Identical pure host',details=supplied_details)
    first_export, second_export = V.export_scorecard(first), V.export_scorecard(second)
    check(name+': identical metric views', first.metrics == second.metrics)
    check(name+': identical HTML', V.render_scorecard_html(first) == V.render_scorecard_html(second))
    check(name+': identical exact exports', first_export == second_export)
    exported = json.loads(first_export)
    check(name+': exact native JSON card snapshot', exported['snapshot']['card'] == V._plain(card))
    metric_by_key = {metric.key:metric for metric in first.metrics}
    for key, value in card['metrics'].items():
        check(name+': exact metric '+key, metric_by_key[key].value == value)
        authority = SC.metric_definition(address.task,key)
        check(name+': task authority '+key, authority is not None and metric_by_key[key].definition == authority.definition)
        destination = 'scorecard:metric/'+key
        check(name+': valid metric route '+key, V.validate_scorecard_link(destination) == destination)
        check(name+': escaped definition popup '+key,
              escape(authority.definition) in V.render_scorecard_detail(first,key))
    if address.task == SC.T_LABEL:
        check(name+': all-eligible denominator', str(card['counts']['eligible']) in metric_by_key['accuracy'].denominator)
        check(name+': calls-only denominator', str(card['counts']['answered']) in metric_by_key['precision_of_calls'].denominator)
    else:
        check(name+': no label accuracy borrowing', 'precision_of_calls' not in metric_by_key)
    exports[name] = exported
    html_outputs[name] = {'compact':V.render_scorecard_html(first),
                          'expanded':V.render_scorecard_html(first,expanded=True)}
    reports.append({'fixture':name,'task':address.task,'input_identity':source_identity,
        'scope_identity':card['scope_identity'],'records_identity':card['records_identity'],
        'card_identity':digest(card),'displayed_metric_count':len(first.metrics),
        'null_metrics':[key for key,value in card['metrics'].items() if value is None]})
set_view = V.build_scorecard_view(set_card)
interval_view = V.build_scorecard_view(interval_card)
set_nominal = next(metric for metric in set_view.metrics if metric.key=='promised_coverage')
interval_nominal = next(metric for metric in interval_view.metrics if metric.key=='promised_coverage')
check('nominal set coverage has its own definition', 'prediction-set' in set_nominal.label.lower()
      and 'label test' in set_nominal.definition)
check('nominal interval coverage has its own definition', 'interval' in interval_nominal.label.lower()
      and set_nominal.definition != interval_nominal.definition)
check('full binary64 selective accuracy preserved', exports['labels']['snapshot']['card']['metrics']['precision_of_calls'] == 2/3)
pd.DataFrame(reports)
""")
    nb.md('Missing metadata, evidence quality and individual outcomes remain visible.',
        'A source audit is not an inference rate. One right/wrong outcome does not become '
        'a per-gene correctness probability. Unsafe navigation and HTML injection are rejected '
        'or escaped before a desktop host would receive them.')
    nb.code("""
missing = V.build_scorecard_view(cards['labels'])
check('missing source explicitly unavailable', missing.source.sha256 == 'Unavailable / not supplied')
check('missing source version explicitly unavailable', missing.source.version == 'Unavailable / not supplied')
check('missing uncertainty reason retained', 'Biological groups unresolved' in V.render_scorecard_detail(missing,'uncertainty'))
check('missing freshness explicitly unavailable', 'Unavailable / not supplied' in V.render_scorecard_detail(missing,'freshness'))
check('missing baseline explicitly unavailable', 'Unavailable / not supplied' in V.render_scorecard_detail(missing,'baseline'))
untested = V.build_scorecard_view({'task':SC.T_LABEL,'metrics':{}})
check('untested state visible', 'Untested / unavailable' in V.render_scorecard_html(untested))
evidence = V.build_scorecard_view({'evidence':{'source_rows':4,'biological_accuracy':None,
    'biological_admission':False}},kind='evidence',status='Synthetic source handling; performance unavailable')
outcome = V.build_scorecard_view({'entity':'g2','outcome':{'truth':'A','prediction':None,'abstained':True},
    'per_gene_accuracy_probability':None},kind='outcome')
check('evidence does not show performance metrics', not evidence.metrics and evidence.kind=='evidence')
check('individual outcome has no accuracy metric', not outcome.metrics and outcome.kind=='outcome')
check('individual probability remains null', json.loads(V.export_scorecard(outcome))['snapshot']['card']['per_gene_accuracy_probability'] is None)
check('individual probability limitation visible', 'no per-gene accuracy probability' in V.render_scorecard_html(outcome))
injected = V.build_scorecard_view(cards['labels'],title='<script>fixture()</script>',
    details={'rows':'<img src=x onerror=fixture()>'})
check('HTML title escaped', '<script>' not in V.render_scorecard_html(injected)
      and '&lt;script&gt;' in V.render_scorecard_html(injected))
check('HTML metadata escaped', '<img ' not in V.render_scorecard_detail(injected,'rows'))
for bad in ['javascript:fixture()', 'file:///tmp/fixture', 'https://user:secret@example.org',
            'https://example.org/%0aevil', 'scorecard:metric/accuracy/extra']:
    try:
        V.validate_scorecard_link(bad)
    except ValueError:
        check('refused unsafe route '+bad, True)
    else:
        check('refused unsafe route '+bad, False)
exports['evidence'] = json.loads(V.export_scorecard(evidence))
exports['outcome'] = json.loads(V.export_scorecard(outcome))
exports['untested'] = json.loads(V.export_scorecard(untested))
html_outputs['evidence'] = {'compact':V.render_scorecard_html(evidence),'expanded':V.render_scorecard_html(evidence,expanded=True)}
html_outputs['outcome'] = {'compact':V.render_scorecard_html(outcome),'expanded':V.render_scorecard_html(outcome,expanded=True)}
for name, expected in code_receipts.items():
    check('unchanged code '+name, hashlib.sha256((Path(root_path)/name).read_bytes()).hexdigest()==expected)
summary = {'schema_version':1,'status':'verified_software_presentation_only','task_types':list(SC.TASKS),
    'aggregate_cards':len(cards),'check_count':len(checks),'all_checks_passed':all(item['passed'] for item in checks),
    'biological_benchmarks_admitted':0,'calibrated_confidence_estimates_created':0,
    'desktop_hosts_tested':False,'new_sources_acquired':0,'models_fitted':0,
    'runtime':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},
    'limitations':['Synthetic fixtures verify renderer/export behavior only',
        'Two identical hosts are pure consumers; actual Qt host integration is a root acceptance gate',
        'Source/split/calibration/task admission remains unavailable; whole action 64.17 stays open',
        'Authored controls are not independent biological null experiments'],
    'cards':reports,'code_sha256':code_receipts}
summary
""")
    ns = nb.ns
    out.mkdir(parents=True)
    payloads = {'inputs.json': ns['inputs'], 'cards.json': ns['cards'],
                'controls.json': ns['control_cards'], 'checks.json': ns['checks'],
                'summary.json': ns['summary'], 'exports.json': ns['exports']}
    for name, value in payloads.items():
        (out / name).write_text(json.dumps(value, sort_keys=True, ensure_ascii=False,
            indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for name, html in ns['html_outputs'].items():
        (out / (name + '.html')).write_text('<!doctype html><meta charset="utf-8">' + html['expanded'], encoding='utf-8')
    nb.write(str(out / 'audit.ipynb'))
    readme = (
        '# Shared scorecard task views — software acceptance\n\n'
        f"**{ns['summary']['check_count']} checks passed** across all six task contracts, ten synthetic aggregate cards, "
        'prediction-set and numeric-interval metadata, explicit controls, positive-only truth, unavailable states, '
        'evidence quality and individual outcomes.\n\n'
        'The executed `audit.ipynb` computes every recorded input/card/check in its shared namespace. '
        '`exports.json` retains the exact native cards; both independently constructed pure consumers agree '
        'on all float/null values, HTML and JSON. Metric routes resolve their task-specific definitions, '
        'including distinct nominal set and interval coverage. Code/card/record/input identities and '
        'runtime versions are retained in `summary.json` and byte receipts in `SHA256SUMS`.\n\n'
        'These are authored software fixtures, with no source acquisition, estimator fitting, biological '
        'benchmark admission or calibrated confidence. No Qt application ran. Actual desktop host checks '
        'and scientific source/split/calibration/task-admission gates remain separate. **64.17 remains open.**\n\n'
        'Run from the repository with the existing environment and a 400 MB external memory cap:\n\n'
        '```bash\n'
        'systemd-run --user --scope -p MemoryMax=400M /home/carruthers/anaconda3/envs/starplast/bin/python '
        'scripts/audit_scorecard_task_views.py\n'
        '```\n\n'
        'The script refuses an existing output directory. Use `--out` with a new destination to reproduce; '
        'notebook cell timings may differ, while canonical cards and displayed values must agree.\n'
    )
    (out / 'README.md').write_text(readme, encoding='utf-8')
    files = sorted(path for path in out.iterdir() if path.is_file())
    (out / 'SHA256SUMS').write_text(''.join(hashlib.sha256(path.read_bytes()).hexdigest()
        + '  ' + path.name + '\n' for path in files), encoding='utf-8')
    print(json.dumps({'out': str(out), 'checks': ns['summary']['check_count'],
                      'status': ns['summary']['status']}, sort_keys=True))


if __name__ == '__main__':
    main()
