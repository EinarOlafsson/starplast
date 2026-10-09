"""Paired recorded-class comparisons with explicit closed-source semantics."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace

from . import functional_agreement as A, functional_results as F
from .scorecard_view import DetailView

LIMITS = ('Membership is derived from recorded complete source profiles. Recorded '
          'nonmembership is not biological absence. Abstention is unknown, not a '
          'negative call. All metrics use the full frozen known-source cohort. '
          'Biological accuracy, calibrated confidence and class-specific group '
          'intervals are unavailable; pooled profile intervals do not apply. '
          'Shared-source methods are not independent votes. Training support '
          'flags refer to the original full profile, not individual class support.')


def outcome_title(row):
    labels={'both_abstain':'Both abstained','left_only_correct':'kNN matches; forest abstained',
        'left_only_wrong':'kNN differs; forest abstained','right_only_correct':'Forest matches; kNN abstained',
        'right_only_wrong':'Forest differs; kNN abstained','conflict_left_correct':'kNN matches; forest differs',
        'conflict_right_correct':'Forest matches; kNN differs','conflict_both_wrong':'Both differ from source'}
    if row['class_outcome']=='agree_correct':return 'Both match recorded membership' if row['recorded_member'] else 'Both match recorded nonmembership'
    if row['class_outcome']=='agree_wrong':return 'Both miss recorded membership' if row['recorded_member'] else 'Both call unrecorded membership'
    return labels[row['class_outcome']]


def addresses(pair):
    """Return distinct recorded class addresses, including classes seen by one method."""
    major = sorted({c['class'] for b in pair for c in b.major_class_cards})
    profiles = sorted({c['class'] for b in pair for c in b.profile_class_cards})
    return ['major:' + value for value in major] + ['profile:' + value for value in profiles]


def derive(report, address):
    """Count source-class agreements across every gene before filtering display rows."""
    A.build_scorecard(report)
    if not isinstance(address, str) or ':' not in address:
        raise ValueError('Explicit major or complete-profile address required')
    kind, value = address.split(':', 1)
    if kind == 'major':
        if value not in '1234567' or len(value) != 1 or report['namespace'] != 'ec_major':
            raise ValueError('Unsupported recorded major-class mapping')
        contains = lambda profile: value in F.profile_classes(profile, report['namespace'])
        mapping = 'Presence of the selected major class in the complete recorded profile'
    elif kind == 'profile':
        F.profile_classes(value, report['namespace'])
        if value not in {r['truth'] for r in report['rows']} | {
                r[key] for r in report['rows'] for key in ('left_prediction','right_prediction')}:
            raise ValueError('Profile absent from the frozen comparison records')
        contains = lambda profile: profile == value
        mapping = 'Exact complete-profile equality; individual members are not interchangeable'
    else:
        raise ValueError('Explicit major or complete-profile address required')
    rows = []; binary = [[], []]
    for row in report['rows']:
        actual = contains(row['truth'])
        decisions = [None if row[key] is None else contains(row[key])
                     for key in ('left_prediction','right_prediction')]
        for side, prediction in enumerate(decisions):
            binary[side].append({**{k: row[k] for k in ('entity','group','training_supported')},
                'truth': 'member' if actual else 'nonmember',
                'prediction': None if prediction is None else 'member' if prediction else 'nonmember',
                'abstained': prediction is None, 'calibrated_confidence': None})
        rows.append({**deepcopy(row), 'recorded_member': actual,
                     'left_member': decisions[0], 'right_member': decisions[1]})
    projected = A.paired_rows(*binary)
    for row, projected_row in zip(rows, projected):
        row['class_outcome'] = projected_row['outcome']
    summary = A._summary(projected)
    counts = summary['counts']
    counts['full_profile_unsupported']=counts.pop('unsupported')
    counts.update(known_positive_genes=sum(r['recorded_member'] for r in rows),
        agreed_presence=sum(r['left_member'] is True and r['right_member'] is True for r in rows),
        agreed_nonmembership=sum(r['left_member'] is False and r['right_member'] is False for r in rows),
        jointly_called_true_positive=sum(r['recorded_member'] and r['left_member'] is True and r['right_member'] is True for r in rows))
    methods = []
    for side, name in enumerate(('left_member','right_member')):
        tp=sum(r['recorded_member'] and r[name] is True for r in rows)
        fp=sum(not r['recorded_member'] and r[name] is True for r in rows)
        fn=sum(r['recorded_member'] and r[name] is not True for r in rows)
        methods.append({'strategy':report['methods'][side]['strategy'], 'eligible':len(rows),
            'positive_calls':tp+fp,'true_positive':tp,'false_positive':fp,'false_negative':fn,
            'source_precision':tp/(tp+fp) if tp+fp else None,
            'source_recall':tp/(tp+fn) if tp+fn else None})
    displayed=[r for r in rows if r['recorded_member'] or r['left_member'] is True or r['right_member'] is True]
    return {'address':address,'kind':kind,'value':value,'mapping':mapping,'limits':LIMITS,
        'counts':counts,'rates':summary['rates'],'methods':methods,
        'joint_positive_precision':counts['jointly_called_true_positive']/counts['agreed_presence'] if counts['agreed_presence'] else None,
        'rows':rows,'displayed_rows':displayed,'displayed_count':len(displayed),
        'parent_methods':deepcopy(report['methods']),'parent_scope':deepcopy(report['scope']),
        'parent_target':report['target'],'parent_protocol':report['protocol'],
        'parent_overlap':deepcopy(report['overlap']),
        'biological_admission':False,'calibrated_confidence':None,'class_uncertainty':'unavailable'}


def build(pair, report, address):
    """Present counted class metrics and unchanged native per-method class records."""
    if len(pair)!=2 or any(b.metadata['source_artifact_identity']!=m['artifact_identity']
                          or b.metadata['strategy']!=m['strategy']
                          for b,m in zip(pair,report['methods'])):
        raise ValueError('Class method identities differ from the paired report')
    if address not in addresses(pair):
        raise ValueError('Class is outside the recorded method cards')
    result=derive(report,address)
    native=[]
    for benchmark,method in zip(pair,result['methods']):
        cards=benchmark.major_class_cards if result['kind']=='major' else benchmark.profile_class_cards
        card=next((c for c in cards if c['class']==result['value']),None)
        if card is not None:
            observed=card if result['kind']=='major' else card['class_metrics']
            if any(observed[key]!=method[key] for key in ('true_positive','false_positive','false_negative')):
                raise ValueError('Native class counts differ from paired projection')
        native.append({'strategy':method['strategy'],'card':deepcopy(card),
            'status':'recorded' if card is not None else 'unavailable for this method'})
    result['native_method_class_cards']=native
    view=A.build_scorecard(report)
    evidence={'Full known-source cohort':len(result['rows']),
        'Displayed reference members or positive calls':result['displayed_count'],
        'Recorded members':result['counts']['known_positive_genes'],
        'Both predict membership':result['counts']['agreed_presence'],
        'Both predict nonmembership':result['counts']['agreed_nonmembership'],
        'Source precision among jointly predicted members':result['joint_positive_precision']}
    for method in result['methods']:
        def rate(n,d):return f'{n} / {d} ({n/d:.1%})' if d else 'Unavailable: no denominator'
        name=method['strategy'].replace('_',' ')
        evidence[name+' source precision']=rate(method['true_positive'],method['positive_calls'])
        evidence[name+' source recall']=rate(method['true_positive'],method['true_positive']+method['false_negative'])
    details=(DetailView('evidence','Recorded membership comparison',F._canonical(evidence)),
        DetailView('sizes','Full cohort and displayed rows',F._canonical({'eligible':len(result['rows']),'displayed':result['displayed_count']})),
        DetailView('uncertainty','Class-specific uncertainty',F._canonical({'status':'unavailable','reason':'No class-specific group interval computed; pooled profile intervals do not apply'})),
        DetailView('calibration','Calibration',F._canonical({'status':'unavailable','reason':LIMITS})),
        DetailView('baseline','Class baseline',F._canonical({'status':'unavailable','reason':'No matched class-specific consensus control'})),
        DetailView('freshness','Original frozen method identities',F._canonical(report['methods'])),
        DetailView('rows','All full-cohort class outcomes',F._canonical(result['rows'])),
        DetailView('native_class_cards','Original per-method class cards',F._canonical(native)),
        DetailView('gaps','Mapping and negative semantics',F._canonical({'mapping':result['mapping'],'limits':LIMITS})))
    view=replace(view,title='Paired recorded '+('class: ' if result['kind']=='major' else 'complete profile: ')+F.class_title(result['value'],report['namespace']),
        scope=tuple((key,'recorded membership '+address) if key=='target' else (key,value)
                    for key,value in view.scope),
        source=replace(view.source,lineage=LIMITS,negative_semantics=LIMITS),
        counts=tuple(result['counts'].items()),metrics=(),details=details,snapshot_json=F._canonical(result))
    return result,view


def gene_scorecard(report, result, entity):
    """Present one observed class outcome, retaining its original profile calls."""
    rows=[r for r in result['rows'] if r['entity']==entity]
    if len(rows)!=1:raise ValueError('Gene outside the frozen class comparison')
    view=A.build_gene_scorecard(report,entity)
    snapshot={'query':result['address'],'mapping':result['mapping'],'outcome':rows[0],
              'limits':LIMITS,'parent_methods':result['parent_methods'],'calibrated_confidence':None}
    details=tuple(replace(d,text=F._canonical(rows[0])) if d.key=='outcome' else d for d in view.details)
    return replace(view,title='Observed recorded-class outcome: '+entity,
                   scope=tuple((key,'recorded membership '+result['address']) if key=='target' else (key,value)
                               for key,value in view.scope),
                   source=replace(view.source,lineage=LIMITS),details=details,snapshot_json=F._canonical(snapshot))
