"""Source-class membership preserves false calls, abstentions and full denominators."""
from copy import deepcopy
import json
import pytest
from starplast import functional_pair_class as C, functional_agreement as A, functional_results as F
from starplast.scorecard_view import export_scorecard,render_scorecard_html


def fixture():
    values=[('["2","3"]','["2"]','["2","3"]'),('["3"]','["2"]','["3"]'),
            ('["2"]',None,'["3"]'),('["3"]',None,None)]
    base=[]
    for i,(truth,left,right) in enumerate(values):
        base.append({'entity':str(i),'truth':truth,'group':str(i),'training_supported':i!=2,
            'left_prediction':left,'right_prediction':right})
    pairs=[[{**{k:r[k] for k in ('entity','truth','group','training_supported')},
             'prediction':r[side+'_prediction'],'abstained':r[side+'_prediction'] is None,
             'calibrated_confidence':None} for r in base] for side in ('left','right')]
    rows=A.paired_rows(*pairs)
    return {'rows':rows,**A._summary(rows),'methods':[{'strategy':'feature_knn'},{'strategy':'random_forest'}],
        'organism':'synthetic','target':'recorded profiles','namespace':'ec_major','protocol':'synthetic',
        'scope':{'truth_grade':'synthetic_control','target_identity':'fixture','context':'fixture',
                 'negative_semantics':'Recorded nonmembership only'},'overlap':{'independent_evidence':False},
        'group_uncertainty':{'status':'unavailable','reason':'Synthetic fixture only'},
        'biological_admission':False,'biological_accuracy':None,'calibrated_confidence':None,'ensemble_selection':False}


def test_major_mapping_differs_from_whole_profile_equality_and_keeps_false_calls():
    report=fixture();result=C.derive(report,'major:2')
    assert result['counts']['eligible']==4 and result['displayed_count']==3
    assert result['counts']['full_profile_unsupported']==1
    assert result['rows'][0]['outcome']=='conflict_right_correct'
    assert result['rows'][0]['class_outcome']=='agree_correct'
    assert result['rows'][1] in result['displayed_rows'] and not result['rows'][1]['recorded_member']
    assert result['rows'][2]['left_member'] is None
    assert result['methods'][0]['true_positive']==1 and result['methods'][0]['false_positive']==1
    assert result['methods'][0]['source_precision']==.5 and result['methods'][0]['source_recall']==.5
    assert result['methods'][1]['false_negative']==1
    assert result['counts']['both_abstain']==1


def test_profile_mapping_requires_complete_equality_and_empty_positives_are_unavailable():
    result=C.derive(fixture(),'profile:["2","3"]')
    assert result['rows'][0]['recorded_member'] and result['rows'][0]['left_member'] is False
    assert result['rows'][0]['right_member'] is True
    result=C.derive(fixture(),'major:7')
    assert result['counts']['known_positive_genes']==0 and result['joint_positive_precision'] is None
    assert all(m['source_precision'] is None and m['source_recall'] is None for m in result['methods'])


@pytest.mark.parametrize('address',['major:8','major:23','profile:["7"]','bad:2','major','profile:bad'])
def test_unknown_or_implicit_class_queries_are_refused(address):
    with pytest.raises(ValueError):C.derive(fixture(),address)


def test_native_pair_all_class_cards_and_individual_outcome_keep_original_counts():
    benchmarks,reason=F.shipped('Pf');assert not reason
    pair=tuple(next(b for b in benchmarks if b.metadata['strategy']==s) for s in ('feature_knn','random_forest'))
    report=A.compare(*pair)
    for address in C.addresses(pair):
        result,view=C.build(pair,report,address)
        assert result['counts']['eligible']==152 and len(result['rows'])==152
        assert result['counts']['full_profile_unsupported']==4
        assert result['displayed_count']<=152
        assert json.loads(export_scorecard(view))['snapshot']==result
        assert 'not biological absence' in render_scorecard_html(view,expanded=True)
    result,view=C.build(pair,report,'major:2')
    individual=C.gene_scorecard(report,result,result['displayed_rows'][0]['entity'])
    assert individual.kind=='outcome' and not individual.metrics
    assert json.loads(export_scorecard(individual))['snapshot']['query']=='major:2'
    with pytest.raises(ValueError,match='identities'):C.build(tuple(reversed(pair)),report,'major:2')
    broken=deepcopy(pair);broken[0].payloads['major_class_cards.json'][1]['true_positive']+=1
    with pytest.raises(ValueError,match='changed|counts differ'):C.build(broken,report,'major:2')
