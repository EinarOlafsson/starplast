"""Known-truth nomenclature cases and overlapping/unknown functional profiles."""
import pandas as pd
import pytest

from starplast import functional_ontology as F, organisms as O, strategies as S


def entries():
    return F.parse_enzyme('''ID   1.1.1.1
DE   First activity.
//
ID   2.7.1.1
DE   Kinase activity.
//
ID   3.6.3.1
DE   Transferred entry: 7.6.2.1.
//
ID   7.6.2.1
DE   Transport activity.
//
ID   4.1.1.1
DE   Transferred entry: 1.1.1.1 and 2.7.1.1.
//
ID   5.1.1.1
DE   Deleted entry.
//
ID   1.1.1.n1
DE   Preliminary activity.
//
''')


def test_only_unique_active_ec_replacements_become_canonical_memberships():
    ontology=entries()
    resolved=F.resolve_ec('3.6.3.1',ontology)
    assert resolved['canonical_ec']=='7.6.2.1'
    assert resolved['chain']==['3.6.3.1','7.6.2.1']
    assert resolved['status']=='unique_transfer'
    split=F.resolve_ec('4.1.1.1',ontology)
    assert split['status']=='ambiguous_transfer' and split['canonical_ec'] is None
    assert split['alternatives']==['1.1.1.1','2.7.1.1']
    for term,status in [('5.1.1.1','deleted'),('1.1.1.n1','preliminary'),('6.1.1.1','missing')]:
        result=F.resolve_ec(term,ontology)
        assert result['status']==status and result['canonical_ec'] is None


def test_replacement_cycles_and_missing_destinations_remain_unresolved():
    ontology=entries()
    ontology['1.1.1.1']=F.EnzymeEntry('1.1.1.1','Transfer','transferred',('2.7.1.1',))
    ontology['2.7.1.1']=F.EnzymeEntry('2.7.1.1','Transfer','transferred',('1.1.1.1',))
    assert F.resolve_ec('1.1.1.1',ontology)['status']=='cycle'
    ontology.pop('2.7.1.1')
    assert F.resolve_ec('1.1.1.1',ontology)['status']=='missing'


@pytest.mark.parametrize('text',['ID   1.1.1.1\nDE   No terminator.','ID   1.1.1.1\nDE   A.\nID   2.2.2.2\nDE   B.\n//',
    'ID   1.1.1.1\nDE   A.\n//\nID   1.1.1.1\nDE   B.\n//','<html>foreign payload</html>'])
def test_invalid_or_truncated_ontology_is_refused(text):
    with pytest.raises(ValueError):F.parse_enzyme(text)


def test_complete_profiles_keep_all_major_classes_and_unknowns_never_become_negative():
    values=['1.1.1.1 (first);2.7.1.1 (kinase)','3.6.3.1 (Transferred entry: 7.6.2.1)',
        None,'1.1.1.1 (first);4.1.1.1 (ambiguous)','5.1.1.1 (deleted)','malformed',
        '1.1.1.1 (first); bad assignment']
    nodes=pd.DataFrame({'gene_id':[f'TGME49_{100001+i}' for i in range(len(values))],'ec_number':values})
    ctx=S.Context(nodes,graph={},organism=O.TOXOPLASMA)
    original=ctx.nodes.copy(deep=True)
    profiles,resolutions=F.ec_profiles(ctx,'ec_number',entries())
    assert profiles.profile.iloc[:2].tolist()==['["1","2"]','["7"]']
    assert profiles.profile.iloc[2:].isna().all()
    assert profiles.eligible.tolist()==[True,True,False,False,False,False,False]
    assert profiles.status.iloc[2]=='unannotated'
    assert profiles.major_classes.iloc[3]==[], 'Ambiguous source cannot become partial truth'
    assert resolutions[resolutions.gene_id.eq(ctx.gene_ids[3])].status.tolist()==['active','ambiguous_transfer']
    pd.testing.assert_frame_equal(ctx.nodes,original,check_exact=True)


def test_major_class_card_counts_overlapping_truth_and_abstention_without_biological_precision():
    rows=pd.DataFrame({'truth':['["1","2"]','["2"]','["3"]','["1"]'],
        'prediction':['["2"]','["1","2"]',None,'["1"]']})
    cards={r['class']:r for r in F.member_class_cards(rows)}
    assert (cards['1']['true_positive'],cards['1']['false_positive'],cards['1']['false_negative'])==(1,1,1)
    assert cards['1']['precision']==.5 and cards['1']['recall']==.5
    assert cards['2']['precision']==1 and cards['2']['recall']==1
    assert cards['3']['recall']==0 and cards['3']['precision'] is None
    assert cards['7']['recall'] is None
    assert all(r['biological_precision'] is None and r['coverage']==.75 for r in cards.values())


@pytest.mark.parametrize('profile',['[]','["8"]','["2","1"]','["1","1"]','[1]','unknown'])
def test_noncanonical_profile_labels_are_refused(profile):
    with pytest.raises(ValueError):F.profile_members(profile)


@pytest.mark.parametrize('source',['1.1.1.1garbage','1.1.1.1 (open;2.7.1.1 (nested)',
    '1.1.1.1 (closed))','1.1.1.1 (closed) garbage (extra)','1.1.1.1 garbage'])
def test_valid_ec_prefix_cannot_make_malformed_source_a_complete_profile(source):
    assert not F.valid_ec_source(source)
    ctx=S.Context(pd.DataFrame({'gene_id':['TGME49_100001'],'ec_number':[source]}),graph={},organism=O.TOXOPLASMA)
    profiles,_=F.ec_profiles(ctx,'ec_number',entries())
    assert not profiles.eligible.iloc[0] and pd.isna(profiles.profile.iloc[0])
