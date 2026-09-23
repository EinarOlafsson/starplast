"""Observation identity, missingness and review gates survive portable storage."""
from dataclasses import replace
import numpy as np
import pandas as pd
import pytest

from starplast.evidence import (Observation, LiteratureAssertion, write_observations,
                               read_observations, from_gene_table, propose_assertions,
                               save_assertions, load_assertions)


def observation(**kwargs):
    return Observation('g1', 'abundance', kwargs.pop('value', 0), 'study1', 'example', **kwargs)


def test_roundtrip_preserves_types_missingness_and_contradictions(tmp_path):
    records = [observation(value=False), observation(value=0, uncertainty=1),
               observation(value=3.5, source_version='v2'),
               observation(value=None, missing_state='not_assayed'),
               observation(value={'site': 12, 'detected': True}, entity_type='residue')]
    path = tmp_path/'observations.parquet'
    digest = write_observations(records + [records[0]], path)
    actual = read_observations(path)
    assert len(digest) == 64 and len(actual) == 5
    assert {r.record_id for r in actual} == {r.record_id for r in records}
    assert {type(r.value) for r in actual} == {bool, int, float, type(None), dict}


def test_content_tampering_is_detected(tmp_path):
    path = tmp_path/'observations.parquet'
    write_observations([observation()], path)
    frame = pd.read_parquet(path); frame.loc[0,'value_json']='9'
    frame.to_parquet(path,index=False)
    with pytest.raises(ValueError, match='identity'): read_observations(path)


@pytest.mark.parametrize('kwargs', [dict(value=None), dict(value=0,missing_state='not_assayed'),
                                  dict(uncertainty=-1), dict(value=np.nan)])
def test_invalid_or_ambiguous_observations_are_refused(kwargs):
    with pytest.raises(ValueError): observation(**kwargs)


def test_legacy_summary_does_not_invent_assay_status():
    nodes = pd.DataFrame({'gene_id':['g1','g2'], 'custom':[0.,np.nan]})
    records = list(from_gene_table(nodes, 'Tg', 'snapshot1', include_missing=True))
    assert records[0].value == 0 and records[0].evidence_status == 'unclassified'
    assert records[1].value is None and records[1].missing_state == 'unknown'
    assert records[0].context['resolution'] == 'existing_gene_summary'


def test_literature_claims_need_review_and_retain_negation(tmp_path):
    class Index:
        def find(self, text): return [('g1', 'symbol')]
    claims = propose_assertions([{'text':'ABC is not localized to the nucleus.',
                                 'source_id':'paper', 'source_url':'https://example.org'}], Index())
    assert len(claims)==1 and claims[0].negated
    with pytest.raises(ValueError, match='review'): claims[0].observation('example')
    with pytest.raises(ValueError, match='reviewer'): claims[0].review(True, '')
    accepted = claims[0].review(True,'Einar')
    record = accepted.observation('example')
    assert record.value['negated'] and record.evidence_status=='computed'
    assert record.context['reviewer']=='Einar'
    path = tmp_path/'claims.jsonl'
    save_assertions([accepted,claims[0].review(False,'Einar')],path)
    assert [a.review_status for a in load_assertions(path)]==['accepted','rejected']
