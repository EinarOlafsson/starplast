"""Host source retention is evidence, with explicit non-admission and no accuracy."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from starplast import host_foundation_scorecard as H,organisms as O
from starplast.scorecard_view import export_scorecard,render_scorecard_detail,render_scorecard_html


def _summary():
    return json.loads((Path(__file__).resolve().parents[1]/'results/human_gene_space_foundation_2026_10_08/summary.json').read_text())


def test_reviewed_source_counts_and_every_gap_survive_shared_evidence_card():
    summary=_summary();view=H.build(summary)
    assert view.kind=='evidence' and view.task is None and view.metrics==()
    assert dict(view.scope)['organism']==O.HUMAN
    assert dict(view.counts)=={field:summary[field] for field in H.COUNT_FIELDS}
    assert view.source.sha256==summary['source_sha256']
    assert view.source.grade=='unresolved' and 'Unavailable' in view.source.version
    assert 'unavailable' in view.status.lower()
    exported=json.loads(export_scorecard(view))
    evidence=exported['snapshot']['card']['evidence']
    assert evidence['biological_accuracy'] is None and evidence['calibrated_confidence'] is None
    assert evidence['admission']['benchmark_admitted'] is False
    assert evidence['admission']['installed_space_registered'] is False
    assert evidence['admission']['distributable_pack_built'] is False
    assert evidence['admission']['redistribution']=='unresolved'
    assert evidence['protein_measurements_projected']==0
    assert evidence['reviewed_context_and_admission_gaps']==summary['gaps']
    for gap in summary['gaps']:assert gap in render_scorecard_detail(view,'gaps')
    html=render_scorecard_html(view,expanded=True)
    assert 'Cultured adult fibroblasts are not validated HFF measurements' in html
    assert 'Result quality / coverage' not in html and 'scorecard:metric/accuracy' not in html
    assert 'Cultured adult fibroblasts are not validated HFF measurements' in render_scorecard_html(view)


def test_card_is_detached_from_mutable_input_and_missing_provenance_stays_unknown():
    summary=_summary();view=H.build(summary)
    summary['gaps'].clear();summary['canonical_genes']=0
    assert dict(view.counts)['canonical_genes']==58988
    assert 'Cultured adult fibroblasts' in render_scorecard_detail(view,'gaps')
    freshness=render_scorecard_detail(view,'freshness')
    assert 'null' in freshness and 'unresolved' in freshness


@pytest.mark.parametrize('mutation', ['wrong_species','wrong_unit','record_gap','mapping_gap','negative_count',
    'boolean_count','benchmark_admission','space_admission','pack_admission','redistribution',
    'projected_measurements','source_loss','missing_gaps','invalid_sha'])
def test_unreviewed_admission_scope_or_counts_are_refused(mutation):
    summary=deepcopy(_summary())
    if mutation=='wrong_species':summary['organism']=O.MOUSE
    elif mutation=='wrong_unit':summary['unit']='protein'
    elif mutation=='record_gap':summary['source_records']+=1
    elif mutation=='mapping_gap':summary['unmapped_genes']-=1
    elif mutation=='negative_count':summary['qualified_PAR_Y_records']=-1
    elif mutation=='boolean_count':summary['multiple_protein_genes']=False
    elif mutation=='benchmark_admission':summary['benchmark_admitted']=True
    elif mutation=='space_admission':summary['installed_space_registered']=True
    elif mutation=='pack_admission':summary['distributable_pack_built']=True
    elif mutation=='redistribution':summary['redistribution']='approved'
    elif mutation=='projected_measurements':summary['protein_measurements_projected']=10
    elif mutation=='source_loss':summary['source_records_lost']=1
    elif mutation=='missing_gaps':summary['gaps']=[]
    elif mutation=='invalid_sha':summary['source_sha256']='unknown'
    with pytest.raises(ValueError):H.build(summary)


def test_summary_file_identity_is_retained_and_bad_pin_refused(tmp_path):
    path=tmp_path/'foundation.json';path.write_text(json.dumps(_summary()))
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    view=H.load(path,expected_sha256=digest)
    assert digest in render_scorecard_detail(view,'freshness')
    with pytest.raises(ValueError,match='checksum'):H.load(path,expected_sha256='0'*64)
    with pytest.raises(ValueError,match='Summary provenance'):H.build(_summary(),summary_sha256='unknown')
