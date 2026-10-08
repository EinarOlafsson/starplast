"""Pinned functional recovery views refuse altered lineage and preserve user navigation."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from starplast import functional_results as F, organisms as O


def _document():
    path=Path(F.__file__).with_name('data')/'functional_results.json'
    return json.loads(path.read_text())


def _write(tmp_path,document):
    path=tmp_path/'bundle.json'
    path.write_text(json.dumps(document,allow_nan=False))
    return path,hashlib.sha256(path.read_bytes()).hexdigest()


def _repair_source_hashes(document):
    entry=document['benchmarks'][0];manifest=entry['metadata']['source_manifest']
    for name,payload in entry['payloads'].items():
        encoded=(F._canonical(payload)+'\n').encode()
        manifest['files'][name]={'bytes':len(encoded),'sha256':hashlib.sha256(encoded).hexdigest()}
    manifest['key']=F._hash(manifest['spec'])
    manifest['identity']=F._hash({key:value for key,value in manifest.items() if key!='identity'})
    entry['metadata']['source_artifact_identity']=manifest['identity']
    entry['metadata']['summary']['artifact_identity']=manifest['identity']


def test_shipped_view_retains_exact_rows_cards_and_null_biological_accuracy():
    benchmarks,reason=F.shipped(O.TOXOPLASMA)
    assert len(benchmarks)==1 and not reason
    benchmark=benchmarks[0]
    assert benchmark.source_matches('ec_number') and not benchmark.source_matches('interpro_id')
    assert benchmark.metadata['target']=='ec_major_classes'
    assert len(benchmark.rows)==benchmark.card['counts']['eligible']
    assert benchmark.card['counts']['correct']==sum(row['prediction']==row['truth'] for row in benchmark.rows)
    assert [card['class'] for card in benchmark.major_class_cards]==list('1234567')
    assert all(card['biological_precision'] is None and card['biological_recall'] is None for card in benchmark.major_class_cards)
    assert all(row['calibrated_confidence'] is None for row in benchmark.rows)
    assert F.shipped(O.FALCIPARUM)[0]==[]
    source=Path(__file__).resolve().parents[1]/benchmark.metadata['source_directory']/'held_out'
    for name,payload in benchmark.payloads.items():assert payload==json.loads((source/name).read_text())


def test_bundle_checksum_is_external_and_changed_payloads_are_refused(tmp_path):
    document=_document();path,digest=_write(tmp_path,document)
    with pytest.raises(ValueError,match='checksum'):F.load(path,expected_sha256='0'*64)
    with pytest.raises(ValueError,match='checksum'):F.load(path,expected_sha256='')
    document['benchmarks'][0]['payloads']['rows.json'][0]['truth']='["7"]'
    path,digest=_write(tmp_path,document)
    with pytest.raises(ValueError,match='payload changed'):F.load(path,expected_sha256=digest)


@pytest.mark.parametrize('mutation,expected',[
    ('confidence','confidence'),('duplicate','population'),('missing_class','seven'),
    ('wrong_class_count','metrics'),('unknown_truth','Unknown'),('training_overlap','training'),
    ('wrong_organism','organism'),('biological','biological'),
    ('missing_profile','reference-profile'),('wrong_profile_count','profile-class'),
    ('missing_baseline','training-only'),('summary_mismatch','summary')])
def test_structural_semantic_refusals_survive_self_consistent_envelope(tmp_path,mutation,expected):
    document=deepcopy(_document());entry=document['benchmarks'][0];payload=entry['payloads']
    if mutation=='confidence':payload['rows.json'][0]['calibrated_confidence']=.99
    elif mutation=='duplicate':payload['rows.json'][1]['entity']=payload['rows.json'][0]['entity']
    elif mutation=='missing_class':payload['major_class_cards.json'].pop()
    elif mutation=='wrong_class_count':payload['major_class_cards.json'][0]['recall']=.999
    elif mutation=='unknown_truth':payload['rows.json'][0]['truth']=None
    elif mutation=='training_overlap':entry['metadata']['source_manifest']['spec']['fit_entities'].append(payload['rows.json'][0]['entity'])
    elif mutation=='wrong_organism':entry['metadata']['organism']=O.FALCIPARUM
    elif mutation=='biological':entry['metadata']['biological_admission']=True
    elif mutation=='missing_profile':payload['profile_class_cards.json'].pop()
    elif mutation=='wrong_profile_count':payload['profile_class_cards.json'][0]['class_metrics']['true_positive']+=1
    elif mutation=='missing_baseline':payload['baseline_cards.json'].pop('majority')
    elif mutation=='summary_mismatch':entry['metadata']['summary']['test']+=1
    _repair_source_hashes(document)
    path,digest=_write(tmp_path,document)
    with pytest.raises(ValueError,match=expected):F.load(path,expected_sha256=digest)


def test_functional_profile_validation_preserves_overlapping_classes_and_unknown():
    assert F.profile_classes('["2","3"]')==('2','3')
    assert F.profile_classes(None) is None
    for invalid in ('[]','["2","2"]','["3","2"]','["8"]'):
        with pytest.raises(ValueError):F.profile_classes(invalid)


def test_panel_functional_source_links_class_cards_abstentions_and_gene_signal(monkeypatch):
    from starplast import claims as C, strategies as S, track_record as T
    from starplast.discoveries_panel import DiscoveriesPanel
    benchmarks,_=F.shipped(O.TOXOPLASMA);benchmark=benchmarks[0]
    monkeypatch.setattr(C,'shipped',lambda *a:pd.DataFrame())
    monkeypatch.setattr(C,'recipes',lambda *a:pd.DataFrame())
    monkeypatch.setattr(T,'shipped',lambda *a:pd.DataFrame())
    panel=DiscoveriesPanel(O.TOXOPLASMA,context=S.Context.shipped(O.TOXOPLASMA))
    try:
        panel.label.setCurrentIndex(panel.label.findData('ec_number'))
        assert panel.annotation_functional.isEnabled() and not panel.annotation_record.isEnabled()
        panel.annotation_functional.click()
        assert panel.tabs.currentWidget()==panel.functional_page
        assert panel.functional_rows.rowCount()==len(benchmark.rows)
        assert 'Independent biological activity accuracy: not evaluated' in panel.functional_note.text()
        assert 'Calibrated confidence: unavailable' in panel.functional_card.toPlainText()
        opened=[];panel.gene_chosen.connect(opened.append)
        panel._functional_gene_clicked(0,0)
        assert opened==[benchmark.rows[0]['entity']]
        panel.functional_view.setCurrentIndex(panel.functional_view.findData('major:2'))
        expected=[r for r in benchmark.rows if '2' in F.profile_classes(r['truth']) or '2' in (F.profile_classes(r['prediction']) or ())]
        assert panel.functional_shown_rows==expected
        assert 'reference prevalence' in panel.functional_card.toPlainText().casefold()
        major=next(card for card in benchmark.major_class_cards if card['class']=='2')
        outcomes=[panel.functional_rows.item(i,4).text() for i in range(panel.functional_rows.rowCount())]
        assert outcomes.count('recovered reference class')==major['true_positive']
        assert outcomes.count('unexpected reference-class call')==major['false_positive']
        assert sum('class missed' in outcome for outcome in outcomes)==major['false_negative']
        panel.functional_view.setCurrentIndex(panel.functional_view.findData('baseline:majority'))
        assert panel.functional_rows.rowCount()==0, 'Control cards must not display native predictions as control outcomes'
        panel.functional_view.setCurrentIndex(panel.functional_view.findData('strategy:'))
        assert panel.functional_rows.rowCount()==len(benchmark.rows)
        assert panel.shown.empty
    finally:panel.close()


def test_unavailable_functional_bundle_remains_explicit_in_panel(monkeypatch):
    from starplast import claims as C, strategies as S, track_record as T
    from starplast.discoveries_panel import DiscoveriesPanel
    monkeypatch.setattr(F,'shipped',lambda organism:([],'Functional results unavailable: checksum mismatch'))
    monkeypatch.setattr(C,'shipped',lambda *a:pd.DataFrame());monkeypatch.setattr(C,'recipes',lambda *a:pd.DataFrame())
    monkeypatch.setattr(T,'shipped',lambda *a:pd.DataFrame())
    ctx=S.Context(pd.DataFrame({'gene_id':['TGME49_100001'],'ec_number':['2.7.11.1']}),graph={},organism=O.TOXOPLASMA)
    panel=DiscoveriesPanel(O.TOXOPLASMA,context=ctx)
    try:
        panel.label.setCurrentIndex(panel.label.findData('ec_number'))
        assert not panel.annotation_functional.isEnabled()
        assert 'checksum mismatch' in panel.functional_note.text()
        assert panel.functional_rows.rowCount()==0
    finally:panel.close()


def test_frozen_functional_test_opens_real_installed_gene_view():
    from starplast.app import Window
    window=Window()
    try:
        panel=window.discoveries
        panel.label.setCurrentIndex(panel.label.findData('ec_number'))
        assert panel.annotation_functional.isEnabled()
        panel.annotation_functional.click()
        assert panel.tabs.currentWidget()==panel.functional_page
        gene=panel.functional_shown_rows[0]['entity']
        panel._functional_gene_clicked(0,0)
        assert gene in window.detail.toHtml()
        assert panel.functional_choice.currentData().metadata['organism']==O.TOXOPLASMA
        assert panel.functional_choice.currentData().metadata['biological_admission'] is False
        panel.label.setCurrentIndex(panel.label.findData('ec_number'))
        available=[i for i,row in enumerate(panel.coverage_shown)
            if row['reference_recovery']['status']=='verified_reference_recovery']
        assert len(available)==1
        panel._coverage_selected(available[0],0)
        assert panel.coverage_open.isEnabled()
        panel.coverage_open.click()
        assert panel.tabs.currentWidget()==panel.functional_page
        assert panel.functional_card.view is not None
        assert dict(panel.functional_card.view.counts)==panel.functional_choice.currentData().card['counts']
    finally:window.close()


def test_selected_domain_provenance_keeps_source_description_and_unknown_membership(monkeypatch):
    from PyQt6 import QtCore, QtGui
    from starplast import claims as C, discovery_labels as D, strategies as S, track_record as T
    from starplast.discoveries_panel import DiscoveriesPanel

    class Lookup:
        sources={'InterPro':{'metadata_url':'https://example.org/domain-source'}}

        def lookup(self,value):
            return {'name':'Current kinase nomenclature','source':'InterPro',
                'status':'current_metadata' if value=='IPR000001' else 'missing_current_metadata',
                'ontology_release':'110.0','accession_version_status':'original_version_not_recorded',
                'metadata_license_url':'https://example.org/metadata-license',
                'metadata_source_url':'https://example.org/domain-source'}

    monkeypatch.setattr(D.FD,'shipped',lambda:Lookup())
    monkeypatch.setattr(F,'shipped',lambda organism:([],'No frozen functional recovery results are packaged.'))
    monkeypatch.setattr(C,'shipped',lambda *a:pd.DataFrame());monkeypatch.setattr(C,'recipes',lambda *a:pd.DataFrame())
    monkeypatch.setattr(T,'shipped',lambda *a:pd.DataFrame())
    ctx=S.Context(pd.DataFrame({'gene_id':['TGME49_100001','TGME49_100002'],
        'interpro_id':['IPR000001','IPR000002'],'interpro_desc':['Original source kinase',None]}),graph={},organism=O.TOXOPLASMA)
    panel=DiscoveriesPanel(O.TOXOPLASMA,context=ctx)
    try:
        panel.label.setCurrentIndex(panel.label.findData('interpro_id'))
        index=panel.annotation_class.findData('IPR000001')
        assert 'Original source kinase' in panel.annotation_class.itemText(index)
        assert 'Current kinase nomenclature' not in panel.annotation_class.itemText(index)
        tooltip=panel.annotation_class.itemData(index,QtCore.Qt.ItemDataRole.ToolTipRole)
        assert 'Original source description' in tooltip and 'Original source kinase' in tooltip
        assert 'https://example.org/metadata-license' in tooltip and 'https://example.org/domain-source' in tooltip
        panel.annotation_class.setCurrentIndex(index)
        assert 'Displayed name source: Original source annotation' in panel.annotation_note.text()
        assert 'release: 110.0' in panel.annotation_note.text()
        assert panel.shown_members.gene_id.tolist()==['TGME49_100001']
        document=QtGui.QTextDocument()
        document.setHtml(panel.member_table.item(0,2).toolTip())
        assert 'does not verify gene function' in ' '.join(document.toPlainText().split())
        panel.annotation_class.setCurrentIndex(panel.annotation_class.findData('IPR000002'))
        assert 'missing current metadata' in panel.annotation_note.text()
        assert 'does not establish domain absence' in panel.annotation_note.text()
        assert panel.shown_members.gene_id.tolist()==['TGME49_100002']
        panel.annotation_class.setCurrentIndex(0)
        assert 'Domain IPR000002' not in panel.annotation_note.text()
        assert panel.annotation_class.toolTip()==''
    finally:panel.close()
