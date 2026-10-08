"""Functional class membership, missing inference coverage and real navigation."""
import numpy as np
import pandas as pd
import pytest

from starplast import claims as C, discovery_labels as D, organisms as O, strategies as S, track_record as T


def context():
    return S.Context(pd.DataFrame({'gene_id':[f'TGME49_{100001+i}' for i in range(4)],
        'product':['A','B','C','D'],'compartment':['nucleus','cytosol',None,None],
        'function_label':['A','A','B',None],'constant_label':['known']*4,
        'interpro_id':['IPR000001;IPR000002','IPR000002;IPR000002',None,'malformed'],
        'interpro_desc':['Kinase;Binding','Binding;Binding',None,'unpaired'],
        'pfam_id':['PF00069;PF00069',None,'PF00001.2',None],
        'ec_number':['2.7.11.1 (kinase); 3.1.1.- (hydrolase)',None,'2.7.11.1 (kinase)',None],
        'has_domain':[1,1,0,0],'has_ec':[1.,0.,1.,np.nan],
        'numerical_measurement':[1.,2.,3.,4.]}),graph={},organism=O.TOXOPLASMA)


def test_multi_valued_function_membership_deduplicates_terms_and_retains_unknown():
    members = D.annotation_members(context(),'interpro_id')
    assert len(members)==3
    assert members.groupby('value').gene_id.nunique().to_dict()=={'IPR000001':1,'IPR000002':2}
    assert members.description.tolist()==['Kinase','Binding','Binding']
    assert members.organism.eq(O.TOXOPLASMA).all()
    enzyme = D.annotation_members(context(),'ec_number')
    assert enzyme.groupby('value').gene_id.nunique().to_dict()=={'2.7.11.1':2,'3.1.1.-':1}
    assert set(enzyme.gene_id)=={'TGME49_100001','TGME49_100003'}
    pfam = D.annotation_members(context(),'pfam_id')
    assert pfam.value.tolist()==['PF00069','PF00001.2']
    assert pfam.description.eq('').all(), 'InterPro descriptions must never name Pfam IDs'


def test_ec_replacement_mentions_are_not_assigned_classes_and_descriptions_stay_paired():
    ctx = context()
    ctx.nodes.loc[0,'ec_number'] = '3.6.3.1 (Transferred entry: 7.6.2.1; 7.6.2.2 is related);2.7.11.1 (kinase)'
    rows = D.annotation_members(ctx,'ec_number')
    own = rows[rows.gene_id.eq(ctx.gene_ids[0])]
    assert own.value.tolist()==['3.6.3.1','2.7.11.1']
    assert 'kinase' not in own.description.iloc[0]
    inventory,_ = D.catalogue(ctx,pd.DataFrame(),pd.DataFrame())
    assert 'replacement' in inventory.set_index('target').loc['ec_number','annotation_warning']


def test_domain_names_enrich_only_missing_descriptions_without_changing_membership():
    from starplast import functional_domains as FD
    ctx = context()
    original = D.annotation_members(ctx,'pfam_id')
    named = D.annotation_members(ctx,'pfam_id',domain_lookup=FD.shipped())
    pd.testing.assert_frame_equal(named[original.columns].drop(columns='description'),
        original.drop(columns='description'),check_exact=True)
    kinase = named[named.value.eq('PF00069')].iloc[0]
    assert 'kinase' in kinase.description.lower()
    assert kinase.description_source=='Pfam current nomenclature'
    assert kinase.source_description==''
    assert kinase.metadata_source_url.startswith('https://ftp.ebi.ac.uk/')
    classes = D.class_summary(named,'pfam_id').set_index('value')
    assert classes.loc['PF00069','ontology_status']=='current_metadata'
    paired = D.annotation_members(ctx,'interpro_id',domain_lookup=FD.shipped())
    assert paired.description.tolist()==['Kinase','Binding','Binding']
    assert paired.description_source.eq('Original source annotation').all()
    assert paired.source_description.tolist()==['Kinase','Binding','Binding']


def test_missing_or_corrupt_nomenclature_keeps_annotation_browser_available(monkeypatch):
    def unavailable():raise ValueError('snapshot checksum mismatch')
    monkeypatch.setattr(D.FD,'shipped',unavailable)
    inventory,members = D.catalogue(context(),pd.DataFrame(),pd.DataFrame())
    assert members[members.target.eq('pfam_id')].value.tolist()==['PF00069','PF00001.2']
    assert 'unavailable' in inventory.set_index('target').loc['pfam_id','annotation_warning']


def test_installed_falciparum_missing_domain_descriptions_gain_searchable_names():
    ctx = S.Context.shipped(O.FALCIPARUM)
    inventory,members = D.catalogue(ctx,pd.DataFrame(),pd.DataFrame())
    original = D.annotation_members(ctx,'pfam_ids')
    named = members[members.target.eq('pfam_ids')].reset_index(drop=True)
    pd.testing.assert_frame_equal(named[original.columns].drop(columns='description'),
        original.drop(columns='description'),check_exact=True)
    assert named.description.str.contains('kinase',case=False,na=False).any()
    assert 'not verified gene function' in inventory.set_index('target').loc['pfam_ids','annotation_warning']


def test_unannotated_gene_table_retains_an_empty_catalogue_schema():
    ctx = S.Context(pd.DataFrame({'gene_id':['TGME49_100001']}),graph={},organism=O.TOXOPLASMA)
    inventory,members = D.catalogue(ctx,pd.DataFrame(),pd.DataFrame())
    assert inventory.empty and members.empty
    assert 'target' in inventory and 'target' in members


def test_all_annotation_recipe_and_claim_variables_remain_visible_without_calibration():
    recipes = pd.DataFrame({'target':['not_in_nodes','function_label'],'proven':[False,True]})
    claims = pd.DataFrame({'target':['claim_only'],'status':['untested']})
    catalogue,members = D.catalogue(context(),claims,recipes)
    assert {'interpro_id','ec_number','has_domain','has_ec','compartment','function_label','constant_label','not_in_nodes','claim_only'}<=set(catalogue.target)
    assert not {'product','gene_id','numerical_measurement'} & set(catalogue.target)
    assert catalogue.family.iloc[0]=='Function'
    functional = catalogue.set_index('target').loc['interpro_id']
    assert functional.annotated_genes==2 and functional.unannotated_genes==2
    assert functional.source_annotated_genes==3 and functional.unparsed_annotation_genes==1
    assert functional.claims==0 and functional.evaluation_status=='Inference not evaluated'
    assert catalogue.set_index('target').loc['not_in_nodes','evaluation_status']=='Recipe calibration incomplete'
    assert members[members.target.eq('interpro_id')].gene_id.nunique()==2
    classes = D.class_summary(members,'interpro_id')
    assert classes.precision.isna().all() and classes.recall.isna().all()


@pytest.mark.parametrize('organism',[O.TOXOPLASMA,O.FALCIPARUM])
def test_installed_function_catalogue_and_all_existing_claims_recipes_are_accessible(organism):
    ctx = S.Context.shipped(organism)
    claims,recipes = C.shipped(organism),C.recipes(organism)
    catalogue,members = D.catalogue(ctx,claims,recipes)
    assert set(claims.target.astype(str))|set(recipes.target.astype(str))<=set(catalogue.target)
    assert {'ec_number','interpro_id' if organism==O.TOXOPLASMA else 'interpro_ids'}<=set(catalogue.target)
    assert members.organism.eq(organism).all()
    assert set(members.gene_id)<=set(ctx.gene_ids)
    assert not members.duplicated(['organism','target','value','gene_id']).any()
    assert catalogue.claims.sum()==len(claims)


def test_function_search_class_members_and_gene_navigation_are_distinct_from_claims(monkeypatch):
    from PyQt6 import QtWidgets
    from starplast.discoveries_panel import DiscoveriesPanel
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    monkeypatch.setattr(C,'shipped',lambda *a:pd.DataFrame())
    monkeypatch.setattr(C,'recipes',lambda *a:pd.DataFrame())
    monkeypatch.setattr(T,'shipped',lambda *a:pd.DataFrame())
    panel = DiscoveriesPanel(O.TOXOPLASMA,context=context())
    try:
        assert panel.tabs.currentWidget()==panel.annotations_page
        panel.annotation_search.setText('kinase')
        assert set(panel.browsed_labels.target)=={'interpro_id','ec_number','pfam_id'}
        row = int(panel.browsed_labels.index[panel.browsed_labels.target.eq('interpro_id')][0])
        panel._select_annotation_label(row,0)
        assert panel.label.currentData()=='interpro_id'
        assert panel.annotation_class.findData('IPR000001')>=0
        panel.annotation_class.setCurrentIndex(panel.annotation_class.findData('IPR000001'))
        assert panel.shown_members.gene_id.tolist()==['TGME49_100001']
        opened = [];panel.gene_chosen.connect(opened.append)
        panel._member_clicked(0,0)
        assert opened==['TGME49_100001'] and panel.shown.empty
        assert 'not evaluated' in panel.annotation_note.text().lower()
        assert not panel.annotation_record.isEnabled()
        panel.annotation_search.setText('nothing matches this search')
        assert panel.browsed_labels.empty
    finally:panel.close()


def test_installed_class_opens_gene_and_existing_scorecard_in_real_window():
    from PyQt6 import QtWidgets
    from starplast.app import Window
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    window = Window()
    try:
        panel = window.discoveries
        panel.label.setCurrentIndex(panel.label.findData('interpro_id'))
        panel.annotation_class.setCurrentIndex(1)
        gene = panel.shown_members.gene_id.iloc[0]
        panel._member_clicked(0,0)
        assert gene in window.detail.toHtml()
        assert not panel.annotation_record.isEnabled()
        panel.label.setCurrentIndex(panel.label.findData('compartment'))
        panel.annotation_class.setCurrentIndex(1)
        assert panel.annotation_record.isEnabled()
        panel.annotation_record.click()
        html = window.detail.toHtml()
        assert 'starplast://target/compartment' in html and 'rate [95%]' in window.detail.toPlainText()
        # Verifier column names contain spaces and must survive table rendering.
        if len(panel.shown):
            row = panel.shown.iloc[0]
            said = [str(row[c]) for c in row.index if str(c).endswith(' verdict') and pd.notna(row[c]) and str(row[c])!='silent']
            expected = f"{said.count('agrees')} agree, {said.count('disagrees')} disagree" if said else '—'
            assert panel.table.item(0,6).text()==expected
    finally:window.close()
