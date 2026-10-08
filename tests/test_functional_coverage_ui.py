"""Coverage navigation cannot attach frozen results to another source context."""
import pandas as pd

from starplast import claims as C, functional_results as F, organisms as O, strategies as S, track_record as T
from starplast.discoveries_panel import DiscoveriesPanel


def _empty_claims(monkeypatch):
    monkeypatch.setattr(C,'shipped',lambda *args:pd.DataFrame())
    monkeypatch.setattr(C,'recipes',lambda *args:pd.DataFrame())
    monkeypatch.setattr(T,'shipped',lambda *args:pd.DataFrame())


def test_custom_context_keeps_annotations_and_refuses_borrowed_recovery(monkeypatch):
    _empty_claims(monkeypatch)
    ctx=S.Context(pd.DataFrame({'gene_id':['g1','g2'],'ec_number':['2.7.11.1',None]}),
        graph={},organism=O.TOXOPLASMA)
    panel=DiscoveriesPanel(O.TOXOPLASMA,context=ctx)
    try:
        panel.label.setCurrentIndex(panel.label.findData('ec_number'))
        panel.annotation_coverage.click()
        assert panel.tabs.currentWidget()==panel.coverage_page
        assert panel.coverage_shown and 'does not match this context' in panel.coverage_note.text()
        assert all(row['reference_recovery']['status']=='unavailable' for row in panel.coverage_shown)
        assert not panel.functional_benchmarks and not panel.annotation_functional.isEnabled()
        assert panel.functional_rows.rowCount()==0
        panel._coverage_selected(0,0)
        assert not panel.coverage_open.isEnabled()
        assert 'Missing annotation' in panel.coverage_card.toPlainText()
        assert panel.shown_members.gene_id.tolist()==['g1']
        assert panel.host_card.view.kind=='evidence' and panel.host_card.view.task is None
        assert dict(panel.host_card.view.scope)['organism']==O.HUMAN
        assert dict(panel.host_card.view.counts)['canonical_genes']==58988
        assert 'not validated HFF' in panel.host_card.toPlainText()
        assert not panel.host_card.view.metrics
    finally:panel.close()


def test_equal_size_altered_source_cannot_bind_archived_accuracy():
    original=S.Context.shipped(O.TOXOPLASMA)
    nodes=original.nodes.copy(deep=True)
    nodes.loc[nodes.index[0],'gene_id']='altered_fixture_gene'
    ctx=S.Context(nodes,graph={},organism=O.TOXOPLASMA)
    benchmark=F.shipped(O.TOXOPLASMA)[0][0]
    import pytest
    with pytest.raises(ValueError,match='values or entities differ'):
        F.require_context(benchmark,ctx)
    assert F.require_context(benchmark,original)


def test_missing_organism_result_never_inherits_other_organism_accuracy(monkeypatch):
    _empty_claims(monkeypatch)
    monkeypatch.setattr(F,'shipped',lambda organism:([],'No matching frozen results'))
    ctx=S.Context(pd.DataFrame({'gene_id':['g1'],'pfam_ids':['PF00001']}),graph={},organism=O.FALCIPARUM)
    panel=DiscoveriesPanel(O.FALCIPARUM,context=ctx)
    try:
        panel.label.setCurrentIndex(panel.label.findData('pfam_ids'))
        assert panel.coverage_shown
        assert all(row['organism']==O.FALCIPARUM for row in panel.coverage_shown)
        panel._coverage_selected(0,0)
        assert not panel.coverage_open.isEnabled()
        assert 'requires multivalued target adapter' in panel.coverage_card.toPlainText()
    finally:panel.close()
