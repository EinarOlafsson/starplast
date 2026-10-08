"""Sentinel leakage checks for opt-in source closure; no fitting or scientific claims."""
from dataclasses import replace

import pandas as pd
import pytest

from starplast import functional_exclusions as F, organisms as O
from starplast.splits import make_split, read_split, write_split


def fixture():
    ids=[f'fixture_{i}' for i in range(32)]
    split=make_split(ids,ids,organism=O.TOXOPLASMA,benchmark_id='functional_software_fixture')
    train=list(split.entities('train'))
    nodes=pd.DataFrame({'gene_id':train,'pfam_id':['PF00001']*len(train),
        'recorded_profile':['["PF00001"]']*len(train),
        'n_publications':[i%3 for i in range(len(train))],
        'n_fulltext_total':[i%4 for i in range(len(train))],
        'pfam_annotation_digest':[float(i%2) for i in range(len(train))],
        'raw_measurement':[float(i*i%11) for i in range(len(train))]})
    return nodes,split


def closure(nodes,split,**kwargs):
    return F.make_exclusions(nodes,'recorded_profile',source_targets=('pfam_id',),
        benchmark_id=split.benchmark_id,split=split,**kwargs)


def test_functional_family_and_attention_sentinels_are_refused():
    nodes,split=fixture();excluded=closure(nodes,split)
    assert excluded.fit_entities==split.entities('train')
    for column in ('pfam_id','recorded_profile','pfam_annotation_digest','n_publications','n_fulltext_total'):
        with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(columns=[column])
    for layer in ('domain','orthogroup','comention','comention_ft','unwritten_interaction','structural_hole'):
        with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(layers=[layer])
    excluded.guard_inputs(columns=['raw_measurement'],layers=['independent_measured_edges'])


@pytest.mark.parametrize('column', [
    'citation', 'citations', 'n_citations', 'citation_rate', 'paper', 'papers',
    'n_papers', 'paper_count', 'publication', 'publications', 'n_publication',
    'publication_year', 'fulltext', 'fulltexts', 'n_fulltexts', 'full_text_count',
    'n_full_texts', 'CITATIONS_total',
])
def test_attention_aliases_and_their_declared_derivatives_are_withheld(column):
    nodes, split = fixture()
    nodes[column] = nodes.raw_measurement
    nodes['raw_paperweight_measurement'] = nodes.raw_measurement
    lineage = (F.DerivedInput('column', 'opaque_attention_vector', (('column', column),)),)
    excluded = closure(nodes, split, derived_inputs=lineage)
    for field in (column, 'opaque_attention_vector'):
        with pytest.raises(ValueError, match='contamination'):
            excluded.guard_inputs(columns=[field])
    excluded.guard_inputs(columns=['raw_measurement', 'raw_paperweight_measurement'])


def test_declared_graph_embedding_and_second_generation_sources_close_transitively():
    nodes,split=fixture()
    lineage=(F.DerivedInput('layer','second_operator',(('column','encoded_embedding'),)),
        F.DerivedInput('column','encoded_embedding',(('layer','domain'),)),
        F.DerivedInput('column','scaled_embedding',(('layer','second_operator'),)),
        F.DerivedInput('column','independent_transform',(('column','raw_measurement'),)))
    excluded=closure(nodes,split,derived_inputs=lineage)
    for column in ('encoded_embedding','scaled_embedding'):
        with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(columns=[column])
    with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(layers=['second_operator'])
    excluded.guard_inputs(columns=['independent_transform'])


def test_caller_lineage_cannot_remove_known_native_derived_layer_sources():
    nodes,split=fixture()
    declared=(F.DerivedInput('layer','unwritten_interaction',(('column','raw_measurement'),)),)
    excluded=closure(nodes,split,derived_inputs=declared)
    # A caller cannot override the known native source lineage to unban a layer.
    with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(layers=['unwritten_interaction'])


@pytest.mark.parametrize('role',['tune','calibration','test'])
def test_hidden_labels_cannot_select_functional_exclusions(role):
    nodes,split=fixture()
    hidden=nodes.iloc[:1].assign(gene_id=split.entities(role)[0])
    with pytest.raises(ValueError,match='Forbidden'):closure(hidden,split)


def test_source_and_manifest_roundtrip_keep_the_declared_split_and_targets(tmp_path):
    nodes,split=fixture();excluded=closure(nodes,split)
    path=tmp_path/'split.json';write_split(path,split,excluded)
    assert read_split(path)==(split,excluded)
    assert set(excluded.targets)=={'pfam_id','recorded_profile'}
    with pytest.raises(ValueError,match='different benchmark'):excluded.guard_inputs(benchmark_id='other')


@pytest.mark.parametrize('source',[(),('unknown_annotation',)])
def test_unregistered_functional_source_is_refused(source):
    nodes,split=fixture()
    with pytest.raises(ValueError,match='actual functional'):
        F.make_exclusions(nodes,'recorded_profile',source_targets=source,benchmark_id=split.benchmark_id,split=split)


@pytest.mark.parametrize('item',[
    ('unknown','name',(('column','x'),)),('column','',(('column','x'),)),
    ('column','x',()),('column','x',(['column','x'],)),
    ('layer','x',(('column',' x'),)),('layer','x',(('column','x'),('column','x')))])
def test_malformed_lineage_addresses_are_refused(item):
    with pytest.raises(ValueError):F.DerivedInput(*item)


def test_duplicate_derived_input_addresses_are_refused():
    nodes,split=fixture();item=F.DerivedInput('column','x',(('column','raw_measurement'),))
    with pytest.raises(ValueError,match='Duplicate'):closure(nodes,split,derived_inputs=(item,item))


def test_supplemental_source_bans_reclose_registered_derived_columns(monkeypatch):
    from starplast import datasets
    nodes,split=fixture()
    original=datasets.derived_dependents
    def registered(columns,organism=None):
        return (*original(columns,organism),*(['attention_embedding'] if 'n_publications' in columns else []))
    monkeypatch.setattr(datasets,'derived_dependents',registered)
    excluded=closure(nodes,split)
    with pytest.raises(ValueError,match='contamination'):excluded.guard_inputs(columns=['attention_embedding'])


@pytest.mark.parametrize('lineage',[
    (F.DerivedInput('column','a',(('column','a'),)),),
    (F.DerivedInput('column','a',(('column','b'),)),F.DerivedInput('column','b',(('column','a'),))),
    (F.DerivedInput('column','a',(('column','unreviewed_unknown_source'),)),)])
def test_cyclic_and_unresolved_representation_lineage_is_refused(lineage):
    nodes,split=fixture()
    with pytest.raises(ValueError,match='Cyclic|Unresolved'):closure(nodes,split,derived_inputs=lineage)


def test_encoded_question_must_be_declared():
    nodes,split=fixture()
    with pytest.raises(ValueError,match='nonempty encoded'):
        F.make_exclusions(nodes,(),source_targets=('pfam_id',),benchmark_id=split.benchmark_id,split=split)
