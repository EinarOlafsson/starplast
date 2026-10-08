"""Execute exact source replay and functional Discoveries navigation against frozen data."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import pandas as pd
from starplast import claims as C, discovery_labels as D, functional_results as F, organisms as O, strategies as S, track_record as T
from scripts.audit_discovery_labels import _source_pairs


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(output):
    """Verify unchanged sources, complete memberships and actual per-class recovery views."""
    from PyQt6 import QtWidgets
    from starplast.discoveries_panel import DiscoveriesPanel
    output = Path(output)
    if output.exists():raise ValueError('Use a new immutable functional UI audit directory')
    output.mkdir(parents=True)
    old = json.loads((ROOT/'results/discoveries_label_coverage_2026_10_08_v4/manifest.json').read_text())
    for path,digest in old['input_sha256'].items():
        if Path(path).suffix=='.parquet':assert _sha(path)==digest,path
    code = [Path(__file__),ROOT/'scripts/audit_discovery_labels.py',
        *(ROOT/'starplast'/name for name in ('discovery_labels.py','discoveries_panel.py','functional_domains.py','functional_results.py'))]
    inputs = [ROOT/'starplast/data'/name for name in ('nodes.parquet','pf_nodes.parquet','claims.parquet',
        'claim_recipes.parquet','track_record.parquet','functional_domain_names.json','functional_results.json')]+code
    hashes = {str(path):_sha(path) for path in inputs}
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    summaries = []
    for organism in (O.TOXOPLASMA,O.FALCIPARUM):
        ctx = S.Context.shipped(organism)
        inventory,members = D.catalogue(ctx,C.shipped(organism),C.recipes(organism),T.shipped(organism))
        for target in inventory.target:
            group = members[members.target.eq(target)]
            if target in D.FUNCTION_FIELDS:expected = _source_pairs(ctx,target)
            elif target in ctx.nodes:expected = {(str(ctx.gene_ids[i]),str(v)) for i,v in ctx.truth(target).dropna().items()}
            else:expected = set()
            assert set(zip(group.gene_id,group.value))==expected,target
        for target in ('interpro_id','interpro_ids','pfam_id','pfam_ids'):
            if target not in ctx.nodes:continue
            original = D.annotation_members(ctx,target)
            named = members[members.target.eq(target)].reset_index(drop=True)
            pd.testing.assert_frame_equal(named[['organism','target','value','gene_id']],
                original[['organism','target','value','gene_id']],check_exact=True)
            assert named.source_description.tolist()==original.description.tolist()
            described = original.description.str.strip().ne('')
            assert named.loc[described,'description'].tolist()==original.loc[described,'description'].tolist()
        directory = output/organism;directory.mkdir()
        inventory.to_parquet(directory/'labels.parquet',index=False)
        members.to_parquet(directory/'memberships.parquet',index=False)
        panel = DiscoveriesPanel(organism,context=ctx)
        try:
            panel.resize(1300,900)
            panel.annotation_search.setText('kinase')
            domain_targets = set(panel.browsed_labels.target)&{'interpro_id','interpro_ids','pfam_id','pfam_ids'}
            assert domain_targets,organism
            panel.label.setCurrentIndex(panel.label.findData('pfam_id' if organism==O.TOXOPLASMA else 'pfam_ids'))
            assert panel.annotation_class.count()>1
            panel.annotation_class.setCurrentIndex(1)
            assert len(panel.shown_members)>0
            assert 'current nomenclature' in panel.annotation_note.text().lower()
            panel.show();app.processEvents()
            assert panel.grab().save(str(directory/'domain_provenance.png'))
            panel.annotation_search.clear()
            panel.label.setCurrentIndex(panel.label.findData('ec_number'))
            benchmarks,reason = F.shipped(organism)
            if benchmarks:
                panel.annotation_functional.click()
                benchmark = benchmarks[0]
                assert panel.tabs.currentWidget()==panel.functional_page
                assert panel.functional_shown_rows==benchmark.rows
                for card in benchmark.major_class_cards:
                    index = panel.functional_view.findData('major:'+card['class'])
                    assert index>=0
                    panel.functional_view.setCurrentIndex(index)
                    outcomes = [panel.functional_rows.item(i,4).text() for i in range(panel.functional_rows.rowCount())]
                    assert outcomes.count('recovered reference class')==card['true_positive']
                    assert outcomes.count('unexpected reference-class call')==card['false_positive']
                    assert sum('class missed' in text for text in outcomes)==card['false_negative']
                panel.functional_view.setCurrentIndex(panel.functional_view.findData('strategy:'))
                gene = panel.functional_shown_rows[0]['entity']
                chosen = [];panel.gene_chosen.connect(chosen.append)
                panel._functional_gene_clicked(0,0)
                assert chosen==[gene]
                app.processEvents();assert panel.grab().save(str(directory/'functional_scorecards.png'))
            else:
                assert not panel.annotation_functional.isEnabled()
                assert reason and panel.functional_rows.rowCount()==0
            summaries.append({'organism':organism,'labels':len(inventory),
                'functional_memberships':int(members.target.isin(D.FUNCTION_FIELDS).sum()),
                'domain_function_search_targets':sorted(domain_targets),
                'known_membership_and_original_descriptions':'exact source replay',
                'functional_held_out_rows':sum(len(b.rows) for b in benchmarks),
                'functional_artifact_identities':[b.metadata['source_artifact_identity'] for b in benchmarks],
                'independent_biological_benchmarks_added':0,'new_deployment_claims':0,
                'ui':'domain provenance, functional search, seven class outcomes and gene navigation verified' if benchmarks else
                    'domain provenance/function search verified; functional benchmark unavailable explicitly'})
        finally:panel.close()
    assert hashes=={str(path):_sha(path) for path in inputs}
    (output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (output/'code').mkdir()
    for path in code:(output/'code'/path.name).write_bytes(path.read_bytes())
    (output/'manifest.json').write_text(json.dumps({'input_sha256':hashes,
        'output_sha256':{str(path.relative_to(output)):_sha(path) for path in output.rglob('*') if path.is_file()},
        'interpretation':'Annotation recovery only; unknown independent activity accuracy, no calibrated probabilities or deployment claims'},indent=2)+'\n')
    return summaries


def main():
    """Save executed UI/source evidence, including failures in separate immutable diagnostics."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);args = parser.parse_args()
    notebook = ExecutedNotebook('Functional Discoveries: source membership, names and recovery scorecards')
    notebook.md('Replay all available labels against their original source memberships. Additional nomenclature names do not replace source descriptions, change gene assignments or establish measured activity. Exercise actual domain search/provenance and all seven EC class outcome denominators in the UI. Original source tables, claims, recipes and legacy held-out records retain their previous hashes.')
    try:
        notebook.code('from scripts.audit_functional_discoveries import audit',f'summary=audit({str(args.out)!r})','summary')
    except Exception as exc:
        notebook.md('Diagnostic: '+type(exc).__name__+': '+str(exc))
        notebook.write(str(args.out/'diagnostic.ipynb'));raise
    notebook.write(str(args.out/'audit.ipynb'))
    print(json.dumps(notebook.ns['summary'],indent=2))


if __name__=='__main__':main()
