"""Freeze and execute the Discoveries label expansion and source-membership audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

import pandas as pd
from starplast import claims as C, discovery_labels as D, organisms as O, strategies as S, track_record as T


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _source_pairs(ctx,target):
    # Independent source replay uses delimited source tokens, not the browser's
    # regex extraction or membership function. Reject unexpected formats.
    pairs = set()
    for i,value in ctx.truth(target).dropna().items():
        for token in str(value).split(';'):
            token = token.strip()
            if target.startswith('ec_number'):
                term = token.split()[0]
                parts = term.split('.')
                assert len(parts)==4 and parts[0] in {'1','2','3','4','5','6','7','-'},(target,token)
                assert all(p=='-' or p.isdigit() for p in parts[1:]),(target,token)
            else:
                term = token
                pattern = r'PF[0-9]{5}(?:\.[0-9]+)?' if target.startswith('pfam') else r'IPR[0-9]{6}'
                assert re.fullmatch(pattern,term),(target,token)
            pairs.add((str(ctx.gene_ids[i]),term))
    return pairs


def freeze(output):
    """Audit every available label and functional term against pinned installed inputs."""
    from PyQt6 import QtWidgets
    from starplast.discoveries_panel import DiscoveriesPanel
    output = Path(output)
    if output.exists():raise ValueError('Use a new immutable audit directory')
    output.mkdir(parents=True)
    inputs = [Path(O.nodes_path(o)) for o in (O.TOXOPLASMA,O.FALCIPARUM)]
    inputs += [ROOT/'starplast/data'/name for name in ('claims.parquet','claim_recipes.parquet','track_record.parquet')]
    code = [Path(__file__),ROOT/'starplast/discovery_labels.py',ROOT/'starplast/discoveries_panel.py']
    inputs += code
    hashes = {str(p):_sha(p) for p in inputs}
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    summaries = []
    for organism in (O.TOXOPLASMA,O.FALCIPARUM):
        ctx = S.Context.shipped(organism)
        claims,recipes,records = C.shipped(organism),C.recipes(organism),T.shipped(organism)
        inventory,members = D.catalogue(ctx,claims,recipes,records)
        assert set(claims.target)|set(recipes.target)|set(records.target)<=set(inventory.target)
        assert inventory.claims.sum()==len(claims)
        assert inventory.held_out_rows.sum()==len(records)
        assert not members.duplicated(['organism','target','value','gene_id']).any()
        checks = {}
        for target in inventory.target:
            group = members[members.target.eq(target)]
            if target in D.FUNCTION_FIELDS:
                expected = _source_pairs(ctx,target)
            elif target in ctx.nodes:
                truth = ctx.truth(target).dropna()
                expected = {(str(ctx.gene_ids[i]),str(v)) for i,v in truth.items()}
            else:expected = set()
            actual = set(zip(group.gene_id,group.value))
            assert actual==expected,target
            checks[target] = {'membership_rows':len(actual),'known_genes':group.gene_id.nunique(),
                'classes':group.value.nunique(),'source_membership':'exact independent replay'}
        old_default = C.discoveries(organism,min_confidence=.8,min_lift=2.)
        directory = output/organism
        directory.mkdir()
        inventory.to_parquet(directory/'labels.parquet',index=False)
        members.to_parquet(directory/'memberships.parquet',index=False)
        (directory/'source_replay.json').write_text(json.dumps(checks,indent=2)+'\n')
        # Exercise the real widget, both a named function and direct class/gene
        # navigation. Synthetic precision is never attached to known members.
        panel = DiscoveriesPanel(organism,context=ctx)
        try:
            assert set(panel.inventory.target)==set(inventory.target)
            assert len(panel.shown)==len(old_default)
            panel.annotation_search.setText('kinase')
            assert len(panel.browsed_labels)>0
            assert panel.browsed_labels.target.str.startswith(('interpro','ec_number')).any()
            target = 'interpro_id' if organism==O.TOXOPLASMA else 'ec_number'
            panel.label.setCurrentIndex(panel.label.findData(target))
            assert panel.annotation_class.count()>1
            panel.annotation_class.setCurrentIndex(1)
            assert len(panel.shown_members)>0
            chosen = [];panel.gene_chosen.connect(chosen.append)
            panel._member_clicked(0,0)
            assert chosen==[panel.shown_members.gene_id.iloc[0]]
            assert not panel.annotation_record.isEnabled()
            assert 'not evaluated' in panel.annotation_note.text().lower()
            panel.resize(1200,850);panel.show();app.processEvents()
            assert panel.grab().save(str(directory/'function_browser.png'))
        finally:panel.close()
        functional = inventory[inventory.target.isin(D.FUNCTION_FIELDS)]
        summaries.append({'organism':organism,'genes':ctx.n,
            'before_menu_labels':int(claims.target.nunique()),'before_default_visible_claims':len(old_default),
            'before_default_visible_targets':sorted(old_default.target.unique().tolist()),
            'after_available_labels':len(inventory),'families':inventory.groupby('family').size().to_dict(),
            'functional_variables':functional.target.tolist(),'functional_classes':int(functional.classes.sum()),
            'functional_memberships':int(members.target.isin(D.FUNCTION_FIELDS).sum()),
            'functional_inferred_claims':int(functional.claims.sum()),'functional_held_out_rows':int(functional.held_out_rows.sum()),
            'legacy_claims_unchanged':len(claims),'legacy_held_out_rows_unchanged':len(records),
            'independent_biological_benchmarks_added':0,'source_membership':'exact for every available variable',
            'ui':'function search, class membership, gene signal and missing scorecards verified'})
    assert hashes=={str(p):_sha(p) for p in inputs},'Inputs changed during audit'
    (output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (output/'code').mkdir()
    for p in code:(output/'code'/p.name).write_bytes(p.read_bytes())
    manifest = {'code_version':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'+discovery-browser-worktree',
        'input_sha256':hashes,'output_sha256':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*') if p.is_file()},
        'interpretation':'Known annotations and legacy records; no new functional predictions or biological accuracy admission'}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return summaries


def verify(output):
    """Re-read stored tables and verify every immutable input/output hash and count."""
    output = Path(output)
    manifest = json.loads((output/'manifest.json').read_text())
    for p,digest in manifest['input_sha256'].items():assert _sha(p)==digest,p
    for p,digest in manifest['output_sha256'].items():assert _sha(output/p)==digest,p
    for summary in json.loads((output/'summary.json').read_text()):
        directory = output/summary['organism']
        inventory = pd.read_parquet(directory/'labels.parquet')
        members = pd.read_parquet(directory/'memberships.parquet')
        assert len(inventory)==summary['after_available_labels']
        for target,row in inventory.set_index('target').iterrows():
            group = members[members.target.eq(target)]
            assert group.gene_id.nunique()==row.annotated_genes
            assert group.value.nunique()==row.classes
            assert row.annotated_genes+row.unannotated_genes==summary['genes']
        assert not members.duplicated(['organism','target','value','gene_id']).any()
    return {'source_and_code_hashes':'verified','stored_membership_counts':'exact','new_functional_predictions':0,
        'independent_biological_benchmarks_added':0}


def main():
    """Execute an annotated before/after audit, retaining diagnostic failures separately."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Discoveries across functions and labels: actual coverage and gaps')
    nb.md('68.01 expands browsing of existing annotations and legacy inference/test coverage. Function-first ordering, multi-valued InterPro/Pfam and EC classes, phenotypes, stages and annotation flags. Known membership is not a new prediction; unannotated and false annotation flags do not establish biological absence. Source IDs remain organism-specific; absent provenance is unresolved. Legacy scorecards retain their historical truth scope. 68.02 remains responsible for functional benchmarks and precomputed inference.')
    try:
        nb.code('from scripts.audit_discovery_labels import freeze, verify',f'summary=freeze({str(args.out)!r})','summary')
        nb.md('Verify immutable source/code/output hashes and every stored gene/class count. Original nodes, claims, recipes and held-out records must remain unchanged. Functional biological accuracy stays unknown.')
        nb.code(f'verification=verify({str(args.out)!r})','verification')
    except Exception as exc:
        nb.md('Diagnostic: '+type(exc).__name__+': '+str(exc))
        nb.write(str(args.out/'diagnostic.ipynb'))
        raise
    nb.write(str(args.out/'audit.ipynb'))
    print(json.dumps(nb.ns['summary'],indent=2))
    print(nb.ns['verification'])


if __name__=='__main__':main()
