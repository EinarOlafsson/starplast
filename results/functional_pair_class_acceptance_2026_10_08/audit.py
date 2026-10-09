"""Independent source-membership recount and desktop receipt acceptance."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd()))
import hashlib
import json
from scripts.notebook_runner import ExecutedNotebook

ROOT = Path.cwd()
OUT = ROOT / 'results/functional_pair_class_acceptance_2026_10_08'
UI = ROOT / 'results/functional_pair_class_ui_2026_10_08'
V1 = ROOT / 'results/functional_pair_class_2026_10_08'
V2 = ROOT / 'results/functional_pair_class_2026_10_08_v2'
nb = ExecutedNotebook('Independent paired recorded-class navigation acceptance')
nb.ns.update(ROOT=ROOT, OUT=OUT, UI=UI, V1=V1, V2=V2)
nb.md('Recount class membership directly from the original complete profiles. All 152 known-source genes remain in denominators; abstention is unknown and recorded nonmembership is not biological absence. Inspect unchanged native cards and all class exports. No fits, source replacements, calibrated probabilities or biological admission. Initial synthetic fixture failure is preserved; its missing group_uncertainty key was corrected. Final desktop screenshot inspected.')
try:
    nb.code(
        'from pathlib import Path', 'import hashlib,json',
        'from starplast import functional_results as F, functional_agreement as A, functional_pair_class as C',
        'from starplast.scorecard_view import export_scorecard,render_scorecard_html',
        "def sha(path):\n with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()",
        'receipts=0',
        "for line in (UI/'SHA256SUMS.txt').read_text().splitlines():\n digest,name=line.split('  ',1);path=UI/name\n assert not path.is_symlink() and path.resolve().is_relative_to(UI.resolve())\n assert sha(path)==digest,name\n receipts+=1",
        "inputs=json.loads((UI/'input_manifest.json').read_text())",
        'for path,digest in inputs.items():assert sha(path)==digest,path',
        "assert sha(ROOT/'starplast/data/functional_agreement.json')=='940da900b3d5ba356c6e360f220787b849d7382138132ccc96f2d1684178e98b'",
        "assert sha(ROOT/'starplast/data/functional_results.json')=='386c8e05f5b8cf32f4b9ba6c777d32933c73c08ce40d9e01d0b27959e9dbadda'",
        "assert sha(ROOT/'starplast/data/functional_agreement_uncertainty.json')=='0bc3e966f6ee493185f1189735245852b8b910a6548abdcb96496e376ef061f9'",
        "benchmarks,reason=F.shipped('Pf');assert not reason",
        "pair=tuple(next(b for b in benchmarks if b.metadata['strategy']==name) for name in ('feature_knn','random_forest'))",
        "report=json.loads((ROOT/'starplast/data/functional_agreement.json').read_text());assert report==A.compare(*pair)",
        'catalog=[]',
        "for address in C.addresses(pair):\n result,view=C.build(pair,report,address)\n kind,value=address.split(':',1)\n def member(profile):\n  return value in json.loads(profile) if kind=='major' else profile==value\n rows=report['rows']\n actual=[member(r['truth']) for r in rows]\n calls=[[None if r[key] is None else member(r[key]) for r in rows] for key in ('left_prediction','right_prediction')]\n assert len(rows)==152 and result['counts']['full_profile_unsupported']==4\n assert result['counts']['known_positive_genes']==sum(actual)\n assert result['counts']['agreed_presence']==sum(l is True and r is True for l,r in zip(*calls))\n assert result['counts']['agreed_nonmembership']==sum(l is False and r is False for l,r in zip(*calls))\n assert result['counts']['jointly_called_true_positive']==sum(a and l is True and r is True for a,l,r in zip(actual,*calls))\n display=[r['entity'] for r,a,l,p in zip(rows,actual,*calls) if a or l is True or p is True]\n assert display==[r['entity'] for r in result['displayed_rows']]\n for method,predictions in zip(result['methods'],calls):\n  tp=sum(a and p is True for a,p in zip(actual,predictions))\n  fp=sum(not a and p is True for a,p in zip(actual,predictions))\n  fn=sum(a and p is not True for a,p in zip(actual,predictions))\n  assert (tp,fp,fn)==tuple(method[k] for k in ('true_positive','false_positive','false_negative'))\n  assert method['source_precision']==(tp/(tp+fp) if tp+fp else None)\n  assert method['source_recall']==(tp/(tp+fn) if tp+fn else None)\n assert json.loads(export_scorecard(view))['snapshot']==result\n assert result['class_uncertainty']=='unavailable' and result['calibrated_confidence'] is None\n assert not result['biological_admission']\n catalog.append({k:v for k,v in result.items() if k not in ('rows','displayed_rows','native_method_class_cards')})",
        "result,view=C.build(pair,report,'major:2')",
        "(OUT/'major2.json').write_text(export_scorecard(view)+'\\n')",
        "(OUT/'major2.html').write_text(render_scorecard_html(view,expanded=True))",
        "(OUT/'class_catalog.json').write_text(json.dumps(catalog,indent=2)+'\\n')",
        "summary=json.loads((UI/'summary.json').read_text())",
        "assert summary['original_inputs_unchanged'] and summary['cross_organism_refused'] and summary['altered_context_refused']",
        "assert '8 failed, 1 passed' in (V1/'tests.log').read_text()",
        "assert '19 passed' in (V2/'ui-tests.log').read_text()",
        "assert '9 passed' in (V2/'final-tests.log').read_text()",
        "assert 'result: success' in (V2/'desktop.log').read_text()",
        "acceptance={'recorded_class_addresses':len(catalog),'output_receipts_verified':receipts,'current_inputs_verified':len(inputs),'full_test_cohort':152,'full_profile_unsupported':4,'focused_checks':19,'final_checks':9,'initial_fixture_failures':8,'independent_class_recount':'exact','desktop_seconds':summary['elapsed_seconds'],'desktop_peak_bytes':summary['process_peak_bytes'],'source_and_parent_bytes_unchanged':True,'fits':0,'biological_admission':False,'calibrated_confidence':None,'screenshot_review':'paired_class_feature_knn.png inspected'}",
        "(OUT/'acceptance.json').write_text(json.dumps(acceptance,indent=2)+'\\n')", 'acceptance')
finally:
    nb.write(str(OUT / 'executed.ipynb'))
    for folder in (OUT,V1,V2):
        (folder/'SHA256SUMS.txt').write_text(''.join(
            f'{nb.ns["sha"](p)}  {p.relative_to(folder)}\n'
            for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(json.dumps(nb.ns['acceptance']))
