from pathlib import Path
import hashlib
import json
from scripts.notebook_runner import ExecutedNotebook

ROOT=Path.cwd()
OUT=ROOT/'results/pf_functional_agreement_ui_acceptance_2026_10_08'
DIAG=ROOT/'results/pf_functional_agreement_ui_diagnostics_2026_10_08'
UI=ROOT/'results/pf_functional_agreement_ui_2026_10_08_v2'

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

nb=ExecutedNotebook('Independent paired desktop integration receipts')
nb.ns.update(ROOT=ROOT,OUT=OUT,DIAG=DIAG,UI=UI,Path=Path,json=json,sha=sha)
nb.md('Verify immutable installation and both desktop packets, the externally pinned shipped report, all current canonical desktop inputs and recorded test results. Original node/claim/benchmark sources are unchanged. Screenshot paired_feature_knn.png was inspected: rates have fractions/percentages and named method outcomes. No new fitting or biological admission.')
try:
 nb.code("receipts=0", "for folder in (ROOT/'results/pf_functional_agreement_install_2026_10_08',ROOT/'results/pf_functional_agreement_ui_2026_10_08',UI):\n for line in (folder/'SHA256SUMS.txt').read_text().splitlines():\n  digest,name=line.split('  ',1);path=folder/name\n  assert not path.is_symlink() and path.resolve().is_relative_to(folder.resolve())\n  assert sha(path)==digest,name\n  receipts+=1", "inputs=json.loads((UI/'input_manifest.json').read_text())", "for name,digest in inputs.items():assert sha(name)==digest,name", "report_path=ROOT/'starplast/data/functional_agreement.json'", "assert sha(report_path)=='940da900b3d5ba356c6e360f220787b849d7382138132ccc96f2d1684178e98b'", "assert report_path.read_bytes()==(ROOT/'results/pf_functional_agreement_2026_10_08_v2/report.json').read_bytes()", "report=json.loads(report_path.read_text())", "summary=json.loads((UI/'summary.json').read_text())", "assert summary['original_inputs_unchanged'] is True and summary['cross_organism_refused'] is True and summary['altered_context_refused'] is True", "assert len(summary['benchmarks'])==2", "for item in summary['benchmarks']:\n assert item['paired_counts']==report['counts']\n assert item['paired_counts']['eligible']==152 and item['paired_counts']['unsupported']==4", "first=(DIAG/'ui-tests.log').read_text()", "assert '34 passed, 1 error' in first and 'QDialog has been deleted' in first", "broad=(DIAG/'ui-tests-v2.log').read_text()", "final=(DIAG/'ui-tests-v3.log').read_text()", "assert '68 passed' in broad and 'result: success' in broad", "assert '31 passed' in final and 'result: success' in final", "result={'output_receipts_verified':receipts,'current_desktop_inputs_verified':len(inputs),'shipped_report_matches_canonical_bytes':True,'focused_checks':68,'final_presentation_checks':31,'test_genes':152,'unsupported_genes':4,'screenshot_review':'inspected paired_feature_knn.png; fractions and named methods visible','desktop_elapsed_seconds':summary['elapsed_seconds'],'desktop_process_peak_bytes':summary['process_peak_bytes'],'initial_failure':'test teardown attempted to close a deleted dialog; fixed with Qt lifetime check','biological_admission':False,'calibrated_confidence':None,'fits_performed':0}", "(OUT/'acceptance.json').write_text(json.dumps(result,indent=2)+'\\n')", "result")
finally:
 nb.write(str(OUT/'executed.ipynb'))
 (OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
 (DIAG/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(DIAG)}\n' for p in sorted(DIAG.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(json.dumps(nb.ns['result']))
