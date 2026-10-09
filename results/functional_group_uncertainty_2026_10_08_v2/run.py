from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()))
import hashlib
import json
import shutil
from scripts.notebook_runner import ExecutedNotebook
ROOT=Path.cwd();OUT=ROOT/'results/functional_group_uncertainty_2026_10_08_v2'
nb=ExecutedNotebook('Frozen paired source-profile whole-group uncertainty')
nb.ns.update(ROOT=ROOT,OUT=OUT)
nb.md('Predeclared: all152 genes/four unsupported/145 recorded groups; 500 whole-group draws, seed20261008, 2.5/97.5 percentiles with linear interpolation. Existing point estimates remain unchanged. Conditional on supplied groups; source/biological independence and per-gene calibration remain unresolved. No fitting, tuning or model selection.')
try:
 nb.code("from pathlib import Path", "import hashlib,json,numpy as np", "from starplast import group_ratios as G,functional_agreement as A", "from starplast.scorecard_view import export_scorecard,render_scorecard_html", "def sha(path):\n with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()", "source=ROOT/'starplast/data/functional_agreement.json'", "assert sha(source)==A.REPORT_SHA256", "report=json.loads(source.read_text())", "A.build_scorecard(report)", "paths=[source,ROOT/'starplast/group_ratios.py',ROOT/'tests/test_group_ratios.py',ROOT/'starplast/functional_agreement.py',ROOT/'starplast/scorecard_view.py',ROOT/'scripts/notebook_runner.py']", "inputs={str(path):sha(path) for path in paths}", "(OUT/'input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\\n')", "state,arrays=G.estimate(*G.paired_terms(report),seed=20261008,draws=500)", "assert state['entities']==152 and state['recorded_groups']==145", "assert report['counts']['unsupported']==4", "np.savez_compressed(OUT/'draws.npz',**arrays)", "reloaded=np.load(OUT/'draws.npz')", "replay,replay_arrays=G.estimate(*G.paired_terms(report),seed=20261008,draws=500,selections=reloaded['selections'])", "assert replay==state", "for name in arrays:np.testing.assert_array_equal(arrays[name],reloaded[name]);np.testing.assert_array_equal(arrays[name],replay_arrays[name])", "view=G.paired_scorecard(report,state)", "(OUT/'uncertainty.json').write_text(json.dumps(state,sort_keys=True,indent=2,allow_nan=False)+'\\n')", "(OUT/'scorecard.json').write_text(export_scorecard(view))", "(OUT/'scorecard.html').write_text(render_scorecard_html(view,expanded=True))", "assert json.loads((OUT/'scorecard.json').read_text())['snapshot']['group_uncertainty']==state", "assert all(sha(path)==digest for path,digest in inputs.items())", "{'groups':state['recorded_groups'],'draws':state['draws_evaluated'],'rates':state['rates']}")
finally:
 nb.write(str(OUT/'executed.ipynb'))
 (OUT/'SHA256SUMS.txt').write_text(''.join(f'{nb.ns["sha"](p)}  {p.relative_to(OUT)}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(json.dumps(nb.ns['state']['rates'],indent=2))
