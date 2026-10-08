"""Execute exact matrix/call diagnostics without selecting or changing a model."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
from scripts import freeze_functional_profile_controls as P

OUT = Path(__file__).resolve().parent
nb = ExecutedNotebook('Complete-Pfam native parity failure diagnostic')
nb.ns.update(ROOT=ROOT, OUT=OUT, P=P)
nb.md('Original failing script/tests were copied before execution. This bounded diagnostic compares full native inputs, orders, outputs and exact float representations. No settings or exclusions are selected from outcomes.')
try:
    nb.code('import numpy as np, pandas as pd, json, shutil',
        'from scripts import freeze_functional_profile_knn as K',
        'from starplast import label_records as L, strategies as S',
        'contract, source_scope, receipts = P.load_snapshot(P.SNAPSHOT)',
        'original = json.loads((P.SNAPSHOT / "input_manifest.json").read_text())',
        'node_paths = [path for path,digest in original.items() if digest == source_scope["installed_source_sha256"]]',
        'assert len(node_paths) == 1 and P._sha(node_paths[0]) == source_scope["installed_source_sha256"]',
        'receipts[node_paths[0]] = source_scope["installed_source_sha256"]',
        'code = [ROOT/"scripts/freeze_functional_profile_knn.py",ROOT/"scripts/freeze_functional_profile_controls.py",ROOT/"starplast/label_records.py",ROOT/"starplast/strategies.py",ROOT/"starplast/functional_exclusions.py",ROOT/"starplast/datasets.py",OUT/"run_diagnostic.py"]',
        '(OUT/"code").mkdir()',
        'for path in code:',
        '    receipts[str(path)] = P._sha(path)',
        '    shutil.copyfile(path,OUT/"code"/path.name)',
        'P._write(OUT/"input_manifest.json",receipts)',
        'from starplast.splits import ExclusionManifest',
        'previous = ROOT/"results/functional_profile_knn_diagnostic_2026_10_08"',
        'for name in ("features.parquet","exclusions.json","feature_provenance.json"):',
        '    receipts[str(previous/name)] = P._sha(previous/name)',
        'features = pd.read_parquet(previous/"features.parquet")',
        'exclusion_data = json.loads((previous/"exclusions.json").read_text())',
        'for key in ("targets","columns","layers","fit_entities"):',
        '    exclusion_data[key] = tuple(exclusion_data[key])',
        'exclusions = ExclusionManifest(**exclusion_data)',
        'provenance = json.loads((previous/"feature_provenance.json").read_text())',
        'assert tuple(features.index) == tuple(a.entity for a in contract.split.assignments)',
        'assert exclusions.split_identity == contract.split.identity and exclusions.fit_entities == contract.split.entities("train")',
        'assert features.columns.tolist() == provenance["selected_columns"]',
        'P._write(OUT/"input_manifest.json",receipts)',
        'labels = contract.training_labels()',
        'features.to_parquet(OUT/"features.parquet")',
        'labels.to_frame().to_parquet(OUT/"training_labels.parquet")',
        'P._write(OUT/"exclusions.json",K._json(__import__("dataclasses").asdict(exclusions)))',
        'P._write(OUT/"feature_provenance.json",provenance)',
        '{"whole_genes":len(contract.rows),"eligible":len(features),"train":len(labels),"features":len(features.columns),"source":"Pinned prepared inputs from interrupted diagnostic; no whole-node frame retained"}')
    nb.code('captured = []',
        'native_vote = S.knn_vote',
        'def capture(matrix,visible,k,query=None):',
        '    captured.append((matrix.copy(),visible.copy(),np.array(query).copy()))',
        '    return native_vote(matrix,visible,k,query=query)',
        'S.knn_vote = capture',
        'try:',
        '    batch = L.feature_knn(features,labels,split=contract.split,exclusions=exclusions,k=15,min_share=0.3)',
        '    calls,support,scores,votes = K.native_replay(features,labels,contract.split,batch.model_state)',
        'finally:',
        '    S.knn_vote = native_vote',
        'assert len(captured) == 2',
        'P._write(OUT/"model_state.json",K._json(batch.model_state))',
        'batch.rows.to_parquet(OUT/"adapter_rows.parquet")',
        'pd.DataFrame({"entity":list(contract.split.entities("test")),"prediction":calls.tolist(),"support":support.tolist()}).to_parquet(OUT/"replay_rows.parquet",index=False)',
        'batch.class_scores.to_parquet(OUT/"adapter_class_scores.parquet")',
        'scores.to_parquet(OUT/"replay_class_scores.parquet")',
        'for name,data in zip(("adapter","replay"),captured):',
        '    np.save(OUT/(name+"_matrix.npy"),data[0],allow_pickle=False)',
        '    np.save(OUT/(name+"_query.npy"),data[2],allow_pickle=False)',
        '    data[1].to_frame("visible_label").to_parquet(OUT/(name+"_visible_labels.parquet"))',
        '{"captured_matrices":len(captured),"test_rows":len(batch.rows)}')
    nb.code('a,b = captured[0][0],captured[1][0]',
        'different = np.argwhere(a != b)',
        'changed_calls = [i for i,(x,y) in enumerate(zip(batch.rows.prediction.tolist(),calls.tolist())) if x != y]',
        'normalized_calls = batch.rows.prediction.astype(object).where(batch.rows.prediction.notna(),None)',
        'normalized_differences = [i for i,(x,y) in enumerate(zip(normalized_calls.tolist(),calls.tolist())) if x != y]',
        'support_a,support_b = batch.rows.support.to_numpy(dtype=float),support.to_numpy(dtype=float)',
        'support_different = np.flatnonzero(~((support_a == support_b)|(np.isnan(support_a)&np.isnan(support_b))))',
        'score_a,score_b = batch.class_scores.to_numpy(),scores.to_numpy()',
        'score_different = np.argwhere(~((score_a == score_b)|(np.isnan(score_a)&np.isnan(score_b))))',
        'summary = {"matrix_shape":list(a.shape),"matrix_differences":len(different),"query_order_exact":captured[0][2].tolist()==captured[1][2].tolist(),"visible_labels_exact":captured[0][1].equals(captured[1][1]),"prediction_representation_differences":len(changed_calls),"prediction_differences":len(normalized_differences),"support_differences":len(support_different),"score_differences":len(score_different),"score_columns_exact":batch.class_scores.columns.tolist()==scores.columns.tolist(),"settings_changed":False,"biological_admission":False}',
        'P._write(OUT/"matrix_differences.json",[{"row":int(i),"column":int(j),"entity":features.index[i],"feature":batch.model_state["kept_columns"][j],"adapter":float(a[i,j]),"replay":float(b[i,j]),"adapter_hex":float(a[i,j]).hex(),"replay_hex":float(b[i,j]).hex()} for i,j in different])',
        'P._write(OUT/"prediction_differences.json",[{"row":i,"entity":batch.rows.entity.iloc[i],"adapter":normalized_calls.iloc[i],"replay":calls.iloc[i],"adapter_type":type(batch.rows.prediction.iloc[i]).__name__,"replay_type":type(calls.iloc[i]).__name__,"adapter_repr":repr(batch.rows.prediction.iloc[i]),"replay_repr":repr(calls.iloc[i])} for i in changed_calls])',
        'P._write(OUT/"support_differences.json",[{"row":int(i),"adapter_hex":float(support_a[i]).hex(),"replay_hex":float(support_b[i]).hex()} for i in support_different])',
        'P._write(OUT/"score_differences.json",[{"row":int(i),"column":int(j),"adapter_hex":float(score_a[i,j]).hex(),"replay_hex":float(score_b[i,j]).hex()} for i,j in score_different])',
        'P._write(OUT/"summary.json",summary)',
        'for path,digest in receipts.items():',
        '    assert P._sha(path) == digest',
        'summary')
    print(json.dumps(nb.ns['summary'],indent=2))
except Exception as exc:
    P._write(OUT/'failure.json',{'error':type(exc).__name__,'detail':str(exc)})
    nb.md('Diagnostic execution failed: '+type(exc).__name__+': '+str(exc))
    raise
finally:
    nb.write(str(OUT/'executed.ipynb'))
    (OUT/'SHA256SUMS.txt').write_text('\n'.join(f'{P._sha(path)}  {path.relative_to(OUT)}'
        for path in sorted(OUT.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt')+'\n')
