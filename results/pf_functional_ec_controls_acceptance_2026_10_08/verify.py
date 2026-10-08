"""Independently verify the frozen Pf control packet and exact population arithmetic."""
import hashlib, json, sys
from pathlib import Path
ROOT=Path('/media/carruthers/mnt3/claude/repo/starplast')
sys.path.insert(0,str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
OUT=Path(__file__).resolve().parent
PACKET=ROOT/'results/pf_functional_ec_controls_2026_10_08'


def sha(path):
    with Path(path).open('rb') as file:return hashlib.file_digest(file,'sha256').hexdigest()


def verify():
    import pandas as pd
    lines=(PACKET/'SHA256SUMS.txt').read_text().splitlines()
    for line in lines:
        digest,name=line.split('  ',1)
        assert sha(PACKET/name)==digest,name
    inputs=json.loads((PACKET/'input_manifest.json').read_text())
    for path,digest in inputs.items():assert sha(path)==digest,path
    summary=json.loads((PACKET/'summary.json').read_text())
    rows=pd.read_parquet(PACKET/'control_rows.parquet')
    assert len(rows)==rows.entity.nunique()==152
    assert (~rows.training_supported).sum()==4
    assert rows.majority.eq(rows.truth).sum()==43
    assert rows.prevalence_call.eq(rows.truth).sum()==35
    assert summary['controls']['training_majority']['metrics']['accuracy']==43/152
    assert summary['controls']['training_prevalence']['metrics']['accuracy']==35/152
    assert summary['roles']=={'train':600,'tune':151,'calibration':147,'test':152}
    direct=pd.read_parquet(PACKET/'direct_profiles.parquet')
    assert len(direct)==5720 and int(direct.eligible.sum())==1050
    assert direct.status.value_counts().to_dict()=={'unannotated':4500,'complete':1050,'unresolved':170}
    assert summary['missing_groups']==summary['fallback_collision_count']==0
    assert summary['biological_admission'] is False
    return {'outputs':len(lines),'inputs':len(inputs),'test_genes':152,'unsupported_genes':4,
        'majority_correct':43,'prevalence_correct':35,'exact_population_arithmetic':True}


nb=ExecutedNotebook('Independent exact acceptance of Plasmodium EC profile controls')
nb.ns.update(verify=verify)
nb.md('Verify every frozen input/output byte, then recount profile/control outcomes. This is source-annotation arithmetic, without inference-strategy fitting or independent biological accuracy.')
nb.code('result=verify()','result')
(OUT/'summary.json').write_text(json.dumps(nb.ns['result'],indent=2)+'\n')
nb.write(str(OUT/'executed.ipynb'))
(OUT/'SHA256SUMS.txt').write_text('\n'.join(f'{sha(path)}  {path.relative_to(OUT)}' for path in sorted(OUT.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt')+'\n')
print(json.dumps(nb.ns['result']))
