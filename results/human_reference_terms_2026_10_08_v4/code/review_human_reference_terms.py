"""Capture bounded authoritative reference and public-data terms without pack admission.

Every response is retained with its official URL, UTC date, byte count and hash.
The existing human foundation, original acquisition inputs and installed files
are verified before and after this documentary review. A JavaScript shell alone
cannot establish GTEx redistribution permission; only recovered page-component
text can support that finding. Failures remain gaps while independent pages run.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import hashlib
import html
import json
from pathlib import Path
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
FOUNDATION=ROOT/'results/human_gene_space_foundation_2026_10_08'
MAX_BYTES=8*1024*1024
MAX_REQUESTS=16
URLS={'gtex_license':'https://gtexportal.org/home/license',
    'uniprot_license':'https://rest.uniprot.org/help/license',
    'gencode39':'https://www.gencodegenes.org/human/release_39.html',
    'gencode47':'https://www.gencodegenes.org/human/release_47.html'}


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def _write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def foundation_receipts():
    """Verify all frozen foundation outputs and unchanged non-code original inputs."""
    manifest=json.loads((FOUNDATION/'manifest.json').read_text());receipts=[]
    for relative,digest in manifest['outputs'].items():
        path=FOUNDATION/relative
        if _sha(path)!=digest:raise ValueError('Foundation output changed: '+relative)
        receipts.append({'path':str(path),'sha256':digest,'bytes':path.stat().st_size,'scope':'foundation_output'})
    for original,digest in manifest['input_sha256'].items():
        path=Path(original)
        # Current code evolves; its original bytes are already verified above.
        if path.suffix=='.py':continue
        if _sha(path)!=digest:raise ValueError('Original source or installed input changed: '+original)
        receipts.append({'path':original,'sha256':digest,'bytes':path.stat().st_size,'scope':'original_noncode_input'})
    path=FOUNDATION/'manifest.json'
    receipts.append({'path':str(path),'sha256':_sha(path),'bytes':path.stat().st_size,'scope':'foundation_manifest'})
    return receipts


def _capture(output,key,url,receipts):
    if len(receipts)>=MAX_REQUESTS:raise ValueError('Documentary request count exhausted')
    now=datetime.now(timezone.utc).isoformat()
    receipt={'key':key,'requested_url':url,'retrieved_utc':now,'max_response_bytes':MAX_BYTES}
    try:
        request=urllib.request.Request(url,headers={'User-Agent':'starplast-documentary-review/1','Accept-Encoding':'identity'})
        with urllib.request.urlopen(request,timeout=25) as response:
            chunks=[];size=0
            while True:
                block=response.read(min(65536,MAX_BYTES+1-size))
                if not block:break
                size+=len(block)
                if size>MAX_BYTES:raise ValueError('Response exceeds 8 MiB documentary bound')
                chunks.append(block)
            data=b''.join(chunks)
            if urllib.parse.urlparse(response.url).hostname not in {'gtexportal.org','www.gtexportal.org','rest.uniprot.org','www.gencodegenes.org'}:
                raise ValueError('Response redirected outside the declared official sites')
            receipt.update(final_url=response.url,http_status=response.status,content_type=response.headers.get('Content-Type',''))
        path=output/'captures'/(key+'.txt');path.write_bytes(data)
        receipt.update(status='captured',path=str(path.relative_to(output)),bytes=len(data),sha256=_sha(path))
        text=data.decode('utf-8',errors='replace')
    except Exception as error:
        receipt.update(status='unavailable',reason=type(error).__name__+': '+str(error));text=''
    receipts.append(receipt);_write(output/'request_receipts.json',receipts)
    return text


def capture(output, reuse=None):
    """Freeze the documentary scope, verify original bytes and capture independent pages."""
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable documentary directory')
    output.mkdir(parents=True);(output/'captures').mkdir()
    _write(output/'scope.json',{'review':'64.21 existing human reference and public redistribution terms',
        'urls':URLS,'response_limit_bytes':MAX_BYTES,'request_limit':MAX_REQUESTS,'runtime_memory_cap_bytes':400*1024*1024,
        'biological_admission':False,'pack_promoted':False,'source_replacement':False,
        'excluded':['Restricted activity truth','Biological fitting','Genome-wide completeness','New measurement acquisition']})
    before=foundation_receipts();_write(output/'foundation_receipts.json',before)
    receipts=[];texts={}
    if reuse is not None:
        reuse=Path(reuse);manifest=json.loads((reuse/'manifest.json').read_text())
        old=json.loads((reuse/'request_receipts.json').read_text())
        for key,url in URLS.items():
            receipt=next(row for row in old if row['key']==key and row.get('status')=='captured' and row['requested_url']==url)
            path=reuse/receipt['path']
            if _sha(path)!=receipt['sha256'] or manifest['files'][receipt['path']]['sha256']!=receipt['sha256']:
                raise ValueError('Reused authoritative capture changed')
            destination=output/receipt['path'];destination.write_bytes(path.read_bytes())
            receipts.append({**receipt,'reused_from':str(reuse),'reused_capture_sha256':receipt['sha256']})
            texts[key]=destination.read_text()
    else:texts={key:_capture(output,key,url,receipts) for key,url in URLS.items()}
    shell=texts['gtex_license'];queue=[]
    for relative in re.findall(r'<script\b[^>]*\bsrc=["\']([^"\']+)["\']',shell,re.I):
        url=urllib.parse.urljoin(URLS['gtex_license'],html.unescape(relative))
        parsed=urllib.parse.urlparse(url)
        if (parsed.hostname in {'gtexportal.org','www.gtexportal.org'} and parsed.path.startswith('/js/')
                and '.js' in parsed.path and 'vendors' not in parsed.path and 'legacy' not in parsed.path):queue.append(url)
    seen=set();bundles=[]
    while queue and len(receipts)<MAX_REQUESTS:
        url=queue.pop(0)
        if url in seen:continue
        seen.add(url)
        old_bundle=next((receipt for receipt in old if receipt.get('status')=='captured' and receipt['requested_url']==url),None) if reuse is not None else None
        if old_bundle is not None:
            path=reuse/old_bundle['path']
            if _sha(path)!=old_bundle['sha256'] or manifest['files'][old_bundle['path']]['sha256']!=old_bundle['sha256']:
                raise ValueError('Reused frontend capture changed')
            destination=output/old_bundle['path'];destination.write_bytes(path.read_bytes())
            receipts.append({**old_bundle,'reused_from':str(reuse),'reused_capture_sha256':old_bundle['sha256']})
            text=destination.read_text()
        else:text=_capture(output,'gtex_bundle_'+str(len(bundles)+1),url,receipts)
        bundles.append({'url':url,'text':text})
        # Direct named chunks/imports only, bounded to the same official origin.
        for relative in re.findall(r'["\']([^"\'\s]+\.js)["\']',text):
            if 'license' not in relative.lower():continue
            candidate=urllib.parse.urljoin(url,relative)
            if urllib.parse.urlparse(candidate).hostname in {'gtexportal.org','www.gtexportal.org'}:queue.append(candidate)
    _write(output/'request_receipts.json',receipts)
    return {'output':output,'texts':texts,'bundles':bundles,'receipts':receipts,'foundation_before':before}


def analyze(session):
    """Retain exact documentary text and narrow conclusions with explicit remaining gates."""
    output=session['output'];texts=session['texts'];findings={};gaps=[]
    uni=texts['uniprot_license']
    try:
        document=json.loads(uni);content=document['content']
        if not all(term in content for term in ('Creative Commons Attribution 4.0','copyrightable','patents')):
            raise ValueError('Expected official license and rights disclaimer not found')
        (output/'uniprot_license_text.md').write_text(content)
        findings['uniprot']={'status':'documented','license':'CC BY 4.0 for copyrightable database parts',
            'attribution_required':True,'other_rights_disclaimer':True,'last_modified':document.get('lastModified'),
            'source_url':URLS['uniprot_license'],'text_path':'uniprot_license_text.md'}
    except Exception as error:gaps.append('UniProt license: '+str(error))
    for key,release,assembly in (('gencode39',39,'GRCh38.p13'),('gencode47',47,'GRCh38.p14')):
        text=html.unescape(re.sub(r'<[^>]+>',' ',texts[key]));text=' '.join(text.split())
        statement=f'Release {release} ({assembly})'
        if statement in text:
            findings[key]={'status':'documented','release':release,'assembly':assembly,'statement':statement,
                'source_url':URLS[key],'interpretation':'GENCODE release identity only; source-specific GTEx linkage still requires its release metadata'}
        else:gaps.append(key+': exact release/assembly statement not found')
    candidates=[]
    for bundle in session['bundles']:
        text=bundle['text']
        for match in re.finditer(r'redistribut',text,re.I):
            start=max(0,match.start()-3500);end=min(len(text),match.end()+7500)
            window=text[start:end]
            if 'acknowledg' in window.lower():
                candidates.append({'source_url':bundle['url'],'start_character':start,'end_character':end,'raw_frontend_text':window})
    _write(output/'gtex_license_candidates.json',candidates)
    for bundle in session['bundles']:
        text=bundle['text']
        matched=re.search(r"\('<div class=\"col-sm-10\"><h2> GTEx Portal Data License.*?</div>',1\)",text,re.S)
        route=re.search(r'name:"licensePage",component:(\w+),path:"/license"',text)
        if matched is None or route is None:continue
        alias=re.search(r'var '+re.escape(route[1])+r'=(\w+);',text)
        wrapper=re.search(re.escape(alias[1])+r'=\(0,\w+\.A\)\(\w+,\[\["render",(\w+)\]\]\)',text) if alias else None
        render=re.search(r'function '+re.escape(wrapper[1])+r'\(e,a\)',text) if wrapper else None
        if render is None or not 0<matched.start()-render.start()<500:
            gaps.append('GTEx captured license text could not be associated with the /license render component');continue
        body=matched[0][2:-4]
        if not all(term in body for term in ('Users who republish or redistribute','acknowledge the GTEx Portal',
                'not indicate or imply','maintain the most current version','make known in a clear and conspicuous manner',
                'do not reflect the most current/accurate data')):
            gaps.append('GTEx component lacks an expected condition; manual review required');continue
        (output/'gtex_license_component.html').write_text(body)
        (output/'gtex_license_text.txt').write_text(html.unescape(re.sub(r'<[^>]+>',' ',body)))
        receipt=next(row for row in session['receipts'] if row['requested_url']==bundle['url'])
        findings['gtex']={'status':'documented_conditional_public_redistribution_terms',
            'source_url':bundle['url'],'page_url':URLS['gtex_license'],'source_sha256':receipt['sha256'],
            'retrieved_utc':receipt['retrieved_utc'],'source_bytes':receipt['bytes'],
            'component_byte_start':len(text[:matched.start()].encode('utf-8')),
            'component_byte_end':len(text[:matched.end()].encode('utf-8')),
            'route_component':route[1],'render_function':wrapper[1],'text_path':'gtex_license_text.txt',
            'conditions':['Acknowledge GTEx as source, including acquisition date and/or dbGaP accession',
                'Do not imply GTEx endorsement',
                'Keep distributed data current OR clearly disclose that it is not the most current/accurate data',
                'Retain supplied warranty/liability disclaimers'],
            'interpretation':'Public GTEx Portal data terms; no access or redistribution authorization inferred for controlled individual-level records',
            'downstream_pack_requirement':'Pinned v10 must carry a clear stale-version notice because v11 is publicly available; retain source/date and acknowledgements'}
        break
    if 'gtex' not in findings:gaps.append('GTEx license component text not yet recovered; JS shell does not establish redistribution permission')
    if foundation_receipts()!=session['foundation_before']:raise ValueError('Original foundation or inputs changed during review')
    summary={'findings':findings,'gaps':gaps,'verified_foundation_inputs_outputs':len(session['foundation_before']),
        'foundation_and_sources_unchanged':True,'pack_promoted':False,'biological_admission':False,'v11_promoted':False,
        'remaining':['GTEx release-to-publication lineage','Source-specific reference annotation/assembly linkage',
            'Pack attribution/stale-version notices and license implementation','Mapping/reference/graph/opening and benchmark gates'],
        'interpretation':'Documentary evidence only; original foundation admission and redistribution fields are unchanged'}
    _write(output/'summary.json',summary)
    return summary


def main():
    """Write actually executed capture/review cells and preserve every failed-site outcome."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--reuse',type=Path)
    args=parser.parse_args();started=time.monotonic();nb=ExecutedNotebook('Human foundation: bounded public terms and reference documentation')
    nb.md('Scope registered under 64.21 before requests. Verify existing foundation/source bytes, capture official GTEx license frontend text, UniProt JSON and GENCODE39/47 release pages. Each response is bounded to 8 MiB, at most 16 requests, external MemoryMax400M. Missing pages do not stop independent pages; no pack, benchmark, genome-wide or v11 promotion.')
    try:
        reuse=str(args.reuse.resolve()) if args.reuse else None
        nb.code('from scripts.review_human_reference_terms import capture, analyze',f'session=capture({str(args.out.resolve())!r}, reuse={reuse!r})',"session['receipts']")
        nb.code('summary=analyze(session)','summary')
    except Exception as error:
        nb.md('Diagnostic: '+type(error).__name__+': '+str(error));nb.write(str(args.out/'diagnostic.ipynb'));raise
    nb.write(str(args.out/'review.ipynb'))
    code=args.out/'code';code.mkdir();(code/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    (code/'notebook_runner.py').write_bytes((ROOT/'scripts/notebook_runner.py').read_bytes())
    peak=next(int(row.split()[1])*1024 for row in Path('/proc/self/status').read_text().splitlines() if row.startswith('VmHWM:'))
    _write(args.out/'runtime.json',{'seconds':time.monotonic()-started,'peak_rss_bytes':peak,
        'measurement_scope':'Linux executed-process VmHWM, excluding inherited launcher high-water marks',
        'external_memory_cap_bytes':400*1024*1024})
    _write(args.out/'manifest.json',{'files':{str(path.relative_to(args.out)):{'sha256':_sha(path),'bytes':path.stat().st_size}
        for path in args.out.rglob('*') if path.is_file()}})
    print(json.dumps(nb.ns['summary'],indent=2))


if __name__=='__main__':main()
