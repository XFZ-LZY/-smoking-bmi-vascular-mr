from pathlib import Path
import urllib.request,csv,concurrent.futures,re,json,hashlib,sys
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'audit';EV=OUT/'evidence';SRC=ROOT/'bmc_reproducibility/data_archive'
groups=[('laa','GCST005840','GCST006907'),('svs','GCST005841','GCST006909'),('ces','GCST005842','GCST006910'),('ais','GCST005843','GCST006908')]
def reliable_lines(url):
 start=0;carry=b'';total=None
 while total is None or start<total:
  end=start+8_000_000-1
  req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'})
  with urllib.request.urlopen(req,timeout=90) as r:
   assert r.status==206,'Range not honoured'
   cr=r.headers['Content-Range'];part=r.read();total=int(cr.split('/')[-1]);assert len(part)==int(r.headers['Content-Length'])
  start+=len(part);pieces=(carry+part).split(b'\n');carry=pieces.pop()
  for line in pieces:yield line
 if carry:yield carry
def work(task):
 k,g,ancestry=task
 frozen={r['proxy_snp.outcome'] if r.get('proxy.outcome','').upper()=='TRUE' else r['SNP']:r for r in csv.DictReader(open(SRC/(k+'_outcome_associations_raw.csv'),encoding='utf-8-sig'))}
 if '--proxies' in sys.argv:frozen={s:r for s,r in frozen.items() if r.get('proxy.outcome','').upper()=='TRUE'}
 band='GCST005001-GCST006000' if g.startswith('GCST005') else 'GCST006001-GCST007000'
 base=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{band}/{g}/'
 index=urllib.request.urlopen(base,timeout=40).read().decode();(EV/(g+'_index.html')).write_text(index)
 files=re.findall('href="([^"?]+\\.out)"',index);assert len(files)==1,files
 url=base+files[0];matches=[];nbytes=0;lines=0
 if True:
  r=iter(reliable_lines(url)) if '--proxies' in sys.argv else urllib.request.urlopen(url,timeout=90)
  header=next(r).decode().strip();print(k,ancestry,'HEADER',header,flush=True)
  names=header.split()
  for line in r:
   nbytes+=len(line);lines+=1
   s=line.decode().strip().split()
   if not s:continue
   if s[0] in frozen:
    row=dict(zip(names,s));f=frozen[s[0]]
    a=row['Allele1'].upper();a2=row['Allele2'].upper();fa=f['effect_allele.outcome'].upper();fo=f['other_allele.outcome'].upper()
    if f.get('proxy.outcome','').upper()=='TRUE':
     if fa==f['target_a1.outcome'].upper():fa,fo=f['proxy_a1.outcome'].upper(),f['proxy_a2.outcome'].upper()
     else:fa,fo=f['proxy_a2.outcome'].upper(),f['proxy_a1.outcome'].upper()
    sign=1 if (a,a2)==(fa,fo) else -1 if (a2,a)==(fa,fo) else None
    b=float(row['Effect']);se=float(row['StdErr']);fb=float(f['beta.outcome']);fs=float(f['se.outcome'])
    matches.append({'key':k,'accession':g,'file_ancestry':ancestry,'source_file':url,'SNP':f['SNP'],'source_SNP':s[0],'effect_allele_source':a,'other_allele_source':a2,'effect_allele_frozen':fa,'other_allele_frozen':fo,'source_beta':b,'source_se':se,'frozen_beta':fb,'frozen_se':fs,'aligned_sign':sign,'abs_beta_diff':abs(sign*b-fb) if sign else None,'abs_se_diff':abs(se-fs)})
   if len(matches)==len(frozen):break
   if '--full' not in sys.argv and (len(matches)>=8 or nbytes>=80_000_000):break
 (EV/(g+'_fingerprint_meta.json')).write_text(json.dumps({'url':url,'header':header,'bytes_read':nbytes,'lines_read':lines,'matches':len(matches)}))
 print(k,ancestry,len(matches),'matches',nbytes,'bytes',flush=True);return matches
tasks=[(k,g,a) for k,t,e in groups for g,a in [(t,'TRANS'),(e,'EUR')] if ('--full' not in sys.argv and '--proxies' not in sys.argv) or a=='TRANS']
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for rr in pool.map(work,tasks):rows+=rr
name='stroke_proxy_verification.csv' if '--proxies' in sys.argv else 'stroke_effect_full_verification.csv' if '--full' in sys.argv else 'stroke_effect_fingerprints.csv'
with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('DONE',len(rows),flush=True)
