from pathlib import Path
import csv,json,urllib.request,urllib.error,time,concurrent.futures,hashlib,re,datetime
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'audit'
RAW=OUT/'evidence'; RAW.mkdir(parents=True,exist_ok=True)
SRC=ROOT/'bmc_reproducibility'
def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def get(url,path):
 if path.exists():return json.loads(path.read_text(encoding='utf-8'))
 for attempt in range(5):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MR-reproducibility-audit/1.0','Accept':'application/json'}),timeout=50) as r:data=json.load(r)
   path.write_text(json.dumps(data,ensure_ascii=False),encoding='utf-8');return data
  except Exception:
   if attempt==4:raise
   time.sleep(2*(attempt+1))
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Post-hoc instrument-independence sensitivity audit; not a preregistered primary analysis','source':'GWAS Catalog REST API v2 curated associations','threshold':5e-8,'lookup':'All exact rsIDs present in frozen univariable, bidirectional and multivariable instrument sets. No outcome-based instrument selection. No LD proxies in this exact-variant screen.','primary_screen':{'education_cognition':'educational|education attainment|years of schooling|intelligence|cognitive ability|cognitive performance','socioeconomic':'household income|deprivation|socioeconomic|socio-economic','alcohol':'alcohol|drinks per week|drinking behavior','physical_activity':'physical activity|sedentary|accelerometer|exercise'},'expanded_screen':{'blood_pressure':'blood pressure|hypertension','lipids':'cholesterol|triglyceride|apolipoprotein|lipoprotein','glycemia':'diabetes|glycated hemoglobin|glycated haemoglobin|fasting glucose|fasting insulin|insulin resistance'},'exclusions':'Study vascular endpoints never trigger exclusions. Smoking and adiposity phenotypes retained. Primary screen flags potential upstream/shared pathways, expanded cardiometabolic screen may also remove vertical pathways. Association alone does not establish horizontal pleiotropy or independence violation.','models':'Retain all original primary models; compare exact-rsID primary screen and expanded screen in six IVW models, six two-exposure MVMR models and two bidirectional models. Recompute FDR within each original six-outcome family and each sensitivity scheme; compute conditional F after removal.'}
pf=OUT/'independence_protocol.json'
if not pf.exists():pf.write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
files=list((SRC/'data_archive').glob('*harmonised.csv'))
targets={r['SNP'] for p in files for r in read(p) if r.get('mr_keep','TRUE').upper()=='TRUE'}
proxy_map={}
for p in (SRC/'data_archive').glob('*_outcome_associations_raw.csv'):
 for r in read(p):
  if r.get('proxy.outcome','').upper()=='TRUE':proxy_map.setdefault(r['SNP'],set()).add(r['proxy_snp.outcome'])
snps=sorted(targets|{v for vv in proxy_map.values() for v in vv})
protocol=json.loads(pf.read_text(encoding='utf-8'))
protocol['proxy_amendment']={'reason':'Source identity verification identified actual proxies retained under target rsID labels. Add all archived actual proxies to exact-variant lookup, without changing phenotype thresholds.','date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'mapping':{k:sorted(v) for k,v in proxy_map.items()},'screening_rule':'A target is flagged if it or any of its actually used proxies meets the frozen screen. Apply the same exclusion set across outcomes. No genome-wide LD expansion.'}
pf.write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'queried_rsids.txt').write_text('\n'.join(snps),encoding='utf-8')
def query(snp):
 rows=[];page=0
 while True:
  url=f'https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={snp}&size=100&page={page}'
  d=get(url,RAW/f'{snp}_{page}.json');rows+=d.get('_embedded',{}).get('associations',[])
  if page+1>=d.get('page',{}).get('totalPages',1):break
  page+=1
 return snp,rows
allrows=[];status=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 futs={pool.submit(query,s):s for s in snps}
 for i,f in enumerate(concurrent.futures.as_completed(futs),1):
  s=futs[f]
  try:
   s,rows=f.result();status.append({'SNP':s,'status':'ok','n_associations':len(rows)})
   for r in rows:
    if not any(z['rs_id']==s for z in r.get('snp_allele',[])):raise ValueError('API returned different rsID')
    allrows.append({'SNP':s,'association_id':r['association_id'],'reported_trait':'; '.join(r.get('reported_trait',[])),'efo_trait':'; '.join(z['efo_trait'] for z in r.get('efo_traits',[])),'p':r.get('p_value'),'study':r.get('accession_id'),'pmid':r.get('pubmed_id'),'effect':r.get('beta',''),'ancestry':'See linked study; Catalog curated multi-ancestry hits are not restricted to EUR'})
  except Exception as e:status.append({'SNP':s,'status':str(e),'n_associations':0})
  if i%25==0:print('queried',i,'of',len(snps),flush=True)
def write(name,rows):
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write('query_status.csv',sorted(status,key=lambda x:x['SNP']))
write('all_catalog_associations.csv',allrows)
protocol=json.loads(pf.read_text(encoding='utf-8'))
for r in allrows:
 text=(r['reported_trait']+' '+r['efo_trait']).lower();sig=r['p'] is not None and float(r['p'])<5e-8
 r['primary_flags']=';'.join(k for k,v in protocol['primary_screen'].items() if sig and re.search(v,text))
 r['expanded_flags']=';'.join(k for k,v in protocol['expanded_screen'].items() if sig and re.search(v,text))
 # Cardiovascular endpoint records are retained in the ledger, never used for selection.
 if re.search(r'stroke|coronary|peripheral artery|peripheral arterial|myocardial infarction',text):r['primary_flags']='';r['expanded_flags']=''
write('catalog_classified.csv',allrows)
flags=[]
for s in sorted(targets):
 rr=[r for r in allrows if r['SNP'] in ({s}|proxy_map.get(s,set()))];a=sorted({v for r in rr for v in r['primary_flags'].split(';') if v});b=sorted({v for r in rr for v in r['expanded_flags'].split(';') if v})
 flags.append({'SNP':s,'actual_proxies':';'.join(sorted(proxy_map.get(s,set()))),'primary_flags':';'.join(a),'expanded_flags':';'.join(b),'exclude_primary':bool(a),'exclude_expanded':bool(a or b)})
write('instrument_flags.csv',flags)
for gcst in ['GCST005840','GCST005841','GCST005842']:
 d=get('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/'+gcst,RAW/(gcst+'.json'));print(gcst,d['initial_sample_size'],flush=True)
print('DONE',len(snps),'SNPs',len(allrows),'associations',sum(r['exclude_primary'] for r in flags),'primary flags',sum(r['exclude_expanded'] for r in flags),'expanded flags',flush=True)
